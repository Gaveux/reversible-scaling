!-------------------------------------------------------------------
! calculation of action along path
!-------------------------------------------------------------------
module actions

    use molecule_specs
    use pimc_structures
    implicit none

    type action_vals
        real(kind=8) :: act, act_old
    end type action_vals

    interface eval_action
      module procedure evaluate_action
    end interface eval_action

    contains

    subroutine evaluate_action(sys,pimc,Beads,this)
        type (molsysdat), intent(in) :: sys
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this

        !=====================================================================================
        !=====================================================================================
        !Block allowing for compile time or runtime selection of the action and the trial moves
        !=====================================================================================

        if(pimc%act%act_type.eq.0) then
            if(pimc%move%move_type.eq.0) then
                call prim_action(sys,pimc,Beads,this)
            else if (pimc%move%move_type.eq.1) then
                call prim_action_staging(pimc,Beads,this)
            endif
        else if (pimc%act%act_type.eq.1) then
            if(pimc%move%move_type.eq.0) then
                call ti_action(sys,pimc,Beads,this)
            else if (pimc%move%move_type.eq.1) then
                call ti_action_staging(sys,pimc,Beads,this)
            endif
        else if (pimc%act%act_type.eq.2) then
            if(pimc%move%move_type.eq.0)then
                call chin_action(sys,pimc,Beads,this)
            else if(pimc%move%move_type.eq.1) then
                call chin_action_staging(sys,pimc,Beads,this)
            endif
        endif
        !=====================================================================================
        !=====================================================================================

    end subroutine evaluate_action

    !Primitive action
    subroutine prim_action(sys,pimc,Beads,this)

        type (molsysdat), intent(in) :: sys
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this
        real(kind=8) :: ke, pe, dsq, tdsq
        integer :: i,j,k

        !----------------------------------------------------------------------
        ! main loop to calculate the action
        !----------------------------------------------------------------------
        ke = 0.0
        pe = 0.0

        do i=1,pimc%NumBeads
            dsq = 0.0
            do j=1,sys%natom
                tdsq=0.0
                do k=1,sys%dimen
                    tdsq = tdsq + (Beads(i)%x(k,j)-Beads(i+1)%x(k,j))**2
                enddo
                dsq = dsq+tdsq*sys%mass(j)
            enddo
            ke = ke + dsq
            pe = pe + Beads(i)%Vcurr
        enddo

        ke = (ke*pimc%NumBeads) * (0.5*pimc%invBeta**2)
        pe = pe*pimc%invNumBeads
        this%act = ke + pe
        return
    end subroutine prim_action

    !The staging trial move version of the primitive action
    subroutine prim_action_staging(pimc,Beads,this)
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this
        real(kind=8) :: pe
        integer :: i

        pe = 0.0
        do i=1,pimc%NumBeads
            pe = pe + Beads(i)%Vcurr
        enddo
        pe = pe*pimc%invNumBeads
        this%act = pe

        return
    end subroutine prim_action_staging

    !The takahashi-imada action
    subroutine ti_action(sys,pimc,Beads,this)
        type (molsysdat), intent(in) :: sys
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this
        real(kind=8) :: ke, pe, frc, dsq, tdsq, f2, tempf2
        integer :: i,j,k

        ke = 0.0d0
        pe = 0.0d0
        frc = 0.0d0

        do i=1,pimc%NumBeads
            dsq = 0.0d0
            f2 = 0.0d0
            do j=1,sys%natom
                tempf2=0.0d0
                tdsq=0.0d0
                do k=1,sys%dimen
                    tdsq = tdsq + (Beads(i)%x(k,j)-Beads(i+1)%x(k,j))**2
                    tempf2 = tempf2 + Beads(i)%dVdx(k,j)**2
                enddo
                f2 = f2 + tempf2/sys%mass(j)
                dsq=dsq+tdsq*sys%mass(j)
            enddo
            ke = ke + dsq
            pe = pe + Beads(i)%Vcurr
            frc = frc + f2
        enddo
        frc=frc*(pimc%Beta**2)*(pimc%invNumBeads**3)/24.0
        ke = (ke*pimc%NumBeads)* (0.5*pimc%invBeta**2)
        pe = pe*pimc%invNumBeads
        this%act = ke + pe + frc
        return
    end subroutine ti_action

    !The staging trial move version of the takahashi-imada action
    subroutine ti_action_staging(sys,pimc,Beads,this)
        type (molsysdat), intent(in) :: sys
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this
        real(kind=8) :: pe, frc, f2, tempf2
        integer :: i,j,k

        pe = 0.0d0
        frc = 0.0d0

        do i=1,pimc%NumBeads
            f2 = 0.0d0
            do j=1,sys%natom
                tempf2=0.0d0
                do k=1,sys%dimen
                    tempf2 = tempf2 + Beads(i)%dVdx(k,j)**2
                enddo
                f2 = f2 + tempf2/sys%mass(j)
            enddo
            pe = pe + Beads(i)%Vcurr
            frc = frc + f2
        enddo

        frc=frc*pimc%Beta**2*(pimc%invNumBeads**3)/24.0
        pe = pe*pimc%invNumBeads
        this%act = pe + frc
        return
    end subroutine ti_action_staging

    !Chin action
    subroutine chin_action(sys,pimc,Beads,this)
        type (molsysdat), intent(in) :: sys
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this
        real(kind=8) :: ke, pe, frc, dsq, tdsq, temp_frc, ttemp_frc
        real(kind=8) :: dsq_a, dsq_b, dsq_c, frc_a, frc_b, frc_c, pe_a, pe_b, pe_c
        integer :: i,j,k

        ke = 0.0
        pe = 0.0
        frc = 0.0
        dsq_a=0.0
        dsq_b=0.0
        dsq_c=0.0
        frc_a=0.0
        frc_b=0.0
        frc_c=0.0
        pe_a=0.0
        pe_b=0.0
        pe_c=0.0

        do i=0,pimc%NumBeads-1
            dsq = 0.0
            temp_frc=0.0

            do j=1,sys%natom
                tdsq=0.0d0
                ttemp_frc=0.0d0
                do k=1,sys%dimen
                    !index as given because fortran is one indexed
                    dsq_a = pimc%act%t1inv*(Beads(3*i+1)%x(k,j)-Beads(3*i+2)%x(k,j))**2
                    dsq_b = pimc%act%t1inv*(Beads(3*i+2)%x(k,j)-Beads(3*i+3)%x(k,j))**2
                    dsq_c = 0.5*pimc%act%t0inv*(Beads(3*i+3)%x(k,j)-Beads(3*i+4)%x(k,j))**2

                    tdsq = tdsq + (dsq_a+dsq_b+dsq_c)

                    frc_a = pimc%act%a1*Beads(3*i+1)%dVdx(k,j)**2
                    frc_b= (1.0-2.0*pimc%act%a1)*Beads(3*i+2)%dVdx(k,j)**2
                    frc_c= pimc%act%a1*Beads(3*i+3)%dVdx(k,j)**2

                    ttemp_frc = ttemp_frc + (frc_a+frc_b+frc_c)

                enddo
                dsq = dsq + tdsq*sys%mass(j)
                temp_frc = temp_frc + ttemp_frc/sys%mass(j)
            enddo

            ke = ke + dsq
            frc = frc + temp_frc
            pe_a=Beads(3*i+1)%VCurr*pimc%act%v1
            pe_b=Beads(3*i+2)%VCurr*pimc%act%v2
            pe_c=Beads(3*i+3)%VCurr*pimc%act%v1
            pe = pe + pe_a + pe_b + pe_c

        enddo

        pe = pe*pimc%invNumBeads
        ke = ke*pimc%NumBeads*0.5
        frc = frc*pimc%act%u0*(pimc%invNumBeads**3)

        this%act = pe + frc*pimc%Beta*pimc%Beta + ke*(pimc%invBeta**2)
    end subroutine chin_action

    !Staging version of the Chin action
    subroutine chin_action_staging(sys,pimc,Beads,this)
        type (molsysdat), intent(in) :: sys
        type (pimc_par), intent(in) :: pimc
        type (pimc_particle), dimension(:), pointer :: Beads
        type (action_vals) :: this
        real(kind=8) :: pe, frc, temp_frc, ttemp_frc, frc_a, frc_b, frc_c, pe_a, pe_b, pe_c
        integer :: i,j,k

        pe = 0.0
        frc = 0.0
        pe_a=0.0
        pe_b=0.0
        pe_c=0.0

        do i=0,pimc%NumBeads-1
            temp_frc=0.0
            do j=1,sys%natom
                ttemp_frc=0.0d0
                do k=1,sys%dimen
                    !index as given because fortran is one indexed
                    frc_a = pimc%act%a1*Beads(3*i+1)%dVdx(k,j)**2
                    frc_b= (1.0-2.0*pimc%act%a1)*Beads(3*i+2)%dVdx(k,j)**2
                    frc_c= pimc%act%a1*Beads(3*i+3)%dVdx(k,j)**2

                    ttemp_frc = ttemp_frc + (frc_a+frc_b+frc_c)
                enddo
                temp_frc = temp_frc + ttemp_frc/sys%mass(j)
            enddo

            frc = frc + temp_frc
            pe_a=Beads(3*i+1)%VCurr*pimc%act%v1
            pe_b=Beads(3*i+2)%VCurr*pimc%act%v2
            pe_c=Beads(3*i+3)%VCurr*pimc%act%v1
            pe = pe + pe_a + pe_b + pe_c
        enddo

        pe = pe*pimc%invNumBeads
        frc = frc*pimc%act%u0*(pimc%invNumBeads**3)

        this%act = pe + frc*pimc%Beta*pimc%Beta
    end subroutine chin_action_staging

end module actions
