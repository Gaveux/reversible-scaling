module Estimator_class

    use vars_class
    use molecule_specs
    use pimc_structures

    implicit none

    type estimator
        !Variable array - stores all variables that are averaged over the course of the simulation
        type (vars), dimension(9) :: vars
        
        !Projected variables array - stores the variables that are calculated by projecting averaged quantities
        !type (project), dimension(4) :: proj

    end type estimator

    interface new
        module procedure estimator_init
    end interface

    interface update
        module procedure estimator_step
    end interface

    interface update_block
        module procedure estimator_block
    end interface

    interface end_sim
        module procedure estimator_end
    end interface

    contains
        !Initialises all of the variables that are to be estimated
        subroutine estimator_init(this)
            type (estimator), intent(out) :: this

            !Initialise the energies
            call new(this%vars(1))!, 'Energy Thermo: ')
            call new(this%vars(2))!, 'Energy Virial: ')
            call new(this%vars(3))!, 'Energy Centroid Virial: ')
            call new(this%vars(4))!, 'Kinetic Energy Thermo: ')
            call new(this%vars(5))!, 'Kinetic Energy Virial: ')
            call new(this%vars(6))!, 'Kinetic Energy Centroid Virial: ' )
            call new(this%vars(7))!, 'Potential Energy Thermo: ')
            call new(this%vars(8))!, 'Potential Energy Virial: ')
            call new(this%vars(9))!, 'Potential Energy Centroid Virial: ')

        end subroutine estimator_init

        subroutine estimator_step(sys,pimc,Beads,this)
            type (molsysdat), intent(in) :: sys
            type (pimc_par), intent(in) :: pimc
            type (pimc_particle), dimension(:), pointer :: Beads
            type (estimator), intent(inout) :: this

            real(kind=8), dimension(3) :: results

            !=====================================================================================
            !=====================================================================================
            !Block allowing for compile time or runtime selection of the action and the trial moves
            !=====================================================================================

                !If we are using the primitive action evaluate the primitive action estimators
                if(pimc%act%act_type.eq.0) then
                    call energy_thermo_pa(sys,pimc,Beads,results)
                    call update(this%vars(1),results(1))
                    call update(this%vars(4),results(2))
                    call update(this%vars(7),results(3))

                    call energy_virial_pa(sys,pimc,Beads,results)
                    call update(this%vars(2),results(1))
                    call update(this%vars(5),results(2))
                    call update(this%vars(8),results(3))

                    call energy_centvirial_pa(sys,pimc,Beads,results)
                    call update(this%vars(3),results(1))
                    call update(this%vars(6),results(2))
                    call update(this%vars(9),results(3))

                !Else if we are using the takahashi-imada action use the takahashi-imada action estimators
                else if (pimc%act%act_type.eq.1) then
                    call energy_thermo_ta(sys,pimc,Beads,results)
                    call update(this%vars(1),results(1))
                    call update(this%vars(4),results(2))
                    call update(this%vars(7),results(3))

                !Else if we are using the chin action use the chin action estimators
                else if (pimc%act%act_type.eq.2) then
                    call energy_thermo_ca(sys,pimc,Beads,results)
                    call update(this%vars(1),results(1))
                    call update(this%vars(4),results(2))
                    call update(this%vars(7),results(3))
                endif
            !#endif
            !=====================================================================================
            !=====================================================================================

        end subroutine estimator_step

        subroutine estimator_block(this)
            type (estimator), intent(inout) :: this

            integer :: i

            do i=1,9
                if(i.eq.1) then
                    write(*,*) 'Energy'
                
                else if(i.eq.4) then
                    write(*,*) 'Kinetic Energy'
                
                else if (i.eq.7) then
                    write(*,*) 'Potential Energy'
                endif
                
                call update_block(this%vars(i))
                call print_block(this%vars(i))
            enddo


        end subroutine estimator_block

        subroutine estimator_end(this)
            type (estimator), intent(inout) :: this
            integer :: i
            do i=1,9
                if(i.eq.1) then
                    write(*,*) 'Energy'
                
                else if(i.eq.4) then
                    write(*,*) 'Kinetic Energy'
                
                else if (i.eq.7) then
                    write(*,*) 'Potential Energy'

                endif
                call print_end(this%vars(i))
            enddo
             
            write(*,*) 'Total average energy:', this%vars(1)%mean_tot
        end subroutine estimator_end
        
        !==========================================================================
        !==========================================================================
        !Estimators for the primitive action - note lots of redundancy in
        !these estimators. If all estimators are to be calculated use 
        !the energy_est_pa and pressure_est_pa subroutines as they calculate
        !==========================================================================

        !Calculates the thermodynamic energy, kinetic energy and potential energy estimators
        subroutine energy_thermo_pa(sys,pimc,Beads,results)
            type (molsysdat), intent(in) :: sys
            type (pimc_par), intent(in) :: pimc
            type (pimc_particle), dimension(:), pointer :: Beads
            real(kind=8), dimension(3), intent(out) :: results

            integer :: i,j,k
            real(kind=8) :: ke_t
            real(kind=8) :: pe
            real(kind=8) :: dsq,tdsq
            real(kind=8) :: n2b
            real(kind=8) :: kin_t
            real(kind=8) :: en_t

            !initialise all different energies to zero
            ke_t=0.0d0

            pe=0.0d0

            do i=1,pimc%NumBeads
                !Initialise the terms contributing to the kinetic term for that bead to zero
                dsq=0.0d0
                !Loop over all atoms and dimension and calculate the distance squared multiplied by
                !the atom mass of Bead(i) and Bead(i+1) 
                do j=1,sys%natom
                    tdsq=0.0d0
                    do k=1,sys%dimen
                        !Compute the distance squared of neighbouring beads(thermo)
                        tdsq=tdsq+(Beads(i)%x(k,j)-Beads(i+1)%x(k,j))**2
                    enddo
                    dsq = dsq + tdsq*sys%mass(j)
                enddo

                !Update the kinetic energies
                ke_t=ke_t+dsq

                pe=pe+Beads(i)%VCurr
            enddo

            !Multiply each term by the required constants
            ke_t=ke_t*pimc%NumBeads*(0.5*pimc%invBeta**2)

            pe=pe*pimc%invNumBeads


            n2b=(0.5*dble(sys%dimen))*sys%natom*pimc%NumBeads*pimc%invBeta

            en_t=n2b-ke_t+pe
            kin_t=n2b-ke_t

            results(1) = en_t
            results(2) = kin_t
            results(3) = pe
            return

            
        end subroutine energy_thermo_pa

        !Calculates the generalised virial energy, kinetic energy and potential energy estimators
        subroutine energy_virial_pa(sys,pimc,Beads,results)
            type (molsysdat), intent(in) :: sys
            type (pimc_par), intent(in) :: pimc
            type (pimc_particle), dimension(:), pointer :: Beads
            real(kind=8), dimension(3), intent(out) :: results

            integer :: i,j,k
            real(kind=8) :: ke_v
            real(kind=8) :: pe
            real(kind=8) :: rdr_v
            real(kind=8) :: n2b
            real(kind=8) :: kin_v
            real(kind=8) :: en_v

            ke_v=0.0d0

            pe=0.0d0
            n2b=0.5*dble(sys%dimen)*sys%natom*pimc%invBeta

            do i=1,pimc%NumBeads
                !Initialise the terms contributing to the kinetic term for that bead to zero
                rdr_v=0.0d0
                do j=1,sys%natom
                    do k=1,sys%dimen
                        !Compute the corresponding virial term
                        rdr_v = rdr_v+(Beads(i)%x(k,j)-Beads(pimc%NumBeads)%x(k,j))*Beads(i)%dVdx(k,j)
                    enddo
                enddo

                !Update the kinetic energies
                ke_v = ke_v + rdr_v
                
                pe=pe+Beads(i)%VCurr
            enddo

            !Multiply each term by the required constants
            ke_v =  n2b + 0.5*ke_v*pimc%invNumBeads
            
            pe=pe*pimc%invNumBeads
            en_v=ke_v+pe
            kin_v=ke_v

            results(1) = en_v
            results(2) = kin_v
            results(3) = pe
            return

        end subroutine energy_virial_pa

        !Calculates the centroid virial energy, kinetic energy and potential energy estimators
        subroutine energy_centvirial_pa(sys,pimc,Beads,results)
            type (molsysdat), intent(in) :: sys
            type (pimc_par), intent(in) :: pimc
            type (pimc_particle), dimension(:), pointer :: Beads
            real(kind=8), dimension(3), intent(out) :: results

            integer :: i,j,k
            real(kind=8) :: ke_cv
            real(kind=8) :: pe
            real(kind=8) :: en_cv
            real(kind=8) :: rdr_cv
            real(kind=8) :: n2b
            real(kind=8) :: kin_cv

            real(kind=8), dimension(sys%dimen,sys%natom) :: cv

            !Compute the centroid
            do i=1,sys%natom
                do j=1,sys%dimen
                    cv(j,i)=0.0d0
                    do k=1,pimc%NumBeads
                        cv(j,i)=cv(j,i)+Beads(k)%x(j,i)
                    enddo
                    cv(j,i)=cv(j,i)*pimc%invNumBeads
                    
                enddo
            enddo

            ke_cv=0.0d0

            pe=0.0d0

            n2b=0.5*dble(sys%dimen)*sys%natom*pimc%invBeta

            do i=1,pimc%NumBeads
                !Initialise the terms contributing to the kinetic term for that bead to zero
                rdr_cv=0.0d0
                do j=1,sys%natom
                    do k=1,sys%dimen
                        !Compute the correspondin centroid virial term
                        rdr_cv = rdr_cv+(Beads(i)%x(k,j)-cv(k,j))*Beads(i)%dVdx(k,j)
                        
                    enddo
                enddo

                !Update the kinetic energies
                ke_cv = ke_cv + rdr_cv
                
                pe=pe+Beads(i)%VCurr
            enddo

            !Multiply each term by the required constants
            ke_cv= n2b +0.5*ke_cv*pimc%invNumBeads
            
            pe=pe*pimc%invNumBeads

            en_cv=ke_cv+pe
            kin_cv=ke_cv

            !Update the energies
            results(1)=en_cv
            results(2)=kin_cv
            results(3)=pe
            
            return

        end subroutine energy_centvirial_pa

        !==========================================================================
        !==========================================================================








        !==========================================================================
        !==========================================================================
        !Estimators for the takahashi-imada action
        !==========================================================================

        subroutine energy_thermo_ta(sys,pimc,Beads,results)
            type (molsysdat), intent(in) :: sys
            type (pimc_par), intent(in) :: pimc
            type (pimc_particle), dimension(:), pointer :: Beads
            real(kind=8), dimension(3), intent(out) :: results

            real(kind=8) :: e, ke, pe, frc, dsq, tdsq, f2, tempf2
            real(kind=8) :: n2b
            
            
            integer :: i,j,k
            
            ke = 0.0
            pe = 0.0
            frc = 0.0

            n2b=0.5*dble(sys%dimen)*sys%natom*pimc%NumBeads*pimc%invBeta
            
            do i=1,pimc%NumBeads
            
                dsq = 0.0
                f2 = 0.0
    
                do j=1,sys%natom
                    tempf2=0.0
                    tdsq=0.0
                    do k=1,sys%dimen
                        tdsq = tdsq + (Beads(i)%x(k,j)-Beads(i+1)%x(k,j))**2
                        tempf2 = tempf2 + Beads(i)%dVdx(k,j)**2

                    enddo
                    dsq = dsq + tdsq*sys%mass(j)
                    f2=f2+tempf2/sys%mass(j)

                enddo
                
                ke = ke + dsq
                pe = pe + Beads(i)%Vcurr
                frc = frc + f2

            enddo
            !print *, frc*pimc%Beta**2/(24.0*dble(pimc%NumBeads)**3), frc, pimc%NumBeads**3, pimc%Beta**2
            frc=frc*pimc%Beta**2*(pimc%invNumBeads**3)/(24.0)
            ke = (ke*pimc%NumBeads) * (0.5*pimc%invBeta**2)
            pe = pe*pimc%invNumBeads
            
            ke=n2b-ke+frc
            pe=pe+2.0*frc
            e=ke+pe
            
            !Update the energies
            results(1)=e
            results(2)=ke
            results(3)=pe
            
            return
        end subroutine energy_thermo_ta

        !==========================================================================
        !==========================================================================





        !==========================================================================
        !==========================================================================
        !Estimators for the chin action
        !==========================================================================

        subroutine energy_thermo_ca(sys,pimc,Beads,results)
            type (molsysdat), intent(in) :: sys
            type (pimc_par), intent(in) :: pimc
            type (pimc_particle), dimension(:), pointer :: Beads
            real(kind=8), dimension(3), intent(out) :: results
            
            real(kind=8) :: e, ke, pe, frc, dsq, tdsq, ttemp_frc, temp_frc
            real(kind=8) :: dsq_a, dsq_b, dsq_c, frc_a, frc_b, frc_c, pe_a, pe_b, pe_c
            real(kind=8) :: n2b

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
                    tdsq=0.0
                    ttemp_frc=0.0
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
            ke = ke*pimc%NumBeads * (0.5*pimc%invBeta**2)
            frc = frc*pimc%act%u0*pimc%Beta**2*(pimc%invNumBeads**3)

            !n2b = 3NDP/2beta for this case we have hardcoded 
            n2b = 1.5*dble(sys%dimen)*dble(sys%natom)*dble(pimc%NumBeads)*pimc%invBeta

            ke = n2b - ke + frc
            pe = pe + 2.0*frc
            e = ke + pe 
 
            !Update the energies
            results(1)=e
            results(2)=ke
            results(3)=pe
            
            return
        end subroutine energy_thermo_ca

end module Estimator_class
