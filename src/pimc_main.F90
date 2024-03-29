
program pimc90  

    !> @file
    !Description:
    !!@details
    !>program pimc90: The core path integral monte carlo program. Sets up the core data structures
    !>                from input files specified as command line arguments and runs the core monte
    !>                carlo subroutine

    use molecule_specs
    use pimc_structures
    use path_integral_monte_carlo

#if POT == 0
    use potential_msi
#elif POT == 1
    use potential_h2o
#elif POT == 2
    use potential_nh3
#elif POT == 3
    use potential_hcn
#endif

    implicit none

    type (molsysdat) :: sys 
    type (pimc_par) :: pimc
    type (pimc_particle), dimension(:), pointer :: Beads
    type (pimc_particle), dimension(:), pointer :: OldBeads

    real :: t1,t2,hours,minutes,seconds
    integer :: i,iargs
 
    character(len=80) :: OUT_DIRECTORY, IN_PIMC, IN_SYSTEM, IN_ISEED, IN_BINNING

#if POT == 0
    type (msi_params) :: pot
    character(len=80) :: POT_FILE, IN_INTERP, IN_ATOMPERM
#endif

    include 'pimc_setup.int'

    !Get the number of command line arguments
    iargs = IARGC()
  
    !If the number of command line arguments specified is not the number required
    !quit execution and print an error message

#if POT == 0
    if(iargs.ne.8) then
        write(*,*) 'Error: incorrect number of command line arguments'
        write(*,*) iargs
        write(*,*) 'Correct Usage:'
        write(*,*) '    pimc90 pimc.in system.in iseed.in binning.in interp.in pot.in atomperms.in output_dir/'
        stop
    endif
#else

    if(iargs.ne.5) then
        write(*,*) 'Error: incorrect number of command line arguments'
        write(*,*) iargs
        write(*,*) 'Correct Usage:'
        write(*,*) '    pimc90 pimc.in system.in iseed.in binning.in output_dir/'
        stop
    endif
#endif

    !read the command line arguments in

    call getarg(1,IN_PIMC)
    call getarg(2,IN_SYSTEM)
    call getarg(3,IN_ISEED)
    call getarg(4,IN_BINNING)

#if POT == 0
    call getarg(5,IN_INTERP)
    call getarg(6,POT_FILE)
    call getarg(7,IN_ATOMPERM)
    call getarg(8,OUT_DIRECTORY)
#else
    call getarg(5,OUT_DIRECTORY)
#endif

    !print which type of potential energy surface is being used

#ifdef FREE_ENERGY
    write(*,*) "Free Energy Simulation Mode is On"
#endif

#if POT == 0
    write(*,*) "Modified Shepard Potential" 
#elif POT == 1
    write(*,*) "PJT2 Water Potential" 
#elif POT == 2
    write(*,*) "AMMPOT4 Ammonia Potential" 
#elif POT == 3
    write(*,*) "Murrell, Carter, Halonen Hydrogen Cyanide Potential" 
#endif

    !---------------------------------------------------------------
    ! Start timing CPU time used
    !---------------------------------------------------------------
  
    call cpu_time(t1)

    !make the output directory
    call system('mkdir -p '//trim(OUT_DIRECTORY)) 

    !read the input files
    call read_system_data(sys,IN_SYSTEM)
    call read_iseed(sys,IN_ISEED)
    call read_pimc(sys,pimc,IN_PIMC)
    
    !setup the simulation from the read in values
    call pimc_setup(sys,pimc,Beads,OldBeads)
 
#if POT == 0
    !Initialise the modified shepard potential energy surface
    call MSI_INIT(pot,sys,IN_INTERP,POT_FILE,IN_ATOMPERM,pimc%numBeadsEff)
#endif


    ! calculate potential energies for the initial geometry
    do i=1,pimc%NumBeadsEff
#if POT == 0
      !MSI potential energy surfaces
      call potential(i,pot,Beads(i)%x,Beads(i)%r,Beads(i)%VCurr,Beads(i)%dVdx)
#else
      !Analytic potential energy surfaces
      call potential(sys,Beads(i)%x,Beads(i)%r,Beads(i)%VCurr,Beads(i)%dVdx)
#endif
    enddo
 
    ! copy the initial geometry to the end of the array (NumGeoms + 1)
    ! to make summing the action over the whole closed system easier
 
    call copy(Beads(1),Beads(pimc%NumBeads+1))
    
    ! main routine - move beads, evaluate action and monte carlo step
    !Create file target for writing variables at end of simulation
#if POT == 0
    call pimc_monte(sys,pimc,Beads,OldBeads,pot, OUT_DIRECTORY, IN_BINNING)
#else 
    call pimc_monte(sys,pimc,Beads,OldBeads,OUT_DIRECTORY, IN_BINNING)
#endif
    ! stop timer and print CPU time used
    call cpu_time(t2)
    print *,' '
    hours = (t2-t1)/3600.0
    minutes = (hours - int(hours))*60.0
    seconds =  (minutes - int(minutes))*60.0
    print *,''//achar(27)//'[31m CPU time: '//achar(27)//'[0m',int(hours),&
    &       'hour(s)',int(minutes),'minute(s)',int(seconds),'second(s)'

end
