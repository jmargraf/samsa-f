module module_scf
  implicit none

! Convergence        
  double precision                            :: dE       = 0.0d0
  logical                                     :: SCFconv  = .false. 

contains

!#############################################
!#              Run SCF Loop        
!#############################################
subroutine run_SCF(DoPrint)
  use module_data,                      only : MaxSCF,Damp,doSCGKT,dim_1e,Fock
  use module_energy,                    only : Eold,calc_Energy,Etot,Eelec
  use module_wavefun,                   only : ort_Fock,dia_Fock,deort_Fock
  use module_wavefun,                   only : calc_Dens,calc_Fock
  use module_wavefun,                   only : Drmsd
  use module_wavefun,                   only : Fock_to_MO,Fock_to_AO,Semicanonical
  use module_occupy,                    only : run_corrF,calc_GKT
  use module_trans,                     only : trans_full
  use module_io,                        only : print_Mat
  implicit none
  integer                  :: iSCF
  real                     :: time, starttime
  logical, intent(IN)      :: DoPrint

  starttime=secnds(0.0)

  if(DoPrint)then
    write(*,*) ""
    write(*,*) "    starting SCF ..."
    write(*,*) ""
    write(*,*) "          #           Etot                 dE                  Drmsd      Damp"
  endif

  SCFconv  = .false.
  dE       = 0.0d0

  do iSCF=1,MaxSCF
    Eold = Etot

    call ort_Fock()

    call dia_Fock()

    call deort_Fock()

    call calc_Dens()

    call calc_Energy()

!    call calc_DynDamp()

    dE = Etot - Eold     

    !  iSCF   Etot  dE   Drmsd
    if(DoPrint)then
      write(*,'("  SCF  ",I5,"  ",3(F18.10,"  "),F5.2)') iSCF, Etot, dE, Drmsd, Damp
    endif

    call check_Conv(DoPrint)

    if(SCFconv .or. (iSCF==MaxSCF-1) ) exit

    call calc_Fock()

    if(doSCGKT .and. abs(dE)<1.0d-10)then
      call Fock_to_MO()
      call SemiCanonical()
      call Fock_to_AO()

      call trans_full(.false.) 
      call calc_GKT(.false.)
      call run_corrF(.true.)

    endif
  enddo

  if(.not.SCFconv)then
    write(*,*) "  (!) Warning: The SCF seems to be lost in a kafkaesque nightmare and has NOT converged (!)"
  endif

  if(DoPrint)then
    write(*,*) "         Eelec =",Eelec
    write(*,*) "    Final Etot =",Etot
  endif

  time=secnds(starttime)
  if(DoPrint)then
    write(*,*) "      SCF done in ",time," s"
  endif

  call Fock_to_MO()
!  call print_Mat(Fock(:,:,1),dim_1e,7,"MO Fock")
!  call print_Mat(Fock(:,:,2),dim_1e,7,"MO Fock")
  call SemiCanonical()
!  call print_Mat(Fock(:,:,1),dim_1e,7,"MO Fock")
!  call print_Mat(Fock(:,:,2),dim_1e,7,"MO Fock")
  call Fock_to_AO()

end subroutine run_SCF

!#############################################
!#          Recalculate WF non-SCF    
!#############################################
subroutine recalc_WF()
  use module_data,                      only : MaxSCF 
  use module_energy,                    only : Eold,calc_Energy,Etot,Eelec
  use module_wavefun,                   only : ort_Fock,dia_Fock,deort_Fock
  use module_wavefun,                   only : calc_Dens,calc_Fock
  implicit none

  write(*,*) ""
  write(*,*) "    non-SCF recalculation of WF ..."
  write(*,*) ""

  Eold = Etot

  call calc_Dens()

  call calc_Fock()

  call ort_Fock()

  call dia_Fock()

  call deort_Fock()

  call calc_Dens()

  call calc_Energy()

  dE = Etot - Eold

  !  iSCF   Etot  dE   Drmsd
  write(*,'("    Modified WF:  ","  ",2(F18.10,"  "))') Etot, dE 

! write(*,*) "         Eelec=",Eelec
! write(*,*) "    Final Etot=",Etot

end subroutine recalc_WF


!#############################################
!#             Check Convergence       
!#############################################
subroutine check_Conv(DoPrint)
  use module_data,      only : Econv 
  use module_data,      only : Dconv 
  use module_data,      only : DoDamp, Damp,DampTol
  use module_energy,    only : Eelec, Etot
  use module_wavefun,   only : Drmsd
  implicit none
  logical, intent(IN)       :: DoPrint

  if(abs(dE) < 1.0d-2 .and. Damp /= 1.0d0  & 
                      .and. Damp /= 0.7d0  &
                      .and. DoDamp)then
    Damp = 0.7d0
    if(DoPrint)then
      write(*,*) "    Damping set to 0.7"
    endif
  endif

  if(abs(dE) < DampTol .and. Damp /= 1.0d0  & 
                      .and. DoDamp)then
    Damp = 1.0d0
    DoDamp =.false.
    if(DoPrint)then
      write(*,*) "    Damping turned off"
    endif
  endif

  if(Drmsd<Dconv .and. abs(dE)<Econv)then
    if(DoPrint)then
      write(*,*) " "
      write(*,*) "    SCF converged (hurra)!"
      write(*,*) " "
    endif
!    write(*,*) "         Eelec=",Eelec
!    write(*,*) "    Final Etot=",Etot
    SCFconv = .true.
  endif

end subroutine check_Conv

end module module_scf

