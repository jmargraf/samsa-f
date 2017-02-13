module module_occupy
  implicit none
  double precision, allocatable :: dEV(:,:)
! Energies           

contains

!#############################################
!#       Calculate MBPT2 EV corrections
!#############################################
subroutine calc_GKT(doprint)
  use module_data,          only : doMBPT2,doDCPT2
  use module_data,          only : Occ,nOccA,nOccB,dim_1e,Eps,Spins
  use module_energy,        only : calc_Embpt2,Embpt2,Embpt2f 
  use module_energy,        only : calc_Edcpt2,Edcpt2,Edcpt2f
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i, iSpin
  double precision            :: corrEV

  if(.not.allocated(dEV))then
    allocate(dEV(dim_1e,Spins))
  endif

  if(doMBPT2)then
  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating MBPT(2) correlation corrections to eigenvalues (Ha)"
    write(*,*) "    "
    write(*,*) "    alpha channel:"
    write(*,*) "         iOrb      dE(PT2)      Eps(HF)       Eps(PT2)  "
  endif

! Only uses Alpha Channel
! IPs:
  do i=1,nOccA
    Occ(i,1) = 0.0d0
    call calc_Embpt2()
    dEV(i,1) = Embpt2 - Embpt2f
    corrEV = Eps(i,1) + dEV(i,1)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,1), Eps(i,1), corrEV
    endif
    Occ(i,1) = 1.0d0
  enddo

  if(doprint)then
    write(*,*) "        ***********************************************"
  endif

! EAs:
  do i=nOccA+1,dim_1e
    Occ(i,1) = 1.0d0
    call calc_Embpt2()
    dEV(i,1) = Embpt2 - Embpt2f
    corrEV = Eps(i,1) + dEV(i,1)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,1), Eps(i,1), corrEV
    endif
    Occ(i,1) = 0.0d0
  enddo

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    beta channel:"
    write(*,*) "         iOrb      dE(PT2)      Eps(HF)       Eps(PT2)  "
  endif

! Only uses Beta Channel
! IPs:
  do i=1,nOccB
    Occ(i,2) = 0.0d0
    call calc_Embpt2()
    dEV(i,2) = Embpt2 - Embpt2f
    corrEV = Eps(i,2) + dEV(i,2)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,2), Eps(i,2), corrEV
    endif
    Occ(i,2) = 1.0d0
  enddo

  if(doprint)then
    write(*,*) "        ***********************************************"
  endif

! EAs:
  do i=nOccB+1,dim_1e
    Occ(i,2) = 1.0d0
    call calc_Embpt2()
    dEV(i,2) = Embpt2 - Embpt2f
    corrEV = Eps(i,2) + dEV(i,2)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,2), Eps(i,2), corrEV
    endif
    Occ(i,2) = 0.0d0
  enddo
  endif

  if(doDCPT2)then
  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating DCPT(2) correlation corrections to eigenvalues (Ha)"
    write(*,*) "    "
    write(*,*) "    alpha channel:"
    write(*,*) "         iOrb      dE(PT2)      Eps(HF)       Eps(PT2)  "
  endif

! Only uses Alpha Channel
! IPs:
  do i=1,nOccA
    Occ(i,1) = 0.0d0
    call calc_Edcpt2()
    dEV(i,1) = Edcpt2 - Edcpt2f
    corrEV = Eps(i,1) + dEV(i,1)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,1), Eps(i,1), corrEV
    endif
    Occ(i,1) = 1.0d0
  enddo

  if(doprint)then
    write(*,*) "        ***********************************************"
  endif

! EAs:
  do i=nOccA+1,dim_1e
    Occ(i,1) = 1.0d0
    call calc_Edcpt2()
    dEV(i,1) = Edcpt2 - Edcpt2f
    corrEV = Eps(i,1) + dEV(i,1)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,1), Eps(i,1), corrEV
    endif
    Occ(i,1) = 0.0d0
  enddo

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    beta channel:"
    write(*,*) "         iOrb      dE(PT2)      Eps(HF)       Eps(PT2)  "
  endif

! Only uses Beta Channel
! IPs:
  do i=1,nOccA
    Occ(i,2) = 0.0d0
    call calc_Edcpt2()
    dEV(i,2) = Edcpt2 - Edcpt2f
    corrEV = Eps(i,2) + dEV(i,2)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,2), Eps(i,2), corrEV
    endif
    Occ(i,2) = 1.0d0
  enddo

  if(doprint)then
    write(*,*) "        ***********************************************"
  endif

! EAs:
  do i=nOccA+1,dim_1e
    Occ(i,2) = 1.0d0
    call calc_Edcpt2()
    dEV(i,2) = Edcpt2 - Edcpt2f
    corrEV = Eps(i,2) + dEV(i,2)
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,2), Eps(i,2), corrEV
    endif
    Occ(i,2) = 0.0d0
  enddo
  endif

end subroutine calc_GKT

!#############################################
!#       Correlation Correction to Fock
!#############################################
subroutine run_corrF(doprint)
  use module_data,          only : Occ,nOccA,nOccB,dim_1e,Eps,Spins,Fock
  use module_energy,        only : calc_Energy, Etot,E_1
  use module_wavefun,       only : Fock_to_MO,Fock_to_AO
  use module_io,            only : print_Mat
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i,j,iSpin
  double precision            :: corrEV

!  call print_Mat(Fock(:,:,1),dim_1e,7,"AO Fock")

  call Fock_to_MO()

  do iSpin=1,Spins
    do i=1,dim_1e
      if(Occ(i,iSpin)==0.0d0)cycle
      Fock(i,i,iSpin) = Fock(i,i,iSpin) + dEV(i,iSpin)*Occ(i,iSpin)
    enddo
  enddo

  call Fock_to_AO()

  call calc_Energy()

  if(doprint)then
    write(*,*) "    Etot, E_1 = ",Etot,E_1
  endif

end subroutine run_corrF


end module module_occupy


