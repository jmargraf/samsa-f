module module_occupy
  implicit none

! Energies           

contains

!#############################################
!#       Calculate MP2 EV corrections
!#############################################
subroutine calc_EVPT2()
  use module_data, only      : Occ,nOccA,nOccB,dim_1e,Eps
  use module_energy, only    : calc_Emp2,Emp2,Emp2f 
  use module_energy, only    : calc_Edcpt2,Edcpt2,Edcpt2f
  implicit none
  integer                   :: i
  double precision          :: dEV,corrEV

  write(*,*) "    "
  write(*,*) "    Calculating MBPT(2) correlation corrections to eigenvalues (Ha)"
  write(*,*) "    "
  write(*,*) "         iOrb      dE(PT2)      Eps(HF)       Eps(PT2)  "
!                 MO:     1      0.00120     -11.20700     -11.20580 

! Only uses Alpha Channel
! IPs:
  do i=1,nOccA
    Occ(i,1) = 0.0d0
    call calc_Emp2()
    dEV = Emp2 - Emp2f
    corrEV = Eps(i,1) + dEV
    write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV, Eps(i,1), corrEV
    Occ(i,1) = 1.0d0
  enddo

  write(*,*) "        ***********************************************"

! EAs:
  do i=nOccA+1,dim_1e
    Occ(i,1) = 1.0d0
    call calc_Emp2()
    dEV = Emp2 - Emp2f
    corrEV = Eps(i,1) + dEV
    write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV, Eps(i,1), corrEV
    Occ(i,1) = 0.0d0
  enddo


  write(*,*) "    "
  write(*,*) "    Calculating DCPT(2) correlation corrections to eigenvalues (Ha)"
  write(*,*) "    "
  write(*,*) "         iOrb      dE(PT2)      Eps(HF)       Eps(PT2)  "

! Only uses Alpha Channel
! IPs:
  do i=1,nOccA
    Occ(i,1) = 0.0d0
    call calc_Edcpt2()
    dEV = Edcpt2 - Edcpt2f
    corrEV = Eps(i,1) + dEV
    write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV, Eps(i,1), corrEV
    Occ(i,1) = 1.0d0
  enddo

  write(*,*) "        ***********************************************"

! EAs:
  do i=nOccA+1,dim_1e
    Occ(i,1) = 1.0d0
    call calc_Edcpt2()
    dEV = Edcpt2 - Edcpt2f
    corrEV = Eps(i,1) + dEV
    write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV, Eps(i,1), corrEV
    Occ(i,1) = 0.0d0
  enddo


end subroutine calc_EVPT2


end module module_occupy


