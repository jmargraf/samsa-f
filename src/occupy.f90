module module_occupy
  implicit none

! Energies           
!  double precision                            :: Emp2frac = 0.0d0

contains

!#############################################
!#       Calculate MP2 EV corrections
!#############################################
subroutine calc_EVPT2()
  use module_data, only      : Occ,nOccA,nOccB,dim_1e,Eps
  use module_energy, only    : calc_Emp2,Emp2,Emp2frac 
  implicit none
  integer                   :: i
  double precision          :: dEV,corrEV

  write(*,*) "    "
  write(*,*) "    Calculating MBPT(2) correlation corrections to eigenvalues"
  write(*,*) "    "


! Only uses Alpha Channel
! IPs:
  do i=1,nOccA
    Occ(i,1) = 0.0d0
    call calc_Emp2()
    dEV = Emp2 - Emp2frac
    corrEV = Eps(i,1) + dEV
    write(*,*) "MO:  ", i, dEV, Eps(i,1), corrEV
    Occ(i,1) = 1.0d0
  enddo

  write(*,*) "****************"

! EAs:
  do i=nOccA+1,dim_1e
    Occ(i,1) = 1.0d0
    call calc_Emp2()
    dEV = Emp2 - Emp2frac
    corrEV = Eps(i,1) + dEV
    write(*,*) "MO:  ", i, dEV, Eps(i,1), corrEV
    Occ(i,1) = 0.0d0
  enddo

end subroutine calc_EVPT2


end module module_occupy


