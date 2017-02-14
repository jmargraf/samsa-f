module module_dscf
  implicit none

contains

!#############################################
!#       Calculate dSCF EV corrections
!#############################################
subroutine calc_dSCF(doprint)
  use module_data,          only : Occ,nOccA,nOccB,dim_1e,Eps,Spins,DropMO
  use module_energy,        only : Etot
  use module_scf,           only : run_SCF
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i, iSpin
  double precision            :: corrEV, E0, dE01

  E0 = Etot

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating dSCF orbital relaxation corrections to eigenvalues (Ha)"
    write(*,*) "    "
    write(*,*) "    alpha channel:"
    write(*,*) "         iOrb      dE(SCF)      Eps(HF)       Eps(dSCF)  "
  endif

! Only uses Alpha Channel
! IPs:
  do i=nOccA,nOccA
    Occ(i,1) = 0.0d0
    call run_scf(.true.)
    dE01 = E0 - Etot
    corrEV = Eps(i,1) + dE01
    if(doprint)then
      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dE01, Eps(i,1), corrEV
    endif
    Occ(i,1) = 1.0d0
  enddo

  if(doprint)then
    write(*,*) "        ***********************************************"
  endif

! EAs:
!  do i=nOccA+1,dim_1e
!    Occ(i,1) = 1.0d0
!    call calc_Embpt2()
!    dEV(i,1) = Embpt2 - Embpt2f
!    corrEV = Eps(i,1) + dEV(i,1)
!    if(doprint)then
!      write(*,'("    MO: ",I5,3(" ",F12.5," "))') i, dEV(i,1), Eps(i,1), corrEV
!    endif
!    Occ(i,1) = 0.0d0
!  enddo

end subroutine calc_dSCF

end module module_dscf
