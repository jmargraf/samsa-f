module module_dscf
  implicit none
  double precision, allocatable                 :: Fock0(:,:,:)

contains

!#############################################
!#       Calculate fractional SCF curve
!#############################################
subroutine calc_dfrac(doprint)
  use module_data,          only : Occ,nOccA,nOccB,dim_1e,Eps,Spins,DropMO,Damp,Fock
  use module_energy,        only : Etot,calc_Embpt2,Embpt2f,Embpt2
  use module_energy,        only : E_1f
  use module_energy,        only : E_AAxf,E_AAcf,E_BBxf,E_BBcf
  use module_energy,        only : E_OSf
  use module_scf,           only : run_SCF
  use module_trans,         only : trans_full
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i, iSpin
  double precision            :: E2ex, E2c, ESS, EOS
!  E0 = Etot

  allocate(Fock0(dim_1e,dim_1e,spins),)
  Fock0 = Fock

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating fractional occupation curve HOMO (Ha)"
    write(*,*) "    "
    write(*,*) "    alpha channel:"
    write(*,*) &
"                Occ         E(SCF)         E(MBPT2)       E2ex         E2c          E2SS          E2OS          Esing"
  endif

! Only uses Alpha Channel
! IPs:
  do i=0,20
    Fock = Fock0
!    Damp=0.05
    Occ(nOccA,1) = Occ(nOccA,1)-0.05d0*dble(i)
    call run_scf(.false.)
    call trans_full(.false.)
    call calc_Embpt2()
    E2ex = E_AAxf + E_BBxf
    E2c  = E_OSf + E_AAcf + E_BBcf
    ESS = E_AAxf + E_BBxf + E_AAcf + E_BBcf
    EOS = E_OSf
    if(doprint)then
      write(*,'("    Frac: ",8(" ",F12.5," "))') Occ(nOccA,1), Etot, Embpt2f, E2ex, E2c, ESS, EOS,E_1f
    endif
    Occ(nOccA,1) = 1.0d0
  enddo

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating fractional occupation curve LUMO (Ha)"
    write(*,*) "    "
    write(*,*) "    alpha channel:"
    write(*,*) &
"                Occ         E(SCF)         E(MBPT2)       E2ex         E2c           E2SS          E2OS          Esing "
  endif

! Only uses alpha Channel
! EAs:
  do i=0,20
    Fock = Fock0
!    Damp=0.05
    Occ(nOccA+1,1) = Occ(nOccA+1,1)+0.05d0*dble(i)
    call run_scf(.false.)
    call trans_full(.false.)
    call calc_Embpt2()
    E2ex = E_AAxf + E_BBxf
    E2c  = E_OSf + E_AAcf + E_BBcf
    ESS = E_AAxf + E_BBxf + E_AAcf + E_BBcf
    EOS = E_OSf
    if(doprint)then
      write(*,'("    Frac: ",8(" ",F12.5," "))') Occ(nOccA+1,1), Etot, Embpt2f, E2ex, E2c, ESS, EOS,E_1f
    endif
    Occ(nOccA+1,1) = 0.0d0
  enddo


end subroutine calc_dfrac



!#############################################
!#       Calculate dSCF EV corrections
!#############################################
subroutine calc_dSCF(doprint)
  use module_data,          only : Occ,nOccA,nOccB,dim_1e,Eps,Spins,DropMO,Damp,Fock
  use module_energy,        only : Etot,calc_Embpt2,Embpt2f,Embpt2
  use module_scf,           only : run_SCF
  use module_trans,         only : trans_full
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i, iSpin
  double precision            :: corrEV, E0, IPdSCF, relax,E02
  double precision,allocatable:: Eps0(:,:)

  E0 = Etot

  allocate(Fock0(dim_1e,dim_1e,spins),Eps0(dim_1e,Spins))
  Fock0 = Fock
  Eps0  = Eps
  E02   = Embpt2

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating dSCF orbital relaxation corrections to eigenvalues (Ha)"
    write(*,*) "    "
    write(*,*) "    alpha channel:"
    write(*,*) "         iOrb      E+(SCF)   E0(SCF)   Eps(SCF)   IP(dSCF)   relax   E+(2)   E0(2)"
  endif

! Only uses Alpha Channel
! IPs:
  do i=nOccA,nOccA,-1
    Fock = Fock0
    Damp=0.05
    Occ(i,1) = 0.0d0
    call run_scf(.false.)
    call trans_full(.false.)
    call calc_Embpt2()
    IPdSCF = E0 - Etot
    relax  = IPdSCF - Eps0(i,1) 
    if(doprint)then
      write(*,'("    IP: ",I5,7(" ",F12.5," "))') i, Etot, E0, Eps0(i,1), IPdSCF,  relax, Embpt2f, E02
    endif
    Occ(i,1) = 1.0d0
  enddo

!  if(doprint)then
!    write(*,*) "        ***********************************************"
!  endif

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

!#############################################
!#       Calculate fractional CC curve
!#############################################
subroutine calc_fracCC(doprint)
  use module_data,          only : Occ,nOccA,nOccB,dim_1e,Eps,Spins,DropMO,Damp,Fock
  use module_energy,        only : Etot,calc_Embpt2,Embpt2f,Embpt2
  use module_scf,           only : run_SCF
  use module_trans,         only : trans_full,trans_ucc
  use module_cc,            only : Eccd,calc_Eccd
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i, iSpin
  double precision            :: E2ex, E2c, ESS, EOS
!  E0 = Etot

  allocate(Fock0(dim_1e,dim_1e,spins),)
  Fock0 = Fock

!  if(doprint)then
!    write(*,*) "    "
!    write(*,*) "    Calculating fractional occupation curve HOMO (Ha)"
!    write(*,*) "    "
!    write(*,*) "    alpha channel:"
!    write(*,*) &
!"                Occ         E(SCF)         E(MBPT2)       ELCCD  "
!  endif

! Only uses Alpha Channel
! IPs:
!  do i=0,20
!    Fock = Fock0
!    Occ(nOccA,1) = Occ(nOccA,1)-0.05d0*dble(i)
!    call run_scf(.false.)
!    call trans_full(.false.)
!    call trans_ucc(.false.)
!    call calc_Eccd(.false.,.true.)
!    call calc_Embpt2()

!    if(doprint)then
!      write(*,'("    Frac: ",4(" ",F12.5," "))') Occ(nOccA,1), Etot, Embpt2f, Eccd
!    endif
!    Occ(nOccA,1) = 1.0d0
!  enddo
!  if(doprint)then
!    write(*,*) "    "
!    write(*,*) "    Calculating fractional occupation curve LUMO (Ha)"
!    write(*,*) "    "
!    write(*,*) "    alpha channel:"
!    write(*,*) &
!"                Occ         E(SCF)         E(MBPT2)       ELCCD  "
!  endif

! Only uses alpha Channel
! EAs:
  do i=0,20
    Fock = Fock0
    Occ(nOccA+1,1) = Occ(nOccA+1,1)+0.05d0*dble(i)
    call run_scf(.false.)
    call trans_full(.false.)
    call trans_ucc(.false.)
    call calc_Eccd(.false.,.true.)
    call calc_Embpt2()

    if(doprint)then
      write(*,'("    Frac: ",4(" ",F12.5," "))') Occ(nOccA+1,1), Etot, Embpt2f, Eccd
    endif
    Occ(nOccA+1,1) = 0.0d0
  enddo


end subroutine calc_fracCC

end module module_dscf
