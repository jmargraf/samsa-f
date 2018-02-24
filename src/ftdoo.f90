module module_ftdoo
  implicit none

contains

!#############################################
!#   Calculate T dependence of MBPT energy
!#############################################
subroutine calc_ftdmp2(doprint)
  use module_data,          only : Tel, doFTSCF, doGKSCF,ResetOcc
  use module_energy,        only : Etot,calc_Embpt2,Embpt2
  use module_energy,        only : E_1
  use module_scf,           only : run_SCF
  use module_trans,         only : trans_full
  implicit none
  logical, intent(IN)         :: doprint
  integer                     :: i, iSpin
  double precision            :: E2ex, E2c, ESS, EOS
  double precision            :: Tlist(1:22) = (/ &
                                                                  1.0d0, &
                                                                  1.0d1, &
                                                                  1000.0d0, &
                                                                  10000.0d0, &
                                                                  20000.0d0, &
                                                                  30000.0d0, &
                                                                  40000.0d0, &
                                                                  45000.0d0, &
                                                                  50000.0d0, &             
                                                                  55000.0d0, &
                                                                  60000.0d0, &
                                                                  65000.0d0, &
                                                                  70000.0d0, &
                                                                  75000.0d0, &
                                                                  80000.0d0, &
                                                                  85000.0d0, &
                                                                  90000.0d0, &
                                                                  95000.0d0, &
                                                                 100000.0d0, &
                                                                 105000.0d0, &
                                                                 110000.0d0, &
                                                                 115000.0d0  &
                                                                             /)
  double precision            :: Ezero,Eplus,Eminus,dEdT
  double precision            :: Tstep = 10.0d0
  double precision            :: Eold = 0.0d0
  double precision            :: dE = 0.0d0
  double precision            :: alpha = 2.0d10
  double precision            :: Econv = 1.0d-5

!  E0 = Etot

  doFTSCF = .true. 
  doGKSCF = .false.
  ResetOcc = .true.

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Calculating FTD-MP2 energy"
    write(*,*) "    "
    write(*,*) &
"                T_el        E(SCF)         E(MBPT2)        E(sing)        Etot"
  endif
  if(.true.)then
    do i=1,22
      Tel = Tlist(i)
      call run_scf(.false.)
      call trans_full(.false.)
      call calc_Embpt2()
      if(doprint)then
        write(*,'("    FTD-MP2: ",5(" ",F12.5," "))') Tel, Etot, Embpt2,E_1,Etot+Embpt2
      endif
    enddo
  else
    write(*,*)  "                T_el        E(SCF)         E(MBPT2)        E(sing)        Etot        Etot(-Singles)"

    Tel = 20000.0d0
    
    do i=1,1000
      call run_scf(.false.)
      call trans_full(.false.)
      call calc_Embpt2()
      Ezero = Etot+Embpt2-E_1
      dE = Ezero-Eold
      write(*,'("    FTD-MP2: ",6(" ",F14.5," "))') Tel, Etot, Embpt2,E_1,Etot+Embpt2,Ezero
      if(abs(Ezero-Eold)<Econv)then
        write(*,*)  '    ...converged'
        exit
      endif
      Tel = Tel + Tstep
      call run_scf(.false.)
      call trans_full(.false.)
      call calc_Embpt2()
      Eplus = Etot+Embpt2-E_1

      Tel = Tel - 2.0d0*Tstep
      call run_scf(.false.)
      call trans_full(.false.)
      call calc_Embpt2()
      Eminus = Etot+Embpt2-E_1
      Tel = Tel + Tstep

      dEdT = (Eplus-Eminus)/(2.0d0*Tstep)
      !write(*,*) dEdT
      Tel = Tel - alpha*dEdT
      Eold = Ezero
 
    enddo
  endif

end subroutine calc_ftdmp2

end module module_ftdoo


