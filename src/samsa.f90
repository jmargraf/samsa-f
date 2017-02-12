program samsa
  use module_data,   only : mol_name,exe_path
  use module_data,   only : read_input
  use module_data,   only : dimensions
  use module_data,   only : allocate_SCFmat 
  use module_data,   only : doMP2,Fract,doEVPT2
  use module_energy, only : calc_Enuc,calc_Emp2
  use module_ints,   only : calc_Ints,get_Occ
  use module_scf,    only : do_guess,dia_S,run_SCF
  use module_props,  only : print_Eigen
  use module_trans,  only : trans_full,transdone
  use module_occupy, only : calc_EVPT2
  implicit none

! read input and calculate dimensions
  call get_command_argument(0, exe_path)
  call get_command_argument(1, mol_name)
  write(*,*) ""
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(*,*) "  SAMSA - a Semiempirical and Ab Initio Molecular Simulation Application  "
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(*,*) ""
  write(*,*) "  (c) 2017 Johannes Margraf "
  write(*,*) ""
  write(*,*) "  starting calculation: ", trim(mol_name)
  write(*,*) "  ", trim(exe_path)
  write(*,*) ""

  call read_input()
  call dimensions()
  call allocate_SCFmat()
  call calc_Enuc()

! using fractional occupation numbers?
  if(Fract==1)then
    call get_Occ()
  endif

! calc ints
  call calc_Ints()

! calc transformation matrix
  call dia_S() 
  
! initial guess
  call do_guess()

! run scf
  call run_SCF()

! run properties
  call print_Eigen()

! run occupation number schemes
  if(doEVPT2)then
    if(.not.transdone)then
      call trans_full()
    endif
    call calc_EVPT2()
  endif

! using fractional occupation numbers in post-HF?
  if(Fract==2)then
    call get_Occ()
  endif

! run post-scf
  if(doMP2)then
    if(.not.transdone)then
      call trans_full()
    endif
    call run_MP2()
    call run_DCPT2()
  endif

end program samsa



subroutine run_MP2()
  use module_energy, only : calc_Emp2
  use module_energy, only : Emp2     
  use module_energy, only : Emp2frac 
  use module_energy, only : E_OS,E_SS,E_SSx,E_SSc
  use module_energy, only : E_AAx,E_AAc,E_BBx,E_BBc
  use module_data,   only : Spins
  implicit none

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-MBPT(2) correlation energy"
    write(*,*) "    "
    call calc_Emp2()
    write(*,*) "    "
    write(*,'("      E_OS   = ",F12.6)') E_OS
    write(*,'("      E_SS   = ",F12.6)') E_SS
    write(*,'("      E_SSx  = ",F12.6)') E_SSx
    write(*,'("      E_SSc  = ",F12.6)') E_SSc
    write(*,*) "    "
    write(*,'("      EMP2   = ",F12.6)') Emp2
    write(*,*) "    "
  elseif(spins==2)then
    write(*,*) "    "
    write(*,*) "    Calculating UHF-MBPT(2) correlation energy"
    write(*,*) "    "
    call calc_Emp2()
    write(*,*) "    "
    write(*,'("      E_OS   = ",F12.6)') E_OS
    write(*,'("      E_AAx  = ",F12.6)') E_AAx
    write(*,'("      E_AAc  = ",F12.6)') E_AAc
    write(*,'("      E_BBx  = ",F12.6)') E_BBx
    write(*,'("      E_BBc  = ",F12.6)') E_BBc
    write(*,*) "    "
    write(*,'("      EMP2   = ",F12.6)') Emp2
    write(*,*) "    "
    write(*,*) "    "
    write(*,*) "    Calculating UHF-MBPT(2) correlation energy (nOcc version)"
    write(*,*) "    "
    write(*,*) "    "
    write(*,'("      E_OS   = ",F12.6)') E_OS
    write(*,'("      E_AAx  = ",F12.6)') E_AAx
    write(*,'("      E_AAc  = ",F12.6)') E_AAc
    write(*,'("      E_BBx  = ",F12.6)') E_BBx
    write(*,'("      E_BBc  = ",F12.6)') E_BBc
    write(*,*) "    "
    write(*,'("      EMP2   = ",F12.6)') Emp2frac
    write(*,*) "    "
  endif

end subroutine run_MP2

subroutine run_DCPT2()
  use module_energy, only : calc_Edcpt2
  use module_energy, only : Edcpt2
  use module_energy, only : Edcpt2fr
  use module_energy, only : E_OS,E_SS,E_SSx,E_SSc
  use module_energy, only : E_AAx,E_AAc,E_BBx,E_BBc
  use module_data,   only : Spins
  implicit none

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-DCPT(2) correlation energy"
    write(*,*) "    "
    call calc_Edcpt2()
    write(*,*) "    "
    write(*,'("      E_OS  = ",F12.6)') E_OS
    write(*,'("      E_SS  = ",F12.6)') E_SS
    write(*,'("      E_SSx = ",F12.6)') E_SSx
    write(*,'("      E_SSc = ",F12.6)') E_SSc
    write(*,*) "    "
    write(*,'("      EDCPT2 = ",F12.6)') Edcpt2
    write(*,*) "    "
  elseif(spins==2)then
    write(*,*) "    "
    write(*,*) "    Calculating UHF-DCPT(2) correlation energy"
    write(*,*) "    "
    call calc_Edcpt2()
    write(*,*) "    "
    write(*,'("      E_OS  = ",F12.6)') E_OS
    write(*,'("      E_AAx = ",F12.6)') E_AAx
    write(*,'("      E_AAc = ",F12.6)') E_AAc
    write(*,'("      E_BBx = ",F12.6)') E_BBx
    write(*,'("      E_BBc = ",F12.6)') E_BBc
    write(*,*) "    "
    write(*,'("      EDCPT2 = ",F12.6)') Edcpt2

!    write(*,*) "    "
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-MBPT(2) correlation energy (nOcc version)"
!    write(*,*) "    "
!    write(*,*) "    "
!    write(*,*) "    "
!    write(*,*) "    "
!    write(*,'("      E_OS  = ",F12.6)') E_OS
!    write(*,'("      E_AAx = ",F12.6)') E_AAx
!    write(*,'("      E_AAc = ",F12.6)') E_AAc
!    write(*,'("      E_BBx = ",F12.6)') E_BBx
!    write(*,'("      E_BBc = ",F12.6)') E_BBc
!    write(*,*) "    "
!    write(*,'("      EMP2  = ",F12.6)') Emp2frac
!    write(*,*) "    "
  endif

end subroutine run_DCPT2







