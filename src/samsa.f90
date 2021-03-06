program samsa
  use module_data,    only : mol_name,exe_path
  use module_data,    only : read_input
  use module_data,    only : dimensions
  use module_data,    only : allocate_SCFmat 
  use module_data,    only : doMBPT2,doDCPT2,Fract,doGKT
  use module_data,    only : doLCCD,doDSCF,doDFrac,DoRDMFT,doCCD,doCIS,doFTDMP2
  use module_data,    only : spin_homo,spin_lumo,basis_set,doCCDeriv
  use module_energy,  only : calc_Enuc,calc_Embpt2
  use module_ints,    only : calc_Ints,get_Occ
  use module_scf,     only : run_SCF
  use module_wavefun, only : do_guess,dia_S,calc_natorbs
  use module_props,   only : print_Eigen,pop_Mulliken,pop_Loewdin
  use module_trans,   only : trans_full,transdone,trans_ucc
  use module_occupy,  only : calc_GKT,run_corrF
  use module_cc,      only : calc_Eccd
  use module_cc_chen, only : calc_Eccd_chen
  use module_dscf,    only : calc_dSCF,calc_dfrac,calc_fracCC,calc_dEcc
  use module_rdmft,   only : run_rdmft
  use module_cis,     only : calc_CIS
  use module_ftdoo,   only : calc_ftdmp2
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
  call run_SCF(.true.)

! run properties
  call print_Eigen()
  if(basis_set /= 'tzd')then
    call pop_Mulliken()
    call pop_Loewdin()
  endif

! calc NOs?
  call calc_natorbs()
  if(DoRDMFT)then
    call run_rdmft()
  endif

! run occupation number schemes
  if(doGKT)then
    if(.not.transdone)then
      call trans_full(.true.)
    endif
    call calc_GKT(.true.)
!  elseif(doDSCF)then
!    call calc_dSCF(.true.)
   endif

! correct Fock matrix
!  call run_corrF(.true.)

! using fractional occupation numbers in post-HF?
  if(Fract==2)then
    call get_Occ()
  endif

! run post-scf
  if(doMBPT2)then
    if(.not.transdone)then
      call trans_full(.true.)
    endif
    call run_MBPT2()
  endif

  if(doDSCF)then
    call calc_dSCF(.true.)
  endif

  if(doDFrac)then
    !call calc_dfrac(.true.)
    call calc_fracCC(.true.,spin_homo,spin_lumo)
  endif

  if(doLCCD)then
    if(.not.transdone)then
      call trans_full(.true.)
    endif
    call trans_ucc(.true.)
    if(Fract==0)then
      call calc_Eccd(.false.,.false.)
    elseif(Fract>0)then
      call calc_Eccd(.false.,.true.)
      call calc_Eccd_chen(.false.,.true.)
    endif
  endif

  if(doCCD)then
    if(.not.transdone)then
      call trans_full(.true.)
    endif
    call trans_ucc(.true.)
    if(Fract==0)then
      call calc_Eccd(.true.,.false.)
    elseif(Fract>0)then
      call calc_Eccd(.true.,.true.)
      call calc_Eccd_chen(.true.,.true.)
    endif
  endif

  if(doCCDeriv)then
     call calc_dEcc(.true.)
  endif

  if(.false.)then
    call run_LMPT2()
  endif

  if(doCIS)then
    if(.not.transdone)then
      call trans_full(.true.)
    endif
    call trans_ucc(.true.)
    call calc_CIS()
  endif

  if(doDCPT2)then
    if(.not.transdone)then
      call trans_full(.true.)
    endif
    call run_DCPT2()
  endif

  if(doFTDMP2)then
    call calc_ftdmp2(.true.)
  endif

  write(*,*) "    "
  write(*,*) "    Calculation finished successfully! "
  write(*,*) "    "
  write(*,*) "    "

end program samsa

subroutine run_LMPT2()
  use module_energy, only : calc_Elmpt2
  use module_energy, only : Elmpt2
  use module_energy, only : E_OS,E_SS,E_SSx,E_SSc,Etot
!  use module_energy, only : E_AAx,E_AAc,E_BBx,E_BBc
  use module_data,   only : Spins
  implicit none

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-LMPT(2) correlation energy"
    write(*,*) "    "
    call calc_Elmpt2()
    write(*,*) "    "
    write(*,'("      E_OS   = ",F12.6)') E_OS
    write(*,'("      E_SS   = ",F12.6)') E_SS
    write(*,'("      E_SSx  = ",F12.6)') E_SSx
    write(*,'("      E_SSc  = ",F12.6)') E_SSc
    write(*,*) "    "
    write(*,'("      ELMPT2 = ",F12.6)') Elmpt2
    write(*,*) "    "
    write(*,'("      Etotal = ",F12.6)') Elmpt2+Etot
  endif

end subroutine run_LMPT2

subroutine run_MBPT2()
  use module_energy, only : calc_Embpt2
  use module_energy, only : Embpt2,Etot     
  use module_energy, only : Embpt2f 
  use module_energy, only : E_OS,E_SS,E_SSx,E_SSc,E_1,E_1f
  use module_energy, only : E_AAx,E_AAc,E_BBx,E_BBc
  use module_energy, only : E_AAxf,E_AAcf,E_BBxf,E_BBcf
  use module_energy, only : E_OSf
  use module_data,   only : Spins
  implicit none

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-MBPT(2) correlation energy"
    write(*,*) "    "
    call calc_Embpt2()
    write(*,*) "    "
    write(*,'("      E_1    = ",F12.6)') E_1
    write(*,'("      E_OS   = ",F12.6)') E_OS
    write(*,'("      E_SS   = ",F12.6)') E_SS
    write(*,'("      E_SSx  = ",F12.6)') E_SSx
    write(*,'("      E_SSc  = ",F12.6)') E_SSc
    write(*,*) "    "
    write(*,'("      EMBPT2 = ",F12.6)') Embpt2
    write(*,'("      Etotal = ",F12.6)') Embpt2+Etot
    write(*,*) "    "
  elseif(spins==2)then
    write(*,*) "    "
    write(*,*) "    Calculating UHF-MBPT(2) correlation energy"
    write(*,*) "    "
    call calc_Embpt2()
    write(*,*) "    "
    write(*,'("      E_1    = ",F12.6)') E_1
    write(*,'("      E_OS   = ",F12.6)') E_OS
    write(*,'("      E_AAx  = ",F12.6)') E_AAx
    write(*,'("      E_AAc  = ",F12.6)') E_AAc
    write(*,'("      E_BBx  = ",F12.6)') E_BBx
    write(*,'("      E_BBc  = ",F12.6)') E_BBc
    write(*,*) "    "
    write(*,'("      EMBPT2 = ",F12.6)') Embpt2
    write(*,*) "    "
    write(*,*) "    "
    write(*,*) "    Calculating UHF-MBPT(2) correlation energy (nOcc version)"
    write(*,*) "    "
    write(*,*) "    "
    write(*,'("      E_1    = ",F12.6)') E_1f
    write(*,'("      E_OS   = ",F12.6)') E_OSf
    write(*,'("      E_AAx  = ",F12.6)') E_AAxf
    write(*,'("      E_AAc  = ",F12.6)') E_AAcf
    write(*,'("      E_BBx  = ",F12.6)') E_BBxf
    write(*,'("      E_BBc  = ",F12.6)') E_BBcf
    write(*,*) "    "
    write(*,'("      EMBPT2 = ",F12.6)') Embpt2f
    write(*,*) "    "
  endif

end subroutine run_MBPT2

subroutine run_DCPT2()
  use module_energy, only : calc_Edcpt2
  use module_energy, only : Edcpt2
  use module_energy, only : Edcpt2f
  use module_energy, only : E_OS,E_SS,E_SSx,E_SSc
  use module_energy, only : E_AAx,E_AAc,E_BBx,E_BBc
  use module_energy, only : E_AAxf,E_AAcf,E_BBxf,E_BBcf
  use module_energy, only : E_OSf
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

    write(*,*) "    "
    write(*,*) "    "
    write(*,*) "    Calculating UHF-DCPT(2) correlation energy (nOcc version)"
    write(*,*) "    "
    write(*,*) "    "
    write(*,'("      E_OS  = ",F12.6)') E_OSf
    write(*,'("      E_AAx = ",F12.6)') E_AAxf
    write(*,'("      E_AAc = ",F12.6)') E_AAcf
    write(*,'("      E_BBx = ",F12.6)') E_BBxf
    write(*,'("      E_BBc = ",F12.6)') E_BBcf
    write(*,*) "    "
    write(*,'("      EDCPT2 = ",F12.6)') Edcpt2f
    write(*,*) "    "
  endif

end subroutine run_DCPT2







