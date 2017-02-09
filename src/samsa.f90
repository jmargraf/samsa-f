program samsa
  use module_data,   only : mol_name,exe_path
  use module_data,   only : read_input
  use module_data,   only : dimensions
  use module_data,   only : allocate_SCFmat 
  use module_energy, only : calc_Enuc,calc_Emp2
  use module_ints,   only : calc_Ints
  use module_scf,    only : do_guess,dia_S,run_SCF
  use module_props,  only : print_Eigen
  use module_trans,  only : trans_full
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

! run post-scf
  call trans_full()
  call calc_Emp2()

end program samsa


