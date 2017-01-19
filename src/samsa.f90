program samsa
  use module_data,   only : mol_name,exe_path
  use module_data,   only : read_input
  use module_data,   only : dimensions
  use module_data,   only : allocate_SCFmat 
  use module_energy, only : calc_Enuc
  use module_ints,   only : calc_Ints
  use module_scf,    only : guess,dia_S
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
  call dia_S() 
 
! initial guess
  call guess()

! run scf

! run properties

! run post-scf


end program samsa


