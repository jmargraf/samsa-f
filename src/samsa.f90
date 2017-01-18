program samsa
  use module_data,   only : mol_name
  use module_data,   only : read_input
  use module_data,   only : dimensions
  use module_data,   only : allocate_SCFmat 
  use module_energy, only : calc_Enuc
  implicit none

! read input and calculate dimensions
  call get_command_argument(1, mol_name)
  write(*,*) "Start calculation: ", trim(mol_name)," !"

  call read_input()
  call dimensions()
  call allocate_SCFmat()
  call calc_Enuc()

! calc ints

! initial guess

! run scf

! run properties

! run post-scf


end program samsa


