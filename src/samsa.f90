program samsa
  use module_data, only : mol_name,read_input 
  implicit none

! read input
  call get_command_argument(1, mol_name)
  write(*,*) "Start calculation: ", trim(mol_name)," !"

  call read_input()

! calc ints

! initial guess

! run scf

! run properties

! run post-scf


end program samsa


