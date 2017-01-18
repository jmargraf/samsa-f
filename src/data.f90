module module_data
  use module_constants, only : ang2bohr
  implicit none

! Calcularion details
  character(len=32)                             :: mol_name
  integer                                       :: natoms
  integer                                       :: charge
  integer                                       :: mult
  integer                                       :: dim_1e,dim_2e
  integer                                       :: nel, noccA, noccB
  character(len=2), allocatable, dimension(:)   :: atom
  double precision, allocatable, dimension(:,:) :: xyz
  double precision, allocatable, dimension(:,:) :: rAB
  logical                                       :: uhf
  character(len=3)                              :: calctype = "rhf"

! Basis set details
  integer, allocatable, dimension(:)            :: basis
  character(len=3), allocatable, dimension(:)   :: bastyp
  character(len=3)                              :: basis_set = "min"

! SCF options
  integer                                       :: maxSCF = 1000
  double precision                              :: Econv = 1.0d-9
  double precision                              :: Dconv = 1.0d-9
  double precision                              :: Damp = 0.5
  logical                                       :: DoDamp = .true.
  character                                     :: guess = "core"

contains

!#############################################
!#                 Read Input
!#############################################
subroutine read_input
  implicit none
  character(len=36)        :: filename 
  integer                  :: i

  filename = trim(mol_name) // ".xyz"

  open(10,file=trim(filename),status='old',action='read')
  read(10,*) natoms
  call allocate_geo

  read(10,*) charge, mult

  do i=1,natoms
    read(10,*) atom(i),xyz(i,1),xyz(i,2),xyz(i,3)
  enddo
  xyz = xyz*ang2bohr

  write(*,*) ""
  write(*,*) "  Number of atoms = ", natoms
  write(*,*) "  Charge          = ", charge
  write(*,*) "  Mutiplicity     = ", mult
  write(*,*) ""

  do i=1,natoms
    write(*,*) atom(i),xyz(i,1),xyz(i,2),xyz(i,3)
  enddo
  write(*,*) ""

end subroutine read_input 


!#############################################
!#          Allocate Geometry Arrays
!#############################################
subroutine allocate_geo
  implicit none

  allocate(atom(natoms),xyz(natoms,3))

end subroutine allocate_geo


end module module_data



