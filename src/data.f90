module module_data
  implicit none

! Calcularion details
  character(len=32)                             :: mol_name
  integer                                       :: natoms
  integer                                       :: charge
  integer                                       :: mult
  integer                                       :: dim_1e,dim_2e
  integer                                       :: nel, noccA, noccB
  character(len=2), allocatable                 :: atom(:)
  double precision, allocatable                 :: xyz(:,:)
  double precision, allocatable                 :: rAB(:,:)
  logical                                       :: uhf
  character(len=3)                              :: calctype = "rhf"
  integer                                       :: spins = 1

! Basis set details
  integer, allocatable                          :: basis(:)
  character(len=3), allocatable                 :: bastyp(:)
  character(len=3)                              :: basis_set = "min"

! SCF options
  integer                                       :: maxSCF = 1000
  double precision                              :: Econv = 1.0d-9
  double precision                              :: Dconv = 1.0d-9
  double precision                              :: Damp = 0.5
  logical                                       :: DoDamp = .true.
  character                                     :: guess = "core"

! SCF matrices
  double precision, allocatable                 :: Fock(:,:,:)  
  double precision, allocatable                 :: Coef(:,:,:)  
  double precision, allocatable                 :: Dens(:,:,:)  
  double precision, allocatable                 :: Eps(:,:)           
  double precision, allocatable                 :: Hcore(:,:)       
  double precision, allocatable                 :: Sij(:,:)         
  double precision, allocatable                 :: Vij(:,:)         
  double precision, allocatable                 :: Tij(:,:)         
  double precision, allocatable                 :: ERI(:)      
           
contains

!#############################################
!#                 Read Input
!#############################################
subroutine read_input()
  use module_constants, only : ang2bohr
  implicit none
  character(len=36)        :: filename 
  integer                  :: i,j

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

  write(*,*) ""
  write(*,*) "  Geometry (Bohr):"
  write(*,*) ""
  do i=1,natoms
    write(*,*) " ",atom(i),xyz(i,1),xyz(i,2),xyz(i,3)
  enddo
  write(*,*) ""

end subroutine read_input 


!#############################################
!#          Allocate Geometry Arrays
!#############################################
subroutine allocate_geo()
  implicit none

  allocate(atom(natoms),xyz(natoms,3),rAB(natoms,natoms))

end subroutine allocate_geo

!#############################################
!#          Allocate Matrices
!#############################################
subroutine allocate_SCFmat()
  implicit none

  write(*,*) ""
  write(*,*) "  allocating matrices..."
  write(*,*) ""

  allocate(Fock(dim_1e,dim_1e,spins),  &
           Coef(dim_1e,dim_1e,spins),  &
           Dens(dim_1e,dim_1e,spins),  &
           Eps(dim_1e,spins),          &
           Hcore(dim_1e,dim_1e),       &
           Sij(dim_1e,dim_1e),         &
           Vij(dim_1e,dim_1e),         &
           Tij(dim_1e,dim_1e),         &
           ERI(dim_2e),                &
          )

end subroutine allocate_SCFmat


!#############################################
!#       Get Basis Dimensions and Info
!#############################################
subroutine dimensions()
  use module_constants, only : coreq
  implicit none
  integer                  :: i,j,k,l
  integer                  :: ij,kl,ijkl
  logical                  :: even = .true.
  integer                  :: Z

! Dimension of 1e Matrices
  dim_1e = 0

  if(basis_set == "min")then
    do i=1,natoms
      if((atom(i) == "H") .or. (atom(i) == "He"))then
        dim_1e = dim_1e + 1
      elseif((atom(i) == "Li") .or. &
             (atom(i) == "Be") .or. &
             (atom(i) == "B")  .or. &
             (atom(i) == "C")  .or. &
             (atom(i) == "N")  .or. &
             (atom(i) == "O")  .or. &
             (atom(i) == "F")  .or. &
             (atom(i) == "Ne"))then
        dim_1e = dim_1e + 5
      elseif((atom(i) == "Na") .or. &
             (atom(i) == "Mg") .or. &
             (atom(i) == "Al") .or. &
             (atom(i) == "Si") .or. &
             (atom(i) == "P")  .or. &
             (atom(i) == "S")  .or. &
             (atom(i) == "Cl") .or. &
             (atom(i) == "Ar"))then
        dim_1e = dim_1e + 9
      endif
    enddo
  endif 

! Dimension of 2e Matrices
  dim_2e = 0
  do i=0,(dim_1e-1)
    do j=0,i
      do k=0,(dim_1e-1)
        do l=0,k
          ij = i*(i+1)/2 + j
          kl = k*(k+1)/2 + l
          if(ij>=kl)then
            dim_2e = dim_2e + 1
            ijkl = ij*(ij+1)/2 + kl
          endif
        enddo
      enddo
    enddo
  enddo

! Number of electrons
  nel = 0
  noccA = 0
  noccB = 0
  do i=1,natoms
    call  coreq(atom(i),Z)
    nel = nel + Z
  enddo
  nel = nel - charge

  write(*,*) ""
  write(*,*) "  1e Dimension    = ", dim_1e
  write(*,*) "  2e Dimension    = ", dim_2e
  write(*,*) "  No of electrons = ", nel 
  write(*,*) "  No of electrons = ", nel 
  write(*,*) ""


end subroutine dimensions


end module module_data



