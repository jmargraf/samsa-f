module module_data
  implicit none

! Calculation details
  character(len=1024)                           :: exe_path
  character(len=32)                             :: mol_name
  integer                                       :: natoms
  integer                                       :: dim_1e,dim_2e
  integer                                       :: nel, noccA, noccB
  character(len=2), allocatable                 :: atom(:)
  double precision, allocatable                 :: xyz(:,:)
  double precision, allocatable                 :: rAB(:,:)

! Basis set details
  integer, allocatable                          :: Basis(:)
!  character(len=3), allocatable                 :: bastype(:)
  integer, allocatable                          :: Bastype(:)
  character(len=3)                              :: basis_set = "min"

! SCF options
  integer                                       :: charge = 0
  integer                                       :: mult = 1
  character(len=3)                              :: calctype = "rhf"
  integer                                       :: spins = 1
  integer                                       :: maxSCF = 1000
  double precision                              :: Econv = 1.0d-9
  double precision                              :: Dconv = 1.0d-9
  double precision                              :: Damp = 0.5d0
  logical                                       :: DoDamp = .true.
  logical                                       :: DoMP2 = .false.
  integer                                       :: Fract = 0 
  character(len=4)                              :: guess = "core"
  double precision                              :: Par(4) = 0.0d0
  integer                                       :: scaletype = 1

! SCF matrices
  double precision, allocatable                 :: Fock(:,:,:)  
  double precision, allocatable                 :: Coef(:,:,:)  
  double precision, allocatable                 :: Dens(:,:,:)  
  double precision, allocatable                 :: Eps(:,:)           
  double precision, allocatable                 :: Hcore(:,:)       
  double precision, allocatable                 :: Occ(:,:)

! Integrals
  double precision, allocatable                 :: Sij(:,:)         
  double precision, allocatable                 :: S12(:,:)
  double precision, allocatable                 :: Vij(:,:)         
  double precision, allocatable                 :: Tij(:,:)         
  double precision, allocatable                 :: ERI(:)      
! MO integrals           
  double precision, allocatable                 :: MOI(:)
  double precision, allocatable                 :: SMO(:,:,:,:)
           
contains

!#############################################
!#                 Read Input
!#############################################
subroutine read_input()
  use module_constants, only : ang2bohr
  implicit none
  character(len=36)        :: filename 
  integer                  :: i,j
  character(len=1024)      :: keywords
  character(len=64)        :: argument

  filename = trim(mol_name) // ".xyz"

  open(10,file=trim(filename),status='old',action='read')
  read(10,*) natoms
  call allocate_geo

  write(*,*) ""
  write(*,*) "  Keywords:"
  write(*,*) ""

  read(10,'(A)') keywords
  do
    argument = split_string(keywords)
    if(argument=='')exit
    call parse_option(argument)  
  enddo

  if(calctype == "RHF")then
    spins = 1
  elseif(calctype == "UHF")then
    spins = 2
  endif

  do i=1,natoms
    read(10,*) atom(i),xyz(i,1),xyz(i,2),xyz(i,3)
  enddo
  xyz = xyz*ang2bohr

  close(10)

  write(*,*) ""
  write(*,*) "  Number of atoms = ", natoms
  write(*,*) "  Charge          = ", charge
  write(*,*) "  Mutiplicity     = ", mult
  write(*,*) ""

  write(*,*) ""
  write(*,*) "  Geometry (Bohr):"
  write(*,*) ""
  do i=1,natoms
    write(*,*) "    ",atom(i),xyz(i,1),xyz(i,2),xyz(i,3)
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
           S12(dim_1e,dim_1e),         &
           Sij(dim_1e,dim_1e),         &
           Vij(dim_1e,dim_1e),         &
           Tij(dim_1e,dim_1e),         &
           ERI(dim_2e)                 &
          )

end subroutine allocate_SCFmat


!#############################################
!#       Get Basis Dimensions and Info
!#############################################
subroutine dimensions()
  use module_constants, only : coreq
  use module_io,        only : print_SVec,print_Vec
  implicit none
  integer                  :: i,j,k,l,iSpins
  integer                  :: ij,kl,ijkl
  logical                  :: even = .true.
  integer                  :: Z,iBas

! Dimension of 1e Matrices
  dim_1e = 0

  iBas = 0

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
  elseif(basis_set == "svp")then
    do i=1,natoms
      if((atom(i) == "H"))then
        dim_1e = dim_1e + 2
      elseif((atom(i) == "He"))then
        dim_1e = dim_1e + 5
      elseif((atom(i) == "Li") .or. &
             (atom(i) == "Be"))then
        dim_1e = dim_1e + 9
      elseif((atom(i) == "B")  .or. &
             (atom(i) == "C")  .or. &
             (atom(i) == "N")  .or. &
             (atom(i) == "O")  .or. &
             (atom(i) == "F")  .or. &
             (atom(i) == "Ne"))then
        dim_1e = dim_1e + 14
      elseif((atom(i) == "Na"))then  
        dim_1e = dim_1e + 15
      elseif((atom(i) == "Mg") .or. &
             (atom(i) == "Al") .or. &
             (atom(i) == "Si") .or. &
             (atom(i) == "P")  .or. &
             (atom(i) == "S")  .or. &
             (atom(i) == "Cl") .or. &
             (atom(i) == "Ar"))then
        dim_1e = dim_1e + 18
      endif
    enddo    
  endif 

  allocate(bastype(dim_1e),basis(dim_1e))
  call basis_types()

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

  write(*,*) ""
  write(*,*) "  1e Dimension      = ", dim_1e
  write(*,*) "  2e Dimension      = ", dim_2e

! Number of electrons
  nel = 0
  noccA = 0
  noccB = 0
  do i=1,natoms
    call  coreq(atom(i),Z)
    nel = nel + Z
  enddo
  nel = nel - charge

  if(spins == 1)then
    if(mod(nel,2) == 0)then
      even = .true.
      noccA = nel/2 + (mult-1)/2
      noccB = nel/2 - (mult-1)/2
    else
      write(*,*) "  !Error! RHF calc impossible with odd no of electrons "
      call exit(666)
    endif
  elseif(spins == 2)then
    if(mod(nel,2) == 0)then
      even = .true.
      noccA = nel/2 + (mult-1)/2
      noccB = nel/2 - (mult-1)/2
      if(mod(mult,2) == 0)then
        write(*,*) "  !Error! Even no of electrons and mult ", mult," is impossible"
        call exit(666)
      endif
    else
      even = .false.
      noccA = nel/2 + 1 + (mult-1)/2
      noccB = nel/2 - (mult-1)/2
      if(mod(mult,2) == 1)then
        write(*,*) "  !Error! Odd no of electrons and mult ", mult," is impossible"
        call exit(666)
      endif
    endif
  endif

  allocate(Occ(dim_1e,spins))
  Occ = 0.0d0

  do iSpins=1,spins
    if(iSpins==1)then
      do i=1,nOccA
        Occ(i,iSpins) = 2.0d0 * 1.0d0/dble(spins)
      enddo
    endif
    if(iSpins==2)then
      do i=1,nOccB
        Occ(i,iSpins) = 2.0d0 * 1.0d0/dble(spins)
      enddo
    endif
  enddo

  if(Spins==1)then
    call print_Vec(Occ(:,1),dim_1e,18,"Occupation Numbers")
  elseif(Spins==2)then
    call print_SVec(Occ(:,:),dim_1e,18,"Occupation Numbers")
  endif

  write(*,*) "  No of electrons   = ", nel 
  write(*,*) "  No of occ orbs(A) = ", noccA 
  write(*,*) "  No of occ orbs(B) = ", noccB
  write(*,*) ""

end subroutine dimensions

!#############################################
!#                Basis Types    
!#############################################
subroutine basis_types()
  implicit none
  integer                  :: i,iBas

  iBas = 0
  if(basis_set == "min")then
    do i=1,natoms
      if((atom(i) == "H") .or. (atom(i) == "He"))then
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
      elseif((atom(i) == "Li") .or. &
             (atom(i) == "Be") .or. &
             (atom(i) == "B")  .or. &
             (atom(i) == "C")  .or. &
             (atom(i) == "N")  .or. &
             (atom(i) == "O")  .or. &
             (atom(i) == "F")  .or. &
             (atom(i) == "Ne"))then
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 2
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
      elseif((atom(i) == "Na") .or. &
             (atom(i) == "Mg") .or. &
             (atom(i) == "Al") .or. &
             (atom(i) == "Si") .or. &
             (atom(i) == "P")  .or. &
             (atom(i) == "S")  .or. &
             (atom(i) == "Cl") .or. &
             (atom(i) == "Ar"))then
!       Nothing
      endif
    enddo
  elseif(basis_set == "svp")then
    do i=1,natoms
      if((atom(i) == "H"))then
!       2 1S      
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
      elseif((atom(i) == "He"))then
!       2 1S 1 2P
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 2
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
      elseif((atom(i) == "Li") .or. &
             (atom(i) == "Be"))then
!       1 1S 2 2S 2 2P
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 2
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 2
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
      elseif((atom(i) == "B")  .or. &
             (atom(i) == "C")  .or. &
             (atom(i) == "N")  .or. &
             (atom(i) == "O")  .or. &
             (atom(i) == "F")  .or. &
             (atom(i) == "Ne"))then
!       1 1S 2 2S 2 2P 1 3D
        iBas = iBas + 1
        Bastype(iBas) = 1
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 2
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 2
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 3
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 4
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 4
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 4
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 4
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
        iBas = iBas + 1
        Bastype(iBas) = 4
        Basis(iBas) = i
        write(*,*) "Basisfunction ", iBas, " on atom ", i, " is type ", Bastype(iBas)
      elseif((atom(i) == "Na"))then
!       Nothing
      elseif((atom(i) == "Mg") .or. &
             (atom(i) == "Al") .or. &
             (atom(i) == "Si") .or. &
             (atom(i) == "P")  .or. &
             (atom(i) == "S")  .or. &
             (atom(i) == "Cl") .or. &
             (atom(i) == "Ar"))then
!       Nothing
      endif
    enddo
  endif

end subroutine basis_types



!#############################################
!#              Parse Options
!#############################################
subroutine parse_option(argument)
  implicit none
  character(len=64),intent(in)     :: argument
  character(len=64)                :: uargument

  uargument = uppercase(argument)

  if(uargument(1:7)=='CHARGE=')then
    read(UArgument(8:),*) charge
    write(*,*) '    CHARGE    = ', charge

  elseif(uargument(1:5)=='MULT=')then
    read(UArgument(6:),*) mult
    write(*,*) '    MULT      = ', mult

  elseif(uargument(1:5)=='CALC=')then
    read(UArgument(6:),*) calctype
    write(*,*) '    CALC      = ', calctype

  elseif(uargument(1:7)=='MAXSCF=')then
    read(UArgument(8:),*) maxSCF
    write(*,*) '    MAXSCF    = ', maxSCF

  elseif(uargument(1:6)=='ECONV=')then
    read(UArgument(7:),*) Econv
    write(*,*) '    ECONV     = ', Econv

  elseif(uargument(1:6)=='DCONV=')then
    read(UArgument(7:),*) Dconv
    write(*,*) '    DCONV     = ', Dconv

  elseif(uargument(1:5)=='DAMP=')then
    read(UArgument(6:),*) DAMP
    DoDamp = .true.
    write(*,*) '    DAMP      = ',DAMP

  elseif(uargument(1:7)=='DODAMP=')then
    if(uargument(8:9)=='ON')then
      DoDamp = .true.
      write(*,*) '    Damping turned on'
    elseif(uargument(8:10)=='OFF')then
      DoDamp = .false.
      write(*,*) '    Damping turned off'
    else
      write(*,*) "    Unknown Argument :", uargument
    endif

  elseif(uargument(1:6)=='GUESS=')then
    if(uargument(7:10)=='CORE')then
      guess = "core"
      write(*,*) '    GUESS     = CORE'
    elseif(uargument(7:12)=='HUCKEL')then
      guess = "huck"
      write(*,*) '    GUESS     = HUCK'
    else
      write(*,*) "    Unknown Argument :", uargument
    endif

  elseif(uargument(1:6)=='BASIS=')then
    if(uargument(7:10)=='MIN')then
      basis_set = "min"
      write(*,*) '    BASIS     = MIN'
    elseif(uargument(7:10)=='SVP')then
      basis_set = "svp"
      write(*,*) '    BASIS     = SVP'
    else
      write(*,*) "    Unknown Argument :", uargument
    endif

  elseif(uargument(1:5)=='PAR1=')then
    read(UArgument(6:),*) Par(1)
    write(*,*) '    PAR1      = ',Par(1)

  elseif(uargument(1:5)=='PAR2=')then
    read(UArgument(6:),*) Par(2)
    write(*,*) '    PAR2      = ',Par(2)

  elseif(uargument(1:5)=='PAR3=')then
    read(UArgument(6:),*) Par(3)
    write(*,*) '    PAR3      = ',Par(3)

  elseif(uargument(1:5)=='PAR4=')then
    read(UArgument(6:),*) Par(4)
    write(*,*) '    PAR4      = ',Par(4)

  elseif(uargument(1:5)=='SCAL=')then
    read(UArgument(6:),*) scaletype
    write(*,*) '    SCAL      = ',scaletype

  elseif(uargument(1:5)=='MBPT2')then
    doMP2 = .true.
    write(*,*) '    MBPT(2) '

  elseif(uargument(1:6)=='FRACT=')then
    read(UArgument(7:),*) Fract
    write(*,*) '    FRACT     = ',Fract

  else
    write(*,*) "    Unknown Argument :", uargument
  endif
  
end subroutine parse_option


!#############################################
!#              Split Spring              
!#############################################
function split_string(string) result(word)
  implicit none
  character(len=*), intent(inout) :: string
  character(len=len(string)) :: word
  integer :: i

! Remove trailing blanks
  string=adjustl(string)
! Split first word from string
  i = index(string,' ')
  if(i<1)then
    word   = string
    string = ''
  else
    word   = string(1:i-1)
    string = adjustl(string(i+1:))
  endif

end function split_string


!#############################################
!#              Make Uppercase                   
!#############################################
function uppercase(string) result(upper)
  implicit none
  character(len=*), intent(in) :: string
  character(len=len(string))   :: upper
  integer :: i

! Convert to uppercase
  do i=1,len(string)
    if(string(i:i) >= "a" .and. string(i:i) <= "z") then
      upper(i:i)=achar(iachar(string(i:i))-32)
    else
      upper(i:i)=string(i:i)
    endif
  enddo
end function uppercase


end module module_data



