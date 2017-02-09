module module_ints
  implicit none

contains

!#############################################
!#            Calculate Integrals
!#############################################
subroutine calc_Ints()
  use module_data,      only : mol_name,basis_set,dim_1e
  implicit none
  character(len=1024)       :: command 
  character(len=1024)       :: dir

  write(*,*) ""
  write(*,*) "  Running Psi4 integral package ..."
  write(*,*) ""

  call get_environment_variable("SAMSA_DIR",dir)
 
  write(command,*) "perl ",trim(dir),"/tools/runpsiints.pl " , &
                           trim(mol_name)," "                , &
                           trim(basis_set)," "               , &  
                           dim_1e

  write(*,*) "    ",trim(command)
  write(*,*) ""

  call system(trim(command))

  call calc_Overlap()
  call calc_Kinetic()
  call calc_CoreEl()
  call calc_Hcore()
  call calc_ERI()

end subroutine calc_Ints

!#############################################
!#            Calculate Overlap
!#############################################
subroutine calc_Overlap()
  use module_data,        only : Sij,dim_1e
  implicit none
  integer                     :: i,j,a,b
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch  

  write(*,*) "     ... fetching Sij "

  open(20,file=trim(filename),status='old',action='read')

  do
    read(20,*) scratch
    if(trim(scratch) == "OVERLAP") exit
  enddo

  do i=1,dim_1e
    do j=1,i
      read(20,*) a, b, temp
      Sij(i,j) = temp
      Sij(j,i) = temp
    enddo
  enddo

  close(20)

end subroutine


!#############################################
!#            Calculate Kinetic
!#############################################
subroutine calc_Kinetic()
  use module_data,        only : Tij,dim_1e
  implicit none
  integer                     :: i,j,a,b
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch

  write(*,*) "     ... fetching Tij "

  open(30,file=trim(filename),status='old',action='read')

  do
    read(30,*) scratch
    if(trim(scratch) == "KINETIC") exit
  enddo

  do i=1,dim_1e
    do j=1,i
      read(30,*) a,b,temp
!      write(*,*) temp
      Tij(i,j) = temp
      Tij(j,i) = temp
    enddo
  enddo

  close(30)

end subroutine


!#############################################
!#            Calculate CoreEl
!#############################################
subroutine calc_CoreEl()
  use module_data,        only : Vij,dim_1e
  implicit none
  integer                     :: i,j,a,b
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch

  write(*,*) "     ... fetching Vij "

  open(40,file=trim(filename),status='old',action='read')

  do
    read(40,*) scratch
    if(trim(scratch) == "POTENTIAL") exit
  enddo

  do i=1,dim_1e
    do j=1,i
      read(40,*) a,b,temp
!      write(*,*) temp
      Vij(i,j) = temp
      Vij(j,i) = temp
    enddo
  enddo

  close(40)

end subroutine


!#############################################
!#            Calculate Hcore
!#############################################
subroutine calc_Hcore()
  use module_data,        only : Hcore,Tij,Vij
  implicit none

  write(*,*) "     ... calc'ing Hcore "

  Hcore = Tij + Vij

  deallocate(Tij,Vij)

end subroutine


!#############################################
!#            Calculate Repulsion
!#############################################
subroutine calc_ERI()
  use module_data,        only : ERI,dim_1e
  implicit none
  integer             :: i,j,k,l,a,b,c,d
  integer             :: ij,kl,ijkl
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch

  write(*,*) "     ... fetching ERI "

  open(50,file=trim(filename),status='old',action='read')

  do
    read(50,*) scratch
    if(trim(scratch) == "REPULSION") exit
  enddo

  do i=0,(dim_1e-1)
    do j=0,i
      do k=0,(dim_1e-1)
        do l=0,k
          ij = i*(i+1)/2 + j
          kl = k*(k+1)/2 + l
          if(ij>=kl)then
            ijkl = ij*(ij+1)/2 + kl + 1
            read(50,*) a,b,c,d,temp
            ERI(ijkl) = temp
          endif
        enddo
      enddo
    enddo
  enddo

  close(50)

end subroutine

!#############################################
!#              Calc 2e Indices                 
!#############################################
subroutine Index2e(i,j,ij)
  implicit none
  integer, intent(in)           :: i,j
  integer, intent(out)          :: ij

  if (i>j)then
    ij = i*(i+1)/2 + j
  else
    ij = j*(j+1)/2 + i
  endif

end subroutine Index2e

end module module_ints
