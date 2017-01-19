module module_ints
  implicit none

! Energies           
! double precision                            :: Enuc     = 0.0d0
! double precision                            :: Eelec    = 0.0d0         
! double precision                            :: Etot     = 0.0d0          

contains

!#############################################
!#            Calculate Integrals
!#############################################
subroutine calc_Ints()
  use module_data, only      : mol_name,basis_set,dim_1e
  implicit none
  character(len=1024)       :: command 
! integer                   :: i,j
! integer                   :: q1,q2

  write(*,*) ""
  write(*,*) "  Running Psi4 integral package ..."
  write(*,*) ""

!  command = "perl ../tools/runpsiints.pl " // &
!                    trim(mol_name) // " " //  &
!                   trim(basis_set) // " " // 

  write(command,*) "perl ../tools/runpsiints.pl " , &
                               trim(mol_name)," " , &
                               trim(basis_set)," ", &  
                               dim_1e

  write(*,*) "      ",trim(command)

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
  integer                     :: i,j
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
      read(20,*) temp
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
  integer                     :: i,j
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch

  write(*,*) "     ... fetching Tij "

  open(20,file=trim(filename),status='old',action='read')

  do
    read(20,*) scratch
    if(trim(scratch) == "KINETIC") exit
  enddo

  do i=1,dim_1e
    do j=1,i
      read(20,*) temp
      Tij(i,j) = temp
      Tij(j,i) = temp
    enddo
  enddo

  close(20)

end subroutine


!#############################################
!#            Calculate CoreEl
!#############################################
subroutine calc_CoreEl()
  use module_data,        only : Vij,dim_1e
  implicit none
  integer                     :: i,j
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch

  write(*,*) "     ... fetching Vij "

  open(20,file=trim(filename),status='old',action='read')

  do
    read(20,*) scratch
    if(trim(scratch) == "POTENTIAL") exit
  enddo

  do i=1,dim_1e
    do j=1,i
      read(20,*) temp
      Vij(i,j) = temp
      Vij(j,i) = temp
    enddo
  enddo

  close(20)

end subroutine


!#############################################
!#            Calculate Hcore
!#############################################
subroutine calc_Hcore()
  use module_data,        only : Hcore,Tij,Vij
  implicit none

  write(*,*) "     ... calc'ing Hcore "

  Hcore = Tij + Vij

end subroutine


!#############################################
!#            Calculate Repulsion
!#############################################
subroutine calc_ERI()
  use module_data,        only : ERI,dim_1e
  implicit none
  integer             :: i,j,k,l
  integer             :: ij,kl,ijkl
  character(len=4)            :: filename="ints"
  double precision            :: temp
  character(len=32)           :: scratch

  write(*,*) "     ... fetching ERI "

  open(20,file=trim(filename),status='old',action='read')

  do
    read(20,*) scratch
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
            read(20,*) temp
            ERI(ijkl) = temp
          endif
        enddo
      enddo
    enddo
  enddo

  close(20)

end subroutine


!#############################################
!#            Test LAPACK        
!#############################################
subroutine test_lapack()
  implicit none

  double precision :: a(4,4), d(4), e(3), t(3), ew(4)
  integer          :: inf, tw(8), z

!  external dgeprt
!  external dsteqr, dsytrd

  data a/1, 2, 3, 4, &
        2, 2, 6, 4,  & 
        3, 6, 5, 6,  &
        4, 4, 6, 6/

!  call dgeprt(4,5,a,'a=')
  call dsytrd('U', 4, a, 4, d, e, t, tw, 4, inf)

  if (inf .eq. 0) then
     write (*,*) 'successful tridiagonal reduction'
  else if (inf .lt. 0) then
     write (*,*) 'illegal value at: %d', -inf
     stop
  else
     write (*,*) 'unknown result (can''t happen!)'
     stop
  end if

end subroutine


end module module_ints

