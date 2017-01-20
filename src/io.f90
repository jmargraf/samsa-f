module module_io
  implicit none

! Energies           

contains

!#############################################
!#               Print Matrix       
!#############################################
subroutine print_Mat(Mat,n,nname,MatName)
  implicit none
  double precision, intent(in)     :: Mat(:,:)
  integer, intent(in)              :: n
  integer, intent(in)              :: nname
  character(len=nname), intent(in) :: MatName
  integer                          :: i,j

  write(*,*) ""
  write(*,*) "  Matrix: ", trim(MatName)
  write(*,*) ""

  do i=1,n
    do j=1,n
      write(*,'(" ",F10.5," ")',advance="no") Mat(i,j)
!      write(*,*)  Mat(i,j)
    enddo
    write(*,*) ""
  enddo

end subroutine print_Mat


!#############################################
!#                Print Vector       
!#############################################
subroutine print_Vec(Vec,n,nname,VecName)
  use module_data,            only : noccA
  implicit none
  double precision, intent(in)     :: Vec(:)
  integer, intent(in)              :: n
  integer, intent(in)              :: nname
  character(len=nname), intent(in) :: VecName
  integer                          :: i

  write(*,*) ""
  write(*,*) "  Vector: ", trim(VecName)
  write(*,*) ""

  do i=1,n
    write(*,'(" ",F10.5," ")',advance="no") Vec(i)
    if(i==noccA) write(*,'(A)',advance="no") " HOMO"
    write(*,*) ""
  enddo

end subroutine print_Vec


!#############################################
!#             Print Spin-Vectors       
!#############################################
subroutine print_SVec(Vec,n,nname,VecName)
  use module_data,            only : noccA,noccB
  implicit none
  double precision, intent(in)     :: Vec(:,:)
  integer, intent(in)              :: n
  integer, intent(in)              :: nname
  character(len=nname), intent(in) :: VecName
  integer                          :: i,j

  write(*,*) ""
  write(*,*) "  Vector: ", trim(VecName)
  write(*,*) ""

  do i=1,n
    do j=1,2
      write(*,'(" ",F10.5," ")',advance="no") Vec(i,j)
    enddo
    if(i==noccA) write(*,'(A)',advance="no") " SOMOa"
    if(i==noccB) write(*,'(A)',advance="no") " SOMOb"
    write(*,*) ""
  enddo

end subroutine print_SVec


end module module_io

