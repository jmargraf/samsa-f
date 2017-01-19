module module_io
  implicit none

! Energies           

contains

!#############################################
!#               Print Matrix       
!#############################################
subroutine print_Mat(Mat,n,MatName)
  implicit none
  double precision, intent(in)  :: Mat(:,:)
  integer, intent(in)           :: n
  character(len=32), intent(in) :: MatName
  integer                       :: i,j

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


end module module_io

