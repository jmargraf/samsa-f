module module_scf
  implicit none

! Energies           
!  double precision                            :: Enuc     = 0.0d0
!  double precision                            :: Eelec    = 0.0d0         
!  double precision                            :: Etot     = 0.0d0          

contains

!#############################################
!#              Initial Guess       
!#############################################
subroutine guess()
  use module_data, only      : Fock,Hcore,Spins,dim_1e
  use module_io, only        : print_Mat
  implicit none
  integer                   :: i

  write(*,*) ""
  write(*,*) "    Initial guess..."
  write(*,*) ""

  do i=1,Spins
    Fock(:,:,i) = Hcore(:,:)
    ! add random perturbation
  enddo

  call print_Mat(Hcore(:,:),dim_1e,"Hcore")

  call print_Mat(Fock(:,:,1),dim_1e,"Fock")

  call dia_Fock()

end subroutine guess

!#############################################
!#          Diagonalize Fockian     
!#############################################
subroutine dia_fock()
  use module_data, only      : Fock,Coef,Eps,dim_1e,Spins
  implicit none
  integer                         :: i,j,lwork,inf
  double precision, allocatable   :: work(:)

  lwork = dim_1e*(3+dim_1e/2)
  allocate(work(lwork))

  write(*,*) ""
  write(*,*) "    Diagonalizing..."
  write(*,*) ""

  do i=1,Spins
    Coef(:,:,i) = Fock(:,:,i)
    call dsyev('V','U',dim_1e,Coef(1:dim_1e,1:dim_1e,i),dim_1e,Eps(1:dim_1e,i),work,lwork,inf)
  enddo

  do i=1,dim_1e
    do j=1,Spins
      write(*,*) Eps(i,j)
    enddo
  enddo

  deallocate(work)

end subroutine dia_fock

!#############################################
!#          Diagonalize Overlap     
!#############################################
subroutine dia_S()
  use module_data, only      : S12,Sij,dim_1e
  implicit none
  integer                         :: i,j,lwork,inf
  double precision, allocatable   :: work(:)
  double precision, allocatable   :: S_vec(:,:)
  double precision, allocatable   :: S_val(:)
  double precision, allocatable   :: S_mat(:,:)
  double precision, allocatable   :: temp(:,:)

  lwork = dim_1e*(3+dim_1e/2)
  allocate(work(lwork),           &
           S_vec(dim_1e,dim_1e),  &
           S_mat(dim_1e,dim_1e),  &
           S_val(dim_1e),         &
           temp(dim_1e,dim_1e))

  write(*,*) ""
  write(*,*) "    Diagonalizing S..."
  write(*,*) ""

  S_vec = Sij
  call dsyev('V','U',dim_1e,S_vec(1:dim_1e,1:dim_1e),dim_1e,S_val(1:dim_1e),work,lwork,inf)

  do i=1,dim_1e
    write(*,*) S_val(i)
  enddo

  deallocate(work)

  S_val = S_val**(-0.5d0)
! stretch matrix
  do i=1,dim_1e
    do j=1,dim_1e
      if(i == j)then
        S_mat(i,j) = S_val(i)
      else
        S_mat(i,j) = 0.0d0
      endif
    enddo
  enddo

! S12 = S_vec x S_mat x transpose(S_vec)
! temp = S_mat x transpose(S_vec) 
  call dgemm('N','T',dim_1e,dim_1e,dim_1e,1.0d0,S_mat,dim_1e,S_vec,dim_1e,0.0d0,temp,dim_1e)
! S12 = S_vec x temp 
  call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,S_vec,dim_1e,temp,dim_1e,0.0d0,S12,dim_1e)

end subroutine dia_S


end module module_scf

