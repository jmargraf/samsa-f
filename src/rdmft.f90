module module_rdmft
  implicit none
  double precision, allocatable               :: Jmat(:,:,:)  ! Coulomb operator matrix
  double precision, allocatable               :: Kmat(:,:,:)  ! Exchange operator matrix
  double precision                            :: Erdmft = 0.0d0

contains

!#############################################
!#                 Run 1-RDMFT
!#############################################
subroutine run_rdmft()
  use module_data, only      : dim_1e,Hcore,NOEps,NOCoef,Coef,Fock
  use module_wavefun, only   : Fock_to_MO,Fock_to_AO
  implicit none
  integer                   :: i,j
  double precision          :: E1e,E2ec,E2ex,E2e

  allocate(Jmat(1:dim_1e,1:dim_1e,1),Kmat(1:dim_1e,1:dim_1e,1))

  Jmat = 0.0d0
  Kmat = 0.0d0

  call AOMat_to_MO()

  Coef = NOCoef
  call Fock_to_AO()

!  Coef = NOCoef
  call Fock_to_MO()

  Erdmft = 0.0d0
  E1e  = 0.0d0
  E2e  = 0.0d0
  E2ec = 0.0d0
  E2ex = 0.0d0

!  do iSpin=1,Spins
!    do i=1,dim_1e
!      do j=1,dim_1e
!        E1e = E1e + Dens(i,j,iSpin)*Hcore(i,j)*(2.0d0/dble(Spins))
!        E2e = E2e + Dens(i,j,iSpin)*(Fock(i,j,iSpin)-Hcore(i,j))*(1.0d0/dble(Spins))
!      enddo
!    enddo
!  enddo

  do i=1,dim_1e
      E1e = E1e + 1.0d0*NOEps(i,1)*Hcore(i,i)
    do j=1,dim_1e
!      E2ec = E2ec + 1.0d0*NOEps(i,1)*NOEps(j,1)*Jmat(i,j,1)
!      E2ex = E2ex - 0.5d0*NOEps(i,1)*NOEps(j,1)*Kmat(i,j,1)
      E2e = E2e + NOEps(i,1)*NOEps(j,1)*(Fock(i,j,1)-Hcore(i,j))      
    enddo
  enddo

!  E2e = E2ex + E2ec
  Erdmft = E1e + E2e!+ E2ec + E2ex

  write(*,*) "         E1e = ", E1e   
  write(*,*) "         E2ex = ", E2ex   
  write(*,*) "         E2ec = ", E2ec
  write(*,*) "         E2e  = ", E2e
  write(*,*) "      Erdmft = ", Erdmft

  deallocate(Jmat,Kmat)

end subroutine run_rdmft


!#############################################
!#     Transform J,K and Hcore to MO Basis
!#############################################
subroutine AOMat_to_MO()
  use module_data, only      : NOCoef,dim_1e,Spins,Hcore,Sm12,Coef
  use module_io, only        : print_Mat
  use module_wavefun, only   : calc_Dens
  implicit none
  integer                         :: i,j
  double precision, allocatable   :: temp(:,:)

!  write(*,*) ""
!  write(*,*) "    transform Fockian to MO-basis ..."
!  write(*,*) ""

  allocate(temp(dim_1e,dim_1e))

  temp = 0.0d0

! deorthogonalize NO-det
! do i=1,Spins
!  $::Coeff = $::Sm12 x $::Coeff;
    temp(1:dim_1e,1:dim_1e) = NOCoef(1:dim_1e,1:dim_1e,1)
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,        &
                Sm12(1:dim_1e,1:dim_1e),dim_1e,           &
                temp(1:dim_1e,1:dim_1e),dim_1e,           &
                0.0d0,NOCoef(1:dim_1e,1:dim_1e,1),dim_1e)
!    call print_Mat(Coef(:,:,i),dim_1e,8,"deort WF")
! enddo

!  do i=1,dim_1e
!    do j=1,dim_1e
!      Coef(dim_1e+1-i,j,1) = NOCoef(i,j,1)
!    enddo
!  enddo

!  call calc_Dens()

  call build_JK()

  temp = 0.0d0

! Coulomb
! do i=1,Spins
!  $::Jmat = transpose($::Coef) x $::Jmat x $::Coef;
!   call print_Mat(Jmat(:,:,i),dim_1e,7,"AO Jmat")
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                Jmat(1:dim_1e,1:dim_1e,1),dim_1e,     &
                NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,   &
                0.0d0,temp(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,   &
                temp(1:dim_1e,1:dim_1e),dim_1e,       &
                0.0d0,Jmat(1:dim_1e,1:dim_1e,1),dim_1e)
!   call print_Mat(Jmat(:,:,i),dim_1e,7,"MO Jmat")
! enddo

  temp = 0.0d0

! Exchange
! do i=1,Spins
!  $::Kmat = transpose($::Coef) x $::Kmat x $::Coef;
!   call print_Mat(Kmat(:,:,i),dim_1e,7,"AO Kmat")
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                Kmat(1:dim_1e,1:dim_1e,1),dim_1e,     &
                NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,   &
                0.0d0,temp(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,   &
                temp(1:dim_1e,1:dim_1e),dim_1e,       &
                0.0d0,Kmat(1:dim_1e,1:dim_1e,1),dim_1e)
!   call print_Mat(Kmat(:,:,i),dim_1e,7,"MO Kmat")
! enddo

  temp = 0.0d0

! Hcore
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                Hcore(1:dim_1e,1:dim_1e),dim_1e,      &
                NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,   &
                0.0d0,temp(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,   &
                temp(1:dim_1e,1:dim_1e),dim_1e,       &
                0.0d0,Hcore(1:dim_1e,1:dim_1e),dim_1e)
!   call print_Mat(Hcore(:,:),dim_1e,8,"MO Hcore")
  deallocate(temp)

end subroutine AOMat_to_MO


!#############################################
!#           Build J and K in AO basis
!#############################################
subroutine build_JK()
  use module_data, only      : dim_1e,ERI,Dens,S12
  implicit none
  integer                       :: i,j,k,l
  integer                       :: ij,kl
  integer                       :: ik,jl
  integer                       :: ijkl,ikjl
!  double precision, allocatable :: temp(:,:),OrtDens(:,:)

!  allocate(temp(1:dim_1e,1:dim_1e),OrtDens(1:dim_1e,1:dim_1e))

! Dens = S12 x D x S12
!  call print_Mat(S12(:,:),dim_1e,3,"S12")
!  call print_Mat(temp(:,:),dim_1e,4,"temp")
! D x S12
!  call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,S12(:,:),dim_1e,Dens(:,:,1),dim_1e,0.0d0,temp(:,:),dim_1e)
!  call print_Mat(NOCoef(:,:,1),dim_1e,7,"NO Coef")
! S12 x temp
!  call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,temp,dim_1e,S12,dim_1e,0.0d0,OrtDens(:,:),dim_1e)
!  call print_Mat(NOCoef(:,:,1),dim_1e,7,"NO Coef")

!  deallocate(temp)

  Jmat = 0.0d0
  Kmat = 0.0d0

  do i=0,dim_1e-1
    do j=0,dim_1e-1
      do k=0,dim_1e-1
        do l=0,dim_1e-1
          ij = i*(i+1)/2 + j
          kl = k*(k+1)/2 + l
          ik = i*(i+1)/2 + k
          jl = j*(j+1)/2 + l
          if(j>i) ij = j*(j+1)/2 + i
          if(l>k) kl = l*(l+1)/2 + k
          if(k>i) ik = k*(k+1)/2 + i
          if(l>j) jl = l*(l+1)/2 + j
          ijkl = ij*(ij+1)/2 + kl +1
          ikjl = ik*(ik+1)/2 + jl +1
          if(kl>ij) ijkl = kl*(kl+1)/2 + ij +1
          if(jl>ik) ikjl = jl*(jl+1)/2 + ik +1

          Jmat(i+1,j+1,1)       = Jmat(i+1,j+1,1)                                 &
                                + Dens(k+1,l+1,1)*ERI(ijkl)
!                                + OrtDens(k+1,l+1)*ERI(ijkl)

          Kmat(i+1,j+1,1)       = Kmat(i+1,j+1,1)                                 &
                                + Dens(k+1,l+1,1)*ERI(ikjl)
!                                + OrtDens(k+1,l+1)*ERI(ikjl)  

        enddo
      enddo
    enddo
  enddo

!  deallocate(OrtDens)

end subroutine build_JK

end module module_rdmft
