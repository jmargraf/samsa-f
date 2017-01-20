module module_scf
  implicit none

! Convergence        
  double precision                            :: dE       = 0.0d0
  double precision                            :: Drmsd    = 0.0d0         
  logical                                     :: SCFconv  = .false. 

contains

!#############################################
!#              Run SCF Loop        
!#############################################
subroutine run_SCF
  use module_data, only     : MaxSCF, Damp
  use module_energy, only   : Eold,calc_Energy,Etot
  implicit none
  integer                  :: iSCF

  write(*,*) ""
  write(*,*) "   Starting SCF ..."
  write(*,*) ""
  write(*,*) "        #           Etot                 dE                  Drmsd      Damp"

  do iSCF=1,MaxSCF

    Eold = Etot

    call ort_Fock()

    call dia_Fock()

    call deort_Fock()

    call calc_Dens()

    call calc_Energy()

    dE = Etot - Eold     

    !  iSCF   Etot  dE   Drmsd
    write(*,'("SCF  ",I5,"  ",3(F18.10,"  "),F5.2)') iSCF, Etot, dE, Drmsd, Damp

    call check_Conv()

    if(SCFconv) exit

    call calc_Fock()

  enddo

end subroutine run_SCF

!#############################################
!#              Initial Guess       
!#############################################
subroutine guess()
  use module_data, only      : Fock,Hcore,Spins,dim_1e
  use module_io, only        : print_Mat
  implicit none
  integer                   :: i

!  write(*,*) ""
!  write(*,*) "    Initial guess..."
!  write(*,*) ""

  do i=1,Spins
    Fock(:,:,i) = Hcore(:,:)
    ! add random perturbation
  enddo

!  call print_Mat(Hcore(:,:),dim_1e,5,"Hcore")
!  call print_Mat(Fock(:,:,1),dim_1e,4,"Fock")

end subroutine guess


!#############################################
!#           Calculate Fock Matrix 
!#############################################
subroutine calc_Fock()
  use module_data, only      : Dens,Fock,Spins,dim_1e,ERI,Hcore
  use module_data, only      : DoDamp, Damp
  use module_io, only        : print_Mat
  implicit none
  integer                         :: iSpin,i,j,k,l
  integer                         :: ij,kl,ijkl
  integer                         :: ik,jl,ikjl
  double precision, allocatable   :: Fockold(:,:,:)

!  write(*,*) ""
!  write(*,*) "    Build new Fock..."
!  write(*,*) ""

  allocate(Fockold(dim_1e,dim_1e,Spins))

  Fockold = Fock

  do iSpin=1,Spins
    Fock(1:dim_1e,1:dim_1e,iSpin) = Hcore(1:dim_1e,1:dim_1e)
  enddo

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
          !UHF
          if(Spins==2)then
            do iSpin=1,Spins
              Fock(i+1,j+1,iSpin) = Fock(i+1,j+1,iSpin)        &
                                  + Dens(k+1,l+1,1)*ERI(ijkl)  &
                                  + Dens(k+1,l+1,2)*ERI(ijkl)  &
                                  - Dens(k+1,l+1,iSpin)*ERI(ikjl)                                                         
            enddo
          elseif(Spins==1)then
            Fock(i+1,j+1,1)       = Fock(i+1,j+1,1)        &
                                  + Dens(k+1,l+1,1)*(2.0d0*ERI(ijkl) - ERI(ikjl))
          endif
        enddo
      enddo

    enddo
  enddo

  ! Damping    
  if(DoDamp)then
    Fock = Damp*Fock + (1.0d0-Damp)*Fockold
  endif

  do iSpin=1,Spins
!    call print_Mat(Fock(:,:,iSpin),dim_1e,4,"Fock")
  enddo

  deallocate(Fockold)

end subroutine calc_Fock


!#############################################
!#         Calculate Density Matrix 
!#############################################
subroutine calc_Dens()
  use module_data, only      : Dens,Coef,Spins,dim_1e
  use module_data, only      : noccA,noccB
  use module_io, only        : print_Mat
  implicit none
  integer                         :: iSpin,i,j,k,l
  double precision, allocatable   :: Densold(:,:,:)

  allocate(Densold(dim_1e,dim_1e,Spins))
  Drmsd = 0.0d0
  Densold = Dens

!  write(*,*) ""
!  write(*,*) "    Calc density matrix"
!  write(*,*) ""

  do i=1,dim_1e
    do j=1,i
      Dens(i,j,1) = 0.0d0
      do k=1,nOccA
        Dens(i,j,1) = Dens(i,j,1) + Coef(i,k,1)*Coef(j,k,1)
      enddo
      Dens(j,i,1) = Dens(i,j,1)
      Drmsd = Drmsd + (Dens(i,j,1)-Densold(i,j,1))**2
    enddo
  enddo
!  call print_Mat(Dens(:,:,1),dim_1e,5,"DensA")

  
  if(Spins==2)then
    do i=1,dim_1e
      do j=1,i
        Dens(i,j,2) = 0.0d0
        do k=1,nOccB
          Dens(i,j,2) = Dens(i,j,2) + Coef(i,k,2)*Coef(j,k,2)
        enddo
        Dens(j,i,2) = Dens(i,j,2)
        Drmsd = Drmsd + (Dens(i,j,2)-Densold(i,j,2))**2
      enddo
    enddo
!    call print_Mat(Dens(:,:,2),dim_1e,5,"DensB")
  endif

  Drmsd = sqrt(Drmsd/dble(dim_1e*dim_1e*Spins))


end subroutine calc_Dens



!#############################################
!#          Orthogonalize Fock      
!#############################################
subroutine ort_fock()
  use module_data, only      : Fock,Coef,S12,dim_1e,Spins
  use module_io, only        : print_Mat
  implicit none
  integer                         :: i
  double precision, allocatable   :: temp(:,:)

!  write(*,*) ""
!  write(*,*) "    orthogonalizing Fock ..."
!  write(*,*) ""

  allocate(temp(dim_1e,dim_1e))

  do i=1,Spins
!  $::Coeff = transpose($::S12) x $::Fock x $::S12;
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                Fock(1:dim_1e,1:dim_1e,i),dim_1e,   &
                S12(1:dim_1e,1:dim_1e),dim_1e,      &
                0.0d0,temp(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                S12(1:dim_1e,1:dim_1e),dim_1e,      &
                temp(1:dim_1e,1:dim_1e),dim_1e,     &
                0.0d0,Coef(1:dim_1e,1:dim_1e,i),dim_1e)
!    call print_Mat(Coef(:,:,i),dim_1e,8,"ort Fock")
  enddo


  deallocate(temp)

end subroutine ort_fock


!#############################################
!#            Deorthogonalize Fock      
!#############################################
subroutine deort_fock()
  use module_data, only      : Coef,S12,dim_1e,Spins
  use module_io, only        : print_Mat
  implicit none
  integer                         :: i
  double precision, allocatable   :: temp(:,:)

!  write(*,*) ""
!  write(*,*) "    deorthogonalizing Fock ..."
!  write(*,*) ""

  allocate(temp(dim_1e,dim_1e))

  do i=1,Spins
!  $::Coeff = $::S12 x $::Coeff;
    temp(1:dim_1e,1:dim_1e) = coef(1:dim_1e,1:dim_1e,i) 
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                S12(1:dim_1e,1:dim_1e),dim_1e,      &
                temp(1:dim_1e,1:dim_1e),dim_1e,   &
                0.0d0,Coef(1:dim_1e,1:dim_1e,i),dim_1e)
!    call print_Mat(Coef(:,:,i),dim_1e,8,"deort WF")
  enddo

  deallocate(temp)

end subroutine deort_fock



!#############################################
!#          Diagonalize Fockian     
!#############################################
subroutine dia_fock()
  use module_data, only      : Fock,Coef,Eps,dim_1e,Spins
  use module_io, only        : print_Mat,print_Vec
  implicit none
  integer                         :: i,j,lwork,inf
  double precision, allocatable   :: work(:)

  lwork = dim_1e*(3+dim_1e/2)
  allocate(work(lwork))

!  write(*,*) ""
!  write(*,*) "    Diagonalizing..."
!  write(*,*) ""

  do i=1,Spins
    call dsyev('V','U',dim_1e,Coef(1:dim_1e,1:dim_1e,i),dim_1e,Eps(1:dim_1e,i),work,lwork,inf)
!    call print_Mat(Coef(:,:,i),dim_1e,6,"ort WF")
  enddo

  do i=1,Spins
!    call print_Vec(Eps(:,i),dim_1e,11,"Eigenvalues")
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

!  write(*,*) ""
!  write(*,*) "    Diagonalizing S..."
!  write(*,*) ""

  S_vec = Sij
  call dsyev('V','U',dim_1e,S_vec(1:dim_1e,1:dim_1e),dim_1e,S_val(1:dim_1e),work,lwork,inf)

!  call print_Vec(S_val(:),dim_1e,11,"Ovrlpvalues")

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

  deallocate(temp,S_vec,S_mat,S_val)

end subroutine dia_S

!#############################################
!#             Check Convergence       
!#############################################
subroutine check_Conv()
  use module_data, only      : Econv 
  use module_data, only      : Dconv 
  use module_data, only      : DoDamp, Damp
  implicit none

  if(abs(dE) < 1.0d-2 .and. Damp /= 1.0d0  & 
                      .and. Damp /= 0.7d0  &
                      .and. DoDamp)then
    Damp = 0.7d0
    write(*,*) "    Damping set to 0.7"
  endif

  if(abs(dE) < 1.0d-5 .and. Damp /= 1.0d0  & 
                      .and. DoDamp)then
    Damp = 1.0d0
    DoDamp =.false.
    write(*,*) "    Damping turned off"
  endif

  if(Drmsd<Dconv .and. abs(dE)<Econv)then
    write(*,*) " "
    write(*,*) "    SCF converged (hurra)!"
    write(*,*) " "

    SCFconv = .true.
  endif

end subroutine check_Conv

end module module_scf

