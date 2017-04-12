module module_wavefun
  implicit none
  double precision                            :: Drmsd    = 0.0d0

! Convergence        
contains

!#############################################
!#        Semicanonical Transformation
!#############################################
subroutine SemiCanonical()
  use module_data, only      : Fock,dim_1e,Spins,nOccA,nOccB,Eps
  use module_io, only        : print_Mat
  implicit none
  integer                         :: iSpin,nOcc,nVirt,lWork,i,j,inf
  double precision, allocatable   :: temp(:,:)
  double precision, allocatable   :: Vec(:)
  double precision, allocatable   :: work(:)

  write(*,*) ""
  write(*,*) "    Semicanonical Transformation ..."
  write(*,*) ""

! Occupied-Occupied Block
  do iSpin=1,Spins
    if(iSpin==1)then
      nOcc=nOccA
    else
      nOcc=nOccB
    endif
    lwork = nOcc*(3+nOcc/2)
    allocate(work(lwork))
    allocate(temp(1:nOcc,1:nOcc),Vec(1:nOcc))

    temp(1:nOcc,1:nOcc) = Fock(1:nOcc,1:nOcc,iSpin)

    call dsyev('V','U',nOcc,temp(1:nOcc,1:nOcc),nOcc,Vec(1:nOcc),work,lwork,inf)

    temp = 0.0d0
    do i=1,nOcc
      temp(i,i) = Vec(i)
      Eps(i,iSpin) = Vec(i)
    enddo
    Fock(1:nOcc,1:nOcc,iSpin) = temp(1:nOcc,1:nOcc)

    deallocate(temp,vec,work)
  enddo

! Virtual-Virtual Block
  do iSpin=1,Spins
    if(iSpin==1)then
      nOcc  = nOccA
      nVirt = dim_1e-nOcc
    else
      nOcc=nOccB
      nVirt = dim_1e-nOcc
    endif
    lwork = nVirt*(3+nVirt/2)
    allocate(work(lwork))
    allocate(temp(1:nVirt,1:nVirt),Vec(1:nVirt))

    temp(1:nVirt,1:nVirt) = Fock(nOcc+1:dim_1e,nOcc+1:dim_1e,iSpin)

    call dsyev('V','U',nVirt,temp(1:nVirt,1:nVirt),nVirt,Vec(1:nVirt),work,lwork,inf)

    temp = 0.0d0
    do i=1,nVirt
      temp(i,i) = Vec(i)
      Eps(nOcc+i,iSpin) = Vec(i)
    enddo
    Fock(nOcc+1:dim_1e,nOcc+1:dim_1e,iSpin) = temp(1:nVirt,1:nVirt)

    deallocate(temp,vec,work)
  enddo


end subroutine SemiCanonical


!#############################################
!#       Transform Fockian to MO Basis
!#############################################
subroutine Fock_to_MO()
  use module_data, only      : Fock,Coef,dim_1e,Spins
  use module_io, only        : print_Mat
  implicit none
  integer                         :: i
  double precision, allocatable   :: temp(:,:)

!  write(*,*) ""
!  write(*,*) "    transform Fockian to MO-basis ..."
!  write(*,*) ""

  allocate(temp(dim_1e,dim_1e))

  do i=1,Spins
!  $::Fock = transpose($::Coef) x $::Fock x $::Coef;
!    call print_Mat(Fock(:,:,i),dim_1e,7,"AO Fock")
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                Fock(1:dim_1e,1:dim_1e,i),dim_1e,   &
                Coef(1:dim_1e,1:dim_1e,i),dim_1e,   &
                0.0d0,temp(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                Coef(1:dim_1e,1:dim_1e,i),dim_1e,   &
                temp(1:dim_1e,1:dim_1e),dim_1e,     &
                0.0d0,Fock(1:dim_1e,1:dim_1e,i),dim_1e)
!    call print_Mat(Fock(:,:,i),dim_1e,7,"MO Fock")
  enddo

  deallocate(temp)

end subroutine Fock_to_MO


!#############################################
!#       Transform Fockian to AO Basis
!#############################################
subroutine Fock_to_AO()
  use module_data, only      : Fock,Coef,Sij,dim_1e,Spins
  use module_io, only        : print_Mat
  implicit none
  integer                         :: i
  double precision, allocatable   :: temp1(:,:)
  double precision, allocatable   :: temp2(:,:)
  double precision, allocatable   :: temp3(:,:)

!  write(*,*) ""
!  write(*,*) "    transform Fockian to AO-basis ..."
!  write(*,*) ""

  allocate(temp1(dim_1e,dim_1e))
  allocate(temp2(dim_1e,dim_1e))
  allocate(temp3(dim_1e,dim_1e))

  do i=1,Spins
!  $::Fock = $::Sij x $::Coeff x $::Fock x transpose($::Coeff) x $::Sij;
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                Sij(1:dim_1e,1:dim_1e),dim_1e,        &
                Coef(1:dim_1e,1:dim_1e,i),dim_1e,     &
                0.0d0,temp1(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                Coef(1:dim_1e,1:dim_1e,i),dim_1e,     &
                Sij(1:dim_1e,1:dim_1e),dim_1e,        &
                0.0d0,temp2(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                Fock(1:dim_1e,1:dim_1e,i),dim_1e,     &
                temp2(1:dim_1e,1:dim_1e),dim_1e,      &
                0.0d0,temp3(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,    &
                temp1(1:dim_1e,1:dim_1e),dim_1e,      &
                temp3(1:dim_1e,1:dim_1e),dim_1e,      &
                0.0d0,Fock(1:dim_1e,1:dim_1e,i),dim_1e)
!    call print_Mat(Fock(:,:,i),dim_1e,7,"AO Fock")
  enddo
  
  deallocate(temp1,temp2,temp3)

end subroutine Fock_to_AO



!#############################################
!#              Initial Guess       
!#############################################
subroutine do_guess()
  use module_data,      only : Fock,Hcore,Spins,dim_1e
  use module_data,      only : guess
  use module_io,        only : print_Mat
  implicit none
  integer                   :: i

!  write(*,*) ""
!  write(*,*) "    Initial guess..."
!  write(*,*) ""

!  write(*,*) "guess = ", guess

  if(guess == "core")then
    write(*,*) ""
    write(*,*) "    core initial guess..."
    write(*,*) ""
    do i=1,Spins
      Fock(:,:,i) = Hcore(:,:)
      ! add random perturbation
    enddo
  elseif(Guess == "huck")then
    call guess_huckel()
  endif

!  call print_Mat(Hcore(:,:),dim_1e,5,"Hcore")
!  call print_Mat(Fock(:,:,1),dim_1e,4,"Fock")

end subroutine do_guess


!#############################################
!#               Hückel Guess       
!#############################################
subroutine guess_huckel()
  use module_data,          only : Fock,Hcore,Spins,dim_1e,Sij
  use module_data,          only : Guess,atom,Bastype,Basis
  use module_io,            only : print_Mat
  use module_constants,     only : CoreQ,eV2Ha
  implicit none
  integer                       :: iSpin,i,j,Zi,Zj
  double precision              :: K = 1.75d0
!  double precision              :: Hij
  double precision              :: Par_ii(8,3)

  write(*,*) ""
  write(*,*) "    Extended Hückel guess..."
  write(*,*) ""

!  allocate(Hii(dim_1e))
  Par_ii(1,1) = -13.06d0
  Par_ii(1,2) = 0.0d0
  Par_ii(1,3) = 0.0d0
  Par_ii(2,1) = -23.836940d0
  Par_ii(2,2) = 0.0d0
  Par_ii(2,3) = 0.0d0
  Par_ii(3,1) = -64.465135d0
  Par_ii(3,2) = -5.41d0
  Par_ii(3,3) = -3.61d0
  Par_ii(4,1) = -122.009435d0
  Par_ii(4,2) = -9.33d0
  Par_ii(4,3) = -5.88d0
  Par_ii(5,1) = -197.703359d0
  Par_ii(5,2) = -14.0d0
  Par_ii(5,3) = -8.24d0
  Par_ii(6,1) = -297.452179d0
  Par_ii(6,2) = -19.42d0
  Par_ii(6,3) = -10.7d0
  Par_ii(7,1) = -417.519956d0
  Par_ii(7,2) = -25.58d0
  Par_ii(7,3) = -13.25d0
  Par_ii(8,1) = -552.668469d0
  Par_ii(8,2) = -32.49d0
  Par_ii(8,3) = -15.88d0

  Par_ii = Par_ii*ev2Ha

  do iSpin=1,Spins
    do i=1,dim_1e
      do j=1,dim_1e
!        Fock(i,j,iSpin) = Hcore(i,j)
        Fock(i,j,iSpin) = 0.0d0
        call CoreQ(atom(Basis(i)),Zi)
        call CoreQ(atom(Basis(j)),Zj)

!        if(Bastype(i) > 1 .and. Bastype(j) > 1) then
!          Fock(i,j,iSpin)  = K * Sij(i,j)    &
!                           * (Par_ii(Zi,Bastype(i)) + Par_ii(Zj,Bastype(j)))/2.0d0
!        endif
        if(i == j) then
          Fock(i,j,iSpin)  = K * Sij(i,j)    &
                           * (Par_ii(Zi,Bastype(i)) + Par_ii(Zj,Bastype(j)))/2.0d0
        endif
      enddo
    enddo
  enddo

!  call print_Mat(Hcore(:,:),dim_1e,5,"Hcore")
!  call print_Mat(Fock(:,:,1),dim_1e,4,"Fock")

end subroutine guess_huckel



!#############################################
!#           Calculate Fock Matrix 
!#############################################
subroutine calc_Fock()
  use module_data, only      : Dens,Fock,Spins,dim_1e,ERI,Hcore
  use module_data, only      : Sij,Par,Bastype,scaletype,Fockold
  use module_io, only        : print_Mat
  implicit none
  integer                         :: iSpin,i,j,k,l
  integer                         :: ij,kl,ijkl
  integer                         :: ik,jl,ikjl
!  double precision, allocatable   :: Fockold(:,:,:)
  double precision                :: ScaleFactor=1.0d0
  double precision                :: AddFactor=0.0d0
  double precision                :: TempPar

!  write(*,*) ""
!  write(*,*) "    Build new Fock..."
!  write(*,*) ""

!  allocate(Fockold(dim_1e,dim_1e,Spins))

  Fockold = Fock

  do iSpin=1,Spins
    Fock(1:dim_1e,1:dim_1e,iSpin) = Hcore(1:dim_1e,1:dim_1e)
  enddo

!$OMP PARALLEL PRIVATE(j,k,l,ij,kl,ik,jl,ijkl,ikjl,iSpin, &
!$OMP                  TempPar,ScaleFactor,AddFactor )
!$OMP DO
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

          ScaleFactor=1.0d0
          AddFactor=0.0d0

          if(scaletype==1)then
!            ScaleFactor = 1.0d0 - ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0)*Sij(i+1,j+1)
            if(Bastype(i+1)==4 .or. Bastype(j+1)==4)then
              TempPar = 0.0d0        
            elseif(Bastype(i+1)==3 .or. Bastype(j+1)==3)then
              TempPar = Par(3)
            elseif(Bastype(i+1)==2 .or. Bastype(j+1)==2)then
              TempPar = Par(2)
            elseif(Bastype(i+1)==1 .or. Bastype(j+1)==1)then
              TempPar = Par(1)
            endif
            ScaleFactor = 1.0d0 - (TempPar)*Sij(i+1,j+1)

          elseif(scaletype==2)then
            ScaleFactor = 1.0d0 -                                               &
                                  ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0) &
                                 *Sij(i+1,j+1) *Sij(i+1,k+1) *Sij(i+1,l+1)      &
                                 *Sij(j+1,k+1) *Sij(j+1,l+1) *Sij(k+1,l+1)
          elseif(scaletype==3)then
            ScaleFactor = 1.0d0 -                                                  &
                                  ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0)    &
                                 *(Sij(i+1,j+1) +Sij(i+1,k+1) +Sij(i+1,l+1)        &
                                 + Sij(j+1,k+1) +Sij(j+1,l+1) +Sij(k+1,l+1))/6.0d0
          elseif(scaletype==4)then
            ScaleFactor = 1.0d0 -                                                &
                                  ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0)  &
                              *abs(Sij(i+1,j+1) *Sij(i+1,k+1) *Sij(i+1,l+1)      &
                                 * Sij(j+1,k+1) *Sij(j+1,l+1) *Sij(k+1,l+1))
          elseif(scaletype==5)then
            ScaleFactor = 1.0d0 -                                                 &
                                  ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0)   &
                              *abs(Sij(i+1,j+1) *Sij(i+1,k+1) *Sij(i+1,l+1)       &
                                  *Sij(j+1,k+1) *Sij(j+1,l+1) *Sij(k+1,l+1))**(1.0d0/6.0d0)
          elseif(scaletype==6)then
            do iSpin=1,Spins
              ScaleFactor = 1.0d0 -                                                 &
                                    ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0)   &
                                   *(Dens(i+1,j+1,iSpin) *Dens(i+1,k+1,iSpin)       &
                                   * Dens(i+1,l+1,iSpin) *Dens(j+1,k+1,iSpin)       &
                                   * Dens(j+1,l+1,iSpin) *Dens(k+1,l+1,iSpin))
            enddo
          elseif(scaletype==7)then
            do iSpin=1,Spins
              ScaleFactor = 1.0d0 -                                                 &
                                    ((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0)   &
                                   *(Dens(i+1,j+1,iSpin) *Dens(i+1,k+1,iSpin)       &
                                   * Dens(i+1,l+1,iSpin) *Dens(j+1,k+1,iSpin)       & 
                                   * Dens(j+1,l+1,iSpin) *Dens(k+1,l+1,iSpin))      &
                                   *(Sij(i+1,j+1) *Sij(i+1,k+1) *Sij(i+1,l+1)       &
                                   * Sij(j+1,k+1) *Sij(j+1,l+1) *Sij(k+1,l+1))
            enddo
!          elseif(scaletype==8)then
!            ScaleFactor = 1.0d0 - (Par(1))*Sij(i+1,j+1)
!            AddFactor   = Par(2)*Sij(i+1,j+1)*Sij(k+1,l+1) + Par(3)
          elseif(scaletype==8)then
            ScaleFactor = 1.0d0 - (((Par(Bastype(i+1))+Par(Bastype(j+1)))/2.0d0))*Sij(i+1,j+1)*Sij(k+1,l+1)
!           AddFactor   = ((Par(Bastype(i+1))+Par(Bastype(i+1)))/2.0d0)*Sij(i+1,j+1)*Sij(k+1,l+1) 
          elseif(scaletype==9)then
!           AddFactor   =   (Par(1))*Sij(i+1,j+1)*Sij(k+1,l+1)
!           AddFactor   =   ((Par(Bastype(i+1))+Par(Bastype(i+1)))/2.0d0)*Sij(i+1,j+1)*Sij(k+1,l+1) 
!           write(*,*)  "",
!           AddFactor   =   -sqrt(abs(Par(Bastype(i+1))*Par(Bastype(j+1))))*Sij(i+1,j+1)*Sij(k+1,l+1)
!           AddFactor   = - Par(Bastype(i+1))*Par(Bastype(j+1))*Sij(i+1,j+1)*Sij(k+1,l+1)
           AddFactor   = -1.0d0*(Par(Bastype(i+1))*Par(Bastype(j+1)))*Sij(i+1,j+1)

          endif

          !UHF
          if(Spins==2)then
            do iSpin=1,Spins
              Fock(i+1,j+1,iSpin) = Fock(i+1,j+1,iSpin)            &
                                  +(Dens(k+1,l+1,1)*ERI(ijkl)      &
                                  + Dens(k+1,l+1,2)*ERI(ijkl)      &
                                  - Dens(k+1,l+1,iSpin)*ERI(ikjl)) &
                                  * ScaleFactor + AddFactor
            enddo
          !RHF
          elseif(Spins==1)then
            Fock(i+1,j+1,1)       = Fock(i+1,j+1,1)                                 &
                                  +(Dens(k+1,l+1,1)*(2.0d0*ERI(ijkl) - ERI(ikjl)))  &
                                  * ScaleFactor + AddFactor
          endif
        enddo
      enddo

    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  ! Damping    
!  if(DoDamp)then
!    Fock = Damp*Fock + (1.0d0-Damp)*Fockold
!  endif

!  do iSpin=1,Spins
!    call print_Mat(Fock(:,:,iSpin),dim_1e,4,"Fock")
!  enddo

!  deallocate(Fockold)

end subroutine calc_Fock


!#############################################
!#         Calculate Density Matrix 
!#############################################
subroutine calc_Dens()
  use module_data, only      : Dens,Coef,Spins,dim_1e,Densold
  use module_data, only      : noccA,noccB,Occ,DoDamp,Damp,DynDamp
  use module_io, only        : print_Mat
  implicit none
  integer                         :: iSpin,i,j,k,l
!  double precision, allocatable   :: Densold(:,:,:)

!  allocate(Densold(dim_1e,dim_1e,Spins))
  Drmsd = 0.0d0
  Densold = Dens

!  write(*,*) ""
!  write(*,*) "    Calc density matrix"
!  write(*,*) ""

  do iSpin=1,Spins
    do i=1,dim_1e
      do j=1,i
        Dens(i,j,iSpin) = 0.0d0
        do k=1,dim_1e
          Dens(i,j,iSpin) = Dens(i,j,iSpin) + Occ(k,iSpin)*Coef(i,k,iSpin)*Coef(j,k,iSpin)*0.5d0*Spins
        enddo
        Dens(j,i,iSpin) = Dens(i,j,iSpin)
        Drmsd = Drmsd + (Dens(i,j,iSpin)-Densold(i,j,iSpin))**2
      enddo
    enddo
!    call print_Mat(Dens(:,:,iSpin),dim_1e,5,"DensX")
  enddo

  ! Damping    
  if(DoDamp)then
    Dens = Damp*Dens + (1.0d0-Damp)*Densold
!    call calc_DynDamp()
  elseif(DynDamp)then
    call calc_DynDamp()
    Dens = Damp*Dens + (1.0d0-Damp)*Densold
  endif

  Drmsd = sqrt(Drmsd/dble(dim_1e*dim_1e*Spins))

end subroutine calc_Dens


!#############################################
!#       Calculate Dynamic Damping Factor
!#############################################
subroutine calc_DynDamp()
  use module_data, only      : Fock,Coef,dim_1e,Spins,Occ
  use module_data, only      : Dens,Densold,Hcore,Fockold,Damp
  implicit none
  double precision          :: E1a, E2a, E1b, E2b, E2ab, temp,tempa,tempb
  integer                   :: i,j,iSpin

  E1a  = 0.0d0
  E1b  = 0.0d0
  E2a  = 0.0d0
  E2b  = 0.0d0
  E2ab = 0.0d0

!        Eelec = Eelec + Dens(i,j,iSpin)                &
!                      *(Hcore(i,j) + Fock(i,j,iSpin))  &
!                      *(1.0d0/dble(Spins))

  do iSpin=1,Spins
    do i=1,dim_1e
      do j=1,dim_1e
        E1a   = E1a + Dens(i,j,iSpin)*(Hcore(i,j)*(1.0d0/dble(Spins)))
        E2a   = E2a + Dens(i,j,iSpin)*((Fock(i,j,iSpin))*(1.0d0/dble(Spins)))

        E1b   = E1b + Densold(i,j,iSpin)*(Hcore(i,j)*(1.0d0/dble(Spins)))
        E2b   = E2b + Densold(i,j,iSpin)*((Fockold(i,j,iSpin))*(1.0d0/dble(Spins)))

        E2ab  = E2ab + Dens(i,j,iSpin)*((Fockold(i,j,iSpin))*(1.0d0/dble(Spins)))

      enddo
    enddo
  enddo


  tempa = E1a + E2a
  tempb = E1b + E2b
  temp  = (E1a-E1b+E2ab-2.0d0*E2b)/(4.0d0*E2ab-2.0d0*E2a-2.0d0*E2b)
  damp = min(abs(temp),1.0d0)

!  write(*,*) E1a,"-",E1b,"+",E2ab,"-2",E2b
!  write(*,*) "4",E2ab,"-2",E2a,"-2",E2b
!  write(*,*) tempa,tempb,E2ab,temp 

end subroutine calc_DynDamp


!#############################################
!#          Orthogonalize Fock      
!#############################################
subroutine ort_fock()
  use module_data, only      : Fock,Coef,Sm12,dim_1e,Spins
  use module_io, only        : print_Mat
  implicit none
  integer                         :: i
  double precision, allocatable   :: temp(:,:)

!  write(*,*) ""
!  write(*,*) "    orthogonalizing Fock ..."
!  write(*,*) ""

  allocate(temp(dim_1e,dim_1e))

  do i=1,Spins
!  $::Coeff = transpose($::Sm12) x $::Fock x $::Sm12;
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                Fock(1:dim_1e,1:dim_1e,i),dim_1e,   &
                Sm12(1:dim_1e,1:dim_1e),dim_1e,      &
                0.0d0,temp(1:dim_1e,1:dim_1e),dim_1e)
    call dgemm('T','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                Sm12(1:dim_1e,1:dim_1e),dim_1e,      &
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
  use module_data, only      : Coef,Sm12,dim_1e,Spins
  use module_io, only        : print_Mat
  implicit none
  integer                         :: i
  double precision, allocatable   :: temp(:,:)

!  write(*,*) ""
!  write(*,*) "    deorthogonalizing Fock ..."
!  write(*,*) ""

  allocate(temp(dim_1e,dim_1e))

  do i=1,Spins
!  $::Coeff = $::Sm12 x $::Coeff;
    temp(1:dim_1e,1:dim_1e) = coef(1:dim_1e,1:dim_1e,i) 
    call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,  &
                Sm12(1:dim_1e,1:dim_1e),dim_1e,      &
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
!#              Diagonalize 1-RDM       
!#############################################
subroutine calc_natorbs()
! see Pulay + Hamilton JCP 88, 4926 (1988)
  use module_data, only      : Dens,NOCoef,NOEps,dim_1e,Spins,S12,nel,NoccA,NoccB
  use module_io, only        : print_Mat,print_Vec
  implicit none
  integer                         :: i,j,lwork,inf
  double precision, allocatable   :: work(:),temp(:,:)
  double precision                :: SpinQN,N,Na,Nb,sumocc2

  lwork = dim_1e*(3+dim_1e/2)
  allocate(work(lwork))
  allocate(NOCoef(1:dim_1e,1:dim_1e,Spins),NOEps(1:dim_1e,Spins),temp(1:dim_1e,1:dim_1e))

  NOCoef = 0.0d0
  NOEps  = 0.0d0
  temp   = 0.0d0

  write(*,*) ""
  write(*,*) "    Diagonalizing 1-RDM..."
  write(*,*) ""

  do i=1,Spins
    temp = temp + Dens(1:dim_1e,1:dim_1e,i)*2.0d0/dble(Spins) 
  enddo

  call build_S12()

! S12 = S12 x D x S12
!  call print_Mat(S12(:,:),dim_1e,3,"S12")
!  call print_Mat(temp(:,:),dim_1e,4,"temp")
! D x S12
  call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,temp,dim_1e,S12,dim_1e,0.0d0,NOCoef(:,:,1),dim_1e)
!  call print_Mat(NOCoef(:,:,1),dim_1e,7,"NO Coef")
  temp = NOCoef(:,:,1)
! S12 x temp
  call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,S12,dim_1e,temp,dim_1e,0.0d0,NOCoef(:,:,1),dim_1e)
!  call print_Mat(NOCoef(:,:,1),dim_1e,7,"NO Coef")

  deallocate(temp)

!  do i=1,Spins
!  call print_Mat(Dens(:,:,1),dim_1e,5,"1-RDM")
  call dsyev('V','U',dim_1e,NOCoef(1:dim_1e,1:dim_1e,1),dim_1e,NOEps(1:dim_1e,1),work,lwork,inf)
!  call print_Mat(NOCoef(:,:,1),dim_1e,7,"NO Coef")
!  enddo

!  do i=1,Spins
  call print_Vec(NOEps(:,1),dim_1e,15,"(U)NO occ. num.")
!  enddo

  sumocc2 = 0.0d0
  do i=1,dim_1e
    sumocc2 = sumocc2 + NOEps(i,1)*NOEps(i,1)
  enddo

  N = dble(nel)
  Na = dble(NoccA)
  Nb = dble(NoccB)

  SpinQN = N*(N+4.0d0)/4.0d0 - Na*Nb - 0.5d0*sumocc2

  write(*,*) ""
  write(*,*) "    <S**2> = ",SpinQN  
  write(*,*) ""

  deallocate(work)

end subroutine calc_natorbs


!#############################################
!#              Build S**1/2               
!#############################################
subroutine build_S12()
  use module_data, only      : S12,Sij,dim_1e
  implicit none
  integer                         :: i,j,lwork,inf
  double precision, allocatable   :: work(:)
  double precision, allocatable   :: S_vec(:,:)
  double precision, allocatable   :: S_val(:)
  double precision, allocatable   :: S_mat(:,:)
  double precision, allocatable   :: temp(:,:)

  lwork = dim_1e*(3+dim_1e/2)
  allocate(work(lwork),            &
           S_vec(dim_1e,dim_1e),   &
           S_mat(dim_1e,dim_1e),   &
           S_val(dim_1e),          &
           S12(1:dim_1e,1:dim_1e), &
           temp(dim_1e,dim_1e))

!  write(*,*) ""
!  write(*,*) "    Diagonalizing S..."
!  write(*,*) ""

  S_vec = Sij
  call dsyev('V','U',dim_1e,S_vec(1:dim_1e,1:dim_1e),dim_1e,S_val(1:dim_1e),work,lwork,inf)

!  call print_Vec(S_val(:),dim_1e,11,"Ovrlpvalues")

  deallocate(work)

  S_val = S_val**(0.5d0)
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

end subroutine build_S12


!#############################################
!#          Diagonalize Overlap     
!#############################################
subroutine dia_S()
  use module_data, only      : Sm12,Sij,dim_1e
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

! Sm12 = S_vec x S_mat x transpose(S_vec)
! temp = S_mat x transpose(S_vec) 
  call dgemm('N','T',dim_1e,dim_1e,dim_1e,1.0d0,S_mat,dim_1e,S_vec,dim_1e,0.0d0,temp,dim_1e)
! S12 = S_vec x temp 
  call dgemm('N','N',dim_1e,dim_1e,dim_1e,1.0d0,S_vec,dim_1e,temp,dim_1e,0.0d0,Sm12,dim_1e)

  deallocate(temp,S_vec,S_mat,S_val)

end subroutine dia_S


end module module_wavefun

