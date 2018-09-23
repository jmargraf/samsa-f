module module_trans
  implicit none
  logical                       :: transdone = .false. 

contains

!#############################################
!#         Full Integral Transformation
!#############################################
subroutine trans_full(doprint) 
  use module_ints, only          : Index2e
  use module_data, only          : dim_1e,dim_2e,ERI,Coef,MOI,spins,SMO
  implicit none
  integer                       :: i,j,k,l,iSpin
  integer                       :: ij,kl,ijkl
  integer                       :: p,q,r,s,pq,rs,pqrs
  integer                       :: dimAO,dimMO,ab

  double precision, allocatable :: in_pqrl(:,:,:,:)
  double precision, allocatable :: in_pqkl(:,:,:,:)
  double precision, allocatable :: in_pjkl(:,:,:,:)
  double precision, allocatable :: in_ijkl(:,:,:,:)
  real                          :: starttime,stoptime,time
  logical, intent(IN)           :: doprint

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Performing full integral transformation... "
    write(*,*) "    "
  endif

  starttime = secnds(0.0)
!  call cpu_time(starttime)
!!$  call omp_get_wtime(starttime)

  dimMO = (spins*dim_1e) - 1
  dimAO = dim_1e - 1

  allocate(in_pqrl(0:dimAO,0:dimAO,0:dimAO,0:dimMO))
  in_pqrl = 0.0d0

! pqrs -> pqrl
!!$ write(*,*) "      (in parallel)" !, OMP_NUM_THREADS
  if(doprint)then
    write(*,*) "      (pq|rs) -> (pq|rl)"
  endif
!$OMP PARALLEL PRIVATE(p,q,r,s,pq,rs,pqrs)
!$OMP DO
  do l=0,dimMO
    do p=0,dimAO
      do q=0,dimAO
        do r=0,dimAO
          do s=0,dimAO
            call Index2e(p,q,pq)
            call Index2e(r,s,rs)
            call Index2e(pq,rs,pqrs)
            pqrs = pqrs + 1
            if(spins == 2)then
              if(mod(l,2) == 0)then
                in_pqrl(p,q,r,l) = in_pqrl(p,q,r,l) + Coef(s+1,(l/2+1),1)*ERI(pqrs) 
              else
                in_pqrl(p,q,r,l) = in_pqrl(p,q,r,l) + Coef(s+1,((l-1)/2+1),2)*ERI(pqrs)
              endif
            elseif(spins == 1)then
              in_pqrl(p,q,r,l) = in_pqrl(p,q,r,l) + Coef(s+1,l+1,1)*ERI(pqrs)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL


  allocate(in_pqkl(0:dimAO,0:dimAO,0:dimMO,0:dimMO))
  in_pqkl = 0.0d0

  ! pqrl -> pqkl
  if(doprint)then
    write(*,*) "      (pq|rl) -> (pq|kl)"
  endif
!$OMP PARALLEL PRIVATE(p,q,r,l)
!$OMP DO
  do k=0,dimMO
    do p=0,dimAO
      do q=0,dimAO
        do r=0,dimAO
          do l=0,dimMO
            if(spins == 2)then
              if(mod(k,2) == 0)then
                in_pqkl(p,q,k,l) = in_pqkl(p,q,k,l) + Coef(r+1,(k/2+1),1)     * in_pqrl(p,q,r,l)
              else
                in_pqkl(p,q,k,l) = in_pqkl(p,q,k,l) + Coef(r+1,((k-1)/2+1),2) * in_pqrl(p,q,r,l)
              endif
            elseif(spins == 1)then
                in_pqkl(p,q,k,l) = in_pqkl(p,q,k,l) + Coef(r+1,k+1,1)         * in_pqrl(p,q,r,l)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  deallocate(in_pqrl)

  allocate(in_pjkl(0:dimAO,0:dimMO,0:dimMO,0:dimMO))
  in_pjkl = 0.0d0

  ! pqkl -> pjkl
  if(doprint)then
    write(*,*) "      (pq|kl) -> (pj|kl)"
  endif
!$OMP PARALLEL PRIVATE(p,q,k,l)
!$OMP DO
  do j=0,dimMO
    do p=0,dimAO
      do q=0,dimAO
        do k=0,dimMO
          do l=0,dimMO
            if(spins == 2)then
              if(mod(j,2) == 0)then
                in_pjkl(p,j,k,l) = in_pjkl(p,j,k,l) + Coef(q+1,(j/2+1),1)     * in_pqkl(p,q,k,l)
              else
                in_pjkl(p,j,k,l) = in_pjkl(p,j,k,l) + Coef(q+1,((j-1)/2+1),2) * in_pqkl(p,q,k,l)
              endif
            elseif(spins == 1)then
                in_pjkl(p,j,k,l) = in_pjkl(p,j,k,l) + Coef(q+1,j+1,1)         * in_pqkl(p,q,k,l)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  deallocate(in_pqkl)

  allocate(in_ijkl(0:dimMO,0:dimMO,0:dimMO,0:dimMO))
  in_ijkl = 0.0d0

  ! pjkl -> ijkl
  if(doprint)then
    write(*,*) "      (pj|kl) -> (ij|kl)"
  endif
!$OMP PARALLEL PRIVATE(p,j,k,l)
!$OMP DO
  do i=0,dimMO
    do p=0,dimAO
      do j=0,dimMO
        do k=0,dimMO
          do l=0,dimMO
            if(spins == 2)then
              if(mod(i,2) == 0)then
                in_ijkl(i,j,k,l) = in_ijkl(i,j,k,l) + Coef(p+1,(i/2+1),1)     * in_pjkl(p,j,k,l)
              else
                in_ijkl(i,j,k,l) = in_ijkl(i,j,k,l) + Coef(p+1,((i-1)/2+1),2) * in_pjkl(p,j,k,l)
              endif
            elseif(spins == 1)then
                in_ijkl(i,j,k,l) = in_ijkl(i,j,k,l) + Coef(p+1,i+1,1)         * in_pjkl(p,j,k,l)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  deallocate(in_pjkl)

  if(doprint)then
    write(*,*) "      done ..."
    write(*,*) "    "
  endif

  if(spins == 1)then
    if(.not.allocated(MOI))then
      allocate(MOI(0:dim_2e-1))
    endif
! sort integrals into 1D array

!$OMP PARALLEL PRIVATE(j,k,l,ij,kl,ijkl)
!$OMP DO
    do i=0,dimAO
      do j=0,i
        call Index2e(i,j,ij)
        do k=0,dimAO
          do l=0,k
            call Index2e(k,l,kl)
            call Index2e(ij,kl,ijkl)
!            write(*,*) ijkl
            MOI(ijkl) = in_ijkl(i,j,k,l)       
            !write(*,*) i,j,k,l,MOI(ijkl)
          enddo
        enddo
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

  elseif(spins == 2)then
    dimMO = dim_1e*spins
    if(.not.allocated(smo))then
      allocate(SMO(1:dimMO,1:dimMO,1:dimMO,1:dimMO))
    endif

!$OMP PARALLEL PRIVATE(j,k,l)
!$OMP DO
    do i=1,dimMO
      do j=1,dimMO
        do k=1,dimMO
          do l=1,dimMO
            SMO(i,j,k,l) = in_ijkl(i-1,j-1,k-1,l-1)
            !write(*,*) i,j,k,l,SMO(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

  endif
    
  deallocate(in_ijkl)

  transdone = .true.

  time=secnds(starttime)
!  call cpu_time(stoptime)
!!$  call omp_get_wtime(stoptime)
!  time = stoptime-starttime
  if(doprint)then
    write(*,*) "      Integral transformation done in ",time," s"
  endif

end subroutine trans_full


! translate integrals to spin integrated antisymmetrized form
! the number of spin orbitals nmo is 2*n, with alpha and beta spins alternating
subroutine trans_ucc(doprint)
  use module_ints, only          : Index2e
  use module_data, only          : dim_1e,dim_2e,ERI,Coef,MOI,spins,SMO,AMO
  use module_data, only          : Fock,Eps,Eps_SO,F_SO
  use module_wavefun, only       : Fock_to_MO
  use module_io, only            : print_Vec
  implicit none
  integer                       :: p,q,r,s,pr,qs,prqs,ps,qr,psqr
  integer                       :: nmo,nso

  real                          :: starttime,stoptime,time
  logical, intent(IN)           :: doprint
  double precision              :: value1,value2 

  if(doprint)then
    write(*,*) "    "
    write(*,*) "    Converting integrals to spin MO basis ... "
  endif

  if(spins==1)then

    nmo = dim_1e*2

    if(.not.allocated(AMO))then
      allocate(AMO(1:nmo,1:nmo,1:nmo,1:nmo),Eps_SO(1:nmo),F_SO(1:nmo,1:nmo))
    endif

    AMO = 0.0d0

    call Fock_to_MO()

    do p=0,nmo-1
      Eps_SO(p+1) = Eps(p/2+1,1)
      do q=0,nmo-1
        F_SO(p+1,q+1) = 0.0d0
        if(mod(p,2) == mod(q,2))then
          F_SO(p+1,q+1) =  Fock(p/2+1,q/2+1,1)
        endif
        do r=0,nmo-1
          do s=0,nmo-1
            call Index2e(p/2,r/2,pr)
            call Index2e(q/2,s/2,qs)
            call Index2e(p/2,s/2,ps)
            call Index2e(q/2,r/2,qr)
            call Index2e(pr,qs,prqs)
            call Index2e(ps,qr,psqr)
            value1 = 0.0d0
            value2 = 0.0d0

            if((mod(p,2) == mod(r,2)) .and. (mod(q,2) == mod(s,2)))then
              value1 = MOI(prqs)
            endif

            if((mod(p,2) == mod(s,2)) .and. (mod(q,2) == mod(r,2)))then
              value2 = MOI(psqr)
            endif

            AMO(p+1,q+1,r+1,s+1) = value1 - value2

!            write(*,*) p,q,r,s,value1,value2,AMO(p+1,q+1,r+1,s+1)

          enddo
        enddo
      enddo
    enddo

    call print_Vec(Eps_SO,nmo,6,"Eps_SO")

  elseif(spins==2)then

    nmo = dim_1e*2

    if(.not.allocated(AMO))then
      allocate(AMO(1:nmo,1:nmo,1:nmo,1:nmo),Eps_SO(1:nmo),F_SO(1:nmo,1:nmo))
    endif

    AMO = 0.0d0

    call Fock_to_MO()

    do p=1,nmo
      if((mod(p,2) == 1))then
        Eps_SO(p) = Eps((p+1)/2,1)
      else
        Eps_SO(p) = Eps(p/2,2)
      endif

      do q=1,nmo
        F_SO(p,q) = 0.0d0
        if(mod(p,2) == mod(q,2) .and. mod(p,2) == 1)then
          F_SO(p,q) =  Fock((p+1)/2,(q+1)/2,1)
        elseif(mod(p,2) == mod(q,2) .and. mod(p,2) == 0)then
          F_SO(p,q) =  Fock(p/2,q/2,2)
        endif

        do r=1,nmo
          do s=1,nmo
            value1 = 0.0d0
            value2 = 0.0d0

            if((mod(p,2) == mod(r,2)) .and. (mod(q,2) == mod(s,2)))then
              value1 = SMO(p,r,q,s)
            endif

            if((mod(p,2) == mod(s,2)) .and. (mod(q,2) == mod(r,2)))then
              value2 = SMO(p,s,q,r)
            endif

            AMO(p,q,r,s) = value1 - value2

          enddo
        enddo
      enddo
    enddo

    call print_Vec(Eps_SO,nmo,6,"Eps_SO")

  endif

  if(doprint)then
    write(*,*) "    ... done "
    write(*,*) "    "
  endif

end subroutine trans_ucc


end module module_trans

