module module_trans
  implicit none

contains

!#############################################
!#         Full Integral Transformation
!#############################################
subroutine trans_full 
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

  write(*,*) "    "
  write(*,*) "    Performing full integral transformation... "
  write(*,*) "    "

  dimMO = (spins*dim_1e) - 1
  dimAO = dim_1e - 1

  allocate(in_pqrl(0:dimAO,0:dimAO,0:dimAO,0:dimMO))
  in_pqrl = 0.0d0

  ! pqrs -> pqrl
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

  write(*,*) "      25% ..."

  allocate(in_pqkl(0:dimAO,0:dimAO,0:dimMO,0:dimMO))
  in_pqkl = 0.0d0

  ! pqrl -> pqkl
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

  deallocate(in_pqrl)

  write(*,*) "      50% ..."
  
  allocate(in_pjkl(0:dimAO,0:dimMO,0:dimMO,0:dimMO))
  in_pjkl = 0.0d0

  ! pqkl -> pjkl
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

  deallocate(in_pqkl)

  write(*,*) "      75% ..."

  allocate(in_ijkl(0:dimMO,0:dimMO,0:dimMO,0:dimMO))
  in_ijkl = 0.0d0

  ! pjkl -> ijkl
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

  deallocate(in_pjkl)

  write(*,*) "     100% ..."
  write(*,*) "    "

  if(spins == 1)then
    allocate(MOI(0:dim_2e-1))

    ! sort integrals into 1D array
    do i=0,dimAO
      do j=0,i
        call Index2e(i,j,ij)
        do k=0,dimAO
          do l=0,k
            call Index2e(k,l,kl)
            call Index2e(ij,kl,ijkl)
            MOI(ijkl) = in_ijkl(i,j,k,l)       
          enddo
        enddo
      enddo
    enddo

  elseif(spins == 2)then
    dimMO = dim_1e*spins
    allocate(SMO(1:dimMO,1:dimMO,1:dimMO,1:dimMO))

    do i=1,dimMO
      do j=1,dimMO
        do k=1,dimMO
          do l=1,dimMO
            SMO(i,j,k,l) = in_ijkl(i-1,j-1,k-1,l-1)
          enddo
        enddo
      enddo
    enddo

  endif
    
  deallocate(in_ijkl)

end subroutine trans_full

end module module_trans

