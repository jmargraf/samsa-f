module module_cis
  implicit none

! Energies           
  double precision, allocatable               :: Eexcite(:)   
  double precision, allocatable               :: matA(:,:)

contains

!#############################################
!#      Calculate CIS excitation energies
!#############################################
subroutine calc_CIS
  use module_data,      only : Spins,dim_1e,nOccA,nOccB,Eps,Occ
  use module_data,      only : DropMO,DoDrop,DoSingles,dim_1e,AMO
  use module_data,      only : Eps_SO,F_SO
  implicit none
  logical                   :: CCD
  integer                   :: nOcc,nmo,dimCIS
  integer                   :: iSpin,MOZero,iT2
  integer                   :: i,j,a,b
  integer                   :: ia,jb


  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-CIS excitation energies"
    write(*,*) "    "

    nmo = dim_1e*2
    MOZero = 1

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif

    nOcc =  nOccA*2
    dimCIS = (nOcc-DropMO*2)*(nmo-nOcc)


    allocate(matA(1:dimCIS,1:dimCIS),  &
       Eexcite(1:dimCIS))

    matA = 0.0d0

    ia = 0
    jb = 0

    write(*,*) '  building CIS matrix '

    ! Build CIS matrix
    do i = MOZero,nOcc
      do a = nOcc+1,nmo
        ia = ia + 1
        do j = MOZero,nOcc
          do b = nOcc+1,nmo
            jb = jb + 1
            if(ia == jb)then
              matA(ia,jb) = Eps_SO(a)- Eps_SO(i) 
            endif
            matA(ia,jb) = matA(ia,jb) + AMO(a,j,i,b)
          enddo
        enddo
        jb = 0
      enddo
    enddo

    write(*,*) '  done building CIS matrix '

    call dia_CI(dimCIS)

    ia = 0
    do i = MOZero,nOcc
      do a = nOcc+1,nmo
        ia = ia + 1
        write(*,*) ia, Eexcite(ia)
      enddo
    enddo

  endif

end subroutine calc_CIS


!#############################################
!#          Diagonalize CI Matrix
!#############################################
subroutine dia_CI(dimCIS)
  use module_io, only        : print_Mat,print_Vec
  implicit none
  integer                         :: i,j,lwork,inf,dimCIS
  double precision, allocatable   :: work(:)

  lwork = dimCIS*(3+dimCIS/2)
  allocate(work(lwork))

  write(*,*) ""
  write(*,*) "    Diagonalizing CI Matrix"
  write(*,*) ""

  call dsyev('V','U',dimCIS,matA(1:dimCIS,1:dimCIS),dimCIS,Eexcite(1:dimCIS),work,lwork,inf)
  !call print_Mat(matA(:,:),dimCIS,6,"CIS WF")

  !call print_Vec(Eexcite(:),dimCIS,11,"Eigenvalues")

  deallocate(work)

end subroutine dia_CI


end module module_cis 
