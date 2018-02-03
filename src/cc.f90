module module_cc
  implicit none

! Energies           
  double precision                            :: Elccd    = 0.0d0
  double precision, allocatable               :: t2(:,:,:,:),t2old(:,:,:,:)
  double precision, allocatable               :: SO_Occ(:)

!  integer                                     :: In4

contains

!#############################################
!#         Calculate Linearized CC
!#############################################
subroutine calc_Elccd
  use module_data,      only : MOI,Spins,dim_1e,nOccA,nOccB,Eps,SMO,Occ
  use module_data,      only : DropMO,DoDrop,Fock,DoSingles,dim_1e,AMO
  use module_data,      only : Eps_SO,F_SO
  use module_ints,      only : Index2e
  use module_wavefun,   only : Fock_to_MO,Fock_to_AO
  implicit none
  integer                   :: iSpin,MOZero,iT2
  integer                   :: p3,p4,h1,h2,h5,h6,p5,p6
  integer                   :: i,j,a,b 
  integer                   :: k,l,c,d


  integer                   :: nOcc,nmo
  double precision          :: occ_factor,Dijab,E_2,E_2OS,E_2SS
  double precision          :: DampCC = 0.5d0
  double precision          :: in1

  Elccd  = 0.0d0

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-LCCD correlation energy"
    write(*,*) "    "

    nmo = dim_1e*2
    MOZero = 1

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif
    nOcc =  nOccA*2

    allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
          t2old(1:nmo,1:nmo,1:nmo,1:nmo))

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0

  do iT2=1,30 ! T amplitude loop
    Elccd  = 0.0d0
    t2old  = t2 !R0.5d0*(t2+t2old) ! uses mixing
    !T2 equations
    do i=MOZero,nOcc
      do j=MOZero,nOcc 
        do a=nOcc+1,nmo
          do b=nOcc+1,nmo
            Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j)
            !write(*,*) i,j,a,b,in1

!           [ - 1.0 + 1.0 * P( p4 p3 h1 h2 => p3 p4 h1 h2 ) ] * Sum ( p5 ) * f ( p4 p5 ) * t ( p5 p3 h1 h2 ) done
            !do d = nOcc+1,nmo
            !  if ((d == a) .or. (d == b)) then
            !      cycle 
            !  endif
            !  in1 = in1 - 1.0d0 * F_SO(b,d) * t2old(a,d,i,j)
            !  in1 = in1 + 1.0d0 * F_SO(a,d) * t2old(d,b,i,j)
            !enddo

!           [ - 1.0 + 1.0 * P( p3 p4 h1 h2 => p3 p4 h2 h1 ) ] * Sum ( h5 ) * f ( h5 h1 ) * t ( p3 p4 h5 h2 ) done
            !do k = MOZero,nOcc
            !  if ((k == i) .or. (k == j)) then
            !      cycle
            !  endif
            !  in1 = in1 - 1.0d0 * F_SO(k,j) * t2old(a,b,i,k)
            !  in1 = in1 + 1.0d0 * F_SO(k,i) * t2old(a,b,k,j)
            !enddo

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do c = nOcc+1,nmo
              do d = nOcc+1,c-1
                in1 = in1 + 1.0d0 * t2old(c,d,i,j) * AMO(a,b,c,d)
              enddo
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do k = MOZero,nOcc
              do l = MOZero,k-1
                in1 = in1 + 1.0d0 * t2old(a,b,k,l) * AMO(k,l,i,j)
              enddo
            enddo

!          [ + 1.0 
!            - 1.0 * P( p3 p4 h2 h1 => p4 p3 h2 h1 ) 
!            - 1.0 * P( p3 p4 h2 h1 => p3 p4 h1 h2 ) 
!            + 1.0 * P( p3 p4 h2 h1 => p4 p3 h1 h2 ) ] 
!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do c = nOcc+1,nmo
              do k = MOZero,nOcc
                in1 = in1 - 1.0d0 * t2old(a,c,i,k) * AMO(k,b,j,c)
                in1 = in1 + 1.0d0 * t2old(b,c,i,k) * AMO(k,a,j,c)
                in1 = in1 + 1.0d0 * t2old(a,c,j,k) * AMO(k,b,i,c)
                in1 = in1 - 1.0d0 * t2old(b,c,j,k) * AMO(k,a,i,c)
              enddo
            enddo

!            write(*,*) "t2   ",p3,p4,h1,h2,AMO(p3,p4,h1,h2),AMO(h1,h2,p3,p4)
!            Dijab = (Eps_SO(p3) + Eps_SO(p4) - Eps_SO(h1) - Eps_SO(h2))
            t2(a,b,i,j) =  1.0d0*in1/Dijab

          enddo
        enddo
      enddo
    enddo

    Elccd = 0.0d0
!    E_2 = 0.0d0
!    E_2OS = 0.0d0
!    E_2SS = 0.0d0

    do i=MOZero,nOcc
      do j=MOZero,i-1
        do a=nOcc+1,nmo
          do b=nOcc+1,a-1
            Elccd = Elccd + 1.0d0*t2(a,b,i,j)*AMO(i,j,a,b)
          enddo
        enddo
      enddo
    enddo

    write(*,*) "      LCCD: ",iT2,Elccd 
 
  enddo ! t2 loop

  elseif(spins == 2)then
    write(*,*) "    "
    write(*,*) "    Calculating UHF-LCCD correlation energy"
    write(*,*) "    "

    allocate(SO_Occ(1:dim_1e*2))

    nmo = dim_1e*2
    MOZero = 1

    do i=2,nmo,2
      SO_Occ(i-1) = Occ(i/2,1)
      SO_Occ(i) = Occ(i/2,2)
    enddo

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif
    nOcc =  nOccA*2

    allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
          t2old(1:nmo,1:nmo,1:nmo,1:nmo))

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0

  do iT2=1,30 ! T amplitude loop
    Elccd  = 0.0d0
    t2old  = t2 !R0.5d0*(t2+t2old) ! uses mixing

    !T2 equations
    !Nocc version
    do i=MOZero,nmo
      do j=MOZero,nmo
        do a=MOZero,nmo
          do b=MOZero,nmo
            if((SO_Occ(a) /= 0.0) .or. (SO_Occ(b) /= 0.0))then
              cycle
            elseif((SO_Occ(i) /= 1.0) .or. (SO_Occ(j) /= 1.0))then
              cycle
            endif

            Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j)

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do c = MOZero,nmo
              do d = MOZero,c-1
                if((SO_Occ(c) /= 0.0) .or. (SO_Occ(d) /= 0.0))then
                  cycle
                endif
                in1 = in1 + 1.0d0 * t2old(c,d,i,j) * AMO(a,b,c,d)
              enddo
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do k = MOZero,nmo
              do l = MOZero,k-1
                if((SO_Occ(k) /= 1.0) .or. (SO_Occ(l) /= 1.0))then
                  cycle
                endif
                in1 = in1 + 1.0d0 * t2old(a,b,k,l) * AMO(k,l,i,j)
              enddo
            enddo

!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do c = MOZero,nmo
              do k = MOZero,nmo
                if((SO_Occ(c) /= 0.0) .or. (SO_Occ(k) /= 1.0))then
                  cycle
                endif
                in1 = in1 - 1.0d0 * t2old(a,c,i,k) * AMO(k,b,j,c)
                in1 = in1 + 1.0d0 * t2old(b,c,i,k) * AMO(k,a,j,c)
                in1 = in1 + 1.0d0 * t2old(a,c,j,k) * AMO(k,b,i,c)
                in1 = in1 - 1.0d0 * t2old(b,c,j,k) * AMO(k,a,i,c)
              enddo
            enddo

            t2(a,b,i,j) =  1.0d0*in1/Dijab

          enddo
        enddo
      enddo
    enddo

    Elccd = 0.0d0
!    E_2 = 0.0d0
!    E_2OS = 0.0d0
!    E_2SS = 0.0d0

    do i=MOZero,nOcc
      do j=MOZero,i-1
        do a=nOcc+1,nmo
          do b=nOcc+1,a-1
            Elccd = Elccd + 1.0d0*t2(a,b,i,j)*AMO(i,j,a,b)
          enddo
        enddo
      enddo
    enddo

    write(*,*) "      LCCD: ",iT2,Elccd

  enddo ! t2 loop
  endif


contains

!#############################################
!#              Calc 2e Indices                 
!#############################################
function In4(i,j,k,l) 
  implicit none
  integer, intent(in)           :: i,j,k,l
  integer                       :: In4
  integer                       :: ij,kl

  if (i>j)then
    ij = i*(i+1)/2 + j
  else
    ij = j*(j+1)/2 + i
  endif

  if (k>l)then
    kl = k*(k+1)/2 + l
  else
    kl = l*(l+1)/2 + k
  endif

  if (ij>kl)then
    In4 = ij*(ij+1)/2 + kl
  else
    In4 = kl*(kl+1)/2 + ij
  endif

end function In4


end subroutine calc_Elccd


end module module_cc
