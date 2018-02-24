module module_cc
  implicit none

! Energies           
  double precision                            :: Eccd    = 0.0d0
  double precision, allocatable               :: t2(:,:,:,:),t2old(:,:,:,:)
  double precision, allocatable               :: SO_Occ(:)

contains

!#############################################
!#         Calculate CC doubles
!#############################################
subroutine calc_Eccd(CCD,frac_occ)
  use module_data,      only : Spins,dim_1e,nOccA,nOccB,Eps,Occ
  use module_data,      only : DropMO,DoDrop,DoSingles,dim_1e,AMO
  use module_data,      only : Eps_SO,F_SO
  implicit none
  logical                   :: CCD 
  logical                   :: frac_occ
  logical                   :: MBPT3 = .false.
  integer                   :: iSpin,MOZero,iT2
  integer                   :: i,j,a,b 
  integer                   :: k,l,c,d

  integer                   :: nOcc,nmo
  integer                   :: maxT = 150
  double precision          :: occ_factor,Dijab
  double precision          :: DampCC = 0.5d0
  double precision          :: Ecctol = 1.0d-7
  double precision          :: Eold   = 0.0d0
  double precision          :: deltaE = 99.0d0
  double precision          :: in1 ! intermediate during T Amp equations

  Eccd  = 0.0d0

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-(L)CCD correlation energy"
    write(*,*) "    "

    nmo = dim_1e*2
    MOZero = 1

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif
    nOcc =  nOccA*2

    if(.not.allocated(t2))then
      allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
            t2old(1:nmo,1:nmo,1:nmo,1:nmo))
    endif

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0
    Eccd = 0.0d0


  do iT2=1,maxT ! T amplitude loop
    Eold  = Eccd
    t2old  = t2 !R0.5d0*(t2+t2old) ! uses mixing

    !T2 equations
!$OMP PARALLEL PRIVATE(i,j,a,b,k,l,c,d,Dijab,in1)
!$OMP DO
    do i=MOZero,nOcc
      do j=MOZero,nOcc 
        do a=nOcc+1,nmo
          do b=nOcc+1,nmo
            Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j)
            !write(*,*) i,j,a,b,in1

!           [ - 1.0 + 1.0 * P( p4 p3 h1 h2 => p3 p4 h1 h2 ) ] * Sum ( p5 ) * f ( p4 p5 ) * t ( p5 p3 h1 h2 ) done
            do d = nOcc+1,nmo
              if ((d == a) .or. (d == b)) then
                  cycle 
              endif
              in1 = in1 - 1.0d0 * F_SO(b,d) * t2old(a,d,i,j)
              in1 = in1 + 1.0d0 * F_SO(a,d) * t2old(d,b,i,j)
            enddo

!           [ - 1.0 + 1.0 * P( p3 p4 h1 h2 => p3 p4 h2 h1 ) ] * Sum ( h5 ) * f ( h5 h1 ) * t ( p3 p4 h5 h2 ) done
            do k = MOZero,nOcc
              if ((k == i) .or. (k == j)) then
                  cycle
              endif
              in1 = in1 - 1.0d0 * F_SO(k,j) * t2old(a,b,i,k)
              in1 = in1 + 1.0d0 * F_SO(k,i) * t2old(a,b,k,j)
            enddo

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

            if(CCD)then ! do full CCD
            !Terms that are quadratic in T2
              do c = nOcc+1,nmo
                do d = nOcc+1,nmo
                  do k = MOZero,nOcc
                    do l = MOZero,nOcc

                      in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(b,c,i,k)*t2old(d,a,l,j)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,j,k)*t2old(d,a,l,i)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(a,c,j,k)*t2old(d,b,l,i)   * AMO(k,l,c,d)

                      in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)

                      in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,i,j)*t2old(a,d,k,l)   * AMO(k,l,c,d)

                      in1 = in1 - 0.5d0 *  t2old(a,b,i,k)*t2old(c,d,j,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(a,b,j,k)*t2old(c,d,i,l)   * AMO(k,l,c,d)

                    enddo
                  enddo
                enddo
              enddo
            endif


!            write(*,*) "t2   ",p3,p4,h1,h2,AMO(p3,p4,h1,h2),AMO(h1,h2,p3,p4)
!            Dijab = (Eps_SO(p3) + Eps_SO(p4) - Eps_SO(h1) - Eps_SO(h2))
            t2(a,b,i,j) =  1.0d0*in1/Dijab

          enddo
        enddo
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    Eccd = 0.0d0

    do i=MOZero,nOcc
      do j=MOZero,i-1
        do a=nOcc+1,nmo
          do b=nOcc+1,a-1
            Eccd = Eccd + 1.0d0*t2(a,b,i,j)*AMO(i,j,a,b)
          enddo
        enddo
      enddo
    enddo

    deltaE = Eccd-Eold
    if(CCD)then
      write(*,*) "      CCD: ",iT2,Eccd,deltaE
    else
      write(*,*) "      LCCD: ",iT2,Eccd,deltaE 
    endif
    if(abs(deltaE)<Ecctol)then
      write(*,*) '    '
      write(*,*) '    CC iterations converged '
      write(*,*) '    '
      exit
    endif

  enddo ! t2 loop

  elseif(spins == 2 .and. .not. frac_occ)then
    write(*,*) "    "
    write(*,*) "    Calculating UHF-(L)CCD correlation energy"
    write(*,*) "    "

    if(.not.allocated(SO_Occ))then
      allocate(SO_Occ(1:dim_1e*2))
    endif

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

    if(.not.allocated(t2))then
      allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
          t2old(1:nmo,1:nmo,1:nmo,1:nmo))
    endif

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0
    Eccd = 0.0d0

  do iT2=1,maxT ! T amplitude loop
    Eold   = Eccd
    t2old  = t2 !R0.5d0*(t2+t2old) ! uses mixing

    !T2 equations
    !Nocc version
!$OMP PARALLEL PRIVATE(i,j,a,b,k,l,c,d,Dijab,in1)
!$OMP DO
    do i=MOZero,nmo
      do j=MOZero,nmo
        do a=MOZero,nmo
          do b=MOZero,nmo
            if((SO_Occ(a) == 1.0) .or. (SO_Occ(b) == 1.0))then
              cycle
            elseif((SO_Occ(i) == 0.0) .or. (SO_Occ(j) == 0.0))then
              cycle
            endif

            Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j)

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do c = MOZero,nmo
              do d = MOZero,c-1
                if((SO_Occ(c) == 1.0) .or. (SO_Occ(d) == 1.0))then
                  cycle
                endif
                in1 = in1 + 1.0d0 * t2old(c,d,i,j) * AMO(a,b,c,d)
              enddo
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do k = MOZero,nmo
              do l = MOZero,k-1
                if((SO_Occ(k) == 0.0) .or. (SO_Occ(l) == 0.0))then
                  cycle
                endif
                in1 = in1 + 1.0d0 * t2old(a,b,k,l) * AMO(k,l,i,j)
              enddo
            enddo

!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do c = MOZero,nmo
              do k = MOZero,nmo
                if((SO_Occ(c) == 1.0) .or. (SO_Occ(k) == 0.0))then
                  cycle
                endif
                in1 = in1 - 1.0d0 * t2old(a,c,i,k) * AMO(k,b,j,c)
                in1 = in1 + 1.0d0 * t2old(b,c,i,k) * AMO(k,a,j,c)
                in1 = in1 + 1.0d0 * t2old(a,c,j,k) * AMO(k,b,i,c)
                in1 = in1 - 1.0d0 * t2old(b,c,j,k) * AMO(k,a,i,c)
              enddo
            enddo

            if(CCD)then ! do full CCD
            !Terms that are quadratic in T2
              do c = MOZero,nmo
                do d = MOZero,nmo
                  do k = MOZero,nmo
                    do l = MOZero,nmo
                      if((SO_Occ(k) == 0.0) .or. (SO_Occ(l) == 0.0))then
                        cycle
                      elseif((SO_Occ(c) == 1.0) .or. (SO_Occ(d) == 1.0))then
                        cycle
                      endif

                      in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(b,c,i,k)*t2old(d,a,l,j)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,j,k)*t2old(d,a,l,i)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(a,c,j,k)*t2old(d,b,l,i)   * AMO(k,l,c,d)

                      in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)

                      in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,i,j)*t2old(a,d,k,l)   * AMO(k,l,c,d)

                      in1 = in1 - 0.5d0 *  t2old(a,b,i,k)*t2old(c,d,j,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(a,b,j,k)*t2old(c,d,i,l)   * AMO(k,l,c,d)

                    enddo
                  enddo
                enddo
              enddo
            endif

            t2(a,b,i,j) =  1.0d0*in1/Dijab

          enddo
        enddo
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    Eccd = 0.0d0

    do i=MOZero,nmo
      do j=MOZero,i-1
        do a=MOZero,nmo
          do b=MOZero,a-1
            if((SO_Occ(a) == 1.0) .or. (SO_Occ(b) == 1.0))then
              cycle
            elseif((SO_Occ(i) == 0.0) .or. (SO_Occ(j) == 0.0))then
              cycle
            endif
            Eccd = Eccd + 1.0d0*t2(a,b,i,j)*AMO(i,j,a,b)
          enddo
        enddo
      enddo
    enddo

    deltaE = Eccd-Eold
    if(CCD)then
      write(*,*) "      CCD: ",iT2,Eccd,deltaE
    else
      write(*,*) "      LCCD: ",iT2,Eccd,deltaE 
    endif
    if(abs(deltaE)<Ecctol)then
      write(*,*) '    '
      write(*,*) '    CC iterations converged '
      write(*,*) '    '
      exit
    endif

  enddo ! t2 loop

  elseif(spins == 2 .and. frac_occ)then
    write(*,*) "    "
    write(*,*) "    Calculating frac. occ. UHF-(L)CCD correlation energy"
    write(*,*) "    "

    if(.not.allocated(SO_Occ))then
      allocate(SO_Occ(1:dim_1e*2))
    endif

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

    if(.not.allocated(t2))then
      allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
          t2old(1:nmo,1:nmo,1:nmo,1:nmo))
    endif

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0
    Eccd = 0.0d0

  do iT2=1,maxT ! T amplitude loop
    Eold   = Eccd
    t2old  = t2 !R0.5d0*(t2+t2old) ! uses mixing

    !T2 equations
    !Nocc version
!$OMP PARALLEL PRIVATE(i,j,a,b,k,l,c,d,Dijab,in1)
!$OMP DO
    do i=MOZero,nmo
      do j=MOZero,nmo
        do a=MOZero,nmo
          do b=MOZero,nmo
            occ_factor = (SO_Occ(i)*(1.0d0-SO_Occ(a))*SO_Occ(j)*(1.0d0-SO_Occ(b)))

            Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))
            if(Dijab == 0.0d0)then
              cycle
            endif

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j)*occ_factor

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do c = MOZero,nmo
              do d = MOZero,c-1
                occ_factor = ((1.0d0-SO_Occ(a))*(1.0d0-SO_Occ(b))*(1.0d0-SO_Occ(c))*(1.0d0-SO_Occ(d)))
                in1 = in1 + 1.0d0 * t2old(c,d,i,j) * AMO(a,b,c,d)*occ_factor
              enddo
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do k = MOZero,nmo
              do l = MOZero,k-1
                occ_factor = (SO_Occ(k)*(SO_Occ(l))*SO_Occ(i)*(SO_Occ(j)))
                in1 = in1 + 1.0d0 * t2old(a,b,k,l) * AMO(k,l,i,j)*occ_factor
              enddo
            enddo

!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do c = MOZero,nmo
              do k = MOZero,nmo
                occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(b))*SO_Occ(j)*(1.0d0-SO_Occ(c)))
                in1 = in1 - 1.0d0 * t2old(a,c,i,k) * AMO(k,b,j,c)*occ_factor
                occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(a))*SO_Occ(j)*(1.0d0-SO_Occ(c)))
                in1 = in1 + 1.0d0 * t2old(b,c,i,k) * AMO(k,a,j,c)*occ_factor
                occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(b))*SO_Occ(i)*(1.0d0-SO_Occ(c)))
                in1 = in1 + 1.0d0 * t2old(a,c,j,k) * AMO(k,b,i,c)*occ_factor
                occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(a))*SO_Occ(i)*(1.0d0-SO_Occ(c)))
                in1 = in1 - 1.0d0 * t2old(b,c,j,k) * AMO(k,a,i,c)*occ_factor
              enddo
            enddo

            if(CCD)then ! do full CCD
            !Terms that are quadratic in T2
              do c = MOZero,nmo
                do d = MOZero,nmo
                  do k = MOZero,nmo
                    do l = MOZero,nmo
                      occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(c))*SO_Occ(l)*(1.0d0-SO_Occ(d)))
                      in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d) *occ_factor
                      in1 = in1 - 0.5d0 *  t2old(b,c,i,k)*t2old(d,a,l,j)   * AMO(k,l,c,d) *occ_factor
                      in1 = in1 + 0.5d0 *  t2old(b,c,j,k)*t2old(d,a,l,i)   * AMO(k,l,c,d) *occ_factor
                      in1 = in1 - 0.5d0 *  t2old(a,c,j,k)*t2old(d,b,l,i)   * AMO(k,l,c,d) *occ_factor

                      in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d) *occ_factor

                      in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d) *occ_factor
                      in1 = in1 + 0.5d0 *  t2old(b,c,i,j)*t2old(a,d,k,l)   * AMO(k,l,c,d) *occ_factor

                      in1 = in1 - 0.5d0 *  t2old(a,b,i,k)*t2old(c,d,j,l)   * AMO(k,l,c,d) *occ_factor
                      in1 = in1 + 0.5d0 *  t2old(a,b,j,k)*t2old(c,d,i,l)   * AMO(k,l,c,d) *occ_factor

                    enddo
                  enddo
                enddo
              enddo
            endif

            t2(a,b,i,j) =  1.0d0*in1/Dijab

          enddo
        enddo
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    Eccd = 0.0d0

    do i=MOZero,nmo
      do j=MOZero,i-1
        do a=MOZero,nmo
          do b=MOZero,a-1
            !if((SO_Occ(a) == 1.0) .or. (SO_Occ(b) == 1.0))then
            !  cycle
            !elseif((SO_Occ(i) == 0.0) .or. (SO_Occ(j) == 0.0))then
            !  cycle
            !endif
            Eccd = Eccd + 1.0d0*t2(a,b,i,j)*AMO(i,j,a,b)
          enddo
        enddo
      enddo
    enddo

    deltaE = Eccd-Eold
    if(CCD)then
      write(*,*) "      CCD: ",iT2,Eccd,deltaE
    elseif(MBPT3)then
      write(*,*) "      MBPT3: ",iT2,Eccd
      exit
    else
      write(*,*) "      LCCD: ",iT2,Eccd,deltaE
    endif
    if(abs(deltaE)<Ecctol)then
      write(*,*) '    '
      write(*,*) '    CC iterations converged '
      write(*,*) '    '
      exit
    endif

  enddo ! t2 loop

  endif

end subroutine calc_Eccd

end module module_cc

