module module_cc
  implicit none

! Energies           
  double precision                            :: Eccd    = 0.0d0
  double precision                            :: Ecc1    = 0.0d0
  double precision                            :: Ecc2    = 0.0d0
  double precision, allocatable               :: t2(:,:,:,:),t2old(:,:,:,:)
  double precision, allocatable               :: X1(:,:,:,:),X2(:,:,:,:),X3(:,:),X4(:,:)
  double precision, allocatable               :: SO_Occ(:)

contains

!#############################################
!#         Calculate CC doubles
!#############################################
subroutine calc_Eccd(CCD,frac_occ)
  use module_data,      only : Spins,dim_1e,nOccA,nOccB,Eps,Occ
  use module_data,      only : DropMO,DoDrop,DoSingles,dim_1e,AMO
  use module_data,      only : DampCC,DoCCLS,CCMax
  use module_data,      only : Eps_SO,F_SO,DoReNorm,OmegaReg
  use module_energy,    only : Etot
  implicit none
  logical                   :: CCD 
  logical                   :: frac_occ
  logical                   :: MBPT3 = .false.
  logical                   :: factorized = .true.
  integer                   :: iSpin,MOZero,iT2,nfocca,nfoccb,nfocc
  integer                   :: i,j,a,b 
  integer                   :: k,l,c,d

  integer                   :: nOcc,nmo
  integer                   :: maxT = 1000
  double precision          :: occ_factor,Dijab
  double precision          :: Ecctol = 1.0d-6
  double precision          :: LSmult = 1.0d4
  double precision          :: Eold   = 0.0d0
  double precision          :: deltaE = 99.0d0
  double precision          :: in1 ! intermediate during T Amp equations
  double precision          :: Xin1,Xin2,Xin3,Xin4

  maxT = CCMax

  if (.not. DoCCLS) then
      LSmult=1.0d0
  else
      LSmult=1.0d4
  endif

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

    !construct intermediates here
    if(factorized)then
       if(.not.allocated(X1))then
         allocate(X1(MOZero:nOcc,nOcc+1:nmo,nOcc+1:nmo,MOZero:nOcc), &
                  X2(MOZero:nOcc,MOZero:nOcc,MOZero:nOcc,MOZero:nOcc), &
                  X3(nOcc+1:nmo,nOcc+1:nmo), &  
                  X4(MOZero:nOcc,MOZero:nOcc))
       endif       

       do l=MOZero,nOcc
         do d=nOcc+1,nmo 
           do a=nOcc+1,nmo
             do i=MOZero,nOcc
               Xin1 = 0.0d0
               do k=MOZero,nOcc
                 do c=nOcc+1,nmo 
                   !in1 = in1 + 0.5d0 *  t2old(d,b,l,j) * X1(l,d,a,i)
                   !in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                   Xin1 = Xin1 + t2old(a,c,i,k)*AMO(k,l,c,d)
                 enddo
               enddo
               X1(l,d,a,i) = Xin1
             enddo
           enddo
         enddo
       enddo

       do k=MOZero,nOcc
         do l=MOZero,nOcc
           do i=MOZero,nOcc
             do j=MOZero,nOcc
               Xin2 = 0.0d0
               do c=nOcc+1,nmo 
                 do d=nOcc+1,nmo 
                   !in1 = in1 + 0.25d0*  t2old(a,b,k,l) * X2(k,l,i,j)
                   !in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)
                   Xin2 = Xin2 + t2old(c,d,i,j)*AMO(k,l,c,d)
                 enddo
               enddo
               X2(k,l,i,j) = Xin2
             enddo
           enddo
         enddo
       enddo

       do b=nOcc+1,nmo
         do c=nOcc+1,nmo
           Xin3 = 0.0d0
           do k=MOZero,nOcc
             do l=MOZero,nOcc
               do d=nOcc+1,nmo
                 !in1 = in1 - 0.5d0 *  t2old(a,c,i,j) * X3(b,c)
                 !in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                 Xin3 = Xin3 + t2old(b,d,k,l)*AMO(k,l,c,d)
               enddo
             enddo
           enddo
           X3(b,c) = Xin3
         enddo
       enddo


       do k=MOZero,nOcc 
         do j=MOZero,nOcc
           Xin4 = 0.0d0
           do l=MOZero,nOcc
             do c=nOcc+1,nmo 
               do d=nOcc+1,nmo
                 Xin4 = Xin4 + t2old(c,d,j,l)*AMO(k,l,c,d) 
               enddo
             enddo
           enddo
           X4(k,j) = Xin4
         enddo
       enddo

    endif

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

            if(CCD .and. .not. factorized)then ! do full CCD
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
            elseif(CCD .and. factorized)then 
               ! Sum over l,d
                do d = nOcc+1,nmo
                  do l = MOZero,nOcc
                      in1 = in1 + 0.5d0 *  t2old(d,b,l,j) * X1(l,d,a,i)
                     !in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(d,a,l,j) * X1(l,d,b,i)
                     !in1 = in1 - 0.5d0 *  t2old(b,c,i,k)*t2old(d,a,l,j)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(d,a,l,i) * X1(l,d,b,j) 
                     !in1 = in1 + 0.5d0 *  t2old(b,c,j,k)*t2old(d,a,l,i)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(d,b,l,i) * X1(l,d,a,j)
                     !in1 = in1 - 0.5d0 *  t2old(a,c,j,k)*t2old(d,b,l,i)   * AMO(k,l,c,d)
                  enddo
                enddo

               ! Sum over k,l
                do k = MOZero,nOcc
                  do l = MOZero,nOcc
                      in1 = in1 + 0.25d0*  t2old(a,b,k,l) * X2(k,l,i,j)
                     !in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)
                  enddo
                enddo

               ! Sum over c
                do c = nOcc+1,nmo
                      in1 = in1 - 0.5d0 *  t2old(a,c,i,j) * X3(b,c) 
                     !in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,i,j) * X3(a,c)
                     !in1 = in1 + 0.5d0 *  t2old(b,c,i,j)*t2old(a,d,k,l)   * AMO(k,l,c,d)
                enddo

               ! Sum over k
                do k = MOZero,nOcc
                      in1 = in1 - 0.5d0 *  t2old(a,b,i,k) * X4(k,j)
                     !in1 = in1 - 0.5d0 *  t2old(a,b,i,k)*t2old(c,d,j,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(a,b,j,k) * X4(k,i)
                     !in1 = in1 + 0.5d0 *  t2old(a,b,j,k)*t2old(c,d,i,l)   * AMO(k,l,c,d)
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
      write(*,*) '    Total Energy = ', Eccd+Etot, '(',Etot,'+',Eccd,')'
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

    !construct intermediates here
    if(factorized)then
       if(.not.allocated(X1))then
         allocate(X1(MOZero:nmo,MOZero:nmo,MOZero:nmo,MOZero:nmo), &
                  X2(MOZero:nmo,MOZero:nmo,MOZero:nmo,MOZero:nmo), &
                  X3(MOZero:nmo,MOZero:nmo), &
                  X4(MOZero:nmo,MOZero:nmo))
       endif

       do l=MOZero,nmo 
         if (SO_Occ(l) == 0.0) cycle
         do d=MOZero,nmo 
           if (SO_Occ(d) == 1.0) cycle
           do a=MOZero,nmo 
             if (SO_Occ(a) == 1.0) cycle
             do i=MOZero,nmo 
               if (SO_Occ(i) == 0.0) cycle
               Xin1 = 0.0d0
               do k=MOZero,nmo 
                 if (SO_Occ(k) == 0.0) cycle
                 do c=MOZero,nmo 
                   if (SO_Occ(c) == 1.0) cycle
                   !in1 = in1 + 0.5d0 *  t2old(d,b,l,j) * X1(l,d,a,i)
                   !in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                   Xin1 = Xin1 + t2old(a,c,i,k)*AMO(k,l,c,d)
                 enddo
               enddo
               X1(l,d,a,i) = Xin1
             enddo
           enddo
         enddo
       enddo

       do k=MOZero,nmo 
         if (SO_Occ(k) == 0.0) cycle
         do l=MOZero,nmo 
          if (SO_Occ(l) == 0.0) cycle
           do i=MOZero,nmo 
             if (SO_Occ(i) == 0.0) cycle
             do j=MOZero,nmo 
               if (SO_Occ(j) == 0.0) cycle 
               Xin2 = 0.0d0
               do c=MOZero,nmo 
                 if (SO_Occ(c) == 1.0) cycle
                 do d=MOZero,nmo 
                   if (SO_Occ(d) == 1.0) cycle
                   !in1 = in1 + 0.25d0*  t2old(a,b,k,l) * X2(k,l,i,j)
                   !in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)
                   Xin2 = Xin2 + t2old(c,d,i,j)*AMO(k,l,c,d)
                 enddo
               enddo
               X2(k,l,i,j) = Xin2
             enddo
           enddo
         enddo
       enddo

       do b=MOZero,nmo 
         if (SO_Occ(b) == 1.0) cycle
         do c=MOZero,nmo 
           if (SO_Occ(c) == 1.0) cycle
           Xin3 = 0.0d0
           do k=MOZero,nmo 
             if (SO_Occ(k) == 0.0) cycle
             do l=MOZero,nmo 
               if (SO_Occ(l) == 0.0) cycle
               do d=MOZero,nmo 
                 if (SO_Occ(d) == 1.0) cycle
                 !in1 = in1 - 0.5d0 *  t2old(a,c,i,j) * X3(b,c)
                 !in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                 Xin3 = Xin3 + t2old(b,d,k,l)*AMO(k,l,c,d)
               enddo
             enddo
           enddo
           X3(b,c) = Xin3
         enddo
       enddo


       do k=MOZero,nmo 
         if (SO_Occ(k) == 0.0) cycle
         do j=MOZero,nmo 
           if (SO_Occ(j) == 0.0) cycle
           Xin4 = 0.0d0
           do l=MOZero,nmo 
             if (SO_Occ(l) == 0.0) cycle
             do c=MOZero,nmo 
               if (SO_Occ(c) == 1.0) cycle
               do d=MOZero,nmo 
                 if (SO_Occ(d) == 1.0) cycle
                 Xin4 = Xin4 + t2old(c,d,j,l)*AMO(k,l,c,d)
               enddo
             enddo
           enddo
           X4(k,j) = Xin4
         enddo
       enddo

    endif

    !T2 equations
    !Nocc version
!$OMP PARALLEL PRIVATE(i,j,a,b,k,l,c,d,Dijab,in1)
!$OMP DO
    do i=MOZero,nmo
      if (SO_Occ(i) == 0.0) cycle
      do j=MOZero,nmo
        if (SO_Occ(j) == 0.0) cycle
        do a=MOZero,nmo
          if (SO_Occ(a) == 1.0) cycle
          do b=MOZero,nmo
            if (SO_Occ(b) == 1.0) cycle
            !if((SO_Occ(a) == 1.0) .or. (SO_Occ(b) == 1.0))then
            !  cycle
            !elseif((SO_Occ(i) == 0.0) .or. (SO_Occ(j) == 0.0))then
            !  cycle
            !endif

            Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j)

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do c = MOZero,nmo
              if (SO_Occ(c) == 1.0) cycle
              do d = MOZero,c-1
                if (SO_Occ(d) == 1.0) cycle
                !if((SO_Occ(c) == 1.0) .or. (SO_Occ(d) == 1.0))then
                !  cycle
                !endif
                in1 = in1 + 1.0d0 * t2old(c,d,i,j) * AMO(a,b,c,d)
              enddo
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do k = MOZero,nmo
              if (SO_Occ(k) == 0.0) cycle
              do l = MOZero,k-1
                if (SO_Occ(l) == 0.0) cycle
                !if((SO_Occ(k) == 0.0) .or. (SO_Occ(l) == 0.0))then
                !  cycle
                !endif
                in1 = in1 + 1.0d0 * t2old(a,b,k,l) * AMO(k,l,i,j)
              enddo
            enddo

!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do c = MOZero,nmo
              if (SO_Occ(c) == 1.0) cycle
              do k = MOZero,nmo
                if (SO_Occ(k) == 0.0) cycle
                !if((SO_Occ(c) == 1.0) .or. (SO_Occ(k) == 0.0))then
                !  cycle
                !endif
                in1 = in1 - 1.0d0 * t2old(a,c,i,k) * AMO(k,b,j,c)
                in1 = in1 + 1.0d0 * t2old(b,c,i,k) * AMO(k,a,j,c)
                in1 = in1 + 1.0d0 * t2old(a,c,j,k) * AMO(k,b,i,c)
                in1 = in1 - 1.0d0 * t2old(b,c,j,k) * AMO(k,a,i,c)
              enddo
            enddo

            !if(CCD)then ! do full CCD
            if(CCD .and. .not. factorized)then ! do full CCD
            !Terms that are quadratic in T2
              do c = MOZero,nmo
                if (SO_Occ(c) == 1.0) cycle
                do d = MOZero,nmo
                  if (SO_Occ(d) == 1.0) cycle
                  do k = MOZero,nmo
                    if (SO_Occ(k) == 0.0) cycle
                    do l = MOZero,nmo
                      if (SO_Occ(l) == 0.0) cycle
                      !if((SO_Occ(k) == 0.0) .or. (SO_Occ(l) == 0.0))then
                      !  cycle
                      !elseif((SO_Occ(c) == 1.0) .or. (SO_Occ(d) == 1.0))then
                      !  cycle
                      !endif

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
            elseif(CCD .and. factorized)then
               ! Sum over l,d
                do d = MOZero,nmo
                  if (SO_Occ(d) == 1.0) cycle
                  do l = MOZero,nmo 
                      if (SO_Occ(l) == 0.0) cycle
                     !if(SO_Occ(l) == 0.0)then
                      !  cycle
                      !elseif(SO_Occ(d) == 1.0)then
                      !  cycle
                      !endif
                      in1 = in1 + 0.5d0 *  t2old(d,b,l,j) * X1(l,d,a,i)
                     !in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(d,a,l,j) * X1(l,d,b,i)
                     !in1 = in1 - 0.5d0 *  t2old(b,c,i,k)*t2old(d,a,l,j)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(d,a,l,i) * X1(l,d,b,j)
                     !in1 = in1 + 0.5d0 *  t2old(b,c,j,k)*t2old(d,a,l,i)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(d,b,l,i) * X1(l,d,a,j)
                     !in1 = in1 - 0.5d0 *  t2old(a,c,j,k)*t2old(d,b,l,i)   * AMO(k,l,c,d)
                  enddo
                enddo

               ! Sum over k,l
                do k = MOZero,nmo
                  if (SO_Occ(k) == 0.0) cycle
                  do l = MOZero,nmo
                       if (SO_Occ(l) == 0.0) cycle
                      !if((SO_Occ(k) == 0.0) .or. (SO_Occ(l) == 0.0))then
                      !  cycle
                      !endif
                      in1 = in1 + 0.25d0*  t2old(a,b,k,l) * X2(k,l,i,j)
                     !in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)
                  enddo
                enddo

               ! Sum over c
                do c = MOZero,nmo
                      if (SO_Occ(c) == 1.0) cycle
                      in1 = in1 - 0.5d0 *  t2old(a,c,i,j) * X3(b,c)
                     !in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,i,j) * X3(a,c)
                     !in1 = in1 + 0.5d0 *  t2old(b,c,i,j)*t2old(a,d,k,l)   * AMO(k,l,c,d)
                enddo

               ! Sum over k
                do k = MOZero,nmo
                      if (SO_Occ(k) == 0.0) cycle
                      in1 = in1 - 0.5d0 *  t2old(a,b,i,k) * X4(k,j)
                     !in1 = in1 - 0.5d0 *  t2old(a,b,i,k)*t2old(c,d,j,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(a,b,j,k) * X4(k,i)
                     !in1 = in1 + 0.5d0 *  t2old(a,b,j,k)*t2old(c,d,i,l)   * AMO(k,l,c,d)
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
      write(*,*) '    Total Energy = ', Eccd+Etot, '(',Etot,'+',Eccd,')'
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

    nfocca = -1
    nfoccb = -1
    nfocc = 0

    do i=2,nmo,2
      SO_Occ(i-1) = Occ(i/2,1)
      SO_Occ(i) = Occ(i/2,2)
      if ((SO_Occ(i-1) == 0.0d0) .and. (nfocca == -1)) then
          nfocca = i-1
      endif
      if ((SO_Occ(i) == 0.0d0) .and. (nfoccb == -1)) then
          nfoccb = i
      endif
    enddo

    nfocc = max(nfocca,nfoccb)

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif
    nOcc =  nOccA*2

    if(.not.allocated(t2))then
      allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
          t2old(1:nmo,1:nmo,1:nmo,1:nmo))
          !t2 = 0.0d0
          !t2old = 0.0d0
    endif

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0
    Eccd = 0.0d0

  do iT2=1,maxT ! T amplitude loop
    Eold   = Eccd
    t2old  = (1.0d0-DampCC)*t2+DampCC*t2old ! uses mixing

    !construct intermediates here
    if(factorized)then
       if(.not.allocated(X1))then
         allocate(X1(MOZero:nmo,MOZero:nmo,MOZero:nmo,MOZero:nmo), &
                  X2(MOZero:nmo,MOZero:nmo,MOZero:nmo,MOZero:nmo), &
                  X3(MOZero:nmo,MOZero:nmo), &
                  X4(MOZero:nmo,MOZero:nmo))
       endif

       do l=MOZero,nmo
         if (SO_Occ(l) == 0.0d0) cycle
         do d=MOZero,nmo
           if (SO_Occ(d) == 1.0d0) cycle
           do a=MOZero,nmo
             if (SO_Occ(a) == 1.0d0) cycle
             do i=MOZero,nmo
               if (SO_Occ(i) == 0.0d0) cycle
               Xin1 = 0.0d0
               do k=MOZero,nmo
                 if (SO_Occ(k) == 0.0d0) cycle
                 do c=MOZero,nmo
                   if (SO_Occ(c) == 1.0d0) cycle
                   !in1 = in1 + 0.5d0 *  t2old(d,b,l,j) * X1(l,d,a,i)
                   !in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                   occ_factor = (1.0d0-SO_Occ(c))*(SO_Occ(k))      
                   Xin1 = Xin1 + t2old(a,c,i,k)*AMO(k,l,c,d)*occ_factor
                 enddo
               enddo
               X1(l,d,a,i) = Xin1
             enddo
           enddo
         enddo
       enddo

       do k=MOZero,nmo
         if (SO_Occ(k) == 0.0d0) cycle
         do l=MOZero,nmo
          if (SO_Occ(l) == 0.0d0) cycle
           do i=MOZero,nmo
             if (SO_Occ(i) == 0.0d0) cycle
             do j=MOZero,nmo
               if (SO_Occ(j) == 0.0d0) cycle
               Xin2 = 0.0d0
               do c=MOZero,nmo
                 if (SO_Occ(c) == 1.0d0) cycle
                 do d=MOZero,nmo
                   if (SO_Occ(d) == 1.0d0) cycle
                   !in1 = in1 + 0.25d0*  t2old(a,b,k,l) * X2(k,l,i,j)
                   !in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)
                   occ_factor = (1.0d0-SO_Occ(c))*(1.0d0-SO_Occ(d))       
                   Xin2 = Xin2 + t2old(c,d,i,j)*AMO(k,l,c,d)*occ_factor
                 enddo
               enddo
               X2(k,l,i,j) = Xin2
             enddo
           enddo
         enddo
       enddo

       do b=MOZero,nmo
         if (SO_Occ(b) == 1.0) cycle
         do c=MOZero,nmo
           if (SO_Occ(c) == 1.0) cycle
           Xin3 = 0.0d0
           do k=MOZero,nmo
             if (SO_Occ(k) == 0.0) cycle
             do l=MOZero,nmo
               if (SO_Occ(l) == 0.0) cycle
               do d=MOZero,nmo
                 if (SO_Occ(d) == 1.0) cycle
                 !in1 = in1 - 0.5d0 *  t2old(a,c,i,j) * X3(b,c)
                 !in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                 occ_factor = (1.0d0-SO_Occ(d))*(SO_Occ(k))*(SO_Occ(l))       
                 Xin3 = Xin3 + t2old(b,d,k,l)*AMO(k,l,c,d)*occ_factor
               enddo
             enddo
           enddo
           X3(b,c) = Xin3
         enddo
       enddo


       do k=MOZero,nmo
         if (SO_Occ(k) == 0.0) cycle
         do j=MOZero,nmo
           if (SO_Occ(j) == 0.0) cycle
           Xin4 = 0.0d0
           do l=MOZero,nmo
             if (SO_Occ(l) == 0.0) cycle
             do c=MOZero,nmo
               if (SO_Occ(c) == 1.0) cycle
               do d=MOZero,nmo
                 if (SO_Occ(d) == 1.0) cycle
                 occ_factor = (1.0d0-SO_Occ(d))*(1.0d0-SO_Occ(c))*(SO_Occ(l))
                 Xin4 = Xin4 + t2old(c,d,j,l)*AMO(k,l,c,d)*occ_factor
               enddo
             enddo
           enddo
           X4(k,j) = Xin4
         enddo
       enddo
       !write(*,*) 'factorization done'
    endif

    !T2 equations
    !Nocc version
!$OMP PARALLEL PRIVATE(i,j,a,b,k,l,c,d,Dijab,in1)
!$OMP DO
    do i=MOZero,nfocc
      if (SO_Occ(i) == 0.0) cycle
      do j=MOZero,nfocc
        if (SO_Occ(j) == 0.0) cycle
        do a=MOZero,nmo
          if (SO_Occ(a) == 1.0) cycle
          do b=MOZero,nmo
            if (SO_Occ(b) == 1.0) cycle
            !occ_factor = (SO_Occ(i)*(1.0d0-SO_Occ(a))*SO_Occ(j)*(1.0d0-SO_Occ(b)))

            if(DoRenorm)then
              Dijab = (SO_Occ(i)*Eps_SO(i)        + SO_Occ(j)*Eps_SO(j) & 
                    - (1.0d0-SO_Occ(a))*Eps_SO(a) - (1.0d0-SO_Occ(b))*Eps_SO(b))
            else
              Dijab = (Eps_SO(i) + Eps_SO(j) - Eps_SO(a) - Eps_SO(b))
            endif

            Dijab = (Dijab*Dijab+LSmult*OmegaReg)/(Dijab)
            !if(Dijab == 0.0d0)then
            !  cycle
            !endif

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(a,b,i,j) !*occ_factor

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do c = MOZero,nmo
              if (SO_Occ(c) == 1.0) cycle
              do d = MOZero,c-1
                if (SO_Occ(d) == 1.0) cycle
                occ_factor = (1.0d0-SO_Occ(c))*(1.0d0-SO_Occ(d))   !((1.0d0-SO_Occ(a))*(1.0d0-SO_Occ(b))*(1.0d0-SO_Occ(c))*(1.0d0-SO_Occ(d)))
                in1 = in1 + 1.0d0 * t2old(c,d,i,j) * AMO(a,b,c,d)*occ_factor
              enddo
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do k = MOZero,nmo
              if (SO_Occ(k) == 0.0) cycle
              do l = MOZero,k-1
                if (SO_Occ(l) == 0.0) cycle
                occ_factor = (SO_Occ(k))*(SO_Occ(l)) !*SO_Occ(i)*(SO_Occ(j)))
                in1 = in1 + 1.0d0 * t2old(a,b,k,l) * AMO(k,l,i,j)*occ_factor
              enddo
            enddo

!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do c = MOZero,nmo
              if (SO_Occ(c) == 1.0) cycle
              do k = MOZero,nmo
                if (SO_Occ(k) == 0.0) cycle
                occ_factor = (SO_Occ(k))*(1.0d0-SO_Occ(c))

                !occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(b))*SO_Occ(j)*(1.0d0-SO_Occ(c)))
                in1 = in1 - 1.0d0 * t2old(a,c,i,k) * AMO(k,b,j,c)*occ_factor
                !occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(a))*SO_Occ(j)*(1.0d0-SO_Occ(c)))
                in1 = in1 + 1.0d0 * t2old(b,c,i,k) * AMO(k,a,j,c)*occ_factor
                !occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(b))*SO_Occ(i)*(1.0d0-SO_Occ(c)))
                in1 = in1 + 1.0d0 * t2old(a,c,j,k) * AMO(k,b,i,c)*occ_factor
                !occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(a))*SO_Occ(i)*(1.0d0-SO_Occ(c)))
                in1 = in1 - 1.0d0 * t2old(b,c,j,k) * AMO(k,a,i,c)*occ_factor
              enddo
            enddo

            if(CCD .and. .not. factorized)then ! do full CCD
            !Terms that are quadratic in T2
              do c = MOZero,nmo
                if (SO_Occ(c) == 1.0) cycle
                do d = MOZero,nmo
                  if (SO_Occ(d) == 1.0) cycle
                  do k = MOZero,nmo
                    if (SO_Occ(k) == 0.0) cycle
                    do l = MOZero,nmo
                      if (SO_Occ(l) == 0.0) cycle
      
                      occ_factor = (SO_Occ(k)*(1.0d0-SO_Occ(c))*SO_Occ(l)*(1.0d0-SO_Occ(d)))
                      !write(*,*) occ_factor
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
            elseif(CCD .and. factorized)then
               ! Sum over l,d
                do d = MOZero,nmo
                  if (SO_Occ(d) == 1.0) cycle
                  do l = MOZero,nmo
                      if (SO_Occ(l) == 0.0) cycle
                      occ_factor = (1.0d0-SO_Occ(d))*(SO_Occ(l))
                      in1 = in1 + 0.5d0 *  t2old(d,b,l,j) * X1(l,d,a,i) *occ_factor
                     !in1 = in1 + 0.5d0 *  t2old(a,c,i,k)*t2old(d,b,l,j)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(d,a,l,j) * X1(l,d,b,i) *occ_factor
                     !in1 = in1 - 0.5d0 *  t2old(b,c,i,k)*t2old(d,a,l,j)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(d,a,l,i) * X1(l,d,b,j) *occ_factor
                     !in1 = in1 + 0.5d0 *  t2old(b,c,j,k)*t2old(d,a,l,i)   * AMO(k,l,c,d)
                      in1 = in1 - 0.5d0 *  t2old(d,b,l,i) * X1(l,d,a,j) *occ_factor
                     !in1 = in1 - 0.5d0 *  t2old(a,c,j,k)*t2old(d,b,l,i)   * AMO(k,l,c,d)
                  enddo
                enddo

               ! Sum over k,l
                do k = MOZero,nmo
                  if (SO_Occ(k) == 0.0) cycle
                  do l = MOZero,nmo
                      if (SO_Occ(l) == 0.0) cycle
                      occ_factor = (1.0d0-SO_Occ(c))*(SO_Occ(k))
                      in1 = in1 + 0.25d0*  t2old(a,b,k,l) * X2(k,l,i,j) *occ_factor
                     !in1 = in1 + 0.25d0*  t2old(c,d,i,j)*t2old(a,b,k,l)   * AMO(k,l,c,d)
                  enddo
                enddo

               ! Sum over c
                do c = MOZero,nmo
                      if (SO_Occ(c) == 1.0) cycle
                      occ_factor = (1.0d0-SO_Occ(c))
                      in1 = in1 - 0.5d0 *  t2old(a,c,i,j) * X3(b,c) *occ_factor
                     !in1 = in1 - 0.5d0 *  t2old(a,c,i,j)*t2old(b,d,k,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(b,c,i,j) * X3(a,c) *occ_factor
                     !in1 = in1 + 0.5d0 *  t2old(b,c,i,j)*t2old(a,d,k,l)   * AMO(k,l,c,d)
                enddo

               ! Sum over k
                do k = MOZero,nmo
                      if (SO_Occ(k) == 0.0) cycle
                      occ_factor = (SO_Occ(k))
                      in1 = in1 - 0.5d0 *  t2old(a,b,i,k) * X4(k,j) *occ_factor
                     !in1 = in1 - 0.5d0 *  t2old(a,b,i,k)*t2old(c,d,j,l)   * AMO(k,l,c,d)
                      in1 = in1 + 0.5d0 *  t2old(a,b,j,k) * X4(k,i) *occ_factor
                     !in1 = in1 + 0.5d0 *  t2old(a,b,j,k)*t2old(c,d,i,l)   * AMO(k,l,c,d)
                enddo
            endif
 
            in1 = 1.0d0*in1/Dijab

            if(abs(in1)<=1.0d0)then
                t2(a,b,i,j) =  in1
            elseif(in1 < -1.0d0)then
                !t2(a,b,i,j) =  -1.0d0
                write(*,*) '! warning: tamp < -1.0', in1
            elseif(in1 > 1.0d0)then
                !t2(a,b,i,j) =  1.0d0
                write(*,*) '! warning: tamp > 1.0', in1
            endif

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
            occ_factor = (SO_Occ(i)*(1.0d0-SO_Occ(a))*SO_Occ(j)*(1.0d0-SO_Occ(b)))
            !if((SO_Occ(a) == 1.0) .or. (SO_Occ(b) == 1.0))then
            !  cycle
            !elseif((SO_Occ(i) == 0.0) .or. (SO_Occ(j) == 0.0))then
            !  cycle
            !endif
            Eccd = Eccd + 1.0d0*t2(a,b,i,j)*AMO(i,j,a,b)*occ_factor
          enddo
        enddo
      enddo
    enddo

    if(iT2 == 1)then
      Ecc1 = Eccd
    elseif(iT2 == 2)then
      Ecc2 = Eccd
    endif

    deltaE = Eccd-Eold
    if(CCD)then
      write(*,*) "      CCD: ",iT2,Eccd,deltaE,DampCC,LSmult*OmegaReg
    elseif(MBPT3)then
      write(*,*) "      MBPT3: ",iT2,Eccd
      exit
    else
      write(*,*) "      LCCD: ",iT2,Eccd,deltaE,DampCC,LSmult*OmegaReg
    endif
    if(abs(deltaE)<Ecctol*1.0d0 .and. DampCC==0.6d0 .and. LSmult == 1.0d0)then
      DampCC=0.0d0
      write(*,*) '    Turning off damping '
    elseif(abs(deltaE)<Ecctol*1.0d0 .and. DampCC>0.8d0 .and. LSmult == 1.0d0)then
      DampCC=0.8d0
      write(*,*) '    set damping = 0.8 '
    elseif(abs(deltaE)<Ecctol*1.0d0 .and. DampCC>0.6d0 .and. LSmult == 1.0d0)then
      DampCC=0.6d0
      write(*,*) '    set damping = 0.6 '
    elseif(abs(deltaE)<Ecctol .and. LSmult == 1.0d0 .and. DampCC==0.0)then
      write(*,*) '    '
      write(*,*) '    CC iterations converged '
      write(*,*) '    '
      write(*,*) '    Total Energy = ', Eccd+Etot, '(',Etot,'+',Eccd,')' 
      exit
    endif

    if(DoCCLS .and. abs(deltaE)<Ecctol*0.1d0 .and. LSmult <= 10.0d0)then
      LSmult=1.0d0
      write(*,*) '    Turning off level shift '
    elseif(DoCCLS .and. abs(deltaE)<0.1d0*Ecctol .and. LSmult > 1.0d0)then
      LSmult=0.5d0*LSmult
      write(*,*) '    decrease level shift to',LSmult*OmegaReg
    endif

  enddo ! t2 loop

  endif

end subroutine calc_Eccd

end module module_cc

