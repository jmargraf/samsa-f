module module_cc
  implicit none

! Energies           
  double precision                            :: Elccd    = 0.0d0
  double precision, allocatable               :: t2(:,:,:,:),t2old(:,:,:,:)
!  integer                                     :: In4

contains

!#############################################
!#         Calculate Linearized CC
!#############################################
subroutine calc_Elccd
  use module_data,      only : MOI,Spins,dim_1e,nOccA,nOccB,Eps,SMO,Occ
  use module_data,      only : DropMO,DoDrop,Fock,DoSingles,dim_2e
  use module_ints,      only : Index2e
  use module_wavefun,   only : Fock_to_MO,Fock_to_AO
  implicit none
  integer                   :: iSpin,MOZero,iT2
  integer                   :: i,j,a,b,ia,ja,jb,ib,iajb,ibja
  integer                   :: c,d
  integer                   :: k,l

  integer                   :: nOcc,dimMO
  double precision          :: occ_factor,Dijab,E_2,E_2OS,E_2SS
  double precision          :: DampCC = 0.5d0

  Elccd  = 0.0d0
  MOZero = 0

!  call Fock_to_MO()

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-LCCD correlation energy"
    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO
    else
      MOZero = 0
    endif
    nOcc =  nOccA-1

    allocate(t2(0:dim_1e-1,0:dim_1e-1,0:dim_1e-1,0:dim_1e-1),  &
          t2old(0:dim_1e-1,0:dim_1e-1,0:dim_1e-1,0:dim_1e-1))

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0

  do iT2=1,50 ! T2 loop
    Elccd  = 0.0d0
    t2old  = t2
    do i=MOZero,nOcc
      do a=nOcc+1,dim_1e-1
        do j=MOZero,nOcc
          do b=nOcc+1,dim_1e-1

            iajb = In4(i,a,j,b)

            t2(i,a,j,b) = MOI(iajb) 

!                write(*,*) "cd terms"
            do c=nOcc+1,dim_1e-1
              do d=nOcc+1,dim_1e-1
!                if(a == b)then
!                  if(c == a .or. d == a)cycle
!                else
!                  if(c == a .and. d == a)cycle
!                  if(c == b .and. d == b)cycle
!                endif

            !sum:                     0.5d0  *      (ac|bd)   *      t(icjd) +
                t2(i,a,j,b) = t2(i,a,j,b) + 0.5d0*MOI(In4(a,c,b,d)) * t2old(i,c,j,d)
            !sum:                     0.5d0 *       (bc|ad)   *      t(jcid) +
                t2(i,a,j,b) = t2(i,a,j,b) + 0.5d0*MOI(In4(b,c,a,d)) * t2old(j,c,i,d)
              enddo
            enddo

!           write(*,*) "kl terms"
            do k=MoZero,nOcc
              do l=MOZero,nOcc
!                if(i == j)then
!                  if((k == i .or. l == i) .or.    &
!                     (k == j .or. l == j)) cycle
!                 else
!                  if((k == i .and. l == i) .or.    &
!                     (k == j .and. l == j)) cycle
!                 endif

            !sum:                      0.5d0 *       (ki|lj)  *      t(kalb) +
                 t2(i,a,j,b) = t2(i,a,j,b) + 0.5d0*MOI(In4(k,i,l,j))*t2old(k,a,l,b)
            !sum:                      0.5d0 *       (kj|li)  *      t(kbla) +
                 t2(i,a,j,b) = t2(i,a,j,b) + 0.5d0*MOI(In4(k,j,l,i))*t2old(k,b,l,a)

              enddo
            enddo

!           write(*,*) "kc terms"
            do k=MoZero,nOcc
              do c=nOcc+1,dim_1e-1
!               if(a == b .and. c == b)cycle
!               if(i == j)then
!                 if(k == i .or. k == j)cycle
!               endif

            !sum:                     2.0d0 *        (kc|bj)  *       t(iakc) +
                t2(i,a,j,b) = t2(i,a,j,b) + 2.0d0*(MOI(In4(k,c,b,j))*t2old(i,a,k,c))
            !sum:                     2.0d0 *        (kc|ai)  *       t(ijkc) +
                t2(i,a,j,b) = t2(i,a,j,b) + 2.0d0*(MOI(In4(k,c,a,i))*t2old(j,b,k,c))


            !sum:                    -1.0d0 *        (kc|bj)  *       t(icka) +
                t2(i,a,j,b) = t2(i,a,j,b) - 1.0d0*(MOI(In4(k,c,b,j))*t2old(i,c,k,a))
            !sum:                    -1.0d0 *        (kc|ai)  *       t(jckb) +
                t2(i,a,j,b) = t2(i,a,j,b) - 1.0d0*(MOI(In4(k,c,a,i))*t2old(j,c,k,b))


            !sum:                    -1.0d0 *        (ki|bc)  *       t(kajc) +
                t2(i,a,j,b) = t2(i,a,j,b) - 1.0d0*(MOI(In4(k,i,b,c))*t2old(k,a,j,c))
            !sum:                    -1.0d0 *        (kj|ac)  *       t(kbic) +
                t2(i,a,j,b) = t2(i,a,j,b) - 1.0d0*(MOI(In4(k,j,a,c))*t2old(k,b,i,c))


            !sum:                    -1.0d0 *        (kj|bc)  *       t(iakc) +
                t2(i,a,j,b) = t2(i,a,j,b) - 1.0d0*(MOI(In4(k,j,b,c))*t2old(i,a,k,c))
            !sum:                    -1.0d0 *        (ki|ac)  *       t(jbkc) +
                t2(i,a,j,b) = t2(i,a,j,b) - 1.0d0*(MOI(In4(k,i,a,c))*t2old(j,b,k,c))

              enddo
            enddo

          enddo
        enddo
      enddo
    enddo

    Elccd = 0.0d0
    E_2 = 0.0d0
    E_2OS = 0.0d0
    E_2SS = 0.0d0

    do i=MOZero,nOccA-1
      do a=nOccA,dim_1e-1
        call Index2e(i,a,ia)
        do j=MOZero,nOccA-1
          call Index2E(j,a,ja)
          do b=nOccA,dim_1e-1
            call Index2e(j,b,jb)
            call Index2e(i,b,ib)
            call Index2e(ia,jb,iajb)
            call Index2e(ib,ja,ibja)

            Dijab = Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1)
            Elccd = Elccd - t2(i,a,j,b)*MOI(iajb)/Dijab
            Elccd = Elccd - t2(i,a,j,b)*(MOI(iajb)-MOI(ibja))/Dijab

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

    dimMO = dim_1e*spins
    allocate(t2(1:dimMO,1:dimMO,1:dimMO,1:dimMO),t2old(1:dimMO,1:dimMO,1:dimMO,1:dimMO))

    if(DoDrop)then
      MOZero = DropMO*2 + 1
    else
      MOZero = 1
    endif

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0

    do iT2=1,3 ! T2 loop
      Elccd  = 0.0d0
      t2old  = t2

    ! T2     
    do i=MOZero,nOccA*2-1,2
      do a=nOccA*2+1,dim_1e*2-1,2
        do j=MOZero,nOccB*2,2
          do b=nOccB*2+1,dim_1e*2-1,2

            ! aa
              t2(i,a,j,b)         = t2(i,a,j,b) + (SMO(i,a,j,b) - SMO(i,b,j,a))
            ! bb
              t2(i+1,a+1,j+1,b+1) = t2(i+1,a+1,j+1,b+1) + & 
                                  (SMO(i+1,a+1,j+1,b+1) - &
                                   SMO(i+1,b+1,j+1,a+1)) 
            ! ab
              t2(i,a,j+1,b+1)     = t2(i,a,j+1,b+1) + SMO(i,a,j+1,b+1)

            enddo
          enddo
        enddo
      enddo

      Elccd = 0.0d0

    ! Energy 
    do i=MOZero,nOccA*2-1,2
      do a=nOccA*2+1,dim_1e*2-1,2
        do j=MOZero,nOccB*2,2
          do b=nOccB*2+1,dim_1e*2-1,2

            ! aa
              Dijab = Eps(a,1)-Eps(b,1)-Eps(i,1)-Eps(j,1)
              Elccd = Elccd - 0.5d0*(SMO(i,a,j,b) - SMO(i,b,j,a)) * t2(i,a,j,b)/Dijab
            ! bb
              Dijab = Eps(a,2)-Eps(b,2)-Eps(i,2)-Eps(j,2)
              Elccd = Elccd - 0.5d0*(SMO(i+1,a+1,j+1,b+1) - SMO(i+1,b+1,j+1,a+1)) * t2(i+1,a+1,j+1,b+1)/Dijab
            ! ab
              Dijab = Eps(a,1)-Eps(b,2)-Eps(i,1)-Eps(j,2)
              Elccd = Elccd -  SMO(i,a,j+1,b+1) * t2(i,a,j+1,b+1)/Dijab

            enddo
          enddo
        enddo
      enddo

      write(*,*) "      LCCD: ",iT2,Elccd

    enddo ! End T2 loop

  endif

!  call Fock_to_AO()

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
