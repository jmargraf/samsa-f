module module_cc
  implicit none

! Energies           
  double precision                            :: Elccd    = 0.0d0
  double precision, allocatable               :: t2(:),t2old(:)

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
  integer                   :: c,d,ac,bd,ad,bc,acbd,adbc,ic,jd,icjd
  integer                   :: k,l,ki,lj,kj,li,kilj,kjli,ka,lb,kalb
  integer                   :: kc,bj,kjbc,kcbj,iakc,bjkc
  integer                   :: aj,kjac,kcaj,jc,ibjc,id
  integer                   :: bi,kibc,kcbi,jakc,bcad,jcid
  integer                   :: ai,kiac,kcai,jbkc,la,kbla
  integer                   :: icka,kb,jckb,kajc,kbjc,aikc
  integer                   :: kbic

  integer                   :: nOcc
  double precision          :: occ_factor,Dijab,E_2,E_2OS,E_2SS
  double precision          :: DampCC = 0.5d0

  Elccd  = 0.0d0
  MOZero = 0

!  call Fock_to_MO()

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-LCCD correlation energy"
    write(*,*) "    "
            !t2(iajb) = (ia|jb) +
            !sum:       0.5d0*(ac|bd)*t(icjd) +
            !sum:       0.5d0*(bc|ad)*t(jcid) +

            !sum:       0.5d0*(ki|lj)*t(kalb) +
            !sum:       0.5d0*(kj|li)*t(kbla) +

            !sum:       2.0d0*(kc|bj)*t(iakc) +
            !sum:       2.0d0*(kc|ai)*t(jbkc) +

            !sum:      -1.0d0*(kc|bj)*t(icka) +
            !sum:      -1.0d0*(kc|ai)*t(jckb) +

            !sum:      -1.0d0*(ki|bc)*t(kajc) +
            !sum:      -1.0d0*(kj|ac)*t(kbic) +

            !sum:      -1.0d0*(kj|bc)*t(iakc) +
            !sum:      -1.0d0*(ki|ac)*t(jbkc) +


    if(DoDrop)then
      MOZero = DropMO
    else
      MOZero = 0
    endif
    nOcc =  nOccA-1

    allocate(t2(0:dim_2e-1),t2old(0:dim_2e-1))

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0

  do iT2=1,150 ! T2 loop
    Elccd  = 0.0d0
    t2old = t2
    do i=MOZero,nOcc
      do a=nOcc+1,dim_1e-1
        call Index2e(i,a,ia)
        do j=MOZero,nOcc
          call Index2E(j,a,ja)
          do b=nOcc+1,dim_1e-1
            call Index2e(j,b,jb)
            call Index2e(i,b,ib)
            call Index2e(ia,jb,iajb)
            call Index2e(ib,ja,ibja)

            t2(iajb) = MOI(iajb) 

            !sum:       0.5d0*(ac|bd)*t(icjd) +
            !sum:       0.5d0*(ic|jd)*t(acbd) +
!                write(*,*) "cd terms"
            do c=nOcc+1,dim_1e-1
              do d=nOcc+1,dim_1e-1
!                if(a == b)then
!                  if(c == a .or. d == a)cycle
!                else
!                  if(c == a .and. d == a)cycle
!                  if(c == b .and. d == b)cycle
!                endif
                call Index2e(a,c,ac)
                call Index2e(b,d,bd)
                call Index2e(ac,bd,acbd)
                call Index2e(i,c,ic)
                call Index2e(j,d,jd)
                call Index2e(ic,jd,icjd)
                call Index2e(b,c,bc)
                call Index2e(a,d,ad)
                call Index2e(bc,ad,bcad)
                call Index2e(j,c,jc)
                call Index2e(i,d,id)
                call Index2e(jc,id,jcid)

            !sum:                     0.5d0 *(ac|bd) * t(icjd) +
                t2(iajb) = t2(iajb) + 0.5d0*MOI(acbd)*t2old(icjd)
            !sum:                     0.5d0 * (bc|ad) * t(jcid) +
                t2(iajb) = t2(iajb) + 0.5d0*MOI(bcad)*t2old(jcid)
              enddo
            enddo

            !sum:       0.5d0*(ki|lj)*t(kalb) +
            !sum:       0.5d0*(kj|li)*t(kbla) +
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
                call Index2e(k,i,ki)
                call Index2e(l,j,lj)
                call Index2e(ki,lj,kilj)
                call Index2e(k,j,kj)
                call Index2e(l,i,li)
                call Index2e(kj,li,kjli)
                call Index2e(k,a,ka)
                call Index2e(l,b,lb)
                call Index2e(ka,lb,kalb)
                call Index2e(k,b,kb)
                call Index2e(l,a,la)
                call Index2e(kb,la,kbla)

            !sum:                     0.5d0 * (ki|lj) * t(kalb) +
                t2(iajb) = t2(iajb) + 0.5d0*MOI(kilj)*t2old(kalb)
            !sum:                     0.5d0 * (kj|li) * t(kbla) +
                t2(iajb) = t2(iajb) + 0.5d0*MOI(kjli)*t2old(kbla)

              enddo
            enddo

!           write(*,*) "kc terms"
            do k=MoZero,nOcc
              do c=nOcc+1,dim_1e-1
!               if(a == b .and. c == b)cycle
!               if(i == j)then
!                 if(k == i .or. k == j)cycle
!               endif
                call Index2e(k,c,kc)
                call Index2e(b,j,bj)
                call Index2e(kc,bj,kcbj)
                call Index2e(i,a,ia)
                call Index2e(ia,kc,iakc)
                call Index2e(a,i,ai)
                call Index2e(kc,ai,kcai)
                call Index2e(j,b,jb)
                call Index2e(jb,kc,jbkc)

            !                         2.0d0 * (kc|bj)  * t(iakc) +
!                t2(iajb) = t2(iajb) + 2.0d0*(MOI(kcbj))*t2old(iakc)
            !                         2.0d0 * (kc|ai) * t(ijkc) +
!                t2(iajb) = t2(iajb) + 2.0d0*(MOI(kcai))*t2old(jbkc)

                call Index2e(i,c,ic)
                call Index2e(k,a,ka)
                call Index2e(ic,ka,icka)
                call Index2e(j,c,jc)
                call Index2e(k,b,kb)
                call Index2e(jc,kb,jckb)
            !sum:                    -1.0d0 * (kc|bj) *  t(icka) +
!                t2(iajb) = t2(iajb) - 1.0d0*(MOI(kcbj))*t2old(icka)
            !sum:                    -1.0d0 * (kc|ai) *  t(jckb) +
!                t2(iajb) = t2(iajb) - 1.0d0*(MOI(kcai))*t2old(jckb)


                call Index2e(k,i,ki)
                call Index2e(b,c,bc)
                call Index2e(ki,bc,kibc)
                call Index2e(k,a,ka)
                call Index2e(j,c,jc)
                call Index2e(ka,jc,kajc)
                call Index2e(k,j,kj)
                call Index2e(a,c,ac)
                call Index2e(kj,ac,kjac)
                call Index2e(k,b,kb)
                call Index2e(i,c,ic)
                call Index2e(kb,ic,kbic)

            !sum:                    -1.0d0 * (ki|bc) * t(kajc) +
!                t2(iajb) = t2(iajb) - 1.0d0*(MOI(kibc))*t2old(kajc)
            !sum:                    -1.0d0 * (kj|ac) * t(kbic) +
!                t2(iajb) = t2(iajb) - 1.0d0*(MOI(kjac))*t2old(kbic)


                call Index2e(k,j,kj)
                call Index2e(b,c,bc)
                call Index2e(kj,bc,kjbc)
                call Index2e(i,a,ia)
                call Index2e(k,c,kc)
                call Index2e(ia,kc,iakc)
                call Index2e(k,i,ki)
                call Index2e(a,c,ac)
                call Index2e(ki,ac,kiac)
                call Index2e(j,b,jb)
                call Index2e(k,c,kc)
                call Index2e(jb,kc,jbkc)
            !sum:                    -1.0d0 * (kj|bc) * t(iakc) +
!                t2(iajb) = t2(iajb) - 1.0d0*(MOI(kjbc))*t2old(iakc)
            !sum:                    -1.0d0 * (ki|ac) * t(jbkc) +
!                t2(iajb) = t2(iajb) - 1.0d0*(MOI(kiac))*t2old(jbkc)

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
            Elccd = Elccd - t2(iajb)*MOI(iajb)/Dijab
            Elccd = Elccd - t2(iajb)*(MOI(iajb)-MOI(ibja))/Dijab

          enddo
        enddo
      enddo
    enddo

    write(*,*) "      LCCD: ",iT2,Elccd 
 
  enddo ! t2 loop
  endif

!  call Fock_to_AO()

end subroutine calc_Elccd

end module module_cc
