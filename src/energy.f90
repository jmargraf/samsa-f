module module_energy
  implicit none

! Energies           
  double precision                            :: Enuc     = 0.0d0
  double precision                            :: Eelec    = 0.0d0         
  double precision                            :: Etot     = 0.0d0          
  double precision                            :: Eold     = 0.0d0 
  double precision                            :: Embpt2   = 0.0d0
  double precision                            :: Embpt2f  = 0.0d0
  double precision                            :: Edcpt2   = 0.0d0
  double precision                            :: Edcpt2f  = 0.0d0
  double precision                            :: Elmpt2   = 0.0d0
  double precision                            :: E_1,E_1f,E_OS,E_OSf,E_SS,E_SSx,E_SSc
  double precision                            :: E_AAx,E_AAc,E_BBx,E_BBc
  double precision                            :: E_AAxf,E_AAcf,E_BBxf,E_BBcf
  double precision                            :: E1e, E2e
  double precision, allocatable               :: LMmat(:,:,:)

contains

!#############################################
!#            Calculate E_mbpt2
!#############################################
subroutine calc_Embpt2
  use module_data,      only : MOI,Spins,dim_1e,nOccA,nOccB,Eps,SMO,Occ
  use module_data,      only : DropMO,DoDrop,Fock,DoSingles,DoReNorm
  use module_ints,      only : Index2e
  use module_wavefun,   only : Fock_to_MO,Fock_to_AO
  implicit none
  integer                   :: iSpin,i,j,a,b,ia,ja,jb,ib,iajb,ibja,MOZero
  double precision          :: occ_factor,occi,occj,occa,occb

  Embpt2 = 0.0d0
  E_OS   = 0.0d0
  E_SS   = 0.0d0
  E_SSx  = 0.0d0
  E_SSc  = 0.0d0
  E_AAx  = 0.0d0
  E_AAc  = 0.0d0
  E_BBx  = 0.0d0
  E_BBc  = 0.0d0
  E_1    = 0.0d0

  MOZero = 0

  call Fock_to_MO()

  if(spins==1)then
!    write(*,*) "    "
!    write(*,*) "    Calculating RHF-MBPT(2) correlation energy"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO
    else
      MOZero = 0
    endif

    do i=MOZero,nOccA-1
      do a=nOccA,dim_1e-1
        call Index2e(i,a,ia)
        if(doSingles)then
          E_1   = E_1 + Fock(i+1,a+1,1)*Fock(i+1,a+1,1)/(Eps(a+1,1)-Eps(i+1,1))
        endif

        do j=MOZero,nOccA-1
          call Index2E(j,a,ja)
          do b=nOccA,dim_1e-1
            call Index2e(j,b,jb)
            call Index2e(i,b,ib)
            call Index2e(ia,jb,iajb)
            call Index2e(ib,ja,ibja)

            E_OS  = E_OS  - (MOI(iajb)*MOI(iajb))/                 &
                            (Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1))

            E_SS  = E_SS  - (MOI(iajb)*(MOI(iajb)-MOI(ibja)))/     &
                            (Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1))

            E_SSx = E_SSx + (MOI(iajb)*MOI(ibja))/                 &
                            (Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1))

            E_SSc = E_SSc - (MOI(iajb)*MOI(iajb))/                 &
                            (Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1))

          enddo
        enddo
      enddo
    enddo

    Embpt2 = E_1 + E_OS + E_SS

  elseif(spins==2)then
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-MBPT(2) correlation energy"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif

    ! Opposite Spin 
    do i=MOZero,nOccA*2-1,2
      do a=nOccA*2+1,dim_1e*2-1,2
        do j=MOZero+1,nOccB*2,2
          do b=nOccB*2+2,dim_1e*2,2

            E_OS  = E_OS  - (SMO(i,a,j,b)*SMO(i,a,j,b))/            &
                            (Eps((a+1)/2,1)+Eps((b/2),2)            &
                            -Eps((i+1)/2,1)-Eps((j/2),2))
          enddo
        enddo
      enddo
    enddo

    ! Same Spin AA
    do i=MOZero,nOccA*2-1,2
      do a=nOccA*2+1,dim_1e*2-1,2
        if(doSingles)then
          E_1   = E_1 - Fock((i+1)/2,(a+1)/2,1)*Fock((i+1)/2,(a+1)/2,1)*0.5d0/  & 
                        (Eps((a+1)/2,1)-Eps((i+1)/2,1))
        endif

        do j=MOZero,nOccA*2-1,2
          do b=nOccA*2+1,dim_1e*2-1,2

            E_AAc = E_AAc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0/      &
                            (Eps((a+1)/2,1)+Eps((b+1)/2,1)          &
                            -Eps((i+1)/2,1)-Eps((j+1)/2,1))


            E_AAx = E_AAx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0/      &
                            (Eps((a+1)/2,1)+Eps((b+1)/2,1)          &
                            -Eps((i+1)/2,1)-Eps((j+1)/2,1))
          enddo
        enddo
      enddo
    enddo

    ! Same Spin BB
    do i=MOZero+1,nOccB*2,2
      do a=nOccB*2+2,dim_1e*2,2
        if(doSingles)then
          E_1   = E_1 - Fock(i/2,a/2,2)*Fock(i/2,a/2,2)*0.5d0/  & 
                        (Eps(a/2,2)-Eps(i/2,2))
        endif

        do j=MOZero+1,nOccB*2,2
          do b=nOccB*2+2,dim_1e*2,2

            E_BBc = E_BBc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0/      &
                            (Eps((a/2),2)+Eps((b/2),2)              &
                            -Eps((i/2),2)-Eps((j/2),2))


            E_BBx = E_BBx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0/      &
                            (Eps((a/2),2)+Eps((b/2),2)              &
                            -Eps((i/2),2)-Eps((j/2),2))
          enddo
        enddo
      enddo
    enddo

    Embpt2 = E_1 + E_OS + E_AAc + E_AAx + E_BBx + E_BBc

  endif

  if(spins==2)then
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-MBPT(2) correlation energy (nOcc version)"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO*2+1
    endif

    Embpt2f    = 0.0d0
    E_1f       = 0.0d0
    E_OSf      = 0.0d0
    E_AAxf     = 0.0d0
    E_AAcf     = 0.0d0
    E_BBxf     = 0.0d0
    E_BBcf     = 0.0d0
    occi       = 1.0d0
    occj       = 1.0d0
    occa       = 1.0d0
    occb       = 1.0d0

    ! Opposite Spin 
    do i=MOZero,dim_1e*2-1,2
      do a=MOZero,dim_1e*2-1,2
        do j=MOZero+1,dim_1e*2,2
          do b=MOZero+1,dim_1e*2,2
            occ_factor = (Occ((i+1)/2,1)*(1.0d0-Occ((a+1)/2,1))*Occ(j/2,2)*(1.0d0-Occ(b/2,2)))
            if(DoReNorm)then
              occi = Occ((i+1)/2,1)
              occa = (1.0d0 - Occ((a+1)/2,1))
              occj = Occ(j/2,2)
              occb = (1.0d0-Occ(b/2,2))
            endif
            if (occ_factor==0.0d0) cycle
            if (a == b .or. i == j) cycle
            E_OSf = E_OSf - (SMO(i,a,j,b)*SMO(i,a,j,b))*occ_factor/       &
                            (occa*Eps((a+1)/2,1)+occb*Eps((b/2),2)        &
                            -occi*Eps((i+1)/2,1)-occj*Eps((j/2),2))
          enddo
        enddo
      enddo
    enddo

    ! Same Spin AA
    do i=MOZero,dim_1e*2-1,2
      do a=MOZero,dim_1e*2-1,2
        occ_factor = (Occ((i+1)/2,1)*(1.0d0-Occ((a+1)/2,1)))
!        write(*,*) i,a,Eps((a+1)/2,1),Eps((i+1)/2,1)
        if(DoReNorm)then
          occi = Occ((i+1)/2,1)
          occa = (1.0d0 - Occ((a+1)/2,1))
        endif
        if (occ_factor==0.0d0) cycle
        if (doSingles .and. i/=a) then
          E_1f  = E_1f - Fock((i+1)/2,(a+1)/2,1)*Fock((i+1)/2,(a+1)/2,1)*occ_factor*0.5d0/  &
                         (occa*Eps((a+1)/2,1)-occi*Eps((i+1)/2,1))
        endif
        do j=MOZero,dim_1e*2-1,2
          do b=MOZero,dim_1e*2-1,2
            occ_factor = (Occ((i+1)/2,1)*(1.0d0-Occ((a+1)/2,1))*Occ((j+1)/2,1)*(1.0d0-Occ((b+1)/2,1)))
            if(DoReNorm)then
              occi = Occ((i+1)/2,1)
              occa = (1.0d0 - Occ((a+1)/2,1))
              occj = Occ((j+1)/2,1)
              occb = (1.0d0-Occ((b+1)/2,1))
            endif
            if (occ_factor==0.0d0) cycle
            if (a == b .or. i == j) cycle
            E_AAcf = E_AAcf- (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0*occ_factor/  &
                             (occa*Eps((a+1)/2,1)+occb*Eps((b+1)/2,1)       &
                             -occi*Eps((i+1)/2,1)-occj*Eps((j+1)/2,1))


            E_AAxf = E_AAxf+ (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0*occ_factor/  &
                             (occa*Eps((a+1)/2,1)+occb*Eps((b+1)/2,1)       &
                             -occi*Eps((i+1)/2,1)-occj*Eps((j+1)/2,1))
          enddo
        enddo
      enddo
    enddo


    ! Same Spin BB
    do i=MOZero+1,dim_1e*2,2
      do a=MOZero+1,dim_1e*2,2
        occ_factor = (Occ(i/2,2)*(1.0d0-Occ(a/2,2)))
        if(DoReNorm)then
          occi = Occ(i/2,2)
          occa = (1.0d0 - Occ(a/2,2))
        endif
        if (occ_factor==0.0d0) cycle
        if (doSingles .and. i/=a ) then
          E_1f  = E_1f - Fock(i/2,a/2,2)*Fock(i/2,a/2,2)*occ_factor*0.5d0/  &
                       (occa*Eps(a/2,2)-occi*Eps(i/2,2))
        endif
        do j=MOZero+1,dim_1e*2,2
          do b=MOZero+1,dim_1e*2,2
            occ_factor = (Occ(i/2,2)*(1.0d0-Occ(a/2,2))*Occ(j/2,2)*(1.0d0-Occ(b/2,2)))
            if(DoReNorm)then
              occi = Occ(i/2,2)
              occa = (1.0d0 - Occ(a/2,2))
              occj = Occ(j/2,2)
              occb = (1.0d0-Occ(b/2,2))
            endif
            if (occ_factor==0.0d0) cycle
            if (a == b .or. i == j) cycle
            E_BBcf= E_BBcf- (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0*occ_factor/  &
                            (occa*Eps((a/2),2)+occb*Eps((b/2),2)           &
                            -occi*Eps((i/2),2)-occj*Eps((j/2),2))


            E_BBxf= E_BBxf+ (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0*occ_factor/  &
                            (occa*Eps((a/2),2)+occb*Eps((b/2),2)           &
                            -occi*Eps((i/2),2)-occj*Eps((j/2),2))
          enddo
        enddo
      enddo
    enddo

    Embpt2f = E_1f + E_OSf + E_AAcf + E_AAxf + E_BBxf + E_BBcf

  endif

  call Fock_to_AO()

end subroutine calc_Embpt2


!#############################################
!#            Calculate E_lmpt2
!#############################################
subroutine calc_Elmpt2
  use module_data,      only : MOI,Spins,dim_1e,nOccA,nOccB,Eps,Occ,scaleOS,scaleSS
  use module_data,      only : DropMO,DoDrop,Fock,DoSingles,scaleJ,scaleK
  use module_ints,      only : Index2e
  use module_wavefun,   only : Fock_to_MO,Fock_to_AO
  implicit none
  integer                   :: iSpin,i,j,a,b,MOZero
  double precision          :: moJ,moK
!  double precision          :: occ_factor,occi,occj,occa,occb

  Elmpt2 = 0.0d0
  E_OS   = 0.0d0
  E_SS   = 0.0d0
  E_SSx  = 0.0d0
  E_SSc  = 0.0d0

  MOZero = 1

  call calc_lmmat()

  if(spins==1)then
!    write(*,*) "    "
!    write(*,*) "    Calculating RHF-LMPT(2) correlation energy"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO+1
    else
      MOZero = 1
    endif

    do i=MOZero,nOccA
      do a=nOccA+1,dim_1e

        do j=MOZero,nOccA
          do b=nOccA+1,dim_1e
            moJ = 0.0d0
            moK = 0.0d0
            call calc_lmint(i,a,j,b,scaleJ,moJ)
            call calc_lmint(i,b,j,a,scaleK,moK)

            E_OS  = E_OS  - (moJ*moJ)/                 &
                            (Eps(a,1)+Eps(b,1)-Eps(i,1)-Eps(j,1))

            E_SS  = E_SS  - (moJ*(moJ-moK))/           &
                            (Eps(a,1)+Eps(b,1)-Eps(i,1)-Eps(j,1))

            E_SSx = E_SSx + (moJ*moK)/                 &
                            (Eps(a,1)+Eps(b,1)-Eps(i,1)-Eps(j,1))

            E_SSc = E_SSc - (moJ*moJ)/                 &
                            (Eps(a,1)+Eps(b,1)-Eps(i,1)-Eps(j,1))

          enddo
        enddo
      enddo
    enddo

    Elmpt2 = scaleOS*E_OS + scaleSS*E_SS

  endif

  deallocate(LMmat)

end subroutine calc_Elmpt2

!#############################################
!#            Calculate lmint
!#############################################
subroutine calc_lmint(p,q,r,s,alp,integral)
  use module_data,          only : nAtoms,atom
  use module_constants,     only : hardness,coreq
  implicit none
  integer, intent(in)           :: p,q,r,s
  double precision, intent(in)  :: alp
  double precision, intent(out) :: integral
  integer                       :: A,B,ZA,ZB
  double precision              :: eta ! average hardness
  double precision              :: Gam

  integral = 0.0d0

  do A=1,nAtoms
    do B=1,nAtoms
      call coreq(atom(A),ZA)
      call coreq(atom(B),ZB)
      eta = 0.5d0*(hardness(ZA)+hardness(ZB))
      call calc_klopman(A,B,eta,alp,Gam)
      integral = integral + LMmat(p,q,A)*LMmat(r,s,B)*Gam
    enddo
  enddo

end subroutine calc_lmint


!#############################################
!#         Calculate Klopman Damping
!#############################################
subroutine calc_klopman(A,B,eta,alp,Gam)
  use module_data,          only : rAB
  implicit none
  integer, intent(in)           :: A,B
  double precision, intent(in)  :: eta,alp
  double precision, intent(out) :: Gam

!  Gam = 1.0d0/(rAB(A,B)+(1.0d0/eta))
  Gam = (1.0d0/( (rAB(A,B)**alp) + (eta**(-alp)) ))**(1.0d0/alp)

end subroutine calc_klopman


!#############################################
!#    Calculate transition density matrix
!#############################################
subroutine calc_lmmat()
  use module_data,             only : natoms,dim_1e,Spins,Basis,Coef
  use module_data,             only : atom
  use module_wavefun,          only : build_S12,ort_coef,deort_fock
  use module_io,               only : print_Mat
  implicit none
  integer                         :: i,j,k,iSpin,A

  allocate(LMmat(dim_1e,dim_1e,nAtoms))

  call build_S12()
  call ort_coef()


  do iSpin=1,Spins
    do A=1,nAtoms

      do i=1,dim_1e !MO
        do j=1,dim_1e !MO
          LMmat(i,j,A) = 0.0d0
          do k=1,dim_1e !AO
            if(Basis(k)/=A) cycle
            LMmat(i,j,A) = LMmat(i,j,A) + Coef(k,i,iSpin)*Coef(k,j,iSpin) *2.0d0/Spins
          enddo
        enddo
      enddo
      !call print_Mat(LMmat(:,:,A),dim_1e,22,"Loewdin Transition Density Matrix")
    enddo
  enddo

  call deort_fock()

end subroutine calc_lmmat


!#############################################
!#            Calculate E_dcpt2
!#############################################
subroutine calc_Edcpt2
  use module_data,      only : MOI,Spins,dim_1e,nOccA,nOccB,Eps,SMO,Occ
  use module_data,      only : DropMO,DoDrop
  use module_ints,      only : Index2e
  implicit none
  integer                   :: iSpin,i,j,a,b,ia,ja,jb,ib,iajb,ibja,MOZero
  double precision          :: occ_factor
  double precision          :: Dabij,Vijab,Vibja,Viajb

  Edcpt2 = 0.0d0
  E_OS   = 0.0d0
  E_SS   = 0.0d0
  E_SSx  = 0.0d0
  E_SSc  = 0.0d0
  E_AAx  = 0.0d0
  E_AAc  = 0.0d0
  E_BBx  = 0.0d0
  E_BBc  = 0.0d0

  MOZero = 0

  if(spins==1)then
!    write(*,*) "    "
!    write(*,*) "    Calculating RHF-DCPT(2) correlation energy"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO
    else
      MOZero = 0
    endif

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

            Dabij = Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1)
            Viajb = MOI(iajb) 
            Vibja = MOI(ibja)

            E_OS  = E_OS  + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/2.0d0

            E_SSc = E_SSc + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/2.0d0
            E_SSx = E_SSx - (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Vibja))/2.0d0

          enddo
        enddo
      enddo
    enddo

    E_SS   = E_SSx + E_SSc
    Edcpt2 = E_OS  + E_SS

  elseif(spins==2)then
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-DCPT(2) correlation energy"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO*2+1
    else
      MOZero = 1
    endif

    ! Opposite Spin 
    do i=MOZero,nOccA*2-1,2
      do a=nOccA*2+1,dim_1e*2-1,2
        do j=MOZero+1,nOccB*2,2
          do b=nOccB*2+2,dim_1e*2,2

             Dabij = Eps((a+1)/2,1)+Eps((b/2),2)-Eps((i+1)/2,1)-Eps((j/2),2)
             Viajb = SMO(i,a,j,b)

             E_OS  = E_OS + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/2.0d0
          enddo
        enddo
      enddo
    enddo

    ! Same Spin AA
    do i=MOZero,nOccA*2-1,2
      do a=nOccA*2+1,dim_1e*2-1,2
        do j=MOZero,nOccA*2-1,2
          do b=nOccA*2+1,dim_1e*2-1,2

            Dabij = Eps((a+1)/2,1)+Eps((b+1)/2,1)-Eps((i+1)/2,1)-Eps((j+1)/2,1)
            Viajb = SMO(i,a,j,b)
            Vibja = SMO(i,b,j,a)

            E_AAc = E_AAc + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/4.0d0
            E_AAx = E_AAx - (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Vibja))/4.0d0
          enddo
        enddo
      enddo
    enddo


    ! Same Spin BB
    do i=MOZero+1,nOccB*2,2
      do a=nOccB*2+2,dim_1e*2,2
        do j=MOZero+1,nOccB*2,2
          do b=nOccB*2+2,dim_1e*2,2

            Dabij = Eps((a/2),2)+Eps((b/2),2)-Eps((i/2),2)-Eps((j/2),2) 
            Viajb = SMO(i,a,j,b)
            Vibja = SMO(i,b,j,a)

            E_BBc = E_BBc + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/4.0d0
            E_BBx = E_BBx - (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Vibja))/4.0d0
          enddo
        enddo
      enddo
    enddo

    Edcpt2 = E_OS + E_AAc + E_AAx + E_BBx + E_BBc
  endif

  if(spins==2)then
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-DCPT(2) correlation energy (nOcc version)"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO*2+1
    endif

    Edcpt2f  = 0.0d0
    E_OSf    = 0.0d0
    E_AAxf   = 0.0d0
    E_AAcf   = 0.0d0
    E_BBxf   = 0.0d0
    E_BBcf   = 0.0d0

    ! Opposite Spin 
    do i=MOZero,dim_1e*2-1,2
      do a=MOZero,dim_1e*2-1,2
        do j=MOZero+1,dim_1e*2,2
          do b=MOZero+1,dim_1e*2,2
            occ_factor = (Occ((i+1)/2,1)*(1-Occ((a+1)/2,1))*Occ(j/2,2)*(1-Occ(b/2,2)))
            if (occ_factor==0.0d0) cycle
            if (a == b .or. i == j) cycle
!           E_OS  = E_OS  - (SMO(i,a,j,b)*SMO(i,a,j,b))*occ_factor/       &
!                           (Eps((a+1)/2,1)+Eps((b/2),2)                &
!                           -Eps((i+1)/2,1)-Eps((j/2),2))
             Dabij = Eps((a+1)/2,1)+Eps((b/2),2)-Eps((i+1)/2,1)-Eps((j/2),2)
             Viajb = SMO(i,a,j,b)

             E_OSf = E_OSf+ (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb*occ_factor))/2.0d0
          enddo
        enddo
      enddo
    enddo

    ! Same Spin AA
    do i=MOZero,dim_1e*2-1,2
      do a=MOZero,dim_1e*2-1,2
        do j=MOZero,dim_1e*2-1,2
          do b=MOZero,dim_1e*2-1,2
            occ_factor = (Occ((i+1)/2,1)*(1-Occ((a+1)/2,1))*Occ((j+1)/2,1)*(1-Occ((b+1)/2,1)))
            if (occ_factor==0.0d0) cycle
            if (a == b .or. i == j) cycle
            Dabij = Eps((a+1)/2,1)+Eps((b+1)/2,1)-Eps((i+1)/2,1)-Eps((j+1)/2,1)
            Viajb = SMO(i,a,j,b)
            Vibja = SMO(i,b,j,a)

            E_AAcf= E_AAcf+ (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb*occ_factor))/4.0d0
            E_AAxf= E_AAxf- (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Vibja*occ_factor))/4.0d0

!           E_AAc = E_AAc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0*occ_factor/  &
!                           (Eps((a+1)/2,1)+Eps((b+1)/2,1)                 &
!                           -Eps((i+1)/2,1)-Eps((j+1)/2,1))


!           E_AAx = E_AAx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0*occ_factor/  &
!                           (Eps((a+1)/2,1)+Eps((b+1)/2,1)                 &
!                           -Eps((i+1)/2,1)-Eps((j+1)/2,1))

          enddo
        enddo
      enddo
    enddo


    ! Same Spin BB
    do i=MOZero+1,dim_1e*2,2
      do a=MOZero+1,dim_1e*2,2
        do j=MOZero+1,dim_1e*2,2
          do b=MOZero+1,dim_1e*2,2
            occ_factor = (Occ(i/2,2)*(1-Occ(a/2,2))*Occ(j/2,2)*(1-Occ(b/2,2)))
            if (occ_factor==0.0d0) cycle
            if (a == b .or. i == j) cycle
            Dabij = Eps((a/2),2)+Eps((b/2),2)-Eps((i/2),2)-Eps((j/2),2)
            Viajb = SMO(i,a,j,b)
            Vibja = SMO(i,b,j,a)

            E_BBcf= E_BBcf+ (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb*occ_factor))/4.0d0
            E_BBxf= E_BBxf- (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Vibja*occ_factor))/4.0d0

!           E_BBc = E_BBc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0*occ_factor/  &
!                           (Eps((a/2),2)+Eps((b/2),2)                     &
!                           -Eps((i/2),2)-Eps((j/2),2))


!           E_BBx = E_BBx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0*occ_factor/  &
!                           (Eps((a/2),2)+Eps((b/2),2)                     &
!                           -Eps((i/2),2)-Eps((j/2),2))

          enddo
        enddo
      enddo
    enddo

   Edcpt2f  = E_OSf+ E_AAcf+ E_AAxf+ E_BBxf+ E_BBcf

  endif


end subroutine calc_Edcpt2


!#############################################
!#            Calculate E_total
!#############################################
subroutine calc_Energy
  use module_data, only      : Dens,Hcore,Fock,Spins,dim_1e
  implicit none
  integer                   :: iSpin,i,j 

  Eelec = 0.0d0
  Etot  = 0.0d0
  E1e   = 0.0d0
  E2e   = 0.0d0

  do iSpin=1,Spins
    do i=1,dim_1e
      do j=1,dim_1e
        E1e = E1e + Dens(i,j,iSpin)*Hcore(i,j)*(2.0d0/dble(Spins))
        E2e = E2e + Dens(i,j,iSpin)*(Fock(i,j,iSpin)-Hcore(i,j))*(1.0d0/dble(Spins))          
      enddo
    enddo
  enddo

  Eelec = E1e + E2e

  Etot = Eelec + Enuc

!  write(*,*) "Eelec = ", Eelec
!  write(*,*) "Etot  = ", Etot


end subroutine calc_Energy


!#############################################
!#            Calculate E_nuclear
!#############################################
subroutine calc_Enuc()
  use module_constants, only : CoreQ
  use module_data, only      : rAB,atom,xyz,natoms
  implicit none
  integer                   :: i,j
  integer                   :: q1,q2

  write(*,*) ""
  write(*,*) "  Distance matrix:"
  write(*,*) ""

  Enuc = 0.0d0
  do i=1,natoms
    call CoreQ(atom(i),q1)
    rAB(i,i)=0.0d0
    do j=i+1,natoms
      call CoreQ(atom(j),q2)
      rAB(i,j) = sqrt((xyz(i,1)-xyz(j,1))**2+ &
                      (xyz(i,2)-xyz(j,2))**2+ &
                      (xyz(i,3)-xyz(j,3))**2  )
!      write(*,*) "q1 q2: ",q1,q2
      Enuc = Enuc + (dble(q1)*dble(q2))/rAB(i,j)
      write(*,*) "    rAB ",atom(i),i,atom(j),j," = ",rAB(i,j)
      rAB(j,i)=rAB(i,j)
    enddo
  enddo

  write(*,*) "    Enuc = ", Enuc
  write(*,*) ""

end subroutine calc_Enuc


end module module_energy

