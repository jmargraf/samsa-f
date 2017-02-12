module module_energy
  implicit none

! Energies           
  double precision                            :: Enuc     = 0.0d0
  double precision                            :: Eelec    = 0.0d0         
  double precision                            :: Etot     = 0.0d0          
  double precision                            :: Eold     = 0.0d0 
  double precision                            :: Emp2     = 0.0d0
  double precision                            :: Emp2frac = 0.0d0
  double precision                            :: Edcpt2   = 0.0d0
  double precision                            :: Edcpt2fr = 0.0d0
  double precision                            :: E_OS,E_SS,E_SSx,E_SSc
  double precision                            :: E_AAx,E_AAc,E_BBx,E_BBc

contains

!#############################################
!#            Calculate E_mp2
!#############################################
subroutine calc_Emp2
  use module_data,      only : MOI,Spins,dim_1e,nOccA,nOccB,Eps,SMO,Occ
  use module_data,      only : DropMO,DoDrop
  use module_ints,      only : Index2e
  implicit none
  integer                   :: iSpin,i,j,a,b,ia,ja,jb,ib,iajb,ibja,MOZero
  double precision          :: occ_factor

  Emp2  = 0.0d0
  E_OS  = 0.0d0
  E_SS  = 0.0d0
  E_SSx = 0.0d0
  E_SSc = 0.0d0
  E_AAx = 0.0d0
  E_AAc = 0.0d0
  E_BBx = 0.0d0
  E_BBc = 0.0d0

  MOZero = 0

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

    Emp2 = E_OS + E_SS

!    write(*,*) "    "
!    write(*,*) "      E_OS  = ", E_OS
!    write(*,*) "      E_SS  = ", E_SS
!    write(*,*) "      E_SSx = ", E_SSx
!    write(*,*) "      E_SSc = ", E_SSc
!    write(*,*) "    "
!    write(*,*) "      Emp2  = ", Emp2
!    write(*,*) "    "


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

    Emp2 = E_OS + E_AAc + E_AAx + E_BBx + E_BBc

!    write(*,*) "    "
!    write(*,*) "      E_OS  = ", E_OS
!    write(*,*) "      E_AAx = ", E_AAx
!    write(*,*) "      E_AAc = ", E_AAc
!    write(*,*) "      E_BBx = ", E_BBx
!    write(*,*) "      E_BBc = ", E_BBc
!    write(*,*) "    "
!    write(*,*) "      Emp2  = ", Emp2
!    write(*,*) "    "

  endif

  if(spins==2)then
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-MBPT(2) correlation energy (nOcc version)"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO*2+1
    endif

    Emp2frac = 0.0d0
    E_OS     = 0.0d0
    E_SS     = 0.0d0
    E_SSx    = 0.0d0
    E_SSc    = 0.0d0
    E_AAx    = 0.0d0
    E_AAc    = 0.0d0
    E_BBx    = 0.0d0
    E_BBc    = 0.0d0

    ! Opposite Spin 
    do i=MOZero,dim_1e*2-1,2
      do a=MOZero,dim_1e*2-1,2
        do j=MOZero+1,dim_1e*2,2
          do b=MOZero+1,dim_1e*2,2
            occ_factor = (Occ((i+1)/2,1)*(1-Occ((a+1)/2,1))*Occ(j/2,2)*(1-Occ(b/2,2)))
            if (occ_factor==0.0d0) cycle
            E_OS  = E_OS  - (SMO(i,a,j,b)*SMO(i,a,j,b))*occ_factor/       &
                            (Eps((a+1)/2,1)+Eps((b/2),2)                &
                            -Eps((i+1)/2,1)-Eps((j/2),2))


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
            E_AAc = E_AAc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0*occ_factor/  &
                            (Eps((a+1)/2,1)+Eps((b+1)/2,1)                 &
                            -Eps((i+1)/2,1)-Eps((j+1)/2,1))


            E_AAx = E_AAx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0*occ_factor/  &
                            (Eps((a+1)/2,1)+Eps((b+1)/2,1)                 &
                            -Eps((i+1)/2,1)-Eps((j+1)/2,1))

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
            E_BBc = E_BBc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0*occ_factor/  &
                            (Eps((a/2),2)+Eps((b/2),2)                     &
                            -Eps((i/2),2)-Eps((j/2),2))


            E_BBx = E_BBx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0*occ_factor/  &
                            (Eps((a/2),2)+Eps((b/2),2)                     &
                            -Eps((i/2),2)-Eps((j/2),2))

          enddo
        enddo
      enddo
    enddo

    Emp2frac = E_OS + E_AAc + E_AAx + E_BBx + E_BBc

!    write(*,*) "    "
!    write(*,*) "      E_OS  = ", E_OS
!    write(*,*) "      E_AAx = ", E_AAx
!    write(*,*) "      E_AAc = ", E_AAc
!    write(*,*) "      E_BBx = ", E_BBx
!    write(*,*) "      E_BBc = ", E_BBc
!    write(*,*) "    "
!    write(*,*) "      Emp2  = ", Emp2frac
!    write(*,*) "    "

  endif


end subroutine calc_Emp2


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

!            Edcpt2 = Edcpt2 + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Vijab*Vijab))

!            E_OS  = E_OS  - (MOI(iajb)*MOI(iajb))/                 

!            E_SS  = E_SS  - (MOI(iajb)*(MOI(iajb)-MOI(ibja)))/     

!            E_SSx = E_SSx + (MOI(iajb)*MOI(ibja))/                 

!            E_SSc = E_SSc - (MOI(iajb)*MOI(iajb))/                 

            E_OS  = E_OS  + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/2.0d0
!
            E_SSc = E_SSc + (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Viajb))/2.0d0
            E_SSx = E_SSx - (Dabij - sqrt(Dabij*Dabij + 4.0d0*Viajb*Vibja))/2.0d0

          enddo
        enddo
      enddo
    enddo

    E_SS   = E_SSx + E_SSc
    Edcpt2 = E_OS  + E_SS

!    write(*,*) "    "
!    write(*,*) "      E_OS  = ", E_OS
!    write(*,*) "      E_SS  = ", E_SS
!    write(*,*) "      E_SSx = ", E_SSx
!    write(*,*) "      E_SSc = ", E_SSc
!    write(*,*) "    "
!    write(*,*) "      Edcpt2 = ", Edcpt2
!    write(*,*) "    "


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
!            E_OS  = E_OS  - (SMO(i,a,j,b)*SMO(i,a,j,b))/            &
!                            (Eps((a+1)/2,1)+Eps((b/2),2)            &
!                            -Eps((i+1)/2,1)-Eps((j/2),2))
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


!            E_AAc = E_AAc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0/      &
!                            (Eps((a+1)/2,1)+Eps((b+1)/2,1)          &
!                            -Eps((i+1)/2,1)-Eps((j+1)/2,1))


!            E_AAx = E_AAx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0/      &
!                            (Eps((a+1)/2,1)+Eps((b+1)/2,1)          &
!                            -Eps((i+1)/2,1)-Eps((j+1)/2,1))

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

!            E_BBc = E_BBc - (SMO(i,a,j,b)*SMO(i,a,j,b))*0.5d0/      &
!                            (Eps((a/2),2)+Eps((b/2),2)              &
!                            -Eps((i/2),2)-Eps((j/2),2))


!            E_BBx = E_BBx + (SMO(i,a,j,b)*SMO(i,b,j,a))*0.5d0/      &
!                            (Eps((a/2),2)+Eps((b/2),2)              &
!                            -Eps((i/2),2)-Eps((j/2),2))

          enddo
        enddo
      enddo
    enddo

    Edcpt2 = E_OS + E_AAc + E_AAx + E_BBx + E_BBc

!    write(*,*) "    "
!    write(*,*) "      E_OS  = ", E_OS
!    write(*,*) "      E_AAx = ", E_AAx
!    write(*,*) "      E_AAc = ", E_AAc
!    write(*,*) "      E_BBx = ", E_BBx
!    write(*,*) "      E_BBc = ", E_BBc
!    write(*,*) "    "
!    write(*,*) "      Emp2  = ", Emp2
!    write(*,*) "    "

  endif

  if(spins==2)then
!    write(*,*) "    "
!    write(*,*) "    Calculating UHF-DCPT(2) correlation energy (nOcc version)"
!    write(*,*) "    "

    if(DoDrop)then
      MOZero = DropMO*2+1
    endif

    Edcpt2fr = 0.0d0
!   E_OS     = 0.0d0
!   E_SS     = 0.0d0
!   E_SSx    = 0.0d0
!   E_SSc    = 0.0d0
!   E_AAx    = 0.0d0
!   E_AAc    = 0.0d0
!   E_BBx    = 0.0d0
!   E_BBc    = 0.0d0

    ! Opposite Spin 
    do i=MOZero,dim_1e*2-1,2
      do a=MOZero,dim_1e*2-1,2
        do j=MOZero+1,dim_1e*2,2
          do b=MOZero+1,dim_1e*2,2
            occ_factor = (Occ((i+1)/2,1)*(1-Occ((a+1)/2,1))*Occ(j/2,2)*(1-Occ(b/2,2)))
            if (occ_factor==0.0d0) cycle
!           E_OS  = E_OS  - (SMO(i,a,j,b)*SMO(i,a,j,b))*occ_factor/       &
!                           (Eps((a+1)/2,1)+Eps((b/2),2)                &
!                           -Eps((i+1)/2,1)-Eps((j/2),2))


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

!   Edcpt2fr = E_OS + E_AAc + E_AAx + E_BBx + E_BBc

!    write(*,*) "    "
!    write(*,*) "      E_OS  = ", E_OS
!    write(*,*) "      E_AAx = ", E_AAx
!    write(*,*) "      E_AAc = ", E_AAc
!    write(*,*) "      E_BBx = ", E_BBx
!    write(*,*) "      E_BBc = ", E_BBc
!    write(*,*) "    "
!    write(*,*) "      Emp2  = ", Emp2frac
!    write(*,*) "    "

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

  do iSpin=1,Spins
    do i=1,dim_1e
      do j=1,dim_1e
        Eelec = Eelec + Dens(i,j,iSpin)                & 
                      *(Hcore(i,j) + Fock(i,j,iSpin))  &
                      *(1.0d0/dble(Spins))
      enddo
    enddo
  enddo

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
    do j=i+1,natoms
      call CoreQ(atom(j),q2)
      rAB(i,j) = sqrt((xyz(i,1)-xyz(j,1))**2+ &
                      (xyz(i,2)-xyz(j,2))**2+ &
                      (xyz(i,3)-xyz(j,3))**2  )
!      write(*,*) "q1 q2: ",q1,q2
      Enuc = Enuc + (dble(q1)*dble(q2))/rAB(i,j)
      write(*,*) "    rAB ",atom(i),i,atom(j),j," = ",rAB(i,j)
    enddo
  enddo

  write(*,*) "    Enuc = ", Enuc
  write(*,*) ""

end subroutine calc_Enuc


end module module_energy

