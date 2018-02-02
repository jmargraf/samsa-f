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
  use module_data,      only : DropMO,DoDrop,Fock,DoSingles,dim_1e,AMO
  use module_data,      only : Eps_SO,F_SO
  use module_ints,      only : Index2e
  use module_wavefun,   only : Fock_to_MO,Fock_to_AO
  implicit none
  integer                   :: iSpin,MOZero,iT2
  integer                   :: p3,p4,h1,h2,h5,h6,p5,p6

  integer                   :: nOcc,nmo
  double precision          :: occ_factor,Dijab,E_2,E_2OS,E_2SS
  double precision          :: DampCC = 0.5d0
  double precision          :: in1

  Elccd  = 0.0d0

!  call Fock_to_MO()

  if(spins==1)then
    write(*,*) "    "
    write(*,*) "    Calculating RHF-LCCD correlation energy"
    write(*,*) "    "

    nmo = dim_1e
    MOZero = 1


    if(DoDrop)then
      MOZero = DropMO+1
    else
      MOZero = 1
    endif
    nOcc =  nOccA

    allocate(t2(1:nmo,1:nmo,1:nmo,1:nmo),  &
          t2old(1:nmo,1:nmo,1:nmo,1:nmo))

    t2 = 0.0d0
    t2old = 0.0d0
    Dijab = 1.0d0

  do iT2=1,50 ! T2 loop
    Elccd  = 0.0d0
    t2old  = t2
    do h1=MOZero,nOcc
      do h2=MOZero,nOcc    
        do p3=nOcc+1,nmo
          do p4=nOcc+1,nmo

            in1 = 0.0d0

!           [ + 1.0 ] * v ( p3 p4 h1 h2 )
            in1 = 1.0d0 * AMO(h1,h2,p3,p4)

!           [ - 1.0 + 1.0 * P( p4 p3 h1 h2 => p3 p4 h1 h2 ) ] * Sum ( p5 ) * f ( p4 p5 ) * t ( p5 p3 h1 h2 ) done
            do p5 = nOcc+1,nmo
              in1 = in1 + 1.0d0 * F_SO(p4,p5) * t2old(p3,p5,h1,h2)
              in1 = in1 + 1.0d0 * F_SO(p3,p5) * t2old(p4,p5,h2,h1)
            enddo

!           [ - 1.0 + 1.0 * P( p3 p4 h1 h2 => p3 p4 h2 h1 ) ] * Sum ( h5 ) * f ( h5 h1 ) * t ( p3 p4 h5 h2 ) done
            do h5 = MOZero,nOcc
              !write(*,*) F_SO(h5,h2), F_SO(h5,h1)
              in1 = in1 - 1.0d0 * F_SO(h5,h2) * t2old(p3,p4,h1,h5)
              in1 = in1 - 1.0d0 * F_SO(h5,h1) * t2old(p4,p3,h2,h5)
            enddo

!           [ + 0.5 ] * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 ) * v ( h5 h6 h1 h2 ) done
            do h5 = MOZero,nOcc
              do h6 = MOZero,nOcc
                in1 = in1 + 0.5d0 * t2old(p3,p4,h5,h6) * AMO(h5,h6,h1,h2)
                in1 = in1 + 0.5d0 * t2old(p4,p3,h5,h6) * AMO(h5,h6,h2,h1)
              enddo
            enddo

!          [ + 1.0 
!            - 1.0 * P( p3 p4 h2 h1 => p4 p3 h2 h1 ) 
!            - 1.0 * P( p3 p4 h2 h1 => p3 p4 h1 h2 ) 
!            + 1.0 * P( p3 p4 h2 h1 => p4 p3 h1 h2 ) ] 
!          * Sum ( p5 h6 ) * t ( p5 p3 h6 h2 ) * v ( h6 p4 h1 p5 ) done
            do p5 = nOcc+1,nmo
              do h5 = MOZero,nOcc
                in1 = in1 + 2.0d0 * t2old(p3,p5,h1,h5) * AMO(h5,p4,p5,h2)
                in1 = in1 + 2.0d0 * t2old(p4,p5,h2,h5) * AMO(h5,p3,p5,h1)

                in1 = in1 - 1.0d0 * t2old(p5,p4,h1,h5) * AMO(h5,p4,p5,h2)
                in1 = in1 - 1.0d0 * t2old(p5,p3,h2,h5) * AMO(h5,p3,p5,h1)

                in1 = in1 - 1.0d0 * t2old(p3,p5,h5,h2) * AMO(h5,p4,h1,p5)
                in1 = in1 - 1.0d0 * t2old(p4,p5,h5,h1) * AMO(h5,p3,h2,p5)

                in1 = in1 - 1.0d0 * t2old(p3,p5,h1,h5) * AMO(h5,p4,h2,p5)
                in1 = in1 - 1.0d0 * t2old(p4,p5,h2,h5) * AMO(h5,p3,h1,p5)
              enddo
            enddo

!           [ + 0.5 ] * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 ) * v ( p3 p4 p5 p6 ) done
            do p5 = nOcc+1,nmo
              do p6 = nOcc+1,nmo
                in1 = in1 + 0.5d0 * t2old(p5,p6,h1,h2) * AMO(p3,p4,p5,p6)
                in1 = in1 + 0.5d0 * t2old(p5,p6,h2,h1) * AMO(p4,p3,p5,p6)
              enddo
            enddo


!            write(*,*) "t2   ",p3,p4,h1,h2,AMO(p3,p4,h1,h2),AMO(h1,h2,p3,p4)
!            Dijab = (Eps_SO(p3) + Eps_SO(p4) - Eps_SO(h1) - Eps_SO(h2))
            t2(p3,p4,h1,h2) = in1

          enddo
        enddo
      enddo
    enddo

    Elccd = 0.0d0
!    E_2 = 0.0d0
!    E_2OS = 0.0d0
!    E_2SS = 0.0d0

    do h1=MOZero,nOcc
      do h2=MOZero,nOcc
        do p3=nOcc+1,nmo
          do p4=nOcc+1,nmo

!            Dijab = Eps(a+1,1)+Eps(b+1,1)-Eps(i+1,1)-Eps(j+1,1)
!            Elccd = Elccd - t2(i,a,j,b)*MOI(iajb)/Dijab
!            Elccd = Elccd - t2(i,a,j,b)*(MOI(iajb)-MOI(ibja))/Dijab
            Dijab = (Eps_SO(p3) + Eps_SO(p4) - Eps_SO(h1) - Eps_SO(h2))
!            write(*,*) "Dijab",p3,p4,h1,h2,Dijab
            if(t2(p3,p4,h1,h2) /= 0.0d0)then

              Elccd = Elccd - 0.25d0*t2(p3,p4,h1,h2)*AMO(h1,h2,p3,p4)/Dijab
            endif
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
