module module_props
  implicit none

contains

!#############################################
!#            Population Analysis  
!#############################################
subroutine pop_Mulliken()
  use module_data,             only : natoms,dim_1e,Dens,Spins,Basis,Sij
  use module_data,             only : atom
  use module_io,               only : print_Vec
  use module_constants,        only : CoreQ
  implicit none
  integer                         :: i,j,iSpin,intq
  double precision, allocatable   :: PopQ(:)
  double precision                :: q

  allocate(PopQ(nAtoms))

  PopQ = 0.0d0

  do iSpin = 1,Spins
    do i=1,dim_1e
      do j=1,dim_1e
        PopQ(Basis(i)) = PopQ(Basis(i)) - Dens(i,j,iSpin)*Sij(i,j)*1.0d0/dble(Spins)
        PopQ(Basis(j)) = PopQ(Basis(j)) - Dens(i,j,iSpin)*Sij(i,j)*1.0d0/dble(Spins)
      enddo
    enddo
  enddo 

  do i=1,nAtoms
    call CoreQ(atom(i),intq)
    q = dble(intq)
    PopQ(i) = PopQ(i) + q
  enddo

  call print_Vec(PopQ(:),nAtoms,10,"Q Mulliken")

end subroutine pop_Mulliken

!#############################################
!#        Löwdin Population Analysis  
!#############################################
subroutine pop_Loewdin()
  use module_data,             only : natoms,dim_1e,Spins,Basis,Coef,Occ
  use module_data,             only : atom,Loew
  use module_io,               only : print_Vec,print_Mat
  use module_constants,        only : CoreQ
  use module_wavefun,          only : build_S12,ort_coef,deort_fock
  implicit none
!  double precision, intent(IN)    :: Dens(dim_1e,dim_1e)
  integer                         :: i,j,k,iSpin,intq
  double precision, allocatable   :: PopQ(:)
  double precision                :: q

  allocate(PopQ(nAtoms))
  allocate(Loew(dim_1e,dim_1e,Spins))

  call build_S12()
  call ort_coef()

  do iSpin=1,Spins
    do i=1,dim_1e
      do j=1,i
        Loew(i,j,iSpin) = 0.0d0
        do k=1,dim_1e
          Loew(i,j,iSpin) = Loew(i,j,iSpin) + Occ(k,iSpin)*Coef(i,k,iSpin)*Coef(j,k,iSpin)*0.5d0*Spins
        enddo
        Loew(j,i,iSpin) = Loew(i,j,iSpin)
      enddo
    enddo
    !call print_Mat(Loew(:,:,iSpin),dim_1e,22,"Loewdin Density Matrix")
  enddo

  call deort_fock()

  PopQ = 0.0d0

  do iSpin = 1,Spins
    do i=1,dim_1e
!      do j=1,dim_1e
        PopQ(Basis(i)) = PopQ(Basis(i)) - Loew(i,i,iSpin)*2.0d0/dble(Spins)
!        PopQ(Basis(j)) = PopQ(Basis(j)) - Dens(i,j,iSpin)!*Sij(i,j)!*1.0d0/dble(Spins)
!      enddo
    enddo
  enddo

  do i=1,nAtoms
    call CoreQ(atom(i),intq)
    q = dble(intq)
    PopQ(i) = PopQ(i) + q
  enddo

  call print_Vec(PopQ(:),nAtoms,9,"Q Loewdin")

end subroutine pop_Loewdin


!#############################################
!#           Print Orbital Energies
!#############################################
subroutine print_Eigen()
  use module_data,            only : Eps,dim_1e,Spins
  use module_constants,       only : ha2eV
  use module_io,              only : print_SVec,print_Vec
  implicit none
  integer                         :: i,j
  double precision, allocatable   :: Eps_eV(:,:)

  allocate(Eps_eV(dim_1e,Spins))

  Eps_eV = Eps*ha2eV

  if(Spins==1)then
    call print_Vec(Eps(:,1),dim_1e,21,"Orbital Energies / Ha")
    call print_Vec(Eps_eV(:,1),dim_1e,21,"Orbital Energies / eV")
  elseif(Spins==2)then
    call print_SVec(Eps(:,:),dim_1e,21,"Orbital Energies / Ha")
    call print_SVec(Eps_eV(:,:),dim_1e,21,"Orbital Energies / eV")
  endif

end subroutine print_Eigen


end module module_props
