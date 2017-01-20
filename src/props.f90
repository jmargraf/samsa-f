module module_props
  implicit none

contains

!#############################################
!#            Population Analysis  
!#############################################
subroutine pop_Mulliken()
  use module_data,             only : natoms,dim_1e
  implicit none
!  integer                         :: i,j
!  double precision, allocatable   :: PopQ(:)

end subroutine pop_Mulliken



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
