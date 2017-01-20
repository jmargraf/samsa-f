module module_energy
  implicit none

! Energies           
  double precision                            :: Enuc     = 0.0d0
  double precision                            :: Eelec    = 0.0d0         
  double precision                            :: Etot     = 0.0d0          
  double precision                            :: Eold     = 0.0d0 

contains

!#############################################
!#            Calculate E_nuclear
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

