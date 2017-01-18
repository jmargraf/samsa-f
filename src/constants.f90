module module_constants
  implicit none

! Calcularion details
  double precision, parameter                 :: ha2eV = 27.21d0
  double precision, parameter                 :: ang2bohr = 1.889725989d0
  double precision, parameter                 :: pi = 3.14159265359d0

contains

!#############################################
!#             Get Core Charge
!#############################################
subroutine coreq(El,Q)
  implicit none
  character, intent(in)                       :: El
  integer, intent(out)                        :: Q

  if(El == 'H')then
    Q = 1
  elseif(El == 'He')then
    Q = 2
  elseif(El == 'Li')then
    Q = 3
  elseif(El == 'Be')then
    Q = 4
  elseif(El == 'B')then
    Q = 5
  elseif(El == 'C')then
    Q = 6
  elseif(El == 'N')then
    Q = 7
  elseif(El == 'O')then
    Q = 8
  elseif(El == 'F')then
    Q = 9
  elseif(El == 'Ne')then
    Q = 10
  elseif(El == 'Na')then
    Q = 11
  elseif(El == 'Mg')then
    Q = 12
  elseif(El == 'Al')then
    Q = 13
  elseif(El == 'Si')then
    Q = 14
  elseif(El == 'P')then
    Q = 15
  elseif(El == 'S')then
    Q = 16
  elseif(El == 'Cl')then
    Q = 17
  elseif(El == 'Ar')then
    Q = 18
  else
    write(*,*) "  !Error: Element ", El, " not found!"
  endif

end subroutine coreq

!#############################################
!#            Get Atomic Radius 
!#############################################
subroutine arad(El,Rad)
  implicit none
  character, intent(in)                       :: El
  double precision, intent(out)               :: Rad

  if(El == 'H')then
    Rad = 1.0015547741d0
  elseif(El == 'He')then
    Rad = 0.585815056d0
  elseif(El == 'Li')then
    Rad = 3.155842401d0
  elseif(El == 'Be')then
    Rad = 2.116493107d0
  elseif(El == 'B')then
    Rad = 1.64406161d0
  elseif(El == 'C')then
    Rad = 1.266116412d0
  elseif(El == 'N')then
    Rad = 1.058246554d0
  elseif(El == 'O')then
    Rad = 0.907068475d0
  elseif(El == 'F')then
    Rad = 0.793684915d0
  elseif(El == 'Ne')then
    Rad = 0.718095876d0
  elseif(El == 'Na')then
    Rad = 3.5904796544504d0
  elseif(El == 'Mg')then
    Rad = 2.7401028941858d0
  elseif(El == 'Al')then
    Rad = 2.2298768380271d0
  elseif(El == 'Si')then
    Rad = 2.0975960086526d0
  elseif(El == 'P')then
    Rad = 1.8519316112428d0
  elseif(El == 'S')then
    Rad = 1.6629589978507d0
  elseif(El == 'Cl')then
    Rad = 1.4928836457978d0
  elseif(El == 'Ar')then
    Rad = 1.3417055550841d0
  else
    write(*,*) "  !Error: Element ", El, " not found!"
  endif

end subroutine arad


end module module_constants



