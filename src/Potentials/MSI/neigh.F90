
! choose which data points to include in the neighbour list


subroutine neighbour(interp,pot,RawWeight,r,neigh)
  use interpolation
  implicit none

  type (interp_params), intent(in) :: interp
  type (pot_data_point), dimension(:), pointer :: pot
  real(kind=8), dimension(:), intent(out) :: RawWeight
  real(kind=8), dimension(:), intent(in) :: r
  type (neighbour_list), intent(out) :: neigh

  integer :: i
  real(kind=8) totsum,tol

  neigh%numInner=0    ! number of neighbours
  neigh%inner = 0      ! list of (inner) neighbours

  !----------------------------------------------------------
  ! calculate raw weights and totsum in two loops for speed
  !----------------------------------------------------------

  do i=1,interp%ndata
     RawWeight(i) = 1.0/(sum((r-pot(i)%r)**2)**interp%ipow)
  enddo

  totsum = sum(RawWeight)
  
  !----------------------------------------------------------
  !  build the inner neighbour list
  !----------------------------------------------------------

  tol = interp%wtol*totsum
  do i=1,interp%ndata
     if (RawWeight(i) > tol) then
       neigh%numInner = neigh%numInner + 1
       neigh%inner(neigh%numInner) = i
     endif
  enddo

  !----------------------------------------------------------

  return
end subroutine

