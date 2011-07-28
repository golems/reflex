
Subroutine rfx_lqg_observe( n_x, n_u, n_z, A, B, C, x, u, z, K, dx, zwork )
  Implicit None
  integer, intent(in) :: n_x, n_u, n_z                  ! state/input space
  real(8), intent(in), dimension(n_x) :: x
  real(8), intent(in), dimension(n_u) :: u
  real(8), intent(in), dimension(n_z) :: z
  real(8), intent(in), dimension(n_x,n_x) :: A
  real(8), intent(in), dimension(n_x,n_u) :: B
  real(8), intent(in), dimension(n_z,n_x) :: C
  real(8), intent(in), dimension(n_x,n_z) :: K
  real(8), intent(out), dimension(n_x) :: dx
  real(8), intent(out), dimension(n_z) :: zwork
  !! dx = A*x + B*u + K*(z-C*x)
  !! gfortran is pretty dumb

  !dx = matmul(A,x)
  !dx = dx + matmul(B,u)
  !zwork = matmul(C,x)
  !zwork = z - zwork

  dx = matmul(A,x) + matmul(B,u) + matmul(K, z - matmul(C,x))
End Subroutine rfx_lqg_observe
