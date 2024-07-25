program vegetation_resistance
   implicit none
   integer, parameter :: N = 400
   double precision, dimension(N) :: u_values, theta_values, F_D_values, compute_k_values, compute_value_values
   integer :: i
   double precision :: pi, g, rho, C_D, C_DL, A_1, l, d, C_SL
   character(len=50) :: filename
   pi = 4.0d0 * atan(1.0d0)
   filename = 'vegetation_data.csv'

   ! 物理定数とパラメータの定義
   g = 9.8d0
   rho = 1.0d0
   C_D = 1.0d0
   C_DL = 0.095d0
   A_1 = 0.096d0
   l = 0.3d0
   d = 0.01d0
   C_SL = 0.015d0

   ! u_valuesの生成
   do i = 1, N
      u_values(i) = (i-1) * 10.0d0 / (N-1)
   end do

   ! theta_valuesの計算
   do i = 1, N
      call theta_function(u_values(i), theta_values(i))
   end do

   ! F_D_valuesの計算
   do i = 1, N
      call F_D(rho, C_D, C_DL, A_1, l, d, u_values(i), 1.0d0, theta_values(i), C_SL, F_D_values(i))
   end do

   ! compute_k_valuesの計算
   do i = 1, N
      call compute_k(g, rho, u_values(i), 1.0d0, 10, F_D_values(i), compute_k_values(i))
   end do

   ! compute_value_valuesの計算
   do i = 1, N
      call compute_value(g, u_values(i), 1.0d0, compute_k_values(i), compute_value_values(i))
   end do

   ! CSVファイルに結果を書き込む
   call write_to_csv(filename, u_values, compute_value_values)

 contains

   subroutine theta_function(u, theta)
     double precision, intent(in) :: u
     double precision, intent(out) :: theta
     double precision, parameter :: pi = 3.141592653589793d0
     theta = ((pi / 2) - asin(0.2d0)) * ((u - 0.12d0)**2)
     if (theta >= ((pi / 2) - asin(0.2d0))) then
        theta = ((pi / 2) - asin(0.2d0))
     elseif (u <= 0.12d0) then
        theta = 0.0d0
     end if
   end subroutine theta_function

   subroutine F_D(rho, C_D, C_DL, A_1, l, d, u, v, theta, C_SL, F_D_value)
     double precision, intent(in) :: rho, C_D, C_DL, A_1, l, d, u, v, theta, C_SL
     double precision, intent(out) :: F_D_value
     F_D_value = 0.5d0 * rho * (C_D + C_DL * A_1 / (2d0 * l * d)) * (u**2 + v**2) * cos(theta)**2 &
               + 0.5d0 * rho * C_SL * A_1 * (u**2 + v**2) / cos(theta)
   end subroutine F_D

   subroutine compute_k(g, rho, u, v, N, F_D, k)
     double precision, intent(in) :: g, rho, u, v, F_D
     integer, intent(in) :: N
     double precision, intent(out) :: k
     k = (g * rho * (u**2 + v**2) / (N * F_D))**0.5
   end subroutine compute_k

   subroutine compute_value(g, u, v, k, value)
     double precision, intent(in) :: g, u, v, k
     double precision, intent(out) :: value
     if (k**2 * (u**2 + v**2)**0.5 /= 0d0) then
         value = (g * u) / (k**2 * (u**2 + v**2)**0.5)
     else
         value = 0d0
     end if
   end subroutine compute_value

   subroutine write_to_csv(filename, u, values)
     character(len=*), intent(in) :: filename
     double precision, dimension(:), intent(in) :: u, values
     integer :: i, n
     open(unit=10, file=filename, status='replace')
     write(10, '(a,a)') 'u_values', ',compute_value_values'
     n = size(u)
     do i = 1, n
        write(10, '(F8.5,",",F8.5)') u(i), values(i)
     end do
     close(10)
   end subroutine write_to_csv

 end program vegetation_resistance
