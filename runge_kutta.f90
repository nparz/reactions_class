program runge_kutta
implicit none 

real(8),allocatable,dimension(:,:) :: U,Upr
real(8) :: dx,x,phase
real(8) ::energy,RL,k,int1,int2,int3,phases(3)
complex(8) :: deltaL,SL
real(8),parameter :: twomuoverh2 = .0478450d0,pi = 3.14159265359d0
integer :: gridpoints,i,L,Lmax

print*, 'how many grid points? '
read*, gridpoints
print*, 'grid spacing? '
read*, dx
print*, 'How many partial waves? ' 
read*, Lmax

! allocate storage for the wavefunctions
allocate(U(gridpoints,Lmax),Upr(gridpoints,Lmax)) 



open(unit=42,file='phase_v_energy.dat') 


energy = .1d0
do while ( energy < 4.005d0) 
   
   print*, energy
  
   k = sqrt( twomuoverh2 * energy )  ! wave vector

   do L = 0,Lmax-1 
      x = 0.d0 
      call runge_kutta_integrate(gridpoints,x,dx,U(:,L+1),Upr(:,L+1),0.d0,1.d0,L,energy) 
   
      ! R matrix
      Rl = 1.d0/x*U(gridpoints,L+1)/Upr(gridpoints,L+1) 
 
      ! intermediates
      int1 = 1.d0 - (k*x*Rl)**2
      int2 = 2.d0*k*x
      int3 = 1.d0 + (k*x*Rl)**2
   
      ! calculate complex S matrix, phase
      SL = (-1)**L / int3 * complex( int1*Cos(int2) + int2*Rl*Sin(int2) ,  int2*Rl*Cos(int2) - int1*Sin(int2) )
      deltaL = Log( (-1)**L * complex( int1*cos(int2) + int2*Rl*Sin(int2), int2*Rl*Cos(int2) - int1*Sin(int2) ) / int3 )

      phase = aimag(deltaL)/2.d0

      !regularize the phase
      do 
         if ( phase > 0.d0 ) then 
            
            if (phase > pi) then
               phase = phase -pi 
            else
               exit
            end if
      
         else 
            if ( abs(phase) < 1e-2 ) exit
            phase = phase + pi
      
         end if
      end do
   
      ! store the phase in a vector for each L
      phases(L+1) = phase



   end do

   write(42,*) energy, phases


! write the wave functions if you want
   if (abs(energy - 0.1d0) < 1e-6) then 
      
      x = 0.d0
      open(unit=45,file='wfs_e01.dat') 
      do i = 1, gridpoints
         write(45,*) x,U(i,1),U(i,2),U(i,3)
         x = x + dx
      end do
      close(45)
      
   else if ( abs(energy - 3.d0) < 1e-6 ) then 
      
      x = 0.d0
      open(unit=45,file='wfs_e30.dat') 
      do i = 1, gridpoints
         write(45,*) x,U(i,1),U(i,2),U(i,3)
         x = x + dx
      end do
      close(45)
      
   end if 

   
   energy = energy + .01d0 
   
end do


close(42)

end program
!========================================================================
!========================================================================
subroutine runge_kutta_integrate(gridpoints,t0,dt,x,v,x0,v0,L,energy) 
  implicit none 

  integer :: L,gridpoints,i  
  real(8),dimension(gridpoints) :: x,v
  real(8) :: t0,dt,t,energy
  real(8),intent(in) :: x0,v0
  real(8) :: xvd(2),k1x,k2x,k3x,k4x
  real(8) :: k1v,k2v,k3v,k4v

  ! initial conditions

  x(1) = x0
  v(1) = v0   
  t = t0+dt

do i = 2, gridpoints 
   
   call derivatives(t,x(i-1),v(i-1),L,energy,xvd) ! get derivatives
   k1x = xvd(1)
   k1v = xvd(2) 
   
   call derivatives(t+dt/2.d0,x(i-1)+k1x*dt/2.d0,v(i-1)+k1v*dt/2.d0,L,energy,xvd) ! get derivatives
   k2x = xvd(1)
   k2v = xvd(2) 

   call derivatives(t+dt/2.d0,x(i-1)+k2x*dt/2.d0,v(i-1)+k2v*dt/2.d0,L,energy,xvd) ! get derivatives
   k3x = xvd(1)
   k3v = xvd(2) 

   call derivatives(t+dt,x(i-1)+k3x*dt,v(i-1)+k3v*dt,L,energy,xvd) ! get derivatives
   k4x = xvd(1)
   k4v = xvd(2) 

   x(i) = x(i-1) + dt/6.d0 * (k1x + 2*k2x + 2*k3x + k4x ) 
   v(i) = v(i-1) + dt/6.d0 * (k1v + 2*k2v + 2*k3v + k4v )    
   
   t = t + dt
end do
t0 = t - dt


end subroutine 
!========================================================================
!========================================================================
subroutine derivatives(t,x,v,L,E,xvderivs)  
  implicit none 
  
  real(8),parameter :: twomuoverh2 = .0478450d0 , hbarc = 197.32705d0 , V0 = -61.1d0, A = 10.d0 
  real(8),intent(in) :: t,x,v,E
  real(8),dimension(2) :: xvderivs
  integer ,intent(in) :: L
  
  xvderivs(1) = v !dx/dt
  xvderivs(2) = x* (dfloat(L*(L+1))/t**2 + twomuoverh2 * (  V0 /(1 + exp((t - 1.2*A**(.333333d0))/.65)  )-E))  !dv/dt
 
end subroutine 
