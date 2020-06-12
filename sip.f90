!*****************************************************************
!*Strongly Implicite procedure
!*****************************************************************

subroutine fields_comp_SIP (Ra, Pr, Bi, ampl, omega, hx, hy, psi, T,T1,T2, fi,fi1,fi2, time_in, time_out, time, eps_poisson, Max_counter_poisson, relativ, relax, ampl_psi, psi_aver)
                           
!---------------------variable declaration(beg)---------------------------
implicit none
integer, parameter :: Nx=50    ! must be even (четным)
integer, parameter :: Ny=100   ! must be even (четным)
integer i, j, k

real(kind(1d0)) Ra, Pr, Bi, k1, k2, ampl, omega, ampl_psi           ! manag. parameters
real(kind(1d0)) psi_new(0:Nx,0:Ny), psi(0:Nx,0:Ny), T(0:Nx,0:Ny), T1(0:Nx,0:Ny), T2(0:Nx,0:Ny), fi(0:Nx,0:Ny), fi1(0:Nx,0:Ny), fi2(0:Nx,0:Ny), psi_aver(0:Nx,0:Ny)           ! arrays

real(kind(1d0)) time_in, time_out, hx, hy, dt, time, relativ, relax, max_psi, max_psi_old, max_eta, max_eta_old, eps_threshold, v(0:Nx,0:Ny)
real(kind(1d0)) alpha /0.92/
real(kind(1d0)) delta_poisson, sum_poisson, eps_poisson   ! for Poisson equation
integer icounter_poisson, Max_counter_poisson, out_iter_counter, inner_iter_counter 
logical stac_flg, comp_aver/.true./
real(kind(1d0)) Un(0:Nx,0:Ny), Ue(0:Nx,0:Ny), Lp(0:Nx,0:Ny), Ls(0:Nx,0:Ny), Lw(0:Nx,0:Ny), res(0:Nx,0:Ny)
real(kind(1d0)) Ap(0:Nx,0:Ny), Ae(0:Nx,0:Ny), Aw(0:Nx,0:Ny), As(0:Nx,0:Ny), An(0:Nx,0:Ny), Q(0:Nx,0:Ny), resN
real(kind(1d0)) eps_SIP/0.000001/, eps_stac/0.000001/, nonlinearity1/0.848826/, nonlinearity2/0.63661977/ !/0./ 
real(kind(1d0)) coeff1/1.273239544/, coeff2/0.63661977/, time_start_aver/3000./, time_stop_aver/6000./, Nu
real(kind(1d0)) pi, r, v2, v2r_ratio
!---------------------variable declaration(end)---------------------------

dt=relativ*hx*hy  ! time step
time=time_in
stac_flg=.true.

!Coefficents are derived from Galerkin procedure
pi = 4*atan(1.0)
k1 = pi**2/4.0
nonlinearity1 = 8.0/(3.0*pi)
r = 2.0*(2.0+Bi)/pi

v2 = 2.0 + 4.0 * Bi * (Bi + 4.0)/pi**2
k2 = 0.5*pi*Bi*r/v2
nonlinearity2 = 4.*(1.0 + Bi + 8.*Bi**2/(3.*pi**2))/pi/v2
v2r_ratio = r/v2

!r = 1.0
!v2 = 1.0
!k2 = k1
!v2r_ratio = 1.0
!nonlinearity2 = 8.0/(3.0*pi)
!period_counter = 0

! for linear analysis stability problem
!nonlinearity1 = 0.
!nonlinearity2 = 0.

do i=0, Nx
  Un(i,0)=0.; Un(i,Ny)=0.; 
  Ue(i,0)=0.; Ue(i,Ny)=0.;
  Lp(i,0)=0.; Lp(i,Ny)=0.;
  Lw(i,0)=0.; Lw(i,Ny)=0.;
  Ls(i,0)=0.; Ls(i,Ny)=0.;
  res(i,0)=0.; res(i,Ny)=0.; 
enddo
do j=0, Ny
  Un(0,j)=0.;  Un(Nx,j)=0.; 
  Ue(0,j)=0.;  Ue(Nx,j)=0.; 
  Lp(0,j)=0.;  Lp(Nx,j)=0.; 
  Lw(0,j)=0.;  Lw(Nx,j)=0.; 
  Ls(0,j)=0.;  Ls(Nx,j)=0.; 
  res(0,j)=0.;  res(Nx,j)=0.; 
enddo

max_psi_old=0.
max_eta_old=0.
stac_flg=.true.

!print*, 'in the SIP'
!pause

!--------------------------Main cycle on time(beg)-------------------------
do while (time<time_out ) ! .and. stac_flg
time=time+dt
!
!out_iter_counter=0.
!do while (out_iter_counter<20) 
!out_iter_counter=out_iter_counter+1
!!write (*,'(I3, x, F6.3, x, F6.3, x, F6.3, x, I3)') out_iter_counter, psi(20,20), T(20,20), eta(20,20), inner_iter_counter 

!------------------------Vorticity computing(beg)--------------------
 
!---------------moving time layers (beg)----------------------
 do j=0, Ny
  do i=0, Nx
   fi2(i,j)=fi1(i,j)
   fi1(i,j)=fi(i,j)
  enddo
 enddo
!---------------moving time layers (end)----------------------
 
 do j=1, Ny-1
  do i=1, Nx-1
   Ap(i,j) = 1.5/dt + 2./(hx**2) + 2./(hy**2)  + k1
   As(i,j) = -nonlinearity1*(psi(i+1,j)-psi(i-1,j))/(Pr*4.*hx*hy) - 1./(hy**2)
   An(i,j) =  nonlinearity1*(psi(i+1,j)-psi(i-1,j))/(Pr*4.*hx*hy) - 1./(hy**2)
   Ae(i,j) = -nonlinearity1*(psi(i,j+1)-psi(i,j-1))/(Pr*4.*hx*hy) - 1./(hx**2)
   Aw(i,j) =  nonlinearity1*(psi(i,j+1)-psi(i,j-1))/(Pr*4.*hx*hy) - 1./(hx**2)
   Q(i,j) = r*(-Ra+ampl*cos(omega*time))*(T(i+1,j)-T(i-1,j))/(2.*hx)+ 2.*fi1(i,j)/dt - fi2(i,j)/(2.*dt) !+ampl*cos(omega*time)*((T(i,j+1)-T(i,j-1))/(2.*hy)-1.2732395) !vertical/horizontal
  enddo
 enddo
  
 do j=1, Ny-1
  do i=1, Nx-1
   Lw(i,j)=Aw(i,j)/(1.+alpha*Un(i-1,j))
   Ls(i,j)=As(i,j)/(1.+alpha*Ue(i,j-1))
   Lp(i,j)=Ap(i,j) + alpha*(Lw(i,j)*Un(i-1,j) + Ls(i,j)*Ue(i,j-1)) - Lw(i,j)*Ue(i-1,j) - Ls(i,j)*Un(i,j-1)
   Un(i,j)=(An(i,j) - alpha*Lw(i,j)*Un(i-1,j))/Lp(i,j)
   Ue(i,j)=(Ae(i,j) - alpha*Ls(i,j)*Ue(i,j-1))/Lp(i,j)
  enddo
 enddo

resN = 1.
inner_iter_counter=0.

do while (resN>eps_SIP)
inner_iter_counter=inner_iter_counter+1
 
 resN=0.
 do j=1, Ny-1
  do i=1, Nx-1
   res(i,j)=Q(i,j)-Ae(i,j)*fi(i+1,j)-Aw(i,j)*fi(i-1,j)-An(i,j)*fi(i,j+1)-As(i,j)*fi(i,j-1)-Ap(i,j)*fi(i,j)
   resN=resN+ABS(res(i,j))
   enddo
 enddo
 
 do j=1, Ny-1
  do i=1, Nx-1
   res(i,j)=(res(i,j) - Lw(i,j)*res(i-1,j) - Ls(i,j)*res(i,j-1))/Lp(i,j)
  enddo
 enddo

 do j=Ny-1, 1, -1 
  do i=Nx-1, 1, -1 
   res(i,j)=res(i,j) - Un(i,j)*res(i,j+1) - Ue(i,j)*res(i+1,j)
   fi(i,j)=fi(i,j)+res(i,j)
  enddo
 enddo
enddo

!---------------Boundary conditions - the Toma formula(beg)----------
do i=0, Nx
     fi(i,0)=-(8.*psi(i,1)-psi(i,2))/(2.*hy**2)
     fi(i,Ny)=-(8.*psi(i,Ny-1)-psi(i,Ny-2))/(2.*hy**2)
enddo
do j=0, Ny
     fi(0,j)=-(8.*psi(1,j)-psi(2,j))/(2.*hx**2)
     fi(Nx,j)=-(8.*psi(Nx-1,j)-psi(Nx-2,j))/(2.*hx**2)
enddo
!---------------Boundary conditions - the Toma formula(end)----------

!---------------Boundary conditions - zero vorticity (beg)----------
!do i=0, Nx
!     fi(i,0)=0.
!     fi(i,Ny)=0.
!enddo
!do j=0, Ny
!     fi(0,j)=0.
!     fi(Nx,j)=0.
!enddo
!---------------Boundary conditions - zero vorticity (end)----------

!------------------------Vorticity computing(end)--------------------

!------------------------Temperature computing(beg)------------------------
 
!---------------moving time layers (beg)---------------------- 
 do j=0, Ny
  do i=0, Nx
   T2(i,j)=T1(i,j)
   T1(i,j)=T(i,j)
  enddo
 enddo
!---------------moving time layers (end)----------------------
 
 do j=1, Ny-1
  do i=1, Nx-1
   Ap(i,j) = Pr*1.5/dt + 2./(hx**2) + 2./(hy**2) + k2
   As(i,j) = -nonlinearity2*(psi(i+1,j)-psi(i-1,j))/(4.*hx*hy) - 1./(hy**2)
   An(i,j) =  nonlinearity2*(psi(i+1,j)-psi(i-1,j))/(4.*hx*hy) - 1./(hy**2)
   Ae(i,j) = -nonlinearity2*(psi(i,j+1)-psi(i,j-1))/(4.*hx*hy) - 1./(hx**2)
   Aw(i,j) =  nonlinearity2*(psi(i,j+1)-psi(i,j-1))/(4.*hx*hy) - 1./(hx**2)
   Q(i,j) = v2r_ratio*(psi(i+1,j)-psi(i-1,j))/(2.*hx) + 2.*Pr*T1(i,j)/dt - Pr*T2(i,j)/(2.*dt)
  enddo
 enddo
 
 do j=1, Ny-1
  do i=1, Nx-1
   Lw(i,j)=Aw(i,j)/(1.+alpha*Un(i-1,j))
   Ls(i,j)=As(i,j)/(1.+alpha*Ue(i,j-1))
   Lp(i,j)=Ap(i,j) + alpha*(Lw(i,j)*Un(i-1,j) + Ls(i,j)*Ue(i,j-1)) - Lw(i,j)*Ue(i-1,j) - Ls(i,j)*Un(i,j-1)
   Un(i,j)=(An(i,j) - alpha*Lw(i,j)*Un(i-1,j))/Lp(i,j)
   Ue(i,j)=(Ae(i,j) - alpha*Ls(i,j)*Ue(i,j-1))/Lp(i,j)
  enddo
 enddo

resN = 1.
inner_iter_counter=0.

do while (resN>eps_SIP)
inner_iter_counter=inner_iter_counter+1
 
 resN=0.
 do j=1, Ny-1
  do i=1, Nx-1
   res(i,j)=Q(i,j)-Ae(i,j)*T(i+1,j)-Aw(i,j)*T(i-1,j)-An(i,j)*T(i,j+1)-As(i,j)*T(i,j-1)-Ap(i,j)*T(i,j)

   resN=resN+ABS(res(i,j))
   enddo
 enddo
 

 do j=1, Ny-1
  do i=1, Nx-1
   res(i,j)=(res(i,j) - Lw(i,j)*res(i-1,j) - Ls(i,j)*res(i,j-1))/Lp(i,j)
  enddo
 enddo


 do j=Ny-1, 1, -1 
  do i=Nx-1, 1, -1 
   res(i,j)=res(i,j) - Un(i,j)*res(i,j+1) - Ue(i,j)*res(i+1,j)
   T(i,j)=T(i,j)+res(i,j)
  enddo
 enddo
enddo

!---------------Boundary conditions - zero derivative (beg)----------
do j=0, Ny
 T(0,j) = (4.*T(1,j)-T(2,j))/3. !1.0 - j*1.0/Ny 
 T(Nx,j)= (4.*T(Nx-1,j)-T(Nx-2,j))/3. !1.0 - j*1.0/Ny 
enddo
!---------------Boundary conditions - zero derivative (end)----------

!------------------------Temperature computing(end)------------------------

!----------------------Stream function computing(beg)----------------------
icounter_poisson=0                   

do							 
delta_poisson=0.     
sum_poisson=0.              

 do i=1, Nx-1
  do j=1, Ny-1
  
   psi_new(i,j)=(1.-relax)*psi(i,j)+relax*(psi_new(i-1,j)+psi(i+1,j)+(hx/hy)**2*(psi_new(i,j-1)+psi(i,j+1)) + hx**2*fi(i,j))/(2.*(1.+(hx/hy)**2))
 
   delta_poisson=delta_poisson+abs(psi_new(i,j)-psi(i,j))
   sum_poisson=sum_poisson+abs(psi_new(i,j))

  enddo
 enddo
 
 do i=1, Nx-1
  do j=1, Ny-1
   psi(i,j)=psi_new(i,j)     
  enddo
 enddo
 
 do i= Nx-1, 1, -1  
  do j= Ny-1, 1, -1
   
   psi_new(i,j)=(1.-relax)*psi(i,j)+relax*(psi(i-1,j)+psi_new(i+1,j)+(hx/hy)**2*(psi(i,j-1)+psi_new(i,j+1)) + hx**2*fi(i,j))/(2.*(1.+(hx/hy)**2))
 
   delta_poisson=delta_poisson+abs(psi_new(i,j)-psi(i,j))
   sum_poisson=sum_poisson+abs(psi_new(i,j))

  enddo
 enddo

 icounter_poisson=icounter_poisson+1
  
 do i=1, Nx-1
  do j=1, Ny-1
   psi(i,j)=psi_new(i,j)     
  enddo
 enddo

 if ((delta_poisson.LE.sum_poisson*eps_poisson).or.(icounter_poisson.GT.Max_counter_poisson)) exit

!write(6,*) icounter_poisson, delta_poisson, sum_poisson
enddo
!write(*,*) icounter_poisson, delta_poisson, sum_poisson
!----------------------Stream function computing(end)----------------------

!enddo !out_iter_counter

!---------------Finding maximum value of stream function(beg)--------------
max_psi=abs(psi(1,1))

do i=1, Nx-1
 do j=1, Ny-1
  if (abs(psi(i,j))>max_psi) max_psi=abs(psi(i,j))
 enddo
enddo
!---------------Finding maximum value of stream function(end)--------------

if (abs(max_psi_old-max_psi)<eps_stac) stac_flg=.false.
max_psi_old=max_psi

!-----------------Nusselt number calculation(beg)--------------------------
Nu=0
do i=1,Nx-1
 Nu=Nu+(hx/hy)*(T(i,2)-4.*T(i,1)+3.*T(i,0))
enddo
!-----------------Nusselt number calculation(end)--------------------------

write(*,'(x,A,f12.4,x,A,f12.5,x,A,i4,x,A,f8.4)') 'time=', time, 'psi_max=', max_psi , 'iter=', icounter_poisson, 'Nu=', Nu
write(5,'(x,f12.4,x,f10.5,x,f10.5,x,f10.5,x,f10.5)') time, max_psi, Nu, T(int(5./hx),int(10./hy)), T(int(10./hx),int(10./hy))

!---------------Constructing avarege velocity field(beg)-------------------
if ((time>time_start_aver).and.(time<time_stop_aver)) then
 do i=0, Nx
  do j=0, Ny
  psi_aver(i,j)=psi_aver(i,j)+psi(i,j)
  enddo
 enddo
endif

if ((time>time_stop_aver).and.comp_aver) then
 comp_aver=.false.
 do i=0, Nx
  do j=0, Ny
   psi_aver(i,j)=dt*psi_aver(i,j)/(time_stop_aver-time_start_aver)
  enddo
 enddo

 open(13,file='out_aver.dat')
 do i=1, Nx-1
  do j=1, Ny-1
   v(i,j)=sqrt((psi_aver(i+1,j)-psi_aver(i-1,j)/(2*hx))**2+(psi_aver(i,j+1)-psi_aver(i,j-1)/(2*hy))**2)
   write(13,'(x,f10.5,x,f10.5,x,f10.5,x,f10.5)') i*hx,  j*hy, psi_aver(i+1,j), v(i,j)
  enddo
 enddo 
  
 open(14,file='out_velocity_inf.dat')
  do i=0, 8
   write(14,*) i*sqrt(hx**2+hy**2), v(13+i,15+i)
  enddo
endif
!---------------Constructing avarege velocity field(end)-------------------

enddo 
!--------------------------Main cycle on time(end)-------------------------
ampl_psi=max_psi

end subroutine
