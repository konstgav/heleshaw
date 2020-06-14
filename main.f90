!******************************************************************************************
!*Convection in the Hele-Shaw cell with vertical/horizontal vibrations of finite frequency*
!******************************************************************************************
!Parameters:
!	  Ra - Rayleigh number
!	  Pr - Prandtl number
!     k1 - hydrodynamic damping coefficient on the walls
!     k2 - heat damping coefficient on the walls
!   ampl - vibration amplitude
!  omega - vibration frequency
!      L - cell length
!      H - cell height

program Hele_Shaw
!use GVLRG_int

!---------------------variable declaration(beg)---------------------------
implicit none
integer, parameter :: Nx=50    ! must be even (четным)
integer, parameter :: Ny=100    ! must be even (четным)
integer, parameter :: N_vec=2

real(kind(1d0)) Ra, Pr, Bi, Ra1, Ra2, L, H, ampl, omega, ampl_psi           ! manag. parameters
real(kind(1d0)) eps_poisson, eps_threshold     ! accuracy
real(kind(1d0)) psi(0:Nx,0:Ny), T(0:Nx,0:Ny), fi(0:Nx,0:Ny), T1(0:Nx,0:Ny), fi1(0:Nx,0:Ny), T2(0:Nx,0:Ny), fi2(0:Nx,0:Ny)             ! arrays

real(kind(1d0)) time_for_count, hx, hy, dt, time, relativ, relax, period, psi_aver(0:Nx,0:Ny)
real(kind(1d0)) max_psi, max_psi_new, max_psi_old , increment_new, increment, max_eta_new, freq , time_in, time_out, time_out_shift ! inc
integer i, j, k, inc_counter, time_counter, N_half
logical first_time

integer  Max_counter_poisson

real(kind(1d0)) PI/3.14159265/
real(kind(1d0)) x1, x2, z1, z2, y1, y2, y3 !, a_sq, b_sq, c_sq, x1, x2, x3, y1, y2, y3, R_sub

character(len=8) GetFileNameCounter
character(len=8) filename_counter
character(:),allocatable :: filename_psi, filename_T, filename_fi
character(:),allocatable :: parametersFileName

integer period_counter, period_Num_steps
real(kind(1d0)) time_step_period

integer vec_counter, max_mult_num       ! for eigenvalues computations
real(kind(1d0)) Scal_Prod_ful(N_vec,N_vec+1), vec(1:N_vec+1,1:3*(Nx+1)*(Ny+1)), A_vec(N_vec,N_vec), B_vec(N_vec,N_vec)
real(kind(1d0)) beta_vec(N_vec), eval(N_vec), max_mult
complex(kind(1d0)) alpha_vec(N_vec), evalc(N_vec), complx
real(kind(1d0)) delta_R/10.0/, factor/1./, time_inc_comp/10./
real(kind(1d0)) x,y
real(kind(1d0)) incrementOld, incrementNew, deltaTimeInc, espIncRelative /0.0001/, ampl_psi_old
real(kind(1d0)) deltaRa, RaMax, Ra_beg
integer caseNumber, counterBi

character(:),allocatable :: resultsDirName        !   for output
character(:),allocatable :: psiDirName, TDirName, phiDirName
character(len=2) :: month, day, hour, minutes
character(len=4) :: year
integer,dimension(8) :: dateTimeValues
!---------------------variable declaration(end)---------------------------

!------------------creating output directories(beg)-----------------------
call date_and_time(VALUES=dateTimeValues)
write(year, '(I4.4)') dateTimeValues(1)
write(month, '(I2.2)') dateTimeValues(2)
write(day, '(I2.2)') dateTimeValues(3)
write(hour, '(I2.2)') dateTimeValues(5)
write(minutes, '(I2.2)') dateTimeValues(6)
resultsDirName = '../hs_results' // year // '_' // month // '_' // day // '_' // hour // '_' // minutes
call system('mkdir ' // resultsDirName)
psiDirName = resultsDirName // '/psi'
TDirName = resultsDirName // '/T'
phiDirName = resultsDirName // '/phi'
call system('mkdir ' // psiDirName)
call system('mkdir ' // TDirName)
call system('mkdir ' // phiDirName)
call system('mkdir ' // resultsDirName // '/animation')
call system('cp -p ./makeanimation.sh ' // resultsDirName // '/makeanimation.sh')
call system('cp -p ./fieldsplot.plt ' // resultsDirName // '/fieldsplot.plt')
!------------------creating output directories(end)-----------------------

!--------------------Reading data from file(beg)--------------------------
parametersFileName = './input2.dat'
open(1, file = parametersFileName)
read(1,*) Ra
read(1,*) Pr
read(1,*) Bi
read(1,*) ampl
read(1,*) omega
read(1,*) L
read(1,*) H
read(1,*) time_for_count
read(1,*) eps_poisson
read(1,*) Max_counter_poisson
read(1,*) eps_threshold
read(1,*) relativ
read(1,*) relax
!---------------------Reading data from file(end)--------------------------

open(5, file = resultsDirName // '/out_psimax_t.dat')

hx=L/Nx		 	! grid step
hy=H/Ny		 	! grid step

!open(unit = 11, file = resultsDirName // '/Ra_psiMax.dat', &
!            action = "write", status = "replace")

!----------------------Initial and boundary conditions(beg)----------------    
do i=0, Nx
     psi(i,0)=0.
     psi(i,Ny)=0.
     T(i,0)=0.
     T(i,Ny)=0.
     fi(i,0)=0.
     fi(i,Ny)=0.
enddo
do j=0, Ny
     psi(0,j)=0.
     psi(Nx,j)=0.
     T(0,j)=0.
     T(Nx,j)=0.
     fi(0,j)=0.
     fi(Nx,j)=0.
enddo

! Initial conditions
do i=1, Nx-1
 do j=1, Ny-1
    psi(i,j) = sin(2.*PI*i*hx/L)*sin(PI*j*hy/H) !sin(caseNumber*PI*i*hx/L)*sin(PI*j*hy/H)
    fi(i,j) = 0. !sin(2.*PI*i*hx/L)*sin(PI*j*hy/H)
    T(i,j) = cos(2.*PI*i*hx/L)*sin(PI*j*hy/H) !cos(caseNumber*PI*i*hx/L)*sin(PI*j*hy/H)
 enddo
enddo

!!!wrong initializing
fi1=fi
T1=T
fi2=fi
T2=T
psi_aver=0.
!----------------------Initial and boundary conditions(end)----------------

time_in = 0.
time_out = time_for_count

! Initializing computations
call fields_comp_SIP (Ra, Pr, Bi, ampl, omega, hx, hy, psi, T,T1,T2, fi,fi1,fi2, time_in, time_out, time, eps_poisson, Max_counter_poisson, relativ, relax, ampl_psi,psi_aver)

period = 40.
period_counter = 0
time_in = time
period_Num_steps = 20
time_step_period = period/period_Num_steps

incrementOld = 1.
deltaTimeInc = time_step_period
ampl_psi_old = ampl_psi

do while(period_counter <= period_Num_steps)

!-------------------------Write data to file(beg)--------------------------
 filename_counter = GetFileNameCounter(period_counter)
 filename_psi=trim(psiDirName//'/'//filename_counter//'.dat')
 filename_T=trim(TDirName//'/'//filename_counter//'.dat')
 filename_fi=trim(phiDirName//'/'//filename_counter//'.dat')

 open(period_counter+1000,file=filename_psi)
 open(period_counter+2000,file=filename_T)
 open(period_counter+3000,file=filename_fi)

 do i=0, Nx
  do j=0, Ny
   write(period_counter+1000,'(x,f10.5,x,f10.5,x,f10.5)') i*hx,  j*hy, psi(i,j)
   write(period_counter+2000,'(x,f10.5,x,f10.5,x,f10.5)') i*hx,  j*hy, T(i,j)
   write(period_counter+3000,'(x,f10.5,x,f10.5,x,f10.5)') i*hx,  j*hy, fi(i,j)
  enddo
   write(period_counter+1000,*) 
   write(period_counter+2000,*)
   write(period_counter+3000,*)
 enddo
!-------------------------Write data to file(end)--------------------------

time_out=time_in+time_step_period
call fields_comp_SIP (Ra, Pr, Bi, ampl, omega, hx, hy, psi, T,T1,T2, fi,fi1,fi2, time_in, time_out, time, eps_poisson, Max_counter_poisson, relativ, relax, ampl_psi,psi_aver)
deltaTimeInc = time - time_in
time_in=time
period_counter=period_counter+1

enddo

print*, 'end of the program'

call system('ls')
call system('chmod 777 '//resultsDirName//'/makeanimation.sh')
call system('cd '//resultsDirName//' && '//resultsDirName//'/makeanimation.sh')

end program Hele_Shaw

!Convert an integer to string
character(len=8) function GetFileNameCounter(counter)
    integer counter
    write (GetFileNameCounter, '(I8.8)' ) counter
    GetFileNameCounter = trim(GetFileNameCounter)
end function GetFileNameCounter


! read data frrom file
!open(12,file='field_start_psi.dat')
!open(13,file='field_start_fi.dat')
!open(14,file='field_start_T.dat')
! do i=0, Nx
!  do j=0, Ny
!   read(12,*) x,  y, psi(i,j)  !'(x,f10.5,x,f10.5,x,f10.5)'
!   read(13,'(x,f10.5,x,f10.5,x,f10.5)') x,  y, fi(i,j)
!   read(14,'(x,f10.5,x,f10.5,x,f10.5)') x,  y, T(i,j)
!  enddo
! enddo