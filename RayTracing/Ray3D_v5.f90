!****************************************************************************
!
!  PROGRAM: Full3D - Version 5
!
!  PURPOSE:  Ray-tracing code
!
!  Major modification from the last version: 1)Code cleaned and re-estructured.
!                                            2)New module for Ray/Bin veriables.
!                                            2)Removed dependency from MATLAB script.
!                                            3)Calculates all Sources in a single run.
!                                            4)Included modifications to consider Velocity in Y and Z directions.
!                                            5)Support Parallel Processing.
!                                            6)Support Linux Compilers.
!
!  PROGRAMERS: Dr. Carlos R. Ilário da Silva (carlos.ilario@yahoo.com)
!              M.Sc. Pedro Ricardo C. Souza (eng.pedro.ricardo@gmail.com)
!  Date: 2015
!
!
!  COMPILER COMMAND: mpif90 -O3 -g -fno-backtrace -Wall FILE_NAME.f90 -o EXE_NAME
!  RUN COMMAND: mpirun -n 4 ./EXE_NAME
!
!****************************************************************************

! Module globais calls the variables that are used through out the code
MODULE globais
DOUBLE PRECISION :: DIA
DOUBLE PRECISION, ALLOCATABLE :: XX1(:),YY1(:),ZZ1(:),Velxmat(:,:,:),Velymat(:,:,:),Velzmat(:,:,:),Cmat(:,:,:),RHOmat(:,:,:),&
&TKE(:,:,:),EPS(:,:,:),DU(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: xs(:),ys(:),zs(:)
INTEGER :: Nx,Ny,Nz,RHOflag
END MODULE globais

! Module convection calls the variables that are atualized for each ray position
MODULE convection
DOUBLE PRECISION :: Vel(3),OMEGA,RHO,Slown(3),C
DOUBLE PRECISION :: Vxamb,Vyamb,Vzamb,Camb,RHOamb,RHOstart
END MODULE

! Module invariant_blokh calls the variables used for the Blokh constant
MODULE invariant_blokh
DOUBLE PRECISION :: Velstart(3),Cstart,Tstart
END MODULE

MODULE RayVar
DOUBLE PRECISION, ALLOCATABLE :: rayCount(:),rayCountFF(:),RAYBINS(:,:),PHI(:),THETA(:)
INTEGER :: numBins,numRays
END MODULE

PROGRAM Version05
    !------------ Modules --------------------------!
    USE globais
    USE convection
    USE invariant_blokh
    USE RayVar
    USE mpi
    !-----------------------------------------------!
    DOUBLE PRECISION :: T_END
    DOUBLE PRECISION, ALLOCATABLE :: Y(:,:),TIMESTEPS(:)
    INTEGER NSTEP,Ray,Path,Ntot2,k,j,i,nproc,proc_id,fileID_1,ierr

call MPI_INIT ( ierr )
call MPI_COMM_RANK ( MPI_COMM_WORLD , proc_id , ierr )
call MPI_COMM_SIZE ( MPI_COMM_WORLD , nproc , ierr )
    
! Open basic input file

fileID_1=(proc_id+1)*100
OPEN(UNIT=fileID_1,FILE='Ray_input.dat',status='old',IOSTAT=IOS)

IF(IOS/=0) THEN
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  PRINT*,"ERROR #03: Input file: Open Error - Please check the Input File!"
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  call MPI_FINALIZE ( ierr )
  STOP
END IF

REWIND fileID_1 

READ(fileID_1,*)DIA ! Read jet's diamenter, velocity and the sound speed
READ(fileID_1,*)Vxamb,Vyamb,Vzamb,Camb,RHOamb ! Read the ambient conditions
READ(fileID_1,*)T_END,NSTEP ! Read time over which rays should be traced, and number of steps
READ(fileID_1,*)Ray ! Read the number of Rays to launch
READ(fileID_1,*)Path ! Read if the Rays Paths will be saved or not
READ(fileID_1,*)Nx,Ny,Nz ! Grid
    
CLOSE(fileID_1)

if (proc_id==0) then
  PRINT*,"*************************************"
  PRINT*,"    3D RayTracing - V5 (Parallel)   " 
  PRINT*,"*************************************"
  PRINT*,""
  PRINT*,"> Reading data from the CFD input"
  PRINT*,""

  PRINT*,"> Loading Sources File."
  PRINT*,""
end if


OPEN(UNIT=fileID_1,FILE='Sources_New.dat',status='old',IOSTAT=IOS)

IF(IOS/=0) THEN
 PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 PRINT*,"ERROR #03: Input file: Open Error - Please check the Input File!"
 PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 call MPI_FINALIZE ( ierr )
 STOP
END IF

REWIND fileID_1

read(fileID_1,*)Ntot2
ALLOCATE(xs(Ntot2),ys(Ntot2),zs(Ntot2))
    
do k=1,Ntot2
  read(fileID_1,*)xs(k),ys(k),zs(k)
end do

CLOSE(fileID_1)

if(proc_id==0) then
 PRINT*,"> Loading Flow-Field."
 PRINT*,""
end if

OPEN(unit=fileID_1,file='CFD_input.dat',status='old',IOSTAT=IOS)

IF(IOS/=0) THEN
 PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 PRINT*,"Error #01: Input CFD file: Open Error. Check the CFD input file and correct!!"
 PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 call MPI_FINALIZE ( ierr )
 STOP
END IF

REWIND(fileID_1)
    
ALLOCATE(XX1(Nx),YY1(Ny),ZZ1(Nz),Velxmat(Nx,Ny,Nz),Velymat(Nx,Ny,Nz),Velzmat(Nx,Ny,Nz),Cmat(Nx,Ny,Nz),RHOmat(Nx,Ny,Nz),&
&TKE(Nx,Ny,Nz),EPS(Nx,Ny,Nz),DU(Nx,Ny,Nz))

READ(fileID_1,*)!1
READ(fileID_1,*)!2
READ(fileID_1,*)!3
READ(fileID_1,*)!4
READ(fileID_1,*)!5
READ(fileID_1,*)!6
READ(fileID_1,*)!7
READ(fileID_1,*)!8
READ(fileID_1,*)!9
READ(fileID_1,*)!10
READ(fileID_1,*)!11
READ(fileID_1,*)!12
READ(fileID_1,*)!13
READ(fileID_1,*)!14
READ(fileID_1,*)!15
READ(fileID_1,*)!16

do k=1,nz
   do j=1,ny
      do i=1,nx
          read(fileID_1,*)xx1(i),yy1(j),zz1(k),Velxmat(i,j,k),Velymat(i,j,k),Velzmat(i,j,k),rhomat(i,j,k),tke(i,j,k),&
          &eps(i,j,k),cmat(i,j,k)!,du(i,j,k)
      end do
   end do
end do
    
CLOSE(fileID_1)

! ---- Load Farfield bin ----
if (proc_id==0) then
  PRINT*,"> Loading Farfield Bins File."
  PRINT*,""
end if

OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/farfieldBins10k',status='old',IOSTAT=IOS)

IF(IOS/=0) THEN
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  PRINT*,"ERROR #02: Farfield Bins file: Open Error - Please check the Farfield Bins input file!"
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  call MPI_FINALIZE ( ierr )
  STOP
END IF

REWIND fileID_1

! Read number of bins and allocate an array
READ(fileID_1,*)numBins
ALLOCATE(RAYBINS(numBins,2),rayCount(numBins),rayCountFF(numBins))

! Read in polar and azi angles respectively
DO I=1,numBins
  READ(fileID_1,*)RAYBINS(I,1),RAYBINS(I,2)
  rayCount(I) = 0.0D0
  rayCountFF(I) = 0.0D0
END DO

CLOSE(fileID_1)
! ------------------------------

! ---- Load Rays ----

if(proc_id==0) then
  PRINT*,"> Loading Ray File."
  PRINT*,""
end if

    SELECT CASE (Ray)
    
        CASE (1)
            OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/rayLaunch162rays',status='old',IOSTAT=IOS)   
        CASE (2)
            OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/rayLaunchAngles40k',status='old',IOSTAT=IOS)  
        CASE (3)
            OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/rayLaunchAngles163k',status='old',IOSTAT=IOS)
        CASE (4)
            OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/rayLaunchAngles300k',status='old',IOSTAT=IOS)
        CASE (5)
            OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/rayLaunchAngles655k',status='old',IOSTAT=IOS)
        CASE (6)
            OPEN(UNIT=fileID_1,FILE='/home/LRT_Comum_Files/rayLaunchAngles2M',status='old',IOSTAT=IOS)
            
    END SELECT

IF(IOS/=0) THEN
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  PRINT*,"ERROR #04: Ray Launch Angles File: Open Error. Please check the Ray Angle input file"
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  call MPI_FINALIZE ( ierr )
  STOP
END IF
REWIND fileID_1
READ(fileID_1,*)numRays
ALLOCATE(PHI(numRays),THETA(numRays))

DO J=1,numRays    
  READ(fileID_1,*)PHI(J),THETA(J)
end do

CLOSE(fileID_1)
! -------------------
if (proc_id==0) then
  PRINT*,"> Data ready!"
  PRINT*,""
  CALL system('mkdir ./output')
end if

!Main Loop for all Sources

ALLOCATE(Y(6,NSTEP+1),TIMESTEPS(NSTEP+1))
Y=0; TIMESTEPS=0

if(proc_id==0) then
  PRINT*,"-------------------------------------"
  Write(*,fmt="(a)",ADVANCE='NO')" Sources Running: "
end if

CALL MPI_BARRIER ( MPI_COMM_WORLD , ierr )

do k=int((proc_id*(Ntot2/nproc))+1),int((proc_id+1)*Ntot2/nproc)
   if(proc_id==0) Write(*,fmt="(I4,a)",ADVANCE='NO')k," "
   CALL FULLRAYTRC(k,T_END,NSTEP,Y,TIMESTEPS,Ray,Path,proc_id)
end do

if (proc_id==0) then
  print*,""
  print*,""
  print*,"Waiting for other processes to finish ..."
end if

CALL MPI_BARRIER ( MPI_COMM_WORLD , ierr )
if (proc_id==0) then
  print*,""
  print*," ***** Ray-Tracing Finished *****"
end if
call MPI_FINALIZE ( ierr )
    
end PROGRAM Version05

!***************************** Start the main code *********************************!
SUBROUTINE FULLRAYTRC(ks,T_END,NSTEP,Y,TIMESTEPS,Ray,Path,proc_id)
!
! Ray Tracer 3D is a shell code which calls RAY_TRACE
!
   
   !------------ Modules --------------------------!
    USE globais
    USE convection
    USE invariant_blokh
    USE RayVar
    !-----------------------------------------------!
    IMPLICIT NONE

    DOUBLE PRECISION :: YSTART(6),BLOKH
    DOUBLE PRECISION :: R,THETA_END,T_END,PHI_END,FF(3),TEMP,T2,T4,rayDist
    DOUBLE PRECISION :: Y(6,NSTEP+1),TIMESTEPS(NSTEP+1)
    INTEGER :: ks,NSTEP,K,IOS,I,J,Ray,Path,Ntot,Ttot,proc_id,fileID_1,fileID_2,ierr

    character(4) ::source_num
    
    !-------Modification done in September 2011-----!
    character(200),dimension(10)   :: cases
    !-----------------------------------------------
    
    CALL CPU_TIME(T2)

   !-------Modification done in September 2011-----!
    SELECT CASE (Ray)
        CASE (1)
            cases(4) ='rayLaunch162rays'
        CASE (2)
            cases(4) ='rayLaunchAngles40k'
        CASE (3)
            cases(4) ='rayLaunchAngles163k'
        CASE (4)
            cases(4) ='rayLaunchAngles300k'
        CASE (5)
            cases(4) ='rayLaunchAngles655k'
        CASE (6)
            cases(4) ='rayLaunchAngles2M'
    END SELECT

    !-----------------------------------------------!
     
!Total grid points
Ntot=Nx*Ny*Nz 

!Define source Number and Location
WRITE(source_num,fmt="(I4.4)")ks
YSTART(1)=xs(ks)
YSTART(2)=ys(ks)
YSTART(3)=zs(ks)

!Open output file, which we will write rays paths to, if switched the input file has 1 otherwise is not going to print.	

fileID_2=(proc_id+1)*101
SELECT CASE (Path)
    CASE (1) 
!       PRINT*,""
!       PRINT*,"Note: The rays paths will be saved!"
!       PRINT*,""
        OPEN(UNIT=fileID_2,FILE='./output/RaysPath_out_S'//trim(source_num),IOSTAT=IOS)
            If (Ray > 1) then
              if (proc_id==0) then
                PRINT*,""
                PRINT*,"Note #01: RaysPath output file can be extremely big!"
              end if
            end if
            IF(IOS/=0) THEN
              PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              PRINT*,"ERROR #05: Impossible to create output file."
              PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              call MPI_FINALIZE ( ierr )
              STOP
            END IF  
        REWIND fileID_2
        WRITE(fileID_2,*)numRays
    
    CASE (0)
 !      PRINT*,""
 !      PRINT*,"Note: Rays paths will not be saved!"
 !      PRINT*,""
END SELECT

if (proc_id==0) then
  PRINT*,""
  PRINT*,""
  PRINT*,"        *** Tracing rays ***        "
end if

!******************* Call flow parameters to return vel at start postion of ray ****************!
RHOflag=1
CALL CFD_Interpolation(YSTART(1:3),YSTART(1:3))
RHOflag=2
        
!***** Copying the initial values to be used in the invariant of Blokhintzev calculation *******!
do j=1,3                    
  Velstart(j)=Vel(j)
end do
Cstart=C
    
    !************************************* Starting the main LOOP *********************************!
DO J=1,numRays

  if (proc_id==0) then
    if (j==(NINT(0.1*numRays+1))) then
        PRINT*,"Progress: 10% of rays already traced!"
    elseif(j==(NINT(0.25*numRays+1))) then
        PRINT*,"Progress: 25% of Rays already traced!"
    elseif (j==(NINT(0.5*numRays+1))) then
        PRINT*,"Progress: 50% of Rays already traced!"
    elseif (j==(NINT(0.75*numRays+1))) then
        PRINT*,"Progress: 75% of Rays already traced!"
    end if
  end if

 ! Loading the initial values  
 do k=1,3                    
   Vel(k)=Velstart(k)
 end do
 C=Cstart

 CALL LAUNCH_VECTOR(THETA(J),PHI(J),YSTART(4:6))   
 CALL RAY_TRACE(YSTART,T_END,Y,TIMESTEPS,NSTEP)         

 !Loop to write full ray path to main file 'outputRays'
 SELECT CASE (Path)
    CASE (1)
        DO K=1,NSTEP
            WRITE(fileID_2,*)Y(1:3,k)
        END DO
        WRITE(fileID_2,*)""   
 END SELECT

 rayDist=( Y(2,NSTEP)**2 + Y(3,NSTEP)**2)**0.5

 IF (rayDist > 0.8) THEN
   ! Takes last two ray positions and propagates as a straight line to farfield position FF
   CALL FARFIELD(Y(1:3,NSTEP-1),Y(1:3,NSTEP),FF)
   ! Calculates farfield polar angles to the jet nozzle (0,0,0)
   CALL END_ANGLES(R,PHI_END,THETA_END,FF)
   ! Calculates the Blokhintzev invariant along the ray
   CALL BLOKHINTZEV(BLOKH,YSTART(4:6))
   ! Sum rays in ray bins with the jet switched on
   CALL RAYSUM(rayCount,RAYBINS,numBins,PHI_END,THETA_END,BLOKH)
 END IF
             
 ! Sum free field rays in bins - i.e. jet switched off       
 CALL RAYSUMFF(rayCountFF,RAYBINS,numBins,PHI(J),THETA(J))                        
END DO

if (proc_id==0) then
  PRINT*,"Progress: 100% of Rays already traced!"
end if

CLOSE(fileID_2)

!****************** Open output file, which we will write ray bin count to ******************!

fileID_1=(proc_id+1)*100
OPEN(UNIT=fileID_1,FILE='./output/DeltaSPL_out_S'//trim(source_num),IOSTAT=IOS)

IF(IOS/=0) THEN
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  PRINT*,"ERROR #05: Impossible to create output file."
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  call MPI_FINALIZE ( ierr )
  STOP
END IF

! Write ray bins total to file
WRITE(fileID_1,*)numBins

DO I=1,numBins
  IF(rayCount(I)==0 .OR. rayCountFF(I)==0)THEN
    WRITE(fileID_1,*)-200.0
  ELSE
    TEMP=rayCount(I)/rayCountFF(I)
    WRITE(fileID_1,*)10*LOG10(TEMP)
  END IF
END DO
   
CLOSE(fileID_1)
  
CALL CPU_TIME(T4)
Ttot=NINT(T4-T2)


if (proc_id==0) then
  WRITE(*,fmt="(a,I4)")" CPU time in seconds is ",Ttot
end if

if (proc_id==0) then
  print*,""
  print*,""
  PRINT*,"-------------------------------------"
  Write(*,fmt="(a)",ADVANCE='NO')" Sources Running: "
end if

return
    
END SUBROUTINE FULLRAYTRC

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE RAY_TRACE(YSTART,T_END,Y,TIMESTEPS,NSTEP)        
USE convection
! This routine traces a single ray over a fixed time
! through the medium defined by FLOW_PARAM.
!	
!	
!	
!   VEL = { v_x v_y v_z }
!	
!   X = { x y z }
!	
!	S = { s_x s_y s_z }
!	
! 	OMEGA = 1 - V . S
!	
! 	RHO is density
! 		
!   C is speed of sound
!	
!	Velocity component partial derivatives 
!	
!	DV = 	 { dv_x 	dv_x 	dv_x }
!		 { ----		---- 	---- }
!		 {  dx 		 dy  	 dz  }
!		 {			     }
!		 { dv_y 	dv_y 	dv_y }
!		 { ----		---- 	---- }
!		 {  dx 		 dy  	 dz  }
!		 {			     }
!		 { dv_z 	dv_z 	dv_z }
!		 { ----		---- 	---- }
!		 {  dx 		 dy  	 dz  }
!	
!		
!  Sound speed partial derivatives	
!	
!	DC = {  dc 	   dc 	  dc  }
!	     { ----	  ----   ---- }
!	     {  dx 	   dy  	  dz  }
!	
!---------------------------------------------------------------------------

DOUBLE PRECISION, INTENT(IN) :: T_END,YSTART(6)
    INTEGER, INTENT(IN) :: NSTEP
    DOUBLE PRECISION, INTENT(INOUT) :: Y(6,NSTEP+1),TIMESTEPS(NSTEP+1)
    
    DOUBLE PRECISION :: DC(3),DV(3,3),T0
    INTEGER :: nvar

    nvar=6
    T0=0

OMEGA=1-(Vel(1)*YSTART(4)+Vel(2)*YSTART(5)+Vel(3)*YSTART(6))
    
    CALL DERIVS_NUMERIC(YSTART,DC,DV)

CALL RKDUMB(YSTART,nvar,T0,T_END,NSTEP,Y,TIMESTEPS)


END SUBROUTINE RAY_TRACE

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE RKDUMB(vstart,nvar,x1,x2,nstep,y,xx) 

!	Fourth order Runge Kutta method, based on numerical recipes in fortran algorithm


INTEGER nstep,nvar

!PARAMETER(NMAX=6,NSTPMX=1000) 								! Maximum number of functions and
DOUBLE PRECISION x1,x2,vstart(nvar),xx(nstep+1),y(6,nstep+1) ! maximum number of values to be stored


!COMMON/path/xx,y            ! Storage ofresults. 

!	USES rk4
!	Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
!	advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
!	evaluates derivatives. Results are stored in the common block path.Be sure to dimension
!	the common block appropriately.


INTEGER i,k 
DOUBLE PRECISION h,x,dv(6),v(6) 

  
do  i=1,nvar  ! Load starting values.
v(i)=vstart(i) 
y(i,1)=v(i) 
end do  

xx(1)=x1 
x=x1 
h=(x2-x1)/ nstep 

do  k=1,nstep ! Take nstep steps. 
call derivs(v,dv) 
 call rk4(v,dv,nvar,x,h,v) 
if(x+h.eq.x)then
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  PRINT*,"ERROR: Step size not significant in  rkdumb"
  PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  call MPI_FINALIZE ( ierr )
  STOP
  end if
x=x+h 
!print*,"Max=",nstep+1
!print*,"k=",k,"i=",i
xx(k+1)=x    ! Store intermediate steps. 

do  i=1,nvar 
y(i,k+1)=v(i) 
end do
      
end do

return 

END SUBROUTINE rkdumb

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE rk4(y,dydx,n,x,h,yout) 

INTEGER n,NMAX 

DOUBLE PRECISION h,x,dydx(n),y(n),yout(n) 


PARAMETER(NMAX=50)! Set to the maximum number of functions. 

!	Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
! 	the fourth-order Runge-Kutta method to advance the solution over an interval h and return
!	the incremented variables as yout(1:n), which need not be a distinct array from y. The
!	user supplies the subroutine derivs(x,y,dydx) ,which returns derivatives dydx at x.
 

INTEGER i 

DOUBLE PRECISION h6,hh,xh,dym(6),dyt(6),yt(6) 
hh=h*0.5 
h6=h/6. 
xh=x+hh 


do i=1,n ! Firststep. 
yt(i)=y(i)+hh*dydx(i) 
enddo

call derivs(yt,dyt)  

do i=1,n ! Secondstep.

yt(i)=y(i)+hh*dyt(i) 

enddo

call derivs(yt,dym) 

do i=1,n ! Thirdstep. 
yt(i)=y(i)+h*dym(i) 
dym(i)=dyt(i)+dym(i)
enddo

call derivs(yt,dyt) ! Fourth step. 

do i=1,n ! Accumulate increments with proper weights. 
yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i)) 
enddo

return 

END SUBROUTINE rk4

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE derivs(y,dy)
Use convection
!	This routine contains the system of ODE's which rk4 call to integrate over.
!
!	t is time
!	y is intial conditions vector
!	dy is pde eqn vector
!
DOUBLE PRECISION, INTENT(IN) :: y(6)
DOUBLE PRECISION, INTENT(OUT) :: dy(6)

DOUBLE PRECISION :: DC(3),DV(3,3)

    CALL CFD_Interpolation(y(1:3),y(4:6))
    
    CALL DERIVS_NUMERIC(y(1:3),DC,DV)

    dy(1)= vel(1) + y(4)*(c**2)/OMEGA
    dy(2)= vel(2) + y(5)*(c**2)/OMEGA
    dy(3)= vel(3) + y(6)*(c**2)/OMEGA
    
    dy(4)= - (OMEGA/C)*DC(1) - (y(4)*DV(1,1) + y(5)*DV(2,1) + y(6)*DV(3,1))
    dy(5)= - (OMEGA/C)*DC(2) - (y(4)*DV(1,2) + y(5)*DV(2,2) + y(6)*DV(3,2))
    dy(6)= - (OMEGA/C)*DC(3) - (y(4)*DV(1,3) + y(5)*DV(2,3) + y(6)*DV(3,3))  

END SUBROUTINE derivs

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE LAUNCH_VECTOR(THETA,PHI,S)
USE convection
!
! Launch vector return the slowness vector of a ray launched at given angles
! 
! It solves the eional equation to determine the slowness vector for a known 
! wavefront normal vector at the initial position of the ray.
!

    
DOUBLE PRECISION, INTENT(INOUT) :: S(3),PHI
    DOUBLE PRECISION, INTENT(IN) :: THETA
    DOUBLE PRECISION :: A,B,C1,A1,A2,A3,TT,TT2,TP,TP2,C2,MINUS,PI
   
   PI=3.141592653589793D0

    TT=TAN(THETA)  
    TP=TAN(PHI)
    TT2=(TAN(THETA))**2
    TP2=(TAN(PHI))**2
    C2=C**2

A1= C2*(1+TT2)/TP2 + C2 + C2*TT2
    A2= -(VEL(1)**2)*(1 + TT2)/TP2 - 2*VEL(1)*VEL(2)*((1+TT2)**0.5)/TP - 2*VEL(1)*VEL(3)*((1 + TT2)**0.5)*TT/TP
    A3= -VEL(2)**2 - 2*VEL(2)*VEL(3)*TT - (VEL(3)**2)*TT2

A = A1 + A2 + A3

    B = 2*VEL(1)*((1 + TT2)**0.5)/TP + 2*VEL(2) + 2*VEL(3)*TT

    C1 = -1

  MINUS=1.0
IF (THETA>PI/2 .AND. THETA<3*PI/2) THEN    
   MINUS=-1.0      
    END IF

    S(2)=MINUS*(-B + (B**2 - 4*A*C1)**0.5 )/(2*A)

S(3)=S(2)*TT

    S(1)=((S(2)**2 + S(3)**2)**0.5)/TP
    
    Slown(1)=S(1)
    
    Slown(2)=S(2)
    
    Slown(3)=S(3)
 

END SUBROUTINE LAUNCH_VECTOR

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE END_ANGLES(R,PHI,THETA,X)

! This routine returns the farfield angles from the origin
! to the point the ray ends
!

DOUBLE PRECISION, INTENT(IN) :: X(3)
DOUBLE PRECISION, INTENT(OUT) :: THETA,PHI

DOUBLE PRECISION :: R,PI

PI=3.141592653589793D0
    R=(X(1)**2 + X(2)**2 + X(3)**2)**0.5

   PHI = ACOS(X(1)/R)
    THETA = ACOS(X(2)/(R*SIN(PHI)))

    IF ( X(3) < 0) THEN
      THETA=2*PI-THETA
    END IF
   
END SUBROUTINE END_ANGLES

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE RAYSUM(rayCount,RAYBINS,numBins,PHI_END,THETA_END,BLOKH)
!
!	This routine sums the number of rays turning up each ray bin.
!	after being traced through the jet.
!

INTEGER, INTENT(IN) :: numBins
    DOUBLE PRECISION, INTENT(IN) :: PHI_END,THETA_END,RAYBINS(numBins,2),BLOKH
    DOUBLE PRECISION, INTENT(INOUT) :: rayCount(numBins)

    INTEGER :: I
    DOUBLE PRECISION :: TEMP,SMALLEST


SMALLEST=100
N=1
    DO I=1,numBins
TEMP= ABS( RAYBINS(I,1) - PHI_END) + ABS( RAYBINS(I,2) - THETA_END)
        IF (TEMP<SMALLEST) THEN
          SMALLEST=TEMP
          N=I
        END IF
END DO
    
rayCount(N)=rayCount(N) + 1.0*BLOKH

END SUBROUTINE RAYSUM

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE RAYSUMFF(rayCount,RAYBINS,numBins,PHI_END,THETA_END)
!
!	This routine sums the number of rays turning up each ray bin in a free field
!

INTEGER, INTENT(IN) :: numBins
    DOUBLE PRECISION, INTENT(IN) :: PHI_END,THETA_END,RAYBINS(numBins,2)
    DOUBLE PRECISION, INTENT(INOUT) :: rayCount(numBins)

    INTEGER :: I
    DOUBLE PRECISION :: TEMP,SMALLEST


SMALLEST=100
N=1
    DO I=1,numBins
TEMP= ABS( RAYBINS(I,1) - PHI_END) + ABS( RAYBINS(I,2) - THETA_END)
        IF (TEMP<SMALLEST) THEN
          SMALLEST=TEMP
          N=I
        END IF
END DO
    
 
rayCount(N)=rayCount(N) + 1.0 

END SUBROUTINE RAYSUMFF

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------
    
SUBROUTINE FARFIELD(P1,P2,FF)
!
! Determines the farfield position a ray would intersect a sphere radius 1000m
! if proagated as straight line.
!

    USE globais

DOUBLE PRECISION, INTENT(IN) :: P1(3), P2(3)
    DOUBLE PRECISION, INTENT(INOUT) :: FF(3)

    DOUBLE PRECISION :: x1,y1,z1,x2,y2,z2,x3,y3,z3,r,a,b,c,u

    ! FF lines a sphere 1000m from origin
    ! method from
    ! http://local.wasp.uwa.edu.au/~pbourke/geometry/sphereline/

    x1=P1(1)
    y1=P1(2)
    z1=P1(3)

    x2=P2(1)
    y2=P2(2)
    z2=P2(3)

    ! centre of farfield sphere
    x3=0.0
    y3=0.0
    z3=0.0

    ! sphere radius     
    r=1000              

    a = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2

b = 2*( (x2 - x1)*(x1 - x3) + (y2 - y1)*(y1 - y3) + (z2 - z1)*(z1 - z3) )

c = x3**2 + y3**2 + z3**2 + x1**2 + y1**2 + z1**2 - 2*(x3*x1 + y3*y1 + z3*z1) - r**2 

    u = (-b + (b**2 - 4*a*C)**0.5)/(2*a)

    FF(1) = x1 + u*( x2 - x1 )
    FF(2) = y1 + u*( y2 - y1 )
    FF(3) = z1 + u*( z2 - z1 )
    
END SUBROUTINE FARFIELD

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE BLOKHINTZEV(BLOKH,SStart) 
! dermines the blokhintzev invariant
USE invariant_blokh
USE convection
implicit none
DOUBLE PRECISION, INTENT(IN) :: SStart(3)
    DOUBLE PRECISION, INTENT(INOUT) :: BLOKH

    DOUBLE PRECISION :: CEnd,Term,Vrend,Vrstart,Vratio,Omgstart,Vr_x,Vr_y,Vr_z

    CEnd=Camb   
          
        Omgstart=1-(VelStart(1)*SStart(1)+VelStart(2)*SStart(2)+VelStart(3)*SStart(3)) !Omgstart=1-VelStart(1)*SStart(1)
        Term=(RHOstart*(Cstart**2))/(RHOamb*(CEnd**2)) !Term for temperature effects on the Blokhintzev constant
        
        Vr_x=VelStart(1)+(SStart(1)/Omgstart)*Cstart**2
        
        
        Vr_y=VelStart(2)+(SStart(2)/Omgstart)*Cstart**2 !Vr_y=(SStart(2)/Omgstart)*Cstart**2
        Vr_z=VelStart(3)+(SStart(3)/Omgstart)*Cstart**2 !Vr_z=(SStart(3)/Omgstart)*Cstart**2
        
        Vrstart=SQRT((Vr_x**2)+(Vr_y**2)+(Vr_z**2)) !SQRT((SStart(1)*Cend**2)**2)
        Vrend=Camb !Vrstart=SQRT((VelStart(1)+((SStart(1)*Cstart**2)/Omgstart))**2)
                                    
        Vratio=Vrend/Vrstart
        BLOKH=Omgstart*Term*Vratio



END SUBROUTINE BLOKHINTZEV

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!--------------------------------------------------------------------------

SUBROUTINE DERIVS_NUMERIC(X,DC,DV)
USE convection

! This routine calculates the derivatives of specific 
! fluid parameters at the postion X.
! Velocity component partial derivatives 
!
!	DV = { dv_x 	dv_x 	dv_x }
!		 { ----		---- 	---- }
!		 {  dx 		 dy  	 dz  }
!		 {			     }	
!		 { dv_y 	dv_y 	dv_y }
!		 { ----		---- 	---- }
!		 {  dx 		 dy  	 dz  }
!		 {			     }	
!		 { dv_z 	dv_z 	dv_z }
!		 { ----		---- 	---- }
!		 {  dx 		 dy  	 dz  }
!	
!		
!  Sound speed partial derivatives	
!	
!	DC = {  dc     dc 	  dc  }
!	     { ----	  ----   ---- }
!	     {  dx 	   dy  	  dz  }
!	

DOUBLE PRECISION, INTENT(IN) :: X(3)
    DOUBLE PRECISION, INTENT(OUT):: DC(3),DV(3,3)

DOUBLE PRECISION :: XD(3),S(3),Ctemp,Omegatemp,DELTA
DOUBLE PRECISION :: Velx1,Velx2,Velx3
  
    ! Saving the values before marching in the directions to calculate the derivatives
    Velx1=Vel(1)     
    Velx2=Vel(2)     
    Velx3=Vel(3)     
    Ctemp=C          
    Omegatemp=OMEGA  
    ! ***************************************************
    
    S(1)=Slown(1)
    S(2)=Slown(2)
    S(3)=Slown(3)
    
    DELTA = 0.00000001
    
! Small change in the x direction
XD(1)=X(1) + DELTA
    XD(2)=X(2)
    XD(3)=X(3)
    
    CALL CFD_Interpolation(XD,S) ! Com a modificaçao de X será calculado a nova matriz de velocidade Vel(1:3) e utilizada abaixo para o cálculo da derivada.
  
  DC(1)= (C - Ctemp)/DELTA

    DV(1,1)= (Vel(1) - Velx1)/ DELTA    
DV(2,1)= (Vel(2) - Velx2)/ DELTA
    DV(3,1)= (Vel(3) - Velx3)/ DELTA


! Small change in the y direction
XD(1)=X(1) 
    XD(2)=X(2) + DELTA
    XD(3)=X(3)

    CALL CFD_Interpolation(XD,S)
  
DC(2)= (C - Ctemp)/DELTA

    DV(1,2)= (Vel(1) - Velx1)/ DELTA    
DV(2,2)= (Vel(2) - Velx2)/ DELTA
    DV(3,2)= (Vel(3) - Velx3)/ DELTA
    
! Small change in the z direction
XD(1)=X(1) 
    XD(2)=X(2) 
    XD(3)=X(3) + DELTA
    
    CALL CFD_Interpolation(XD,S)   
    
    DC(3)= (C - Ctemp)/DELTA
    
DV(1,3)= (Vel(1) - Velx1)/ DELTA    
DV(2,3)= (Vel(2) - Velx2)/ DELTA
    DV(3,3)= (Vel(3) - Velx3)/ DELTA
    
    ! Retornando o valor original de Vel antes do cálculo das derivadas
    Vel(1)=Velx1    !gravando a velocidade antes do avanço em X de DELTA
    Vel(2)=Velx2    !gravando a velocidade antes do avanço em X de DELTA
    Vel(3)=Velx3    !gravando a velocidade antes do avanço em X de DELTA
    ! Retornando o valor da velocidade do som antes do cálculo das derivadas
    C=Ctemp
    ! Retornando o valor de omega antes do cálculo das derivadas
    OMEGA=Omegatemp

END SUBROUTINE DERIVS_NUMERIC

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE CFD_Interpolation(XP,S)
!
! This subroutine checks if the ray point is inside the computational domain and if it is it calls the trilinear interpolation.
!
    USE globais
    USE convection
    
    DOUBLE PRECISION :: XP(3),S(3),XtempP(3)
    
    XtempP(1)=XP(1)
    XtempP(2)=XP(2)
    XtempP(3)=XP(3)
    
   if (XP(2)>YY1(Ny) .OR. XP(2)<YY1(1)) then
            Vel(1)=Vxamb
            Vel(2)=Vyamb
            Vel(3)=Vzamb
            C=Camb
     
    elseif (XP(3)>ZZ1(Nz) .OR. XP(3)<ZZ1(1)) then
            Vel(1)=Vxamb
            Vel(2)=Vyamb
            Vel(3)=Vzamb
            C=Camb
  
    else
        if (XtempP(1)<=XX1(1))then
            XtempP(1)=XX1(1)
   
        elseif (Xtempp(1)>=XX1(Nx)) then
            XtempP(1)=XX1(Nx)
   
        else
            XtempP(1)=XP(1)
 
        end if

        call Trilinear_Interpolation(XtempP(1:3))        
    end if
        
    RHO=1
OMEGA=1-(VEL(1)*S(1)+VEL(2)*S(2)+VEL(3)*S(3))

    return
    
END SUBROUTINE CFD_Interpolation

!---------------------------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!---------------------------------------------------------------------------

SUBROUTINE Trilinear_Interpolation(XP)
!
! This subroutine calculates by a trilinear interpolation the variables on the point of interest
!
    USE Convection
    USE Globais
    
    IMPLICIT NONE
    
    INTEGER :: i,nx1,nx2,ny1,ny2,nz1,nz2,nt
    DOUBLE PRECISION :: dif,dift,XP(3)
    DOUBLE PRECISION :: v1,v2,v3,v4,v5,v6,v7,v8,c1,c2,c3,c4,c5,c6,c7,c8
    DOUBLE PRECISION :: Lxcel,Lycel,Lzcel
    DOUBLE PRECISION :: U00,U10,U01,U11,U0,U1
    
   ! ***************************** Finding the cell that surrounds the point **************************!
    dif=XP(1)-XX1(1)
    nx1=1
    do i=1,Nx
        dift=abs(XX1(i)-XP(1))
        if (dift<abs(dif)) then
            dif=dift
            nx1=i
        end if
    end do
    
    if ((XP(1)-XX1(nx1)).lt.0.0) then
        nt=nx1-1
        nx2=nx1
        nx1=nt
    else
        nx2=nx1+1
    end if
    
    dif=XP(2)-YY1(1)
    ny1=1
    do i=1,Ny
        dift=abs(YY1(i)-XP(2))
        if (dift<abs(dif)) then
            dif=dift
            ny1=i
        end if
    end do
   
    if ((XP(2)-YY1(ny1)).lt.0.0) then
        nt=ny1-1
        ny2=ny1
        ny1=nt
    else
        ny2=ny1+1
    end if
    
    dif=XP(3)-ZZ1(1)
    nz1=1
    do i=1,Nz
        dift=abs(ZZ1(i)-XP(3))
        if (dift<abs(dif)) then
            dif=dift
            nz1=i
        end if
    end do
    
    if ((XP(3)-zz1(nz1)).lt.0.0) then
        nt=nz1-1
        nz2=nz1
        nz1=nt
    else
        nz2=nz1+1
    end if
    
    ! Correction for the X direction - if outside the X domain it copies the values from the boundaries 
    If (nx1>=Nx) then      
        nx2=nx1
        Lxcel=1
    elseif (nx1<1) then
        nx1=1
        nx2=1
        Lxcel=1
    else
        Lxcel=abs(XX1(nx1)-XX1(nx2))    !Calculates the lenght of the cell in the X direction
    end if
    
    
    ! Correction for the Y direction - if the point is in the Y boundary 
    If (ny1.eq.Ny) then      
        ny2=ny1
        Lycel=1
    elseif (ny1.eq.1) then
        ny1=1
        ny2=1
        Lycel=1
    else
        Lycel=abs(YY1(ny1)-YY1(ny2))    !Calculates the lenght of the cell in the X direction
    end if
    
    ! Correction for the Z direction - if the point is in the Z boundary 
    If (nz1.eq.Nz) then      
        nz2=nz1
        Lzcel=1
    elseif (nz1.eq.1) then
        nz1=1
        nz2=1
        Lzcel=1
    else
        Lzcel=abs(ZZ1(nz1)-ZZ1(nz2))    !Calculates the lenght of the cell in the X direction
    end if
    
    
          
    !******************************** VELOCITY in the X direction ******************************!
    ! Writing the velocity at each point of the cell that contains the point
    v1=Velxmat(nx1,ny1,nz1)
    v2=Velxmat(nx1,ny1,nz2)
    v3=Velxmat(nx1,ny2,nz1)
    v4=Velxmat(nx1,ny2,nz2)
    v5=Velxmat(nx2,ny1,nz1)
    v6=Velxmat(nx2,ny1,nz2)
    v7=Velxmat(nx2,ny2,nz1)
    v8=Velxmat(nx2,ny2,nz2)
    
    ! Trilinear interpolation as described in http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/index.html
    U00=v1+(XP(3)-ZZ1(nz1))*(v2-v1)/Lzcel
    U10=v5+(XP(3)-ZZ1(nz1))*(v6-v5)/Lzcel
    U01=v3+(XP(3)-ZZ1(nz1))*(v4-v3)/Lzcel
    U11=v7+(XP(3)-ZZ1(nz1))*(v8-v7)/Lzcel
    U0=U00+(XP(1)-XX1(nx1))*(U10-U00)/Lxcel
    U1=U01+(XP(1)-XX1(nx1))*(U11-U01)/Lxcel
    
    Vel(1)=U0+(XP(2)-YY1(ny1))*(U1-U0)/Lycel !Calculates the velocity on the desired point
    
    
    !************************************ VELOCITY in the Y direction ******************************!
    ! Writing the velocity at each point of the cell that contains the point
    v1=Velymat(nx1,ny1,nz1)
    v2=Velymat(nx1,ny1,nz2)
    v3=Velymat(nx1,ny2,nz1)
    v4=Velymat(nx1,ny2,nz2)
    v5=Velymat(nx2,ny1,nz1)
    v6=Velymat(nx2,ny1,nz2)
    v7=Velymat(nx2,ny2,nz1)
    v8=Velymat(nx2,ny2,nz2)
    
    ! Trilinear interpolation as described in http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/index.html
    U00=v1+(XP(3)-ZZ1(nz1))*(v2-v1)/Lzcel
    U10=v5+(XP(3)-ZZ1(nz1))*(v6-v5)/Lzcel
    U01=v3+(XP(3)-ZZ1(nz1))*(v4-v3)/Lzcel
    U11=v7+(XP(3)-ZZ1(nz1))*(v8-v7)/Lzcel
    U0=U00+(XP(1)-XX1(nx1))*(U10-U00)/Lxcel
    U1=U01+(XP(1)-XX1(nx1))*(U11-U01)/Lxcel
    
    Vel(2)=U0+(XP(2)-YY1(ny1))*(U1-U0)/Lycel !Calculates the velocity on the desired point
        
    !************************************ VELOCITY in the Z direction ******************************!
    ! Writing the velocity at each point of the cell that contains the point
    v1=Velzmat(nx1,ny1,nz1)
    v2=Velzmat(nx1,ny1,nz2)
    v3=Velzmat(nx1,ny2,nz1)
    v4=Velzmat(nx1,ny2,nz2)
    v5=Velzmat(nx2,ny1,nz1)
    v6=Velzmat(nx2,ny1,nz2)
    v7=Velzmat(nx2,ny2,nz1)
    v8=Velzmat(nx2,ny2,nz2)
    
    ! Trilinear interpolation as described in http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/index.html
    U00=v1+(XP(3)-ZZ1(nz1))*(v2-v1)/Lzcel
    U10=v5+(XP(3)-ZZ1(nz1))*(v6-v5)/Lzcel
    U01=v3+(XP(3)-ZZ1(nz1))*(v4-v3)/Lzcel
    U11=v7+(XP(3)-ZZ1(nz1))*(v8-v7)/Lzcel
    U0=U00+(XP(1)-XX1(nx1))*(U10-U00)/Lxcel
    U1=U01+(XP(1)-XX1(nx1))*(U11-U01)/Lxcel
    
    Vel(3)=U0+(XP(2)-YY1(ny1))*(U1-U0)/Lycel !Calculates the velocity on the desired point
    
    
    !******************************** SOUND SPEED ******************************!
    ! Writing the sound speed at each point of the cell that contains the point
    c1=Cmat(nx1,ny1,nz1)
    c2=Cmat(nx1,ny1,nz2)
    c3=Cmat(nx1,ny2,nz1)
    c4=Cmat(nx1,ny2,nz2)
    c5=Cmat(nx2,ny1,nz1)
    c6=Cmat(nx2,ny1,nz2)
    c7=Cmat(nx2,ny2,nz1)
    c8=Cmat(nx2,ny2,nz2)
    
    ! Trilinear interpolation
    U00=c1+(XP(3)-ZZ1(nz1))*(c2-c1)/Lzcel
    U10=c5+(XP(3)-ZZ1(nz1))*(c6-c5)/Lzcel
    U01=c3+(XP(3)-ZZ1(nz1))*(c4-c3)/Lzcel
    U11=c7+(XP(3)-ZZ1(nz1))*(c8-c7)/Lzcel
    U0=U00+(XP(1)-XX1(nx1))*(U10-U00)/Lxcel
    U1=U01+(XP(1)-XX1(nx1))*(U11-U01)/Lxcel
    
    C=U0+(XP(2)-YY1(ny1))*(U1-U0)/Lycel !Calculates the speed of sound on the desired point
    
   !******************************** Density ******************************! 
    if (RhoFlag==1) then
        ! Writing the density value at each point of the cell that contains the point
        c1=RHOmat(nx1,ny1,nz1)
        c2=RHOmat(nx1,ny1,nz2)
        c3=RHOmat(nx1,ny2,nz1)
        c4=RHOmat(nx1,ny2,nz2)
        c5=RHOmat(nx2,ny1,nz1)
        c6=RHOmat(nx2,ny1,nz2)
        c7=RHOmat(nx2,ny2,nz1)
        c8=RHOmat(nx2,ny2,nz2)
        
        ! Trilinear interpolation
        U00=c1+(XP(3)-ZZ1(nz1))*(c2-c1)/Lzcel
        U10=c5+(XP(3)-ZZ1(nz1))*(c6-c5)/Lzcel
        U01=c3+(XP(3)-ZZ1(nz1))*(c4-c3)/Lzcel
        U11=c7+(XP(3)-ZZ1(nz1))*(c8-c7)/Lzcel
        U0=U00+(XP(1)-XX1(nx1))*(U10-U00)/Lxcel
        U1=U01+(XP(1)-XX1(nx1))*(U11-U01)/Lxcel
        
        RHOstart=U0+(XP(2)-YY1(ny1))*(U1-U0)/Lycel !Calculates the speed of sound on the desired point
    end if
    return
    
END SUBROUTINE Trilinear_Interpolation
