!****************************************************************************
!
!  PROGRAM: GridMaker Version 5
!
!  PURPOSE: Mesh generation application for the 3D RayTracing methodology
!
!  Major modification from the last version: 1)Input file changed.
!                                            2)Mesh stretching corrected.
!                                            3)Domain positining in Dj added.
!                                            4)Loading bar added.
!
!  PROGRAMERS: Dr. Carlos R. Ilário da Silva (carlos.ilario@yahoo.com)
!              M.Sc. Pedro Ricardo C. Souza (eng.pedro.ricardo@gmail.com)
!  Date: 2015
!
!
!  COMPILER COMMAND: gfortran -O3 -g -fno-backtrace -Wall FILE_NAME.f90 -o EXE_NAME
!  RUN COMMAND: ./EXE_NAME
!
!****************************************************************************

    module globais
        INTEGER :: NX,NY,NZ
    end module
    
    program GridMaker
    Use globais

    IMPLICIT NONE
    DOUBLE PRECISION, ALLOCATABLE :: xx(:,:),yy(:,:),zz(:,:)
    REAL :: LenX,LenY,LenZ,Dj
    REAL :: Ntemp,xx0,yy0,zz0,xx1,yy1,zz1,T1,T2
   
    INTEGER :: NY05,NZ05,i,j,k,Ntot,y2D,z2D,IOS,Ttot,Tec2D
    REAL :: XSF,YSF,ZSF,deltx,delty,deltz,dxsf,dysf,dzsf
    
    PRINT*,"**********************"
    PRINT*,"     GridMaker V5       "
    PRINT*,"**********************"
    PRINT*,""
    PRINT*,"1. Preparing to generate the grid."
    PRINT*,""
    
    CALL CPU_TIME(T1)
    
    OPEN(unit=20,file='GridMaker_input.dat',status='old',IOSTAT=IOS)
    IF(IOS/=0) THEN
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      PRINT*,"GridMaker - Error #01: Input file: Open Error. Check the input file and correct!!"
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      STOP
    END IF
    
    READ(20,*)Dj                                        !Jet Diameter  
    READ(20,*)NX,NY,NZ                                  !Number of points in the x,y,z direction, respectively
    
    ALLOCATE(xx(1,NX),yy(1,NY),zz(1,NZ))
    
    Ntot=NX*NY*NZ                                       !Total number of points in the grid
    
    NZ05=INT(NZ/2)
    NY05=INT(NY/2)
    
    PRINT*,"1.1 Total number of nodes=",Ntot
    PRINT*,""
   
    READ(20,*)xx0,yy0,zz0                               !Origin coordinates (center of the jet)  
    READ(20,*)xx1,yy1,zz1                               !First point coordinate on the mesh in Jet Dj
    xx1=Dj*xx1; yy1=Dj*yy1; zz1=Dj*zz1
    READ(20,*)LenX,LenY,LenZ                            !Lenght of the retangular domain (Lx,Ly,Lz) 
    LenX=Dj*LenX; LenY=Dj*LenY; LenZ=Dj*LenZ
    READ(20,*)XSF,YSF,ZSF                               !Stretching factor in the x,y,z direction
    READ(20,*)dxsf,dysf,dzsf
    READ(20,*)Tec2D                                     !Flag to generate the 2D grid
    
    CLOSE(20)
    
    PRINT*,"1.2 Generating mesh..."
    PRINT*,""
    
   
    !************** Calculates the space of each point in function of the number of points *****************
    deltx=LenX/(NX-1)                           
    delty=LenY/(NY-1)
    deltz=LenZ/(NZ-1)
    
    !*************** Writing the arrays for x, y, z, respectively, in function of the delta and the streching factor    
        
    if (XSF.eq.1.0d0) then !Regular uniform mesh
            xx(1,1)=xx0+xx1                                   !First x point of the domain 
            xx(1,Nx)=xx0+LenX+xx1                             !Last x point of the domain
        do i=2,Nx-1
            xx(1,i)=xx(1,i-1)+deltx
        end do
    else
        xx(1,1)=xx0+xx1                                   !First x point of the domain 
        xx(1,2)=xx(1,1)+dxsf
        xx(1,Nx)=xx0+LenX+xx1                             !Last x point of the domain
        Nx=int( log((dxsf - xx(1,Nx) + XSF*xx(1,Nx))/(XSF*dxsf))/log(XSF) + 2 )
        DEALLOCATE(xx)
        PRINT*,"---------------WARNING----------------"
        PRINT*,"Nx changed due to Streatching Factor"
        PRINT*,"New Nx value: ",Nx
        PRINT*,"--------------------------------------"
        PRINT*,""
        ALLOCATE(xx(1,Nx))
        xx(1,1)=xx0+xx1                                   
        xx(1,2)=xx(1,1)+dxsf
            
        do i=3,Nx
            xx(1,i)=xx(1,2)+(xx(1,i-1)-xx(1,1))*XSF
        end do
    end if
    if (YSF.eq.1.0d0) then
         yy(1,1)=yy0+yy1                          !First y point of the domain
         yy(1,Ny)=yy0+LenY+yy1                         !Last y point of the domain
        do i=2,Ny-1
            yy(1,i)=yy(1,i-1)+delty
        end do
    else
        yy(1,1)=yy0+yy1                                   !First x point of the domain 
        yy(1,2)=yy(1,1)+dysf
        yy(1,Ny)=yy0+LenY+yy1                             !Last x point of the domain
        Ny=int( log((dysf - yy(1,Ny) + YSF*yy(1,Ny))/(YSF*dysf))/log(YSF) + 2 )
        DEALLOCATE(yy)
        PRINT*,"---------------WARNING----------------"
        PRINT*,"Ny changed due to Streatching Factor"
        PRINT*,"New Ny value: ",Ny
        PRINT*,"--------------------------------------"
        PRINT*,""
        ALLOCATE(yy(1,Ny))
        yy(1,1)=yy0+yy1                                   
        yy(1,2)=yy(1,1)+dysf
            
        do i=3,Ny
            yy(1,i)=yy(1,2)+(yy(1,i-1)-yy(1,1))*YSF
        end do
    end if
    if (ZSF.eq.1.0d0) then
             zz(1,1)=zz0+zz1                          !First z point of the domain
             zz(1,Nz)=zz0+LenZ+zz1                         !First z point of the domain
        do i=2,Nz-1
            zz(1,i)=zz(1,i-1)+deltz
        end do
    else                                                 !Non-uniform mesh
        zz(1,1)=zz0+zz1                                   !First x point of the domain 
        zz(1,2)=zz(1,1)+dzsf
        zz(1,Nz)=zz0+LenZ+zz1                             !Last x point of the domain
        Nz=int( log((dzsf - zz(1,Nz) + ZSF*zz(1,Nz))/(ZSF*dzsf))/log(ZSF) + 2 )
        DEALLOCATE(zz)
        PRINT*,"---------------WARNING----------------"
        PRINT*,"Nz changed due to Streatching Factor"
        PRINT*,"New Nz value: ",Nz
        PRINT*,"--------------------------------------"
        PRINT*,""
        ALLOCATE(zz(1,Nz))
        zz(1,1)=zz0+zz1                                   
        zz(1,2)=zz(1,1)+dzsf
            
        do i=3,Nz
            zz(1,i)=zz(1,2)+(zz(1,i-1)-zz(1,1))*ZSF
        end do
   end if 
   
   if((ZSF.ne.1.0d0).or.(YSF.ne.1.0d0).or.(XSF.ne.1.0d0)) then
    Ntot=NX*NY*NZ
    PRINT*,"New total number of nodes= ",Ntot
    WRITE(*,fmt="(a,i5,a,i5,a,i5,a)")" New Mesh resolution [ ",Nx," ; ",Ny," ; ",Nz," ]"
    print*,""
   end if
    
    PRINT*,"Mesh ready!"
    PRINT*,""
    PRINT*,"1.3 Writing the TecPlot output file..."
    PRINT*,""    
        
!***************************** Writing the solution file to TecPlot *************************************
    OPEN(unit=11,file='TecPlot_out_Interp.dat',IOSTAT=IOS)
      IF(IOS/=0) THEN
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      PRINT*,"GridMaker - Error #02: Impossible to create output file. Please verify permissions."
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      STOP
    END IF
    REWIND(11)    
    WRITE(11,*)'Variables = "X", "Y", "Z"'
    WRITE(11,*)'ZONE I=',Nx,'J=',Ny,'K=',Nz,'ZONETYPE=Ordered, DATAPACKING=POINT'
    
    
    do i=1,Nz
        do j=1,Ny
            do k=1,Nx
                WRITE(11,fmt="(6(f13.5,1X))",ADVANCE='NO')xx(1,k),yy(1,j),zz(1,i)
                WRITE(11,*)" "
            end do
       end do
       call progress(i)             !********************** Modificação PEDRO RICARDO - Loading bar                                
    end do
    
    CLOSE(11)
    
    !*************** 2 Dimensional grid ********************!
    Select case (Tec2D)
        Case (1)
            Ntemp=Ny/2.0
            y2D=NINT(Ntemp)          
            Ntemp=Nz/2.0
            z2D=NINT(Ntemp)          
            
            OPEN(unit=30,file='Flight_out.dat',IOSTAT=IOS)
            IF(IOS/=0) THEN
              PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              PRINT*,"GridMaker - Error #02: Impossible to create output file. Please verify permissions."
              PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               STOP
            END IF
            REWIND(30)
            WRITE(30,*)'Variables = "X", "Y"'
            WRITE(30,*)'ZONE I=',Nx,'J=',Ny-y2D+1,'ZONETYPE=Ordered, DATAPACKING=POINT'
            
            do i=1,Ny !y2D
                do k=1,Nx
                    WRITE(30,fmt="(2(f16.8,1X))",ADVANCE='NO')xx(1,k),yy(1,i)
                    Write(30,*)" "
                end do
            end do
            
            CLOSE(30)
           
    end select 
    
    CALL CPU_TIME(T2)
         
    Ttot=NINT(T2-T1)                             !Calculates the time of the numerical simulation
    PRINT*,""
    PRINT*,""
    PRINT*,"----------------------------------------"
    PRINT*,"Files ready to interpolate the CFD data!"
    PRINT*,"----------------------------------------"
    PRINT*,""
    PRINT*,"Time of the execution:",Ttot,"seconds"
    PRINT*,""    
   
   end program GridMaker
   
 
!!********************************** Modificação PEDRO RICARDO - Loading Bar 
subroutine progress(j)
  Use globais
  implicit none
  integer(kind=4)::j
  character(len=4)::bar="???%"
  write(unit=bar(1:3),fmt="(i3)") 100*j/NZ
  write(unit=6,fmt="(a1,a4,$)") ,char(13), bar  ! print the progress bar.
  return
end subroutine progress
!***********************************  
