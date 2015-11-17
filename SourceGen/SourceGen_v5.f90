!****************************************************************************
!
!  PROGRAM: SourceGen Version 5
!
!  PURPOSE: Source generation for the 3D RayTracing methodology.
!
!  Major modification from the last version: 1)Input file changed.
!                                            2)Code restructured.
!                                            3)Sources location generalized for non symmetric jets (TKE based).
!                                            4)Mesh based source spacing.
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
    DOUBLE PRECISION, ALLOCATABLE :: X(:),Y(:),Z(:),U(:,:,:),C(:,:,:),V(:,:,:),W(:,:,:),TKE(:,:,:),EPS(:,:,:),RHO(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Xpl(:),Ypl(:),Zpl(:),ipl(:),jpl(:),kpl(:),Vpl(:)
    DOUBLE PRECISION :: Umax, Tkemax
    REAL :: Vflight
    INTEGER :: Nx,Ny,Nsx,imax,jmax,kmax,Nz,Ny05,Nz05
end module globais

program coaxial
    Use globais
    implicit none
    
    DOUBLE PRECISION, ALLOCATABLE :: Sx(:,:),Vel(:,:),xs(:,:),ys(:,:),zs(:,:)
    DOUBLE PRECISION :: Vx,xmax,ymax,zmax,Mach,Dsx,TKEtemp
    DOUBLE PRECISION :: diaP,diaS,xx0,yy0,zz0,Ldjet,Xend
    REAL :: T1,T2,Ttot
    INTEGER :: IOS
    INTEGER :: i,j,k,posx,NtotSources,test,LdJetInt,ndjet,Nslmax,l,kl,jl,kk,dSz,dSy
    INTEGER, ALLOCATABLE :: cont(:)
    
    PRINT*,"*******************************"
    PRINT*,"         SourceGen V5          "
    PRINT*,"*******************************"
    PRINT*," Version: Any jet - TKE 10%"
    
    CALL CPU_TIME(T1)
    PRINT*,""
    PRINT*,"2. Reading the CFD grid information..."
    PRINT*,""
    
    !Input file for Source Generation code
    OPEN(unit=2,file='SourceGen_input.dat',status='old',IOSTAT=IOS)
      IF(IOS/=0) THEN
       PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       PRINT*,"SourceGen - Error #01: Input file: Open Error. Check the input file and correct!!"
       PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       STOP
      END IF
    READ(2,*)diaP, diaS       !Jet diameter
    READ(2,*)xx0,yy0,zz0      !Initial jet position -> need to be [0,0,0]
    READ(2,*)Ldjet,ndjet      !Size downstream the jet in function of jet's diameter that is required to find source positions
    READ(2,*)dSy,dSz          !Source Spacing Y,Z
    READ(2,*)Nx,Ny,Nz         !Grid
    READ(2,*)Vflight          !Velocity flight - included on 15th June 2010.  
    
    CLOSE(2)
    
   
   !!************Modificação para bocal simples - SMC simulations Abril 2012
    If (diaP==diaS) then
        diaP=0.5*diaS
    end if
   !*********************************
    
    Nz05=NINT(Nz/2.0)   !Calculates the next integer of half domain for the Z and Y direction
    Ny05=NINT(Ny/2.0)
    
    ALLOCATE(X(Nx),Y(Ny),Z(Nz),U(Nx,Ny,Nz),Vel(Nx,Ny),C(Nx,Ny,Nz),RHO(Nx,Ny,Nz),TKE(Nx,Ny,Nz),&
    &EPS(Nx,Ny,Nz),V(Nx,Ny,Nz),W(Nx,Ny,Nz),Xpl(Nx),Ypl(Nx),Zpl(Nx),ipl(Nx),jpl(Nx),kpl(Nx),Vpl(Nx))
    X=0; Y=0; Z=0; U=0; Vel=0; C=0; RHO=0; TKE=0; EPS=0; Xpl=0; Ypl=0; Zpl=0; ipl=0; jpl=0; kpl=0; Vpl=0;
    
    OPEN(unit=1,file='CFD_input.dat',status='old',IOSTAT=IOS)
    IF(IOS/=0) THEN
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      PRINT*,"SourceGen - Error #01: Input file: Open Error. Check the input file and correct!!"
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      STOP
    END IF
   
   Umax=0.0d0
   Tkemax=0.0d0

   READ(1,*)!1
   READ(1,*)!2
   READ(1,*)!3
   READ(1,*)!4
   READ(1,*)!5
   READ(1,*)!6
   READ(1,*)!7
   READ(1,*)!8
   READ(1,*)!9
   READ(1,*)!10
   READ(1,*)!11
   READ(1,*)!12
   READ(1,*)!13
   READ(1,*)!14
   READ(1,*)!15
   READ(1,*)!16
     
   do k=1,Nz
        do j=1,Ny
           do i=1,Nx

           !READ(1,*)X(i),Y(j),Z(k),U(i,j,k),RHO(i,j,k),TKE(i,j,k),EPS(i,j,k),C(i,j,k),DU(i,j,k)
            READ(1,*)X(i),Y(j),Z(k),U(i,j,k),V(i,j,k),W(i,j,k),RHO(i,j,k),TKE(i,j,k),EPS(i,j,k),C(i,j,k)

            TKEtemp=TKE(i,j,k)
            Vx=U(i,j,k)
                If (Vx>Umax) then
                   !Maximum velocity on the domain
                    Umax=Vx         
                   !Position of the Umax
                    xmax=X(i)       
                    ymax=Y(j) 
                    zmax=Z(k)          
                   !Node localization of the Umax
                    imax=i          
                    jmax=j        
                    kmax=k  
                end if

                !********* TKE máxima
                If (TKEtemp>Tkemax) then
                    Tkemax=TKEtemp     
                end if
        end do
    end do 
  end do
  
 CLOSE(1)
  !********* Bloco adicionado PEDRO RICARDO 2015 - Cálculo da posição de velocidade máxima em cada plano
  k=Nz05            !Plano de simetria
  do i=1,Nx
  Vpl(i)=0.0
    do j=1,Ny
        Vx=U(i,j,k)
        If (Vx>Vpl(i)) then
           !Velocidade máxima no plano
            Vpl(i)=Vx         
           !Posição da Velocidade maxima no plano
            Xpl(i)=X(i)       
            Ypl(i)=Y(j) 
            Zpl(i)=Z(k)          
           !Nó da velocidade máxima no plano
            ipl(i)=i          
            jpl(i)=j        
            kpl(i)=k  
         end if
    end do
end do      

    
    PRINT*,"-----------"  
    PRINT*,"Data ready!"
    PRINT*,"-----------"
    PRINT*,""
    
    Write(*,fmt="(a,f8.2)")"The maximum velocity found on the flow is ",Umax
    Write(*,fmt="(a,E8.2,a,E8.2,a,E8.2,a)")"located at P=[ ",xmax," ; ",ymax," ; ",zmax," ]"
    PRINT*,""    
    PRINT*,"2.1 Generating the sources..."
    
    LdjetInt=10*NINT(Ldjet)
    
    Xend=Ldjet*diaS              !The end point position dowstream the jet to insert sources
    Dsx=Xend/ndjet               !Espaçamento das fontes 
    
    If (Xend>X(Nx)) then
        PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        PRINT*,"SourceGen - Error #03: The endpoint for source location is outside"
        PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        STOP
    end if
    
    
    ALLOCATE(Sx(LdjetInt+1,1))
    Sx=0 !Inicializa todo o Array como ZERO
    
    Sx(1,1)=0.01*diaS+xx0      !Position of the first plane defined to be 1% of the jet diameter in the downstream position
    Mach=U(imax,jmax,kmax)/C(imax,jmax,kmax) !Determines the local Mach of the first plane position. 
                                             !Considered to be the first node near the jet center position (xx0,yy0,zz0)
    
    do i=2,LdjetInt
        if (Mach<0.1) then
            Sx(i,1)=Sx(i-1,1)+2*diaS
        else
           !Sx(i,1)=Sx(i-1,1)+diaS*(-1.6667*Mach+2.5) ! Jato ISVR
           !Sx(i,1)=Sx(i-1,1)+diaS/2.0d0 ! Crossflow: espaçamento fixo
            Sx(i,1)=Sx(i-1,1)+Dsx ! Crossflow: espaçamento dinâmico de acordo com Ldjet do input
        end if
        if (Sx(i,1)>Xend) then
            Nsx=i-1
            exit
        end if
        test=1
        CALL Near(X(1),Sx(i,1),posx,test,Mach)
        Nsx=i
   end do
      
    Nslmax=Nz*Ny
    
    ALLOCATE(xs(Nslmax,Nsx),ys(Nslmax,Nsx),zs(Nslmax,Nsx),cont(Nsx+1))
    xs=0; ys=0; zs=0;
    
    NtotSources=0
    cont=0
    kl=0
    jl=0
    k=1
    ! Principal LOOP
    do l=1,NSx
        test=0
        CALL Near(X(1),Sx(l,1),posx,test,Mach)
        
        if (l==1) then
        do k=1,Nz
            
            do j=1,Ny
                TKEtemp=TKE(posx,j,k)
                
                if (TKEtemp>(Tkemax*0.1d0)) then
                    
                    cont(l)=cont(l)+1
                    ys(cont(l),l)=Y(j)
                    zs(cont(l),l)=Z(k)
                    xs(cont(l),l)=Xpl(posx)
                    NtotSources=NtotSources+1
                    
                end if   
            end do
        end do
        k=1
        else
        do kk=1,int(Nz/dSz)
            
            do j=1,Ny
                TKEtemp=TKE(posx,j,k)
                
                if (TKEtemp>Tkemax*0.1d0) then
                    
                    if (abs(j-jl)>dSy .or. cont(l)==0) then
                    cont(l)=cont(l)+1
                    ys(cont(l),l)=Y(j)
                    zs(cont(l),l)=Z(k)
                    xs(cont(l),l)=Xpl(posx)
                    NtotSources=NtotSources+1
                    kl=k; jl=j;
                    end if
                    
                end if
            k=kk*dSz    
            end do
        end do
        end if        
    end do      ! Fim do loop principal
    !************
    
    PRINT*,""
    Write(*,fmt="(a,i5)")"Number of sources that will be generated: ", NtotSources
    Write(*,fmt="(a,i5)")"Number of source planes: ", Nsx
    PRINT*,""
    
    PRINT*,"----------------------------------------------"
    PRINT*,"Simulation finished! Writing the output files."
    PRINT*,"----------------------------------------------"
    
    OPEN(unit=3,file="Sources_New.dat",IOSTAT=IOS)
         IF(IOS/=0) THEN
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      PRINT*,"SourceGen - Error #02: Impossible to create output file. Please verify permissions."
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      STOP
     END IF

     OPEN(unit=5,file="Sources_tecplot.dat",IOSTAT=IOS) ! Arquivo p/ visualização das Sources
     IF(IOS/=0) THEN
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      PRINT*,"SourceGen - Error #02: Impossible to create output file. Please verify permissions."
      PRINT*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      STOP
     END IF
    WRITE(5,*)"VARIABLES = X Y Z"    
    WRITE(3,*)NtotSources
    
    do i=1,Nsx
        do j=1,cont(i)
            WRITE(3,fmt="(3(f11.6,1X))",ADVANCE='NO')xs(j,i),ys(j,i),zs(j,i)
            WRITE(5,fmt="(3(f11.6,1X))",ADVANCE='NO')xs(j,i),ys(j,i),zs(j,i)
            Write(3,*)" "
            Write(5,*)" "
        end do

    end do
    
    CLOSE(3)
    CLOSE(5)
    
    CALL CPU_TIME(T2)
    Ttot=T2-T1
   
    PRINT*,""
    PRINT*,"Sources's coordinates are ready! Program finished."
    WRITE(*,fmt="(a,f8.2)")"Execution time in second was ",Ttot
    PRINT*,"..."
    
    
    DEALLOCATE (xs,ys,zs,Sx,X,Y,Z,U,Vel,C,RHO,TKE,EPS,Xpl,Ypl,Zpl,ipl,jpl,kpl,Vpl)
    end program Coaxial
    
!****************************************************************************************
!****************************************************************************************
    
SUBROUTINE Near(x1,Sx,posx,test,Mach)
! This subroutine finds the nearest grid point in the x direction of the source position. This
! information will be used to calculates the approximate length of the jet plume in the Y direction.
!
  Use globais
  
  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: Sx,x1
  DOUBLE PRECISION, INTENT(INOUT) :: Mach
  INTEGER, INTENT(OUT) :: posx,test
  DOUBLE PRECISION :: dif,dift
  INTEGER :: i
  
  dif=abs(Sx-x1)
  posx=1
  do i=1,Nx
        dift=abs(Sx-X(i))
        if (dift<dif) then
            dif=dift
            posx=i
        end if
  end do
  
  if (test==1)then
    CALL LinearInterp(Sx,posx,Mach)
  else
    RETURN
  end if
  test=0

END SUBROUTINE Near

SUBROUTINE LinearInterp(Sx,posx,Mach)
     
      USE globais
     
     IMPLICIT NONE
     
     DOUBLE PRECISION, INTENT(IN) :: Sx
     DOUBLE PRECISION, INTENT(INOUT) :: Mach
     INTEGER, INTENT(IN) :: posx
     DOUBLE PRECISION :: Utemp,Ctemp,Lx,dif,U1,U2,C1,C2
     INTEGER :: nod2,nod1
     
     dif=Sx-X(posx)
     If (dif>0.0) then
        nod1=posx
        nod2=nod1+1
     else
        nod1=posx
        nod2=nod1-1
     end if
     
     if (nod1==1) then
        nod2=nod1
        Lx=1
     elseif (nod1>=Nx)then
        nod2=nod1
        Lx=1
     else
        Lx=abs(X(nod1)-X(nod2))
     end if
     
     If (Vpl(nod1)<0) then
        U1=Umax
        C1=C(imax,jmax,kmax)
     else
        U1=Vpl(nod1)
        C1=C(nod1,jmax,kmax)
     end if
     If (Vpl(nod2)<0) then
        U2=Umax
        C2=C(imax,jmax,kmax)
     else
        U2=Vpl(nod2)
        C2=C(nod2,jmax,kmax)
     end if
         
     Utemp=(U1*abs(Sx-X(nod1))/Lx)+(U2*abs(Sx-X(nod2))/Lx)
     
     Ctemp=(C1*abs(Sx-X(nod1))/Lx)+(C2*abs(Sx-X(nod2))/Lx)
     
     Mach=Utemp/Ctemp
     
END SUBROUTINE LinearInterp
