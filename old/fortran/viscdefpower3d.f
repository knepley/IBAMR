cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        Filename : viscdefpower3d.f
c        written by amneet bhalla on Taylor@mech.northwestern.edu
c        created on 07/14/2011
c        It calculates pointwise viscous power dissipated on a 3d patch due to deformation kinematics.

c        power_def = 2*mu*D(u_def):D(u_def), where D(x) = 0.5( \nabla x + (\nabla x)' ) and ":"
c        represents the inner tensor product.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




       subroutine calculateviscdefpower3d(
     &           TAG,
     &           VEL,
     &           VISC_DEF_POWER,
     &           mu,
     &           gcw,
     &           ilower0, iupper0,
     &           ilower1, iupper1,
     &           ilower2, iupper2,
     &           DX )

       implicit none

c      INPUT
       
       integer gcw,ilower0,iupper0, ilower1, iupper1,ilower2,
     &         iupper2

       integer TAG( ilower0-gcw:iupper0+gcw, ilower1-gcw:
     &              iupper1+gcw, ilower2-gcw:iupper2+gcw )

       double precision VEL( ilower0-gcw:iupper0+gcw, ilower1-gcw:
     &                       iupper1+gcw, ilower2-gcw:iupper2+gcw, 0:2)

       double precision   DX(0:2)
       double precision   mu 

c      OUTPUT
       double precision VISC_DEF_POWER(ilower0:iupper0,ilower1:iupper1, 
     &                                 ilower2:iupper2, 0:2)

c      LOCAL VARIABLES

       integer i,j,k, East,West,North,South,Top,Bottom, Center
       double precision dudx,dudy,dudz, dvdx, dvdy, dvdz, dwdx, dwdy,
     &  dwdz
  

c      COMPUTE DERIVATIVES

       do k = ilower2, iupper2,1
        do j = ilower1, iupper1,1
         do i = ilower0, iupper0, 1
   
          Center = TAG(i,j,k)
          if(Center .NE. 0 ) then
            
            East   = TAG(i+1,j,k)
            West   = TAG(i-1,j,k)  
            North  = TAG(i,j+1,k)
            South  = TAG(i,j-1,k)
            Top    = TAG(i,j,k+1)
            Bottom = TAG(i,j,k-1)

c	    COMPUTE 'X' DERIVATIVES 

c           central difference in x 
            if( (East .NE. 0) .AND. (West .NE. 0) ) then     
              dudx = (VEL(i+1,j,k,0) - VEL(i-1,j,k,0))/(2d0*DX(0)) 
              dvdx = (VEL(i+1,j,k,1) - VEL(i-1,j,k,1))/(2d0*DX(0)) 
              dwdx = (VEL(i+1,j,k,2) - VEL(i-1,j,k,2))/(2d0*DX(0))  
            endif

c           forward difference in x   
            if( (East .NE. 0) .AND. (West .EQ. 0) ) then         
               dudx = ( VEL(i+1,j,k,0) - VEL(i,j,k,0) ) / DX(0) 
               dvdx = ( VEL(i+1,j,k,1) - VEL(i,j,k,1) ) / DX(0) 
               dwdx = ( VEL(i+1,j,k,2) - VEL(i,j,k,2) ) / DX(0) 
             endif

c            backward difference in x 
             if( (East .EQ. 0) .AND. (West .NE. 0) ) then         
               dudx = ( VEL(i,j,k,0) - VEL(i-1,j,k,0) ) / DX(0)  
               dvdx = ( VEL(i,j,k,1) - VEL(i-1,j,k,1) ) / DX(0)
               dwdx = ( VEL(i,j,k,2) - VEL(i-1,j,k,2) ) / DX(0)
             endif
         
c            nothing right && left is present in x direction
             if( (East .EQ. 0) .AND. (West .EQ. 0) ) then        
               dudx = 0.0d0                 
               dvdx = 0.0d0 
               dwdx = 0.0d0
             endif


c            COMPUTE 'Y' DERIVATIVES 

c           central difference in y 
            if( (North .NE. 0) .AND. (South .NE. 0) ) then  
               dudy = (VEL(i,j+1,k,0) - VEL(i,j-1,k,0))/(2d0*DX(1)) 
               dvdy = (VEL(i,j+1,k,1) - VEL(i,j-1,k,1))/(2d0*DX(1))  
               dwdy = (VEL(i,j+1,k,2) - VEL(i,j-1,k,2))/(2d0*DX(1))
            endif

c           forward difference in y 
            if( (North .NE. 0) .AND. (South .EQ. 0) ) then       
               dudy = ( VEL(i,j+1,k,0) - VEL(i,j,k,0) ) / DX(1) 
               dvdy = ( VEL(i,j+1,k,1) - VEL(i,j,k,1) ) / DX(1) 
               dwdy = ( VEL(i,j+1,k,2) - VEL(i,j,k,2) ) / DX(1) 
            endif

c           backward difference in y 
            if( (North .EQ. 0) .AND. (South .NE. 0) ) then       
               dudy = ( VEL(i,j,k,0) - VEL(i,j-1,k,0) ) / DX(1) 
               dvdy = ( VEL(i,j,k,1) - VEL(i,j-1,k,1) ) / DX(1) 
               dwdy = ( VEL(i,j,k,2) - VEL(i,j-1,k,2) ) / DX(1) 
            endif
         
c           nothing above && below is present in y direction
            if( (North .EQ.0) .AND. (South .EQ.0) ) then       
               dudy = 0.0d0                 
               dvdy = 0.0d0 
               dwdy = 0.0d0
            endif


c           COMPUTE 'Z' DERIVATIVES 

c           central difference in z 
            if( (Top .NE. 0) .AND. (Bottom .NE.0) ) then       
               dudz = (VEL(i,j,k+1,0) - VEL(i,j,k-1,0))/(2d0*DX(2))  
               dvdz = (VEL(i,j,k+1,1) - VEL(i,j,k-1,1))/(2d0*DX(2))  
               dwdz = (VEL(i,j,k+1,2) - VEL(i,j,k-1,2))/(2d0*DX(2))
            endif

c           forward difference in z 
            if( (Top .NE. 0) .AND. (Bottom .EQ. 0) ) then       
               dudz = ( VEL(i,j,k+1,0) - VEL(i,j,k,0) ) / DX(2)
               dvdz = ( VEL(i,j,k+1,1) - VEL(i,j,k,1) ) / DX(2) 
               dwdz = ( VEL(i,j,k+1,2) - VEL(i,j,k,2) ) / DX(2) 
            endif

c           backward difference in z 
            if( (Top .EQ.0) .AND. (Bottom .NE.0) ) then         
               dudz = ( VEL(i,j,k,0) - VEL(i,j,k-1,0) ) / DX(2)  
               dvdz = ( VEL(i,j,k,1) - VEL(i,j,k-1,1) ) / DX(2) 
               dwdz = ( VEL(i,j,k,2) - VEL(i,j,k-1,2) ) / DX(2) 
            endif

c           nothing above && below is present in z direction
            if( (Top .EQ.0) .AND. (Bottom .EQ. 0) ) then     
               dudz = 0.0d0                 
               dvdz = 0.0d0 
               dwdz = 0.0d0
            endif


c           COMPUTE VISCOUS POWER

            VISC_DEF_POWER(i,j,k,0) = 2d0*mu*( dudx**2 + (0.5d0 * 
     &      (dudy + dvdx) )**2 + (0.5d0 * (dudz + dwdx) )**2  )

            VISC_DEF_POWER(i,j,k,1) = 2d0*mu*( (0.5d0 * (dudy + 
     &      dvdx) )**2 + dvdy**2 + (0.5d0 * (dvdz + dwdy) )**2  )

            VISC_DEF_POWER(i,j,k,2) = 2d0*mu*( (0.5d0 * (dudz + 
     &      dwdx) )**2 + (0.5d0 * (dvdz + dwdy) )**2 + dwdz**2  )

             
          endif ! Center != 0
              
         enddo
        enddo
       enddo

       return
       end



