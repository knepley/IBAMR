
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        Filename : viscdefpower2d.f
c        written by amneet bhalla on Taylor@mech.northwestern.edu
c        created on 07/14/2011
c        It calculates pointwise viscous power dissipated on a 2d patch due to deformation kinematics.

c        power_def = 2*mu*D(u_def):D(u_def), where D(x) = 0.5( \nabla x + (\nabla x)' ) and ":"
c        represents the inner tensor product.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       
       subroutine calculateviscdefpower2d(
     &           TAG,
     &           VEL,
     &           VISC_DEF_POWER,
     &           mu,
     &           gcw,
     &           ilower0, iupper0,
     &           ilower1, iupper1,
     &           DX )

       implicit none

c      INPUT

       integer  gcw,ilower0,iupper0, ilower1, iupper1
       integer  TAG( ilower0-gcw:iupper0+gcw, ilower1-gcw:
     &               iupper1+gcw )

       double precision VEL( ilower0-gcw:iupper0+gcw, ilower1-gcw:
     &                       iupper1+gcw,0:1)

       double precision DX(0:1)
       double precision mu 

c      OUTPUT
       double precision   VISC_DEF_POWER(ilower0:iupper0,ilower1:
     &                                   iupper1,0:1)

c      LOCAL VARIABLES

       integer i,j, East,West,North,South,Center
       double precision dudx,dudy,dvdx, dvdy

       do j = ilower1, iupper1,1
        do i = ilower0, iupper0, 1
        
         Center = TAG(i,j) 
         if( Center .NE. 0 ) then   

           East   = TAG(i+1,j)
           West   = TAG(i-1,j)
           North  = TAG(i,j+1)
           South  = TAG(i,j-1)

c	   COMPUTE 'x' DERIVATIVES FIRST

c          central difference in x
           if( (East .NE. 0) .AND. (West .NE. 0) ) then   
              dudx = ( VEL(i+1,j,0) - VEL(i-1,j,0) ) / (2d0 * DX(0)) 
              dvdx = ( VEL(i+1,j,1) - VEL(i-1,j,1) ) / (2d0 * DX(0))  
           endif

c          forward difference in x
           if((East .NE. 0) .AND. (West .EQ. 0) ) then  
              dudx = ( VEL(i+1,j,0) - VEL(i,j,0) ) / DX(0) 
              dvdx = ( VEL(i+1,j,1) - VEL(i,j,1) ) / DX(0) 
           endif

c          backward difference in x
           if((East .EQ. 0) .AND. (West .NE. 0) ) then  
              dudx = ( VEL(i,j,0) - VEL(i-1,j,0) ) / DX(0) 
              dvdx = ( VEL(i,j,1) - VEL(i-1,j,1) ) / DX(0) 
           endif
         
c          nothing right && left is present in x direction
           if((East .EQ. 0) .AND. (West .EQ. 0) ) then
              dudx = 0.0d0                 
              dvdx = 0.0d0 
           endif


c          COMPUTE 'y' DERIVATIVES

c          central difference in y 
           if((North .NE. 0) .AND. (South .NE. 0) ) then 
               dudy = ( VEL(i,j+1,0) - VEL(i,j-1,0) ) / (2d0 * DX(1))
               dvdy = ( VEL(i,j+1,1) - VEL(i,j-1,1) ) / (2d0 * DX(1))  
           endif

c          forward difference in y
           if( (North .NE. 0) .AND. (South .EQ. 0) ) then 
               dudy = ( VEL(i,j+1,0) - VEL(i,j,0) ) / DX(1) 
               dvdy = ( VEL(i,j+1,1) - VEL(i,j,1) ) / DX(1) 
           endif

c           backward difference in y 
            if( (North .EQ. 0) .AND. (South .NE. 0 ) ) then 
               dudy = ( VEL(i,j,0) - VEL(i,j-1,0) ) / DX(1)
               dvdy = ( VEL(i,j,1) - VEL(i,j-1,1) ) / DX(1) 
            endif

c           nothing above && below is present in y direction        
            if((North .EQ. 0) .AND. (South .EQ. 0) ) then 
               dudy = 0.0d0  
               dvdy = 0.0d0 
            endif



c           COMPUTE VISCOUS POWER

            VISC_DEF_POWER(i,j,0) = 2d0*mu*( dudx**2 + (0.5d0 * (dudy 
     &                                                    + dvdx) )**2 )

            VISC_DEF_POWER(i,j,1) = 2d0*mu*( dvdy**2 + (0.5d0 * (dudy 
     &                                                    + dvdx) )**2 )

         endif ! Center /= 0

        enddo
       enddo

       return
       end



