module simple_ocean

!------------------------------------------------------------
! Purpose:
!
! A collection of routines used to specify fixed 
! or compute interactive SSTs, like slab-ocean model, etc.
!
! Author: Marat Khairoutdinov
! Based on dynamic ocean impelemntation from the UW version of SAM.
!------------------------------------------------------------

use grid
implicit none

public set_sst     ! set SST 
public sst_evolve ! evolve SST according to a model set by the ocean_type

CONTAINS


SUBROUTINE set_sst()

 use vars, only: sstxy,t00
 use params, only: tabs_s, delta_sst, &
      ocean_type,bubble_y0,bubble_x0,bubble_radius_hor,bubble_dtemp

! parameters of the sinusoidal SST destribution 
! along the X for Walker-type simulatons( ocean-type = 1):

 real(8) tmpx(nx), pii, lx,ly, yy
 integer i,j, it,jt, ihot, jhot, nhot
real xxx,yyy, ttt

 select case (ocean_type)

   case(0) ! fixed constant SST

      sstxy = tabs_s - t00

   case(1) ! Sinusoidal distribution along the x-direction:

     lx = float(nx_gl)*dx
     do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
     end do
     pii = atan2(0.d0,-1.d0)
     do j=1,ny
       do i=1,nx
         sstxy(i,j) = tabs_s-delta_sst*cos(2.*pii*tmpx(i)/lx) - t00
       end do
     end do
   
   case(2) ! Sinusoidal distribution along the y-direction:
     
     call task_rank_to_index(rank,it,jt)
     
     pii = atan2(0.d0,-1.d0)
     lx = float(ny_gl)*dy
     do j=1,ny
        yy = dy*(j+jt-(ny_gl+YES3D-1)/2-1)
       do i=1,nx
         sstxy(i,j) = tabs_s+delta_sst*(2.*cos(pii*yy/lx)-1.) - t00
       end do
     end do 

   case(3) ! Hot spot

     call task_rank_to_index(rank,it,jt)
     ly = float(ny_gl)*dy
     lx = float(nx_gl)*dx
     pii = atan2(0.d0,-1.d0)
     yy=pii*bubble_radius_hor**2/(lx*ly) 
     do j=1,ny
      yyy = dy*(j+jt)
       do i=1,nx
          xxx = dx*(i+it)
          if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2) then
            !rrr =
            !atan((bubble_radius_hor-sqrt((xxx-bubble_x0)**2+(yyy-bubble_y0)**2))/(2.*dx))/atan(bubble_radius_hor/(2.*dx))
            sstxy(i,j) =tabs_s + bubble_dtemp*(1-yy) - t00  !bubble_dtemp*rrr*sin(2.*acos(-1.)*time/86400.)
          else
            sstxy(i,j)=tabs_s - bubble_dtemp*yy -t00
          end if
       end do
      end do


   case(5) !CMuller: Hot spots

      call task_rank_to_index(rank,it,jt)
      ly = float(ny_gl)*dy
      lx = float(nx_gl)*dx
      pii = atan2(0.d0,-1.d0)
      nhot = 8 
      yy=(float(nhot))**2*pii*bubble_radius_hor**2/(lx*ly)

     do j=1,ny
      yyy = dy*(j+jt)
      do i=1,nx
       xxx = dx*(i+it)
       sstxy(i,j)=tabs_s - bubble_dtemp*yy -t00
       do ihot=1,nhot
         bubble_x0 = 0.5*lx/float(nhot) + float(ihot-1)*lx/float(nhot)
         do jhot=1,nhot
          bubble_y0 = 0.5*ly/float(nhot) + float(jhot-1)*ly/float(nhot)
          if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2) then
            sstxy(i,j) = tabs_s + bubble_dtemp*(1-yy) -t00  
          end if
         end do
       end do
      end do
     end do
   !CMuller 


   case default

     if(masterproc) then
         print*, 'unknown ocean type in set_sst. Exitting...'
         call task_abort
     end if

 end select

end subroutine set_sst



SUBROUTINE sst_evolve
 use vars, only: sstxy, t00, fluxbt, fluxbq, rhow,qocean_xy
 use params, only: cp, lcond, tabs_s, ocean_type, dossthomo, &
                   depth_slab_ocean, Szero, deltaS, timesimpleocean
 use rad, only: swnsxy, lwnsxy

 real, parameter :: rhor = 1000. ! density of water (kg/m3)
 real, parameter :: cw = 4187.   ! Liquid Water heat capacity = 4187 J/kg/K
 real factor_cp, factor_lc, qoceanxy
 real tmpx(nx), lx
 real(8) sss,ssss
 integer i,j

      if(time.lt.timesimpleocean) return

      lx = float(nx_gl)*dx
      do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
      end do

      ! Define weight factors for the mixed layer heating due to
      ! the model's sensible and latent heat flux.
      factor_cp = rhow(1)*cp
      factor_lc = rhow(1)*lcond

      ! Use forward Euler to integrate the differential equation
      ! for the ocean mixed layer temperature: dT/dt = S - E.
      ! The source: CPT?GCSS WG4 idealized Walker-circulation 
      ! RCE Intercomparison proposed by C. Bretherton.
      do j=1,ny
         do i=1,nx
           qoceanxy = Szero + deltaS*abs(2.*tmpx(i)/lx - 1)
	   qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy

            sstxy(i,j) = sstxy(i,j) &
                 + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
                 - lwnsxy(i,j)               & ! LW Radiative Heating
                 - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                 - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                 + qoceanxy)            & ! Ocean Heating
                 /(rhor*cw*depth_slab_ocean)        ! Convert W/m^2 Heating to K/s
         end do
      end do

     if(dossthomo) then
        sss = 0.
        do j=1,ny
         do i=1,nx
           sss = sss + sstxy(i,j)
         end do
        end do
        sss = sss / dble(nx*ny)
        if(dompi) then
            call task_sum_real8(sss,ssss,1)
            sss = ssss /float(nsubdomains)
        end if ! dompi
        if(ocean_type.eq.2) then
            tabs_s = sss + t00
            call set_sst()
        else
           sstxy(:,:) = sss
        end if
     end if

end subroutine sst_evolve


end module simple_ocean
