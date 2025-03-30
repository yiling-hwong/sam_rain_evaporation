module microphysics

! module for original SAM bulk microphysics
! Marat Khairoutdinov, 2006

use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information 
use params, only: doprecip, docloud
use micro_params
implicit none

!----------------------------------------------------------------------
!!! required definitions:

integer, parameter :: nmicro_fields = 3   ! total number of prognostic water vars

!!! microphysics prognostic variables are storred in this array:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

integer, parameter :: flag_wmass(nmicro_fields) = (/1,1,0/)
integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
integer, parameter :: index_cloud_ice = -1   ! index for cloud ice (sedimentation)
integer, parameter :: flag_precip(nmicro_fields) = (/0,1,1/)

! both variables correspond to mass, not number
integer, parameter :: flag_number(nmicro_fields) = (/0,0,1/)

! SAM1MOM 3D microphysical fields are output by default.
integer, parameter :: flag_micro3Dout(nmicro_fields) = (/0,0,1/)

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface flux of tracers
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary flux of tracers

!!! these arrays are needed for output statistics:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mklsadv(nz,1:nmicro_fields)  ! tendency due to large-scale vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion

!======================================================================
! UW ADDITIONS

!bloss: arrays with names/units for microphysical outputs in statistics.
character*3, dimension(nmicro_fields) :: mkname
character*80, dimension(nmicro_fields) :: mklongname
character*10, dimension(nmicro_fields) :: mkunits
real, dimension(nmicro_fields) :: mkoutputscale

! END UW ADDITIONS
!======================================================================

!------------------------------------------------------------------
! Optional (internal) definitions)

! make aliases for prognostic variables:
! note that the aliases should be local to microphysics

real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! total nonprecipitating water
real qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! total precipitating water
real conp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! precipitating water number concentration
equivalence (q(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (qp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))
equivalence (conp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,3))

real qn(nx,ny,nzm)  ! cloud condensate (liquid + ice)

real qpsrc(nz)  ! source of precipitation microphysical processes
real qpfall(nz) ! source of precipitating water due to fall out in a given level
real qpevp(nz)  ! sink of precipitating water due to evaporation

real vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
subroutine micro_setparm()
  use vars
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder
  
  !======================================================================
  ! UW ADDITION
  NAMELIST /MICRO_DRIZZLE/ &
       Nc0 ! (cm-3) Prescribed cloud drop concentration 

  !bloss: Create dummy namelist, so that we can figure out error code
  !       for a mising namelist.  This lets us differentiate between
  !       missing namelists and those with an error within the namelist.
  NAMELIST /BNCUIODSBJCB/ place_holder

  Nc0 = 40. ! default cloud drop number concentration

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 

  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_DRIZZLE namelist
  read (55,MICRO_DRIZZLE,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_DRIZZLE namelist'
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '***** No MICRO_DRIZZLE namelist in prm file ********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

  if(masterproc) then
     write(*,*) 'Cloud droplet number concentration = ', Nc0, '/cm3'
  end if
  ! END UW ADDITION
  !======================================================================

end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:

subroutine micro_init()

  use vars, only: q0
  use grid, only: nrestart
  integer k
  
  if(doprecip) call precip_init() 

  if(nrestart.eq.0) then

     micro_field = 0.
     do k=1,nzm
      q(:,:,k) = q0(k)
     end do
     qn = 0.
     fluxbmk = 0.
     fluxtmk = 0.

     if(docloud) then
       call cloud()
       call micro_diagnose()
     end if

  end if

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.

  qpsrc = 0.
  qpevp = 0.

  mkname(1) = 'QT'
  mklongname(1) = 'TOTAL WATER: VAPOR + CONDENSATE'
  mkunits(1) = 'g/kg'
  mkoutputscale(1) = 1.e3

  mkname(2) = 'QR'
  mklongname(2) = 'DRIZZLE MASS MIXING RATIO'
  mkunits(2) = 'g/kg'
  mkoutputscale(2) = 1.e3

  mkname(3) = 'CONP'
  mklongname(3) = 'DRIZZLE DROP CONCENTRATION'
  mkunits(3) = 'cm-3'
  mkoutputscale(3) = 1.e-6

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
subroutine micro_flux()

  use vars, only: fluxbq, fluxtq

  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (bayond advection and SGS diffusion):
!
subroutine micro_proc()

   use grid, only: icycle
   integer k

   ! Update bulk coefficient
   if(doprecip.and.icycle.eq.1) call precip_init() 

   if(docloud) then
     call cloud()
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine micro_diagnose()
 
   use vars
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = q(i,j,k) - qn(i,j,k)
       qcl(i,j,k) = qn(i,j,k)
       qci(i,j,k) = 0.
       qpl(i,j,k) = qp(i,j,k)
       qpi(i,j,k) = 0.
     end do
    end do
   end do

end subroutine micro_diagnose

!----------------------------------------------------------------------
!!! function to compute terminal velocity for precipitating variables:

real function term_vel_qp(i,j,k,ind)
  
  integer, intent(in) ::  i,j,k ! current indexes
  real rvdr
  integer ind ! placeholder dummy variable

  term_vel_qp = 0.
  if(qp(i,j,k).gt.qp_threshold) then
      conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))
      rvdr = (coefrv*qp(i,j,k)/conp(i,j,k))**0.333
      term_vel_qp= max(0.,1.2e4*rvdr-0.2)
  else
      qp(i,j,k)=0.
      conp(i,j,k) = 0.
  endif

end function term_vel_qp

real function term_vel_conp(i,j,k,ind)
  
  integer, intent(in) ::  i,j,k ! current indexes
  real rvdr
  integer ind ! placeholder dummy variable

  term_vel_conp = 0.
  if(qp(i,j,k).gt.qp_threshold) then
      conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))
      rvdr = (coefrv*qp(i,j,k)/conp(i,j,k))**0.333
      term_vel_conp= max(0.,0.7e4*rvdr-0.1)
  else
      qp(i,j,k)=0.
      conp(i,j,k) = 0.
  endif

end function term_vel_conp

!----------------------------------------------------------------------
!!! compute sedimentation 
!
subroutine micro_precip_fall()
  
  use vars
  use params, only : pi

  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real f0(nzm),df0(nzm)
  real omega(nx,ny,nzm)
  integer ind
  integer i,j,k

! Initialize arrays that accumulate surface precipitation flux

 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if

 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
 end do

 call precip_fall(qp, term_vel_qp, 0, omega, ind)
 call precip_fall(conp, term_vel_conp, 3, omega, ind)

! keep reasonable values:

do k=1,nzm
  do j=1,ny
      do i=1,nx
         if(qp(i,j,k).gt.qp_threshold) then
             conp(i,j,k) = max(qp(i,j,k)*coefconpmin, conp(i,j,k))
         else
             q(i,j,k) = q(i,j,k) +qp(i,j,k)
             qp(i,j,k)=0.
             conp(i,j,k) = 0.
         endif
      end do
  end do
end do

 if(dostatis) then
   do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
          df(i,j,k) = t(i,j,k)
       end do
     end do
   end do
   call stat_varscalar(t,df,f0,df0,t2leprec)
   call setvalue(twleprec,nzm,0.)
   call stat_sw2(t,df,twleprec)
 endif

end subroutine micro_precip_fall

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine micro_statistics()
  
  use vars
  use hbuffer, only: hbuf_avg_put, hbuf_put
  use params, only : lcond

  real tmp(2), factor_xy 
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm), nrz(nzm), omg
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

  call t_startf ('statistics')

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,1) = mkwsb(k,1) * tmp(1) * rhow(k) * lcond
      mkwle(k,1) = mkwle(k,1)*tmp(2)*rhow(k)*lcond + mkwsb(k,1)
      if(docloud.and.doprecip) then
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
  end do

  call hbuf_put('QTFLUX',mkwle(:,1),factor_xy)
  call hbuf_put('QTFLUXS',mkwsb(:,1),factor_xy)
  call hbuf_put('QPFLUX',mkwle(:,2),factor_xy)
  call hbuf_put('QPFLUXS',mkwsb(:,2),factor_xy)

  do k=1,nzm
    qcz(k) = 0.
    qiz(k) = 0.
    qrz(k) = 0.
    qsz(k) = 0.
    qgz(k) = 0.
    do j=1,ny
    do i=1,nx
      qcz(k)=qcz(k)+qcl(i,j,k)
      qiz(k)=0.
      qrz(k)=qrz(k)+qpl(i,j,k)
      qsz(k)=0.
      qgz(k)=0.
    end do
    end do
  end do

  call hbuf_put('QC',qcz,1.e3*factor_xy)
  call hbuf_put('QR',qrz,1.e3*factor_xy)

  call hbuf_avg_put('CONP',conp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.e-6)

  call hbuf_put('QTADV',mkadv(:,1)+qifall,factor_xy*86400000./dtn)
  call hbuf_put('QTDIFF',mkdiff(:,1),factor_xy*86400000./dtn)
  call hbuf_put('QTSINK',qpsrc,-factor_xy*86400000./dtn)
  call hbuf_put('QTSRC',qpevp,-factor_xy*86400000./dtn)
  call hbuf_put('QPADV',mkadv(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPDIFF',mkdiff(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPFALL',qpfall,factor_xy*86400000./dtn)
  call hbuf_put('QPSRC',qpsrc,factor_xy*86400000./dtn)
  call hbuf_put('QPEVP',qpevp,factor_xy*86400000./dtn)

  do n = 1,nmicro_fields
     if(flag_wmass(n).lt.1) then
        ! remove factor of rho from number concentrations
        mklsadv(1:nzm,n) = mklsadv(1:nzm,n)*rho(:)
     end if
     call hbuf_put(trim(mkname(n))//'LSADV', &
          mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
  end do

  do n = 1,ncondavg

     do k=1,nzm
        qcz(k) = 0.
        qrz(k) = 0.
        nrz(k) = 0.
        do j=1,ny
           do i=1,nx
              qcz(k)=qcz(k)+qcl(i,j,k)*condavg_mask(i,j,k,n)
              qrz(k)=qrz(k)+qpl(i,j,k)*condavg_mask(i,j,k,n)
              ! drizzle number concentration
              nrz(k)=nrz(k)+rho(k)*micro_field(i,j,k,3)*condavg_mask(i,j,k,n)
           end do
        end do
     end do

     call hbuf_put('QC' // TRIM(condavgname(n)),qcz,1.e3)
     call hbuf_put('QR' // TRIM(condavgname(n)),qrz,1.e3)
     call hbuf_put('CONP' // TRIM(condavgname(n)),qrz,1.e-6)
  end do

  call t_stopf ('statistics')

end subroutine micro_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
  call fminmax_print('conp:',conp,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics 
!
subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

  use vars


   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr, n


   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUX'
   deflist(count) = 'Nonprecipitating water flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUXS'
   deflist(count) = 'Nonprecipitating-water flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUX'
   deflist(count) = 'Precipitating-water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUXS'
   deflist(count) = 'Precipitating-water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QC'
   deflist(count) = 'Liquid cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'CONP'
   deflist(count) = 'Drizzle drop concentration'
   unitlist(count) = '1/cm3'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QR'
   deflist(count) = 'Rain water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   do n = 1,nmicro_fields
      count = count + 1
      trcount = trcount + 1
      namelist(count) = TRIM(mkname(n))//'LSADV'
      deflist(count) = 'Source of '//TRIM(mklongname(n))//' due to large-scale vertical advection'
      unitlist(count) = TRIM(mkunits(n))//'day'
      status(count) = 1    
      average_type(count) = 0
   end do

  !bloss: setup to add an arbitrary number of conditional statistics
   do n = 1,ncondavg

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QC' // TRIM(condavgname(n))
      deflist(count) = 'Mean Liquid cloud water in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QR' // TRIM(condavgname(n))
      deflist(count) = 'Mean Drizzle Mixing ratio in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'CONP' // TRIM(condavgname(n))
      deflist(count) = 'Mean drizzle drop concentration in ' // TRIM(condavglongname(n))
      unitlist(count) = '/cm3'
      status(count) = 1    
      average_type(count) = n

   end do

end subroutine micro_hbuf_init


!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

! -------------------------------------------------------------------------------
! Dummy effective radius functions:

function Get_reffc() ! liquid water
  real, pointer, dimension(:,:,:) :: Get_reffc
end function Get_reffc

function Get_reffi() ! ice
  real, pointer, dimension(:,:,:) :: Get_reffi
end function Get_reffi


end module microphysics



