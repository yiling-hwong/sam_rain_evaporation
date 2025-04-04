  
subroutine precip_proc

use vars
use microphysics
use micro_params
use params

implicit none

integer i,j,k
real autor, autos, accrr, accris, accrcs, accrig, accrcg
real dq, omn, omp, omg, qsatt
real pows1, pows2, powg1, powg2, powr1, powr2, tmp
real qii, qcc, qrr, qss, qgg
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)

powr1 = (3 + b_rain) / 4.
powr2 = (5 + b_rain) / 8.
pows1 = (3 + b_snow) / 4.
pows2 = (5 + b_snow) / 8.
powg1 = (3 + b_grau) / 4.
powg2 = (5 + b_grau) / 8.
      
call t_startf ('precip_proc')
     
if(dostatis) then
        
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = q(i,j,k)
     end do
    end do
  end do
         
endif


 
do k=1,nzm
 qpsrc(k)=0.
 qpevp(k)=0.
 do j=1,ny
  do i=1,nx	  
	  
!-------     Autoconversion/accretion 

   if(qn(i,j,k)+qp(i,j,k).gt.0.) then


         omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
         omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))

	 if(qn(i,j,k).gt.0.) then
     
           qcc = qn(i,j,k) * omn
           qii = qn(i,j,k) * (1.-omn)

           if(qcc .gt. qcw0) then
            autor = alphaelq
           else
            autor = 0.
           endif 

           if(qii .gt. qci0) then
            autos = betaelq*coefice(k)
           else
            autos = 0.
           endif 

           accrr = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             accrr = accrrc(k) * qrr ** powr1
           end if
           accrcs = 0.
           accris = 0. 
           if(omp.lt.0.999.and.omg.lt.0.999) then
             qss = qp(i,j,k) * (1.-omp)*(1.-omg)
             tmp = qss ** pows1
             accrcs = accrsc(k) * tmp
             accris = accrsi(k) * tmp 
           end if
           accrcg = 0.
           accrig = 0. 
           if(omp.lt.0.999.and.omg.gt.0.001) then
             qgg = qp(i,j,k) * (1.-omp)*omg
             tmp = qgg ** powg1
             accrcg = accrgc(k) * tmp
             accrig = accrgi(k) * tmp 
           endif
           qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
           qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
           dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+ &
             (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))
           dq = min(dq,qn(i,j,k))
           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
           qn(i,j,k) = qn(i,j,k) - dq
	   qpsrc(k) = qpsrc(k) + dq

!-------    CMuller 03/2013: Comment out the evaporation (partial or total) of precipitating condensates
         elseif(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k).eq.0.) then

           qsatt = 0.
           if(omn.gt.0.001) qsatt = qsatt + omn*qsatw(tabs(i,j,k),pres(k))
           if(omn.lt.0.999) qsatt = qsatt + (1.-omn)*qsati(tabs(i,j,k),pres(k))
           dq = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp
             dq = dq + evapr1(k)*sqrt(qrr) + evapr2(k)*qrr**powr2 
           end if
           if(omp.lt.0.999.and.omg.lt.0.999) then
             qss = qp(i,j,k) * (1.-omp)*(1.-omg)
             dq = dq + evaps1(k)*sqrt(qss) + evaps2(k)*qss**pows2 
           end if
           if(omp.lt.0.999.and.omg.gt.0.001) then
             qgg = qp(i,j,k) * (1.-omp)*omg
             dq = dq + evapg1(k)*sqrt(qgg) + evapg2(k)*qgg**powg2
           end if
           dq = dq * dtn * (q(i,j,k) /qsatt-1.) 
           dq = max(-0.5*qp(i,j,k),dq) 

            if(z(k).lt.1000.) then !do not allow rain evaporation below 1km !CMuller
               dq=0.*dq !YLH
            endif !CMuller 

           qp(i,j,k) = qp(i,j,k) + dq
           q(i,j,k) = q(i,j,k) - dq
           qpevp(k) = qpevp(k) + dq

         else
            dq=qp(i,j,k) 

            if(z(k).lt.1000.) then !do not allow rain evaporation below 1km !CMuller
               dq=0.*dq !YLH
            endif !CMuller

            q(i,j,k) = q(i,j,k) + dq
            qpevp(k) = qpevp(k) - dq
            qp(i,j,k) = 0.
 
         endif
 
    endif

    dq = qp(i,j,k)
    qp(i,j,k)=max(0.,qp(i,j,k))
    q(i,j,k) = q(i,j,k) + (dq-qp(i,j,k))

  end do
 enddo
enddo
    


if(dostatis) then
                  
  call stat_varscalar(q,df,f0,df0,q2leprec)
  call setvalue(qwleprec,nzm,0.)
  call stat_sw2(q,df,qwleprec)

endif

call t_stopf ('precip_proc')

end subroutine precip_proc

