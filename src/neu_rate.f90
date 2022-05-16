module neu_rate

use const_mod
use gauss_integral

implicit none
public:: &
     pair, &
     urca, &
     func1, &
     func2, &
     func3, &
     brem, &
     plasmon
contains



!c..This routine calculates the neutrino emissivity of the 
!c..electron-positron annihilation.
!c..Input:
!c..      xkt temperature in unit of MeV, and
!c..      xetae chemical potential of electron in unit of xkt.
!c..Output:
!c..      ev the total rate in unit of ergs /s/cm^3,
!c..      ev1 the emission rate of electron neutrino,
!c..      ev2 that of mu neutrino
!c..      ev3 that of tau neutrino
!c..

	subroutine pair(xkt,xetae,ev,ev1,ev2,ev3)
	
	double precision xetae,betae,ev,ev1,ev2,ev3,f12,f32,f52,f72, &
      fdeta,fdtheta,cv1e,cv2e,cv1m,cv2m,cv1t,cv2t,QC,xkt, &
      xetapos,u_1,u0,u1,u2,fai_1,fai0,fai1,fai2,cv,ca
	parameter(QC=1.023d23)

	betae=xkt/0.510998910d0
	cv=2.0*0.23d0+0.5d0
	ca=0.5d0
       	cv1e=cv*cv+ca*ca
  	cv2e=cv*cv-ca*ca
	cv=2.0*0.23d0-0.5d0
	ca=-0.5d0
     	cv1m=cv*cv+ca*ca
	cv2m=cv*cv-ca*ca
	cv=2.0*0.23d0-0.5d0
	ca=-0.5d0
     	cv1t=cv*cv+ca*ca
	cv2t=cv*cv-ca*ca
!	print*,'c+=',CV1E+CV1M+CV1T
!	print*,'c-=',CV2E+CV2M+CV2T

	call dfermi(0.5d0,xetae,betae,f12,fdeta,fdtheta)
	call dfermi(1.5d0,xetae,betae,f32,fdeta,fdtheta)
	call dfermi(2.5d0,xetae,betae,f52,fdeta,fdtheta)
	call dfermi(3.5d0,xetae,betae,f72,fdeta,fdtheta)

	u_1=f12
	u0=f12+betae*f32
	u1=f12+2.0*betae*f32+betae*betae*f52
	u2=f12+3.0*betae*f32+3.0*betae*betae*f52+betae**3*f72

	xetapos=-xetae-2.0/betae
	call dfermi(0.5d0,xetapos,betae,f12,fdeta,fdtheta)
	call dfermi(1.5d0,xetapos,betae,f32,fdeta,fdtheta)
	call dfermi(2.5d0,xetapos,betae,f52,fdeta,fdtheta)
	call dfermi(3.5d0,xetapos,betae,f72,fdeta,fdtheta)

	fai_1=f12
	fai0=f12+betae*f32
	fai1=f12+2.0*betae*f32+betae*betae*f52
	fai2=f12+3.0*betae*f32+3.0*betae*betae*f52+betae**3*f72

	ev1=QC/(36.0*PI)*(2.0/PI**4)*betae**3* &
     	  (cv1e*(8.0*(fai1*u2+fai2*u1)-2.0*(fai_1*u2+fai2*u_1)+ &
         7.0*(fai0*u1+fai1*u0)+5.0*(fai0*u_1+fai_1*u0))+ &
         9.0*cv2e*(fai0*(u1+u_1)+(fai_1+fai1)*u0))

	ev2=QC/(36.0*PI)*(2.0/PI**4)*betae**3* &
     	  (cv1m*(8.0*(fai1*u2+fai2*u1)-2.0*(fai_1*u2+fai2*u_1)+ &
         7.0*(fai0*u1+fai1*u0)+5.0*(fai0*u_1+fai_1*u0))+ &
         9.0*cv2m*(fai0*(u1+u_1)+(fai_1+fai1)*u0))

	ev3=QC/(36.0*PI)*(2.0/PI**4)*betae**3* &
     	  (cv1t*(8.0*(fai1*u2+fai2*u1)-2.0*(fai_1*u2+fai2*u_1)+ &
         7.0*(fai0*u1+fai1*u0)+5.0*(fai0*u_1+fai_1*u0))+ &
         9.0*cv2t*(fai0*(u1+u_1)+(fai_1+fai1)*u0))

	ev=ev1+ev2+ev3

	return
end subroutine pair
 
!c..This routine calculates the neutrino emissivity of the 
!c..Urca processes.
!c..Input:
!c..      b1 dimensionless factor which characteristic 
!c..         the effect of neutrino trapping,
!c..      xkt temperature in unit of MeV, 
!c..      xnn the number density of free neutron,
!c..      xnp the number density of free proton,
!c..      xetae chemical potential of electron in unit of xkt,
!c..      xetaep chemical potential of positron in unit of xkt,
!c..Output:
!c..      evep_k the rate of electron capture in unit of ergs /s/cm^3,
!c..      even_k the rate of positron capture in unit of ergs /s/cm^3,
!c..      evn_k the rate of neutron decay in unit of ergs /s/cm^3,
!c..
subroutine urca(b1,xkt,xnn,xnp,xetae,xetaep,evep_k,even_k,evn_k)
	
!      external func1,func2,func3
      integer n
!      double precision func1,func2,func3
      double precision xkt,par(3),rate,a,b,fdxyz1,fdxyz2,fdxyz3,b1, &
       t11,xnn,xnp,q,xetae,xetaep,evn_k,evep_k,even_k

      b1=0.0d0

      a=xme/xkt
      b=(xmn-xmp)/xkt

      n=3
      t11=xkt*1.16059d-1 
      par(1)=xkt
      par(2)=xetae
      par(3)=b1

!c..n--->p+e+v~

      call gauss_legendre(func1,a,b,fdxyz1,par,n)

!c      evn_k=2.03813046d35*t11**6*xnn*fdxyz1
      evn_k=1.14473926d37*t11**6*xnn*fdxyz1

!c..p+e--->n+v

      a=(xmn-xmp)/xkt
      call gauss_laguerre(func2,a,1.0d0,fdxyz1,par,n)
      evep_k=1.14473926d37*t11**6*xnp*fdxyz1

!c..n+e--->p+v~

      a=xme/xkt
      par(2)=xetaep
      call gauss_laguerre(func3,a,1.0d0,fdxyz1,par,n)
      even_k=1.14473926d37*t11**6*xnn*fdxyz1

      return
    end subroutine urca


!c..........................................................
!c..n--->p+e+v~

    subroutine func1(x,par,n,fdyz)
      
      integer n
      double precision x,par(n),fdyz,xetae,xkt,betae,q,xe,factor1,b1,factor2
       

      xkt=par(1)
      xetae=par(2)
      b1=par(3)

      betae=xkt/xme
      xe=xetae-(x-1.0/betae)
!      if(xe.lt.5.0d0) then
             factor1=1.0/(exp(xe)+1.0)
!      else
!        factor1=exp(-xe)
!      endif


      q=(xmn-xmp)/xkt
      factor2=abs(1.0-b1/(exp(q-x)+1.0))

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(q-x)**3*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func1

!c..p+e--->n+v
      subroutine func2(x,par,n,fdyz)
      
      integer n
      double precision x,par(n),fdyz,xetae,xkt,betae,q,xe,factor1,factor2,b1
      
      xkt=par(1)
      xetae=par(2)
      b1=par(3)

      betae=xkt/xme
      xe=-xetae+(x-1.0/betae)
!      if(xe.lt.5.0d1) then
             factor1=1.0/(exp(xe)+1.0)
!      else
!        factor1=exp(-xe)
!      endif

      q=(xmn-xmp)/xkt
      factor2=abs(1.0-b1/(exp(x-q)+1.0))
      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(x-q)**3*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func2

!c..n+e--->p+v~
      subroutine func3(x,par,n,fdyz)
      
      integer n
      double precision x,par(n),fdyz,xetae,xkt,betae,q,xe,factor1,factor2,b1
      
      xkt=par(1)
      xetae=par(2)
      b1=par(3)

      betae=xkt/xme
      xe=(x-1.0/betae)-xetae
!      if(xe.lt.5.0d1) then
             factor1=1.0/(exp(xe)+1.0)
!      else
!        factor1=exp(-xe)
!      endif


      q=(xmn-xmp)/xkt
      factor2=abs(1.0-b1/(exp(q+x)+1.0))
      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(x+q)**3*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func3


!----------------------------------------------------------
    subroutine brem(nb,xkt,rate)
      
      double precision nb,xkt,rho10,t11,rate
      rho10=nb*1.66d5
      t11=xkt*1.16059d-1 
      rate=3.35d27*rho10**2.0*t11**5.5 
      return 
    end subroutine brem

!----------------------------------------------------------
    subroutine plasmon(xkt,xetae,rate)
      
      double precision xkt,xetae,rate,xmue,t11,etae,gp
      t11=xkt*1.16059d-1
    
      etae=xetae+0.510998910d0/xkt
      gp=5.565d-2*dsqrt((PI**2+3.0*etae**2)/3.0) 
      rate=1.5d32*t11**9*gp**6*exp(-gp)*(1.0+gp)*(2.0+gp**2/(1.0+gp))
      return
    end subroutine plasmon


  end module
