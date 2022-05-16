
module accretion 

use number_den
use const_mod
use neu_rate

implicit none


public::  &
     xnuc, &
     h_dep, &
     opt_dep, &
     b_neu, &
     c_neu

contains

       

!
!The mass fraction of nucleons (Popham et al. 1999, P.360
!see also Qian & Woosley 1996)
!


!-------------------------------------------------------------------------------
      double precision function xnuc(nb,xkt)
      
        double precision nb,xkt,rho10,t10,t11
        rho10=nb*1.66d5 ! density /10^10 [g/cm^3]
        t10=xkt*1.16059d0 ! T 
        t11=xkt*1.16059d-1 ! T 
        ! agnes      xnuc=30.97/rho10**(0.75d0)*t10**1.125d0*exp(-6.096d0/t10)
!        xnuc=29.55/rho10**(0.75d0)*t10**1.125d0*exp(-8.209d0/t10)
      	xnuc= 295.5d0/rho10**(0.75d0)*t11**1.125d0*exp(-0.8209d0/t11)

        if(xnuc.gt.1.0d0) xnuc=1.0d0
!       if(xnuc.lt.1.0d-50) xnuc=0.0d0
        if(xnuc.lt.0.0d0) xnuc=0.0d0
        return
      end function xnuc


! disk thickness in [cm] from model, need to calculate that separatelly
!------------------------------------------------------------------------------
      double precision function h_dep(HdFac,nb,p,e)
        
        double precision HdFac,nb,p,e,rho,cs,CC
!        h_dep=1.93d5*dsqrt(mdot/(6.0*PI*alpha*dsqrt(rho*p)))

        CC=2.99792458d10
        rho=nb*1.66d15
        p=p*1.6021772d33
        e=e*1.6021772d33

!        cs=dsqrt(p/rho)
        cs=CC*dsqrt(p/(rho*CC*CC + e))
        if(cs.gt.CC) cs=CC
        h_dep=HdFac*cs

        return
      end function h_dep

! optical depth for neutrino absorption, and neutrino scattering on n and p      
! we want that in harm
!------------------------------------------------------------------------------
      subroutine opt_dep(hd,b1,xkt,nnucl,xnn,xnp,xetae,xetaep,ta1,ta2,ts)
      
        double precision hd,xkt,nnucl,xnn,xnp,xetae,xetaep,ta1,ta2,ts, &
             t11,ev,ev1,ev2,ev3,evep_k,even_k,evn_k, &
             evbrem,evplasmon,qe,qo,xl,rho10,b1, &
             csn,csp,cs0,xme,cv,alph,tsn,tsp
        rho10=nnucl*1.66d5
        t11=xkt*1.16059d-1
        b1=0.0d0
        call pair(xkt,xetae,ev,ev1,ev2,ev3)
        call urca(b1,xkt,xnn,xnp,xetae,xetaep,evep_k,even_k,evn_k)
        call brem(nnucl,xkt,evbrem)
        call plasmon(xkt,xetae,evplasmon)
        qe=ev1+evep_k+even_k+evn_k+evbrem/3.0+evplasmon
!        qo=ev2*2.0+evbrem*2.0/3.0     ! contribution from mu and tau neutrinos
        qo=ev2+evbrem/3.0     ! contribution from mu and tau neutrinos
!       xl=5.67d39*4.0*7.0/8.0*t11*4  ! 4.0*7.0/8.0*sigma*T^4
        xl=5.67051d39*4.0*7.0/8.0*t11**4  ! 4.0*7.0/8.0*sigma*T^4
        ta1=qe*hd/xl
        ta2=qo*hd/xl
        
        cv=0.5+2.0*0.23
        alph=1.25d0
        xme=0.510998910d0
        csn=(1.0d0+5.0*alph**2)/24.0
        csp=(4.0d0*(cv-1.0)**2+5.0*alph**2)/24.0 
        
        tsn=csn*1.76d-5*xnn*hd*13.8*(xkt/xme)**2
        tsp=csp*1.76d-5*xnp*hd*13.8*(xkt/xme)**2
        ts=tsn+tsp
        return
      end subroutine opt_dep

! calculate bloking factors for taon and electron neutrinos 
!--------------------------------------------------------------------

      subroutine b_neu(ta1,ta2,ts,b1,b2)
        
        double precision ta1,ta2,ts,b1,b2,r3,te,tt
        te=ta1+ts
        tt=ta2+ts
        r3=dsqrt(3.0d0)
        b1=(te/2.0+1.0/r3)/(te/2.0+1.0/r3+1./3./ta1)
        b2=(tt/2.0+1.0/r3)/(tt/2.0+1.0/r3+1./3./ta2)
        if(ta1.le.0.0d0 .or. ta2.le.0.0d0) then
          b1=0.0d0
          b2=0.0d0
        endif
        return
      end subroutine b_neu

! 
!----------------------------------------------------
      subroutine c_neu(ta1,ta2,ts,c1,c2)
      
        double precision ta1,ta2,ts,c1,c2,r3,te,tt
        te=ta1+ts
        tt=ta2+ts
        r3=dsqrt(3.0d0)
        c1=4.0d0/3.0d0/(te/2.0+1.0d0/r3+1./3./ta1)
        c2=4.0d0/3.0d0/(tt/2.0+1.0d0/r3+1./3./ta2)
        if(ta1.le.0.0d0 .or. ta2.le.0.0d0) then
          c1=0.0d0
          c2=0.0d0
        endif
        return
      end subroutine c_neu


    end module accretion
