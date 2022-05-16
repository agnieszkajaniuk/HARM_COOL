
module chemistry


use const_mod
use gauss_integral
use accretion 
use number_den
use eos
use bisection



implicit none


public:: &
     chem_eq, &
     eqnucl, &
     nuclear, &
     funcnb, &
     neutrino, &
     func1k, &
     func1r, &
     func2k, &
     func2r, &
     func3k, &
     func3r
     

contains



!this routine computes the fractions of the particles under
!chemical equilibrium. 
!Input: 
!       nb the baryon number density (in unit of fm^{-3}), and
!       xkt the tmpertature (in unit of MeV). 
!Output: 
!       xnn the number density of neutron, 
!       xnp the number density of proton, 
!       xnee the number density of electron, 
!       xnep the number density of positron, 
!       xetan the chemical potential of neutron 
!            (mu/xkt, mu does not contain the rest mass), 
!       xetap the chemical potential of proton, and 
!       xetae the chemical potential of electron.
 subroutine chem_eq(xnb,xkt,xnn,xnp,xnee,xnep,xetan,xetap,xetae,hdf)

   double precision xnb,xnbl,xkt,xnn,xnp,xnee,xnep,xetan,xetap,xetae,hdf
   
   call eqnucl(xkt,xnb,xnn,xnp,xnee,xnep,xetan,xetap,xetae,hdf) ! solves some chemical equlibrium equation

   if(xnep.lt.1.0d-60) xnep=0.0d0
   return

 end subroutine chem_eq
 


!-----------------------------------------------------------------
!This routine dertermines the chemcial equilibrium above the neutron
!drip.
!-----------------------------------------------------------------
 subroutine eqnucl(xkt,nb,xnn,xnp,xnee,xnep,xetan,xetap,xetae,hdf)
        
   double precision xkt,nb,xp,xnnucl,y1,y2,yacc, &
        xnn,xnp,xnee,xnep,xetan,xetap,xetae,hd, &
        q,nl,betae,xne,xnuc0,xne0,xkt1,hdf,HdFac,nb0

   common /nucl/xnnucl,xkt1,xne0,HdFac,nb0


   y1=0.0d0
   y2=1.0d0
   xnuc0=xnuc(nb,xkt)

   nb0=nb
   xnnucl=nb*xnuc0
   xne0=dabs(1.0d0-xnuc0)*nb/2.0d0

   HdFac=hdf
   xkt1=xkt
   betae=xkt/0.510998910d0
   yacc=1.0d-10

!   if(xnnucl.le.0.0) then
!      xp=0.0d0
!      xnn=0.0d0
!      xnp=0.0d0
!   else
      ! for a given temperature and density
      ! bisection method to calculate the n and p ratio, using reaction balance (Eq. 17)-> xp=np/nb  
!      xp=rtbis(funcnb,y1,y2,yacc)    
      xp=rtbis_funcnb(y1,y2,yacc)    
      call nuclear(xkt,xnnucl,xp,xetan,xetap) ! calculate potentials for np and ne
      xnp=xnnucl*xp
      xnn=xnnucl-xnp
      xne=xnp+xne0
!   endif

   xetae=eeta(xkt,xne)
   call lepton(xkt,xetae,nl,xnee,xnep)  ! calculate number of leptions e+ e-
   return

 end subroutine eqnucl


!-----------------------------------------------------------------
 subroutine nuclear(xkt,xnb,xp,xetan,xetap)
   double precision xkt,xnb,xp,xetan,xetap,xnp,xnn
   
   
   xnp=xnb*xp
   xnn=xnb-xnp
   if(xnp.le.0.0d0) then
!      xetap=-1.0d5
      xnp=0.0d0
!   else
!      xetap=peta(xkt,xnp)
   endif
   xetap=peta(xkt,xnp)
   
   if(xnn.le.0.0d0) then
!      xetan=-1.0d5
      xnn=0.0d0
!   else
!      xetan=neta(xkt,xnn)
   endif
   xetan=neta(xkt,xnn)
   
 end subroutine nuclear



! used to calculate the xp in bisection method, balance beteen 6 reactions eq. 11-16
!------------------------------------------------------------------------------

 function funcnb(xp)

   integer i
   double precision funcnb,xp,xnnucl,xkt,q,xetan,xetap,xetae,nl, &
        xnee,xnep,betae,difrate,xne,xnn,xnp,xetaep,pneu,sneu, &
        e,p,s,mdot0,alpha0,hd,b1,ta1,ta2,ts,b2,er,pr,sr,eneu, &
        e1,p1,xne0,nhe,ehe,phe,she,rho,HdFac,nb0
 
   common /nucl/xnnucl,xkt,xne0,HdFac,nb0

   q=(-xmn+xmp+2.0*xme)/xkt ! some neutrality conditions
   betae=xkt/xme
   call nuclear(xkt,xnnucl,xp,xetan,xetap)
   nhe=xne0/2.0
   xnp=xnnucl*xp
   xnn=xnnucl-xnp
   xne=xnnucl*xp+xne0
   xetae=eeta(xkt,xne)
   xetaep=-xetae-2.0/betae
   
   b1=0.0d0
   b2=0.0d0
!.ij???
   call lepton(xkt,xetae,nl,xnee,xnep)
   
   call eos_npe(xkt,xetan,xetap,xetae,xnn,xnp,xnee,xnep,e,p,s)
   call eos_rad(xkt,er,pr,sr)
   call eos_he(xkt,nhe,ehe,phe,she)
!   rho=nhe*3755.6740+xnn*939.565+xnp*938.272
   p1=p+pr+phe
   e1=e+er+ehe

   do i=1,10 ! to get b1 and b2 iteration
      call eos_neu(xkt,b1,b2,eneu,pneu,sneu)
      p=p1+pneu
      e=e1+eneu
      hd=h_dep(HdFac,nb0,p,e)
      call opt_dep(hd,b1,xkt,xnnucl,xnn,xnp,xetae,xetaep,ta1,ta2,ts)
      call b_neu(ta1,ta2,ts,b1,b2)
   enddo

   !        print*,'b1,b2=',b1,b2
   call neutrino(b1,xkt,xnn,xnp,xetae,difrate) ! calculate balance between 6 (3-3) reactions, n-p balance
   funcnb=difrate
   return
 end function funcnb





!-----------------------------------------------------------------

 subroutine neutrino(b1,xkt,xnn,xnp,xetae,difrate)
   integer n
   double precision xkt,par(3),rate,a,b,fdxyz1,fdxyz2,fdxyz3, &
        t11,xnn,xnp,q,xetae,xetaep,evn_k,evep_k,even_k, &
        difrate,b1,fdxyz1r,fdxyz2r,fdxyz3r
   
   xetaep=-xetae-2.0*xme/xkt

   a=xme/xkt
   b=(xmn-xmp)/xkt
   
   n=3
   t11=xkt*1.16059d-1 
   par(1)=xkt
   par(2)=xetae
   par(3)=b1

!.ij
   if(xnn.lt.0.0d0) then
      xnn=0.0d0
   endif
   if(xnp.lt.0.0d0) then
      xnp=0.0d0
   endif
	  
! Gamma in paper A1-A6
!n--->p+e+v~

   call gauss_legendre(func1k,a,b,fdxyz1,par,n)
   call gauss_legendre(func1r,a,b,fdxyz1r,par,n)

   evn_k=2.03813046d35*t11**5*xnn*(fdxyz1+fdxyz1r)

!p+e--->n+v

   a=(xmn-xmp)/xkt
   call gauss_laguerre(func2k,a,1.0d0,fdxyz1,par,n)
   call gauss_laguerre(func2r,a,1.0d0,fdxyz1r,par,n)
   evep_k=2.03813046d35*t11**5*xnp*(fdxyz1+fdxyz1r)

!n+e--->p+v~

   a=xme/xkt
   par(2)=xetaep
   call gauss_laguerre(func3k,a,1.0d0,fdxyz1,par,n)
   call gauss_laguerre(func3r,a,1.0d0,fdxyz1r,par,n)
   even_k=2.03813046d35*t11**5*xnn*(fdxyz1+fdxyz1r)
   
   if (abs(evep_k).lt.1.0d-80) then
     difrate=(evep_k-even_k-evn_k)/(abs(evep_k)+1.0d-80)
   else
     difrate=(evep_k-even_k-evn_k)/abs(evep_k)
   endif
   return

 end subroutine neutrino

!-----------------------------------------------------------------
! funkcje podcalkowe dla roznych reakcji 
!-----------------------------------------------------------------
!n--->p+e+v~
 subroutine func1k(x,par,n,fdyz)
   
   integer n
   double precision x,par(n),fdyz,xetae,xkt,betae,q,xe,factor1,factor2,b1


   xkt=par(1)
   xetae=par(2)
   b1=par(3)
   
   betae=xkt/xme
      xe=xetae-(x-1.0/betae)
!      if(xe.lt.5.0d1) then
         factor1=1.0/(exp(xe)+1.0)
!      else
!        factor1=exp(-xe)
!      endif

      q=(xmn-xmp)/xkt

      factor2=dabs(1.0-b1/(exp(q-x)+1.0))

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(q-x)**2*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
      end subroutine func1k
!-----------------------------------------------------------------
!p+e+v~--->n
      subroutine func1r(x,par,n,fdyz)
      
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

      factor2=abs(b1/(exp(q-x)+1.0))

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(q-x)**2*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func1r
!-----------------------------------------------------------------
!p+e--->n+v
      subroutine func2k(x,par,n,fdyz)
      
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

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(x-q)**2*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func2k

!-----------------------------------------------------------------
!n+v--->p+e
      subroutine func2r(x,par,n,fdyz)
      
      integer n
      double precision x,par(n),fdyz,xetae,xkt,betae,q,xe,factor1,factor2,b1


      xkt=par(1)
      xetae=par(2)
      b1=par(3)

      betae=xkt/xme
      xe=+xetae-(x-1.0/betae)
!      if(xe.lt.5.0d1) then
             factor1=1.0/(exp(xe)+1.0)
!      else
!        factor1=exp(-xe)
!      endif

      q=(xmn-xmp)/xkt

      factor2=abs(b1/(exp(x-q)+1.0))

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(x-q)**2*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func2r

!-----------------------------------------------------------------
!p+v~--->n+e
      subroutine func3r(x,par,n,fdyz)
      
      integer n
      double precision x,par(n),fdyz,xetae,xkt,betae,q,xe,factor1,factor2,b1


      xkt=par(1)
      xetae=par(2)
      b1=par(3)

      betae=xkt/xme
      xe=xetae-(x-1.0/betae)
!      if(xe.lt.5.0d1) then
             factor1=1.0/(exp(xe)+1.0)
!      else
!        factor1=exp(-xe)
!      endif

      q=(xmn-xmp)/xkt
    
      factor2=abs(b1/(exp(x+q)+1.0))

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(x+q)**2*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func3r

!-----------------------------------------------------------------
!n+e--->p+v~
      subroutine func3k(x,par,n,fdyz)

      
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
    
      factor2=abs(1.0-b1/(exp(x+q)+1.0))

      fdyz=x*dsqrt(x*x-1.0/(betae*betae))*(x+q)**2*factor1*factor2
      if(fdyz.lt.0.0d0) print*,'fdyz=',fdyz
      
      return
    end subroutine func3k


subroutine usrfun_rtbis_funcnb(x,n,np,fjac,fvec,par,npar) 
!  use chemistry, only: funcnb

  integer n,np,npar
  double precision x(n),par(npar) 
  double precision fjac(np,np),fvec(np) 
  double precision reqv
  double precision xmid,xmidp,xmidm
  double precision nn,nnp,nnm
  double precision f,fp,fm
  double precision xacc

  reqv=par(1)

  xacc=dabs(x(1))*1.0d-8
  if(xacc.lt.1.0d-50) then
     xacc=1.0d-50
  endif
  xmid=x(1)
  xmidp=xmid+xacc
  xmidm=xmid-xacc

  nn=funcnb(xmid)
  nnp=funcnb(xmidp)
  nnm=funcnb(xmidm)

  f=(reqv-nn)
  fp=(reqv-nnp)
  fm=(reqv-nnm)

  fjac(1,1)=(fp-fm)/(2.0d0*xacc);
  fvec(1)=-f

end subroutine usrfun_rtbis_funcnb

!-----------------------------------------------------------------
double precision function rtbis_funcnb(X1,X2,XACC)
  use mnewtf

  double precision X1,X2,XACC
  double precision x(1),par(1),errf,tol
  double precision xnnucl,xkt,xne0,HdFac,nb0
  common /nucl/xnnucl,xkt,xne0,HdFac,nb0

  x(1)=(X2-X1)/2.0
  par(1)=0.0
  tol=1.0d-8*xnnucl
  if(tol.lt.1.0d-50) tol=1.0d-50
  call mnewt(200,x,1,1.0d-20,tol,usrfun_rtbis_funcnb,par,1,errf) 
  
  if(x(1).lt.X1) then
    x(1)=X1
  endif
  if(x(1).gt.X2) then
    x(1)=X2
  endif
  
  rtbis_funcnb=x(1)

end function rtbis_funcnb



  end module chemistry
