module number_den

use const_mod
use gauss_integral
use bisection
use mnewtf

implicit none


public:: &
     fxnn, &
     fxnp, &
     neta, &
     peta, &
     eeta, &
     lepton

contains

! number densities of fermions eq. (20): e+,e-,p,n


! number density of neutrons eq (20)
!-----------------------------------------------------------------
  double precision function fxnn(xkt,xetan)
    double precision xkt,xetan,betan,cnn,fs(3),fdeta,fdtheta,xmn
    
    xmn=939.56536d0
    betan=xkt/xmn
    cnn=dsqrt(2.0d0)/PI**2*(xmn/HC)**3*betan**1.5
    
    call dfermi(0.5d0,xetan,betan,fs(1),fdeta,fdtheta)
    call dfermi(1.5d0,xetan,betan,fs(2),fdeta,fdtheta)
    fxnn=cnn*(fs(1)+betan*fs(2))
    
    return
  end function fxnn
  

! number density of protons eq(20)
!-----------------------------------------------------------------
double precision function fxnp(xkt,xetap)
  double precision xkt,xetap,betap,cnp,fs(3),fdeta,fdtheta,xmp
  
  xmp=938.272013d0
  betap=xkt/xmp
  cnp=dsqrt(2.0d0)/PI**2*(xmp/HC)**3*betap**1.5
  
  call dfermi(0.5d0,xetap,betap,fs(1),fdeta,fdtheta)
  call dfermi(1.5d0,xetap,betap,fs(2),fdeta,fdtheta)
  fxnp=cnp*(fs(1)+betap*fs(2))

  return
end function fxnp





subroutine usrfun_neta(x,n,np,fjac,fvec,par,npar) 
  integer n,np,npar
  double precision x(n),par(npar) 
  double precision fjac(np,np),fvec(np) 
  double precision xkt,xnn
  double precision xetan,xetanp,xetanm
  double precision nn,nnp,nnm
  double precision f,fp,fm
  double precision xacc

  xkt=par(1)
  xnn=par(2)

  xacc=dabs(x(1))*1.0d-8
  if(xacc.lt.1.0d-50) then
    xacc=1.0d-50
  endif
  xetan=x(1)
  xetanp=xetan+xacc
  xetanm=xetan-xacc

  nn=fxnn(xkt,xetan)
  nnp=fxnn(xkt,xetanp)
  nnm=fxnn(xkt,xetanm)

  f=(xnn-nn)
  fp=(xnn-nnp)
  fm=(xnn-nnm)

  fjac(1,1)=(fp-fm)/(2.0d0*xacc);
  fvec(1)=-f

end subroutine usrfun_neta

!-----------------------------------------------------------------
double precision function neta(xkt,xnn)

  double precision xkt,xnn,tol
  double precision x(1),par(2),errf
!  double precision oneta

  x(1)=-10.0d0
  par(1)=xkt
  par(2)=xnn
  tol=1.0d-8*xnn
  if(tol.lt.1.0d-50) tol=1.0d-50
  call mnewt(200,x,1,1.0d-20,tol,usrfun_neta,par,2,errf) 

  
  neta=x(1)
!  oneta = neta_old(xkt,xnn)

!  print*,'neta mnewt, neta=', neta,'errf=',errf,'oneta=',oneta

end function neta


! czemu iteracjynie
!-----------------------------------------------------------------
double precision function neta_old(xkt,xnn)
  double precision xkt,xetan,xnn,nn
  double precision xetan1,xetan2,deta
  integer j
  double precision dx,xmid,fmid,xacc
  double precision rtbis1
  xetan1=-100.0d0
  deta=1.0d0
  
  xetan=xetan1
  nn=fxnn(xkt,xetan)
  do while(nn.le.xnn)
     nn=fxnn(xkt,xetan)
     xetan=xetan+deta
  enddo
  
  rtbis1=xetan-2.0*deta
  dx=deta
  xacc=dx*1.0d-10
  
  do j=1,40
     dx=dx*.5
     xmid=rtbis1+dx
     nn=fxnn(xkt,xmid)
     fmid=nn-xnn
     if(fmid.lt.0.) rtbis1=xmid
     if(dabs(dx).lt.xacc .or. fmid.eq.0.) exit
  enddo
  
  neta_old=rtbis1
  
end function neta_old




!-----------------------------------------------------------------
subroutine usrfun_peta(x,n,np,fjac,fvec,par,npar) 
  integer n,np,npar
  double precision x(n),par(npar) 
  double precision fjac(np,np),fvec(np) 
  double precision xkt,xnp
  double precision xetap,xetapp,xetapm
  double precision npt,npp,npm
  double precision f,fp,fm
  double precision xacc

  xkt=par(1)
  xnp=par(2)

  xacc=dabs(x(1))*1.0d-8
  if(xacc.lt.1.0d-50) then
    xacc=1.0d-50
  endif
  xetap=x(1)
  xetapp=xetap+xacc
  xetapm=xetap-xacc

  npt=fxnp(xkt,xetap)
  npp=fxnp(xkt,xetapp)
  npm=fxnp(xkt,xetapm)

  f=(xnp-npt)
  fp=(xnp-npp)
  fm=(xnp-npm)

  fjac(1,1)=(fp-fm)/(2.0d0*xacc);
  fvec(1)=-f

end subroutine usrfun_peta

!-----------------------------------------------------------------
double precision function peta(xkt,xnp)

  double precision xkt,xnp,tol
  double precision x(1),par(2),errf
!  double precision oneta

  x(1)=-10.0d0
  par(1)=xkt
  par(2)=xnp
  tol=1.0d-8*xnp
  if(tol.lt.1.0d-50) tol=1.0d-50
  call mnewt(200,x,1,1.0d-20,tol,usrfun_peta,par,2,errf) 

  
  peta=x(1)
!  opeta = neta_old(xkt,xnn)

!  print*,'peta mnewt, peta=', peta,'errf=',errf,'opeta=',opeta

end function peta



double precision function peta_old(xkt,xnp)
  double precision xkt,xetap,xnp,np
  double precision xetap1,xetap2,deta
  integer j
  double precision rtbis,dx,xmid,fmid,xacc
  
  xetap1=-100.0d0
  deta=1.0d0
  
  xetap=xetap1
  np=fxnp(xkt,xetap)
  do while(np.le.xnp)
     np=fxnp(xkt,xetap)
     xetap=xetap+deta
  enddo
  
  rtbis=xetap-2.0*deta
  dx=deta
  xacc=dx*1.0d-10
  
  do j=1,40
     dx=dx*.5
     xmid=rtbis+dx
     np=fxnp(xkt,xmid)
     fmid=np-xnp
     if(fmid.lt.0.) rtbis=xmid
     if(dabs(dx).lt.xacc .or. fmid.eq.0.) exit
  enddo
        
  peta_old=rtbis

end function peta_old

subroutine usrfun_eeta(x,n,np,fjac,fvec,par,npar) 
  integer n,np,npar
  double precision x(n),par(npar) 
  double precision fjac(np,np),fvec(np) 
  double precision xkt,xne,xnee,xnep
  double precision xetae,xetaep,xetaem
  double precision nl,nlp,nlm
  double precision f,fp,fm
  double precision xacc

  xkt=par(1)
  xne=par(2)

  xacc=dabs(x(1))*1.0d-8
  if(xacc.lt.1.0d-50) then
    xacc=1.0d-50
  endif
  xetae=x(1)
  xetaep=xetae+xacc
  xetaem=xetae-xacc

  call lepton(xkt,xetae,nl,xnee,xnep)
  call lepton(xkt,xetaep,nlp,xnee,xnep)
  call lepton(xkt,xetaem,nlm,xnee,xnep)

  f=(xne-nl)
  fp=(xne-nlp)
  fm=(xne-nlm)


  fjac(1,1)=(fp-fm)/(2.0d0*xacc);
  fvec(1)=-f

end subroutine usrfun_eeta


!-----------------------------------------------------------------
double precision function eeta(xkt,xne)

  double precision xkt,xne,tol
  double precision x(1),par(2),errf
!interface
!         subroutine usrfun_eeta(x,n,np,fjac,fvec,par,npar) 
!         integer n,np,npar
!         double precision x(n),par(npar) 
!         double precision fjac(np,np),fvec(np) 
!         end subroutine usrfun_eeta
!end interface 
 
!  print*,'eeta xkt=',xkt,'xne=',xne
  
  x(1)=0.510998910d0/xkt
  par(1)=xkt
  par(2)=xne
  tol=1.0d-8*xne
  if(tol.lt.1.0d-50) tol=1.0d-50
  call mnewt(200,x,1,1.0d-20,tol,usrfun_eeta,par,2,errf)
  eeta=x(1)

  

  if(eeta.gt.0.0d0) then
	eeta = 0.0d0
  endif


!  print*,'eeta mnewt, eeata=', eeta,'errf=',errf
!  eeta = eeta_old(xkt,xne)

end function eeta

!
!-----------------------------------------------------------------
double precision function eeta_old(xkt,xne)

  double precision xkt,xetae,nl,xnee,xnep
  double precision xetae1,xetae2,deta,xne
  integer j
  double precision rtbis,dx,xmid,fmid,xacc
  double precision x(1),par(2),errf
 
!  print*,'eeta xkt=',xkt,'xne=',xne
  
  xetae1=-0.511d0/xkt
  deta=0.1d0
  
  xetae=xetae1
  nl=0.0d0
!   print*,'eeta 1...'

  ! iteracja do liczenia xnee i xnep
  do while(nl.le.xne)
!   print*,'eeta 1.x...', xne, nl,'xetae=',xetae,xnee,xnep
     call lepton(xkt,xetae,nl,xnee,xnep)
     xetae=xetae+deta
  enddo
!  print*,'eeta 2...'

  rtbis=xetae1
  dx=xetae-xetae1
  xacc=dx*1.0d-10
 

! interacja do liczenia xnee i xnep
  do j=1,40
     dx=dx*.5
     xmid=rtbis+dx
     call lepton(xkt,xmid,nl,xnee,xnep)
     fmid=nl-xne
     if(fmid.lt.0.) rtbis=xmid
     if(dabs(dx).lt.xacc .or. fmid.eq.0.) exit
  enddo
 
  print*,'eeta end,   rtbis=', rtbis

  eeta_old=rtbis

end function eeta_old



! oblicza liczbe leptonow:e+ e-
!-----------------------------------------------------------------
subroutine lepton(xkt,xetae,nl,xnee,xnep)
        
  double precision xkt,xetae,nl,betae,cne,fs(3),fdeta,fdtheta, &
       xnee,xnep,etapos
  
  betae=xkt/xme
  cne=dsqrt(2.0d0)/PI**2*(xme/HC)**3*betae**1.5
  
  call dfermi(0.5d0,xetae,betae,fs(1),fdeta,fdtheta)
  call dfermi(1.5d0,xetae,betae,fs(2),fdeta,fdtheta)
  xnee=cne*(fs(1)+betae*fs(2)) ! equation (20) for electrons
  
  etapos=-xetae-2.0/betae

!  if(etapos.lt.-1.0d2) then
!     xnep=0.0d0
!  else
    call dfermi(0.5d0,etapos,betae,fs(1),fdeta,fdtheta)
    call dfermi(1.5d0,etapos,betae,fs(2),fdeta,fdtheta)
    xnep=cne*(fs(1)+betae*fs(2))   ! equation (20) for positrons
!  endif
  
  nl=xnee-xnep ! nominator of eq (21)
  return
end subroutine lepton

end module number_den
