module eos

use const_mod
use gauss_integral

implicit none

public:: &
     eos_npe, &
     eos_rad, &
     eos_neu, &
     eos_he

contains
!c..
!c..This routine computes the EOS of the free neutron, proton, 
!c..and electron-positron pair gas in beta equilibrium.
!c..Input:  xkt teperature in unit of MeV,
!c..        xetan chemical potential of neutron in unit of xkt, 
!c..              which does not contain the rest energy of the particle,
!c..        xetap chemical potential of proton,
!c..        xetae chemical potential of electron. The chemical potential
!c..                 of positron is -xetae-2.0/betae.
!c..        xnn the number density of neutron,
!c..        xnp the number density of proton,
!c..        xnee the number density of electron, and
!c..        xnep the number density of positron.
!c..Output: e the energy density in unit of MeV fm^{-3},
!c..        p the pressure in unit of MeV fm^{-3}, and
!c..        s the entropy density in unit of k_B*fm^{-3}
!c..
        subroutine eos_npe(xkt,xetan,xetap,xetae,xnn,xnp,xnee,xnep,e,p,s)
        
	double precision xkt,xetan,xetap,xetae,xnn,xnp,xnee,xnep,e,p,s, & 
         q,betae,betap,betan,cne,cnp,cnn,cpe,cpp,cpn, &
         cee,cep,cen,pee,pep,pn,pp,eee,eep,en,ep,sn,sp,fs(3),fdeta, &
         fdtheta,see,sep,xetaep
	

	q=(xmn-xmp-xme)/xkt

	betae=xkt/xme
	betap=xkt/xmp
	betan=xkt/xmn

	cne=dsqrt(2.0d0)/PI**2*(xme/HC)**3*betae**1.5
	cnp=dsqrt(2.0d0)/PI**2*(xmp/HC)**3*betap**1.5
	cnn=dsqrt(2.0d0)/PI**2*(xmn/HC)**3*betan**1.5

	cpe=2.0/3.0*dsqrt(2.0d0)/PI**2*(xme/HC)**3*xme*betae**2.5
	cpp=2.0/3.0*dsqrt(2.0d0)/PI**2*(xmp/HC)**3*xmp*betap**2.5
	cpn=2.0/3.0*dsqrt(2.0d0)/PI**2*(xmn/HC)**3*xmn*betan**2.5

	cee=3.0/2.0*cpe
	cep=3.0/2.0*cpp
	cen=3.0/2.0*cpn

!----------Electrons---------------------------------
	call dfermi(1.5d0,xetae,betae,fs(2),fdeta,fdtheta)
	call dfermi(2.5d0,xetae,betae,fs(3),fdeta,fdtheta)
        pee=cpe*(fs(2)+0.5*betae*fs(3))
        eee=cee*(fs(2)+betae*fs(3))
	see=(eee+pee)/xkt-xnee*(xetae+1.0/betae)
	if(see.le.0.0d0) see=0.0d0

!----------Positrons---------------------------------
	if(xnee.le.0.0d0) then 
	  pep=0.0d0
	  eep=0.0d0
	  sep=0.0d0
	else
	  xetaep=-xetae-2.0/betae
	  call dfermi(1.5d0,xetaep,betae,fs(2),fdeta,fdtheta)
	  call dfermi(2.5d0,xetaep,betae,fs(3),fdeta,fdtheta)
          pep=cpe*(fs(2)+0.5*betae*fs(3))
          eep=cee*(fs(2)+betae*fs(3))
	  sep=(eep+pep)/xkt-xnep*(xetaep+1.0/betae)
	endif
	if(sep.le.0.0d0) sep=0.0d0

!----------Protons---------------------------------
	call dfermi(1.5d0,xetap,betap,fs(2),fdeta,fdtheta)
	call dfermi(2.5d0,xetap,betap,fs(3),fdeta,fdtheta)
        pp=cpp*(fs(2)+0.5*betap*fs(3))
        ep=cep*(fs(2)+betap*fs(3))
	sp=(ep+pp)/xkt-xnp*(xetap+1.0/betap)
	if(sp.le.0.0d0) sp=0.0d0

!----------Neutrons---------------------------------
	call dfermi(1.5d0,xetan,betan,fs(2),fdeta,fdtheta)
	call dfermi(2.5d0,xetan,betan,fs(3),fdeta,fdtheta)
        pn=cpn*(fs(2)+0.5*betan*fs(3))
        en=cen*(fs(2)+betan*fs(3))
	sn=(en+pn)/xkt-xnn*(xetan+1.0/betan)
	if(sn.le.0.0d0) sn=0.0d0
	
!----------Total------------------------------------
	p=pee+pep+pp+pn
	e=eee+eep+ep+en
	s=see+sep+sp+sn

!  print*,'xetae=',xetae
!  print*,'betae=',betae

!  print*,'eee=',eee*1.6021772d33
!  print*,'eep=',eep*1.6021772d33
!  print*,'ep=',ep*1.6021772d33
!  print*,'en=',en*1.6021772d33

end subroutine eos_npe

!c..
!c..This routine computes the EOS of the radiation field
!c..Input:  xkt teperature in unit of MeV
!c..Output: e the energy density in unit of MeV fm^{-3},
!c..        p the pressure in unit of MeV fm^{-3}, and
!c..        s the entropy density in unit of k_B*fm^{-3}
!.

subroutine eos_rad(xkt,e,p,s)
       
       double precision xkt,e,p,s
       
       p=PI*PI/15.0*xkt*(xkt/HC)**3/3.0
       e=3.0d0*p
       s=4.0d0*p/xkt

       return 
     end subroutine eos_rad

!c..
!c..This routine computes the EOS of the trapped neutrinos
!c..Input:  xkt teperature in unit of MeV
!c..        b1 (without unit)
!c..        b2 (without unit)
!c..Output: e the energy density in unit of MeV fm^{-3},
!c..        p the pressure in unit of MeV fm^{-3}, and
!c..        s the entropy density in unit of k_B*fm^{-3}
!c..
       subroutine eos_neu(xkt,b1,b2,e,p,s)
       
       double precision xkt,e,p,s,b1,b2
     
       
       p=7.0d0/8.0d0*(b1+b2)*PI*PI/15.0*xkt*(xkt/HC)**3/3.0
       e=3.0d0*p
       s=4.0d0*p/xkt
       
       return 
     end subroutine eos_neu

!
!This routine computes the EOS of Helium
!Input:  xkt teperature in unit of MeV
!        nhe number density of helium in unit of fm^{-3}
!Output: e the energy density in unit of MeV fm^{-3},
!        p the pressure in unit of MeV fm^{-3}, and
!        s the entropy density in unit of k_B*fm^{-3}
!
       subroutine eos_he(xkt,nhe,e,p,s)
       
       double precision xkt,nhe,e,p,s,mhe
       
       mhe=3755.6740d0
       if(nhe.le.0.0d0) then
          e=0.0d0
          p=0.0d0
          s=0.0d0
       else
          e=1.5d0*nhe*xkt
          p=nhe*xkt
          s=nhe*(2.5d0+1.5d0*log(mhe*xkt/HC**2/2.0/PI)-log(nhe))
!          s=abs(s)
         if(s.le.0.0d0) then
           e=0.0d0
           p=0.0d0
           s=0.0d0
         endif
       endif
       return
     end subroutine eos_he
       
   end module eos
