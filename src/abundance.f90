! Calculates the abundances for a given rho, T and hd
! calculates neutriono cooling rates using Janiuk(2007) code
! modified Agnieszka code

! subroutine abundance(rho, T, hd, Qnu, ta1, ta2, ts)
 subroutine abundance(rho, T, HdFac, e, p, Qnu, tau, xnee, xnep)

  use const_mod
  use chemistry
  use accretion
  use number_den
  use eos
  use gauss_integral

  implicit none

  integer i
!  integer i,j,k
!  data xetan/-0.3145212d2/
!  SAVE xetan
!  integer, parameter :: N=60   ! N=30 in old version
!  double precision, parameter :: logT_min=8 !9 for 3 Msun
!  double precision, parameter :: logT_max=15 !14 for 3 Msun
!  double precision, parameter :: logrho_min=7 !9 for 3 Msun
!  double precision, parameter :: logrho_max=13 !14 for 3 Msun
!  double precision, parameter :: MBH = 10. !in MSUN units

!  integer, parameter :: N=30   ! N=30 in old version
!  double precision, parameter :: logT_min=9 !9 for 3 Msun
!  double precision, parameter :: logT_max=14 !14 for 3 Msun
!  double precision, parameter :: logrho_min=9 !9 for 3 Msun
!  double precision, parameter :: logrho_max=14 !14 for 3 Msun
!  double precision, parameter :: MBH = 3. !in MSUN units

  double precision xnb,xt,xnn,xnp,xnee,xnep,xetan,xetap,xetae
  double precision b1,b2,e,p,s,nhe,nnuc,nl,xne,xnuc0,xkt
  double precision cne,cnp,cnn,betae,betap,betan,q
  double precision x1,x2,x3,x4,x5,rad,ye,mdot,alpha,neurate,prad2
  double precision mdot0,alpha0,enucl,pnucl,snucl,prad,erad,srad,pneu,eneu,sneu
  double precision phe,ehe,she,p1,e1,rho0,hd,xetaep,ta1,ta2,ts,c1,c2,tau
  
!  double precision T,T11,rho,r
  double precision T,T11,rho
  double precision ev,ev1,ev2,ev3,evep_k,even_k,evn_k,evbrem,evplasmon,Qnu
!  double precision Lunit
  double precision dlogT,dlogrho
  double precision HdFac	
  character*40 output

!  print*,'Calculating EOS table ...'
!  print*,'abundance 1...'

  ! set the lenght scale for MBH 
!  Lunit=G*MBH*MSUN/CL/CL

!  dlogT=(logT_max-logT_min)/N
!  dlogrho=(logrho_max-logrho_min)/N

!100  format(3X,A8,18X,A6,20X,A8,17X,A3)


  ! open file to write down a table
!  output="cooling_tab.dat"
!  open(1,FILE=output,status='new')
!  write(1,100) 'log(rho)','log(T)','log(Qnu)','n_n'  
  
!  do j=1,N,1
!     do k=1,N,1
     !P=(gam-1)u
     !calculate T from the code using: P m_avg /rho_0 k where m_avg=1/5*(2*m_e+mp+mn+mhe) is aveaged mass 
  !in cgs units parameters lets assume, make table
!     T=10**(logT_min+k*dlogT)
     !T=1.e11     ! K  
!     rho=10**(logrho_min+j*dlogrho)
     !rho=1e12   ! g/cm^3
     
!     r=10.*Lunit ! cm where?
     xt=T
!     hd=0.5*r   ! cm, lets assume for simplicity, should be somehow calculated from the model
 
 !    hd=0.1*r ! for testing thickness
    

  ! unit conversion
     xnb=rho/(1.66d15)  ! number density of barions:p,n,e,He (number per fm^3) [1fm=1e-13 cm]
     t11=xt*1.0d-11     ! temperature in T11 units
     xkt=xt/1.16059d10  ! kT in [MeV] units !
 
     q=(-xmn+xmp+2.0d0*xme)/xkt  ! what is it?

!this is in common , used later in !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     betae=xkt/xme     ! dimentionless temperature of electrons
     ! not used here...
     betap=xkt/xmp     ! dimentionless temperature of protons
     betan=xkt/xmn     ! dimentionless temperature of neutrons
     
     cne=dsqrt(2.0d0)/PI**2*(xme/HC)**3*betae**1.5
     cnp=dsqrt(2.0d0)/PI**2*(xmp/HC)**3*betap**1.5
     cnn=dsqrt(2.0d0)/PI**2*(xmn/HC)**3*betan**1.5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
     xnuc0=xnuc(xnb,xkt)  ! calculate Xnuc-from temperature and rho eq (7)
     nnuc=xnb*xnuc0  
     nhe=(xnb-nnuc)/4.0d0 ! number density of He nucleons eq (6)
     
  ! calculate number densities and chemical potentials of barions
  ! using equalibrium reactions
  ! will get n_n, n+p, n_e-, n_e+ and eta_e-, eta_p, eta_n, eta_e+

!  write(*,*) 'initial Xnuc_0',xnuc0
!  print*,'abundance 2...'

!     if(xnuc0.gt.1.0d-9) then         
        ! input: number density xnb of barions xnb and xkt [MeV]
        ! output: xnn,xnp,xnee,xnep,xetan,xetap,xetae
        call chem_eq(xnb,xkt,xnn,xnp,xnee,xnep,xetan,xetap,xetae,HdFac)
!     else
!        xnn=0.0d0 ! xnuc small than no free nucleons p and n
!        xnp=0.0d0
!        xetan=-1.0d3
!        xetap=-1.0d3
!        xne=xnb/2.0d0
!        xetae=eeta(xkt,xne)
!     endif
!  print*,'abundance 4...'

     ! chemical potential for positrons
     xetaep=-xetae-2.0d0/betae
  
     ! what is that?
     ! input xkt, xetae,
     call lepton(xkt,xetae,nl,xnee,xnep) ! what does it do??? XXXXXXXXXXXXXX

!  print*,'abundance 5...'

  ! blocking factors for two neutrino types, initially zero, 
     ! need to calculate that later when all number densities are known
     b1=0.0d0
     b2=0.0d0
  
  ! calculates eos (pressure,internal energy and entropy density -I dont need all of that)
  ! for different spcies
    call eos_npe(xkt,xetan,xetap,xetae,xnn,xnp,xnee,xnep,enucl,pnucl,snucl)
    call eos_rad(xkt,erad,prad,srad)
    call eos_he(xkt,nhe,ehe,phe,she)
  
     ! input xkt, xetae, output cooling rates ->table: ev,ev1,ev2,ev3
  !   call pair(xkt,xetae,ev,ev1,ev2,ev3)
     
     ! input b1,xkt,xnn,xnp,xetae,xetaep, output3 cooling rates->table:evep_k,even_k,evn_k
  !   call urca(b1,xkt,xnn,xnp,xetae,xetaep,evep_k,even_k,evn_k)

     ! to samo mozna policzyc w harmie
  !   call brem(nnuc,xkt,evbrem)

     ! to tez bedzie w harmie
  !   call plasmon(xkt,xetae,evplasmon)
  
     ! cooling rates
     !  write(*,*) evplasmon,evbrem,evep_k,even_k,evn_k,ev
     
     ! table structure
     ! rho, T, ev, evep_k+even_k+evn_k, xnp,xnn,xnee,xnep,xhe, b1, b2
     

!I need to code, up to this points,rest of it calculated in harm
!*****************************************************************
!*****************************************************************
!*****************************************************************

  ! p1- partial pressure without neutrions
    p1=pnucl+prad+phe
    e1=enucl+erad+ehe
  ! new density? - what units
  ! nHe*mHe+nn*mn+np*mp
!     rho0=nhe*3755.6740+xnn*939.565+xnp*938.272
 
   ! iteration to get correct 
     do i=1,10
        call eos_neu(xkt,b1,b2,eneu,pneu,sneu)  
        p=p1+pneu
        e=e1+eneu
!        hd=0.5*r 
     !   hd=0.1*r ! testing thickness

        hd=h_dep(HdFac,xnb,p,e) ! disk thickness
        call opt_dep(hd,b1,xkt,nnuc,xnn,xnp,xetae,xetaep,ta1,ta2,ts) ! optical depts for scattering and absorption
        call b_neu(ta1,ta2,ts,b1,b2) ! traping parameter, lets put it into table too
     enddo
     
  ! call eos for neutrinos once again to get the neutrono e,p and s - I dont need that
     call eos_neu(xkt,b1,b2,eneu,pneu,sneu)

!   write(*,*) ta1,ta2,ts
!   write(*,*) 'n_n=  ',xnn
!   write(*,*) 'n_p=  ',xnp
!   write(*,*) 'n_e-= ',xnee
!   write(*,*) 'n_e+= ',xnep
!   write(*,*) 'n_He= ',nhe
!   write(*,*) 'n_b=',xnb,nhe+xnep+xnee+xnp+xnn

! chlodzenie zaklada ze dysk jest w przyblizeniu H=0.5r
     Qnu=7.0d0/8.0d0*sigmasb*T**4/0.75*( &
          1./(0.5*ta1+0.5*ts+0.557+0.3/ta1) + &
          1./(0.5*ta2+0.5*ts+0.557+0.3/ta2))

     tau=ts+ta1+ta2
       
!   write(1,*) log10(rho), log10(T), log10(Qnu) 
 
   !  write(*,*) 'b1,b2',b1,b2

   !  write(*,*) 'tau1,tau2,ts',ta1,ta2,ts

! dont need that below
! total: pressure, internal energy and entropy with corrected neutrino pressure
  p=pnucl+prad+pneu+phe
  e=enucl+erad+eneu+ehe
!!  s=snucl+srad+sneu+she

  ! convert units from SI to cgs?
  p=p*1.6021772d33
  e=e*1.6021772d33
!!  s=s/1.16059d10*1.6021772d33

!  print*,'enucl=',enucl*1.6021772d33
!  print*,'erad=',erad*1.6021772d33  
!  print*,'eneu=',eneu*1.6021772d33
!  print*,'ehe=',ehe*1.6021772d33

! does not know what is this, capture by neutrons??? dont know
!!  call c_neu(ta1,ta2,ts,c1,c2)
!!  neurate=(c1+c2)*5.670373d39*7.0d0/8.0d0*t11**4
!!  if(neurate.le.1.0d-50) neurate=0.0d0
  
!  prad2=prad*1.60218d33
!enddo
!enddo

!close(1)
!  print*,'abundance end.'

end subroutine abundance
