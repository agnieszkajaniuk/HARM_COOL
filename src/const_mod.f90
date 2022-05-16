! constants c.g.s.

module const_mod 
  integer, parameter :: single = 4
  integer, parameter :: si = single
  integer, parameter :: double = 8
  integer, parameter :: dbl = double
  integer, parameter :: quad = 8
  integer, parameter :: qd = quad
  integer, parameter :: oct = 8
  integer, parameter :: oc = oct
  integer, parameter :: logic = 1
  integer, parameter :: lg = logic

  public
  real(8), parameter :: pi=3.1415926535897932
  real(8), parameter :: pi__2=1.570796327
  real(8), parameter :: pi_3__2=4.71238898
  real(8), parameter :: pi_4=12.56637061
  real(8), parameter :: pi_2=6.28318531
  real(8), parameter :: kb=1.380658E-16
  real(8), parameter :: CL=29979245800.0
  real(8), parameter :: CL2=9.0e20
  real(8), parameter :: h=6.6260755E-27
  real(8), parameter :: sigmasb=5.6e-5
  real(8), parameter :: sigmath=6.65248e-25
  real(8), parameter :: miu=0.5
  real(8), parameter :: miue=1.14
  real(8), parameter :: miup=1.23
  real(8), parameter :: miui=1.23
  real(8), parameter :: ee=4.8e-10
  real(8), parameter :: ee2=23.e-20
  real(8), parameter :: h__mc2=8.09e-21
  real(8), parameter :: kb__mc2=1.69e-10
  real(8), parameter :: G=6.67e-8
  real(8), parameter :: MSUN=1.99e33
  real(8), parameter :: alphaf=0.00729
  real(8), parameter :: HC=197.3269718e0
  real(8), parameter :: me=9.1093897E-28
  real(8), parameter :: mp=1.6726231E-24 
  real(8), parameter :: mn=1.6749272E-24 
  real(8), parameter :: mh=1.6726231E-24
  real(8), parameter :: xme=0.510998910d0   ! MeV
  real(8), parameter :: xmp=938.272013d0 ! MeV
  real(8), parameter :: xmn=939.56536d0 ! MeV  
  real(8), parameter :: SMALL=1e-30  



end module const_mod
