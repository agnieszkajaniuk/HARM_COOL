

!routine dfermi gets the fermi-dirac functions and their derivaties
!routine fdfunc1 forms the integrand of the fermi-dirac functions
!routine fdfunc2 same as fdfunc but with the change of variable z**2=x
!routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
!routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
!routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
!routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
!routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
!routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
!routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
!routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy





      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta)

!this routine computes the fermi-dirac integrals of 
!index dk, with degeneracy parameter eta and relativity parameter theta.
!input is dk the double precision index of the fermi-dirac function,
!eta the degeneracy parameter, and theta the relativity parameter.
!the output is fd is computed by applying three 10-point 
!gauss-legendre and one 10-point gauss-laguerre rules over
!four appropriate subintervals. the derivative with respect to eta is
!output in fdeta, and the derivative with respct to theta is in fdtheta.
!within each subinterval the fd kernel.
!
!this routine delivers at least 9, and up to 32, figures of accuracy
!
!reference: j.m. aparicio, apjs 117, 632 1998


!declare
      external         fdfunc1,fdfunc2
      double precision dk,eta,theta,fd,fdeta,fdtheta,
                      d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3, &
                      eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3), &
                      res1,dres1,ddres1,res2,dres2,ddres2, &
                      res3,dres3,ddres3,res4,dres4,ddres4 &


!   parameters defining the location of the breakpoints for the
!   subintervals of integration:
      data d   / 3.3609d 0 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d 0 /
      data b1  / 1.1418d 0 /
      data c1  / 2.9826d 0 /
      data a2  / 3.7601d 0 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d 1 /
      data e2  / 1.0056d 0 /
      data a3  / 7.5669d 0 /
      data b3  / 1.1695d 0 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d 0 /
      data e3  /-1.2819d-1 /


!   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


!   definition of xi:
      eta1=sg*(eta-d)
      if (eta1.le.5.d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

!   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2)
     +  /(1.d0+c1*xi)
      x2=(a2  +b2*xi+c2*d2*xi2)
     +  /(1.d0+e2*xi+c2*   xi2)
      x3=(a3  +b3*xi+c3*d3*xi2)
     +  /(1.d0+e3*xi+c3*   xi2)

!   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=dsqrt(s1)

!   quadrature integrations: 

! 9 significant figure accuracy
      call dqleg010(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
      call dqleg010(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
      call dqleg010(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
      call dqlag010(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)


! 14 significant figure accuracy
      call dqleg020(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
      call dqleg020(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
      call dqleg020(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
      call dqlag020(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)


! 18 significant figure accuracy
!      call dqleg040(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
!      call dqleg040(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
!     call dqleg040(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
!     call dqlag040(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

! 32 significant figure accuracy
!      call dqleg080(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
!      call dqleg080(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
!      call dqleg080(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
!      call dqlag080(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)


!sum the contributions
      fd      = res1 + res2 + res3 + res4
      fdeta   = dres1 + dres2 + dres3 + dres4
      fdtheta = ddres1 + ddres2 + ddres3 + ddres4
      return
    end subroutine dfermi


!--------------------------------------------------------------------------

      subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta)
!
!forms the fermi-dirac integrand and its derivatives with eta and theta.
!on input x is the integration variable, par(1) is the double precision 
!index, par(2) is the degeneravy parameter, and par(3) is the relativity
!parameter. on output fd is the integrand, fdeta is the derivative
!with respect to eta, and fdtheta is the derivative with respect to theta.
!
!declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta, &
           factor,dxst,denom,denom2,xdk,xdkp1

!initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xdkp1 = x * xdk
      dxst  = dsqrt(1.0d0 + 0.5d0*x*theta)

!   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor  = exp(x-eta)
       denom   = factor + 1.0d0
       fd      = xdk * dxst / denom
       fdeta   = fd * factor / denom 
       denom2  = 4.0d0 * dxst * denom
       fdtheta = xdkp1 / denom2

      else
       factor   = exp(eta-x)
       fd       = xdk * dxst * factor
       fdeta    = fd
       denom2   = 4.0d0 * dxst
       fdtheta  = xdkp1/denom2 * factor
      endif

      return
    end subroutine fdfunc1


!-----------------------------------------------------------

      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta)
!
!forms the fermi-dirac integrand and its derivatives with eta and theta,
!when the z**2=x variable change has been made.
!on input x is the integration variable, par(1) is the double precision 
!index, par(2) is the degeneravy parameter, and par(3) is the relativity
!parameter. on output fd is the integrand, fdeta is the derivative
!with respect to eta, and fdtheta is the derivative with respect to theta.
!
!declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta, &
                      factor,dxst,denom,denom2,xdk,xdkp1,xsq

      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xdkp1 = xsq * xdk
      dxst  = dsqrt(1.0d0 + 0.5d0 * xsq * theta)

!   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.d2) then
       factor  = exp(xsq - eta)
       denom   = factor + 1.0d0
       fd      = 2.0d0 * xdk * dxst/denom 
       fdeta   = fd * factor/denom
       denom2  = 4.0d0 * dxst * denom
       fdtheta = 2.0d0 * xdkp1/denom2

      else
       factor  = exp(eta - xsq)
       fd      = 2.0d0 * xdk * dxst * factor
       fdeta   = fd 
       denom2  = 4.0d0 * dxst
       fdtheta = 2.0d0 * xdkp1/denom2 * factor
      endif

      return
    end subroutine fdfunc2



!-------------------------------------------------------------------

      subroutine dqleg010(f,a,b,result,dresult,ddresult,par,n)
!
!10 point gauss-legendre rule for the fermi-dirac function and
!its derivatives with respect to eta and theta.
!on input f is the name of the subroutine containing the integrand,
!a is the lower end point of the interval, b is the higher end point,
!par is an array of constant parameters to be passed to subroutine f,
!and n is the length of the par array. on output result is the 
!approximation from applying the 10-point gauss-legendre rule,
!dresult is the derivative with respect to eta, and ddresult is the
!derivative with respect to theta.
!
!note: since the number of nodes is even, zero is not an abscissa.
!
!declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n), &
           absc1,absc2,center,hlfrun,wg(5),xg(5), &
           fval1,dfval1,ddfval1,fval2,dfval2,ddfval2

! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 20-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 20-point rule
! wg     - weights of the 20-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.48874338981631210884826001129719984d-1 /
      data xg (  2) /   4.33395394129247190799265943165784162d-1 /
      data xg (  3) /   6.79409568299024406234327365114873575d-1 /
      data xg (  4) /   8.65063366688984510732096688423493048d-1 /
      data xg (  5) /   9.73906528517171720077964012084452053d-1 /

      data wg (  1) /   2.95524224714752870173892994651338329d-1 /
      data wg (  2) /   2.69266719309996355091226921569469352d-1 /
      data wg (  3) /   2.19086362515982043995534934228163192d-1 /
      data wg (  4) /   1.49451349150580593145776339657697332d-1 /
      data wg (  5) /   6.66713443086881375935688098933317928d-2 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 10-point gauss formula

      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,5
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
    end subroutine dqleg010


!------------------------------------------------------------------

      subroutine dqleg020(f,a,b,result,dresult,ddresult,par,n)
!
!20 point gauss-legendre rule for the fermi-dirac function and
!its derivatives with respect to eta and theta.
!on input f is the name of the subroutine containing the integrand,
!a is the lower end point of the interval, b is the higher end point,
!par is an array of constant parameters to be passed to subroutine f,
!and n is the length of the par array. on output result is the 
!approximation from applying the 20-point gauss-legendre rule,
!dresult is the derivative with respect to eta, and ddresult is the
!derivative with respect to theta.
!
!note: since the number of nodes is even, zero is not an abscissa.
!
!declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n), &
                      absc1,absc2,center,hlfrun,wg(10),xg(10), &
                      fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


! the abscissae and weights are given for the interval (-1,1).
! xg     - abscissae of the 20-point gauss-legendre rule
!          for half of the usual run (-1,1), i.e.
!          the positive nodes of the 20-point rule
! wg     - weights of the 20-point gauss rule.
!
! abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.65265211334973337546404093988382110d-2 /
      data xg (  2) /   2.27785851141645078080496195368574624d-1 /
      data xg (  3) /   3.73706088715419560672548177024927237d-1 /
      data xg (  4) /   5.10867001950827098004364050955250998d-1 /
      data xg (  5) /   6.36053680726515025452836696226285936d-1 /
      data xg (  6) /   7.46331906460150792614305070355641590d-1 /
      data xg (  7) /   8.39116971822218823394529061701520685d-1 /
      data xg (  8) /   9.12234428251325905867752441203298113d-1 /
      data xg (  9) /   9.63971927277913791267666131197277221d-1 /
      data xg ( 10) /   9.93128599185094924786122388471320278d-1 /

      data wg (  1) /   1.52753387130725850698084331955097593d-1 /
      data wg (  2) /   1.49172986472603746787828737001969436d-1 /
      data wg (  3) /   1.42096109318382051329298325067164933d-1 /
      data wg (  4) /   1.31688638449176626898494499748163134d-1 /
      data wg (  5) /   1.18194531961518417312377377711382287d-1 /
      data wg (  6) /   1.01930119817240435036750135480349876d-1 /
      data wg (  7) /   8.32767415767047487247581432220462061d-2 /
      data wg (  8) /   6.26720483341090635695065351870416063d-2 /
      data wg (  9) /   4.06014298003869413310399522749321098d-2 /
      data wg ( 10) /   1.76140071391521183118619623518528163d-2 /


!           list of major variables
!           -----------------------
!
!           absc   - abscissa
!           fval*  - function value
!           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
    end function value



!------------------------------------------------------------------

      subroutine dqleg040(f,a,b,result,dresult,ddresult,par,n)
c..
c..40 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 40-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(20),xg(20),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 40-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 40-point rule
c wg     - weights of the 40-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.8772417506 0508219331 9344402462 32946 d -2 /
      data xg (  2) /   1.1608407067 5255208483 4512844080 24113 d -1 /
      data xg (  3) /   1.9269758070 1371099715 5168520651 49894 d -1 /
      data xg (  4) /   2.6815218500 7253681141 1843448085 96183 d -1 /
      data xg (  5) /   3.4199409082 5758473007 4924811791 94310 d -1 /
      data xg (  6) /   4.1377920437 1605001524 8797458037 13682 d -1 /
      data xg (  7) /   4.8307580168 6178712908 5665742448 23004 d -1 /
      data xg (  8) /   5.4946712509 5128202075 9313055295 17970 d -1 /
      data xg (  9) /   6.1255388966 7980237952 6124502306 94877 d -1 /
      data xg ( 10) /   6.7195668461 4179548379 3545149614 94109 d -1 /
      data xg ( 11) /   7.2731825518 9927103280 9964517549 30548 d -1 /
      data xg ( 12) /   7.7830565142 6519387694 9715455064 94848 d -1 /
      data xg ( 13) /   8.2461223083 3311663196 3202306660 98773 d -1 /
      data xg ( 14) /   8.6595950321 2259503820 7818083546 19963 d -1 /
      data xg ( 15) /   9.0209880696 8874296728 2533308684 93103 d -1 /
      data xg ( 16) /   9.3281280827 8676533360 8521668452 05716 d -1 /
      data xg ( 17) /   9.5791681921 3791655804 5409994527 59285 d -1 /
      data xg ( 18) /   9.7725994998 3774262663 3702837129 03806 d -1 /
      data xg ( 19) /   9.9072623869 9457006453 0543522213 72154 d -1 /
      data xg ( 20) /   9.9823770971 0559200349 6227024205 86492 d -1 /

      data wg (  1) /   7.7505947978 4248112637 2396295832 63269 d -2 /
      data wg (  2) /   7.7039818164 2479655883 0753428381 02485 d -2 /
      data wg (  3) /   7.6110361900 6262423715 5807592249 48230 d -2 /
      data wg (  4) /   7.4723169057 9682642001 8933626132 46731 d -2 /
      data wg (  5) /   7.2886582395 8040590605 1068344251 78358 d -2 /
      data wg (  6) /   7.0611647391 2867796954 8363085528 68323 d -2 /
      data wg (  7) /   6.7912045815 2339038256 9010823192 39859 d -2 /
      data wg (  8) /   6.4804013456 6010380745 5452956675 27300 d -2 /
      data wg (  9) /   6.1306242492 9289391665 3799640839 85959 d -2 /
      data wg ( 10) /   5.7439769099 3915513666 1773091042 59856 d -2 /
      data wg ( 11) /   5.3227846983 9368243549 9647977226 05045 d -2 /
      data wg ( 12) /   4.8695807635 0722320614 3416044814 63880 d -2 /
      data wg ( 13) /   4.3870908185 6732719916 7468604171 54958 d -2 /
      data wg ( 14) /   3.8782167974 4720176399 7203129044 61622 d -2 /
      data wg ( 15) /   3.3460195282 5478473926 7818308641 08489 d -2 /
      data wg ( 16) /   2.7937006980 0234010984 8915750772 10773 d -2 /
      data wg ( 17) /   2.2245849194 1669572615 0432418420 85732 d -2 /
      data wg ( 18) /   1.6421058381 9078887128 6348488236 39272 d -2 /
      data wg ( 19) /   1.0498284531 1528136147 4217106727 96523 d -2 /
      data wg ( 20) /   4.5212770985 3319125847 1732878185 33272 d -3 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end





      subroutine dqleg080(f,a,b,result,dresult,ddresult,par,n)
c..
c..80 point gauss-legendre rule for the fermi-dirac function and
c..its derivatives with respect to eta and theta.
c..on input f is the name of the subroutine containing the integrand,
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to subroutine f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 80-point gauss-legendre rule,
c..dresult is the derivative with respect to eta, and ddresult is the
c..derivative with respect to theta.
c..
c..note: since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc1,absc2,center,hlfrun,wg(40),xg(40),
     2                 fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


c the abscissae and weights are given for the interval (-1,1).
c xg     - abscissae of the 80-point gauss-legendre rule
c          for half of the usual run (-1,1), i.e.
c          the positive nodes of the 80-point rule
c wg     - weights of the 80-point gauss rule.
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   1.9511383256 7939976543 5123410745 45479 d -2 /
      data xg (  2) /   5.8504437152 4206686289 9332188341 77944 d -2 /
      data xg (  3) /   9.7408398441 5845990632 7845010493 69020 d -2 /
      data xg (  4) /   1.3616402280 9143886559 2410780007 17067 d -1 /
      data xg (  5) /   1.7471229183 2646812559 3390480112 86195 d -1 /
      data xg (  6) /   2.1299450285 7666132572 3885386663 21823 d -1 /
      data xg (  7) /   2.5095235839 2272120493 1588160350 04797 d -1 /
      data xg (  8) /   2.8852805488 4511853109 1393014347 13898 d -1 /
      data xg (  9) /   3.2566437074 7701914619 1129436273 58695 d -1 /
      data xg ( 10) /   3.6230475349 9487315619 0432863589 63588 d -1 /
      data xg ( 11) /   3.9839340588 1969227024 3796425175 33757 d -1 /
      data xg ( 12) /   4.3387537083 1756093062 3867003631 81958 d -1 /
      data xg ( 13) /   4.6869661517 0544477036 0783649358 08657 d -1 /
      data xg ( 14) /   5.0280411188 8784987593 6727503675 68003 d -1 /
      data xg ( 15) /   5.3614592089 7131932019 8572531254 00904 d -1 /
      data xg ( 16) /   5.6867126812 2709784725 4857866248 27158 d -1 /
      data xg ( 17) /   6.0033062282 9751743154 7462991640 06848 d -1 /
      data xg ( 18) /   6.3107577304 6871966247 9283872893 36863 d -1 /
      data xg ( 19) /   6.6085989898 6119801735 9671228443 17234 d -1 /
      data xg ( 20) /   6.8963764434 2027600771 2076124389 35266 d -1 /
      data xg ( 21) /   7.1736518536 2099880254 0682582938 15278 d -1 /
      data xg ( 22) /   7.4400029758 3597272316 5405279309 13673 d -1 /
      data xg ( 23) /   7.6950242013 5041373865 6160687490 26083 d -1 /
      data xg ( 24) /   7.9383271750 4605449948 6393117384 54358 d -1 /
      data xg ( 25) /   8.1695413868 1463470371 1249940122 95707 d -1 /
      data xg ( 26) /   8.3883147358 0255275616 6230439028 67064 d -1 /
      data xg ( 27) /   8.5943140666 3111096977 1921234916 56492 d -1 /
      data xg ( 28) /   8.7872256767 8213828703 7733436391 24407 d -1 /
      data xg ( 29) /   8.9667557943 8770683194 3240719673 95986 d -1 /
      data xg ( 30) /   9.1326310257 1757654164 7336561509 47478 d -1 /
      data xg ( 31) /   9.2845987717 2445795953 0459590754 53133 d -1 /
      data xg ( 32) /   9.4224276130 9872674752 2660045000 01735 d -1 /
      data xg ( 33) /   9.5459076634 3634905493 4815170210 29508 d -1 /
      data xg ( 34) /   9.6548508904 3799251452 2731556714 54998 d -1 /
      data xg ( 35) /   9.7490914058 5727793385 6452300691 36276 d -1 /
      data xg ( 36) /   9.8284857273 8629070418 2880277091 16473 d -1 /
      data xg ( 37) /   9.8929130249 9755531026 5031671366 31385 d -1 /
      data xg ( 38) /   9.9422754096 5688277892 0635036649 11698 d -1 /
      data xg ( 39) /   9.9764986439 8237688899 4942081831 22985 d -1 /
      data xg ( 40) /   9.9955382265 1630629880 0804990945 67184 d -1 /

      data wg (  1) /   3.9017813656 3066548112 8043925275 40483 d -2 /
      data wg (  2) /   3.8958395962 7695311986 2552477226 08223 d -2 /
      data wg (  3) /   3.8839651059 0519689317 7418266878 71658 d -2 /
      data wg (  4) /   3.8661759774 0764633270 7711026715 66912 d -2 /
      data wg (  5) /   3.8424993006 9594231852 1243632949 01384 d -2 /
      data wg (  6) /   3.8129711314 4776383442 0679156573 62019 d -2 /
      data wg (  7) /   3.7776364362 0013974897 7497642632 10547 d -2 /
      data wg (  8) /   3.7365490238 7304900267 0537705783 86691 d -2 /
      data wg (  9) /   3.6897714638 2760088391 5099657340 52192 d -2 /
      data wg ( 10) /   3.6373749905 8359780439 6499104652 28136 d -2 /
      data wg ( 11) /   3.5794393953 4160546028 6158881615 44542 d -2 /
      data wg ( 12) /   3.5160529044 7475934955 2659238869 68812 d -2 /
      data wg ( 13) /   3.4473120451 7539287943 6422673102 98320 d -2 /
      data wg ( 14) /   3.3733214984 6115228166 7516306423 87284 d -2 /
      data wg ( 15) /   3.2941939397 6454013828 3618090195 95361 d -2 /
      data wg ( 16) /   3.2100498673 4877731480 5649028725 06960 d -2 /
      data wg ( 17) /   3.1210174188 1147016424 4286672060 35518 d -2 /
      data wg ( 18) /   3.0272321759 5579806612 2001009090 11747 d -2 /
      data wg ( 19) /   2.9288369583 2678476927 6758601957 91396 d -2 /
      data wg ( 20) /   2.8259816057 2768623967 5319796501 45302 d -2 /
      data wg ( 21) /   2.7188227500 4863806744 1870668054 42598 d -2 /
      data wg ( 22) /   2.6075235767 5651179029 6874360026 92871 d -2 /
      data wg ( 23) /   2.4922535764 1154911051 1784700321 98023 d -2 /
      data wg ( 24) /   2.3731882865 9301012931 9252461356 84162 d -2 /
      data wg ( 25) /   2.2505090246 3324619262 2158968616 87390 d -2 /
      data wg ( 26) /   2.1244026115 7820063887 1073725061 31285 d -2 /
      data wg ( 27) /   1.9950610878 1419989288 9192871511 35633 d -2 /
      data wg ( 28) /   1.8626814208 2990314287 3541415215 72090 d -2 /
      data wg ( 29) /   1.7274652056 2693063585 8420713129 09998 d -2 /
      data wg ( 30) /   1.5896183583 7256880449 0290922917 85257 d -2 /
      data wg ( 31) /   1.4493508040 5090761169 6207458346 05500 d -2 /
      data wg ( 32) /   1.3068761592 4013392937 8682589705 63403 d -2 /
      data wg ( 33) /   1.1624114120 7978269164 6676999543 26348 d -2 /
      data wg ( 34) /   1.0161766041 1030645208 3185035240 69436 d -2 /
      data wg ( 35) /   8.6839452692 6085842640 9452204034 28135 d -3 /
      data wg ( 36) /   7.1929047681 1731275267 5570867956 50747 d -3 /
      data wg ( 37) /   5.6909224514 0319864926 9107117162 01847 d -3 /
      data wg ( 38) /   4.1803131246 9489523673 9304201681 35132 d -3 /
      data wg ( 39) /   2.6635335895 1268166929 3535831668 45546 d -3 /
      data wg ( 40) /   1.1449500031 8694153454 4171941315 63611 d -3 /


c           list of major variables
c           -----------------------
c
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end





      subroutine dqlag010(f,a,b,result,dresult,ddresult,par,n)
c..
c..10 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=exp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 10-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(10),xg(10),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 10-point gauss-laguerre rule
c wg     - weights of the 10-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.3779347054 0492430830 7725056527 11188 d -1 /
      data xg (  2) /   7.2945454950 3170498160 3731216760 78781 d -1 /
      data xg (  3) /   1.8083429017 4031604823 2920075750 60883 d  0 /
      data xg (  4) /   3.4014336978 5489951448 2532221408 39067 d  0 /
      data xg (  5) /   5.5524961400 6380363241 7558486868 76285 d  0 /
      data xg (  6) /   8.3301527467 6449670023 8767197274 52218 d  0 /
      data xg (  7) /   1.1843785837 9000655649 1853891914 16139 d  1 /
      data xg (  8) /   1.6279257831 3781020995 3265393583 36223 d  1 /
      data xg (  9) /   2.1996585811 9807619512 7709019559 44939 d  1 /
      data xg ( 10) /   2.9920697012 2738915599 0879334079 91951 d  1 /

      data wg (  1) /   3.5400973860 6996308762 2268914420 67608 d -1 /
      data wg (  2) /   8.3190230104 3580738109 8296581278 49577 d -1 /
      data wg (  3) /   1.3302885617 4932817875 2792194393 99369 d  0 /
      data wg (  4) /   1.8630639031 1113098976 3988735482 46693 d  0 /
      data wg (  5) /   2.4502555580 8301016607 2693731657 52256 d  0 /
      data wg (  6) /   3.1227641551 3518249615 0818263314 55472 d  0 /
      data wg (  7) /   3.9341526955 6152109865 5812459248 23077 d  0 /
      data wg (  8) /   4.9924148721 9302310201 1485652433 15445 d  0 /
      data wg (  9) /   6.5722024851 3080297518 7668710376 11234 d  0 /
      data wg ( 10) /   9.7846958403 7463069477 0086638718 59813 d  0 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 10-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end





      subroutine dqlag020(f,a,b,result,dresult,ddresult,par,n)
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(20),xg(20),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.0539889691 9887533666 8900458421 50958 d -2 /
      data xg (  2) /   3.7212681800 1611443794 2413887611 46636 d -1 /
      data xg (  3) /   9.1658210248 3273564667 7162770741 83187 d -1 /
      data xg (  4) /   1.7073065310 2834388068 7689667413 05070 d  0 /
      data xg (  5) /   2.7491992553 0943212964 5030460494 81338 d  0 /
      data xg (  6) /   4.0489253138 5088692237 4953369133 33219 d  0 /
      data xg (  7) /   5.6151749708 6161651410 4539885651 89234 d  0 /
      data xg (  8) /   7.4590174536 7106330976 8860218371 81759 d  0 /
      data xg (  9) /   9.5943928695 8109677247 3672734282 79837 d  0 /
      data xg ( 10) /   1.2038802546 9643163096 2340929886 55158 d  1 /
      data xg ( 11) /   1.4814293442 6307399785 1267971004 79756 d  1 /
      data xg ( 12) /   1.7948895520 5193760173 6579099261 25096 d  1 /
      data xg ( 13) /   2.1478788240 2850109757 3517036959 46692 d  1 /
      data xg ( 14) /   2.5451702793 1869055035 1867748464 15418 d  1 /
      data xg ( 15) /   2.9932554631 7006120067 1365613516 58232 d  1 /
      data xg ( 16) /   3.5013434240 4790000062 8493590668 81395 d  1 /
      data xg ( 17) /   4.0833057056 7285710620 2956770780 75526 d  1 /
      data xg ( 18) /   4.7619994047 3465021399 4162715285 11211 d  1 /
      data xg ( 19) /   5.5810795750 0638988907 5077344449 72356 d  1 /
      data xg ( 20) /   6.6524416525 6157538186 4031879146 06659 d  1 /

      data wg (  1) /   1.8108006241 8989255451 6754059131 10644 d -1 /
      data wg (  2) /   4.2255676787 8563974520 3441725664 58197 d -1 /
      data wg (  3) /   6.6690954670 1848150373 4821149925 15927 d -1 /
      data wg (  4) /   9.1535237278 3073672670 6046847718 68067 d -1 /
      data wg (  5) /   1.1695397071 9554597380 1478222395 77476 d  0 /
      data wg (  6) /   1.4313549859 2820598636 8449948915 14331 d  0 /
      data wg (  7) /   1.7029811379 8502272402 5332616332 06720 d  0 /
      data wg (  8) /   1.9870158907 9274721410 9218392751 29020 d  0 /
      data wg (  9) /   2.2866357812 5343078546 2228546814 95651 d  0 /
      data wg ( 10) /   2.6058347275 5383333269 4989509540 33323 d  0 /
      data wg ( 11) /   2.9497837342 1395086600 2354168272 85951 d  0 /
      data wg ( 12) /   3.3253957820 0931955236 9519374217 51118 d  0 /
      data wg ( 13) /   3.7422554705 8981092111 7072932653 77811 d  0 /
      data wg ( 14) /   4.2142367102 5188041986 8080637824 78746 d  0 /
      data wg ( 15) /   4.7625184614 9020929695 2921978390 96371 d  0 /
      data wg ( 16) /   5.4217260442 4557430380 3082979899 81779 d  0 /
      data wg ( 17) /   6.2540123569 3242129289 5184903007 07542 d  0 /
      data wg ( 18) /   7.3873143890 5443455194 0300191964 64791 d  0 /
      data wg ( 19) /   9.1513287309 8747960794 3482425529 50528 d  0 /
      data wg ( 20) /   1.2893388645 9399966710 2628712874 85278 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end




      subroutine dqlag040(f,a,b,result,dresult,ddresult,par,n)
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(40),xg(40),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.5700394308 8883851220 8447128660 08554 d -2 /
      data xg (  2) /   1.8816228315 8698516003 5893462190 95913 d -1 /
      data xg (  3) /   4.6269428131 4576453564 9375245611 90364 d -1 /
      data xg (  4) /   8.5977296397 2934922257 2722246887 22412 d -1 /
      data xg (  5) /   1.3800108205 2733718649 8000329595 26559 d  0 /
      data xg (  6) /   2.0242091359 2282673344 2066002800 13075 d  0 /
      data xg (  7) /   2.7933693535 0681645765 3514486026 64039 d  0 /
      data xg (  8) /   3.6887026779 0827020959 1526351908 68698 d  0 /
      data xg (  9) /   4.7116411465 5497269361 8722836277 47369 d  0 /
      data xg ( 10) /   5.8638508783 4371811427 3164237995 82987 d  0 /
      data xg ( 11) /   7.1472479081 0228825068 5691951979 42362 d  0 /
      data xg ( 12) /   8.5640170175 8616376271 8522042088 13232 d  0 /
      data xg ( 13) /   1.0116634048 4519394068 4962965639 52448 d  1 /
      data xg ( 14) /   1.1807892294 0045848428 4158670436 06304 d  1 /
      data xg ( 15) /   1.3640933712 5370872283 7167636065 01202 d  1 /
      data xg ( 16) /   1.5619285893 3390738372 0196365218 80145 d  1 /
      data xg ( 17) /   1.7746905950 0956630425 7387749542 43772 d  1 /
      data xg ( 18) /   2.0028232834 5748905296 1261481017 51172 d  1 /
      data xg ( 19) /   2.2468249983 4984183513 7178622899 45366 d  1 /
      data xg ( 20) /   2.5072560772 4262037943 9608620940 09769 d  1 /
      data xg ( 21) /   2.7847480009 1688627207 5170414045 57997 d  1 /
      data xg ( 22) /   3.0800145739 4454627007 5438519619 11114 d  1 /
      data xg ( 23) /   3.3938657084 9137196090 9885858628 19990 d  1 /
      data xg ( 24) /   3.7272245880 4760043283 2076099060 74207 d  1 /
      data xg ( 25) /   4.0811492823 8869204661 5567558160 06426 d  1 /
      data xg ( 26) /   4.4568603175 3344627071 2302063449 83559 d  1 /
      data xg ( 27) /   4.8557763533 0599922809 6204880670 67936 d  1 /
      data xg ( 28) /   5.2795611187 2169329693 5202113739 17638 d  1 /
      data xg ( 29) /   5.7301863323 3936274950 3374699589 21651 d  1 /
      data xg ( 30) /   6.2100179072 7751116121 6819905789 89921 d  1 /
      data xg ( 31) /   6.7219370927 1269987990 8027755188 87054 d  1 /
      data xg ( 32) /   7.2695158847 6124621175 2192772426 19385 d  1 /
      data xg ( 33) /   7.8572802911 5713092805 4389683348 12596 d  1 /
      data xg ( 34) /   8.4911231135 7049845427 0156470966 63186 d  1 /
      data xg ( 35) /   9.1789874671 2363769923 3719348062 73153 d  1 /
      data xg ( 36) /   9.9320808717 4468082501 0905416548 68123 d  1 /
      data xg ( 37) /   1.0767244063 9388272520 7967676113 22664 d  2 /
      data xg ( 38) /   1.1712230951 2690688807 6506441235 50702 d  2 /
      data xg ( 39) /   1.2820184198 8255651192 5411043896 31263 d  2 /
      data xg ( 40) /   1.4228004446 9159997888 3488353595 41764 d  2 /

      data wg (  1) /   9.1625471157 4598973115 1169808013 74830 d -2 /
      data wg (  2) /   2.1342058490 5012080007 1933671215 12341 d -1 /
      data wg (  3) /   3.3571811668 0284673880 5107016162 92191 d -1 /
      data wg (  4) /   4.5854093503 3497560385 4323803764 52497 d -1 /
      data wg (  5) /   5.8206816577 9105168990 9963654015 43283 d -1 /
      data wg (  6) /   7.0649521636 7219392989 8300156730 16682 d -1 /
      data wg (  7) /   8.3202690300 3485238099 1129479783 49523 d -1 /
      data wg (  8) /   9.5887819879 4443111448 1226796760 28906 d -1 /
      data wg (  9) /   1.0872761620 3054971575 3869333172 02661 d  0 /
      data wg ( 10) /   1.2174623279 7778097895 4277850665 60948 d  0 /
      data wg ( 11) /   1.3496954913 5676530792 3938594423 94519 d  0 /
      data wg ( 12) /   1.4842549297 7684671120 5611786129 78719 d  0 /
      data wg ( 13) /   1.6214441628 1182197802 3168843164 54527 d  0 /
      data wg ( 14) /   1.7615953746 7676961118 4242204209 81598 d  0 /
      data wg ( 15) /   1.9050746658 9479967668 2993205972 79371 d  0 /
      data wg ( 16) /   2.0522883472 6171671760 1995822729 47454 d  0 /
      data wg ( 17) /   2.2036905532 4509588909 8283443281 40570 d  0 /
      data wg ( 18) /   2.3597925385 2320332354 0373753789 01497 d  0 /
      data wg ( 19) /   2.5211741403 7643299165 3136902874 22820 d  0 /
      data wg ( 20) /   2.6884980554 0884226415 9505447063 74659 d  0 /
      data wg ( 21) /   2.8625278132 1044881203 4763959831 04311 d  0 /
      data wg ( 22) /   3.0441506653 1151710041 0439679543 33670 d  0 /
      data wg ( 23) /   3.2344070972 6353194177 4902394288 67111 d  0 /
      data wg ( 24) /   3.4345293984 2774809220 3984818916 02464 d  0 /
      data wg ( 25) /   3.6459928249 9408907238 9656466994 90434 d  0 /
      data wg ( 26) /   3.8705845972 1651656808 4753202134 44338 d  0 /
      data wg ( 27) /   4.1104986804 3282265583 5822472639 51577 d  0 /
      data wg ( 28) /   4.3684687232 5406347450 8083382729 45025 d  0 /
      data wg ( 29) /   4.6479589840 7446688299 3033998838 83991 d  0 /
      data wg ( 30) /   4.9534461124 0989326218 6961507855 62721 d  0 /
      data wg ( 31) /   5.2908484059 0073657468 7373657188 58968 d  0 /
      data wg ( 32) /   5.6682046090 3297677000 7305290232 63795 d  0 /
      data wg ( 33) /   6.0967964147 4342030593 3760108591 98806 d  0 /
      data wg ( 34) /   6.5931088610 3999953794 4296642062 94899 d  0 /
      data wg ( 35) /   7.1824959955 3689315064 4298016266 99574 d  0 /
      data wg ( 36) /   7.9066663113 8422877369 3107423105 86595 d  0 /
      data wg ( 37) /   8.8408924928 1034652079 1255950630 26792 d  0 /
      data wg ( 38) /   1.0140899265 6211694839 0946003069 40468 d  1 /
      data wg ( 39) /   1.2210021299 2046038985 2264858758 81108 d  1 /
      data wg ( 40) /   1.6705520642 0242974052 4687743985 73553 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula


      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end





      subroutine dqlag080(f,a,b,result,dresult,ddresult,par,n)
c..
c..20 point gauss-laguerre rule for the fermi-dirac function.
c..on input f is the external function defining the integrand
c..f(x)=g(x)*w(x), where w(x) is the gaussian weight 
c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function, 
c..a is the lower end point of the interval, b is the higher end point,
c..par is an array of constant parameters to be passed to the function f,
c..and n is the length of the par array. on output result is the 
c..approximation from applying the 20-point gauss-laguerre rule. 
c..since the number of nodes is even, zero is not an abscissa.
c..
c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),
     1                 absc,wg(80),xg(80),fval,dfval,ddfval


c the abscissae and weights are given for the interval (0,+inf).
c xg     - abscissae of the 20-point gauss-laguerre rule
c wg     - weights of the 20-point gauss rule. since f yet
c          includes the weight function, the values in wg 
c          are actually exp(xg) times the standard
c          gauss-laguerre weights
c
c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.7960423300 6983655540 1031924740 16803 d -2 /
      data xg (  2) /   9.4639912994 3539888113 9027246521 72943 d -2 /
      data xg (  3) /   2.3262286812 5867569207 7061572163 49831 d -1 /
      data xg (  4) /   4.3199254780 2387480255 7861724977 70411 d -1 /
      data xg (  5) /   6.9282886135 2021839905 7022136354 46867 d -1 /
      data xg (  6) /   1.0152325561 8947143744 6254368599 35350 d  0 /
      data xg (  7) /   1.3993276878 4287277414 4190514309 78382 d  0 /
      data xg (  8) /   1.8452623038 3584513811 1771177695 99966 d  0 /
      data xg (  9) /   2.3532088716 0926152447 2447080161 40181 d  0 /
      data xg ( 10) /   2.9233646865 5542632483 6912342597 32862 d  0 /
      data xg ( 11) /   3.5559523140 4613405944 9673083246 38370 d  0 /
      data xg ( 12) /   4.2512200823 0987808316 4857664485 77637 d  0 /
      data xg ( 13) /   5.0094426336 2016477243 3677068182 06389 d  0 /
      data xg ( 14) /   5.8309215386 0871901982 1271132956 05083 d  0 /
      data xg ( 15) /   6.7159859778 5131711156 5500876351 99430 d  0 /
      data xg ( 16) /   7.6649934948 9177306073 4189090478 23480 d  0 /
      data xg ( 17) /   8.6783308251 6770109543 4422555426 61083 d  0 /
      data xg ( 18) /   9.7564148057 4293071316 5509446173 66591 d  0 /
      data xg ( 19) /   1.0899693371 2878553774 3610010214 89406 d  1 /
      data xg ( 20) /   1.2108646642 3656999007 0548486983 15593 d  1 /
      data xg ( 21) /   1.3383788112 7786473701 6298406038 33297 d  1 /
      data xg ( 22) /   1.4725665943 5085855393 3580762618 38437 d  1 /
      data xg ( 23) /   1.6134864371 6624665791 6585454289 90907 d  1 /
      data xg ( 24) /   1.7612005243 8144378598 6356869435 86520 d  1 /
      data xg ( 25) /   1.9157749684 2412479221 7299702056 74985 d  1 /
      data xg ( 26) /   2.0772799909 7920960924 4193790104 89579 d  1 /
      data xg ( 27) /   2.2457901204 5404583114 0959169508 77516 d  1 /
      data xg ( 28) /   2.4213844068 9586473771 9224693924 47092 d  1 /
      data xg ( 29) /   2.6041466560 1655866929 3900535654 35682 d  1 /
      data xg ( 30) /   2.7941656841 8594655558 2330692936 92111 d  1 /
      data xg ( 31) /   2.9915355964 9009855011 2704121157 37715 d  1 /
      data xg ( 32) /   3.1963560902 2089207107 7488875426 36533 d  1 /
      data xg ( 33) /   3.4087327864 7261898749 8349473428 60505 d  1 /
      data xg ( 34) /   3.6287775928 7814544588 0319884362 16948 d  1 /
      data xg ( 35) /   3.8566091009 2922104582 5630521729 08535 d  1 /
      data xg ( 36) /   4.0923530218 0312671999 0958505955 44326 d  1 /
      data xg ( 37) /   4.3361426651 7312302957 8267604682 19500 d  1 /
      data xg ( 38) /   4.5881194661 2788863456 2664899748 78378 d  1 /
      data xg ( 39) /   4.8484335660 8331891358 7372733535 63006 d  1 /
      data xg ( 40) /   5.1172444544 6070105959 8894323349 07144 d  1 /
      data xg ( 41) /   5.3947216789 5544471206 2102787875 72430 d  1 /
      data xg ( 42) /   5.6810456334 6362231341 2485032441 02122 d  1 /
      data xg ( 43) /   5.9764084342 1099549427 2959612774 71927 d  1 /
      data xg ( 44) /   6.2810148963 9264772036 2729175902 88682 d  1 /
      data xg ( 45) /   6.5950836257 4560573434 6406271607 92248 d  1 /
      data xg ( 46) /   6.9188482420 2362773741 9802886482 37373 d  1 /
      data xg ( 47) /   7.2525587544 2633453588 3896526165 68450 d  1 /
      data xg ( 48) /   7.5964831127 8641748269 4497949747 96502 d  1 /
      data xg ( 49) /   7.9509089629 0888369620 5728262599 80809 d  1 /
      data xg ( 50) /   8.3161456401 0536896630 4295068758 48705 d  1 /
      data xg ( 51) /   8.6925264419 6156234481 1659260404 48396 d  1 /
      data xg ( 52) /   9.0804112300 9407559518 4117278203 18427 d  1 /
      data xg ( 53) /   9.4801894215 9474332072 0718891387 35302 d  1 /
      data xg ( 54) /   9.8922834446 9405791648 0193727380 36790 d  1 /
      data xg ( 55) /   1.0317152750 8039130233 0470941673 45654 d  2 /
      data xg ( 56) /   1.0755298497 7539906327 6078907989 75954 d  2 /
      data xg ( 57) /   1.1207269048 4128333623 9300461662 11013 d  2 /
      data xg ( 58) /   1.1673666467 3503666318 1578881308 01099 d  2 /
      data xg ( 59) /   1.2155154249 0952625566 8638957521 10813 d  2 /
      data xg ( 60) /   1.2652466579 6515540341 5702654316 53573 d  2 /
      data xg ( 61) /   1.3166419525 2120310870 0890863080 06192 d  2 /
      data xg ( 62) /   1.3697924668 6936973947 5706372894 63788 d  2 /
      data xg ( 63) /   1.4248005891 2161601930 8265692004 55232 d  2 /
      data xg ( 64) /   1.4817820245 5004441818 6523848360 07732 d  2 /
      data xg ( 65) /   1.5408684228 1798697859 4174252655 96259 d  2 /
      data xg ( 66) /   1.6022107287 0095715935 2684168930 10646 d  2 /
      data xg ( 67) /   1.6659835193 4053918744 5211797337 12213 d  2 /
      data xg ( 68) /   1.7323907133 4249503830 9065037750 56999 d  2 /
      data xg ( 69) /   1.8016732304 9032317982 4302089977 01523 d  2 /
      data xg ( 70) /   1.8741194967 6963772390 4901345880 21771 d  2 /
      data xg ( 71) /   1.9500802244 1532991450 3904796005 99643 d  2 /
      data xg ( 72) /   2.0299898419 5074937824 8076778237 14777 d  2 /
      data xg ( 73) /   2.1143987049 4836466691 4849046955 42608 d  2 /
      data xg ( 74) /   2.2040236815 1735739654 0442066777 63168 d  2 /
      data xg ( 75) /   2.2998320607 5680004348 4109696758 44754 d  2 /
      data xg ( 76) /   2.4031908705 5841540417 5974604792 19628 d  2 /
      data xg ( 77) /   2.5161587933 0499611167 4449393109 73194 d  2 /
      data xg ( 78) /   2.6421382388 3199102097 6961086914 35553 d  2 /
      data xg ( 79) /   2.7876673304 6004563652 0141725306 11597 d  2 /
      data xg ( 80) /   2.9696651199 5651345758 8528591557 03581 d  2 /

      data wg (  1) /   4.6093103133 0609664705 2513213955 10083 d -2 /
      data wg (  2) /   1.0731300778 3932752564 1503203043 98860 d -1 /
      data wg (  3) /   1.6866442954 7948111794 2204577827 02406 d -1 /
      data wg (  4) /   2.3008808938 4940054411 2571819781 93282 d -1 /
      data wg (  5) /   2.9160130250 2437964832 1693187729 43752 d -1 /
      data wg (  6) /   3.5322675357 5408236352 7231258056 47046 d -1 /
      data wg (  7) /   4.1498817755 0940466187 1976863112 80092 d -1 /
      data wg (  8) /   4.7690979230 2936241314 7770254185 05661 d -1 /
      data wg (  9) /   5.3901621847 4955374499 5076565223 27912 d -1 /
      data wg ( 10) /   6.0133249744 7190529086 7652488407 39512 d -1 /
      data wg ( 11) /   6.6388413639 6680571849 4422407272 99214 d -1 /
      data wg ( 12) /   7.2669716361 4156688973 5672962491 40514 d -1 /
      data wg ( 13) /   7.8979818942 8428531349 7930783987 88294 d -1 /
      data wg ( 14) /   8.5321447143 8152298354 5981624313 62968 d -1 /
      data wg ( 15) /   9.1697398383 3892698590 3429000315 53302 d -1 /
      data wg ( 16) /   9.8110549100 4005747195 0601559842 18607 d -1 /
      data wg ( 17) /   1.0456386258 0654218147 5684456631 76029 d  0 /
      data wg ( 18) /   1.1106039730 0025890771 1247632597 29371 d  0 /
      data wg ( 19) /   1.1760331584 1226175056 6510765192 08666 d  0 /
      data wg ( 20) /   1.2419589444 9809359279 3517618178 71338 d  0 /
      data wg ( 21) /   1.3084153330 3134064261 1885428459 54645 d  0 /
      data wg ( 22) /   1.3754376757 4892843813 1559170934 90796 d  0 /
      data wg ( 23) /   1.4430627938 7849270398 3124172072 47308 d  0 /
      data wg ( 24) /   1.5113291075 8830693847 6550205599 17703 d  0 /
      data wg ( 25) /   1.5802767765 3099415830 2018787231 21659 d  0 /
      data wg ( 26) /   1.6499478528 0267874116 0120428193 55036 d  0 /
      data wg ( 27) /   1.7203864478 1283277182 0042814522 90770 d  0 /
      data wg ( 28) /   1.7916389147 6093832891 4426205276 88915 d  0 /
      data wg ( 29) /   1.8637540486 4909708435 9257090286 88162 d  0 /
      data wg ( 30) /   1.9367833060 3070923513 9254343278 41646 d  0 /
      data wg ( 31) /   2.0107810470 1134222912 6149881755 55546 d  0 /
      data wg ( 32) /   2.0858048023 8741046429 3039785129 89079 d  0 /
      data wg ( 33) /   2.1619155692 4159897378 3163440488 27763 d  0 /
      data wg ( 34) /   2.2391781388 2364652373 4539974474 45645 d  0 /
      data wg ( 35) /   2.3176614611 4651854068 6060480434 96370 d  0 /
      data wg ( 36) /   2.3974390514 4001430514 1172386388 49980 d  0 /
      data wg ( 37) /   2.4785894444 4973417756 3691644552 22527 d  0 /
      data wg ( 38) /   2.5611967035 7790455335 1155092225 72643 d  0 /
      data wg ( 39) /   2.6453509930 6968892850 4634410003 67534 d  0 /
      data wg ( 40) /   2.7311492228 9915138861 4102871311 69260 d  0 /
      data wg ( 41) /   2.8186957777 5934171703 1418737478 11157 d  0 /
      data wg ( 42) /   2.9081033436 8223018934 5502767774 92687 d  0 /
      data wg ( 43) /   2.9994938483 9685626832 4124518299 68724 d  0 /
      data wg ( 44) /   3.0929995346 9357468116 6951083530 33660 d  0 /
      data wg ( 45) /   3.1887641899 4712376429 3652715016 23466 d  0 /
      data wg ( 46) /   3.2869445597 5337531998 3781070122 16956 d  0 /
      data wg ( 47) /   3.3877119796 0397652334 0549087621 54571 d  0 /
      data wg ( 48) /   3.4912542659 8732012281 7324237827 64895 d  0 /
      data wg ( 49) /   3.5977779176 9613046096 2947301749 02943 d  0 /
      data wg ( 50) /   3.7075106900 1745708341 0271556592 28179 d  0 /
      data wg ( 51) /   3.8207046196 5311695152 0299594304 67622 d  0 /
      data wg ( 52) /   3.9376395977 1430720676 8005406573 30923 d  0 /
      data wg ( 53) /   4.0586276133 8354481597 4201161879 88679 d  0 /
      data wg ( 54) /   4.1840178238 1424031850 6076923345 03121 d  0 /
      data wg ( 55) /   4.3142026492 9613425820 0845732179 87912 d  0 /
      data wg ( 56) /   4.4496251505 3655906604 9828201553 77774 d  0 /
      data wg ( 57) /   4.5907880226 3617511042 9598491489 29810 d  0 /
      data wg ( 58) /   4.7382646459 8929537394 7538735058 38770 d  0 /
      data wg ( 59) /   4.8927127796 6692168696 8869367432 83567 d  0 /
      data wg ( 60) /   5.0548916853 4039512820 5725071351 75938 d  0 /
      data wg ( 61) /   5.2256837559 4272391089 2780101660 22467 d  0 /
      data wg ( 62) /   5.4061221337 9727909853 3235123407 17863 d  0 /
      data wg ( 63) /   5.5974264018 4041404016 5536941589 80053 d  0 /
      data wg ( 64) /   5.8010493213 7643943530 6261624553 94841 d  0 /
      data wg ( 65) /   6.0187389387 8099108768 0151515140 26344 d  0 /
      data wg ( 66) /   6.2526224749 1437403092 9342134800 91928 d  0 /
      data wg ( 67) /   6.5053217351 7668675787 4827196636 96133 d  0 /
      data wg ( 68) /   6.7801152120 0777294201 2873479800 59368 d  0 /
      data wg ( 69) /   7.0811712202 5414518776 1743119167 59402 d  0 /
      data wg ( 70) /   7.4138924461 5305421421 6956062266 87752 d  0 /
      data wg ( 71) /   7.7854415484 1612700386 2327403392 30532 d  0 /
      data wg ( 72) /   8.2055734781 4596472333 9050861009 17119 d  0 /
      data wg ( 73) /   8.6880138399 6161871469 4199580582 55237 d  0 /
      data wg ( 74) /   9.2528697341 5578523923 5565062019 79918 d  0 /
      data wg ( 75) /   9.9311447184 0215736008 3709865340 09772 d  0 /
      data wg ( 76) /   1.0773973641 4646829405 7508435229 90655 d  1 /
      data wg ( 77) /   1.1873891246 5097447081 9508877108 77400 d  1 /
      data wg ( 78) /   1.3422885849 7264236139 7349401540 89734 d  1 /
      data wg ( 79) /   1.5919780161 6897924449 5542522001 85978 d  1 /
      data wg ( 80) /   2.1421454296 4372259537 5210361864 15127 d  1 /


c           list of major variables
c           -----------------------
c           absc   - abscissa
c           fval*  - function value
c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,80
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end
