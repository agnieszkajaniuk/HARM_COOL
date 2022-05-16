module mnewtf

implicit none
public:: &
     mnewt

contains 

SUBROUTINE lubksb(a,n,np,indx,b) 
DOUBLE PRECISION a(np,np),b(n) 
INTEGER i,ii,j,n,np,indx(n) 
DOUBLE PRECISION sum
INTEGER ll

ii=0 
do i=1,n 
    ll=indx(i) 
    sum=b(ll) 
    b(ll)=b(i) 
    if (ii/=0) then 
        do j=ii,i-1 
            sum=sum-a(i,j)*b(j) 
        end do 
    else if (sum/=0.) then 
        ii=i 
    endif 
    b(i)=sum 
end do 
do i=n,1,-1 
    sum=b(i) 
    if(i<n) then 
        do j=i+1,n 
            sum=sum-a(i,j)*b(j) 
        end do 
    endif 
    b(i)=sum/a(i,i) 
end do 
END SUBROUTINE lubksb 

INTEGER FUNCTION ludcmp(a,n,np,indx,d) 
INTEGER nmax,n,np
DOUBLE PRECISION tiny
PARAMETER (nmax=100,tiny=1.0d-20) 
DOUBLE PRECISION a(np,np),vv(nmax)
INTEGER indx(n)
INTEGER i,j,k,imax
DOUBLE PRECISION d,sum,aamax,dum


ludcmp=1
 
d=1. 
do i=1,n 
    aamax=0. 
    do j=1,n 
        if (dabs(a(i,j))>aamax) aamax=dabs(a(i,j)) 
    end do 
    if (aamax==0.) then 
		print*, 'singular matrix.'
		ludcmp=0
		return
	endif
    vv(i)=1./aamax 
end do 
do j=1,n 
    if (j>1) then 
        do i=1,j-1 
            sum=a(i,j) 
            if (i>1) then 
                do k=1,i-1 
                    sum=sum-a(i,k)*a(k,j) 
                end do 
                a(i,j)=sum 
            endif 
        end do 
    endif 
    aamax=0. 
    do i=j,n 
        sum=a(i,j) 
        if (j>1) then 
            do k=1,j-1 
                sum=sum-a(i,k)*a(k,j) 
            end do 
            a(i,j)=sum 
        endif 
        dum=vv(i)*dabs(sum) 
        if (dum>=aamax) then 
            imax=i 
            aamax=dum 
        endif 
    end do 
    if (j/=imax) then 
        do k=1,n 
            dum=a(imax,k) 
            a(imax,k)=a(j,k) 
            a(j,k)=dum 
        end do 
        d=-d 
        vv(imax)=vv(j) 
    endif 
    indx(j)=imax 
    if (j/=n) then 
        if (a(j,j)==0.) a(j,j)=tiny 
        dum=1./a(j,j) 
        do i=j+1,n 
            a(i,j)=a(i,j)*dum 
        end do 
    endif 
end do 
if(a(n,n)==0.) a(n,n)=tiny 
END FUNCTION ludcmp 



SUBROUTINE mnewt(ntrial,x,n,tolx,tolf,usrfun,par,npar,errf) 
INTEGER n,ntrial,NP,npar
DOUBLE PRECISION tolf,tolx,x(n),par(npar) 
PARAMETER (NP=15) 
!USES lubksb,ludcmp,usrfun 
INTEGER i,k,indx(n) 
DOUBLE PRECISION d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)
interface
         subroutine usrfun(x,n,np,fjac,fvec,par,npar) 
         integer n,np,npar
         double precision x(n),par(npar) 
         double precision fjac(np,np),fvec(np) 
         end subroutine usrfun
end interface 


do k=1,ntrial 

  call usrfun(x,n,np,fjac,fvec,par,npar) 

  errf=0. 
  do i=1,n 
    errf=errf+dabs(fvec(i)) 
  end do 
  if(errf<=tolf) return 
  if(ludcmp(fjac,n,np,indx,d)==1) then 
  	call lubksb(fjac,n,np,indx,fvec)
  endif
  errx=0. 
  do i=1,n 
    errx=errx+dabs(fvec(i)) 
    x(i)=x(i)+fvec(i) 
  end do 
  if(errx<=tolx) return 
end do 
END SUBROUTINE mnewt

end module mnewtf
