
module bisection

implicit none
public:: &
     rtbis

contains 
double precision function rtbis(FUNC,X1,X2,XACC)
        
  integer JMAX,J
  PARAMETER (JMAX=40)
  double precision X1,X2,XACC,FMID,F,DX,XMID,FUNC
  !double precision rtbis

  FMID=FUNC(X2)
  F=FUNC(X1)
  !         IF(F*FMID.GE.0.) PAUSE 
  !		 'Root must be bracketed for bisection.'

  IF(F.LT.0.)THEN
     rtbis=X1
     DX=X2-X1
  ELSE
     rtbis=X2
     DX=X1-X2
  ENDIF

  DO J=1,JMAX
     DX=DX*.5
     XMID=rtbis+DX
     FMID=FUNC(XMID)
     IF(FMID.LT.0.)rtbis=XMID
     IF(DABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
  ENDDO

     !        PAUSE 
     !		'too many bisections in rtbis'

END FUNCTION rtbis

end module bisection
