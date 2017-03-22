* $Id$
C
C
C
	FUNCTION ROMG(X)
C		********This gives the eikonal function from a table
C			calculated in the first call
	DIMENSION FR(0:1000)
	DATA I0/0/
	SAVE
	IF(I0.NE.0) GO TO 100
	DO 50 I=1,1001
	XR=(I-1)*0.01
	FR(I-1)=OMG0(XR)
50	CONTINUE
100	I0=1
	IF(X.GE.10.0) THEN
		ROMG=0.0
		RETURN
	ENDIF
	IX=INT(X*100)
	ROMG=(FR(IX)*((IX+1)*0.01-X)+FR(IX+1)*(X-IX*0.01))/0.01
	RETURN
	END
