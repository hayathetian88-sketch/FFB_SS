      SUBROUTINE FNDOSH(N,NEG,NE,LIST,TINY,
     *                 NP,NEX,X,Y,Z,NODE,IFRM,LEFRM,
     *                 XP,YP,ZP,XMINE,YMINE,ZMINE,XMAXE,YMAXE,ZMAXE,
     *                 EN,GP,EP,TP,ERR,IUT0,IERR)
      IMPLICIT NONE
C
***** DEFINE ARGUMENTS *****
      INTEGER N,NEG,NE,NP,NEX,NODE,IFRM,LEFRM(NEG)
      INTEGER LIST(NE) 
      REAL X,Y,Z
      REAL XP,YP,ZP,
     *     XMINE(NEG),YMINE(NEG),ZMINE(NEG),
     *     XMAXE(NEG),YMAXE(NEG),ZMAXE(NEG)
C
      INTEGER EN
      REAL GP,EP,TP,ERR
C
      DIMENSION NODE(N,NEG),NEX(8),X(NP),Y(NP),Z(NP)
C      
***** OBJECTS *****
      INTEGER JHEX,JWED,JPRD,JTET
      INTEGER IFOUND
      REAL XI(3),DXI(3),XPV(3),XQV(3),DX(3),X1,X2,X3
      REAL XE(8,3),NN(8),PSI(8,3),J(3,3),INVJ(3,3)
      INTEGER IUT0,IERR
      INTEGER IN,IE,IE0,IP,IT
      INTEGER ITMAX
      INTEGER II,JJ,KK

      REAL*4 TINY,DD,DDX,DDY,DDZ
      REAL*4 TMP,DETJ
      REAL*4 ERR1,ERR2,ERR3
      REAL*4 XBUF,YBUF,ZBUF
C
      INTEGER NTMP,NHEX,NPRD,NWED,NTET     
C
***** DATA *****
      DATA IFOUND  /0/
      DATA ITMAX  /10/
C
      DATA DD   /1.0E-1/ 
C      
      IFOUND=0
      IERR =0
C
      NTET=NEX(5)
      NPRD=NEX(6)
      NWED=NEX(7)
      NHEX=NEX(8)
C
***** MAKE XPV  *****
      XBUF=XP
      YBUF=YP
      ZBUF=ZP
C      
************************
***** ELEMENT LOOP *****
************************
      DO 1000 IE0=1,NE
C
          IE=LIST(IE0)
C
          XPV(1) = XBUF
          XPV(2) = YBUF
          XPV(3) = ZBUF
C
          IF(LEFRM(IE).NE.IFRM) GOTO 1000 
C
          DDX=XMAXE(IE)-XMINE(IE)
          DDY=YMAXE(IE)-YMINE(IE)
          DDZ=ZMAXE(IE)-ZMINE(IE)
          IF(     (XPV(1).LT.XMINE(IE)-DDX) 
     *       .OR. (XPV(1).GT.XMAXE(IE)+DDX)
     *       .OR. (XPV(2).LT.YMINE(IE)-DDY) 
     *       .OR. (XPV(2).GT.YMAXE(IE)+DDY)
     *       .OR. (XPV(3).LT.ZMINE(IE)-DDZ) 
     *       .OR. (XPV(3).GT.ZMAXE(IE)+DDZ) )
     *    GOTO 1000 
C
***** DETERMIN ELEMENT TYPE *****
          JHEX=1
          XI(1)=0.0
          XI(2)=0.0
          XI(3)=0.0
C
          DO 1100 II=1,8
              XE(II,1) =X(NODE(II,IE))
              XE(II,2) =Y(NODE(II,IE))
              XE(II,3) =Z(NODE(II,IE))
 1100     CONTINUE
C            
************HEX************
          DO 1200 IT=1,ITMAX    
                   NN(1)=(1.-XI(1))*(1.-XI(2))*(1.-XI(3))
                   NN(2)=(1.+XI(1))*(1.-XI(2))*(1.-XI(3))
                   NN(3)=(1.+XI(1))*(1.+XI(2))*(1.-XI(3))
                   NN(4)=(1.-XI(1))*(1.+XI(2))*(1.-XI(3))
                   NN(5)=(1.-XI(1))*(1.-XI(2))*(1.+XI(3))
                   NN(6)=(1.+XI(1))*(1.-XI(2))*(1.+XI(3))
                   NN(7)=(1.+XI(1))*(1.+XI(2))*(1.+XI(3))
                   NN(8)=(1.-XI(1))*(1.+XI(2))*(1.+XI(3))
C
              DO 1300 II=1,3
                  TMP=0.
                   DO 1310 JJ=1,8
                       TMP =TMP +XE(JJ,II)*NN(JJ)*0.125E0
 1310              CONTINUE
                   XQV(II)=TMP
                   DX (II)=XPV(II)-TMP
 1300         CONTINUE     
C
                  PSI(1,1)=-(1.-XI(2))*(1.-XI(3))
                  PSI(2,1)= (1.-XI(2))*(1.-XI(3))
                  PSI(3,1)= (1.+XI(2))*(1.-XI(3))
                  PSI(4,1)=-(1.+XI(2))*(1.-XI(3))
                  PSI(5,1)=-(1.-XI(2))*(1.+XI(3))
                  PSI(6,1)= (1.-XI(2))*(1.+XI(3))
                  PSI(7,1)= (1.+XI(2))*(1.+XI(3))
                  PSI(8,1)=-(1.+XI(2))*(1.+XI(3))
C
                  PSI(1,2)=-(1.-XI(3))*(1.-XI(1))
                  PSI(2,2)=-(1.-XI(3))*(1.+XI(1))
                  PSI(3,2)= (1.-XI(3))*(1.+XI(1))
                  PSI(4,2)= (1.-XI(3))*(1.-XI(1))
                  PSI(5,2)=-(1.+XI(3))*(1.-XI(1))
                  PSI(6,2)=-(1.+XI(3))*(1.+XI(1))
                  PSI(7,2)= (1.+XI(3))*(1.+XI(1))
                  PSI(8,2)= (1.+XI(3))*(1.-XI(1))
C
                  PSI(1,3)=-(1.-XI(1))*(1.-XI(2))
                  PSI(2,3)=-(1.+XI(1))*(1.-XI(2))
                  PSI(3,3)=-(1.+XI(1))*(1.+XI(2))
                  PSI(4,3)=-(1.-XI(1))*(1.+XI(2))
                  PSI(5,3)= (1.-XI(1))*(1.-XI(2))
                  PSI(6,3)= (1.+XI(1))*(1.-XI(2))
                  PSI(7,3)= (1.+XI(1))*(1.+XI(2))
                  PSI(8,3)= (1.-XI(1))*(1.+XI(2))
C
              DO 1400 II=1,3
                  DO 1410 JJ=1,3
                      TMP=0.
                      DO 1420 KK=1,8
                          TMP =TMP +PSI(KK,II)*XE(KK,JJ)*0.125E0
 1420                 CONTINUE
                      J(II,JJ)=TMP
 1410             CONTINUE
 1400         CONTINUE
C
              DETJ =J(1,1)*J(2,2)*J(3,3) -J(1,1)*J(2,3)*J(3,2)
     *             +J(1,2)*J(2,3)*J(3,1) -J(1,2)*J(2,1)*J(3,3)
     *             +J(1,3)*J(2,1)*J(3,2) -J(1,3)*J(2,2)*J(3,1)
              INVJ(1,1) =(J(2,2)*J(3,3)-J(2,3)*J(3,2))/DETJ
              INVJ(1,2) =(J(1,3)*J(3,2)-J(1,2)*J(3,3))/DETJ
              INVJ(1,3) =(J(1,2)*J(2,3)-J(1,3)*J(2,2))/DETJ
              INVJ(2,1) =(J(2,3)*J(3,1)-J(2,1)*J(3,3))/DETJ
              INVJ(2,2) =(J(1,1)*J(3,3)-J(1,3)*J(3,1))/DETJ
              INVJ(2,3) =(J(1,3)*J(2,1)-J(1,1)*J(2,3))/DETJ
              INVJ(3,1) =(J(2,1)*J(3,2)-J(2,2)*J(3,1))/DETJ
              INVJ(3,2) =(J(1,2)*J(3,1)-J(1,1)*J(3,2))/DETJ
              INVJ(3,3) =(J(1,1)*J(2,2)-J(1,2)*J(2,1))/DETJ
C
              DO 1500 II=1,3
                  DXI(II)= INVJ(1,II)*DX(1)
     *                    +INVJ(2,II)*DX(2)
     *                    +INVJ(3,II)*DX(3)
                   XI(II)=XI(II)+DXI(II)
 1500         CONTINUE     
C
              IF(    (XI(1).GE.-1-TINY) 
     *         .AND. (XI(1).LE.+1+TINY)
     *         .AND. (XI(2).GE.-1-TINY) 
     *         .AND. (XI(2).LE.+1+TINY)
     *         .AND. (XI(3).GE.-1-TINY) 
     *         .AND. (XI(3).LE.+1+TINY) ) THEN
                  EN =IE
                  GP =XI(1)
                  EP =XI(2)
                  TP =XI(3)
                  ERR1=ABS(GP)-1.E0
                  ERR2=ABS(EP)-1.E0
                  ERR3=ABS(TP)-1.E0
                  ERR =AMAX1(ERR1,ERR2,ERR3,0.E0)
C
              RETURN
          ENDIF
C            
 1200     CONTINUE
 1000 CONTINUE
C      
*********************
***** NOT FOUND *****
*********************
      EN=0
      ERR=1.0
      RETURN
C      
      RETURN
      END
