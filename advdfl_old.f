      SUBROUTINE ADVDFL_OLD
     * (MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *  N1,N2,NE,NP,MELM,NEX,NODE,
     *  UG,VG,WG,DT,SCT,VISC,
     *  CGAS,CLQD,P,SIGMA,SN,
     *  CEB,CUPWD,CVISL,ICAVI,
     *  E,EX,EY,EZ,EXX,EYY,EZZ,EXY,EYZ,EXZ,
     *  A,X,RHS,AD,
     *  DELTA,IERR)
      IMPLICIT NONE
C
      INTEGER*4 MCOLOR,MCPART
      INTEGER*4 ICOLOR,ICPART
      INTEGER*4 NCOLOR(4),NCPART(MCOLOR,4),LLOOP(MCPART,MCOLOR,4)
C[INTPUT]
      INTEGER*4 N1,N2,NE,NP,MELM,NEX(12)
      INTEGER*4 NODE(N2,NE)
      REAL*4    UG(NE),VG(NE),WG(NE),DT
      REAL*4    VISC(NE),SCT,CEB,X(NP),DELTA(NE)
      REAL(4),INTENT(IN) :: CGAS,CLQD,P(NE),SIGMA,SN(N1,NE)
      REAL*4    E(MELM),EX(MELM),EY(MELM),EZ(MELM),
     *          EXX(MELM),EYY(MELM),EZZ(MELM),
     *          EXY(MELM),EYZ(MELM),EXZ(MELM)
C
C[OUTPUT]
      REAL*4    A(N1,N2,NE),RHS(NP),AD(NP)
C
C[LOCAL]
      INTEGER*4 NETET,NEPRD,NEWED,NEHEX
      INTEGER*4 NTET,NPRD,NWED,NHEX
      INTEGER*4 NSKIP1,NSKIP2,NSKIP3,NSKIP4
      INTEGER*4 IE,IP,I,IERR,
     *          IES1,IEE1,IES2,IEE2,IES3,IEE3,IES4,IEE4,
     *          IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8,IELM0,
     *          IELM1A,IELM2A,IELM3A,IELM4A,
     *          IELM5A,IELM6A,IELM7A,IELM8A,
     *          IELM1B,IELM2B,IELM3B,IELM4B,
     *          IELM5B,IELM6B,IELM7B,IELM8B,
     *          ITRANS,IDIAGV,J
      REAL*4    UU,VV,WW,AT,DTH,UABS,DELTAP,
     *          AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8,
     *          AC1,AC2,AC3,AC4,AC5,AC6,AC7,AC8,
     *          AT1,AT2,AT3,AT4,AT5,AT6,AT7,AT8,
     *          AUPWD1,AUPWD2,AUPWD3,AUPWD4,
     *          AUPWD5,AUPWD6,AUPWD7,AUPWD8,
     *          CRHS1,CRHS2,CRHS3,CRHS4,CRHS5,CRHS6,CRHS7,CRHS8
C 
      REAL*4    DIJ(8,8)
      DATA DIJ /1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &          0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &          0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &          0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     &          0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
     &          0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     &          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
     &          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/
C
      REAL*4 CUPWD,CVISL
      INTEGER*4 ICAVI
      REAL*4 COEF
C
C          SCT         ; SGS SCHMIDT NUMBER FOR VOLUMETRIC FRACTION
C
      NETET =NEX( 1)
      NEPRD =NEX( 2)
      NEWED =NEX( 3)
      NEHEX =NEX( 4)
      NTET  =NEX( 5)
      NPRD  =NEX( 6)
      NWED  =NEX( 7)
      NHEX  =NEX( 8)
      NSKIP1=NEX( 9)
      NSKIP2=NEX(10)
      NSKIP3=NEX(11)
      NSKIP4=NEX(12)
C
C   == TET. ==  
      IES1=1
      IEE1=NETET 
C
C   == PYRAMID ==  
      IES2=NETET+1
      IEE2=NETET+NEPRD
C
C   == WEDGE ==  
      IES3=NETET+NEPRD+1
      IEE3=NETET+NEPRD+NEWED
C
C   == HEX. ==  
      IES4=NETET+NEPRD+NEWED+1
      IEE4=NETET+NEPRD+NEWED+NEHEX 
C
      DO 100 IP=1,NP
          AD (IP)=0.0E0
          RHS(IP)=0.0E0
  100 CONTINUE
C
      DO ICOLOR=1,NCOLOR(1)
!ocl norecurrence(A,AD,RHS)
      DO ICPART=1,NCPART(ICOLOR,1)
          IES1=LLOOP(ICPART  ,ICOLOR,1)
          IEE1=LLOOP(ICPART+1,ICOLOR,1)-1
!ocl nosimd
!ocl noswp
      DO 1000 IE=IES1,IEE1
          DO 1100 I=1,NTET
              DTH = DT*0.5E0
              IP =NODE(I,IE)
              IP1=NODE(1,IE)
              IP2=NODE(2,IE)
              IP3=NODE(3,IE)
              IP4=NODE(4,IE)
C
              UU=UG(IE)
              VV=VG(IE)
              WW=WG(IE)
              UABS=SQRT(UU*UU+VV*VV+WW*WW+1.0E-10)
C
              IELM0 =1+NSKIP1*(IE-1)
              IELM1A=IELM0+(I-1)+NTET*0
              IELM2A=IELM0+(I-1)+NTET*1
              IELM3A=IELM0+(I-1)+NTET*2
              IELM4A=IELM0+(I-1)+NTET*3
              IELM1B=IELM0+0+NTET*(I-1)
              IELM2B=IELM0+1+NTET*(I-1)
              IELM3B=IELM0+2+NTET*(I-1)
              IELM4B=IELM0+3+NTET*(I-1)
C
CC [TET. TIME-TERM]
              AT  = E(IELM1A)+E(IELM2A)+E(IELM3A)+E(IELM4A) 
              AT1 = AT*DIJ(I,1)
              AT2 = AT*DIJ(I,2)
              AT3 = AT*DIJ(I,3)
              AT4 = AT*DIJ(I,4)
C
CC [TET. COVECTION-TERM]
              AC1 = (UU*EX(IELM1A)+VV*EY(IELM1A)+WW*EZ(IELM1A))*DTH
              AC2 = (UU*EX(IELM2A)+VV*EY(IELM2A)+WW*EZ(IELM2A))*DTH
              AC3 = (UU*EX(IELM3A)+VV*EY(IELM3A)+WW*EZ(IELM3A))*DTH
              AC4 = (UU*EX(IELM4A)+VV*EY(IELM4A)+WW*EZ(IELM4A))*DTH
C
CC [TET. UPWIND-TERM]
              AUPWD1 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM1A)+VV*EXY(IELM1A)+WW*EXZ(IELM1A))
     &           +VV*(UU*EXY(IELM1B)+VV*EYY(IELM1A)+WW*EYZ(IELM1A))
     &           +WW*(UU*EXZ(IELM1B)+VV*EYZ(IELM1B)+WW*EZZ(IELM1A)))*DTH
              AUPWD2 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM2A)+VV*EXY(IELM2A)+WW*EXZ(IELM2A))
     &           +VV*(UU*EXY(IELM2B)+VV*EYY(IELM2A)+WW*EYZ(IELM2A))
     &           +WW*(UU*EXZ(IELM2B)+VV*EYZ(IELM2B)+WW*EZZ(IELM2A)))*DTH
              AUPWD3 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM3A)+VV*EXY(IELM3A)+WW*EXZ(IELM3A))
     &           +VV*(UU*EXY(IELM3B)+VV*EYY(IELM3A)+WW*EYZ(IELM3A))
     &           +WW*(UU*EXZ(IELM3B)+VV*EYZ(IELM3B)+WW*EZZ(IELM3A)))*DTH
              AUPWD4 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM4A)+VV*EXY(IELM4A)+WW*EXZ(IELM4A))
     &           +VV*(UU*EXY(IELM4B)+VV*EYY(IELM4A)+WW*EYZ(IELM4A))
     &           +WW*(UU*EXZ(IELM4B)+VV*EYZ(IELM4B)+WW*EZZ(IELM4A)))*DTH
              AUPWD1 = AUPWD1*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD2 = AUPWD2*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD3 = AUPWD3*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD4 = AUPWD4*(DELTA(IE)**0.33333E0)/(UABS*DT)
C
CC [TET. DIFFUSION-TERM]
              COEF= CVISL*VISC(IE)*DTH/SCT 
              AD1 = (EXX(IELM1A)+EYY(IELM1A)+EZZ(IELM1A))*COEF
              AD2 = (EXX(IELM2A)+EYY(IELM2A)+EZZ(IELM2A))*COEF
              AD3 = (EXX(IELM3A)+EYY(IELM3A)+EZZ(IELM3A))*COEF
              AD4 = (EXX(IELM4A)+EYY(IELM4A)+EZZ(IELM4A))*COEF
C
              A(I,1,IE) = AT1+AC1*(1.0E0+CEB)+AD1*(1.0E0+CEB)+AUPWD1
              A(I,2,IE) = AT2+AC2*(1.0E0+CEB)+AD2*(1.0E0+CEB)+AUPWD2
              A(I,3,IE) = AT3+AC3*(1.0E0+CEB)+AD3*(1.0E0+CEB)+AUPWD3
              A(I,4,IE) = AT4+AC4*(1.0E0+CEB)+AD4*(1.0E0+CEB)+AUPWD4
              AD(IP)=AD(IP)+A(I,I,IE)
C
              CRHS1 = AT1-AC1*(1.0E0-CEB)-AD1*(1.0E0-CEB)-AUPWD1
              CRHS2 = AT2-AC2*(1.0E0-CEB)-AD2*(1.0E0-CEB)-AUPWD2
              CRHS3 = AT3-AC3*(1.0E0-CEB)-AD3*(1.0E0-CEB)-AUPWD3
              CRHS4 = AT4-AC4*(1.0E0-CEB)-AD4*(1.0E0-CEB)-AUPWD4
C
              RHS(IP) = RHS(IP)
     &                 +CRHS1*X(IP1)+CRHS2*X(IP2)
     &                 +CRHS3*X(IP3)+CRHS4*X(IP4)
 1100     CONTINUE
 1000 CONTINUE   
      END DO
      END DO
C
C  /// PYRAMID ///
C
      DO ICOLOR=1,NCOLOR(2)
!ocl norecurrence(A,AD,RHS)
      DO ICPART=1,NCPART(ICOLOR,2)
          IES2=LLOOP(ICPART  ,ICOLOR,2)
          IEE2=LLOOP(ICPART+1,ICOLOR,2)-1
!ocl nosimd
!ocl noswp
      DO 2000 IE=IES2,IEE2
          DO 2100 I=1,NPRD
              DTH = DT*0.5E0
              IP =NODE(I,IE)
              IP1=NODE(1,IE)
              IP2=NODE(2,IE)
              IP3=NODE(3,IE)
              IP4=NODE(4,IE)
              IP5=NODE(5,IE)
C
              UU=UG(IE)
              VV=VG(IE)
              WW=WG(IE)
              UABS=SQRT(UU*UU+VV*VV+WW*WW+1.0E-10)
C
              IELM0 = NETET*NSKIP1
     *               +1+NSKIP2*(IE-1-IEE1)
C
              IELM1A=IELM0+(I-1)+NPRD*0
              IELM2A=IELM0+(I-1)+NPRD*1
              IELM3A=IELM0+(I-1)+NPRD*2
              IELM4A=IELM0+(I-1)+NPRD*3
              IELM5A=IELM0+(I-1)+NPRD*4
              IELM1B=IELM0+0+NPRD*(I-1)
              IELM2B=IELM0+1+NPRD*(I-1)
              IELM3B=IELM0+2+NPRD*(I-1)
              IELM4B=IELM0+3+NPRD*(I-1)
              IELM5B=IELM0+4+NPRD*(I-1)
C
CC [PRD. TIME-TERM]
              AT =E(IELM1A)+E(IELM2A)+E(IELM3A)+E(IELM4A) 
     *           +E(IELM5A)
              AT1 = AT*DIJ(I,1)
              AT2 = AT*DIJ(I,2)
              AT3 = AT*DIJ(I,3)
              AT4 = AT*DIJ(I,4)
              AT5 = AT*DIJ(I,5)
C
CC [PRD. COVECTION-TERM]
              AC1 = (UU*EX(IELM1A)+VV*EY(IELM1A)+WW*EZ(IELM1A))*DTH
              AC2 = (UU*EX(IELM2A)+VV*EY(IELM2A)+WW*EZ(IELM2A))*DTH
              AC3 = (UU*EX(IELM3A)+VV*EY(IELM3A)+WW*EZ(IELM3A))*DTH
              AC4 = (UU*EX(IELM4A)+VV*EY(IELM4A)+WW*EZ(IELM4A))*DTH
              AC5 = (UU*EX(IELM5A)+VV*EY(IELM5A)+WW*EZ(IELM5A))*DTH
C
CC [PRD. UPWIND-TERM]
              AUPWD1 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM1A)+VV*EXY(IELM1A)+WW*EXZ(IELM1A))
     &           +VV*(UU*EXY(IELM1B)+VV*EYY(IELM1A)+WW*EYZ(IELM1A))
     &           +WW*(UU*EXZ(IELM1B)+VV*EYZ(IELM1B)+WW*EZZ(IELM1A)))*DTH
              AUPWD2 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM2A)+VV*EXY(IELM2A)+WW*EXZ(IELM2A))
     &           +VV*(UU*EXY(IELM2B)+VV*EYY(IELM2A)+WW*EYZ(IELM2A))
     &           +WW*(UU*EXZ(IELM2B)+VV*EYZ(IELM2B)+WW*EZZ(IELM2A)))*DTH
              AUPWD3 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM3A)+VV*EXY(IELM3A)+WW*EXZ(IELM3A))
     &           +VV*(UU*EXY(IELM3B)+VV*EYY(IELM3A)+WW*EYZ(IELM3A))
     &           +WW*(UU*EXZ(IELM3B)+VV*EYZ(IELM3B)+WW*EZZ(IELM3A)))*DTH
              AUPWD4 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM4A)+VV*EXY(IELM4A)+WW*EXZ(IELM4A))
     &           +VV*(UU*EXY(IELM4B)+VV*EYY(IELM4A)+WW*EYZ(IELM4A))
     &           +WW*(UU*EXZ(IELM4B)+VV*EYZ(IELM4B)+WW*EZZ(IELM4A)))*DTH
              AUPWD5 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM5A)+VV*EXY(IELM5A)+WW*EXZ(IELM5A))
     &           +VV*(UU*EXY(IELM5B)+VV*EYY(IELM5A)+WW*EYZ(IELM5A))
     &           +WW*(UU*EXZ(IELM5B)+VV*EYZ(IELM5B)+WW*EZZ(IELM5A)))*DTH
              AUPWD1 = AUPWD1*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD2 = AUPWD2*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD3 = AUPWD3*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD4 = AUPWD4*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD5 = AUPWD5*(DELTA(IE)**0.33333E0)/(UABS*DT)
C
CC [PRD. DIFFUSION-TERM]
              COEF= CVISL*VISC(IE)*DTH/SCT 
              AD1 = (EXX(IELM1A)+EYY(IELM1A)+EZZ(IELM1A))*COEF
              AD2 = (EXX(IELM2A)+EYY(IELM2A)+EZZ(IELM2A))*COEF
              AD3 = (EXX(IELM3A)+EYY(IELM3A)+EZZ(IELM3A))*COEF
              AD4 = (EXX(IELM4A)+EYY(IELM4A)+EZZ(IELM4A))*COEF
              AD5 = (EXX(IELM5A)+EYY(IELM5A)+EZZ(IELM5A))*COEF
C
              A(I,1,IE) = AT1+AC1*(1.0E0+CEB)+AD1*(1.0E0+CEB)+AUPWD1
              A(I,2,IE) = AT2+AC2*(1.0E0+CEB)+AD2*(1.0E0+CEB)+AUPWD2
              A(I,3,IE) = AT3+AC3*(1.0E0+CEB)+AD3*(1.0E0+CEB)+AUPWD3
              A(I,4,IE) = AT4+AC4*(1.0E0+CEB)+AD4*(1.0E0+CEB)+AUPWD4
              A(I,5,IE) = AT5+AC5*(1.0E0+CEB)+AD5*(1.0E0+CEB)+AUPWD5
              AD(IP)=AD(IP)+A(I,I,IE)
C
              CRHS1 = AT1-AC1*(1.0E0-CEB)-AD1*(1.0E0-CEB)-AUPWD1
              CRHS2 = AT2-AC2*(1.0E0-CEB)-AD2*(1.0E0-CEB)-AUPWD2
              CRHS3 = AT3-AC3*(1.0E0-CEB)-AD3*(1.0E0-CEB)-AUPWD3
              CRHS4 = AT4-AC4*(1.0E0-CEB)-AD4*(1.0E0-CEB)-AUPWD4
              CRHS5 = AT5-AC5*(1.0E0-CEB)-AD5*(1.0E0-CEB)-AUPWD5
C
              RHS(IP)=RHS(IP) 
     &               +CRHS1*X(IP1)+CRHS2*X(IP2)
     &               +CRHS3*X(IP3)+CRHS4*X(IP4)
     &               +CRHS5*X(IP5)
 2100     CONTINUE
 2000 CONTINUE   
      END DO
      END DO
C
C  /// WEDGE ///
C
      DO ICOLOR=1,NCOLOR(3)
!ocl norecurrence(A,AD,RHS)
      DO ICPART=1,NCPART(ICOLOR,3)
          IES3=LLOOP(ICPART  ,ICOLOR,3)
          IEE3=LLOOP(ICPART+1,ICOLOR,3)-1
!ocl nosimd
!ocl noswp
      DO 3000 IE=IES3,IEE3
          DO 3100 I=1,NWED
              DTH = DT*0.5E0
              IP =NODE(I,IE)
              IP1=NODE(1,IE)
              IP2=NODE(2,IE)
              IP3=NODE(3,IE)
              IP4=NODE(4,IE)
              IP5=NODE(5,IE)
              IP6=NODE(6,IE)
C
              UU=UG(IE)
              VV=VG(IE)
              WW=WG(IE)
              UABS=SQRT(UU*UU+VV*VV+WW*WW+1.0E-10)
C
              IELM0 = NETET*NSKIP1
     *               +NEPRD*NSKIP2
     *               +1+NSKIP3*(IE-1-IEE2)
C
              IELM1A=IELM0+(I-1)+NWED*0
              IELM2A=IELM0+(I-1)+NWED*1
              IELM3A=IELM0+(I-1)+NWED*2
              IELM4A=IELM0+(I-1)+NWED*3
              IELM5A=IELM0+(I-1)+NWED*4
              IELM6A=IELM0+(I-1)+NWED*5
              IELM1B=IELM0+0+NWED*(I-1)
              IELM2B=IELM0+1+NWED*(I-1)
              IELM3B=IELM0+2+NWED*(I-1)
              IELM4B=IELM0+3+NWED*(I-1)
              IELM5B=IELM0+4+NWED*(I-1)
              IELM6B=IELM0+5+NWED*(I-1)
C
CC [WED. TIME-TERM]
              AT =E(IELM1A)+E(IELM2A)+E(IELM3A)+E(IELM4A) 
     *           +E(IELM5A)+E(IELM6A)
              AT1 = AT*DIJ(I,1)
              AT2 = AT*DIJ(I,2)
              AT3 = AT*DIJ(I,3)
              AT4 = AT*DIJ(I,4)
              AT5 = AT*DIJ(I,5)
              AT6 = AT*DIJ(I,6)
C
CC [WED. COVECTION-TERM]
              AC1 = (UU*EX(IELM1A)+VV*EY(IELM1A)+WW*EZ(IELM1A))*DTH
              AC2 = (UU*EX(IELM2A)+VV*EY(IELM2A)+WW*EZ(IELM2A))*DTH
              AC3 = (UU*EX(IELM3A)+VV*EY(IELM3A)+WW*EZ(IELM3A))*DTH
              AC4 = (UU*EX(IELM4A)+VV*EY(IELM4A)+WW*EZ(IELM4A))*DTH
              AC5 = (UU*EX(IELM5A)+VV*EY(IELM5A)+WW*EZ(IELM5A))*DTH
              AC6 = (UU*EX(IELM6A)+VV*EY(IELM6A)+WW*EZ(IELM6A))*DTH
C
CC [WED. UPWIND-TERM]
              AUPWD1 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM1A)+VV*EXY(IELM1A)+WW*EXZ(IELM1A))
     &           +VV*(UU*EXY(IELM1B)+VV*EYY(IELM1A)+WW*EYZ(IELM1A))
     &           +WW*(UU*EXZ(IELM1B)+VV*EYZ(IELM1B)+WW*EZZ(IELM1A)))*DTH
              AUPWD2 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM2A)+VV*EXY(IELM2A)+WW*EXZ(IELM2A))
     &           +VV*(UU*EXY(IELM2B)+VV*EYY(IELM2A)+WW*EYZ(IELM2A))
     &           +WW*(UU*EXZ(IELM2B)+VV*EYZ(IELM2B)+WW*EZZ(IELM2A)))*DTH
              AUPWD3 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM3A)+VV*EXY(IELM3A)+WW*EXZ(IELM3A))
     &           +VV*(UU*EXY(IELM3B)+VV*EYY(IELM3A)+WW*EYZ(IELM3A))
     &           +WW*(UU*EXZ(IELM3B)+VV*EYZ(IELM3B)+WW*EZZ(IELM3A)))*DTH
              AUPWD4 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM4A)+VV*EXY(IELM4A)+WW*EXZ(IELM4A))
     &           +VV*(UU*EXY(IELM4B)+VV*EYY(IELM4A)+WW*EYZ(IELM4A))
     &           +WW*(UU*EXZ(IELM4B)+VV*EYZ(IELM4B)+WW*EZZ(IELM4A)))*DTH
              AUPWD5 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM5A)+VV*EXY(IELM5A)+WW*EXZ(IELM5A))
     &           +VV*(UU*EXY(IELM5B)+VV*EYY(IELM5A)+WW*EYZ(IELM5A))
     &           +WW*(UU*EXZ(IELM5B)+VV*EYZ(IELM5B)+WW*EZZ(IELM5A)))*DTH
              AUPWD6 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM6A)+VV*EXY(IELM6A)+WW*EXZ(IELM6A))
     &           +VV*(UU*EXY(IELM6B)+VV*EYY(IELM6A)+WW*EYZ(IELM6A))
     &           +WW*(UU*EXZ(IELM6B)+VV*EYZ(IELM6B)+WW*EZZ(IELM6A)))*DTH
              AUPWD1 = AUPWD1*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD2 = AUPWD2*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD3 = AUPWD3*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD4 = AUPWD4*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD5 = AUPWD5*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD6 = AUPWD6*(DELTA(IE)**0.33333E0)/(UABS*DT)
C
CC [WED. DIFFUSION-TERM]
              COEF= CVISL*VISC(IE)*DTH/SCT 
              AD1 = (EXX(IELM1A)+EYY(IELM1A)+EZZ(IELM1A))*COEF
              AD2 = (EXX(IELM2A)+EYY(IELM2A)+EZZ(IELM2A))*COEF
              AD3 = (EXX(IELM3A)+EYY(IELM3A)+EZZ(IELM3A))*COEF
              AD4 = (EXX(IELM4A)+EYY(IELM4A)+EZZ(IELM4A))*COEF
              AD5 = (EXX(IELM5A)+EYY(IELM5A)+EZZ(IELM5A))*COEF
              AD6 = (EXX(IELM6A)+EYY(IELM6A)+EZZ(IELM6A))*COEF
C
              A(I,1,IE) = AT1+AC1*(1.0E0+CEB)+AD1*(1.0E0+CEB)+AUPWD1
              A(I,2,IE) = AT2+AC2*(1.0E0+CEB)+AD2*(1.0E0+CEB)+AUPWD2
              A(I,3,IE) = AT3+AC3*(1.0E0+CEB)+AD3*(1.0E0+CEB)+AUPWD3
              A(I,4,IE) = AT4+AC4*(1.0E0+CEB)+AD4*(1.0E0+CEB)+AUPWD4
              A(I,5,IE) = AT5+AC5*(1.0E0+CEB)+AD5*(1.0E0+CEB)+AUPWD5
              A(I,6,IE) = AT6+AC6*(1.0E0+CEB)+AD6*(1.0E0+CEB)+AUPWD6
              AD(IP)=AD(IP)+A(I,I,IE)
C
              CRHS1 = AT1-AC1*(1.0E0-CEB)-AD1*(1.0E0-CEB)-AUPWD1
              CRHS2 = AT2-AC2*(1.0E0-CEB)-AD2*(1.0E0-CEB)-AUPWD2
              CRHS3 = AT3-AC3*(1.0E0-CEB)-AD3*(1.0E0-CEB)-AUPWD3
              CRHS4 = AT4-AC4*(1.0E0-CEB)-AD4*(1.0E0-CEB)-AUPWD4
              CRHS5 = AT5-AC5*(1.0E0-CEB)-AD5*(1.0E0-CEB)-AUPWD5
              CRHS6 = AT6-AC6*(1.0E0-CEB)-AD6*(1.0E0-CEB)-AUPWD6
C
              RHS(IP)=RHS(IP)
     &               +CRHS1*X(IP1)+CRHS2*X(IP2)
     &               +CRHS3*X(IP3)+CRHS4*X(IP4)
     &               +CRHS5*X(IP5)+CRHS6*X(IP6)
 3100     CONTINUE
 3000 CONTINUE   
      END DO
      END DO
C
C  /// HEX ///
C
      DO ICOLOR=1,NCOLOR(4)
!ocl norecurrence(A,AD,RHS)
      DO ICPART=1,NCPART(ICOLOR,4)
          IES4=LLOOP(ICPART  ,ICOLOR,4)
          IEE4=LLOOP(ICPART+1,ICOLOR,4)-1
!ocl nosimd
!ocl noswp
      DO 4000 IE=IES4,IEE4
          DO 4100 I=1,NHEX
              DTH = DT*0.5E0
              IP =NODE(I,IE)
              IP1=NODE(1,IE)
              IP2=NODE(2,IE)
              IP3=NODE(3,IE)
              IP4=NODE(4,IE)
              IP5=NODE(5,IE)
              IP6=NODE(6,IE)
              IP7=NODE(7,IE)
              IP8=NODE(8,IE)
C
              UU=UG(IE)
              VV=VG(IE)
              WW=WG(IE)
              UABS=SQRT(UU*UU+VV*VV+WW*WW+1.0E-10)
              DELTAP=P(IE)+0.5E0*SIGMA
C
              IELM0 = NETET*NSKIP1
     *               +NEPRD*NSKIP2
     *               +NEWED*NSKIP3
     *               +1+NSKIP4*(IE-1-IEE3)
C
              IELM1A=IELM0+(I-1)+NHEX*0
              IELM2A=IELM0+(I-1)+NHEX*1
              IELM3A=IELM0+(I-1)+NHEX*2
              IELM4A=IELM0+(I-1)+NHEX*3
              IELM5A=IELM0+(I-1)+NHEX*4
              IELM6A=IELM0+(I-1)+NHEX*5
              IELM7A=IELM0+(I-1)+NHEX*6
              IELM8A=IELM0+(I-1)+NHEX*7
              IELM1B=IELM0+0+NHEX*(I-1)
              IELM2B=IELM0+1+NHEX*(I-1)
              IELM3B=IELM0+2+NHEX*(I-1)
              IELM4B=IELM0+3+NHEX*(I-1)
              IELM5B=IELM0+4+NHEX*(I-1)
              IELM6B=IELM0+5+NHEX*(I-1)
              IELM7B=IELM0+6+NHEX*(I-1)
              IELM8B=IELM0+7+NHEX*(I-1)
C
CC [HEX. TIME-TERM]
              AT =E(IELM1A)+E(IELM2A)+E(IELM3A)+E(IELM4A) 
     *           +E(IELM5A)+E(IELM6A)+E(IELM7A)+E(IELM8A)
              AT1 = AT*DIJ(I,1)
              AT2 = AT*DIJ(I,2)
              AT3 = AT*DIJ(I,3)
              AT4 = AT*DIJ(I,4)
              AT5 = AT*DIJ(I,5)
              AT6 = AT*DIJ(I,6)
              AT7 = AT*DIJ(I,7)
              AT8 = AT*DIJ(I,8)
C
CC [HEX. COVECTION-TERM]
              AC1 =(UU*EX(IELM1A)+VV*EY(IELM1A)+WW*EZ(IELM1A))*DTH
              AC2 =(UU*EX(IELM2A)+VV*EY(IELM2A)+WW*EZ(IELM2A))*DTH
              AC3 =(UU*EX(IELM3A)+VV*EY(IELM3A)+WW*EZ(IELM3A))*DTH
              AC4 =(UU*EX(IELM4A)+VV*EY(IELM4A)+WW*EZ(IELM4A))*DTH
              AC5 =(UU*EX(IELM5A)+VV*EY(IELM5A)+WW*EZ(IELM5A))*DTH
              AC6 =(UU*EX(IELM6A)+VV*EY(IELM6A)+WW*EZ(IELM6A))*DTH
              AC7 =(UU*EX(IELM7A)+VV*EY(IELM7A)+WW*EZ(IELM7A))*DTH
              AC8 =(UU*EX(IELM8A)+VV*EY(IELM8A)+WW*EZ(IELM8A))*DTH
C
CC [HEX. UPWIND-TERM]
              AUPWD1 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM1A)+VV*EXY(IELM1A)+WW*EXZ(IELM1A))
     &           +VV*(UU*EXY(IELM1B)+VV*EYY(IELM1A)+WW*EYZ(IELM1A))
     &           +WW*(UU*EXZ(IELM1B)+VV*EYZ(IELM1B)+WW*EZZ(IELM1A)))*DTH
              AUPWD2 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM2A)+VV*EXY(IELM2A)+WW*EXZ(IELM2A))
     &           +VV*(UU*EXY(IELM2B)+VV*EYY(IELM2A)+WW*EYZ(IELM2A))
     &           +WW*(UU*EXZ(IELM2B)+VV*EYZ(IELM2B)+WW*EZZ(IELM2A)))*DTH
              AUPWD3 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM3A)+VV*EXY(IELM3A)+WW*EXZ(IELM3A))
     &           +VV*(UU*EXY(IELM3B)+VV*EYY(IELM3A)+WW*EYZ(IELM3A))
     &           +WW*(UU*EXZ(IELM3B)+VV*EYZ(IELM3B)+WW*EZZ(IELM3A)))*DTH
              AUPWD4 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM4A)+VV*EXY(IELM4A)+WW*EXZ(IELM4A))
     &           +VV*(UU*EXY(IELM4B)+VV*EYY(IELM4A)+WW*EYZ(IELM4A))
     &           +WW*(UU*EXZ(IELM4B)+VV*EYZ(IELM4B)+WW*EZZ(IELM4A)))*DTH
              AUPWD5 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM5A)+VV*EXY(IELM5A)+WW*EXZ(IELM5A))
     &           +VV*(UU*EXY(IELM5B)+VV*EYY(IELM5A)+WW*EYZ(IELM5A))
     &           +WW*(UU*EXZ(IELM5B)+VV*EYZ(IELM5B)+WW*EZZ(IELM5A)))*DTH
              AUPWD6 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM6A)+VV*EXY(IELM6A)+WW*EXZ(IELM6A))
     &           +VV*(UU*EXY(IELM6B)+VV*EYY(IELM6A)+WW*EYZ(IELM6A))
     &           +WW*(UU*EXZ(IELM6B)+VV*EYZ(IELM6B)+WW*EZZ(IELM6A)))*DTH
              AUPWD7 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM7A)+VV*EXY(IELM7A)+WW*EXZ(IELM7A))
     &           +VV*(UU*EXY(IELM7B)+VV*EYY(IELM7A)+WW*EYZ(IELM7A))
     &           +WW*(UU*EXZ(IELM7B)+VV*EYZ(IELM7B)+WW*EZZ(IELM7A)))*DTH
              AUPWD8 = DTH*CUPWD
     &          *(UU*(UU*EXX(IELM8A)+VV*EXY(IELM8A)+WW*EXZ(IELM8A))
     &           +VV*(UU*EXY(IELM8B)+VV*EYY(IELM8A)+WW*EYZ(IELM8A))
     &           +WW*(UU*EXZ(IELM8B)+VV*EYZ(IELM8B)+WW*EZZ(IELM8A)))*DTH
              AUPWD1 = AUPWD1*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD2 = AUPWD2*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD3 = AUPWD3*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD4 = AUPWD4*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD5 = AUPWD5*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD6 = AUPWD6*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD7 = AUPWD7*(DELTA(IE)**0.33333E0)/(UABS*DT)
              AUPWD8 = AUPWD8*(DELTA(IE)**0.33333E0)/(UABS*DT)
C
CC [HEX. DIFFUSION-TERM]
              COEF= CVISL*VISC(IE)*DTH/SCT 
              AD1 = (EXX(IELM1A)+EYY(IELM1A)+EZZ(IELM1A))*COEF
              AD2 = (EXX(IELM2A)+EYY(IELM2A)+EZZ(IELM2A))*COEF
              AD3 = (EXX(IELM3A)+EYY(IELM3A)+EZZ(IELM3A))*COEF
              AD4 = (EXX(IELM4A)+EYY(IELM4A)+EZZ(IELM4A))*COEF
              AD5 = (EXX(IELM5A)+EYY(IELM5A)+EZZ(IELM5A))*COEF
              AD6 = (EXX(IELM6A)+EYY(IELM6A)+EZZ(IELM6A))*COEF
              AD7 = (EXX(IELM7A)+EYY(IELM7A)+EZZ(IELM7A))*COEF
              AD8 = (EXX(IELM8A)+EYY(IELM8A)+EZZ(IELM8A))*COEF
C
CC [HEX. MODEL-TERM]
C
              A(I,1,IE) = AT1+AC1*(1.0E0+CEB)
     *                       +AD1*(1.0E0+CEB)+AUPWD1
              A(I,2,IE) = AT2+AC2*(1.0E0+CEB)
     *                       +AD2*(1.0E0+CEB)+AUPWD2
              A(I,3,IE) = AT3+AC3*(1.0E0+CEB)
     *                       +AD3*(1.0E0+CEB)+AUPWD3
              A(I,4,IE) = AT4+AC4*(1.0E0+CEB)
     *                       +AD4*(1.0E0+CEB)+AUPWD4
              A(I,5,IE) = AT5+AC5*(1.0E0+CEB)
     *                       +AD5*(1.0E0+CEB)+AUPWD5
              A(I,6,IE) = AT6+AC6*(1.0E0+CEB)
     *                       +AD6*(1.0E0+CEB)+AUPWD6
              A(I,7,IE) = AT7+AC7*(1.0E0+CEB)
     *                       +AD7*(1.0E0+CEB)+AUPWD7
              A(I,8,IE) = AT8+AC8*(1.0E0+CEB)
     *                       +AD8*(1.0E0+CEB)+AUPWD8
              AD(IP)=AD(IP)+A(I,I,IE)
C
              IF(ICAVI.EQ.1) THEN
                  IF(DELTAP.LE.0.0E0) THEN
                      COEF = (CLQD-CGAS)*DT*DELTAP
                  ELSE
                      COEF = (CLQD-0.10E0*CGAS)*DT*DELTAP
                  ENDIF
              ELSEIF(ICAVI.EQ.2) THEN
                  IF(DELTAP.LE.0.0E0) THEN
                      COEF = CLQD*DT*DELTAP
                  ELSE
                      COEF = -CGAS*DT*DELTAP
                  ENDIF
              ELSEIF(ICAVI.EQ.3) THEN
                  COEF = 0.0E0
              ELSEIF(ICAVI.EQ.4)     THEN
                  COEF = 0.0E0
              ENDIF

              CRHS1 = (1.0E0+COEF)*AT1-AC1*(1.0E0-CEB)
     *              - AD1*(1.0E0-CEB)-AUPWD1
              CRHS2 = (1.0E0+COEF)*AT2-AC2*(1.0E0-CEB)
     *              - AD2*(1.0E0-CEB)-AUPWD2
              CRHS3 = (1.0E0+COEF)*AT3-AC3*(1.0E0-CEB)
     *              - AD3*(1.0E0-CEB)-AUPWD3
              CRHS4 = (1.0E0+COEF)*AT4-AC4*(1.0E0-CEB)
     *              - AD4*(1.0E0-CEB)-AUPWD4
              CRHS5 = (1.0E0+COEF)*AT5-AC5*(1.0E0-CEB)
     *              - AD5*(1.0E0-CEB)-AUPWD5
              CRHS6 = (1.0E0+COEF)*AT6-AC6*(1.0E0-CEB)
     *              - AD6*(1.0E0-CEB)-AUPWD6
              CRHS7 = (1.0E0+COEF)*AT7-AC7*(1.0E0-CEB)
     *              - AD7*(1.0E0-CEB)-AUPWD7
              CRHS8 = (1.0E0+COEF)*AT8-AC8*(1.0E0-CEB)
     *              - AD8*(1.0E0-CEB)-AUPWD8
C
              RHS(IP)=RHS(IP) 
     &               +CRHS1*X(IP1)+CRHS2*X(IP2)
     &               +CRHS3*X(IP3)+CRHS4*X(IP4)
     &               +CRHS5*X(IP5)+CRHS6*X(IP6)
     &               +CRHS7*X(IP7)+CRHS8*X(IP8)
 4100     CONTINUE
 4000 CONTINUE   
      END DO
      END DO
C
      RETURN
      END

