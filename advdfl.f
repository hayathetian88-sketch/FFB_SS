      SUBROUTINE ADVDFL(N1,N2,NE,NP,MEP,NEX,
     *                  UG,VG,WG,DT,SCT,VISC,
     *                  CEB,CUPWD,CVISL,
     *                  EAP1,EAP2,EAP3,EBP,
     *                  IENP,JENP,NODP,NEP,MP,
     *                  ICAVI,CGAS,CLQD,P,SIGMA,
     *                  AP,A,X,RHS,A0,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  LP2NPFREE,
     *                  DELTA,IERR)
      IMPLICIT NONE
C
C[INTPUT]
      INTEGER*4 N1,N2,NE,NP,MEP,NEX(12)
      INTEGER*4 ICAVI
      REAL(4),INTENT(IN) :: CGAS,CLQD,P(NE),SIGMA
      INTEGER*4 MP,IENP(MEP,MP),JENP(MEP,MP),
     *          NODP(N2,MEP,NP),NEP(MP)
      REAL*4    UG(NE),VG(NE),WG(NE),DT
      REAL*4    VISC(NE),SCT,CEB,X(NP+1),DELTA(NE)
      REAL*4    EAP1(N2,MEP,NP),EAP2(3,N2,MEP,NP),
     *          EAP3(6,N2,MEP,NP)
      REAL*4    EBP(3,N2,MEP,NP),AP(N2,MEP,NP)
C
      INTEGER*4 FLAG_CAV_XPOS
      INTEGER*4 LP2NPFREE(NP)
      REAL*4    CAV_XPOS,XPOS(NE)
C
C[OUTPUT]
      REAL*4    A(N1,N2,NE),RHS(NP),A0(NP)
C
C[LOCAL]
      INTEGER*4 NETET,NEPRD,NEWED,NEHEX,NN
      INTEGER*4 IE,IP,I,IERR,
     *          IEE1,IEE2,IEE3,IEE4,IPE,K,IPK
      REAL*4    UU,VV,WW,AT,DTH,UABS,
     *          ADI,ACI,ATI,AUPWDI,CRHSI,DELTAP,XE
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
      REAL*4 CUPWD
      REAL*4 CVISL
      REAL*4 COEF,COEFD
C
C          SCT         ; SGS SCHMIDT NUMBER FOR VOLUMETRIC FRACTION
C
      X(NP+1)=0.0E0
C
      NETET =NEX(1)
      NEPRD =NEX(2)
      NEWED =NEX(3)
      NEHEX =NEX(4)
C
C   == TET. ==  
      IEE1=NETET 
C
C   == PYRAMID ==  
      IEE2=NETET+NEPRD
C
C   == WEDGE ==  
      IEE3=NETET+NEPRD+NEWED
C
C   == HEX. ==  
      IEE4=NETET+NEPRD+NEWED+NEHEX 
C
      DO 100 IP=1,NP
        A0 (IP)=0.0E0
        RHS(IP)=0.0E0
  100 CONTINUE
C
      DTH = DT*0.5E0
C
!ocl norecurrence(A0,RHS)
!$omp parallel
!$omp+private(IP,IPE,IE,I,K,IPK,UU,VV,WW)
!$omp+private(COEF,COEFD,UABS,DELTAP)
!$omp+private(AT,ATI,ACI,AUPWDI,ADI,CRHSI,XE)
!$omp+shared(RHS,NEP,IENP,JENP,UG,VG,WG,P,X,DELTA,XPOS)
!$omp+shared(VISC,EAP1,EAP2,EAP3,EBP,DIJ,A0,AP,NODP)
!$omp+shared(LP2NPFREE)
!$omp+firstprivate(DT,DTH,SCT,CEB,CUPWD,CVISL)
!$omp+firstprivate(ICAVI,CGAS,CLQD,SIGMA)
!$omp+firstprivate(NP,FLAG_CAV_XPOS,CAV_XPOS)
!$omp do
      DO 1040 IP=1,NP
          CRHSI = 0.0E0
!ocl noswp
          DO 1030 IPE=1,NEP(IP)
              IE =IENP(IPE,IP)
              I  =JENP(IPE,IP)
              UU=UG(IE)
              VV=VG(IE)
              WW=WG(IE)
              XE=XPOS(IE)
              DELTAP=P(IE)+0.5E0*SIGMA
              UABS=SQRT(UU*UU+VV*VV+WW*WW+1.0E-10)
              COEFD = CVISL*VISC(IE)*DTH/SCT 
C
              AT =EAP1(1,IPE,IP)+EAP1(2,IPE,IP)+EAP1(3,IPE,IP)
     *           +EAP1(4,IPE,IP)+EAP1(5,IPE,IP)+EAP1(6,IPE,IP)
     *           +EAP1(7,IPE,IP)+EAP1(8,IPE,IP)

C
              DO 1020 K=1,8
C
CC [TIME-TERM]
                  ATI = AT*DIJ(K,I)
C
CC [COVECTION-TERM]
                  ACI =(UU*EAP2(1,K,IPE,IP)+VV*EAP2(2,K,IPE,IP)
     *                 +WW*EAP2(3,K,IPE,IP))*DTH
C
CC [UPWIND-TERM]
              AUPWDI = DTH*CUPWD
     *           *(UU*(UU*EAP3(1,K,IPE,IP)+VV*EAP3(4,K,IPE,IP)
     *            +WW*EAP3(6,K,IPE,IP))
     *            +VV*(UU*EBP(1,K,IPE,IP)+VV*EAP3(2,K,IPE,IP)
     *            +WW*EAP3(5,K,IPE,IP))
     *            +WW*(UU*EBP(3,K,IPE,IP)+VV*EBP(2,K,IPE,IP)
     *            +WW*EAP3(3,K,IPE,IP)))*DTH
              AUPWDI = AUPWDI*(DELTA(IE)**0.33333E0)/(UABS*DT)
C
CC [HEX. DIFFUSION-TERM]

              ADI = (EAP3(1,K,IPE,IP)+EAP3(2,K,IPE,IP)
     *              +EAP3(3,K,IPE,IP))*COEFD
C
              AP(K,IPE,IP) = ATI+ACI*(1.0E0+CEB)+ADI*(1.0E0+CEB)
     *                     + AUPWDI
C
              IPK=NODP(K,IPE,IP)
              COEF = 0.0
C
              IF(ICAVI.EQ.1) THEN
                  IF(FLAG_CAV_XPOS .EQ. 1)THEN
                  IF((DELTAP.LE.0).AND.(XE .GE. CAV_XPOS))THEN
                      COEF = (CLQD-CGAS)*DT*DELTAP
                  ELSE
                      COEF = (CLQD-0.10E0*CGAS)*DT*DELTAP
                  ENDIF
                  ELSE
                  IF(DELTAP.LE.0.0E0) THEN
                      COEF = (CLQD-CGAS)*DT*DELTAP
                  ELSE
                      COEF = (CLQD-0.10E0*CGAS)*DT*DELTAP
                  ENDIF
                  ENDIF
              ELSEIF(ICAVI.EQ.2) THEN
                  IF(DELTAP.LE.0.0E0) THEN
                      COEF = CLQD*DT*DELTAP
                  ELSE
                      COEF = -CGAS*DT*DELTAP
                  ENDIF
              ELSEIF(ICAVI.EQ.3) THEN
                  COEF = 0.0E0
              ELSEIF(ICAVI.EQ.4) THEN
                  COEF = 0.0E0
              ENDIF
              IF(LP2NPFREE(IP) .EQ. 1)THEN
                  COEF = 0.0E0
              ENDIF
C
C
              CRHSI = CRHSI 
     *         + ((1.0E0+COEF)*ATI-ACI*(1.0E0-CEB)
     *            -ADI*(1.0E0-CEB)-AUPWDI)*X(IPK)
C
 1020         CONTINUE
              A0(IP)=A0(IP)+AP(I,IPE,IP)
C
 1030     CONTINUE
          RHS(IP) = CRHSI
C
 1040 CONTINUE   
!$omp end do
!$omp end parallel
C
!$omp parallel
!$omp+private(IP,IPE,K,IE,I)
!$omp+firstprivate(NP)
!$omp+shared(NEP,IENP,JENP,A,AP)
!$omp do
      DO IP=1,NP
        DO IPE=1,NEP(IP)
          DO K=1,8 ! NHEX
            IE =IENP(IPE,IP)
            I  =JENP(IPE,IP)
            A(K,I,IE)=AP(K,IPE,IP)
          ENDDO
        ENDDO
      ENDDO
!$omp end do
!$omp end parallel

      !DO IP=1,NP
      !    IF (A0(IP).EQ.0.0E0) A0(IP)=1.0E0
      !ENDDO
      RETURN
      END

