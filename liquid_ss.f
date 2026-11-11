      SUBROUTINE LIQUID_SS(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *                  N1,N2,ME,MELM,NE,NP,MEP,NEX,NODE,NUMIP,
     *                  ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *                  NCRS,NCRS2,NPPMAX,ITPCRS,IUT6,IUT0,JUNROL,
     *                  NPP,IPCRS,LTAB,
     *                  LSTCLR,LSTDGN,NITR,RES,
     *                  EPS,EPSRE,SIGMA,CGAS,CLQD,COEF1,
     *                  FLINLT,FLIMIT,FLMIN,
     *                  EAP1,EAP2,EAP3,EBP,AP,
     *                  MP,IENP,JENP,NODP,NEP,
     *                  SN,DELTA,CM,
     *                  UG,VG,WG,P,DT,VISC,SCT,F,FESRC,FE,
     *                  NPINLT,LPINLT,
     *                  NPSET, LPSET1,LPSET2,LPSET3,
     *                  COVER1,COVER2,COVER3,
     *                  NPSND ,LPSND ,NPTSND,IPSET,IPSRC,
     *                  NPRCV ,LPRCV ,NPTRCV,
     *                  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *                  LFIX3D,LPFIX,RX,RY,MWRK,WRKN,ACRS,A,
     *                  WRK01,WRK02,WRK03,WRK04,WRK05,
     *                  WRK06,WRK07,WRK08,
     *                  RHS,A0,TS,
     *                  TACRS,IERR,ICRS_T,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  NEFREE,LEFREE,LP2NPFREE,
     *                  CEB,CUPWD,CVISL)
      IMPLICIT NONE
C
C[INPUT]
      INTEGER*4 MCOLOR,MCPART,
     *          NCOLOR(4),NCPART(MCOLOR,4),LLOOP(MCPART,MCOLOR,4),
     *          N1,N2,ME,MELM,NE,NP,MEP,NEX(12),NODE(N2,NE),NUMIP(NP),
     *          ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *          NCRS,NPP(NP),IPCRS(NCRS),LTAB(N1,N2,NE),
     *          LSTCLR(NCRS),LSTDGN(NCRS),
     *          NPPMAX,NCRS2,ITPCRS(NCRS2),IUT6,IUT0,JUNROL
      REAL*4    EPS,EPSRE,SIGMA,CGAS,CLQD,COEF1,FLINLT,FLIMIT,FLMIN
C     SS MODEL PARAMETERS (#CAV_MDL  4  CGAS  CLQD  F0 IN INPUT FILE)
C     CGAS  : NUCLEI NUMBER DENSITY N0   (NON-DIM., 1/L**3)
C     CLQD  : INITIAL NUCLEI DIAMETER D0 (NON-DIM., L)
C     COEF1 : COEFFICIENT IN SQRT TERM   (RECEIVED FROM F0)
      REAL*4    UG(NE),VG(NE),WG(NE),P(NE),DT,VISC(NE),SCT
      REAL*4    SN(N1,NE),DELTA(NE),CM(NP)
      REAL*4    EAP1(N2,MEP,NP),EAP2(3,N2,MEP,NP),EAP3(6,N2,MEP,NP)
      REAL*4    EBP(3,N2,MEP,NP),AP(N2,MEP,NP)
      INTEGER*4 MP,IENP(MEP,MP),JENP(MEP,MP),NODP(N2,MEP,NP),NEP(MP)
      INTEGER*4 NPINLT,LPINLT(NPINLT)
      INTEGER*4 IPART,NDOM,MBPDOM,
     *          LDOM(NDOM),NBPDOM(NDOM),
     *          IPSLF(MBPDOM,NDOM),IPSND(MBPDOM,NDOM)
      INTEGER*4 NPSND,NPRCV
      INTEGER*4 NPSET,LPSET1(NPSET),LPSET2(NPSET),LPSET3(NPSET),
     *          LPSND (NDOM),NPTSND(NDOM) ,LPRCV (NDOM) ,NPTRCV(NDOM),
     *          IPSET (MBPDOM,NDOM),IPSRC (MBPDOM,NDOM)
      REAL*4    COVER1(NPSET),COVER2(NPSET),COVER3(NPSET)
C
C[INPUT/OUTPUT]
      REAL*4    F(NP+1)
C
C[OUTPUT]
      INTEGER*4 IERR,NITR
      REAL*4    FESRC(NE),FE(NE)
      REAL*4    RES
C
C[WORK]
      INTEGER*4 LFIX3D(NP),LPFIX(NP),MWRK
      REAL*4    RX(0:N2,ME),RY(0:N2,ME),ACRS(NCRS),A(N1,N2,NE),
     *          WRK01(NP),WRK02(NP),WRK03(NP),WRK04(NP),
     *          WRK05(NP),WRK06(NP),
     *          WRK07(NP),WRK08(NP),WRKN(MWRK,6),
     *          RHS(NP),A0(NP),TS(0:NP),TACRS(NCRS2)
C
C[LOCAL]
      INTEGER*4 MAXBUF,IDUM,I,IE,IP,IB,ISEND,IDIM,NB,NPFIX,IBE
      INTEGER*4 IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8
      REAL*4    COEF,CEB,RESR,SOURCE,DELTAP,XE,VNUC,ANUC,ARGR,RR
      REAL*4    FLQ,GSS,SUSS,SPSS,ADD
C     SS MODEL PARAMETERS (LOCAL COPIES; DUMMY ARGS ARE NOT MODIFIED)
      REAL*4    ANB0,DIAM0
      REAL*4    PI
      PARAMETER (PI=3.14159266E0)

C
      INTEGER*4 IBCGS
      DATA IBCGS  / 0 /
C
      CHARACTER*60 ERMSGC
     & / ' ## SUBROUTINE LIQUID: FATAL      ERROR REPORT   ; RETURNED' /
C 
      INTEGER*4 ICRS,K
      INTEGER*4 ICRS_T(NP)
C
C     DEBUG FOR OLD PARAM
C
      REAL*4 CUPWD
      REAL*4 CVISL
C
      INTEGER*4 FLAG_CAV_XPOS
      REAL*4    CAV_XPOS,XPOS(NE)
      INTEGER*4 NEFREE
      INTEGER*4 LEFREE(2,NEFREE),LP2NPFREE(NP)
C
      IERR=0
C
C     SS MODEL PARAMETERS FROM INPUT (#CAV_MDL)
      ANB0  = CGAS
      DIAM0 = CLQD
      IF(ICAVI.EQ.4 .AND. (ANB0.LE.0.0E0 .OR. DIAM0.LE.0.0E0)) THEN
          WRITE(IUT0,*) ' ## LIQUID_SS: N0 (CGAS) AND D0 (CLQD) MUST BE'
          WRITE(IUT0,*) ' ## POSITIVE. N0, D0 =', ANB0, DIAM0
          IERR=1
          RETURN
      ENDIF
C
      IF(ITIME.EQ.0) GOTO 4000
      CALL USTSTA(124)
C
CCYY [1] PREPARE FIX B.C. LIST
C
      DO 1100 IP = 1 , NP
          LFIX3D(IP) = 0
          LPFIX(IP) = 0
 1100 CONTINUE
C
      DO 1110 IB = 1 , NPINLT
          LFIX3D(LPINLT(IB))=1
 1110 CONTINUE
C
      DO 1120 IB = 1 , NPSET
          ISEND = LPSET3(IB)
          IF(ISEND.GT.0) GO TO 1120
          LFIX3D(LPSET1(IB))=1
 1120 CONTINUE
C
      NPFIX=0
      DO 1200 IP = 1 , NP
          IF(LFIX3D(IP).EQ.0) GO TO 1200
          NPFIX=NPFIX+1
          LPFIX(NPFIX) = IP
 1200 CONTINUE
C
      IERR=0
      MAXBUF = ME*(N2+1)
      CALL USTEND(124)
C
CCYY [2] SS SOURCE TERM: SPLIT INTO EXPLICIT (SU) AND IMPLICIT (SP)
C
C     S = SIGN(DP)*G*F*(1+ANUC-F),  G = 3*RR*SQRT(2*|DP|*COEF1/3)>=0
C     (ANUC: NUCLEI VOLUME FRACTION. IT GIVES NON-ZERO SOURCE AT F=1,
C      I.E. CAVITATION INCEPTION FROM PURE LIQUID)
C       DP < 0 (EVAPORATION) : S = -[G*(1+ANUC-F^n)]*F^{n+1}
C                              SU = 0              , SP = G*(1+ANUC-F^n)
C       DP > 0 (CONDENSATION): S = G*F^n*(1+ANUC-F^{n+1})
C                              SU = G*F^n*(1+ANUC) , SP = G*F^n
C     S^{n+1} = SU - SP*F^{n+1},  SP >= 0  -> DIAGONAL IS ENHANCED
C     SU IS STORED IN FESRC, SP IS TEMPORARILY STORED IN FE
C     (FE IS RECOMPUTED IN [10] BEFORE RETURN)
C
      CALL USTSTA(127)
      DO 2000 IE = 1 , NE
          IP1=NODE(1,IE)
          IP2=NODE(2,IE)
          IP3=NODE(3,IE)
          IP4=NODE(4,IE)
          IP5=NODE(5,IE)
          IP6=NODE(6,IE)
          IP7=NODE(7,IE)
          IP8=NODE(8,IE)
          IF(IP8.NE.0) THEN 
              FLQ = ( F(IP1)+F(IP2)+F(IP3)+F(IP4)
     *               +F(IP5)+F(IP6)+F(IP7)+F(IP8))/8.0E0
          ELSE IF(IP6.NE.0) THEN
              FLQ = ( F(IP1)+F(IP2)+F(IP3)+F(IP4)
     *               +F(IP5)+F(IP6)              )/6.0E0
          ELSE IF(IP5.NE.0) THEN
              FLQ = ( F(IP1)+F(IP2)+F(IP3)+F(IP4)
     *               +F(IP5)                     )/5.0E0
          ELSE
              FLQ = ( F(IP1)+F(IP2)+F(IP3)+F(IP4))/4.0E0
          ENDIF 
          FLQ = MAX(0.0E0, MIN(1.0E0, FLQ))
C
          SUSS = 0.0E0
          SPSS = 0.0E0
          IF(ICAVI.EQ.4) THEN
              DELTAP = P(IE)+0.5E0*SIGMA
              VNUC = (PI/6.0E0)*ANB0*DIAM0*DIAM0*DIAM0
              ANUC = VNUC/(1.0E0+VNUC)
              ARGR = (4.0E0*PI*ANB0*FLQ)
     *             / (3.0E0*(1.0E0+ANUC-FLQ))
              RR   = ARGR**(0.33333333E0)
              GSS  = 3.0E0*RR*SQRT(2.0E0*ABS(DELTAP)*COEF1/3.0E0)
              IF(DELTAP.LT.0.0E0) THEN
                  SPSS = GSS*(1.0E0+ANUC-FLQ)
              ELSE IF(FLQ.LT.FLIMIT) THEN
                  SUSS = GSS*FLQ*(1.0E0+ANUC)
                  SPSS = GSS*FLQ
              ENDIF
          ENDIF
          FESRC(IE) = SUSS
          FE   (IE) = SPSS
 2000 CONTINUE
      IF(NEFREE.GE.1)THEN
      DO IBE = 1 , NEFREE
          IE = LEFREE(1,IBE)
          IF(IE .GE. 1)THEN
              FESRC(IE) = 0.0E0
              FE   (IE) = 0.0E0
          ENDIF
      END DO
      ENDIF
      CALL USTEND(127)
C
CCYY [3] CAL. TIME, ADV., AND DIFF. TERM + IMPLICIT SOURCE (SP)
C
      CALL USTSTA(125)
      CALL ADVDFL(N1,N2,NE,NP,MEP,NEX,
     *            UG,VG,WG,DT,SCT,VISC,
     *            CEB,CUPWD,CVISL,
     *            EAP1,EAP2,EAP3,EBP,
     *            IENP,JENP,NODP,NEP,MP,
     *            ICAVI,CGAS,CLQD,P,SIGMA,
     *            AP,A,F,RHS,A0,
     *            FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *            LP2NPFREE,
     *            DELTA,IERR)
C
C     ADD DT*SP*SN TO ELEMENT DIAGONAL (SAME WEIGHT SN AS NODLEX)
C     (BEFORE DDCOMX SO THAT A0 IS SUMMED OVER SUB-DOMAINS)
      DO 2050 IE = 1 , NE
          IF(FE(IE).LE.0.0E0) GOTO 2050
          DO 2040 I = 1 , N2
              IP = NODE(I,IE)
              IF(IP.LE.0) GOTO 2040
              IF(LP2NPFREE(IP).EQ.1) GOTO 2040
              ADD = DT*FE(IE)*SN(I,IE)
              A(I,I,IE) = A(I,I,IE) + ADD
              A0(IP)    = A0(IP)    + ADD
 2040     CONTINUE
 2050 CONTINUE
      CALL USTEND(125)
C
      CALL USTSTA(126)
      IDUM = 2
      CALL DDCOMX(IPART,IDUM,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *            RHS,A0,RHS,NP,IUT0,IERR,RX,RY,MAXBUF)
      IF(IERR.NE.0) THEN
          WRITE(IUT0,*)
          WRITE(IUT0,*) ERMSGC
          RETURN
      ENDIF
      CALL USTEND(126)
C
      CALL USTSTA(128)
      CALL NODLEX(NODE,ME,NP,NE,NP,N1,N2,NEX,SN,
     *            IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *            FESRC,WRK01,CM,IUT0,IERR,RX,RY,MAXBUF)
C
      DO 2100 IP = 1,NP
          RHS(IP)=RHS(IP)+DT*WRK01(IP)/CM(IP)
 2100 CONTINUE
      CALL USTEND(128)
C
CCYY [4] MATIRIX CONVET FROM ELEMENT-WIZE TO CRS
C
      CALL USTSTA(129)
      CALL E2PMATR(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *             N2,N1,NE,NEX,NCRS,A,ACRS,LTAB,IUT0,IERR)
      IF(IERR.NE.0) THEN
          WRITE(IUT0,*)'ERROR CODE REPORTED FROM E2PMATR'
          WRITE(IUT0,*)ERMSGC
          RETURN
      ENDIF
      CALL USTEND(129)
C
CCYY [5] DIAGONAL SCALING
C
      CALL USTSTA(130)
      IF(ITIME.EQ.1) THEN
         CALL STLCLR(LSTCLR,LSTDGN,NP,NCRS,IPCRS,NPP,LFIX3D)
      ENDIF
      CALL DGNSCL2(ACRS,A0,NP,NCRS,LSTDGN)
C
C     OLD
C      CALL DGNSCL(ACRS,A0,NP,NE,NCRS,IPCRS,NPP,ME)
C
      DO 1300 IP=1, NP
          RHS(IP) = RHS(IP)/A0(IP)
 1300 CONTINUE
      CALL USTEND(130)
C
CCYY [6] CLEAR CRS MATRIX FOR DIRICHLET B.C.
C
      CALL USTSTA(131)
      CALL CLRCRS2(ACRS,NP,NCRS,LSTCLR,NUMIP)
C
C    OLD
C
C      CALL CLRCRS(IPART,ACRS,NP,NCRS,IPCRS,NPP,
C     *            A0,LFIX3D,NUMIP,WRK01)
      CALL USTEND(131)
C
CCYY [7] SET B.C. (1)INLET, (2)OVERSET
C
      IF (JSET.GE.1) THEN
          CALL OVRST1(IPART,NPSET,N1,N2,ME,NE,NP,NEX,NODE,F,
     *                LPSET1,LPSET2,LPSET3,COVER1,COVER2,COVER3,
     *                NDOM,MBPDOM,NPSND,NPRCV,
     *                LPSND,NPTSND,LPRCV,NPTRCV,IPSET,IPSRC,
     *                WRK01,WRK02,RX,RY,IUT0,IERR)
          IF(IERR.NE.0) THEN
              WRITE(IUT0,*)'ERROR CODE REPORTED FROM OVRST1'
              WRITE(IUT0,*)ERMSGC
              RETURN
          ENDIF
      ENDIF
C
      DO 1400 IB=1,NPINLT
          F(LPINLT(IB)) = FLINLT
 1400 CONTINUE   

      DO 1500 IB = 1 , NPFIX
          IP=LPFIX(IB)
          RHS(IP) = F(IP)
 1500 CONTINUE

C
CCYY [8] SOLVE EQUATION
C
      CALL USTSTA(132)
      IF(JUNROL.EQ.1)
     *CALL CRSCVA(NP,NPPMAX,NCRS,NCRS2,NPP,ACRS,TACRS,ICRS_T)


      IF(ISOLT.EQ.1) THEN
          CALL BCGS3X(IBCGS,IPART,NMAX,EPS,ME,N2,NE,NEX,NP,
     *                NPP,NCRS,IPCRS,ACRS,NODE,RHS,F,NITR,RES,
     *                NDOM,MBPDOM,LDOM,NBPDOM,IPSLF,IPSND,NUMIP,
     *                WRK01,WRK02,WRK03,WRK04,WRK05,
     *                WRK06,WRK07,WRK08,RX,RY,
     *                MWRK,WRKN,
     *                IUT0,IERR,A0,LFIX3D,
     *                JUNROL,NPPMAX,NCRS2,TS,TACRS,ITPCRS)
      ELSE
              WRITE(IUT0,*) 'INVALID VALUE FOR ISOLT'  ,ISOLT
              IERR=1
              RETURN 
      ENDIF
C
      IF(IERR.NE.0) THEN
          WRITE(IUT0,*)'ERROR CODE REPORTED FROM BCGSTX'
          WRITE(IUT0,*)ERMSGC
          RETURN
      ENDIF
      CALL USTEND(132)

C
C
CCYY [8]  SET B.C. (1)INLET, (2)OVERSET
C
 4000 CONTINUE
      CALL USTSTA(133)
C
      DO 4100 IP = 1 , NP
          F(IP)  = AMAX1(FLMIN,F(IP))
          F(IP)  = AMIN1(1.E0 ,F(IP))
 4100 CONTINUE
C
      IF (JSET.GE.1) THEN
          CALL OVRST1(IPART,NPSET,N1,N2,ME,NE,NP,NEX,NODE,F,
     *                LPSET1,LPSET2,LPSET3,COVER1,COVER2,COVER3,
     *                NDOM,MBPDOM,NPSND,NPRCV,
     *                LPSND,NPTSND,LPRCV,NPTRCV,IPSET,IPSRC,
     *                WRK01,WRK02,RX,RY,IUT0,IERR)
          IF(IERR.NE.0) THEN
              WRITE(IUT0,*)'ERROR CODE REPORTED FROM OVRST1'
              WRITE(IUT0,*)ERMSGC
              RETURN
          ENDIF
      ENDIF
C
      DO 4200 IB=1,NPINLT
          F(LPINLT(IB)) = FLINLT
 4200 CONTINUE   
C
      DO 4300 IE = 1 , NE
          XE = XPOS(IE)
          IP1=NODE(1,IE)
          IP2=NODE(2,IE)
          IP3=NODE(3,IE)
          IP4=NODE(4,IE)
          IP5=NODE(5,IE)
          IP6=NODE(6,IE)
          IP7=NODE(7,IE)
          IP8=NODE(8,IE)
          IF(IP8.NE.0) THEN 
              FE(IE) = ( F(IP1)+F(IP2)+F(IP3)+F(IP4)
     *                  +F(IP5)+F(IP6)+F(IP7)+F(IP8))/8.0E0
          ELSE IF(IP6.NE.0) THEN
              FE(IE) = ( F(IP1)+F(IP2)+F(IP3)+F(IP4)
     *                  +F(IP5)+F(IP6)              )/6.0E0
          ELSE IF(IP5.NE.0) THEN
              FE(IE) = ( F(IP1)+F(IP2)+F(IP3)+F(IP4)
     *                  +F(IP5)                     )/5.0E0
          ELSE IF(IP4.NE.0) THEN
              FE(IE) = ( F(IP1)+F(IP2)+F(IP3)+F(IP4))/4.0E0
          ENDIF 
C
C         TOTAL SOURCE S = SIGN(DP)*G*F*(1+ANUC-F) FOR PRES3E
C         (SAME FORMULA AS [2], EVALUATED WITH UPDATED F)
          SOURCE = 0.0E0
          DELTAP = P(IE)+0.5E0*SIGMA
          IF(ICAVI.EQ.4) THEN
              FLQ  = MAX(0.0E0, MIN(1.0E0, FE(IE)))
              VNUC = (PI/6.0E0)*ANB0*DIAM0*DIAM0*DIAM0
              ANUC = VNUC/(1.0E0+VNUC)
              ARGR = (4.0E0*PI*ANB0*FLQ)
     *             / (3.0E0*(1.0E0+ANUC-FLQ))
              RR   = ARGR**(0.33333333E0)
              GSS  = 3.0E0*RR*SQRT(2.0E0*ABS(DELTAP)*COEF1/3.0E0)
              SOURCE = SIGN(1.0E0,DELTAP)*GSS*FLQ*(1.0E0+ANUC-FLQ)
          ENDIF
C
          IF(DELTAP.GE.0.0E0 .AND. FE(IE).GE.FLIMIT) SOURCE=0.0E0
C
          FESRC(IE) = SOURCE
C
 4300 CONTINUE
      IF(NEFREE.GE.1)THEN
      DO IBE = 1 , NEFREE
          IE = LEFREE(1,IBE)
          IF(IE .GE. 1) FESRC(IE) = 0.0E0
      END DO
      ENDIF
      CALL USTEND(133)
C
      RETURN
      END
