      SUBROUTINE LIQUID_MP(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *                  N1,N2,ME,MELM,NE,NP,MEP,NEX,NODE,NUMIP,
     *                  ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *                  NCRS,NCRS2,NPPMAX,ITPCRS,IUT6,IUT0,JUNROL,
     *                  NPP,IPCRS,LTAB,
     *					LSTCLR,LSTDGN,NITR,RES,
     *                  EPS,EPSRE,SIGMA,CGAS,CLQD,F00,
     *                  FLINLT,FLIMIT,FLMIN,
     *                  EAP1,EAP2,EAP3,EBP,AP,
     *                  MP,IENP,JENP,NODP,NEP,
     *                  SN,DELTA,CM,
     *                  UG,VG,WG,P,DT,VISC,SCT,F,FESRC,FE,
     *                  F3,F2,F1,F0,
     *                  CAVRAD,CJS,CALPHA,
     *                  NPINLT,LPINLT,
     *                  F3INLT,F2INLT,F1INLT,F0INLT,
     *                  NPSET, LPSET1,LPSET2,LPSET3,
     *                  COVER1,COVER2,COVER3,
     *                  NPSND ,LPSND ,NPTSND,IPSET,IPSRC,
     *                  NPRCV ,LPRCV ,NPTRCV,
     *                  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *                  LFIX3D,LPFIX,RX,RY,MWRK,WRKN,ACRS,A,
     *                  WRK01,WRK02,WRK03,WRK04,WRK05,
     *                  WRK06,WRK07,WRK08,
     *					RHS,A0,TS,
     *					F3X,F3Y,F2X,F2Y,
     *					F1X,F1Y,F0X,F0Y,
     *                  TACRS,IERR,ICRS_T,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  NEFREE,LEFREE,LP2NPFREE,
     *                  CEB,CUPWD,CVISL)

      IMPLICIT NONE
C
C[INPUT]
      INTEGER*4 MCOLOR,MCPART,
     *          NCOLOR(4),NCPART(MCOLOR,4),LLOOP(MCPART,MCOLOR,4),
     *          N,N1,N2,ME,MELM,NE,NP,MEP,NEX(12),NODE(N2,NE),NUMIP(NP),
     *          ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *          NCRS,NPP(NP),IPCRS(NCRS),LTAB(N1,N2,NE),
     *          LSTCLR(NCRS),LSTDGN(NCRS),
     *          NPPMAX,NCRS2,ITPCRS(NCRS2),IUT6,IUT0,JUNROL
      REAL*4    EPS,EPSRE,SIGMA,CGAS,CLQD,F00,FLINLT,FLIMIT,FLMIN,
     *          CAVRAD,CJS,CALPHA,F3INLT,F2INLT,F1INLT,F0INLT
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
      REAL*4    F(NP),F3(NP),F2(NP),F1(NP),F0(NP)
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
     *          RHS(NP),A0(NP),TS(0:NP),TACRS(NCRS2),
     *          F3Y(NP),F2Y(NP),F1Y(NP),F0Y(NP),
     *          F3X(NP),F2X(NP),F1X(NP),F0X(NP),PN(NP)
C
C[LOCAL]
      INTEGER*4 MAXBUF,IDUM,I,IE,IP,IB,ISEND,IDIM,NB,NPFIX,IBE
      INTEGER*4 IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8
      REAL*4    COEF,CEB,RESR,SOURCE,DELTAP,XE,
     *          CALPHA0,CALPHA1,FY,CAV3,AVERAD,EVPCND,FVX,
     *          PIRQ3,F3LIMIT,F3MAX,MAXH0
C
      INTEGER*4 IBCGS,NPMIN
      DATA IBCGS  / 0 /
C
      CHARACTER*60 ERMSGC
     & / ' ## SUBROUTINE LIQUID: FATAL      ERROR REPORT   ; RETURNED' /
C 
      INTEGER*4 ICRS,K
      INTEGER*4 ICRS_T(NP)
C
C
      REAL*4 CUPWD
      REAL*4 CVISL
C
      INTEGER*4 FLAG_CAV_XPOS
      REAL*4    CAV_XPOS,XPOS(NE)
      INTEGER*4 NEFREE
      REAL*4    LEFREE(2,NEFREE),LP2NPFREE(NP)
C
      IERR=0
      MAXBUF = ME*(N2+1)
C

C
C    FROM PER LIQUID VOLUME TO PER CELL VOLUME
          DO IP = 1 , NP
              F3(IP) = F(IP)*F3(IP)
              F2(IP) = F(IP)*F2(IP)
              F1(IP) = F(IP)*F1(IP)
              F0(IP) = F(IP)*F0(IP)
          END DO
C 
      IF(ITIME.EQ.0) GOTO 4000
      CALL USTSTA(124) 
    
C
CCYY [0] PREPARE FIX B.C. LIST
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
      CALL USTEND(124)
C
C
CCYY [1] DEFINE PARAMETERS
      PIRQ3   = 4.1887902E0
      CALPHA0 = CJS*CAVRAD/3.E0
      CALPHA1 = CALPHA0*CAVRAD
      MAXH0 = 0
      FLIMIT = 0.95
C
      CALL USTSTA(127)
C
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
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4)
     *                  +F3(IP5)+F3(IP6)+F3(IP7)+F3(IP8))/8.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4)
     *                  +F2(IP5)+F2(IP6)+F2(IP7)+F2(IP8))/8.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4)
     *                  +F1(IP5)+F1(IP6)+F1(IP7)+F1(IP8))/8.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4)
     *                  +F0(IP5)+F0(IP6)+F0(IP7)+F0(IP8))/8.0E0
		ELSE IF(IP6.NE.0) THEN
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4)
     *                  +F3(IP5)+F3(IP6))/6.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4)
     *                  +F2(IP5)+F2(IP6))/6.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4)
     *                  +F1(IP5)+F1(IP6))/6.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4)
     *                  +F0(IP5)+F0(IP6))/6.0E0
		ELSE IF(IP5.NE.0) THEN
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4)
     *                  +F3(IP5))/5.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4)
     *                  +F2(IP5))/5.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4)
     *                  +F1(IP5))/5.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4)
     *                  +F0(IP5))/5.0E0
		ELSE IF(IP4.NE.0) THEN
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4))/4.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4))/4.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4))/4.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4))/4.0E0
		ENDIF
C
C
C
		FY = 1.E0-PIRQ3*F3Y(IE)
C
C
        DELTAP=P(IE)+0.5E0*SIGMA
C
	IF(DELTAP.LT.0) THEN
		CAV3 = (CGAS-(CGAS-CLQD)*FY)*DELTAP
	ELSE
		CAV3 = (0.1*CGAS-(0.1*CGAS-0.1*CLQD)*FY)*DELTAP
	ENDIF
		CAV3 = CAV3/((CALPHA+FY-FY*FY)*3.E0)
C
C
		AVERAD = (CALPHA1*FY+F1Y(IE))/(CALPHA0*FY+F0Y(IE))
C
C
C
		SOURCE = -CJS*CAVRAD*CAV3 
                  EVPCND = -AVERAD*CAV3
C
C
C    PER UNIT VOLUME OF LIQUID
C
		F3X(IE)=F3INLT*0.3
		F2X(IE)=0.2268E0
C
	IF(F3Y(IE).LE.F3X(IE)) THEN
C
		F0X(IE) = FY*SOURCE
		 SOURCE  = SOURCE*CAVRAD
		F1X(IE) = FY*SOURCE+EVPCND*F0Y(IE)
		 SOURCE  = SOURCE*CAVRAD
		F2X(IE) = FY*SOURCE+2.E0*EVPCND*F1Y(IE)
		 SOURCE  = SOURCE*CAVRAD
		F3X(IE) = FY*SOURCE+3.E0*EVPCND*F2Y(IE)
		IF(F3X(IE).LE.0) F3X(IE) = 0
          FVX = PIRQ3*F3X(IE)
          F3X(IE) = F3X(IE)*FY
          F2X(IE) = F2X(IE)-F2Y(IE)*FVX
          F1X(IE) = F1X(IE)-F1Y(IE)*FVX
          F0X(IE) = F0X(IE)-F0Y(IE)*FVX
C
	ELSEIF(F3Y(IE).GE.F2X(IE)) THEN
C
		F0X(IE) = FY*SOURCE
		 SOURCE  = SOURCE*CAVRAD
		F1X(IE) = FY*SOURCE+EVPCND*F0Y(IE)
		 SOURCE  = SOURCE*CAVRAD
		F2X(IE) = FY*SOURCE+2.E0*EVPCND*F1Y(IE)
		 SOURCE  = SOURCE*CAVRAD
		F3X(IE) = FY*SOURCE+3.E0*EVPCND*F2Y(IE)
		IF(F3X(IE).GE.0) F3X(IE) = 0
          FVX = PIRQ3*F3X(IE)
          F3X(IE) = F3X(IE)*FY
          F2X(IE) = F2X(IE)-F2Y(IE)*FVX
          F1X(IE) = F1X(IE)-F1Y(IE)*FVX
          F0X(IE) = F0X(IE)-F0Y(IE)*FVX
C
	ELSE
C
		F0X(IE) = FY*SOURCE
		 SOURCE  = SOURCE*CAVRAD
		F1X(IE) = FY*SOURCE+EVPCND*F0Y(IE)
		 SOURCE  = SOURCE*CAVRAD
		F2X(IE) = FY*SOURCE+2.E0*EVPCND*F1Y(IE)
		 SOURCE  = SOURCE*CAVRAD
		F3X(IE) = FY*SOURCE+3.E0*EVPCND*F2Y(IE)
          FVX = PIRQ3*F3X(IE)
          F3X(IE) = F3X(IE)*FY
          F2X(IE) = F2X(IE)-F2Y(IE)*FVX
          F1X(IE) = F1X(IE)-F1Y(IE)*FVX
	  F0X(IE) = F0X(IE)-F0Y(IE)*FVX
		ENDIF
C
C          VISC(IE) = VISC0(IE)*RHO0/RHOBUF
C
 2000 CONTINUE
C
      IF(NEFREE.GE.1)THEN
      DO IBE = 1 , NEFREE
          IE = LEFREE(1,IBE)
          IF(IE .GE. 1)THEN
          F3X(IE) = 0.0E0
          F2X(IE) = 0.0E0
          F1X(IE) = 0.0E0
          F0X(IE) = 0.0E0
          ENDIF
      END DO
      ENDIF
      CALL USTEND(127)
C
            CALL CALC_F(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *                  N1,N2,ME,MELM,NE,NP,MEP,NEX,NODE,NUMIP,
     *                  ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *                  NCRS,NCRS2,NPPMAX,ITPCRS,IUT6,IUT0,JUNROL,
     *                  NPP,IPCRS,LTAB,
     *					LSTCLR,LSTDGN,NITR,RES,
     *                  EPS,EPSRE,SIGMA,CGAS,CLQD,F00,
     *                  FLINLT,FLIMIT,FLMIN,
     *                  EAP1,EAP2,EAP3,EBP,AP,
     *                  MP,IENP,JENP,NODP,NEP,
     *                  SN,DELTA,CM,
     *                  UG,VG,WG,P,DT,VISC,SCT,F,FESRC,FE,
     *                  F3,CAVRAD,CJS,CALPHA,
     *                  NPINLT,LPINLT,F3INLT,
     *                  NPSET, LPSET1,LPSET2,LPSET3,
     *                  COVER1,COVER2,COVER3,
     *                  NPSND ,LPSND ,NPTSND,IPSET,IPSRC,
     *                  NPRCV ,LPRCV ,NPTRCV,
     *                  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *                  LFIX3D,LPFIX,RX,RY,MWRK,WRKN,ACRS,A,
     *                  WRK01,WRK02,WRK03,WRK04,WRK05,
     *                  WRK06,WRK07,WRK08,
     *					RHS,A0,TS,F3X,
     *                  TACRS,IERR,ICRS_T,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  NEFREE,LEFREE,LP2NPFREE,
     *                  CEB,CUPWD,CVISL)
C
            CALL CALC_F(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *                  N1,N2,ME,MELM,NE,NP,MEP,NEX,NODE,NUMIP,
     *                  ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *                  NCRS,NCRS2,NPPMAX,ITPCRS,IUT6,IUT0,JUNROL,
     *                  NPP,IPCRS,LTAB,
     *					LSTCLR,LSTDGN,NITR,RES,
     *                  EPS,EPSRE,SIGMA,CGAS,CLQD,F00,
     *                  FLINLT,FLIMIT,FLMIN,
     *                  EAP1,EAP2,EAP3,EBP,AP,
     *                  MP,IENP,JENP,NODP,NEP,
     *                  SN,DELTA,CM,
     *                  UG,VG,WG,P,DT,VISC,SCT,F,FESRC,FE,
     *                  F2,CAVRAD,CJS,CALPHA,
     *                  NPINLT,LPINLT,F2INLT,
     *                  NPSET, LPSET1,LPSET2,LPSET3,
     *                  COVER1,COVER2,COVER3,
     *                  NPSND ,LPSND ,NPTSND,IPSET,IPSRC,
     *                  NPRCV ,LPRCV ,NPTRCV,
     *                  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *                  LFIX3D,LPFIX,RX,RY,MWRK,WRKN,ACRS,A,
     *                  WRK01,WRK02,WRK03,WRK04,WRK05,
     *                  WRK06,WRK07,WRK08,
     *					RHS,A0,TS,F2X,
     *                  TACRS,IERR,ICRS_T,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  NEFREE,LEFREE,LP2NPFREE,
     *                  CEB,CUPWD,CVISL)
C
            CALL CALC_F(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *                  N1,N2,ME,MELM,NE,NP,MEP,NEX,NODE,NUMIP,
     *                  ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *                  NCRS,NCRS2,NPPMAX,ITPCRS,IUT6,IUT0,JUNROL,
     *                  NPP,IPCRS,LTAB,
     *					LSTCLR,LSTDGN,NITR,RES,
     *                  EPS,EPSRE,SIGMA,CGAS,CLQD,F00,
     *                  FLINLT,FLIMIT,FLMIN,
     *                  EAP1,EAP2,EAP3,EBP,AP,
     *                  MP,IENP,JENP,NODP,NEP,
     *                  SN,DELTA,CM,
     *                  UG,VG,WG,P,DT,VISC,SCT,F,FESRC,FE,
     *                  F1,CAVRAD,CJS,CALPHA,
     *                  NPINLT,LPINLT,F1INLT,
     *                  NPSET, LPSET1,LPSET2,LPSET3,
     *                  COVER1,COVER2,COVER3,
     *                  NPSND ,LPSND ,NPTSND,IPSET,IPSRC,
     *                  NPRCV ,LPRCV ,NPTRCV,
     *                  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *                  LFIX3D,LPFIX,RX,RY,MWRK,WRKN,ACRS,A,
     *                  WRK01,WRK02,WRK03,WRK04,WRK05,
     *                  WRK06,WRK07,WRK08,
     *					RHS,A0,TS,F1X,
     *                  TACRS,IERR,ICRS_T,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  NEFREE,LEFREE,LP2NPFREE,
     *                  CEB,CUPWD,CVISL)
C
            CALL CALC_F(MCOLOR,MCPART,NCOLOR,NCPART,LLOOP,
     *                  N1,N2,ME,MELM,NE,NP,MEP,NEX,NODE,NUMIP,
     *                  ICAVI,ISOLT,NS,NL,ITIME,JSET,NMAX,
     *                  NCRS,NCRS2,NPPMAX,ITPCRS,IUT6,IUT0,JUNROL,
     *                  NPP,IPCRS,LTAB,
     *					LSTCLR,LSTDGN,NITR,RES,
     *                  EPS,EPSRE,SIGMA,CGAS,CLQD,F00,
     *                  FLINLT,FLIMIT,FLMIN,
     *                  EAP1,EAP2,EAP3,EBP,AP,
     *                  MP,IENP,JENP,NODP,NEP,
     *                  SN,DELTA,CM,
     *                  UG,VG,WG,P,DT,VISC,SCT,F,FESRC,FE,
     *                  F0,CAVRAD,CJS,CALPHA,
     *                  NPINLT,LPINLT,F0INLT,
     *                  NPSET, LPSET1,LPSET2,LPSET3,
     *                  COVER1,COVER2,COVER3,
     *                  NPSND ,LPSND ,NPTSND,IPSET,IPSRC,
     *                  NPRCV ,LPRCV ,NPTRCV,
     *                  IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *                  LFIX3D,LPFIX,RX,RY,MWRK,WRKN,ACRS,A,
     *                  WRK01,WRK02,WRK03,WRK04,WRK05,
     *                  WRK06,WRK07,WRK08,
     *					RHS,A0,TS,F0X,
     *                  TACRS,IERR,ICRS_T,
     *                  FLAG_CAV_XPOS,CAV_XPOS,XPOS,
     *                  NEFREE,LEFREE,LP2NPFREE,
     *                  CEB,CUPWD,CVISL)
C
 4000 CONTINUE
C
      CALL USTSTA(133)
C
      DO 4100 IP = 1 , NP
          F3(IP)  = AMAX1(0.0E0,F3(IP))
          F2(IP)  = AMAX1(0.0E0,F2(IP))
          F1(IP)  = AMAX1(0.0E0,F1(IP))
          F0(IP)  = AMAX1(0.0E0,F0(IP))               
 4100 CONTINUE
C
      IF (JSET.GE.1) THEN
          CALL OVRST1(IPART,NPSET,N1,ME,NE,NP,NEX,NODE,F3,
     *                LPSET1,LPSET2,LPSET3,COVER1,COVER2,COVER3,
     *                NDOM,MBPDOM,NPSND,NPRCV,
     *                LPSND,NPTSND,LPRCV,NPTRCV,IPSET,IPSRC,
     *                WRK01,WRK02,RX,RY,IUT0,IERR)
          IF(IERR.NE.0) THEN
              WRITE(IUT0,*)'ERROR CODE REPORTED FROM OVRST1'
              WRITE(IUT0,*)ERMSGC
              RETURN
          ENDIF
          CALL OVRST1(IPART,NPSET,N1,ME,NE,NP,NEX,NODE,F2,
     *                LPSET1,LPSET2,LPSET3,COVER1,COVER2,COVER3,
     *                NDOM,MBPDOM,NPSND,NPRCV,
     *                LPSND,NPTSND,LPRCV,NPTRCV,IPSET,IPSRC,
     *                WRK01,WRK02,RX,RY,IUT0,IERR)
          IF(IERR.NE.0) THEN
              WRITE(IUT0,*)'ERROR CODE REPORTED FROM OVRST1'
              WRITE(IUT0,*)ERMSGC
              RETURN
          ENDIF
          CALL OVRST1(IPART,NPSET,N1,ME,NE,NP,NEX,NODE,F1,
     *                LPSET1,LPSET2,LPSET3,COVER1,COVER2,COVER3,
     *                NDOM,MBPDOM,NPSND,NPRCV,
     *                LPSND,NPTSND,LPRCV,NPTRCV,IPSET,IPSRC,
     *                WRK01,WRK02,RX,RY,IUT0,IERR)
          IF(IERR.NE.0) THEN
              WRITE(IUT0,*)'ERROR CODE REPORTED FROM OVRST1'
              WRITE(IUT0,*)ERMSGC
              RETURN
          ENDIF
          CALL OVRST1(IPART,NPSET,N1,ME,NE,NP,NEX,NODE,F0,
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
          F3(LPINLT(IB)) = F3INLT
          F2(LPINLT(IB)) = F2INLT
          F1(LPINLT(IB)) = F1INLT
          F0(LPINLT(IB)) = F0INLT
 4200 CONTINUE   
C
C     <<<<< INTERPOLATE SIGMA FROM ELEMENT TO NODE>>>>>
C
C
      PIRQ3   = 4.1887902E0
      F3LIMIT = (1.E0-0.99)/PIRQ3
      F3MAX   = (1.E0-FLMIN)/PIRQ3
C
      CALL NODLEX(NODE,ME,NP,NE,NP,N1,N2,NEX,SN,
     *            IPART,LDOM,NBPDOM,NDOM,IPSLF,IPSND,MBPDOM,
     *            P,PN,CM,IUT0,IERR,RX,RY,MAXBUF)
C
		DO 4300 IP = 1 , NP
		IF(F3(IP).LE.F3LIMIT) THEN
C			IF(PN(IP).GE.-0.5E0*SIGMA) THEN
C                  F3(IP) = AMAX1(0.E0 ,F3(IP))
C                  F2(IP) = 0.E0
C                  F1(IP) = 0.E0
C                  F0(IP) = 0.E0
C			ELSE
                  F3(IP) = AMAX1(0.E0 ,F3(IP))
                  F2(IP) = AMAX1(0.E0 ,F2(IP))
                  F1(IP) = AMAX1(0.E0 ,F1(IP))
                  F0(IP) = AMAX1(0.E0 ,F0(IP))
C			ENDIF
		ELSE
              F3(IP) = AMIN1(F3MAX,F3(IP))
              F2(IP) = AMAX1(0.E0 ,F2(IP))
              F1(IP) = AMAX1(0.E0 ,F1(IP))
              F0(IP) = AMAX1(0.E0 ,F0(IP))
		ENDIF
		IF(F0(IP).LE.0) THEN
            F3(IP) = 0
            F2(IP) = 0
            F1(IP) = 0
        ELSEIF(F1(IP).LE.0) THEN
            F3(IP) = 0
            F2(IP) = 0              		
        ELSEIF(F2(IP).LE.0) THEN
            F3(IP) = 0
 		ENDIF
C
 4300 CONTINUE
C
C
C
C
C
C
C	  MAXH0 = F(NPMIN)
C       WRITE(IUT6,*) ' FL AT ELEMENT P-MIN =',MAXH0
C
      PIRQ3   = 4.1887902E0
      CALPHA0 = CJS*CAVRAD/3.E0
      CALPHA1 = CALPHA0*CAVRAD
      MAXH0 = 0
C
      DO 4400 IE = 1 , NE
          IP1=NODE(1,IE)
          IP2=NODE(2,IE)
          IP3=NODE(3,IE)
          IP4=NODE(4,IE)
          IP5=NODE(5,IE)
          IP6=NODE(6,IE)
          IP7=NODE(7,IE)
          IP8=NODE(8,IE) 
		IF(IP8.NE.0) THEN 
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4)
     *                  +F3(IP5)+F3(IP6)+F3(IP7)+F3(IP8))/8.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4)
     *                  +F2(IP5)+F2(IP6)+F2(IP7)+F2(IP8))/8.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4)
     *                  +F1(IP5)+F1(IP6)+F1(IP7)+F1(IP8))/8.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4)
     *                  +F0(IP5)+F0(IP6)+F0(IP7)+F0(IP8))/8.0E0
		ELSE IF(IP6.NE.0) THEN
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4)
     *                  +F3(IP5)+F3(IP6))/6.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4)
     *                  +F2(IP5)+F2(IP6))/6.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4)
     *                  +F1(IP5)+F1(IP6))/6.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4)
     *                  +F0(IP5)+F0(IP6))/6.0E0
		ELSE IF(IP5.NE.0) THEN
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4)
     *                  +F3(IP5))/5.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4)
     *                  +F2(IP5))/5.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4)
     *                  +F1(IP5))/5.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4)
     *                  +F0(IP5))/5.0E0
		ELSE IF(IP4.NE.0) THEN
              F3Y(IE) = ( F3(IP1)+F3(IP2)+F3(IP3)+F3(IP4))/4.0E0
              F2Y(IE) = ( F2(IP1)+F2(IP2)+F2(IP3)+F2(IP4))/4.0E0
              F1Y(IE) = ( F1(IP1)+F1(IP2)+F1(IP3)+F1(IP4))/4.0E0
              F0Y(IE) = ( F0(IP1)+F0(IP2)+F0(IP3)+F0(IP4))/4.0E0
		ENDIF
C
		FE(IE) = 1.E0-PIRQ3*F3Y(IE)
C
        DELTAP=P(IE)+0.5E0*SIGMA		
C	
		IF(DELTAP.LE.0) THEN
		CAV3 =  (CGAS-(CGAS-CLQD)*FE(IE))
             	ELSE
		CAV3 =  (1.0E0-FE(IE))*(0.1*CGAS-(0.1*CGAS-0.1*CLQD)*FE(IE))
             	ENDIF
C
		CAV3 = CAV3*DELTAP
		CAV3 = CAV3/((CALPHA+FE(IE)-FE(IE)*FE(IE))*3.E0)
C
C
		AVERAD = (CALPHA1*FE(IE)+F1Y(IE))/(CALPHA0*FE(IE)+F0Y(IE))
C
C
C
		SOURCE = -CJS*CAVRAD*CAV3 

		EVPCND = -AVERAD*CAV3
C
		SOURCE  = SOURCE*CAVRAD*CAVRAD*CAVRAD
C
C    PER UNIT VOLUME OF LIQUID
C
		F3X(IE)=F3INLT*0.3
		F2X(IE)=0.2268E0
C
	IF(F3Y(IE).LE.F3X(IE)) THEN
		F3X(IE) = FE(IE)*SOURCE+3.E0*EVPCND*F2Y(IE)
		IF(F3X(IE).LE.0) F3X(IE) = 0
	ELSEIF(F3Y(IE).GE.F2X(IE)) THEN
		F3X(IE) = FE(IE)*SOURCE+3.E0*EVPCND*F2Y(IE)
		IF(F3X(IE).GE.0) F3X(IE) = 0
	ELSE
		F3X(IE) = FE(IE)*SOURCE+3.E0*EVPCND*F2Y(IE)
	ENDIF
C
C
C
        FESRC(IE) = -PIRQ3*F3X(IE)
C
C
        IF(F3X(IE).LE.0.0E0 .AND. FE(IE).GE.FLIMIT)  FESRC(IE) = 0.0E0
C
C
 4400 CONTINUE
C
C    FROM PER VOLUME TO PER LIQUID VOLUME
              DO 4500 IP = 1 , NP
                  F(IP) = 1.E0-PIRQ3*F3(IP)
C	
         		  F(IP)  = AMAX1(FLMIN,F(IP))
          	      F(IP)  = AMIN1(1.E0 ,F(IP))
C
                  F3(IP) = F3(IP)/F(IP)
                  F2(IP) = F2(IP)/F(IP)
                  F1(IP) = F1(IP)/F(IP)
                  F0(IP) = F0(IP)/F(IP)
 4500         CONTINUE
C
      CALL USTEND(133)
C
      RETURN
      END
