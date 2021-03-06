CSTART OF AS 100
      SUBROUTINE AJV(SNV, JVAL, ITYPE, GAMMA, DELTA, XLAM, XI, IFAULT)
C
C        ALGORITHM AS 100.1  APPL. STATIST. (1976) VOL.25, P.190
C
C        CONVERTS A STANDARD NORMAL VARIATE (SNV) TO A
C        JOHNSON VARIATE (JVAL)
C
      DOUBLE PRECISION SNV, JVAL, GAMMA, DELTA, XLAM, XI, V, W, ZERO, 
     $  HALF, ONE, ZDABS, ZDEXP, ZDSIGN
C
      DATA ZERO, HALF, ONE /0.0, 0.5, 1.0/
C
      ZDABS(W) = DABS(W)
      ZDEXP(W) = DEXP(W)
      ZDSIGN(W, V) = DSIGN(W, V)
C
      JVAL = ZERO
      IFAULT = 1
      IF (ITYPE .LT. 1 .OR. ITYPE .GT. 4) RETURN
      IFAULT = 0
      GOTO (10, 20, 30, 40), ITYPE
C
C        SL DISTRIBUTION
C
   10 JVAL = XLAM * ZDEXP((XLAM * SNV - GAMMA) / DELTA) + XI
      RETURN
C
C        SU DISTRIBUTION
C
   20 W = ZDEXP((SNV - GAMMA) / DELTA)
      W = HALF * (W - ONE / W)
      JVAL = XLAM * W + XI
      RETURN
C
C        SB DISTRIBUTION
C
   30 W = (SNV - GAMMA) / DELTA
      V = ZDEXP(-ZDABS(W))
      V = (ONE - V) / (ONE + V)
      JVAL = HALF * XLAM * (ZDSIGN(V, W) + ONE) + XI
      RETURN
C
C        NORMAL DISTRIBUTION
C
   40 JVAL = (SNV - GAMMA) / DELTA
      RETURN
      END
C
      SUBROUTINE SNV(AJV, NVAL, ITYPE, GAMMA, DELTA, XLAM, XI, IFAULT)
C
C        ALGORITHM AS 100.2  APPL. STATIST. (1976) VOL.25, P.190
C
C        CONVERTS A JOHNSON VARIATE (AJV) TO A
C        STANDARD NORMAL VARIATE (NVAL)
C
      DOUBLE PRECISION AJV, GAMMA, DELTA, XLAM, XI, V, W, C, ZERO, 
     $  HALF, ONE, ZDLOG, ZDSQRT
C
      DATA ZERO, HALF, ONE, C /0.0, 0.5, 1.0, -63.0/
C
      ZDLOG(W) = DLOG(W)
      ZDSQRT(W) = DSQRT(W)
C
      NVAL = ZERO
      IFAULT = 1
      IF (ITYPE .LT. 1 .OR. ITYPE .GT. 4) RETURN
      IFAULT = 0
      GOTO (10, 20, 30, 40), ITYPE
C
C        SL DISTRIBUTION
C
   10 W = XLAM * (AJV - XI)
      IF (W .LE. ZERO) GOTO 15
      NVAL = XLAM * (ZDLOG(W) * DELTA + GAMMA)
      RETURN
   15 IFAULT = 2
      RETURN
C
C        SU DISTRIBUTION
C
   20 W = (AJV - XI) / XLAM
      IF (W .GT. C) GOTO 23
      W = -HALF / W
      GOTO 27
   23 W = ZDSQRT(W * W + ONE) + W
   27 NVAL = ZDLOG(W) * DELTA + GAMMA
      RETURN
C
C        SB DISTRIBUTION
C
   30 W = AJV - XI
      V = XLAM - W
      IF (W .LE. ZERO .OR. V .LE. ZERO) GOTO 35
      NVAL = ZDLOG(W / V) * DELTA + GAMMA
      RETURN
   35 IFAULT = 2
      RETURN
C
C        NORMAL DISTRIBUTION
C
   40 NVAL = DELTA * AJV + GAMMA
      RETURN
      END
CEND OF AS 100

