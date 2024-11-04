C
C	Modified by GPN on 4th Nov 2024 to remove Computed GOTO
C
CSTART OF AS 100
      SUBROUTINE AJV(SNV, JVAL, ITYPE, GAMMA, DELTA, XLAM, XI, IFAULT)
C
C        ALGORITHM AS 100.1  APPL. STATIST. (1976) VOL.25, P.190
C
C        CONVERTS A STANDARD NORMAL VARIATE (SNV) TO A
C        JOHNSON VARIATE (JVAL)
C
      DOUBLE PRECISION SNV, JVAL, GAMMA, DELTA, XLAM, XI, V, W, ZERO, 
     $  HALF, ONE
C
      DATA ZERO, HALF, ONE /0.0, 0.5, 1.0/
C
C	Replaced following statement functions ZDABS, ZDEXP, ZDSIGN
C     ZDABS(W) = DABS(W)
C     ZDEXP(W) = DEXP(W)
C     ZDSIGN(W, V) = DSIGN(W, V)
C
      JVAL = ZERO
      IFAULT = 1
      IF (ITYPE .LT. 1 .OR. ITYPE .GT. 4) RETURN
      IFAULT = 0

C
C	Replaced Computed GOTO With IF - THEN - ELSEIF - END IF block

      IF (ITYPE .EQ. 1) THEN

C		SL DISTRIBUTION (was GOTO 10)
         JVAL = XLAM * DEXP((XLAM * SNV - GAMMA) / DELTA) + XI
         RETURN

      ELSE IF (ITYPE .EQ. 2) THEN

C               SU DISTRIBUTION (was GOTO 20)
         W = DEXP((SNV - GAMMA) / DELTA)
         W = HALF * (W - ONE / W)
         JVAL = XLAM * W + XI
         RETURN

      ELSE IF (ITYPE .EQ. 3) THEN

C              SB DISTRIBUTION (was GOTO 30)
         W = (SNV - GAMMA) / DELTA
         V = DEXP(-DABS(W))
         V = (ONE - V) / (ONE + V)
         JVAL = HALF * XLAM * (DSIGN(V, W) + ONE) + XI
         RETURN

      ELSE

C        NORMAL DISTRIBUTION (was GOTO 40)

         JVAL = (SNV - GAMMA) / DELTA
         RETURN
      END IF
      END
C
C     waveband does not require subroutine SNV so it has been commented out
C
C     SUBROUTINE SNV(AJV, NVAL, ITYPE, GAMMA, DELTA, XLAM, XI, IFAULT)
C
C        ALGORITHM AS 100.2  APPL. STATIST. (1976) VOL.25, P.190
C
C        CONVERTS A JOHNSON VARIATE (AJV) TO A
C        STANDARD NORMAL VARIATE (NVAL)
C
C     DOUBLE PRECISION AJV, GAMMA, DELTA, XLAM, XI, V, W, C, ZERO, 
C    $  HALF, ONE
C
C     DATA ZERO, HALF, ONE, C /0.0, 0.5, 1.0, -63.0/
C
C	Don't need following statement functions ZDLOG, ZDSQRT
C     ZDLOG(W) = DLOG(W)
C     ZDSQRT(W) = DSQRT(W)
C
C     NVAL = ZERO
C     IFAULT = 1
C     IF (ITYPE .LT. 1 .OR. ITYPE .GT. 4) RETURN
C     IFAULT = 0
C     Removing Computed GOTO (10, 20, 30, 40), ITYPE
C
C     IF (ITYPE .EQ. 1) THEN
C
C        SL DISTRIBUTION (was GOTO 10)
C
C        W = XLAM * (AJV - XI)
C
C        IF (W .LE. ZERO) THEN
C           IFAULT = 2
C        ELSE
C           NVAL = XLAM * (DLOG(W) * DELTA + GAMMA)
C        END IF
C        RETURN
C
C     ELSE IF (ITYPE .EQ. 2) THEN
C
C        SU DISTRIBUTION (was GOTO 20)
C
C        W = (AJV - XI) / XLAM
C        IF (W .GT. C) THEN
C           W = DSQRT(W * W + ONE) + W
C        ELSE
C           W = -HALF / W
C        END IF
C        NVAL = DLOG(W) * DELTA + GAMMA
C        RETURN
C
C     ELSE IF (ITYPE .EQ. 3) THEN
C 
C        SB DISTRIBUTION (was GOTO 30)
C
C        W = AJV - XI
C        V = XLAM - W
C        IF (W .LE. ZERO .OR. V .LE. ZERO) THEN
C           IFAULT = 2
C        ELSE
C           NVAL = DLOG(W / V) * DELTA + GAMMA
C        END IF
C        RETURN
C
C      ELSE
C 
C        NORMAL DISTRIBUTION (was GOTO 40)
C
C        NVAL = DELTA * AJV + GAMMA
C        RETURN
C
C     END IF
C     END
CEND OF AS 100
