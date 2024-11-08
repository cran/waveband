CSTART OF AS 99
      SUBROUTINE JNSN(XBAR, SD, RB1, BB2, ITYPE, GAMMA, DELTA,
     $  XLAM, XI, IFAULT)
C
C        ALGORITHM AS 99  APPL. STATIST. (1976) VOL.25, P.180
C
C        FINDS TYPE AND PARAMETERS OF A JOHNSON CURVE
C        WITH GIVEN FIRST FOUR MOMENTS
C
      DOUBLE PRECISION XBAR, SD, RB1, BB2, GAMMA, DELTA, XLAM, XI, TOL,
     $  B1, B2, Y, X, U, W, ZERO, ONE, TWO, THREE, FOUR, HALF,
     $  QUART
      LOGICAL FAULT
C
      DATA TOL /0.01/
      DATA ZERO, QUART, HALF, ONE, TWO, THREE, FOUR
     $     /0.0,  0.25,  0.5, 1.0, 2.0,   3.0,  4.0/
C

C     Don't need statement fns ZDABS, ZDEXP, ZLOG, ZDSIGN, ZDSQRT
C     ZDABS(X) = DABS(X)
C     ZDEXP(X) = DEXP(X)
C     ZLOG(X) = DLOG(X)
C     ZDSIGN(X, Y) = DSIGN(X, Y)
C     ZDSQRT(X) = DSQRT(X)
C
      IFAULT = 1
      IF (SD .LT. ZERO) RETURN
      IFAULT = 0
      XI = ZERO
      XLAM = ZERO
      GAMMA = ZERO
      DELTA = ZERO
      IF (SD .GT. ZERO) GOTO 10
      ITYPE = 5
      XI = XBAR
      RETURN
   10 B1 = RB1 * RB1
      B2 = BB2
      FAULT = .FALSE.
C
C        TEST WHETHER LOGNORMAL (OR NORMAL) REQUESTED
C
      IF (B2 .GE. ZERO) GOTO 30
   20 IF (DABS(RB1) .LE. TOL) GOTO 70
      GOTO 80
C
C        TEST FOR POSITION RELATIVE TO BOUNDARY LINE
C
   30 IF (B2 .GT. B1 + TOL + ONE) GOTO 60
      IF (B2 .LT. B1 + ONE) GOTO 50
C
C        ST DISTRIBUTION
C
   40 ITYPE = 5
      Y = HALF + HALF * DSQRT(ONE - FOUR / (B1 + FOUR))
      IF (RB1 .GT. ZERO) Y = ONE - Y
      X = SD / DSQRT(Y * (ONE - Y))
      XI = XBAR - Y * X
      XLAM = XI + X
      DELTA = Y
      RETURN
   50 IFAULT = 2
      RETURN
   60 IF (DABS(RB1) .GT. TOL .OR. DABS(B2 - THREE) .GT. TOL) GOTO 80
C
C        NORMAL DISTRIBUTION
C
   70 ITYPE = 4
      DELTA = ONE / SD
      GAMMA = -XBAR / SD
      RETURN
C
C        TEST FOR POSITION RELATIVE TO LOGNORMAL LINE
C
   80 X = HALF * B1 + ONE
      Y = DABS(RB1) * DSQRT(QUART * B1 + ONE)
      U = (X + Y) ** (ONE / THREE)
      W = U + ONE / U - ONE
      U = W * W * (THREE + W * (TWO + W)) - THREE
      IF (B2 .LT. ZERO .OR. FAULT) B2 = U
      X = U - B2
      IF (DABS(X) .GT. TOL) GOTO 90
C
C        LOGNORMAL (SL) DISTRIBUTION
C
      ITYPE = 1
      XLAM = DSIGN(ONE, RB1)
      U = XLAM * XBAR
      X = ONE / DSQRT(DLOG(W))
      DELTA = X
      Y = HALF * X * DLOG(W * (W - ONE) / (SD * SD))
      GAMMA = Y
      XI = XLAM * (U - DEXP((HALF / X - Y) / X))
      RETURN
C
C        SB OR SU DISTRIBUTION
C
   90 IF (X .GT. ZERO) GOTO 100
      ITYPE = 2
      CALL SUFIT(XBAR, SD, RB1, B2, GAMMA, DELTA, XLAM, XI)
      RETURN
  100 ITYPE = 3
      CALL SBFIT(XBAR, SD, RB1, B2, GAMMA, DELTA, XLAM, XI, FAULT)
      IF (.NOT. FAULT) RETURN
C
C        FAILURE - TRY TO FIT APPROXIMATE RESULT
C
      IFAULT = 3
      IF (B2 .GT. B1 + TWO) GOTO 20
      GOTO 40
      END
C
      SUBROUTINE SUFIT (XBAR, SD, RB1, B2, GAMMA, DELTA, XLAM, XI)
C
C        ALGORITHM AS 99.1  APPL. STATIST. (1976) VOL.25, P.180
C
C        FINDS PARAMETERS OF JOHNSON SU CURVE WITH
C        GIVEN FIRST FOUR MOMENTS
C
      DOUBLE PRECISION XBAR, SD, RB1, B2, GAMMA, DELTA, XLAM, XI, 
     $  TOL, B1,
     $  B3, W, Y, W1, WM1, Z, V, A, B, X, ZERO, ONE, TWO, THREE,
     $  FOUR, SIX, SEVEN, EIGHT, NINE, TEN, HALF, ONE5, TWO8,
     $  SIXTEN
C
      DATA TOL /0.0001/
      DATA ZERO,  ONE,  TWO,  THREE, FOUR,  SIX, SEVEN,
     $    EIGHT, NINE,  TEN, SIXTEN, HALF, ONE5,  TWO8
     $     /0.0,  1.0,  2.0,    3.0,  4.0,  6.0,   7.0,
     $      8.0,  9.0, 10.0,   16.0,  0.5,  1.5,   2.8/
C
C     Remove statement fns ZDABS, ZDEXP, ZLOG, ZDSIGN, ZDSQRT
C     ZDABS(X) = DABS(X)
C     ZDEXP(X) = DEXP(X)
C     ZLOG(X) = DLOG(X)
C      ZDSIGN(X, Y) = DSIGN(X, Y)
C     ZDSQRT(X) = DSQRT(X)
C
      B1 = RB1 * RB1
      B3 = B2 - THREE
C
C        W IS FIRST ESTIMATE OF DEXP(DELTA ** (-2))
C
      W = DSQRT(TWO * B2 - TWO8 * B1 - TWO)
      W = DSQRT(W-ONE)
      IF (DABS(RB1) .GT. TOL) GOTO 10
C
C        SYMMETRICAL CASE - RESULTS ARE KNOWN
C
      Y = ZERO
      GOTO 20
C
C        JOHNSON ITERATION (USING Y FOR HIS M)
C
C   10 PRINT*,'Entering iteration; W is', W, ' and B3 is ', B3
   10 W1 = W + ONE
      WM1 = W - ONE
      Z = W1 * B3
C      PRINT*,'           Calculated Z is', Z
      V = W * (SIX + W * (THREE + W))
C      PRINT*,'           Calculated V is', V
      A = EIGHT * (WM1 * (THREE + W * (SEVEN + V)) - Z)
C      PRINT*,'           Calculated A is', A
      B = SIXTEN * (WM1 * (SIX + V) - B3)
C      PRINT*,'           Calculated B is', B
      Y = (DSQRT(A * A - TWO * B * (WM1 * (THREE + W *
     $  (NINE + W * (TEN + V))) - TWO * W1 * Z)) - A) / B
C      PRINT*,'           Calculated Y is', Y
      Z = Y * WM1 * (FOUR * (W + TWO) * Y + THREE * W1 * W1) ** 2 /
     $  (TWO * (TWO * Y + W1) ** 3)
C      PRINT*,'           Calculated Z is', Z
      V = W * W
C      PRINT*,'           Calculated V is', V
      W = DSQRT(ONE - TWO * (ONE5 - B2 + (B1 *
     $  (B2 - ONE5 - V * (ONE + HALF * V))) / Z))
C      PRINT*,'           Calculated W is', W
      W = DSQRT(W-ONE)
C      PRINT*,'           Calculated W is', W
C      PRINT*,'In iteration; Z = ', Z, ' and W = ', W
      IF (DABS(B1 - Z) .GT. TOL) GOTO 10
C
C        END OF ITERATION
C
C      PRINT*,'End of iteration; Z = ', Z, ' and W = ', W
C      PRINT*,'(B1 is', B1, ')'
      Y = Y / W
C      PRINT*,'Y = Y / W = ', Y
      Y = DLOG(DSQRT(Y) + DSQRT(Y + ONE))
C      PRINT*, 'Y = DLOG(DSQRT(Y) + DSQRT(Y + ONE)) = ', Y
      IF (RB1 .GT. ZERO) Y = -Y
C      PRINT*, 'IF (RB1 .GT. ZERO) Y = -Y gives Y = ', Y
   20 X = DSQRT(ONE / DLOG(W))
      DELTA = X
      GAMMA = Y * X
      Y = DEXP(Y)
      Z = Y * Y
      X = SD / DSQRT(HALF * (W - ONE) * (HALF * W *
     $  (Z + ONE / Z) + ONE))
      XLAM = X
      XI = (HALF * DSQRT(W) * (Y - ONE / Y)) * X + XBAR
      RETURN
      END
      SUBROUTINE SBFIT(XBAR, SIGMA, RTB1, B2, GAMMA, DELTA, XLAM,
     $  XI, FAULT)
C
C        ALGORITHM AS 99.2  APPL. STATIST. (1976) VOL.25, P.180
C
C        FINDS PARAMETERS OF JOHNSON SB CURVE WITH
C        GIVEN FIRST FOUR MOMENTS
C
      DOUBLE PRECISION HMU(6), DERIV(4), DD(4), XBAR, SIGMA, RTB1, B2, 
     $  GAMMA,
     $  DELTA, XLAM, XI, TT, TOL, RB1, B1, E, U, X, Y, W, F, D,
     $  G, S, H2, T, H2A, H2B, H3, H4, RBET, BET2, ZERO, ONE,
     $  TWO, THREE, FOUR, SIX, HALF, QUART, ONE5, A1, A2, A3,
     $  A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15,
     $  A16, A17, A18, A19, A20, A21, A22
      LOGICAL NEG, FAULT
C
      DATA TT, TOL, LIMIT /1.0E-4, 0.0001, 5000/
      DATA ZERO, ONE, TWO, THREE, FOUR, SIX, HALF, QUART, ONE5
     $     /0.0, 1.0, 2.0,   3.0,  4.0, 6.0,  0.5,  0.25,  1.5/
      DATA     A1,     A2,     A3,     A4,     A5,     A6,
     $         A7,     A8,     A9,    A10,    A11,    A12,
     $        A13,    A14,    A15,    A16,    A17,    A18,
     $        A19,    A20,    A21,    A22
     $    /0.0124, 0.0623, 0.4043,  0.408,  0.479,  0.485,
     $     0.5291, 0.5955,  0.626,   0.64, 0.7077, 0.7466,
     $        0.8, 0.9281, 1.0614,   1.25, 1.7973,    1.8,
     $      2.163,    2.5, 8.5245, 11.346/
C
C     Remove statement fns ZDABS, ZLOG, ZDSQRT
C     ZDABS(X) = DABS(X)
C     ZLOG(X) = DLOG(X)
C     ZDSQRT(X) = DSQRT(X)
C
      RB1 = DABS(RTB1)
      B1 = RB1 * RB1
      NEG = RTB1 .LT. ZERO
C
C        GET D AS FIRST ESTIMATE OF DELTA
C
      E = B1 + ONE
      X = HALF * B1 + ONE
      Y = DABS(RB1) * DSQRT(QUART * B1 + ONE)
      U = (X + Y) ** (ONE / THREE)
      W = U + ONE / U - ONE
      F = W * W * (THREE + W * (TWO + W)) - THREE
      E = (B2 - E) / (F - E)
      IF (DABS(RB1) .GT. TOL) GOTO 5
      F = TWO
      GOTO 20
    5 D = ONE / DSQRT(DLOG(W))
      IF (D .LT. A10) GOTO 10
      F = TWO - A21 / (D * (D * (D - A19) + A22))
      GOTO 20
   10 F = A16 * D
   20 F = E * F + ONE
      IF (F .LT. A18) GOTO 25
      D = (A9 * F - A4) * (THREE - F) ** (-A5)
      GOTO 30
   25 D = A13 * (F - ONE)
C
C        GET G AS FIRST ESTIMATE OF GAMMA
C
   30 G = ZERO
      IF (B1 .LT. TT) GOTO 70
      IF (D .GT. ONE) GOTO 40
      G = (A12 * D ** A17 + A8) * B1 ** A6
      GOTO 70
   40 IF (D .LE. A20) GOTO 50
      U = A1
      Y = A7
      GOTO 60
   50 U = A2
      Y = A3
   60 G = B1 ** (U * D + Y) * (A14 + D * (A15 * D - A11))
   70 M = 0
C
C        MAIN ITERATION STARTS HERE
C
   80 M = M + 1
      FAULT = M .GT. LIMIT
      IF (FAULT) RETURN
C
C        GET FIRST SIX MOMENTS FOR LATEST G AND D VALUES
C
      CALL MOM(G, D, HMU, FAULT)
      IF (FAULT) RETURN
      S = HMU(1) * HMU(1)
      H2 = HMU(2) - S
      FAULT = H2 .LE. ZERO
      IF (FAULT) RETURN
      T = DSQRT(H2)
      H2A = T * H2
      H2B = H2 * H2
      H3 = HMU(3) - HMU(1) * (THREE * HMU(2) - TWO * S)
      RBET = H3 / H2A
      H4 = HMU(4) - HMU(1) * (FOUR * HMU(3) - HMU(1) *
     $  (SIX * HMU(2) - THREE * S))
      BET2 = H4 / H2B
      W = G * D
      U = D * D
C
C        GET DERIVATIVES
C
      DO 120 J = 1, 2
      DO 110 K = 1, 4
      T = K
      IF (J .EQ. 1) GOTO 90
      S = ((W - T) * (HMU(K) - HMU(K + 1)) + (T + ONE) *
     $  (HMU(K + 1) - HMU(K + 2))) / U
      GOTO 100
   90 S = HMU(K + 1) - HMU(K)
  100 DD(K) = T * S / D
  110 CONTINUE
      T = TWO * HMU(1) * DD(1)
      S = HMU(1) * DD(2)
      Y = DD(2) - T
      DERIV(J) = (DD(3) - THREE * (S + HMU(2) * DD(1) - T * HMU(1))
     $  - ONE5 * H3 * Y / H2) / H2A
      DERIV(J + 2) = (DD(4) - FOUR * (DD(3) * HMU(1) + DD(1) * HMU(3))
     $  + SIX * (HMU(2) * T + HMU(1) * (S - T * HMU(1)))
     $  - TWO * H4 * Y / H2) / H2B
  120 CONTINUE
      T = ONE / (DERIV(1) * DERIV(4) - DERIV(2) * DERIV(3))
      U = (DERIV(4) * (RBET - RB1) - DERIV(2) * (BET2 - B2)) * T
      Y = (DERIV(1) * (BET2 - B2) - DERIV(3) * (RBET - RB1)) * T
C
C        FORM NEW ESTIMATES OF G AND D
C
      G = G - U
      IF (B1 .EQ. ZERO .OR. G .LT. ZERO) G = ZERO
      D = D - Y
      IF (DABS(U) .GT. TT .OR. DABS(Y) .GT. TT) GOTO 80
C
C        END OF ITERATION
C
      DELTA = D
      XLAM = SIGMA / DSQRT(H2)
      IF (NEG) GOTO 130
      GAMMA = G
      GOTO 140
  130 GAMMA = -G
      HMU(1) = ONE - HMU(1)
  140 XI = XBAR - XLAM * HMU(1)
      RETURN
      END
C
      SUBROUTINE MOM(G, D, A, FAULT)
C
C        ALGORITHM AS 99.3  APPL. STATIST. (1976) VOL.25, P.180
C
C        EVALUATES FIRST SIX MOMENTS OF A JOHNSON
C        SB DISTRIBUTION, USING GOODWIN METHOD
C
      DOUBLE PRECISION A(6), B(6), C(6), G, D, ZZ, VV, RTTWO, RRTPI, 
     $  W, E, R,
     $  H, T, U, Y, X, V, F, Z, S, P, Q, AA, AB, DEXPA, DEXPB,
     $  ZERO, QUART, HALF, P75, ONE, TWO, THREE, DPI
      LOGICAL L, FAULT
C
      DATA ZZ, VV, LIMIT /1.0E-5, 1.0E-8, 5000/
C
C        RTTWO IS DSQRT(2.0)
C        RRTPI IS RECIPROCAL OF DSQRT(PI)
C        DEXPA IS A VALUE SUCH THAT DEXP(DEXPA) DOES NOT QUITE
C          CAUSE OVERFLOW
C        DEXPB IS A VALUE SUCH THAT 1.0 + DEXP(-DEXPB) MAY BE
C          TAKEN TO BE 1.0
C
C      DATA     RTTWO,        RRTPI, DEXPA, DEXPB
C     $  /1.414213562, 0.5641895835, 80.0, 23.7/
      DATA   DEXPA, DEXPB
     $     /80.0, 23.7/
      DATA ZERO, QUART, HALF,  P75, ONE, TWO, THREE
     $     /0.0,  0.25,  0.5, 0.75, 1.0, 2.0,   3.0/
C
C	Remove statement fn ZDABS, ZDEXP
C     ZDABS(X) = DABS(X)
C     ZDEXP(X) = DEXP(X)
C
      RTTWO = DSQRT(DBLE(2.0))
      DPI = 4.0D+00*ATAN(1.0D+00)
      RRTPI = 1.0D+00 / DSQRT(DPI)
C
      FAULT = .FALSE.
      DO 10 I = 1, 6
        C(I) = ZERO
   10 CONTINUE
      W = G / D
C
C        TRIAL VALUE OF H
C
      IF (W .GT. DEXPA) GOTO 140
      E = DEXP(W) + ONE
      R = RTTWO / D
      H = P75
      IF (D .LT. THREE) H = QUART * D
      K = 1
      GOTO 40
C
C        START OF OUTER LOOP
C
   20 K = K + 1
      IF (K .GT. LIMIT) GOTO 140
      DO 30 I = 1, 6
        C(I) = A(I)
   30 CONTINUE
C
C        NO CONVERGENCE YET - TRY SMALLER H
C
      H = HALF * H
   40 T = W
      U = T
      Y = H * H
      X = TWO * Y
      A(1) = ONE / E
      DO 50 I = 2, 6
      A(I) = A(I - 1) / E
   50 CONTINUE
      V = Y
      F = R * H
      M = 0
C
C        START OF INNER LOOP
C        TO EVALUATE INFINITE SERIES
C
   60 M = M + 1
      IF (M .GT. LIMIT) GOTO 140
      DO 70 I = 1, 6
        B(I) = A(I)
   70 CONTINUE
      U = U - F
      Z = ONE
      IF (U .GT. -DEXPB) Z = DEXP(U) + Z
      T = T + F
      L = T .GT. DEXPB

C
C	GPN added: If next line is not present then a compiler warning can be thrown
C		about S not being initialised. So, we have initialized S=1 here.
C		If (.NOT. L), then S gets set to the correct value in the line after.
C		If (L), then the value of S is ignored - must be, because it was not previously
C		initialized, so it does not matter than we have set it to something.
C		We have verified by a text search that the variable S does not appear later
C
C
      S = 1.0

      IF (.NOT. L) S = DEXP(T) + ONE
      P = DEXP(-V)
      Q = P
      DO 90 I = 1, 6
      AA = A(I)
      P = P / Z
      AB = AA
      AA = AA + P
      IF (AA .EQ. AB) GOTO 100
      IF (L) GOTO 80
      Q = Q / S
      AB = AA
      AA = AA + Q
      L = AA .EQ. AB
   80 A(I) = AA
   90 CONTINUE
  100 Y = Y + X
      V = V + Y
      DO 110 I = 1, 6
      IF (A(I) .EQ. ZERO) GOTO 140
      IF (DABS((A(I) - B(I)) / A(I)) .GT. VV) GOTO 60
  110 CONTINUE
C
C        END OF INNER LOOP
C
      V = RRTPI * H
      DO 120 I = 1, 6
        A(I) = V * A(I)
  120 CONTINUE
      DO 130 I = 1, 6
      IF (A(I) .EQ. ZERO) GOTO 140
      IF (DABS((A(I) - C(I)) / A(I)) .GT. ZZ) GOTO 20
  130 CONTINUE
C
C        END OF OUTER LOOP
C
      RETURN
  140 FAULT =.TRUE.
      RETURN
      END
CEND OF AS 99

