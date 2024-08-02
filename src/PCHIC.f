
      SUBROUTINE CHFEV (X1, X2, F1, F2, D1, D2, NE, XE, FE, NEXT, IERR)
C***BEGIN PROLOGUE  CHFEV
C***PURPOSE  Evaluate a cubic polynomial given in Hermite form at an
C            array of points.  While designed for use by PCHFE, it may
C            be useful directly as an evaluator for a piecewise cubic
C            Hermite function in applications, such as graphing, where
C            the interval is known in advance.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3
C***TYPE      SINGLE PRECISION (CHFEV-S, DCHFEV-D)
C***KEYWORDS  CUBIC HERMITE EVALUATION, CUBIC POLYNOMIAL EVALUATION,
C             PCHIP
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          CHFEV:  Cubic Hermite Function EValuator
C
C     Evaluates the cubic polynomial determined by function values
C     F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
C     XE(J), J=1(1)NE.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  NE, NEXT(2), IERR
C        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
C
C        CALL  CHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
C
C   Parameters:
C
C     X1,X2 -- (input) endpoints of interval of definition of cubic.
C           (Error return if  X1.EQ.X2 .)
C
C     F1,F2 -- (input) values of function at X1 and X2, respectively.
C
C     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real array of points at which the function is to be
C           evaluated.  If any of the XE are outside the interval
C           [X1,X2], a warning error is returned in NEXT.
C
C     FE -- (output) real array of values of the cubic function defined
C           by  X1,X2, F1,F2, D1,D2  at the points  XE.
C
C     NEXT -- (output) integer array indicating number of extrapolation
C           points:
C            NEXT(1) = number of evaluation points to left of interval.
C            NEXT(2) = number of evaluation points to right of interval.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if NE.LT.1 .
C              IERR = -2  if X1.EQ.X2 .
C                (The FE-array has not been changed in either case.)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811019  DATE WRITTEN
C   820803  Minor cosmetic changes for release 1.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890703  Corrected category record.  (WRB)
C   890703  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  CHFEV
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change CHFEV to DCHFEV wherever it occurs,
C        b. Change the real declaration to double precision, and
C        c. Change the constant ZERO to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  NE, NEXT(2), IERR
      REAL  X1, X2, F1, F2, D1, D2, XE(*), FE(*)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I
      REAL  C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA, ZERO
      SAVE ZERO
      DATA  ZERO /0./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  CHFEV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
C
C  INITIALIZE.
C
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = MIN(ZERO, H)
      XMA = MAX(ZERO, H)
C
C  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
C
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
C                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C3 = (DEL1 + DEL2)/H
C                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
C
C  EVALUATION LOOP.
C
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
C          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
C        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -1
      print*, 'error in PCHIC.f'
      print*, 'NUMBER OF EVALUATION POINTS LESS THAN ONE'
      RETURN
C
 5002 CONTINUE
C     X1.EQ.X2 RETURN.
      IERR = -2
      print*, 'error in PCHIC.f'
      print*, 'INTERVAL ENDPOINTS EQUAL'
      RETURN
C------------- LAST LINE OF CHFEV FOLLOWS ------------------------------
      END 


c    !====================================================================


*DECK PCHFE
      SUBROUTINE PCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
C***BEGIN PROLOGUE  PCHFE
C***PURPOSE  Evaluate a piecewise cubic Hermite function at an array of
C            points.  May be used by itself for Hermite interpolation,
C            or as an evaluator for PCHIM or PCHIC.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3
C***TYPE      SINGLE PRECISION (PCHFE-S, DPCHFE-D)
C***KEYWORDS  CUBIC HERMITE EVALUATION, HERMITE INTERPOLATION, PCHIP,
C             PIECEWISE CUBIC EVALUATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          PCHFE:  Piecewise Cubic Hermite Function Evaluator
C
C     Evaluates the cubic Hermite function defined by  N, X, F, D  at
C     the points  XE(J), J=1(1)NE.
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, NE, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
C        LOGICAL  SKIP
C
C        CALL  PCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on normal return.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real array of points at which the function is to be
C           evaluated.
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XE are increasing relative to X;
C              that is,   XE(J) .GE. X(I)
C              implies    XE(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XE are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C     FE -- (output) real array of values of the cubic Hermite function
C           defined by  N, X, F, D  at the points  XE.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NE.LT.1 .
C             (The FE-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHFEV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811020  DATE WRITTEN
C   820803  Minor cosmetic changes for release 1.
C   870707  Minor cosmetic changes to prologue.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  PCHFE
C  Programming notes:
C
C     1. To produce a double precision version, simply:
C        a. Change PCHFE to DPCHFE, and CHFEV to DCHFEV, wherever they
C           occur,
C        b. Change the real declaration to double precision,
C
C     2. Most of the coding between the call to CHFEV and the end of
C        the IR-loop could be eliminated if it were permissible to
C        assume that XE is ordered relative to X.
C
C     3. CHFEV does not assume that X1 is less than X2.  thus, it would
C        be possible to write a version of PCHFE that assumes a strict-
C        ly decreasing X-array by simply running the IR-loop backwards
C        (and reversing the order of appropriate tests).
C
C     4. The present code has a minor bug, which I have decided is not
C        worth the effort that would be required to fix it.
C        If XE contains points in [X(N-1),X(N)], followed by points .LT.
C        X(N-1), followed by points .GT.X(N), the extrapolation points
C        will be counted (at least) twice in the total returned in IERR.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, NE, IERR
      REAL  X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHFE
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
!      SKIP = .TRUE.
C
C  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
C                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
C
C     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
C
         IF (JFIRST .GT. NE)  GO TO 5000
C
C     LOCATE ALL POINTS IN INTERVAL.
C
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
C
C     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
C
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
C
   40    CONTINUE
         NJ = J - JFIRST
C
C     SKIP EVALUATION IF NO POINTS IN INTERVAL.
C
         IF (NJ .EQ. 0)  GO TO 50
C
C     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
C
C       ----------------------------------------------------------------
        CALL CHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR),
     *              NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)
C       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
C
         IF (NEXT(2) .EQ. 0)  GO TO 42
C        IF (NEXT(2) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
C           RIGHT OF X(IR).
C
            IF (IR .LT. N)  GO TO 41
C           IF (IR .EQ. N)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
C           ELSE
C              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
C           ENDIF
C        ENDIF
   42    CONTINUE
C
         IF (NEXT(1) .EQ. 0)  GO TO 49
C        IF (NEXT(1) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
C           LEFT OF X(IR-1).
C
            IF (IR .GT. 2)  GO TO 43
C           IF (IR .EQ. 2)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
C           ELSE
C              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
C              EVALUATION INTERVAL.
C
C              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
C              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
C                     IN CHFEV.
               GO TO 5005
C
   45          CONTINUE
C              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
C
C              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
C              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
C
   47          CONTINUE
C              AT THIS POINT, EITHER  XE(J) .LT. X(1)
C                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
C              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
C              CYCLING.
               IR = MAX(1, I-1)
C           ENDIF
C        ENDIF
   49    CONTINUE
C
         JFIRST = J
C
C     END OF IR-LOOP.
C
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      print*, 'error in PCHIC.f'
      print*, 'NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      print*, 'error in PCHIC.f'
      print*, 'INCREMENT LESS THAN ONE'
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      print*, 'error in PCHIC.f'
      print*, 'X-ARRAY NOT STRICTLY INCREASING'
      RETURN
C
 5004 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -4
      print*, 'error in PCHIC.f'
      print*, 'NUMBER OF EVALUATION POINTS LESS THAN ONE'
      RETURN
C
 5005 CONTINUE
C     ERROR RETURN FROM CHFEV.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      print*, 'error in PCHIC.f'
      print*, 'ERROR RETURN FROM CHFEV -- FATAL'
      RETURN
C------------- LAST LINE OF PCHFE FOLLOWS ------------------------------
      END


c     !====================================================================


*DECK PCHIM
      SUBROUTINE PCHIM (N, X, F, D, INCFD, IERR)
C***BEGIN PROLOGUE  PCHIM
C***PURPOSE  Set derivatives needed to determine a monotone piecewise
C            cubic Hermite interpolant to given data.  Boundary values
C            are provided which are compatible with monotonicity.  The
C            interpolant will have an extremum at each point where mono-
C            tonicity switches direction.  (See PCHIC if user control is
C            desired over boundary or switch conditions.)
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (PCHIM-S, DPCHIM-D)
C***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
C             PCHIP, PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          PCHIM:  Piecewise Cubic Hermite Interpolation to
C                  Monotone data.
C
C     Sets derivatives needed to determine a monotone piecewise cubic
C     Hermite interpolant to the data given in X and F.
C
C     Default boundary conditions are provided which are compatible
C     with monotonicity.  (See PCHIC if user control of boundary con-
C     ditions is desired.)
C
C     If the data are only piecewise monotonic, the interpolant will
C     have an extremum at each point where monotonicity switches direc-
C     tion.  (See PCHIC if user control is desired in such cases.)
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by PCHFE or PCHFD.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N)
C
C        CALL  PCHIM (N, X, F, D, INCFD, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
C           PCHIM is designed for monotonic data, but it will work for
C           any F-array.  It will force extrema at points where mono-
C           tonicity switches direction.  If some other treatment of
C           switch points is desired, PCHIC should be used instead.
C                                     -----
C     D -- (output) real array of derivative values at the data points.
C           If the data are monotonic, these values will determine a
C           a monotone cubic Hermite function.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that IERR switches in the direction
C                 of monotonicity were detected.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F. N. Fritsch and J. Butland, A method for construc-
C                 ting local monotone piecewise cubic interpolants, SIAM
C                 Journal on Scientific and Statistical Computing 5, 2
C                 (June 1984), pp. 300-304.
C               2. F. N. Fritsch and R. E. Carlson, Monotone piecewise
C                 cubic interpolation, SIAM Journal on Numerical Ana-
C                 lysis 17, 2 (April 1980), pp. 238-246.
C***ROUTINES CALLED  PCHST, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811103  DATE WRITTEN
C   820201  1. Introduced  PCHST  to reduce possible over/under-
C             flow problems.
C           2. Rearranged derivative formula for same reason.
C   820602  1. Modified end conditions to be continuous functions
C             of data when monotonicity switches in next interval.
C           2. Modified formulas so end conditions are less prone
C             of over/underflow problems.
C   820803  Minor cosmetic changes for release 1.
C   870813  Updated Reference 1.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890703  Corrected category record.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920429  Revised format and order of references.  (WRB,FNF)
C***END PROLOGUE  PCHIM
C  Programming notes:
C
C     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
C        either argument is zero, +1 if they are of the same sign, and
C        -1 if they are of opposite sign.
C     2. To produce a double precision version, simply:
C        a. Change PCHIM to DPCHIM wherever it occurs,
C        b. Change PCHST to DPCHST wherever it occurs,
C        c. Change all references to the Fortran intrinsics to their
C           double precision equivalents,
C        d. Change the real declarations to double precision, and
C        e. Change the constants ZERO and THREE to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IERR
      REAL  X(*), F(INCFD,*), D(INCFD,*)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, NLESS1
      REAL  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE,
     *      H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
      SAVE ZERO, THREE
      REAL  PCHST
      DATA  ZERO /0./,  THREE /3./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHIM
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
C
C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
C
C  NORMAL CASE  (N .GE. 3).
C
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
C
C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
      ENDIF
C
C  LOOP THROUGH INTERIOR POINTS.
C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
C
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
C
C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
C
         D(1,I) = ZERO
         IF ( PCHST(DEL1,DEL2) )  42, 41, 45
C
C        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
C
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( PCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
C
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = MAX( ABS(DEL1), ABS(DEL2) )
         DMIN = MIN( ABS(DEL1), ABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
C
   50 CONTINUE
C
C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      print*, 'error in PCHIC.f'
      print*, 'NUMBER OF DATA POINTS LESS THAN TWO'
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      print*, 'error in PCHIC.f'
      print*, 'INCREMENT LESS THAN ONE'
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      print*, 'error in PCHIC.f'
      print*, 'X-ARRAY NOT STRICTLY INCREASING'
      RETURN
C------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
      END


c     !====================================================================


*DECK PCHST
      REAL FUNCTION PCHST (ARG1, ARG2)
C***BEGIN PROLOGUE  PCHST
C***SUBSIDIARY
C***PURPOSE  PCHIP Sign-Testing Routine
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (PCHST-S, DPCHST-D)
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C         PCHST:  PCHIP Sign-Testing Routine.
C
C     Returns:
C        -1. if ARG1 and ARG2 are of opposite sign.
C         0. if either argument is zero.
C        +1. if ARG1 and ARG2 are of the same sign.
C
C     The object is to do this without multiplying ARG1*ARG2, to avoid
C     possible over/underflow problems.
C
C  Fortran intrinsics used:  SIGN.
C
C***SEE ALSO  PCHCE, PCHCI, PCHCS, PCHIM
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811103  DATE WRITTEN
C   820805  Converted to SLATEC library version.
C   870813  Minor cosmetic changes.
C   890411  Added SAVE statements (Vers. 3.2).
C   890411  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
C   930503  Improved purpose.  (FNF)
C***END PROLOGUE  PCHST
C
C**End
C
C  DECLARE ARGUMENTS.
C
      REAL  ARG1, ARG2
C
C  DECLARE LOCAL VARIABLES.
C
      REAL  ONE, ZERO
      SAVE ZERO, ONE
      DATA  ZERO /0./,  ONE /1./
C
C  PERFORM THE TEST.
C
C***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO
C
      RETURN
C------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END


c     !====================================================================


