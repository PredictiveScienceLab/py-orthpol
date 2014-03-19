      SUBROUTINE POLY_EVAL( X, PHI, P, ALPHA, BETA, GAMMA )

      INTEGER               P
      DOUBLE PRECISION      X
      DOUBLE PRECISION      ALPHA( P ), BETA( P ), GAMMA( P ),
     $                      PHI( P )
cf2py double precision intent(in) :: x
cf2py double precision intent(out),depend(p),dimension(p) :: phi
cf2py integer intent(hide),depend(alpha) :: p=len(alpha)
cf2py double precision intent(in) :: alpha
cf2py double precision intent(in),depend(p),dimension(p) :: beta
cf2py double precision intent(in),depend(p),dimension(p) :: gamma
      PHI(1) = 1. / GAMMA(1)
      IF (P.GE.1) THEN
        PHI(2) = (X - ALPHA(1)) * (PHI(1) / GAMMA(1))
      END IF
      DO I=3,P
        PHI(I) = ((X - ALPHA(I-1)) * PHI(I-1) - BETA(I-1) * PHI(I-2)) /
     $           GAMMA(I)
      END DO

      RETURN

      END SUBROUTINE

      SUBROUTINE POLY_EVAL_ALL( N, X, PHI, P, ALPHA, BETA, GAMMA )

      INTEGER               N, P
      DOUBLE PRECISION      X( N ), PHI( N, P ), ALPHA( P ), BETA( P ),
     $                      GAMMA( P )
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py integer intent(hide),depend(alpha) :: p=len(alpha)
cf2py double precision intent(in) :: x
cf2py double precision intent(out),depend(n,p),dimension(n,p) :: phi
cf2py double precision intent(in) :: alpha
cf2py double precision intent(in),depend(p),dimension(p) :: beta
cf2py double precision intent(in),depend(p),dimension(p) :: gamma
      DO I=1,N
        CALL POLY_EVAL( X(I), PHI(I,:), P, ALPHA, BETA, GAMMA )
      END DO

      RETURN

      END SUBROUTINE

      SUBROUTINE POLY_DEVAL(X, DPHI, PHI, P, ALPHA, BETA, GAMMA)

      INTEGER               P
      DOUBLE PRECISION      X
      DOUBLE PRECISION      ALPHA( P ), BETA( P ), GAMMA( P ),
     $                      PHI( P ), DPHI( P )
cf2py double precision intent(in) :: x
cf2py double precision intent(out),depend(p),dimension(p) :: dphi
cf2py integer intent(hide),depend(alpha) :: p=len(alpha)
cf2py double precision intent(in), alpha
cf2py double precision intent(in),depend(p),dimension(p) :: beta
cf2py double precision intent(in),depend(p),dimension(p) :: gamma
cf2py double precision intent(hide),depend(p),dimension(p) :: phi

      CALL POLY_EVAL( X, PHI, P, ALPHA, BETA, GAMMA )
      DPHI(1) = 0.
      IF (P.GE.1) THEN
        DPHI(2) = PHI(1) / GAMMA(1)
      END IF
      DO I=3,P
        DPHI(I) = (PHI(I-1) + (X - ALPHA(I-1)) * DPHI(I-1)
     $             - BETA(I-1) * DPHI(I-2)) / GAMMA(I)
      END DO

      RETURN

      END SUBROUTINE

      SUBROUTINE POLY_DEVAL_ALL( N, X, DPHI, PHI, P, ALPHA, BETA,
     $                           GAMMA )

      INTEGER               N, P
      DOUBLE PRECISION      X( N ), PHI( N, P ), ALPHA( P ), BETA( P ),
     $                      GAMMA( P ), DPHI( N, P )
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py integer intent(hide),depend(alpha) :: p=len(alpha)
cf2py double precision intent(in) :: x
cf2py double precision intent(hide),depend(n,p),dimension(n,p) :: phi
cf2py double precision intent(out),depend(n,p),dimension(n,p) :: dphi
cf2py double precision intent(in) :: alpha
cf2py double precision intent(in),depend(p),dimension(p) :: beta
cf2py double precision intent(in),depend(p),dimension(p) :: gamma
      DO I=1,N
        CALL POLY_DEVAL( X(I), DPHI(I,:), PHI(I, :), P, ALPHA, BETA,
     $                   GAMMA )
      END DO

      RETURN

      END SUBROUTINE

      SUBROUTINE POLY_NORMALIZE( P, BETA, GAMMA )

      IMPLICIT NONE

      INTEGER               P
      DOUBLE PRECISION      BETA( P ), GAMMA( P )
cf2py integer intent(in),depend(beta) :: p=len(beta)
cf2py double precision intent(in,out,copy) :: beta
cf2py double precision intent(in,out,copy),depend(p),dimension(p) :: gamma
      INTEGER               I

      BETA(1) = SQRT(BETA(1))
      GAMMA(1) = BETA(1)
      DO I=1,P
        BETA(I) = SQRT(BETA(I) * GAMMA(I))
        GAMMA(I) = BETA(I)
      END DO

      RETURN

      END SUBROUTINE

      INTEGER FUNCTION N_CHOOSE_K( N, K )

      IMPLICIT NONE

      INTEGER               N, K

      INTEGER               NUM, DEN, I, L

      NUM = 1
      DEN = 1
      L = MIN(N-K, K)
      DO I=0,L-1
        NUM = NUM * (N - I)
        DEN = DEN * (I + L)
      END DO

      N_CHOOSE_K = NUM / DEN

      RETURN

      END FUNCTION
