C                  CVM Class Library
C                  http://cvmlib.com
C
C          Copyright Sergei Nikolaev 1992-2008
C Distributed under the Boost Software License, Version 1.0.
C    (See accompanying file LICENSE_1_0.txt or copy at
C          http://www.boost.org/LICENSE_1_0.txt)
C
C
C     Matrix infinity norm routines
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     A   - matrix (real)(input)
C     LDA - leading dimesion of A (int)(input)

      REAL*8 FUNCTION DINFNM (M, N, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DINFNM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      REAL*8    A(LDA*N)
      INTEGER*4 I
      REAL*8    S
      INTEGER*4 IDAMAX

      DINFNM = 0.D0
      IF (M .EQ. LDA) THEN
          DINFNM = DABS (A (IDAMAX (M * N, A, 1)))
      ELSE
          DO 20 I = 0, N-1
              S = DABS (A (IDAMAX (M, A(I*LDA+1), 1)))
              IF (S .GT. DINFNM) DINFNM = S
20        CONTINUE
      END IF
      RETURN
      END !FUNCTION DINFNM


   
