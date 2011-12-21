C                  CVM Class Library
C                  http://cvmlib.com
C
C          Copyright Sergei Nikolaev 1992-2008
C Distributed under the Boost Software License, Version 1.0.
C    (See accompanying file LICENSE_1_0.txt or copy at
C          http://www.boost.org/LICENSE_1_0.txt)
C
C
C     Matrix axpy routines
C
C     The ?axpy routines perform a vector-vector operation defined as
C     y := a*x + y
C
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     A   - multiplier (real)(input)
C     X   - source matrix (real)(input)
C     LDX - leading dimesion of X (int)(input)
C     Y   - destination matrix (real)(output)
C     LDY - leading dimesion of Y (int)(input)

      SUBROUTINE DAXPYM (M, N, A, X, LDX, Y, LDY)
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DAXPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDX, LDY
      REAL*8    A, X(LDX*N), Y(LDY*N)
      INTEGER*4 I

      IF (M .EQ. LDX .AND. M .EQ. LDY) THEN
          CALL DAXPY (M * N, A, X, 1, Y, 1)
      ELSE
          DO 20 I = 0, N-1
              CALL DAXPY (M, A, X(I*LDX+1), 1, Y(I*LDY+1), 1)
20        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE DAXPYM


