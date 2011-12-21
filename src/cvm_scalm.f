C                  CVM Class Library
C                  http://cvmlib.com
C
C          Copyright Sergei Nikolaev 1992-2008
C Distributed under the Boost Software License, Version 1.0.
C    (See accompanying file LICENSE_1_0.txt or copy at
C          http://www.boost.org/LICENSE_1_0.txt)
C
C
C     Matrix scaling routines
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     S   - scale factor (real)(input)
C     A   - matrix to be scaled (real)(input, output)
C     LDA - leading dimesion of A (int)(input)

      SUBROUTINE DSCALM (M, N, S, A, LDA) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DSCALM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA
      REAL*8    S
      REAL*8    A(LDA*N)
      INTEGER*4 I

      IF (M .EQ. LDA) THEN
          CALL DSCAL (M * N, S, A, 1)
      ELSE
          DO 10 I = 0, N-1
              CALL DSCAL (M, S, A(I*LDA+1), 1)
10        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE DSCALM


