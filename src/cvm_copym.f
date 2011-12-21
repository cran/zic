C                  CVM Class Library
C                  http://cvmlib.com
C
C          Copyright Sergei Nikolaev 1992-2008
C Distributed under the Boost Software License, Version 1.0.
C    (See accompanying file LICENSE_1_0.txt or copy at
C          http://www.boost.org/LICENSE_1_0.txt)
C
C
C     Matrix copy routines
C
C     Input/Output parameters:
C
C     M   - rows (int)(input)
C     N   - columns (int)(input)
C     A   - source matrix (real)(input)
C     B   - destination matrix (real)(output)
C     LDA - leading dimesion of A (int)(input)
C     LDB - leading dimesion of B (int)(input)

      SUBROUTINE DCOPYM (M, N, A, LDA, B, LDB) 
CDEC$ IF DEFINED (FTN_EXPORTS)
CDEC$     ATTRIBUTES DLLEXPORT::DCOPYM
CDEC$ ENDIF
      INTEGER*4 M, N, LDA, LDB
      REAL*8    A(LDA*N), B(LDB*N)
      INTEGER*4 I

      IF (M .EQ. LDA .AND. M .EQ. LDB) THEN
          CALL DCOPY (M * N, A, 1, B, 1)
      ELSE
          DO 20 I = 0, N-1
              CALL DCOPY (M, A(I*LDA+1), 1, B(I*LDB+1), 1)
20        CONTINUE
      END IF
      RETURN
      END !SUBROUTINE DCOPYM

   
