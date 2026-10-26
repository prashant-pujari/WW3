!> @file
!> @brief Placeholder source-term package for gravity-capillary waves.

#include "w3macros.h"
!/ ------------------------------------------------------------------- /

!> @brief Placeholder source-term package for future gravity-capillary physics.
MODULE W3SRCGCMD

  IMPLICIT NONE

  PUBLIC :: W3SPRGC, W3SINGC, W3SDSGC

CONTAINS

  !> @brief Calculate generic spectral moments required by source integration.
  !> @details This contains no gravity-capillary source physics.
  SUBROUTINE W3SPRGC (A, CG, WN, EMEAN, FMEAN, WNMEAN, AMAX, FP)

    USE CONSTANTS, ONLY: TPIINV
    USE W3GDATMD, ONLY: NK, NTH, SIG, DTH, DDEN, DSII

    IMPLICIT NONE

    REAL, INTENT(IN)  :: A(NTH,NK), CG(NK), WN(NK)
    REAL, INTENT(OUT) :: EMEAN, FMEAN, WNMEAN, AMAX, FP

    REAL :: EB(NK)
    REAL, PARAMETER :: HSMIN = 0.05

    ! Generic spectral moments are needed by the common source driver.
    EB     = SUM(A,1) * DDEN / CG
    AMAX   = MAXVAL(A)
    EMEAN  = SUM(EB)
    FMEAN  = SUM(EB / SIG(1:NK))
    WNMEAN = SUM(EB / SQRT(WN))

    FMEAN  = TPIINV * EMEAN / MAX(1.0E-7, FMEAN)
    WNMEAN = (EMEAN / MAX(1.0E-7, WNMEAN))**2

    FP = 0.0
    IF (4.0 * SQRT(EMEAN) .GT. HSMIN) THEN
      EB = SUM(A,1) * SIG(1:NK) / CG * DTH
      FP = SUM(SIG(1:NK) * EB**4 * DSII) / MAX(1.0E-10, SUM(EB**4 * DSII))
      FP = FP * TPIINV
    END IF

  END SUBROUTINE W3SPRGC

  !> @brief Return zero placeholder wind input and diagonal derivative.
  SUBROUTINE W3SINGC (A, CG, WN2, UABS, USTAR, USDIR, CD, DAIR, &
       TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S, D)

    USE W3GDATMD, ONLY: NK, NSPEC

    IMPLICIT NONE

    REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN2(NSPEC), UABS, USTAR, USDIR, CD, DAIR
    REAL, INTENT(OUT) :: TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S(NSPEC), D(NSPEC)

    TAUWX  = 0.0
    TAUWY  = 0.0
    TAUNWX = 0.0
    TAUNWY = 0.0
    CHARN  = 0.0
    S      = 0.0
    D      = 0.0

  END SUBROUTINE W3SINGC

  !> @brief Return zero placeholder dissipation and diagonal derivative.
  SUBROUTINE W3SDSGC (A, CG, WN, S, D)

    USE W3GDATMD, ONLY: NK, NSPEC

    IMPLICIT NONE

    REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN(NK)
    REAL, INTENT(OUT) :: S(NSPEC), D(NSPEC)

    S = 0.0
    D = 0.0

  END SUBROUTINE W3SDSGC

END MODULE W3SRCGCMD
