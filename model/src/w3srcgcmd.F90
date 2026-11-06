!> @file
!> @brief Source-term package for short gravity-capillary waves.
!>
!> @author P. Pujari
!> @date   06-Nov-2026
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /

!> @brief Source-term package for short gravity-capillary waves.
!>
!> @details Source-term package developed for the representation of
!>  short gravity-capillary waves in WAVEWATCH III. The package is
!>  intended to include wind input, viscous dissipation, wave-breaking
!>  dissipation, short-wave generation by breaking waves, and parasitic
!>  capillary-wave generation. Nonlinear four-wave and three-wave
!>  interactions are handled separately from this package.
!>
!> @author P. Pujari
!> @date   06-Nov-2026
!>
MODULE W3SRCGCMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
  !/                  |           P. Pujari               |                  |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         06-Nov-2026 |
  !/                  +-----------------------------------+
  !/
  !/    26-Oct-2026 : Origination of STGC source package (w3srcgcmd.F90)         (P. Pujari) 
  !/    06-Nov-2026 : Add molecular viscous dissipation.                         (P. Pujari) 
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !      Source-term package for the representation of short
  !      gravity-capillary waves in WAVEWATCH III.
  !
  !      The package is intended to contain source terms associated
  !      with:
  !
  !        - wind input;
  !        - molecular viscous dissipation;
  !        - wave-breaking dissipation;
  !        - generation of freely propagating short waves by breaking;
  !        - generation of parasitic capillary waves.
  !
  !      Nonlinear four-wave and three-wave interaction source terms
  !      are treated separately from the STGC source package.
  !
  !      References:
  !       Jähne and Riemer            1990: JGR   95(7)     11,531-11546
  !       Makin and Kudryavtsev       1999: JGR   104(4)    7613-7623
  !       Kudryavtsev et al.          1999: JGR   104(4)    7625-7639
  !       Caudal                      2002: JGR   107(10)
  !       Kudryavtsev et al.          2003: JGR   108(3)  
  !       Kudryavtsev and Johannessen 2004: GRL   31  
  !       Kudryavtsev et al.          2005: JGR   110
  !       Donelan et al.              2012: JGR   117
  !       Yurovskaya et al.           2013: JGR   118       4380-4394
  !
  !  2. Variables and types :
  !
  !      To be defined as individual source terms are implemented.
  !
  !  3. Subroutines and functions :
  !
  !      Name       Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SPRGC    Subr. Public   Integral parameter calculations.
  !      W3SINGC    Subr. Public   Gravity-capillary input/source terms.
  !      W3SDSGC    Subr. Public   Gravity-capillary dissipation terms.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !      Name       Type  Module   Description
  !     ----------------------------------------------------------------
  !      STRACE     Subr. W3SERVMD Subroutine tracing.
  !     ----------------------------------------------------------------
  !
  !  5. Remarks :
  !
  !      Source terms are introduced incrementally and documented
  !      individually as the STGC package is developed.
  !
  !  6. Switches :
  !
  !      !/S    Enable subroutine tracing.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
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
