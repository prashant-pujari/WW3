!> @file
!> @brief Source-term package for short gravity-capillary waves.
!>
!> @author P. Pujari
!> @date   24-Nov-2026
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
!> @date   24-Nov-2026
!>
MODULE W3SRCGCMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
  !/                  |           P. Pujari               |                  
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         24-Nov-2026 |
  !/                  +-----------------------------------+
  !/
  !/    26-Oct-2026 : Origination of STGC source package (w3srcgcmd.F90)         (P. Pujari)
  !/    06-Nov-2026 : Update STGC documentation and source-package header        (P. Pujari)
  !/    11-Nov-2026 : Add molecular viscous dissipation.                         (P. Pujari)
  !/    20-Nov-2026 : Add short-gravity wave-breaking dissipation.               (P. Pujari)
  !/    24-Nov-2026 : Add breaker-generated freely propagating short-wave input. (P. Pujari)
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
  !      W3SINGC implements the Yurovskaya et al. (2013) breaker-generated
  !      freely propagating short-wave input. Parent-wave growth follows the
  !      Kudryavtsev et al. (2005) beta formulation using U10DIR and the source is isotropic in
  !      target direction.
  !
  !      Molecular viscous and short-gravity wave-breaking dissipation are
  !      implemented in W3SDSGC.
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

  !> @brief Calculate breaker-generated freely propagating short-wave input.
  !> @details Uses the Yurovskaya et al. (2013) nonlocal generation source.
  !>  Parent waves satisfy k_j < min(k_i/10,k_b), where k_b=k_gamma/4.
  !>  The K2005 beta uses the actual wind direction WINDDIR=U10DIR. The parent
  !>  integral is applied identically to every target direction, and D_bgen=0.
  SUBROUTINE W3SINGC (A, CG, WN2, UABS, USTAR, WINDDIR, Z0, CD, DAIR, &
       TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S, D)

    USE CONSTANTS, ONLY: GRAV, DWAT, KAPPA, PI
    USE W3DISPMD, ONLY: SURFACE_TENSION
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DSII, DTH, ECOS, ESIN

    IMPLICIT NONE

    REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN2(NSPEC), UABS, USTAR, WINDDIR, Z0, CD, DAIR
    REAL, INTENT(OUT) :: TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S(NSPEC), D(NSPEC)

    INTEGER :: IK, JK, ITH, IS, IST
    REAL, PARAMETER :: CBGEN = 4.5E-3
    REAL :: KGAMMA, KBREAK, KI, KJ, KBM, CJ, COSPHI, CBETA, BETA, BGENINT
    REAL :: SBGEN(NSPEC)

    ! Yurovskaya et al. (2013) breaker-generation wavenumber limit.
    KGAMMA = SQRT(DWAT * GRAV / SURFACE_TENSION)
    KBREAK = 0.25 * KGAMMA

    TAUWX  = 0.0
    TAUWY  = 0.0
    TAUNWX = 0.0
    TAUNWY = 0.0
    CHARN  = 0.0
    SBGEN  = 0.0
    D      = 0.0

    ! Integrate over lower-wavenumber parent waves. DSII is the absolute
    ! intrinsic-frequency interval, and the resulting input is isotropic in
    ! target direction.
    DO IK = 1, NK
      IST = 1 + (IK-1) * NTH
      KI = WN2(IST)
      KBM = MIN(KI/10.0, KBREAK)
      BGENINT = 0.0
      DO JK = 1, NK
        DO ITH = 1, NTH
          IS = ITH + (JK-1) * NTH
          KJ = WN2(IS)
          IF (KJ .LT. KBM) THEN
            CJ = SIG(JK) / KJ
            COSPHI = ECOS(IS)*COS(WINDDIR) + ESIN(IS)*SIN(WINDDIR)
            CBETA = 1.5 * DAIR / DWAT * ( LOG(PI/(KJ*Z0))/KAPPA - CJ/USTAR )
            BETA = CBETA * (USTAR/CJ)**2 * COSPHI * ABS(COSPHI)
            BGENINT = BGENINT + BETA * SIG(JK)**2 * KJ**2 / CG(JK) * &
                 A(IS) * DSII(JK) * DTH
          END IF
        END DO
      END DO
      DO ITH = 1, NTH
        IST = ITH + (IK-1) * NTH
        SBGEN(IST) = CBGEN * BGENINT / (KI**3 * SIG(IK))
      END DO
    END DO

    S = SBGEN

  END SUBROUTINE W3SINGC

  !> @brief Calculate molecular viscous and short-gravity breaking dissipation.
  !> @details Implements S_visc = -4 nu k**2 A and the Kudryavtsev et al.
  !>  (1999) short-gravity breaking source using the local wavenumber.
  SUBROUTINE W3SDSGC (A, CG, WN, S, D)

    USE CONSTANTS, ONLY: GRAV, DWAT
    USE W3DISPMD, ONLY: SURFACE_TENSION
    USE W3GDATMD, ONLY: NK, NTH, NSPEC

    IMPLICIT NONE

    REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN(NK)
    REAL, INTENT(OUT) :: S(NSPEC), D(NSPEC)

    INTEGER :: IK, ITH, IS
    REAL, PARAMETER :: NUWATER = 1.0035E-6
    REAL, PARAMETER :: NG = 2.0, ALPHAG = 0.03
    REAL :: KGAMMA, KBREAK, BREAKCOEF
    REAL :: DVISC(NK), SVISC(NSPEC), SBREAK(NSPEC), DBREAK(NSPEC)

    ! Molecular viscous dissipation: D_visc = -4 nu k**2.
    DVISC = -4.0 * NUWATER * WN**2

    ! Short-gravity breaking range and dissipation parameters.
    KGAMMA = SQRT(DWAT * GRAV / SURFACE_TENSION)
    KBREAK = 0.25 * KGAMMA

    ! Apply the frequency-dependent coefficient to every direction.
    D     = 0.0
    SVISC = 0.0
    SBREAK = 0.0
    DBREAK = 0.0
    DO IK = 1, NK
      DO ITH = 1, NTH
        IS = ITH + (IK-1) * NTH
        SVISC(IS) = DVISC(IK) * A(IS)
        IF (WN(IK) .LT. KBREAK .AND. A(IS) .GT. 0.0) THEN
          BREAKCOEF = GRAV**1.5 / ALPHAG**NG * WN(IK)**(3.0*NG + 1.5)
          SBREAK(IS) = -BREAKCOEF * A(IS)**(NG + 1.0)
          DBREAK(IS) = -(NG + 1.0) * BREAKCOEF * A(IS)**NG
        END IF
        D(IS) = DVISC(IK) + DBREAK(IS)
      END DO
    END DO

    S = SVISC + SBREAK

  END SUBROUTINE W3SDSGC

END MODULE W3SRCGCMD
