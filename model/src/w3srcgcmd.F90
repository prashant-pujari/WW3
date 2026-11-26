!> @file
!> @brief Source-term package for short gravity-capillary waves.
!>
!> @author P. Pujari
!> @date   26-Nov-2026
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
!> @date   26-Nov-2026
!>
MODULE W3SRCGCMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
  !/                  |           P. Pujari               |                  
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         26-Nov-2026 |
  !/                  +-----------------------------------+
  !/
  !/    26-Oct-2026 : Origination of STGC source package (w3srcgcmd.F90)         (P. Pujari)
  !/    06-Nov-2026 : Update STGC documentation and source-package header        (P. Pujari)
  !/    11-Nov-2026 : Add molecular viscous dissipation.                         (P. Pujari)
  !/    20-Nov-2026 : Add short-gravity wave-breaking dissipation.               (P. Pujari)
  !/    24-Nov-2026 : Add breaker-generated freely propagating short-wave input. (P. Pujari)
  !/    26-Nov-2026 : Add direct wind input for short gravity waves.             (P. Pujari)
  !/    26-Nov-2026 : Add parasitic-capillary generation source.                 (P. Pujari)
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
  !      W3SINGC implements the Kudryavtsev et al. (2005) direct
  !      wind-input source for k < k_gamma. The wind-growth rate uses
  !      USTAR, Z0, DAIR and the actual wind direction U10DIR.
  !
  !      W3SINGC also implements the Yurovskaya et al. (2013)
  !      breaker-generated freely propagating short-wave input.
  !      Parent-wave growth follows the Kudryavtsev et al. (2005)
  !      beta formulation using U10DIR, and the generated source is
  !      isotropic in target direction.
  !
  !      W3SINGC also implements Yurovskaya et al. (2013)
  !      parasitic-capillary generation from same-direction parent
  !      saturation spectra.
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

  !> @brief Calculate direct, breaker-generated, and parasitic-capillary input.
  !> @details Direct Kudryavtsev et al. (2005) wind input applies for k < k_gamma. The Yurovskaya
  !>  et al. (2013) nonlocal generation source uses the same wind growth rate, beta.
  !>  Parent waves satisfy k_j < min(k_i/10,k_b), where k_b=k_gamma/4.
  !>  The Kudryavtsev et al. (2005) beta uses the actual wind direction WINDDIR=U10DIR. The parent
  !>  integral is applied identically to every target direction, and D_bgen=0.
  !>  Parasitic-capillary input uses same-direction mapped parent saturation and
  !>  has zero local diagonal derivative.
  SUBROUTINE W3SINGC (A, CG, WN2, UABS, USTAR, WINDDIR, Z0, CD, DAIR, &
       TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S, D)

    USE CONSTANTS, ONLY: GRAV, DWAT, KAPPA, PI
    USE W3DISPMD, ONLY: SURFACE_TENSION
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DSII, DTH, ECOS, ESIN

    IMPLICIT NONE

    REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN2(NSPEC), UABS, USTAR, WINDDIR, Z0, CD, DAIR
    REAL, INTENT(OUT) :: TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S(NSPEC), D(NSPEC)

    INTEGER :: IK, JK, ITH, IS, IST, ISLO, ISHI, JLOW
    REAL, PARAMETER :: CBGEN = 4.5E-3, NGPC = 10.0
    REAL, PARAMETER :: APC = 2.0E-3, CBAR = 0.03
    REAL :: KGAMMA, KBREAK, KI, KJ, KBM, CJ, COSPHI, CBETA, BETA, BGENINT
    REAL :: KMIN, KMAX, KPAR, KLO, KHI, WPAR, KPCL, KPCH, XPAR, FPAR
    REAL :: NPC, ALPHAPC, FPCLO, FPCHI, PHIPC, BLO, BHI, BPAR
    REAL :: SWIND(NSPEC), DWIND(NSPEC), SBGEN(NSPEC), SPC(NSPEC)

    ! Yurovskaya et al. (2013) breaker-generation wavenumber limit.
    KGAMMA = SQRT(DWAT * GRAV / SURFACE_TENSION)
    KBREAK = 0.25 * KGAMMA
    KPCL = 1.5 * KGAMMA
    KPCH = KGAMMA**2 / KBREAK
    KMIN = WN2(1)
    KMAX = WN2(1 + (NK-1)*NTH)

    TAUWX  = 0.0
    TAUWY  = 0.0
    TAUNWX = 0.0
    TAUNWY = 0.0
    CHARN  = 0.0
    SWIND  = 0.0
    DWIND  = 0.0
    SBGEN  = 0.0
    SPC    = 0.0

    ! Kudryavtsev et al. (2005) direct wind input in the short-gravity range.
    DO IK = 1, NK
      DO ITH = 1, NTH
        IS = ITH + (IK-1) * NTH
        KJ = WN2(IS)
        IF (KJ .LT. KGAMMA) THEN
          CJ = SIG(IK) / KJ
          COSPHI = ECOS(IS)*COS(WINDDIR) + ESIN(IS)*SIN(WINDDIR)
          CBETA = 1.5 * DAIR / DWAT * ( LOG(PI/(KJ*Z0))/KAPPA - CJ/USTAR )
          BETA = CBETA * (USTAR/CJ)**2 * COSPHI * ABS(COSPHI)
          DWIND(IS) = BETA * SIG(IK)
          SWIND(IS) = DWIND(IS) * A(IS)
        END IF
      END DO
    END DO

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

    ! Yurovskaya et al. (2013) parasitic-capillary generation. The mapped
    ! parent must lie on the resolved wavenumber grid; no extrapolation is used.
    DO IK = 1, NK
      IST = 1 + (IK-1) * NTH
      KI = WN2(IST)
      KPAR = KGAMMA**2 / KI
      JLOW = 0

      IF (KPAR .GE. KMIN .AND. KPAR .LE. KMAX) THEN
        ! Treat the upper endpoint separately to avoid an NK+1 parent access.
        IF (KPAR .EQ. KMAX) THEN
          JLOW = NK
        ELSE
          DO JK = 1, NK-1
            KLO = WN2(1 + (JK-1)*NTH)
            KHI = WN2(1 + JK*NTH)
            IF (KPAR .GE. KLO .AND. KPAR .LE. KHI) THEN
              JLOW = JK
              WPAR = (LOG(KPAR) - LOG(KLO)) / (LOG(KHI) - LOG(KLO))
              EXIT
            END IF
          END DO
        END IF
      END IF

      IF (JLOW .GT. 0) THEN
        XPAR = KPAR / KGAMMA
        FPAR = 0.5 * (1.0 + TANH(2.0*LOG(4.0*XPAR)))
        NPC = NGPC / (1.0 + (NGPC-1.0)*FPAR)
        ALPHAPC = APC * CBAR**(-1.0/NPC)
        FPCLO = 0.5 * (1.0 + TANH(2.0*LOG(4.0*KI/KPCL)))
        FPCHI = 0.5 * (1.0 + TANH(2.0*LOG(4.0*KI/KPCH)))
        PHIPC = FPCLO - FPCHI

        DO ITH = 1, NTH
          IS = ITH + (IK-1) * NTH
          IF (JLOW .EQ. NK) THEN
            ISLO = ITH + (NK-1)*NTH
            BPAR = SIG(NK) * WN2(ISLO)**3 * A(ISLO)
          ELSE
            ISLO = ITH + (JLOW-1)*NTH
            ISHI = ITH + JLOW*NTH
            BLO = SIG(JLOW) * WN2(ISLO)**3 * A(ISLO)
            BHI = SIG(JLOW+1) * WN2(ISHI)**3 * A(ISHI)
            BPAR = (1.0-WPAR)*BLO + WPAR*BHI
          END IF

          ! Guard against a negative numerical parent saturation before
          ! evaluating the generally non-integer power; BPAR=0 gives SPC=0.
          IF (BPAR .GT. 0.0) THEN
            SPC(IS) = PHIPC / KI**3 * BPAR * (BPAR/ALPHAPC)**NPC
          END IF
        END DO
      END IF
    END DO

    S = SWIND + SBGEN + SPC
    D = DWIND

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
