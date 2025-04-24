! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf
!                                     |
!    Module for Ionosphere Model      |
!                                     |

module ModIonosphere

  use ModNumConst
  use IE_ModSize
  use ModUtilities, ONLY: CON_stop, CON_set_do_test

  implicit none
  save

  ! Ionosphere solution parameters
  real, parameter ::                            &
       IONO_TOLER = 5.0e-05,                    &
       IONO_MU = 1.256637e-06,                  &
       IONO_Theta_0 = 0.0001

  ! This served as a map between iModel and what the number means.
  ! Should be deleted once we're done here.
!  integer, parameter :: IONO_Model_No_Hall = 1, &
!       IONO_Model_With_Hall = 2,                &
!       IONO_Model_With_Simple_Aurora = 3,       &
!       IONO_Model_With_Complex_Aurora = 4

  real :: IONO_Bdp,   &
       IONO_Radius=1.0, IONO_Height=1.0, Radius

  real :: cpcp_north=0.0, cpcp_south=0.0

  logical :: DoUseGMPe = .false., DoUseGMPpar = .false., &
          DoUseGMPepar = .false.

  ! Ionosphere Solution on the whole grid
  real, allocatable :: IONO_Phi(:,:)
  real, allocatable :: IONO_IonNumFlux(:,:)
  real, allocatable :: IONO_Joule(:,:)
  real, allocatable :: IONO_Jr(:,:)
  real, allocatable :: IONO_Ave_E(:,:)
  real, allocatable :: IONO_Eflux(:,:)
  real, allocatable :: IONO_SigmaP(:,:)
  real, allocatable :: IONO_SigmaH(:,:)

  ! Ionosphere solution array definitions
  real, allocatable :: IONO_NORTH_Phi(:,:)
  real, allocatable :: IONO_SOUTH_Phi(:,:)
  real, allocatable :: IONO_NORTH_X(:,:)
  real, allocatable :: IONO_NORTH_Y(:,:)
  real, allocatable :: IONO_NORTH_Z(:,:)
  real, allocatable :: IONO_NORTH_Theta(:,:)
  real, allocatable :: IONO_NORTH_Psi(:,:)
  real, allocatable :: IONO_SOUTH_X(:,:)
  real, allocatable :: IONO_SOUTH_Y(:,:)
  real, allocatable :: IONO_SOUTH_Z(:,:)
  real, allocatable :: IONO_SOUTH_Theta(:,:)
  real, allocatable :: IONO_SOUTH_Psi(:,:)
  real, allocatable :: IONO_NORTH_Ex(:,:)
  real, allocatable :: IONO_NORTH_Ey(:,:)
  real, allocatable :: IONO_NORTH_Ez(:,:)
  real, allocatable :: IONO_NORTH_ETh(:,:)
  real, allocatable :: IONO_NORTH_EPs(:,:)
  real, allocatable :: IONO_SOUTH_Ex(:,:)
  real, allocatable :: IONO_SOUTH_Ey(:,:)
  real, allocatable :: IONO_SOUTH_Ez(:,:)
  real, allocatable :: IONO_SOUTH_ETh(:,:)
  real, allocatable :: IONO_SOUTH_EPs(:,:)
  real, allocatable :: IONO_NORTH_Ux(:,:)
  real, allocatable :: IONO_NORTH_Uy(:,:)
  real, allocatable :: IONO_NORTH_Uz(:,:)
  real, allocatable :: IONO_NORTH_UTh(:,:)
  real, allocatable :: IONO_NORTH_UPs(:,:)
  real, allocatable :: IONO_SOUTH_Ux(:,:)
  real, allocatable :: IONO_SOUTH_Uy(:,:)
  real, allocatable :: IONO_SOUTH_Uz(:,:)
  real, allocatable :: IONO_SOUTH_UTh(:,:)
  real, allocatable :: IONO_SOUTH_UPs(:,:)
  real, allocatable :: IONO_NORTH_EFlux(:,:)
  real, allocatable :: IONO_NORTH_Ave_E(:,:)
  real, allocatable :: IONO_SOUTH_EFlux(:,:)
  real, allocatable :: IONO_SOUTH_Ave_E(:,:)
  real, allocatable :: IONO_NORTH_Sigma0(:,:)
  real, allocatable :: IONO_NORTH_SigmaH(:,:)
  real, allocatable :: IONO_NORTH_SigmaP(:,:)
  real, allocatable :: IONO_NORTH_SigmaThTh(:,:)
  real, allocatable :: IONO_NORTH_SigmaThPs(:,:)
  real, allocatable :: IONO_NORTH_SigmaPsPs(:,:)
  real, allocatable :: IONO_SOUTH_Sigma0(:,:)
  real, allocatable :: IONO_SOUTH_SigmaH(:,:)
  real, allocatable :: IONO_SOUTH_SigmaP(:,:)
  real, allocatable :: IONO_SOUTH_SigmaThTh(:,:)
  real, allocatable :: IONO_SOUTH_SigmaThPs(:,:)
  real, allocatable :: IONO_SOUTH_SigmaPsPs(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThTh_dTheta(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThPs_dTheta(:,:)
  real, allocatable :: IONO_NORTH_dSigmaPsPs_dTheta(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThTh_dPsi(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThPs_dPsi(:,:)
  real, allocatable :: IONO_NORTH_dSigmaPsPs_dPsi(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThTh_dTheta(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThPs_dTheta(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaPsPs_dTheta(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThTh_dPsi(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThPs_dPsi(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaPsPs_dPsi(:,:)
  real, allocatable :: IONO_NORTH_Joule(:,:)
  real, allocatable :: IONO_SOUTH_Joule(:,:)
  real, allocatable :: IONO_NORTH_IonNumFlux(:,:)
  real, allocatable :: IONO_SOUTH_IonNumFlux(:,:)
  real, allocatable :: IONO_NORTH_JR(:,:)
  real, allocatable :: IONO_NORTH_JTh(:,:)
  real, allocatable :: IONO_NORTH_JPs(:,:)
  real, allocatable :: IONO_NORTH_Jx(:,:)
  real, allocatable :: IONO_NORTH_Jy(:,:)
  real, allocatable :: IONO_NORTH_Jz(:,:)
  real, allocatable :: IONO_SOUTH_JR(:,:)
  real, allocatable :: IONO_SOUTH_JTh(:,:)
  real, allocatable :: IONO_SOUTH_JPs(:,:)
  real, allocatable :: IONO_SOUTH_Jx(:,:)
  real, allocatable :: IONO_SOUTH_Jy(:,:)
  real, allocatable :: IONO_SOUTH_Jz(:,:)
  real, allocatable :: IONO_NORTH_TGCM_JR(:,:)
  real, allocatable :: IONO_SOUTH_TGCM_JR(:,:)
  real, allocatable :: IONO_NORTH_Fake_JR(:,:)
  real, allocatable :: IONO_SOUTH_Fake_JR(:,:)
  real, allocatable :: iono_north_im_jr(:,:)
  real, allocatable :: iono_south_im_jr(:,:)
  real, allocatable :: iono_north_im_avee(:,:)
  real, allocatable :: iono_south_im_avee(:,:)
  real, allocatable :: iono_north_im_eflux(:,:)
  real, allocatable :: iono_south_im_eflux(:,:)
  ! Sources of Conductances
  real, allocatable :: IONO_NORTH_DIFFI_Ave_E(:,:)
  real, allocatable :: IONO_SOUTH_DIFFI_Ave_E(:,:)
  real, allocatable :: IONO_NORTH_DIFFE_Ave_E(:,:)
  real, allocatable :: IONO_SOUTH_DIFFE_Ave_E(:,:)
  real, allocatable :: IONO_NORTH_MONO_Ave_E(:,:)
  real, allocatable :: IONO_SOUTH_MONO_Ave_E(:,:)
  real, allocatable :: IONO_NORTH_BBND_Ave_E(:,:)
  real, allocatable :: IONO_SOUTH_BBND_Ave_E(:,:)

  real, allocatable :: IONO_NORTH_DIFFI_EFlux(:,:)
  real, allocatable :: IONO_SOUTH_DIFFI_EFlux(:,:)
  real, allocatable :: IONO_NORTH_DIFFE_EFlux(:,:)
  real, allocatable :: IONO_SOUTH_DIFFE_EFlux(:,:)
  real, allocatable :: IONO_NORTH_MONO_EFlux(:,:)
  real, allocatable :: IONO_SOUTH_MONO_EFlux(:,:)
  real, allocatable :: IONO_NORTH_BBND_EFlux(:,:)
  real, allocatable :: IONO_SOUTH_BBND_EFlux(:,:)

  logical, allocatable :: IsFilledWithIm(:,:)

  real, allocatable :: IONO_NORTH_invB(:,:)
  real, allocatable :: IONO_SOUTH_invB(:,:)
  real, allocatable :: IONO_NORTH_rho(:,:)
  real, allocatable :: IONO_SOUTH_rho(:,:)
  real, allocatable :: IONO_NORTH_p(:,:)
  real, allocatable :: IONO_SOUTH_p(:,:)
  real, allocatable :: IONO_NORTH_Ppar(:,:)
  real, allocatable :: IONO_SOUTH_Ppar(:,:)
  real, allocatable :: IONO_NORTH_Pe(:,:)
  real, allocatable :: IONO_SOUTH_Pe(:,:)
  real, allocatable :: IONO_NORTH_Pepar(:,:)
  real, allocatable :: IONO_SOUTH_Pepar(:,:)
  real, allocatable :: IONO_NORTH_t(:,:)
  real, allocatable :: IONO_SOUTH_t(:,:)
  real, allocatable :: IONO_NORTH_dLat(:,:)
  real, allocatable :: IONO_SOUTH_dLat(:,:)
  real, allocatable :: IONO_NORTH_dLon(:,:)
  real, allocatable :: IONO_SOUTH_dLon(:,:)

  ! Save (X,Y,Z) in different coordinate systems.
  real, allocatable :: IONO_NORTH_GEO_XyzD(:, :, :)
  real, allocatable :: IONO_SOUTH_GEO_XyzD(:, :, :)
  real, allocatable :: IONO_NORTH_GSE_XyzD(:, :, :)
  real, allocatable :: IONO_SOUTH_GSE_XyzD(:, :, :)

  ! Pentadiagonal matrix for the Poisson equation
  real, allocatable :: C_A(:,:)
  real, allocatable :: C_B(:,:)
  real, allocatable :: C_C(:,:)
  real, allocatable :: C_D(:,:)
  real, allocatable :: C_E(:,:)
  real, dimension(:), allocatable :: d_I, e_I, f_I, e1_I, f1_I
  logical :: north, DoPrecond
  integer :: nThetaUsed, nX

  real, dimension(IONO_nTheta) :: dTheta_North, dTheta_South
  real, dimension(IONO_nPsi)   :: dPsi_North, dPsi_South

contains
  !============================================================================

  subroutine init_mod_ionosphere

    !--------------------------------------------------------------------------
    if(allocated(IONO_Phi)) RETURN

    ! Initialize these global grid arrays to 0 (for output before solve)
    allocate(IONO_Phi(2*IONO_nTheta-1,IONO_nPsi));        IONO_Phi = 0
    allocate(IONO_IonNumFlux(2*IONO_nTheta-1,IONO_nPsi)); IONO_IonNumFlux = 0
    allocate(IONO_Joule(2*IONO_nTheta-1,IONO_nPsi));      IONO_Joule = 0
    allocate(IONO_Jr(2*IONO_nTheta-1,IONO_nPsi));         IONO_Jr = 0
    allocate(IONO_Ave_E(2*IONO_nTheta-1,IONO_nPsi));      IONO_Ave_E = 0
    allocate(IONO_Eflux(2*IONO_nTheta-1,IONO_nPsi));      IONO_Eflux = 0
    allocate(IONO_SigmaP(2*IONO_nTheta-1,IONO_nPsi));     IONO_SigmaP = 0
    allocate(IONO_SigmaH(2*IONO_nTheta-1,IONO_nPsi));     IONO_SigmaH = 0

    allocate(IONO_NORTH_PHI(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_PHI(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_X(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Y(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Z(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Theta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Psi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_X(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Y(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Z(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Theta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Psi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ex(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ey(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ez(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_ETh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_EPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ex(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ey(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ez(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_ETh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_EPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Uy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Uz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_UTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_UPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Uy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Uz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_UTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_UPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Sigma0(IONO_nTheta,IONO_nPsi)); IONO_NORTH_Sigma0=1000.
    allocate(IONO_NORTH_SigmaH(IONO_nTheta,IONO_nPsi)); IONO_NORTH_SigmaH = 0
    allocate(IONO_NORTH_SigmaP(IONO_nTheta,IONO_nPsi)); IONO_NORTH_SigmaP = 0
    allocate(IONO_NORTH_SigmaThTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_SigmaThPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_SigmaPsPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Sigma0(IONO_nTheta,IONO_nPsi)); IONO_SOUTH_Sigma0=1000.
    allocate(IONO_SOUTH_SigmaH(IONO_nTheta,IONO_nPsi)); IONO_SOUTH_SigmaH = 0
    allocate(IONO_SOUTH_SigmaP(IONO_nTheta,IONO_nPsi)); IONO_SOUTH_SigmaP = 0
    allocate(IONO_SOUTH_SigmaThTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_SigmaThPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_SigmaPsPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThTh_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaPsPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThTh_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaPsPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThTh_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaPsPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThTh_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaPsPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Joule(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Joule(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_IonNumFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_IonNumFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_JTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_JPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Jx(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Jy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Jz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_JTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_JPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Jx(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Jy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Jz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_TGCM_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_TGCM_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Fake_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Fake_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_north_im_jr(IONO_nTheta,IONO_nPsi))
    allocate(IONO_south_im_jr(IONO_nTheta,IONO_nPsi))
    allocate(IONO_north_im_avee(IONO_nTheta,IONO_nPsi))
    allocate(IONO_south_im_avee(IONO_nTheta,IONO_nPsi))
    allocate(IONO_north_im_eflux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_south_im_eflux(IONO_nTheta,IONO_nPsi))
    allocate(IsFilledWithIm(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_invB(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_invB(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_rho(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_rho(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_p(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_p(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ppar(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ppar(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Pe(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Pe(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Pepar(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Pepar(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_t(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_t(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dLat(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dLat(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dLon(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dLon(IONO_nTheta,IONO_nPsi))
    allocate(C_A(IONO_nTheta,IONO_nPsi))
    allocate(C_B(IONO_nTheta,IONO_nPsi))
    allocate(C_C(IONO_nTheta,IONO_nPsi))
    allocate(C_D(IONO_nTheta,IONO_nPsi))
    allocate(C_E(IONO_nTheta,IONO_nPsi))
    ! Sources of Conductances
    allocate(IONO_NORTH_DIFFI_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_DIFFI_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_DIFFE_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_DIFFE_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_MONO_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_MONO_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_BBND_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_BBND_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_DIFFI_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_DIFFI_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_DIFFE_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_DIFFE_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_MONO_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_MONO_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_BBND_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_BBND_EFlux(IONO_nTheta,IONO_nPsi))

    IONO_NORTH_DIFFI_Ave_E = 0.; IONO_NORTH_DIFFI_EFlux = 0.
    IONO_SOUTH_DIFFI_Ave_E = 0.; IONO_SOUTH_DIFFI_EFlux = 0.
    IONO_NORTH_DIFFE_Ave_E = 0.; IONO_NORTH_DIFFE_EFlux = 0.
    IONO_SOUTH_DIFFE_Ave_E = 0.; IONO_SOUTH_DIFFE_EFlux = 0.
    IONO_NORTH_MONO_Ave_E = 0.; IONO_NORTH_MONO_EFlux = 0.
    IONO_SOUTH_MONO_Ave_E = 0.; IONO_SOUTH_MONO_EFlux = 0.
    IONO_NORTH_BBND_Ave_E = 0.; IONO_NORTH_BBND_EFlux = 0.
    IONO_SOUTH_BBND_Ave_E = 0.; IONO_SOUTH_BBND_EFlux = 0.

    allocate(IONO_NORTH_GEO_XyzD(3, IONO_nTheta, IONO_nPsi))
    allocate(IONO_NORTH_GSE_XyzD(3, IONO_nTheta, IONO_nPsi))
    allocate(IONO_SOUTH_GEO_XyzD(3, IONO_nTheta, IONO_nPsi))
    allocate(IONO_SOUTH_GSE_XyzD(3, IONO_nTheta, IONO_nPsi))

    IONO_NORTH_GEO_XyzD = 0.; IONO_NORTH_GSE_XyzD = 0
    IONO_SOUTH_GEO_XyzD = 0.; IONO_SOUTH_GSE_XyzD = 0

  end subroutine init_mod_ionosphere
  !============================================================================

  subroutine clean_mod_ionosphere

    !--------------------------------------------------------------------------
    if(.not.allocated(IONO_Phi)) RETURN

    deallocate(IONO_Phi)
    deallocate(IONO_IonNumFlux)
    deallocate(IONO_Joule)
    deallocate(IONO_Jr)
    deallocate(IONO_Ave_E)
    deallocate(IONO_Eflux)
    deallocate(IONO_SigmaP)
    deallocate(IONO_SigmaH)
    deallocate(IONO_NORTH_PHI)
    deallocate(IONO_SOUTH_PHI)
    deallocate(IONO_NORTH_X)
    deallocate(IONO_NORTH_Y)
    deallocate(IONO_NORTH_Z)
    deallocate(IONO_NORTH_Theta)
    deallocate(IONO_NORTH_Psi)
    deallocate(IONO_SOUTH_X)
    deallocate(IONO_SOUTH_Y)
    deallocate(IONO_SOUTH_Z)
    deallocate(IONO_SOUTH_Theta)
    deallocate(IONO_SOUTH_Psi)
    deallocate(IONO_NORTH_Ex)
    deallocate(IONO_NORTH_Ey)
    deallocate(IONO_NORTH_Ez)
    deallocate(IONO_NORTH_ETh)
    deallocate(IONO_NORTH_EPs)
    deallocate(IONO_SOUTH_Ex)
    deallocate(IONO_SOUTH_Ey)
    deallocate(IONO_SOUTH_Ez)
    deallocate(IONO_SOUTH_ETh)
    deallocate(IONO_SOUTH_EPs)
    deallocate(IONO_NORTH_Ux)
    deallocate(IONO_NORTH_Uy)
    deallocate(IONO_NORTH_Uz)
    deallocate(IONO_NORTH_UTh)
    deallocate(IONO_NORTH_UPs)
    deallocate(IONO_SOUTH_Ux)
    deallocate(IONO_SOUTH_Uy)
    deallocate(IONO_SOUTH_Uz)
    deallocate(IONO_SOUTH_UTh)
    deallocate(IONO_SOUTH_UPs)
    deallocate(IONO_NORTH_EFlux)
    deallocate(IONO_NORTH_Ave_E)
    deallocate(IONO_SOUTH_EFlux)
    deallocate(IONO_SOUTH_Ave_E)
    deallocate(IONO_NORTH_Sigma0)
    deallocate(IONO_NORTH_SigmaH)
    deallocate(IONO_NORTH_SigmaP)
    deallocate(IONO_NORTH_SigmaThTh)
    deallocate(IONO_NORTH_SigmaThPs)
    deallocate(IONO_NORTH_SigmaPsPs)
    deallocate(IONO_SOUTH_Sigma0)
    deallocate(IONO_SOUTH_SigmaH)
    deallocate(IONO_SOUTH_SigmaP)
    deallocate(IONO_SOUTH_SigmaThTh)
    deallocate(IONO_SOUTH_SigmaThPs)
    deallocate(IONO_SOUTH_SigmaPsPs)
    deallocate(IONO_NORTH_dSigmaThTh_dTheta)
    deallocate(IONO_NORTH_dSigmaThPs_dTheta)
    deallocate(IONO_NORTH_dSigmaPsPs_dTheta)
    deallocate(IONO_NORTH_dSigmaThTh_dPsi)
    deallocate(IONO_NORTH_dSigmaThPs_dPsi)
    deallocate(IONO_NORTH_dSigmaPsPs_dPsi)
    deallocate(IONO_SOUTH_dSigmaThTh_dTheta)
    deallocate(IONO_SOUTH_dSigmaThPs_dTheta)
    deallocate(IONO_SOUTH_dSigmaPsPs_dTheta)
    deallocate(IONO_SOUTH_dSigmaThTh_dPsi)
    deallocate(IONO_SOUTH_dSigmaThPs_dPsi)
    deallocate(IONO_SOUTH_dSigmaPsPs_dPsi)
    deallocate(IONO_NORTH_Joule)
    deallocate(IONO_SOUTH_Joule)
    deallocate(IONO_NORTH_IonNumFlux)
    deallocate(IONO_SOUTH_IonNumFlux)
    deallocate(IONO_NORTH_JR)
    deallocate(IONO_NORTH_JTh)
    deallocate(IONO_NORTH_JPs)
    deallocate(IONO_NORTH_Jx)
    deallocate(IONO_NORTH_Jy)
    deallocate(IONO_NORTH_Jz)
    deallocate(IONO_SOUTH_JR)
    deallocate(IONO_SOUTH_JTh)
    deallocate(IONO_SOUTH_JPs)
    deallocate(IONO_SOUTH_Jx)
    deallocate(IONO_SOUTH_Jy)
    deallocate(IONO_SOUTH_Jz)
    deallocate(IONO_NORTH_TGCM_JR)
    deallocate(IONO_SOUTH_TGCM_JR)
    deallocate(IONO_NORTH_Fake_JR)
    deallocate(IONO_SOUTH_Fake_JR)
    deallocate(IONO_north_im_jr)
    deallocate(IONO_south_im_jr)
    deallocate(IONO_north_im_avee)
    deallocate(IONO_south_im_avee)
    deallocate(IONO_north_im_eflux)
    deallocate(IONO_south_im_eflux)
    deallocate(IsFilledWithIm)
    deallocate(IONO_NORTH_invB)
    deallocate(IONO_SOUTH_invB)
    deallocate(IONO_NORTH_rho)
    deallocate(IONO_SOUTH_rho)
    deallocate(IONO_NORTH_p)
    deallocate(IONO_SOUTH_p)
    deallocate(IONO_NORTH_Ppar)
    deallocate(IONO_SOUTH_Ppar)
    deallocate(IONO_NORTH_Pe)
    deallocate(IONO_SOUTH_Pe)
    deallocate(IONO_NORTH_Pepar)
    deallocate(IONO_SOUTH_Pepar)
    deallocate(IONO_NORTH_t)
    deallocate(IONO_SOUTH_t)
    deallocate(IONO_NORTH_dLat)
    deallocate(IONO_SOUTH_dLat)
    deallocate(IONO_NORTH_dLon)
    deallocate(IONO_SOUTH_dLon)
    deallocate(C_A)
    deallocate(C_B)
    deallocate(C_C)
    deallocate(C_D)
    deallocate(C_E)

    deallocate(IONO_NORTH_GEO_XyzD)
    deallocate(IONO_NORTH_GSE_XyzD)
    deallocate(IONO_SOUTH_GEO_XyzD)
    deallocate(IONO_SOUTH_GSE_XyzD)

    ! Sources of Conductances
    deallocate(IONO_NORTH_DIFFI_Ave_E)
    deallocate(IONO_SOUTH_DIFFI_Ave_E)
    deallocate(IONO_NORTH_DIFFE_Ave_E)
    deallocate(IONO_SOUTH_DIFFE_Ave_E)
    deallocate(IONO_NORTH_MONO_Ave_E)
    deallocate(IONO_SOUTH_MONO_Ave_E)
    deallocate(IONO_NORTH_BBND_Ave_E)
    deallocate(IONO_SOUTH_BBND_Ave_E)
    deallocate(IONO_NORTH_DIFFI_EFlux)
    deallocate(IONO_SOUTH_DIFFI_EFlux)
    deallocate(IONO_NORTH_DIFFE_EFlux)
    deallocate(IONO_SOUTH_DIFFE_Eflux)
    deallocate(IONO_NORTH_MONO_EFlux)
    deallocate(IONO_SOUTH_MONO_EFlux)
    deallocate(IONO_NORTH_BBND_EFlux)
    deallocate(IONO_SOUTH_BBND_EFlux)

  end subroutine clean_mod_ionosphere
  !============================================================================

end module ModIonosphere
!==============================================================================
