! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

module ModMagnit

  use ModConst, ONLY: cBoltzmann, cProtonMass, cElectronMass, cKToKEV, cKEVToK
  use ModUtilities, ONLY: CON_stop, CON_set_do_test
  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi

  implicit none
  save

  ! Unit conversion factors and other constants:
  real, parameter :: cPtoKProt = cKToKEV/cBoltzmann
  real, parameter :: cPtoKElec = cElectronMass/cBoltzmann * cKToKEV

  ! Configuration parameters:
  logical :: DoUseGmPe=.false.  ! Use electron pressure?

  ! Diffuse auroral parameters
  real :: ratioPe = 1./6. ! Ratio of electron P to proton P.

  ! Loss cone factors set flux in loss cones for proton & electron diffuse,
  ! monoenergetic, and broadband precip. Factors are for Energy Flux and
  ! Number Flux, each.
  real :: ConeEfluxDife = 0.217, ConeEfluxDifp = 0.207, &
          ConeEfluxMono = 1.000, ConeEfluxBbnd = 2.247, &
          ConeNfluxDife = 0.055, ConeNfluxDifp = 0.038, &
          ConeNfluxMono = 0.741, ConeNfluxBbnd = 0.494

!  MinWidth = 5.0 * cPi / 180.0
!  nHalfSmooth = 5
!     MulFac_Dae = 1.0e22
!     MulFac_Def = 5.0e19
!     MulFac_ef = 0.2e7
!     MulFac_ae = 1.0 / 1.0e11
!  MulFac_Dae = 1.0e22
!  MulFac_Def = 1.0e19
!  MulFac_ef = 0.3e6! * 1e03
!  MulFac_ae = 3.65e-12

  ! reminder: GmRhoFloor = 1E-21, GmPFloor = 1E-13

  contains
  !============================================================================
  subroutine print_magnit_config(iUnitIn)

    integer, intent(in) :: iUnitIn
    !--------------------------------------------------------------------------

    write(iUnitin,'(a)') "MAGNIT Physics-based Aurora"
    write(iUnitin,'(a)') "(Beta-testing Phase)"
    write(iUnitIn,*) '################## CAUTION ##################'
    write(iUnitIn,*) 'MAGNIT is in BETA TESTING and may not be stable.'
    write(iUnitIn,*) 'Proceed with caution. For issues, contact developers.'
    write(iUnitIn,*) '#############################################'
    if(DoUseGmPe)then
      write(iUnitIn,*) "Electron precip obtained directly from Pe."
    else
      write(iUnitIn,*) "Electron precip obtained from LT-flipped protons"
      write(iUnitIn,'(a, f5.3)') "Pe/P ratio set to ", ratioPe
    end if

  end subroutine print_magnit_config
  !============================================================================

  subroutine magnit_gen_fluxes(NameHemiIn, &
      AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
      EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II)

    use ModConst, ONLY: cPi, cKEV
    use ModIonosphere, ONLY: IONO_North_p, IONO_North_rho, &
        IONO_South_p, IONO_South_rho

    ! Given magnetospheric density, pressure, and FACs, calculate diffuse and
    ! monoenergetic precipitating fluxes.

    character(len=*), intent(in) :: NameHemiIn

    ! Set arrays to hold precip values. Magnetospheric pressure and density
    ! from GM (SI units), Average Energy and Energy Flux (units of KeV and
    ! mW/m2 aka ergs/cm2, respectively) for each precip type.
    real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
        AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
        EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

    ! Set arrays to hold magnetospheric values.
    real, dimension(IONO_nTheta, IONO_nPsi) :: &
        MagP_II, MagNp_II, MagPe_II, MagNe_II, MagRho_II, &
        NfluxDiffe_II, NfluxDiffi_II

    integer :: j

    ! Debug variables:
    logical :: DoTestMe, DoTest
    character(len=*), parameter:: NameSub = 'magnit_gen_fluxes'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe) &
         write(*,*)'IE: '//NameSub//' called for hemisphere='//NameHemiIn

    ! Set values based on hemisphere. Calculate number density.
    if(NameHemiIn == 'north')then
        MagP_II = iono_north_p
        MagNp_II = iono_north_rho / cProtonMass
    else if (NameHemiIn == 'south')then
        MagP_II = iono_south_p
        MagRho_II = iono_south_rho / cProtonMass
    else
      call CON_stop(NameSub//' : unrecognized hemisphere - '//NameHemiIn)
    end if

    ! Set default/background values.

    ! If we need OCFLB calculation, we would put it here.

    ! If not using explicit Pe from GM, obtain values from proton values.
    ! flip across noon midnight, etc. to fill MagPe, MagRhoe
    ! Scale Pe using ratioPe
    if(.not. DoUseGmPe) then
      do j=1, Iono_nPsi
        MagPe_II(:, j) = ratioPe * MagP_II(:, IONO_nPsi-j+1)
        MagNe_II(:, j) = MagNp_II(:, IONO_nPsi-j+1)
      end do
    end if

    ! Calculate diffuse precipitation: protons.
    AvgEDiffi_II  = MagP_II / (MagNp_II * cKEV)  ! T = P/nk in keV
    NfluxDiffi_II = ConeNfluxDifp * MagNp_II * (cKEV * AvgEDiffi_II)**0.5 / &
                    sqrt(2 * cPi * cProtonMass)  ! units of mW/m2 or ergs/cm2
    ! units of mW/m2 or ergs/cm2
    EfluxDiffi_II = 2000 * ConeEfluxDifp * MagNp_II * &
                    (cKEV * AvgEDiffi_II)**1.5 / sqrt(2 * cPi * cProtonMass)
    ! Recalc to make consistent with ConeFactors
    AvgEDiffi_II = EfluxDiffi_II / NfluxDiffi_II

    ! Calculate diffuse precipitation: electrons.
    AvgEDiffe_II  = MagPe_II / (MagNe_II * cKEV)  ! T = P/nk in keV
    NfluxDiffe_II = ConeNfluxDife * MagNe_II * (cKEV * AvgEDiffe_II)**0.5 / &
                 sqrt(2 * cPi * cElectronMass)  ! units of mW/m2 or ergs/cm2
    ! units of mW/m2 or ergs/cm2
    EfluxDiffe_II = 2000 * ConeEfluxDife * MagNe_II * &
                    (cKEV * AvgEDiffe_II)**1.5 / sqrt(2 * cPi * cElectronMass)
    ! Recalc to make consistent with ConeFactors
    AvgEDiffe_II = EfluxDiffe_II / NfluxDiffe_II

  end subroutine magnit_gen_fluxes
  !============================================================================

end module ModMagnit!==============================================================================
