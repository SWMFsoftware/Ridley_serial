! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

module ModMagnit

  use ModConst, ONLY: cBoltzmann, cProtonMass, cElectronMass, cElectronCharge, &
          cKToKEV, cKEVToK
  use ModUtilities, ONLY: CON_stop, CON_set_do_test
  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi, DoUseGMPe, DoUseGMPpar, &
          DoUseGMPepar

  implicit none
  save

  ! Unit conversion factors and other constants:
  real, parameter :: cPtoKProt = cKToKEV/cBoltzmann
  real, parameter :: cPtoKElec = cElectronMass/cBoltzmann * cKToKEV
  real, parameter :: cFACFloor = 1.0E-12 ! A/m2

          ! Configuration parameters:
  ! Diffuse auroral parameters
  real :: ratioPe = 1./6. ! Ratio of electron P to proton P.

  ! Loss cone factors set flux in loss cones for proton & electron diffuse,
  ! monoenergetic, and broadband precip. Factors are for Energy Flux and
  ! Number Flux, each. See Mukhopadhyay et al. 2022 for details.
  real :: ConeEfluxDife = 0.217, ConeEfluxDifp = 0.207, &
          ConeEfluxMono = 1.000, ConeEfluxBbnd = 2.247, &
          ConeNfluxDife = 0.055, ConeNfluxDifp = 0.038, &
          ConeNfluxMono = 0.741, ConeNfluxBbnd = 0.494
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
    if (DoUseGMPpar) write(iUnitIn,*) "Using anisotropic ion pressure."
    if (DoUseGMPepar) write(iUnitIn,*) "Using anisotropic electron pressure."

  end subroutine print_magnit_config
  !============================================================================

  subroutine magnit_gen_fluxes(NameHemiIn, &
      AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
      EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II, &
      LatIn_II)

    use ModConst, ONLY: cPi, cKEV
    use ModIonosphere, ONLY: IONO_North_p, IONO_North_rho, &
        IONO_South_p, IONO_South_rho, IONO_NORTH_JR, IONO_SOUTH_JR, &
        IONO_NORTH_invB, IONO_SOUTH_invB, IONO_NORTH_Joule, &
        IONO_SOUTH_Joule, IONO_North_Pe, IONO_South_Pe
    use ModPlanetConst, ONLY: rPlanet_I, IonoHeightPlanet_I, Earth_

    ! Given magnetospheric density, pressure, and FACs, calculate diffuse and
    ! monoenergetic precipitating fluxes.
    ! Resulting units are W/m2 and KeV for energy flux
    ! and average energy, respectively.

    character(len=*), intent(in) :: NameHemiIn

    ! Set arrays to hold precip values. Magnetospheric pressure and density
    ! from GM (SI units), Average Energy and Energy Flux (units of KeV and
    ! W/m2, respectively) for each precip type.
    real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
        AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
        EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

    ! Import Hemispheric Latitudes for Magnetic Field Calculations
    real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: LatIn_II

    ! Set arrays to hold magnetospheric values.
    real, dimension(IONO_nTheta, IONO_nPsi) :: &
        MagP_II, MagNp_II, MagPe_II, MagNe_II, NfluxDiffe_II, NfluxDiffi_II, &
        NfluxMono_II, NfluxBbnd_II, OCFL_II, PrecipRatio_II=0, Potential_II=0, &
        MirrorRatio_II=0, VExponent_II=0, PotentialTerm_II=0, FAC_II=0, &
        Joule_II=0, Poynting_II=0, ElectronTemp_II=0, OCFL_flip_II=0

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
        FAC_II = IONO_NORTH_JR
        OCFL_II = IONO_NORTH_invB
        Joule_II = IONO_NORTH_Joule
        if (DoUseGMPe) MagPe_II = iono_north_pe
    else if (NameHemiIn == 'south')then
        MagP_II = iono_south_p
        MagNp_II = iono_south_rho / cProtonMass
        FAC_II = IONO_SOUTH_JR
        OCFL_II = IONO_SOUTH_invB
        Joule_II = IONO_SOUTH_Joule
        if (DoUseGMPe) MagPe_II = iono_south_pe
    else
      call CON_stop(NameSub//' : unrecognized hemisphere - '//NameHemiIn)
    end if

    ! Set default/background values - TBD.
    ! If we need OCFLB calculation, we would put it here.

    ! If not using explicit Pe from GM, obtain values from proton values.
    ! Flip across noon midnight, etc. to fill MagPe, MagNe
    ! Scale Pe using ratioPe
    if(.not. DoUseGMPe) then
      do j=1, Iono_nPsi
        MagPe_II(:, j) = ratioPe * MagP_II(:, IONO_nPsi-j+1)
        MagNe_II(:, j) = MagNp_II(:, IONO_nPsi-j+1)
        OCFL_flip_II(:, j) = OCFL_II(:, IONO_nPsi-j+1)
      end do
      OCFL_II = OCFL_flip_II
    else
      MagNe_II = MagNp_II
    end if

    ! Calculate diffuse precipitation: protons.
    AvgEDiffi_II  = MagP_II / MagNp_II  ! Temp = P/nk in Joules
    NfluxDiffi_II = ConeNfluxDifp * MagNp_II * AvgEDiffi_II**0.5 / &
                    sqrt(2 * cPi * cProtonMass)  ! units of #/m2/s
    ! units of W/m2
    EfluxDiffi_II = 2 * NfluxDiffi_II * AvgEDiffi_II * &
            ConeEfluxDifp/ConeNfluxDifp
    ! Recalc to make consistent with ConeFactors (and get units of keV)
    AvgEDiffi_II = EfluxDiffi_II / (NfluxDiffi_II * cKEV)

! Calculate diffuse precipitation: electrons.
    ElectronTemp_II  = MagPe_II / MagNe_II  ! T = P/nk in Joules
    NfluxDiffe_II = ConeNfluxDife * MagNe_II * ElectronTemp_II**0.5 / &
                 sqrt(2 * cPi * cElectronMass)  ! units of #/m2/s
    ! units of W/m2
    EfluxDiffe_II = 2 * NfluxDiffe_II * ElectronTemp_II * &
            ConeEfluxDife/ConeNfluxDife
    ! Recalc to make consistent with ConeFactors (and get units of keV)
    AvgEDiffe_II = EfluxDiffe_II / (NfluxDiffe_II * cKEV)

    ! Calculate discrete electron precipitation
    ! Calculate Nflux from current
    Potential_II = 0
    NfluxMono_II = FAC_II / cElectronCharge
    ! Calculate Parallel Potential Drop
    ! Calculate ratio of diffuse to monoenergetic precipitation
    PrecipRatio_II = NfluxMono_II / NfluxDiffe_II

    ! Calculate ratio of ionospheric magnetic field to plasma sheet
    ! magnetic field
    MirrorRatio_II = (rPlanet_I(Earth_) / (sin(LatIn_II)**2 * &
            (rPlanet_I(Earth_) + IonoHeightPlanet_I(Earth_))))**3 * &
            sqrt(1 + 3*cos(LatIn_II)**2)

    ! Potential calculations only valid where 1 <= NumCoefficient <= MirrorRatio
    where(1 <= PrecipRatio_II .and. PrecipRatio_II < MirrorRatio_II .and. &
          OCFL_II > 0)
      ! Put it all together into potential
      Potential_II = ElectronTemp_II / cElectronCharge * (1 - MirrorRatio_II) &
              * LOG((MirrorRatio_II - PrecipRatio_II) / (MirrorRatio_II - 1))
    elsewhere
      NfluxMono_II = NfluxDiffe_II
    end where
    ! Split up large calculation into a few steps
    ! Calculate large potential exponent
    VExponent_II = EXP(-cElectronCharge * Potential_II / ((ElectronTemp_II) * &
            (MirrorRatio_II - 1)))
    ! Calculate product term
    PotentialTerm_II = ((1 - VExponent_II) * ConeEfluxMono / &
         (1 + ((1 - 1/MirrorRatio_II) * VExponent_II)) * cElectronCharge * &
         Potential_II) + AvgEDiffe_II * cKEV

    ! Plug into equation for EFlux
    EfluxMono_II = NfluxMono_II * PotentialTerm_II

    ! Calculate Avg E in keV
    AvgEMono_II = EfluxMono_II / (NfluxMono_II * cKEV)
    !
!    ! Calculate broadband electron precipitation
!    ! Note: Joule Heating is in SI Units = W/m2
     ! Need to figure out where this value came from,
     ! I don't like a hardcoded 17 Sig Fig constant
    EfluxBbnd_II = 0
    AvgEBbnd_II = 0
    where(Joule_II > 0)
      Poynting_II = Joule_II * 0.43395593979143521 ! in W/m2

      ! Using empirical relationships from Zhang et al. 2015
      EfluxBbnd_II = 2e-3 * (ConeEfluxBbnd * Poynting_II) ** 0.5
      NfluxBbnd_II = 3e13 * (ConeNfluxBbnd * Poynting_II) ** 0.47
      AvgEBbnd_II = EfluxBbnd_II / (NfluxBbnd_II * cKEV)
    end where

  end subroutine magnit_gen_fluxes
  !============================================================================

end module ModMagnit
!==============================================================================
