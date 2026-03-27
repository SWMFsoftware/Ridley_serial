! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

module ModImp

  ! Inner Magnetosphere Precipitation model

  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi

  use ModUtilities, ONLY: CON_stop, CON_set_do_test

  implicit none
  save

contains
  !============================================================================
  subroutine imp_gen_fluxes(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II, &
       AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II,&
       EfluxBbnd_II, LatIn_II)

    use ModIonosphere, ONLY: DoUseIMSpectrum
    use ModMagnit, ONLY: monoenergetic_flux, broadband_flux, smooth_polar_cap
    use ModIonosphere, ONLY: IONO_NORTH_JR, IONO_SOUTH_JR, &
         IONO_NORTH_invB, IONO_SOUTH_invB, IONO_NORTH_Poynting, &
         IONO_SOUTH_Poynting, IONO_north_im_boundary, &
         IONO_south_im_boundary, DoPolarCapSmoothing, &
         iono_north_im_aveeElec, iono_south_im_aveeElec, &
         iono_north_im_efluxElec, iono_south_im_eFluxElec, &
         iono_north_im_aveeHydr, iono_south_im_aveeHydr, &
         iono_north_im_efluxHydr, iono_south_im_eFluxHydr, &
         iono_north_im_nElecPrec, iono_south_im_nElecPrec
    use ModConst, ONLY: cKEV

    real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
         AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
         EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

    real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: LatIn_II

    real, dimension(IONO_nTheta, IONO_nPsi) :: &
         FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, Potential_II, &
         Poynting_II, ImBoundary_II

    character(len=*), intent(in) :: NameHemiIn

    character(len=*), parameter:: NameSub = 'imp_gen_fluxes'
    !--------------------------------------------------------------------------
    if (trim(NameHemiIn) == 'south') then
      AvgEDiffe_II = iono_south_im_aveeElec / 1000.0 ! eV to keV
      EfluxDiffe_II = iono_south_im_efluxElec / 1000.0 ! mW/m^2 to W/m^2
      AvgEDiffi_II = iono_south_im_aveeHydr / 1000.0
      EfluxDiffi_II = iono_south_im_efluxHydr / 1000.0
    else if (trim(NameHemiIn) == 'north') then
      AvgEDiffe_II = iono_north_im_aveeElec / 1000.0
      EfluxDiffe_II = iono_north_im_efluxElec / 1000.0
      AvgEDiffi_II = iono_north_im_aveeHydr / 1000.0
      EfluxDiffi_II = iono_north_im_efluxHydr / 1000.0
    else
       call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                NameHemiIn)
    end if

    where(Poynting_II < 0) Poynting_II = 0

    ! Limits AvgE (primarily for low latitudes)
    ! I don't think this is necessary any longer
    where(AvgEMono_II > 100) AvgEMono_II = 100
    where(AvgEDiffe_II > 100)AvgEDiffe_II = 100
    where(AvgEDiffi_II > 200)AvgEDiffi_II= 200

    if(NameHemiIn == 'north')then
      FAC_II = IONO_NORTH_JR
      OCFL_II = IONO_NORTH_invB
      Poynting_II = -IONO_NORTH_Poynting
      ImBoundary_II = IONO_north_im_boundary
    else if (NameHemiIn == 'south')then
      FAC_II = IONO_SOUTH_JR
      OCFL_II = IONO_SOUTH_invB
      Poynting_II = -IONO_SOUTH_Poynting
      ImBoundary_II = IONO_south_im_boundary
    end if

    ! Smooth area between closed and open field lines
    if (DoPolarCapSmoothing) then
      call smooth_polar_cap(AvgEDiffe_II, OCFL_II, NameHemiIn, ImBoundary_II)
      call smooth_polar_cap(AvgEDiffi_II, OCFL_II, NameHemiIn, ImBoundary_II)
      call smooth_polar_cap(EfluxDiffe_II, OCFL_II, NameHemiIn, ImBoundary_II)
      call smooth_polar_cap(EfluxDiffi_II, OCFL_II, NameHemiIn, ImBoundary_II)
    end if

    ! Calculate monoenergetic flux (same as MAGNIT)
    NfluxDiffe_II = EfluxDiffe_II / AvgEDiffe_II / cKEV
    ! Limit for stability, should probably be removed eventually
    ElectronTemp_II = 2.0 * AvgEDiffe_II * cKEV ! kEV to J
    call monoenergetic_flux(FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, &
            AvgEDiffe_II, LatIn_II, EfluxMono_II, AvgEMono_II)

    call broadband_flux(Poynting_II, EfluxBbnd_II, AvgEBbnd_II)

  end subroutine imp_gen_fluxes
  !============================================================================
end module ModImp
!==============================================================================
