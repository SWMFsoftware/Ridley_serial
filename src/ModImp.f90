! Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission
! For more information, see http://csem.engin.umich.edu/tools/swmf

! Inner Magnetosphere Precipitation model

module ModImp

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
      use ModMagnit, ONLY: monoenergetic_flux
      use ModIonosphere, ONLY: IONO_NORTH_JR, IONO_SOUTH_JR, &
              IONO_NORTH_invB, IONO_SOUTH_invB
      use ModConst, ONLY: cKEV

      real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
              AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
              EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

      real, intent(in), dimension(IONO_nTheta, IONO_nPsi) :: LatIn_II

      real, dimension(IONO_nTheta, IONO_nPsi) :: &
              FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II

      character(len=*), intent(in) :: NameHemiIn

    character(len=*), parameter:: NameSub = 'imp_gen_fluxes'
    !--------------------------------------------------------------------------
      ! Call correct IMP model based on filled inputs
      if (DoUseIMSpectrum) then
          call imp_spectral_flux(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II, &
                  AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, &
                  EfluxMono_II, EfluxBbnd_II)
      else
          call imp_integrated_flux(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II, &
                  AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, &
                  EfluxMono_II, EfluxBbnd_II)
      end if

      ! Limits AvgE (primarily for low latitudes)
      where(AvgEMono_II > 100) AvgEMono_II = 100
      where(AvgEDiffe_II > 100)AvgEDiffe_II = 100
      where(AvgEDiffi_II > 200)AvgEDiffi_II= 200

      if(NameHemiIn == 'north')then
          FAC_II = IONO_NORTH_JR
          OCFL_II = IONO_NORTH_invB
      else if (NameHemiIn == 'south')then
          FAC_II = IONO_SOUTH_JR
          OCFL_II = IONO_SOUTH_invB
       end if

      ! Calculate monoenergetic flux (same as MAGNIT)
      NfluxDiffe_II = EfluxDiffe_II / AvgEDiffe_II / cKEV
      ElectronTemp_II = 2.0 * AvgEDiffe_II / cKEV ! kEV to J(????)
      call monoenergetic_flux(FAC_II, OCFL_II, NfluxDiffe_II, ElectronTemp_II, &
              AvgEDiffe_II, LatIn_II, EfluxMono_II, AvgEMono_II)

      contains
    !==========================================================================
      !-----------------------------------------------------------------------
      subroutine imp_integrated_flux(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II,&
              AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, &
              EfluxMono_II, EfluxBbnd_II)

          use ModIonosphere, ONLY: &
                iono_north_im_aveeElec, iono_south_im_aveeElec, &
                iono_north_im_efluxElec, iono_south_im_eFluxElec, &
                iono_north_im_aveeHydr, iono_south_im_aveeHydr, &
                iono_north_im_efluxHydr, iono_south_im_eFluxHydr

          real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
              AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
              EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

          character(len=*), intent(in) :: NameHemiIn
      character(len=*), parameter:: NameSub = 'imp_integrated_flux'
      !------------------------------------------------------------------------
          if (trim(NameHemiIn) == 'south') then
              AvgEDiffe_II = iono_south_im_aveeElec ! in keV
              EfluxDiffe_II = iono_south_im_efluxElec / 1000.0 ! mW/m^2 to W/m^2
              AvgEDiffi_II = iono_south_im_aveeHydr
              EfluxDiffi_II = iono_south_im_efluxHydr / 1000.0
          else if (trim(NameHemiIn) == 'north') then
              AvgEDiffe_II = iono_north_im_aveeElec
              EfluxDiffe_II = iono_north_im_efluxElec / 1000.0
              AvgEDiffi_II = iono_north_im_aveeHydr
              EfluxDiffi_II = iono_north_im_efluxHydr / 1000.0
           else
              call CON_stop(NameSub//' : unrecognized hemisphere - '//&
                      NameHemiIn)
           end if

      end subroutine imp_integrated_flux
    !==========================================================================
      subroutine imp_spectral_flux(NameHemiIn, AvgEDiffe_II, AvgEDiffi_II, &
              AvgEMono_II, AvgEBbnd_II, EfluxDiffe_II, EfluxDiffi_II, &
              EfluxMono_II, EfluxBbnd_II)
          use ModIonosphere, ONLY: &
                iono_north_im_eElecPrec, iono_south_im_eElecPrec, &
                iono_north_im_eHydrPrec, iono_south_im_eHydrPrec, &
                iono_north_im_nElecPrec, iono_south_im_nElecPrec, &
                iono_north_im_nHydrPrec, iono_south_im_nHydrPrec

          real, intent(out), dimension(IONO_nTheta, IONO_nPsi) :: &
                  AvgEDiffe_II, AvgEDiffi_II, AvgEMono_II, AvgEBbnd_II, &
                  EfluxDiffe_II, EfluxDiffi_II, EfluxMono_II, EfluxBbnd_II

          character(len=*), intent(in) :: NameHemiIn
      character(len=*), parameter:: NameSub = 'imp_spectral_flux'
      !------------------------------------------------------------------------
          call CON_stop(NameSub//' not yet implemented!')
      end subroutine imp_spectral_flux
    !==========================================================================
  end subroutine imp_gen_fluxes
  !============================================================================

end module ModImp
!==============================================================================
