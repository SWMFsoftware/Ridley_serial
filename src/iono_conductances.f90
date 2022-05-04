!  Copyright (C) 2002 Regents of the University of Michigan,
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!*************************************************************************
!
! CONDUCTANCE Calculation Routines
!
!*************************************************************************

subroutine FACs_to_fluxes(iModel, iBlock)

  ! The goal here is to convert the ionospheric FAC pattern into a 
  ! particle precipitation pattern, which can then be turned into
  ! a conductance pattern.

  use ModIonosphere
  use IE_ModMain
  use ModNumConst
  use ModConst, ONLY: cEVToK, cElectronCharge
  use ModConductance, ONLY: UseSubOvalCond, UseOval, UseNewOval, &
       smooth_lagrange_polar
  implicit none

  integer, intent(in) :: iModel  ! model number, 
  ! same as conductance_model in IE_ModMain

  integer, intent(in) :: iBlock  ! 1 for northern, 2 for southern hemisphere

  real    :: aurora_peak_lat, eflux_peak_lon, ave_e_peak_lon
  real    :: aurora_peak_width, ave_e_peak
  real    :: distance, mult_fac, ef, efp, mfc
  real    :: hal_a0, hal_a1, hal_a2, ped_a0, ped_a1, ped_a2, hall, ped
  real    :: dlat, dmlt, y1, y2, x1, x2
  real    :: day_colat, dusk_colat, midnight_colat, dawn_colat
  real    :: day_fac, dusk_fac, midnight_fac, dawn_fac
  real    :: noon_mid, dusk_dawn, oval_shift
  real    :: mean_colat, dev_colat, night_width
  real    :: PolarCapHallConductance, PolarCap_AveE, PolarCap_EFlux
  integer :: jlat, imlt
  integer :: i,j, n, nloc, nHalfSmooth
  real, dimension(1:IONO_nPsi) :: Strength_of_Oval,               &
       Loc_of_Oval, Width_of_Oval
  logical :: polarcap, IsPeakFound, IsDone
  real :: Center, Width
  real :: f, MaxP, MaxT, MulFac_ae, MulFac_ef, MinWidth, ThetaOCB, AuroraWidth
  real, dimension(IONO_nTheta,IONO_nPsi) :: nDen, &
       discrete_k, discrete_ae, discrete_ef, diffuse_ae, diffuse_ef, &
       iono_Ne, iono_T, discrete_nf, diffuse_nf, eV, av_mono_value
  real, dimension(IONO_nPsi) :: OCFLB, EquatorwardEdge, Smooth, OCFLB_s, beta
  real :: MulFac_Dae, MulFac_Def, var_rm, var_tot, numerator, denominator, &
       rm, vari, numflux_floor
  !real, dimension(IONO_nTheta,IONO_nPsi) :: nDen, &
  !     discrete_k, discrete_ae, discrete_ef, diffuse_ae, diffuse_ef, &
  !     iono_Ne, iono_T, discrete_nf, diffuse_nf, eV
  !real, dimension(IONO_nPsi) :: OCFLB, OCFLB_s, beta, EquatorwardEdge, Smooth
  !real :: f, MaxP, MaxT, MulFac_ae, MulFac_ef, MinWidth, ThetaOCB, AuroraWidth
  !real :: MulFac_Dae, MulFac_Def
  !real, dimension(IONO_nTheta,IONO_nPsi) :: nDen, &
  !     discrete_k, discrete_ae, discrete_ef, diffuse_ae, diffuse_ef
  !real, dimension(IONO_nPsi) :: OCFLB, EquatorwardEdge, Smooth
  !real :: f, MaxP, MaxT, MulFac_ae, MulFac_ef, MinWidth, ThetaOCB, AuroraWidth
  !real :: MulFac_Dae, MulFac_Def, var_rm, var_tot, numerator, denominator, &
  !     rm, vari, av_mono_value
  !real, dimension(IONO_nTheta,IONO_nPsi) :: nDen, &
  !     discrete_k, discrete_ae, discrete_ef, diffuse_ae, diffuse_ef, &
  !     iono_Ne, iono_T, discrete_nf, diffuse_nf, eV
  !real, dimension(IONO_nPsi) :: OCFLB, OCFLB_s, beta, EquatorwardEdge, Smooth

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='facs_to_fluxes'
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  
  Hall_to_Ped_Ratio = 1.5

  if (PolarCapPedConductance > 0.0) then
     PolarCapHallConductance = Hall_to_Ped_Ratio * PolarCapPedConductance
     PolarCap_AveE = (Hall_to_Ped_Ratio/0.45)**(1.0/0.85)
     PolarCap_EFlux = ((PolarCapPedConductance*(16.0 + PolarCap_AveE**2) / &
          (40.0*PolarCap_AveE))**2)/1000.0
  else
     PolarCap_AveE  = IONO_Min_Ave_E
     PolarCap_EFlux = IONO_Min_EFlux
  endif

  if (iModel.eq.7) then

     write(*,*) "maxval eflux (north) : ",maxval(iono_north_im_eflux)

     iono_north_eflux = iono_north_im_eflux/1000.0
     iono_north_ave_e = iono_north_im_avee

     where (iono_north_ave_e > 20.0) &
          iono_north_ave_e = 20.0

     do j = 1, IONO_nPsi
        i = 1
        do while (iono_north_eflux(i,j) == 0.0 .and. i < Iono_nTheta)
           i = i + 1
        enddo
        if (i > 2) iono_north_eflux(i-1,j) = iono_north_eflux(i,j)* 0.66 
        if (i > 3) iono_north_eflux(i-2,j) = iono_north_eflux(i,j)* 0.33
        if (i > 4) iono_north_eflux(1:i-3,j) = PolarCap_EFlux
        if (i > 2) iono_north_ave_e(i-1,j) = iono_north_ave_e(i,j)
        if (i > 3) iono_north_ave_e(i-2,j) = iono_north_ave_e(i,j)
        if (i > 4) iono_north_ave_e(1:i-3,j) = PolarCap_AveE
     enddo

     do i = 1, IONO_nTheta
        iono_south_eflux(i,:) = iono_north_eflux(Iono_nTheta-i+1,:)
        iono_south_ave_e(i,:) = iono_north_ave_e(Iono_nTheta-i+1,:)
     enddo

     where (iono_south_ave_e < IONO_Min_Ave_E) &
          iono_south_ave_e = IONO_Min_Ave_E
     where (iono_south_eflux < IONO_Min_EFlux) &
          iono_south_eflux = IONO_Min_EFlux

     where (iono_north_ave_e < IONO_Min_Ave_E) &
          iono_north_ave_e = IONO_Min_Ave_E
     where (iono_north_eflux < IONO_Min_EFlux) &
          iono_north_eflux = IONO_Min_EFlux

  endif

  if (iModel.eq.3) then

     mult_fac   = 0.25*0.33e7*2.0e-3
     ave_e_peak = max(3.0, PolarCap_AveE)                        ! keV

     oval_shift = 10.0*cDegToRad

     select case(iBlock)
     case(1)

        ! Do North first -----------------------------------------------------
        if(UseNewOval)then
           call Create_Auroral_Oval(IONO_NORTH_JR, IONO_NORTH_Theta, &
                IONO_NORTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        else
           call Determine_Oval_Characteristics(IONO_NORTH_JR, IONO_NORTH_Theta,&
                IONO_NORTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        end if

        ave_e_peak_lon = 0.0

        do i = 1, IONO_nTheta
           do j = 1, IONO_nPsi

              ! Energy Flux

              distance = IONO_NORTH_Theta(i,j)-(Loc_of_Oval(j)+oval_shift)
              IONO_NORTH_EFlux(i,j) =                                       &
                   Strength_of_Oval(j) * &
                   mult_fac * &
                   exp(-1.0*(distance/Width_of_Oval(j))**2)

              if (distance < 0.0) then
                 IONO_NORTH_EFlux(i,j) = &
                      max(IONO_NORTH_EFlux(i,j), PolarCap_EFlux)
                 polarcap = .true.
              else
                 IONO_NORTH_EFlux(i,j) = &
                      max(IONO_NORTH_EFlux(i,j), IONO_Min_EFlux)
                 polarcap = .false.
              endif

              ! Average Energy

              IONO_NORTH_Ave_E(i,j) =                                      &
                   ave_e_peak*exp(-1.0*(distance/Width_of_Oval(j))**2)

              if (PolarCap_AveE == IONO_Min_Ave_E) then
                 distance = 0.25 +                                          &
                      0.75*(                                                &
                      sin(IONO_NORTH_Psi(i,j)-(ave_e_peak_lon + cHalfPi)) &
                      + 1.0)/2.0
              else
                 distance = 1.0
              endif

              IONO_NORTH_Ave_E(i,j) = IONO_NORTH_Ave_E(i,j)*distance

              if (polarcap) then
                 IONO_NORTH_Ave_E(i,j) = &
                      max(IONO_NORTH_Ave_E(i,j), PolarCap_AveE)
              else
                 IONO_NORTH_Ave_E(i,j) = &
                      max(IONO_NORTH_Ave_E(i,j), IONO_Min_Ave_E)
              endif

           enddo
        enddo

     case(2)

        ! Do South next --------------------------------------------------------
        if(UseNewOval)then
           call Create_Auroral_Oval(IONO_SOUTH_JR, IONO_SOUTH_Theta, &
                IONO_SOUTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        else
           call Determine_Oval_Characteristics(IONO_SOUTH_JR, IONO_SOUTH_Theta,&
                IONO_SOUTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        end if

        Loc_of_Oval = cPi - Loc_of_Oval

        ave_e_peak_lon = 0.0

        do i = 1, IONO_nTheta
           do j = 1, IONO_nPsi

              ! Energy Flux

              distance = IONO_SOUTH_Theta(i,j)-(Loc_of_Oval(j)-oval_shift)
              IONO_SOUTH_EFlux(i,j) =                                            &
                   Strength_of_Oval(j) * &
                   mult_fac * &
                   exp(-1.0*(distance/Width_of_Oval(j))**2)

              if (distance > 0.0) then
                 IONO_SOUTH_EFlux(i,j) = &
                      max(IONO_SOUTH_EFlux(i,j), PolarCap_EFlux)
                 polarcap = .true.
              else
                 IONO_SOUTH_EFlux(i,j) = &
                      max(IONO_SOUTH_EFlux(i,j), IONO_Min_EFlux)
                 polarcap = .false.
              endif

              ! Average Energy

              IONO_SOUTH_Ave_E(i,j) =                                      &
                   ave_e_peak*exp(-1.0*(distance/Width_of_Oval(j))**2)

              if (PolarCap_AveE == IONO_Min_Ave_E) then
                 distance = 0.25 +                                          &
                      0.75*(                                                &
                      sin(IONO_SOUTH_Psi(i,j)-(ave_e_peak_lon + cHalfPi))   &
                      + 1.0)/2.0
              else
                 distance = 1.0
              endif

              IONO_SOUTH_Ave_E(i,j) = IONO_SOUTH_Ave_E(i,j)*distance

              if (polarcap) then
                 IONO_SOUTH_Ave_E(i,j) = &
                      max(IONO_SOUTH_Ave_E(i,j), PolarCap_AveE)
              else
                 IONO_SOUTH_Ave_E(i,j) = &
                      max(IONO_SOUTH_Ave_E(i,j), IONO_Min_Ave_E)
              endif
           enddo
        enddo
     end select
  endif

  if (iModel.eq.11) then ! ADELPHI AMPERE-derived model

     select case(iBlock)
     case(1)
        !  Do North First

        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              ! Robinson ADELPHI model 8/3/2021
              ! iono_north_jr in units of Amps/m^2
              ! IONO_NORTH_Psi is mlt in units of radians 
              
              call ADELPHI_Calc(iono_north_jr(i,j),IONO_NORTH_Psi(i,j),hall,ped)
              
              if (abs(90 - IONO_NORTH_Theta(i,j)*cRadToDeg) <= 50) then
                 hall=0
                 ped=0
              end if
              
              if ((hall.gt.1.0).and.(ped.gt.0.5)) then
                 IONO_NORTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
                 IONO_NORTH_EFlux(i,j) = (ped*(16.0+IONO_NORTH_Ave_E(i,j)**2)/&
                      (40.0*IONO_NORTH_Ave_E(i,j)))**2/1000.0
              else
                 IONO_NORTH_Ave_E(i,j) = IONO_Min_Ave_E
                 IONO_NORTH_EFlux(i,j) = IONO_Min_EFlux
              endif
           end do
        end do
     case(2)
        !  Do South Next
        
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              
              call ADELPHI_Calc(iono_south_jr(i,j),IONO_SOUTH_Psi(i,j),hall,ped)
              
              if (abs(90 - IONO_SOUTH_Theta(i,j)*cRadToDeg) <= 50) then
                 hall=0
                 ped=0
              end if
              
              if ((hall.gt.1.0).and.(ped.gt.0.5)) then
                 IONO_SOUTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
                 IONO_SOUTH_EFlux(i,j) = (ped*(16.0+IONO_SOUTH_Ave_E(i,j)**2)/&
                      (40.0*IONO_SOUTH_Ave_E(i,j)))**2/1000.0
              else
                 IONO_SOUTH_Ave_E(i,j) = IONO_Min_Ave_E
                 IONO_SOUTH_EFlux(i,j) = IONO_Min_EFlux
              endif
              
           enddo
        enddo
     end select
  end if

  if (iModel.eq.4 .or. iModel.eq.5 .or. iModel.eq.9) then

     ! Calculate grid spacing for conductance grid
     dlat = (cond_lats(1) - cond_lats(2))*cDegToRad
     dmlt = (cond_mlts(2) - cond_mlts(1))*cPi/12.0


     select case(iBlock)
     case(1)
        !  Do North First

        if(UseNewOval)then
           call Create_Auroral_Oval(IONO_NORTH_JR, IONO_NORTH_Theta, &
                IONO_NORTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        else
           call Determine_Oval_Characteristics(IONO_NORTH_JR, IONO_NORTH_Theta,&
                IONO_NORTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        end if

        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta

              y1 = IONO_NORTH_Theta(i,j)/dlat + 1.0
              if (y1 > i_cond_nlats-1) then
                 jlat = i_cond_nlats-1
                 y1   = 1.0
              else
                 jlat = y1
                 y1   = 1.0 - (y1 - jlat)
              endif
              y2 = 1.0 - y1

              x1 = mod((IONO_NORTH_Psi(i,j) + cPi), cTwoPi)/dmlt + 1.0
              imlt = x1
              x1   = 1.0 - (x1 - imlt)
              x2   = 1.0 - x1

              if (iono_north_jr(i,j) > 0) then

                 hal_a0 = x1*y1*hal_a0_up(imlt  ,jlat  ) + &
                      x2*y1*hal_a0_up(imlt+1,jlat  ) + &
                      x1*y2*hal_a0_up(imlt  ,jlat+1) + &
                      x2*y2*hal_a0_up(imlt+1,jlat+1)

                 hal_a1 = x1*y1*hal_a1_up(imlt  ,jlat  ) + &
                      x2*y1*hal_a1_up(imlt+1,jlat  ) + &
                      x1*y2*hal_a1_up(imlt  ,jlat+1) + &
                      x2*y2*hal_a1_up(imlt+1,jlat+1)

                 hal_a2 = x1*y1*hal_a2_up(imlt  ,jlat  ) + &
                      x2*y1*hal_a2_up(imlt+1,jlat  ) + &
                      x1*y2*hal_a2_up(imlt  ,jlat+1) + &
                      x2*y2*hal_a2_up(imlt+1,jlat+1)

                 ped_a0 = x1*y1*ped_a0_up(imlt  ,jlat  ) + &
                      x2*y1*ped_a0_up(imlt+1,jlat  ) + &
                      x1*y2*ped_a0_up(imlt  ,jlat+1) + &
                      x2*y2*ped_a0_up(imlt+1,jlat+1)

                 ped_a1 = x1*y1*ped_a1_up(imlt  ,jlat  ) + &
                      x2*y1*ped_a1_up(imlt+1,jlat  ) + &
                      x1*y2*ped_a1_up(imlt  ,jlat+1) + &
                      x2*y2*ped_a1_up(imlt+1,jlat+1)

                 ped_a2 = x1*y1*ped_a2_up(imlt  ,jlat  ) + &
                      x2*y1*ped_a2_up(imlt+1,jlat  ) + &
                      x1*y2*ped_a2_up(imlt  ,jlat+1) + &
                      x2*y2*ped_a2_up(imlt+1,jlat+1)                 

              else

                 hal_a0 = x1*y1*hal_a0_do(imlt  ,jlat  ) + &
                      x2*y1*hal_a0_do(imlt+1,jlat  ) + &
                      x1*y2*hal_a0_do(imlt  ,jlat+1) + &
                      x2*y2*hal_a0_do(imlt+1,jlat+1)

                 hal_a1 = x1*y1*hal_a1_do(imlt  ,jlat  ) + &
                      x2*y1*hal_a1_do(imlt+1,jlat  ) + &
                      x1*y2*hal_a1_do(imlt  ,jlat+1) + &
                      x2*y2*hal_a1_do(imlt+1,jlat+1)

                 hal_a2 = x1*y1*hal_a2_do(imlt  ,jlat  ) + &
                      x2*y1*hal_a2_do(imlt+1,jlat  ) + &
                      x1*y2*hal_a2_do(imlt  ,jlat+1) + &
                      x2*y2*hal_a2_do(imlt+1,jlat+1)

                 ped_a0 = x1*y1*ped_a0_do(imlt  ,jlat  ) + &
                      x2*y1*ped_a0_do(imlt+1,jlat  ) + &
                      x1*y2*ped_a0_do(imlt  ,jlat+1) + &
                      x2*y2*ped_a0_do(imlt+1,jlat+1)

                 ped_a1 = x1*y1*ped_a1_do(imlt  ,jlat  ) + &
                      x2*y1*ped_a1_do(imlt+1,jlat  ) + &
                      x1*y2*ped_a1_do(imlt  ,jlat+1) + &
                      x2*y2*ped_a1_do(imlt+1,jlat+1)

                 ped_a2 = x1*y1*ped_a2_do(imlt  ,jlat  ) + &
                      x2*y1*ped_a2_do(imlt+1,jlat  ) + &
                      x1*y2*ped_a2_do(imlt  ,jlat+1) + &
                      x2*y2*ped_a2_do(imlt+1,jlat+1)

              endif

              ! "Distance" sets if auroral oval is used/relevant.
              ! To turn off auroral oval, set to large negative number.
              ! This will put all locations in polar cap AND reduce
              ! impact of oval conductance to near-zero.
              if(UseOval)then
                 distance = (IONO_NORTH_Theta(i,j) - Loc_of_Oval(j))
              else
                 distance = -1E9
              end if

              if (distance > 0.0) then
                 polarcap = .false.
              else
                 polarcap = .true.
              endif
              
              if (iModel.eq.9) then
                 ! A simple power law relation between FAC and conductance was added by Zihan Wang. 02/26/2021.
                 ! Conductance files in PARAM.IN needs to be changed. The format is the same. However, a2 will not be used.
                
                 hall=exp(hal_a0+hal_a1*log(abs(iono_north_jr(i,j)*1.0e6)))
                 ped=exp(ped_a0+ped_a1*log(abs(iono_north_jr(i,j)*1.0e6)))
                 
              endif
              
              if (iModel.eq.4) then
                 ! Implemented Feb. 7, 2007 as modified version of iModel 5 with
                 !    new narrower fitting of the auroral oval.  DDZ
                 
                 hall=exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      CondFactor*( &
                      hal_a0+(hal_a1-hal_a0)*exp(-abs(iono_north_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      CondFactor*( &
                      ped_a0+(ped_a1-ped_a0)*exp(-abs(iono_north_jr(i,j)*1.0e9)*ped_a2**2))
              endif
                 
              if (iModel.eq.5) then
                 !
                 ! We want minimal conductance lower than the oval
                 !
                 
                 ! decrease conductance equatorward of the center of the auroral oval.
                 ! The decrease is slow down with a factor of 3.
                 if (.not.polarcap .and. .not.UseSubOvalCond) then
                    distance = distance/3.0
                    hal_a0 = hal_a0 * exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)
                    ped_a0 = ped_a0 * exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)
                 end if

                 ! No Conductance Beyond 45 degrees for UseCMEEFitting
                 if (UseCMEEFitting) then
                    ! IONO_NORTH_Theta is in theta, need to switch to latitude 
                    if (abs(90 - IONO_NORTH_Theta(i,j)*cRadToDeg) <= LatNoConductanceSI) then
                       hal_a1 = hal_a0
                       ped_a1 = ped_a0
                    endif
                 end if

                 !
                 ! A sort of correction on the fit
                 !

                 hal_a1 = hal_a0 + (OvalStrengthFactor*hal_a1 - hal_a0)*  &
                      exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)
                 ped_a1 = ped_a0 + (OvalStrengthFactor*ped_a1 - ped_a0)*  &
                      exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)

                 if (UseCMEEFitting) then
                    hal_a1 = min(hal_a0, hal_a1)
                    ped_a1 = min(ped_a0, ped_a1)
                 end if

                 ! Multiply by sqrt(3) to compensate for the 3 times narrower oval
                 hall=CondFactor*( &
                      hal_a0-hal_a1*exp(-abs(iono_north_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =CondFactor*( &
                      ped_a0-ped_a1*exp(-abs(iono_north_jr(i,j)*1.0e9)*ped_a2**2))

                 if (UseCMEEFitting) then
                    hall=hall + CondFactor*( &
                         FactorHallCMEE*exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2))
                    ped =ped  + CondFactor*( &
                         FactorPedCMEE*exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2))
                 endif
              end if

              if ((hall.gt.1.0).and.(ped.gt.0.5)) then
                 IONO_NORTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
                 IONO_NORTH_EFlux(i,j) = (ped*(16.0+IONO_NORTH_Ave_E(i,j)**2)/&
                      (40.0*IONO_NORTH_Ave_E(i,j)))**2/1000.0
              else
                 IONO_NORTH_Ave_E(i,j) = IONO_Min_Ave_E
                 IONO_NORTH_EFlux(i,j) = IONO_Min_EFlux
              endif

              if ((PolarCap_AveE > 0.0).and.(polarcap)) then
                 IONO_NORTH_Ave_E(i,j) = max(IONO_NORTH_Ave_E(i,j), &
                      PolarCap_AveE)
                 IONO_NORTH_EFlux(i,j) = max(IONO_NORTH_EFlux(i,j), &
                      PolarCap_EFlux)
              endif

           enddo
        enddo

     case(2)
        !  Do South Next

        if(UseNewOval)then
           call Create_Auroral_Oval(IONO_SOUTH_JR, IONO_SOUTH_Theta, &
                IONO_SOUTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        else
           call Determine_Oval_Characteristics(IONO_SOUTH_JR, IONO_SOUTH_Theta,&
                IONO_SOUTH_Psi, Loc_of_Oval, Width_of_Oval,Strength_of_Oval)
        end if

        Loc_of_Oval = cPI - Loc_of_Oval

        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              y1 = (cPi - IONO_SOUTH_Theta(i,j))/dlat + 1.0
              if (y1 > i_cond_nlats-1) then
                 jlat = i_cond_nlats-1
                 y1   = 1.0
              else
                 jlat = y1
                 y1   = 1.0 - (y1 - jlat)
              endif
              y2 = 1.0 - y1

              x1 = mod((IONO_SOUTH_Psi(i,j) + cPi), cTwoPi)/dmlt + 1.0
              imlt = x1
              x1   = 1.0 - (x1 - imlt)
              x2   = 1.0 - x1

              if (iono_south_jr(i,j) > 0) then

                 hal_a0 = x1*y1*hal_a0_up(imlt  ,jlat  ) + &
                      x2*y1*hal_a0_up(imlt+1,jlat  ) + &
                      x1*y2*hal_a0_up(imlt  ,jlat+1) + &
                      x2*y2*hal_a0_up(imlt+1,jlat+1)

                 hal_a1 = x1*y1*hal_a1_up(imlt  ,jlat  ) + &
                      x2*y1*hal_a1_up(imlt+1,jlat  ) + &
                      x1*y2*hal_a1_up(imlt  ,jlat+1) + &
                      x2*y2*hal_a1_up(imlt+1,jlat+1)

                 hal_a2 = x1*y1*hal_a2_up(imlt  ,jlat  ) + &
                      x2*y1*hal_a2_up(imlt+1,jlat  ) + &
                      x1*y2*hal_a2_up(imlt  ,jlat+1) + &
                      x2*y2*hal_a2_up(imlt+1,jlat+1)

                 ped_a0 = x1*y1*ped_a0_up(imlt  ,jlat  ) + &
                      x2*y1*ped_a0_up(imlt+1,jlat  ) + &
                      x1*y2*ped_a0_up(imlt  ,jlat+1) + &
                      x2*y2*ped_a0_up(imlt+1,jlat+1)

                 ped_a1 = x1*y1*ped_a1_up(imlt  ,jlat  ) + &
                      x2*y1*ped_a1_up(imlt+1,jlat  ) + &
                      x1*y2*ped_a1_up(imlt  ,jlat+1) + &
                      x2*y2*ped_a1_up(imlt+1,jlat+1)

                 ped_a2 = x1*y1*ped_a2_up(imlt  ,jlat  ) + &
                      x2*y1*ped_a2_up(imlt+1,jlat  ) + &
                      x1*y2*ped_a2_up(imlt  ,jlat+1) + &
                      x2*y2*ped_a2_up(imlt+1,jlat+1)

              else

                 hal_a0 = x1*y1*hal_a0_do(imlt  ,jlat  ) + &
                      x2*y1*hal_a0_do(imlt+1,jlat  ) + &
                      x1*y2*hal_a0_do(imlt  ,jlat+1) + &
                      x2*y2*hal_a0_do(imlt+1,jlat+1)

                 hal_a1 = x1*y1*hal_a1_do(imlt  ,jlat  ) + &
                      x2*y1*hal_a1_do(imlt+1,jlat  ) + &
                      x1*y2*hal_a1_do(imlt  ,jlat+1) + &
                      x2*y2*hal_a1_do(imlt+1,jlat+1)

                 hal_a2 = x1*y1*hal_a2_do(imlt  ,jlat  ) + &
                      x2*y1*hal_a2_do(imlt+1,jlat  ) + &
                      x1*y2*hal_a2_do(imlt  ,jlat+1) + &
                      x2*y2*hal_a2_do(imlt+1,jlat+1)

                 ped_a0 = x1*y1*ped_a0_do(imlt  ,jlat  ) + &
                      x2*y1*ped_a0_do(imlt+1,jlat  ) + &
                      x1*y2*ped_a0_do(imlt  ,jlat+1) + &
                      x2*y2*ped_a0_do(imlt+1,jlat+1)

                 ped_a1 = x1*y1*ped_a1_do(imlt  ,jlat  ) + &
                      x2*y1*ped_a1_do(imlt+1,jlat  ) + &
                      x1*y2*ped_a1_do(imlt  ,jlat+1) + &
                      x2*y2*ped_a1_do(imlt+1,jlat+1)

                 ped_a2 = x1*y1*ped_a2_do(imlt  ,jlat  ) + &
                      x2*y1*ped_a2_do(imlt+1,jlat  ) + &
                      x1*y2*ped_a2_do(imlt  ,jlat+1) + &
                      x2*y2*ped_a2_do(imlt+1,jlat+1)

              endif

              ! Use "distance" to turn oval off/on (see above comment).
              if(UseOval)then
                 distance = (IONO_SOUTH_Theta(i,j) - Loc_of_Oval(j))
              else
                 distance = 1E9 ! Switch sign for southern hemisphere.
              end if

              if (distance < 0.0) then
                 polarcap = .false.
              else
                 polarcap = .true.
              endif

              if (iModel.eq.9) then
                 ! A simple power law relation between FAC and conductance was added by Zihan Wang. 02/26/2021.
                 ! Conductance files in PARAM.IN needs to be changed. The format is the same. However, a2 will not be used.
                 
                 hall=exp(hal_a0+hal_a1*log(abs(iono_south_jr(i,j)*1.0e6)))
                 ped=exp(ped_a0+ped_a1*log(abs(iono_south_jr(i,j)*1.0e6)))

              endif

              if (iModel.eq.4) then
                 ! Implemented Feb. 7, 2007 as modified version of iModel 5 with
                 !    new narrower fitting of the auroral oval.  DDZ
                 
                 hall=exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      CondFactor*( &
                      hal_a0+(hal_a1-hal_a0)*exp(-abs(iono_south_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2) * &
                      CondFactor*( &
                      ped_a0+(ped_a1-ped_a0)*exp(-abs(iono_south_jr(i,j)*1.0e9)*ped_a2**2))
              endif
                 
              if (iModel.eq.5) then
                 ! Restrict FAC-related conductance outside auroral oval.
                 if (.not.polarcap .and. .not.UseSubOvalCond) then
                    distance = distance/3.0
                    hal_a0 = hal_a0 * exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)
                    ped_a0 = ped_a0 * exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)
                 endif

                 ! No Conductance Below 45 degrees
                 if (UseCMEEFitting) then
                    ! IONO_SOUTH_Theta is in theta, need to switch to latitude 
                    if (abs(90 - IONO_SOUTH_Theta(i,j)*cRadToDeg) <= LatNoConductanceSI) then
                       hal_a1 = hal_a0
                       ped_a1 = ped_a0
                    end if
                 end if

                 !
                 ! A sort of correction on the fit
                 !

                 hal_a1 = hal_a0 + (OvalStrengthFactor*hal_a1 - hal_a0)*  &
                      exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)
                 ped_a1 = ped_a0 + (OvalStrengthFactor*ped_a1 - ped_a0)*  &
                      exp(-1.0*(distance/(OvalWidthFactor*Width_of_Oval(j)))**2)

                 ! No Negatives
                 if (UseCMEEFitting) then
                    hal_a1 = min(hal_a0, hal_a1)
                    ped_a1 = min(ped_a0, ped_a1)
                 end if

                 ! Multiply by sqrt(3) to compensate for the 3 times narrower oval
                 hall=CondFactor*( &
                      hal_a0-hal_a1*exp(-abs(iono_south_jr(i,j)*1.0e9)*hal_a2**2))
                 ped =CondFactor*( &
                      ped_a0-ped_a1*exp(-abs(iono_south_jr(i,j)*1.0e9)*ped_a2**2))

                 if (UseCMEEFitting) then
                    hall=hall + CondFactor*( &
                         FactorHallCMEE*exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2))
                    ped =ped  + CondFactor*( &
                         FactorPedCMEE *exp(-(distance/(OvalWidthFactor*Width_of_Oval(j)))**2))
                 end if
              end if

              if ((hall.gt.1.0).and.(ped.gt.0.5)) then
                 IONO_SOUTH_Ave_E(i,j)  = ((hall/ped)/0.45)**(1.0/0.85)
                 IONO_SOUTH_EFlux(i,j) = (ped*(16.0+IONO_SOUTH_Ave_E(i,j)**2)/&
                      (40.0*IONO_SOUTH_Ave_E(i,j)))**2/1000.0
              else
                 IONO_SOUTH_Ave_E(i,j) = IONO_Min_Ave_E
                 IONO_SOUTH_EFlux(i,j) = IONO_Min_EFlux
              endif

              if ((PolarCap_AveE > 0.0).and.(polarcap)) then
                 IONO_SOUTH_Ave_E(i,j) = max(IONO_SOUTH_Ave_E(i,j), &
                      PolarCap_AveE)
                 IONO_SOUTH_EFlux(i,j) = max(IONO_SOUTH_EFlux(i,j),        &
                      PolarCap_EFlux)
              endif

           enddo
        enddo
     end select
  end if

  if (iModel.eq.6) then

     MinWidth = 5.0 * cPi / 180.0
     nHalfSmooth = 5
     !     MulFac_Dae = 1.0e22
     !     MulFac_Def = 5.0e19
     !     MulFac_ef = 0.2e7
     !     MulFac_ae = 1.0 / 1.0e11
     MulFac_Dae = 1.0e22
     MulFac_Def = 1.0e19
     MulFac_ef = 0.3e6
     MulFac_ae = 4.0e-12

     if (iBlock == 1) then 

        iono_north_eflux = 1.0e-6
        iono_north_ave_e = 1.0

        OCFLB = -1.0
        EquatorWardEdge = -1.0

        do j = 1, IONO_nPsi
           IsPeakFound = .false.
           MaxP = max(maxval(iono_north_p(:,IONO_nPsi-j+1)),1.0e-15)
           IsDone = .false.
           i = 1
           do while (.not. IsDone)

              if (iono_north_p(i,IONO_nPsi-j+1) > 0) then

                 if (OCFLB(j) == -1) OCFLB(j) = iono_north_theta(i,j)
                 AuroraWidth = iono_north_theta(i,j) - OCFLB(j)
                 iono_north_eflux(i,j) = MulFac_ef*MaxP

                 if (iono_north_p(i,IONO_nPsi-j+1)==MaxP) IsPeakFound = .true.

                 if (IsPeakFound .and. AuroraWidth >= MinWidth) then
                    EquatorWardEdge(j) = iono_north_theta(i,j)
                    IsDone = .true.
                 endif

                 if (i == IONO_nTheta) IsDone = .true.

              endif

              i = i + 1

           enddo
           if (.not. IsPeakFound) then
              OCFLB(j) = MinWidth
              EquatorwardEdge(j) = MinWidth*2
           endif

        enddo

        do j = 1, IONO_nPsi
           smooth(j) = 0.0
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + OCFLB(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        OCFLB = smooth/(nHalfSmooth*2+1)

        do j = 1, IONO_nPsi
           smooth(j) = 0.0
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + &
                   EquatorWardEdge(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        EquatorWardEdge = smooth/(nHalfSmooth*2+1)

        do j = 1, IONO_nPsi

           Center = (OCFLB(j)*2.0 + EquatorWardEdge(j))/3.0
           Width = abs(EquatorWardEdge(j) - OCFLB(j))/2

           MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)
           iono_north_eflux(:,j) = &
                MulFac_ef*MaxP * &
                exp(-abs(Center-iono_north_theta(:,j))/Width) * &
                (0.375*cos(iono_north_psi(:,j)-cPi)+0.625)

           ! This is the way that seems to work
           iono_north_ave_e(:,j) = &
                iono_north_p(:,IONO_nPsi-j+1) / &
                (iono_north_rho(:,IONO_nPsi-j+1)+1.0e-32) * MulFac_ae

           !           iono_north_ave_e(Poleward-2,j) = iono_north_ave_e(Poleward+1,j)
           !           iono_north_ave_e(Poleward-1,j) = iono_north_ave_e(Poleward+1,j)

           !           do i = Poleward-2, Equatorward+2
           !           do i = Poleward, Equatorward
           !              f = exp(-abs(float(i-Center))/float(width))
           !              iono_north_eflux(i,j) = iono_north_eflux(i,j) * f
           !           enddo

        enddo

        diffuse_ef = iono_north_eflux
        diffuse_ae = iono_north_ave_e

        discrete_ef = 0.0
        discrete_k  = 0.0

        ! This is from Jimmy's Paper on the Knight Relationship
        where (iono_north_p > 0) discrete_k = &
             (iono_north_rho**1.5) / iono_north_p
        ! Mirror the particles, since we are dealing with electrons and not
        ! ions.
        do j = 1, IONO_nPsi
           discrete_k(:,j) = discrete_k(:,IONO_nPsi-j+1)
        enddo
        where (iono_north_jr > 0.0) &
             discrete_ef = (iono_north_jr*1e6)*discrete_k
        discrete_ae = discrete_ef*MulFac_Dae
        discrete_ef = (iono_north_jr*1e6)*discrete_ef*MulFac_Def

        ! Let's add a little conductance on ANY not in the polar cap,
        ! so the code doesn't blow up
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              if (abs(iono_north_jr(i,j)) > 0.0 .and. &
                   discrete_ef(i,j) < 5.0e-3) then
                 discrete_ef(i,j) = MulFac_Def/2.5 * &
                      (iono_north_jr(i,j)*1e6)**2*discrete_k(i,j)
              endif
           enddo
        enddo

        where (diffuse_ae < IONO_Min_Ave_E/2) diffuse_ae = IONO_Min_Ave_E/2
        where (discrete_ae < IONO_Min_Ave_E/2) discrete_ae = IONO_Min_Ave_E/2

        ! Let's weight the average energy by the number flux, which is ef/av
        iono_north_ave_e = &
             (diffuse_ef + discrete_ef) / &
             (diffuse_ef/diffuse_ae + discrete_ef/discrete_ae)

        ! The energy flux should be weighted by the average energy
        iono_north_eflux = &
             (diffuse_ef/diffuse_ae + discrete_ef/discrete_ae) * &
             iono_north_ave_e

        where (iono_north_ave_e < IONO_Min_Ave_E) &
             iono_north_ave_e = IONO_Min_Ave_E

        !        ! Let's add a little conductance on ANY not in the polar cap,
        !        ! so the code doesn't blow up
        !        do j = 1, IONO_nPsi
        !           do i = 1, IONO_nTheta
        !              if (iono_north_p(i,j) <= 0.0) then
        !                 iono_north_ave_e(i,j) = max(iono_north_ave_e(i,j),polarcap_avee)
        !                 iono_north_eflux(i,j) = &
        !                      max(iono_north_eflux(i,j),polarcap_eflux)
        !              endif
        !           enddo
        !        enddo

     else

        iono_south_eflux = 1.0e-6
        iono_south_ave_e = 1.0

        do j = 1, IONO_nPsi
           IsPeakFound = .false.
           MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)

           AuroraWidth = -1.0
           IsDone = .false.
           i = IONO_nTheta
           OCFLB = -1.0
           do while (.not. IsDone)

              if (iono_south_p(i,IONO_nPsi-j+1) > 0) then

                 if (OCFLB(j) == -1) OCFLB(j) = iono_south_theta(i,j)
                 AuroraWidth = OCFLB(j) - iono_south_theta(i,j)
                 MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)
                 iono_south_eflux(i,j) = MulFac_ef*MaxP

                 if (iono_south_p(i,IONO_nPsi-j+1)==MaxP) IsPeakFound = .true.

                 if (IsPeakFound .and. AuroraWidth >= MinWidth) then
                    EquatorWardEdge(j) = iono_south_theta(i,j)
                    IsDone = .true.
                 endif

                 if (i == 1) IsDone = .true.

              endif

              i = i - 1

           enddo

           if (.not. IsPeakFound) then
              OCFLB(j) = cPi - MinWidth
              EquatorwardEdge(j) = cPi - MinWidth*2
           endif

        enddo


        smooth = 0.0
        do j = 1, IONO_nPsi
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + OCFLB(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        OCFLB = smooth/(nHalfSmooth*2+1)

        smooth = 0.0
        do j = 1, IONO_nPsi
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + &
                   EquatorWardEdge(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        EquatorWardEdge = smooth/(nHalfSmooth*2+1)

        do j = 1, IONO_nPsi
           Center = (OCFLB(j)*2.0 + EquatorWardEdge(j))/3.0
           Width = abs(EquatorWardEdge(j) - OCFLB(j))/2

           MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)
           iono_south_eflux(:,j) = &
                MulFac_ef*MaxP * &
                exp(-abs(Center-iono_south_theta(:,j))/Width) * &
                (0.375*cos(iono_south_psi(:,j)-cPi)+0.625)

           iono_south_ave_e(:,j) = &
                iono_south_p(:,IONO_nPsi-j+1) / &
                (iono_south_rho(:,IONO_nPsi-j+1)+1.0e-32) * MulFac_ae

        enddo

        diffuse_ef = iono_south_eflux
        diffuse_ae = iono_south_ave_e

        discrete_ef = 0.0
        discrete_k  = 0.0
        ! This is from Jimmy's Paper on the Knight Relationship
        where (iono_south_p > 0) discrete_k = &
             (iono_south_rho**1.5) / iono_south_p
        ! Reverse the particles again.
        do j = 1, IONO_nPsi
           discrete_k(:,j) = discrete_k(:,IONO_nPsi-j+1)
        enddo

        where (iono_south_jr > 0.5e-7) &
             discrete_ef = (iono_south_jr*1e6)*discrete_k
        discrete_ae = discrete_ef*MulFac_Dae
        discrete_ef = (iono_south_jr*1e6)*discrete_ef*MulFac_Def

        ! Let's add a little conductance on ANY not in the polar cap,
        ! so the code doesn't blow up
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              if (abs(iono_south_jr(i,j)) > 0.0 .and. &
                   discrete_ef(i,j) < 5.0e-3) then
                 discrete_ef(i,j) = MulFac_Def/2.5 * &
                      (iono_south_jr(i,j)*1e6)**2*discrete_k(i,j)
              endif
           enddo
        enddo

        where (diffuse_ae < IONO_Min_Ave_E/2) diffuse_ae = IONO_Min_Ave_E/2
        where (discrete_ae < IONO_Min_Ave_E/2) discrete_ae = IONO_Min_Ave_E/2

        !        write(*,*) "Discrete (south) : ",&
        !             minval(discrete_ef), maxval(discrete_ef),&
        !             minval(discrete_ae), maxval(discrete_ae)

        ! Let's weight the average energy by the number flux, which is ef/av
        iono_south_ave_e = &
             (diffuse_ef + discrete_ef) / &
             (diffuse_ef/diffuse_ae + discrete_ef/discrete_ae)

        ! The energy flux should be weighted by the average energy
        iono_south_eflux = &
             (diffuse_ef/diffuse_ae + discrete_ef/discrete_ae) * &
             iono_south_ave_e

        where (iono_south_ave_e < IONO_Min_Ave_E) &
             iono_south_ave_e = IONO_Min_Ave_E

        write(*,*) "MAGNIT: Done with northern aurora"

        !        do j = 1, IONO_nPsi
        !           do i = 1, IONO_nTheta
        !              if (iono_south_p(i,j) <= 0.0) then
        !                 iono_south_ave_e(i,j) = max(iono_south_ave_e(i,j),polarcap_avee)
        !                 iono_south_eflux(i,j) = &
        !                      max(iono_south_eflux(i,j),polarcap_eflux)
        !              endif
        !           enddo
        !        enddo

     endif

  endif

  if (iModel.eq.8) then

     !=================================================================
     !-----------------------------------------------------------------
     !------------- MAGNetosphere-Ionosphere-Thermosphere -------------
     !----------- (MAGNIT) Physics-based Conductance Model ------------
     !-----------------------------------------------------------------
     !------------ Developed by Agnit Mukhopadhyay (2019) -------------
     !------------   University of Michigan, Ann Arbor    -------------
     !-----------------------------------------------------------------
     !=================================================================
     
     ! Multiplication Factors - For OCFLB Estimation
     MinWidth = 5.0 * cPi / 180.0
     nHalfSmooth = 5
!     MulFac_Dae = 1.0e22
!     MulFac_Def = 5.0e19
!     MulFac_ef = 0.2e7
!     MulFac_ae = 1.0 / 1.0e11
     MulFac_Dae = 1.0e22
     MulFac_Def = 1.0e19
     MulFac_ef = 0.3e6! * 1e03
     MulFac_ae = 3.65e-12

     ! Warning Message
     write(*,*) '################## CAUTION ##################'
     write(*,*) '(Beta Testing) The conductance calculations from MAGNIT are unstable.'
     write(*,*) 'Please proceed with caution. For issues, please contact developers.'
     write(*,*) '#############################################'

     if (iBlock == 1) then ! Northern Hemisphere

        ! First, we set up some basic limits to our flux values
        iono_north_eflux = 1.0e-6
        iono_north_ave_e = 1.0

        ! Next, we describe and calculate the
        ! Open Closed Field Line Boundary or
        ! OCFLB, and the EquatorWard Edge
        OCFLB = -1.0
        EquatorWardEdge = -1.0

        do j = 1, IONO_nPsi ! MLT-wise only
           IsPeakFound = .false.
           MaxP = max(maxval(iono_north_p(:,IONO_nPsi-j+1)),1.0e-15)

           ! OCFLB Calculation was set up in an early version
           ! of Conductance Model No. 6
           IsDone = .false.
           i = 1
           do while (.not. IsDone)
              
              if (iono_north_p(i,IONO_nPsi-j+1) > 0) then

                 if (OCFLB(j) == -1) OCFLB(j) = iono_north_theta(i,j)
                 AuroraWidth = iono_north_theta(i,j) - OCFLB(j)
                 iono_north_eflux(i,j) = MulFac_ef*MaxP

                 if (iono_north_p(i,IONO_nPsi-j+1)==MaxP) IsPeakFound = .true.

                 if (IsPeakFound .and. AuroraWidth >= MinWidth) then
                    EquatorWardEdge(j) = iono_north_theta(i,j)
                    IsDone = .true.
                 endif

                 if (i == IONO_nTheta) IsDone = .true.

              endif

              i = i + 1

           enddo
           if (.not. IsPeakFound) then
              OCFLB(j) = MinWidth
              EquatorwardEdge(j) = MinWidth*2
           endif

        enddo

        do j = 1, IONO_nPsi
           smooth(j) = 0.0
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + OCFLB(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        OCFLB = smooth/(nHalfSmooth*2+1)
           
        do j = 1, IONO_nPsi
           smooth(j) = 0.0
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + &
                   EquatorWardEdge(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        EquatorWardEdge = smooth/(nHalfSmooth*2+1)

        ! Set lower limit for MHD Density
        where (iono_north_rho < 1E-21) &
          iono_north_rho = 1E-21

        ! Calculation of Initial Diffuse Flux
        do j = 1, IONO_nPsi

           Center = (OCFLB(j)*2.0 + EquatorWardEdge(j))/3.0
           Width = abs(EquatorWardEdge(j) - OCFLB(j))/2

           MaxP = max(maxval(iono_north_p(:,IONO_nPsi-j+1)),1.0e-15)

           !write(*,*), "MaxP(j) - N", MaxP, j

           !write(*,*) 'OCFLB, EquatorEdge North', iono_north_psi(1,j)*cRadToDeg, &
           !     iono_north_theta(1,j)*cRadToDeg, &
           !     90. - OCFLB(j)*cRadToDeg, 90. - EquatorWardEdge(j)*cRadToDeg

           ! Energy Flux calculation in Conductance Model No. 6
           !iono_north_eflux(:,j) = &
           !     MulFac_ef*MaxP * &
           !    exp(-abs(Center-iono_north_theta(:,j))/Width) * &
           !     (0.375*cos(iono_north_psi(:,j)-cPi)+0.625)
           
           !write(*,*) 'North - Exp 1', exp(-abs(Center-iono_north_theta(10,j))/Width)

           ! Avg Energy calculation in Conductance Model No. 6
           ! AJR - This is the way that seems to work
           iono_north_ave_e(:,j) = &
                iono_north_p(:,IONO_nPsi-j+1) / &
                (iono_north_rho(:,IONO_nPsi-j+1)+1.0e-32) * MulFac_ae

!           iono_north_ave_e(Poleward-2,j) = iono_north_ave_e(Poleward+1,j)
!           iono_north_ave_e(Poleward-1,j) = iono_north_ave_e(Poleward+1,j)

!           do i = Poleward-2, Equatorward+2
!           do i = Poleward, Equatorward
!              f = exp(-abs(float(i-Center))/float(width))
!              iono_north_eflux(i,j) = iono_north_eflux(i,j) * f
!           enddo

        enddo

        ! ================= SOURCE ==================
        ! Electron Diffuse Precipitation Calculation:
        ! ===========================================
        
        !=================================================================
        !-----------------------------------------------------------------
        !--------------------- MAXWELLIAN DISTRIBUTION -------------------
        !----------- Using same technique as in Zhang et al. 2015, -------
        !----------- Raeder et al. 2001, Wiltberger et al. 2009    -------
        ! Implemented as a test condition --------------------------------
        !-----------------------------------------------------------------
        !=================================================================

        ! We already have a good estimate of the average energy based on
        ! the MHD temperature. In the following section, we find E-flux
        ! again!
        
        ! Number Flux = N_e * ((kB * T)/(2pi * m_elec))**0.5
        ! N_e = Density conversion (amu/cm3) = 1.66E-21 kg/m3
        ! T = Ave_E * 1E03 * 11000 (conversion from keV to Kelvins)
        ! m_elec = 9.1E-31 kg
        ! kB = 1.38E-23 J/K
        ! Multiplication Factor = (kB/(2pi * m_elec))**0.5
        !                       = 1553.5632

        ! Energy Flux = Number Flux * Average Energy

        do j = 1, IONO_nPsi
           ! Conversion to amu/cm3
           iono_Ne(:,j) = iono_north_rho(:,IONO_nPsi-j+1) / 1.66e-21
        enddo

        !write(*,*) 'MaxVal Number Density (in amu/cm3)', maxval(iono_Ne)
        ! Conversion to Kelvins
        iono_T = iono_north_ave_e * 1e03 * 1.1e04! * 1./6.
        ! Energy Flux
        iono_north_eflux = iono_Ne * 1.38e-23 * 1e06 * &
             1553.5632  * iono_T**1.5 !* &
             !0.4 ! Beta = Loss Cone Factor (Can be replaced using same method as in RCM)

        if (DoTest) then
           write(*,*) 'N_Diff: MaxVal Ave_E (in keV)', maxval(iono_north_ave_e)
           !write(*,*) 'MaxVal Temperature (in K)', maxval(iono_T)
           write(*,*) 'N_Diff: MaxVal EFlux (in W/m2)', maxval(iono_north_eflux)
        end if

        ! Loss Cone Factor, as derived in Wolf et al. 1991
        !do j = 1, IONO_nPsi
        !   ! Loss Cone Fraction as a trigonometric relation of MLT
        !   beta = 0.5 + 0.3 * sin(iono_north_psi(:,j))
        !   iono_north_eflux(:,j) = iono_north_eflux(:,j) * beta
        !   !write(*,*) 'LC_frac Psi North', maxval(beta), maxval(iono_north_psi(:,j) * 12.0/cPi)
        !end do

        !==================== Limits ========================
        !-----------------------------------------------------------------
        !-----------------------------------------------------------------
        
        
        ! At this point, we will add EFlux and AveE from IM in addition to
        ! the fluxes calculated from GM. 
        !
        ! Variables:
        ! ==========
        !
        ! IM EFlux              iono_north_im_eflux/1000.0
        ! IM AveE               iono_north_im_avee
        !
        ! GM EFlux              iono_north_eflux
        ! GM AveE               iono_north_ave_e

        ! Limits on IM & GM Precip - Diffuse
        where (iono_north_im_avee > 20.0) &
          iono_north_im_avee = 20.0
        where (iono_north_im_avee < IONO_Min_Ave_E) &
             iono_north_im_avee = IONO_Min_Ave_E
        where (iono_north_im_eflux/1000.0 < IONO_Min_EFlux) &
             iono_north_im_eflux = IONO_Min_EFlux * 1000.0
        where (iono_north_ave_e < IONO_Min_Ave_E) &
             iono_north_ave_e = IONO_Min_Ave_E
        where (iono_north_eflux < IONO_Min_EFlux) &
             iono_north_eflux = IONO_Min_EFlux
        where (iono_north_eflux < IONO_Min_EFlux) &
             iono_north_ave_e = IONO_Min_Ave_E
        where (iono_north_im_eflux/1000.0 < IONO_Min_EFlux) &
             iono_north_im_avee = IONO_Min_Ave_E
        !where (iono_north_im_eflux/1000.0 > IONO_Min_EFlux) &
        !     iono_north_eflux = iono_north_im_eflux/1000.0
        !where (iono_north_im_avee > IONO_Min_Ave_E) &
        !     iono_north_ave_e = iono_north_im_avee
        
        if(iModel.eq.8) then
           ! Precipitation from only GM
           diffuse_ef = iono_north_eflux
           diffuse_ae = iono_north_ave_e

           !do j = 1, IONO_nPsi
           !   do i = 1, IONO_nTheta
           !      if (iono_north_im_eflux(i,j)/1000.0 > IONO_Min_EFlux) then
           !         diffuse_ef(i,j) = iono_north_im_eflux(i,j)/1000.0
           !         diffuse_ae(i,j) = iono_north_im_avee(i,j)
           !      end if
           !      !where (iono_north_im_avee > IONO_Min_Ave_E) &
           !      !     iono_north_ave_e = iono_north_im_avee
           !   end do
           !end do
        !else
           ! Precipitation from both GM and IM
           !diffuse_ae = &!iono_north_im_avee
           
           !diffuse_ef = iono_north_im_eflux
           !diffuse_ae = iono_north_im_avee
           
           !diffuse_ae = (iono_north_eflux + iono_north_im_eflux/1000.0)/ &
           !     (iono_north_eflux/iono_north_ave_e + &
           !     iono_north_im_eflux/(iono_north_im_avee*1000.0))

           !diffuse_ef = iono_north_im_eflux/1000.0
           !diffuse_ef = (iono_north_eflux/iono_north_ave_e + &
           !     iono_north_im_eflux/(iono_north_im_avee*1000.0)) * diffuse_ae
        end if

        !write(*,*) "Max EFlux = ", maxval(diffuse_ef)
        
        !================== Discrete Conductance =========================
        
        !discrete_ef = 0.0
        !discrete_k  = 0.0

        ! =================== SOURCE ====================
        ! Electron Monoenergetic Precipitation (Discrete)
        ! ===============================================
        !
        ! Following similar procedure as in Knight (1973) and
        ! Fridman and Lemaire (1980). Re-done in Zhang et al.
        ! (2015) and Yu et al. (2016)
        !
        
        diffuse_nf = iono_Ne * 1e06 * 1553.5632 * iono_T**0.5 ! * 0.15
        rm = 1.1 !+ 3*SIN(cHalfPi - IONO_NORTH_Theta(i,j))!100

        ! Let's have some default values
        ev = 0.0
        numflux_floor = cFACFloor/cElectronCharge
        discrete_nf = numflux_floor

        ! The Great Loop of our time...
        do j = 1, IONO_nPsi
           ! IONO_nTheta is the equator and the formulas below break down
           do i = 1, IONO_nTheta - 1
              !iono_north_jr(i,j) = iono_north_jr(i,j) * 10.

              ! Trigonometry tells us that sin(pi/2 - angle) = cos(angle)
              ! For the equator this gives rm=1
              rm = sqrt(1 + 3*cos(IONO_NORTH_Theta(i,j))) ! Assume Dipolar
              !----------------- DISCRETE NUMBER FLUX -------------------

              ! Discrete N_flux = J/e
              ! Consider only upward FACs
              if (iono_north_jr(i,j) > 0.) then
                 ! positive numbers don't change with absolute value
                 ! constants like 1.6e-19 should be named constants
                 discrete_nf(i,j) = ABS(iono_north_jr(i,j)) / cElectronCharge
              else
                 discrete_nf(i,j) = numflux_floor
              end if
              !if(i==2 .and. j==1)then
              !   write(*,*)'!!! iono_north_jr(i,j)=',iono_north_jr(i,j)
              !   write(*,*)'!!! discrete_nf(i,j) = ',discrete_nf(i,j)
              !end if

              ! Set a minimum on diffuse number flux
              ! maximum function can be used above to avoid this
              if (diffuse_nf(i,j) < numflux_floor) then
                 diffuse_nf(i,j) = numflux_floor
              end if

              ! Set a minimum on electron temperature = 1 eV.
              if (iono_T(i,j) < cEVtoK) then
                 iono_T(i,j) = cEVtoK
              end if              

              ! Deal with the following portion carefully
              ! Lots of avenues to get non-physical answers
              if (iono_north_jr(i,j) > 0.) then
                 ! Variable must not be greater than
                 ! Why ABS. Why constants.
                 vari = ABS(iono_north_jr(i,j)) / &
                      (1.6E-19* diffuse_nf(i,j))
                 ! Why 5e-5? Why 05? Why then ... endif. Why vari is changed? Why not rm?
                 if ((rm - vari) < 0.) then
                    vari = rm - 5E-05 ! Try
                 end if
                 
                 ! Potential Energy eV
                 ! This formula is later inverted back to (rm-1)/(rm-vari). This one fails from rm=1.
                 eV(i,j)  = 1.38E-23 * iono_T(i,j) * (rm - 1) * &
                      LOG((rm - 1)/(rm - vari))
              else
                 ! This is unneeded if ev is initialized
                 eV(i,j) = 0.
              end if


              ! Grand Calculation for the Energy Flux (in W/m2)
              ! See Knight [1973] and Fridman and Lemaire [1980] for details.

              !write(*,*)'!!! North', i,j,ev(i,j),iono_T(i,j),rm

              ! Dummy vars do not exist. All variables have meaning. 
              ! This is now reverting the above formula, so var_rm = (rm - 1)/(rm - vari) or 1
              var_rm = EXP(-eV(i,j)/(1.38e-23 * iono_T(i,j) * &  !
                   (rm - 1)))                                    !
              numerator   = 1 - var_rm                           ! DUMMY VARS

              denominator = 1 + (1 - 1/rm) * var_rm              !
              if (denominator == 0.) then
                 denominator = 1E-15
              end if
              var_tot = eV(i,j) * numerator/denominator          ! 
              
              ! Hard coded constants are evil. 
              discrete_ef(i,j) = (var_tot + 2*1.38e-23*iono_T(i,j)) * &
                   discrete_nf(i,j)

              !write(*,*) eV(i,j), vari
              
              ! Average Energy can be estimated as,
              ! Ave_E = E_Flux/N_Flux (in eV)

              discrete_ae(i,j) = discrete_ef(i,j)/&!discrete_nf(i,j)
                   (1.6e-19 * discrete_nf(i,j) )
              
              ! Convert EFlux and Ave_E into units for Robinson relation

              ! This line does nothing
              discrete_ef(i,j) = discrete_ef(i,j) !* 10.0! * 1E03 ! Convert to ergs/cm2
              ! This does something
              discrete_ae(i,j) = discrete_ae(i,j) / 1E03 ! Convert to keV

              !iono_north_jr(i,j) = iono_north_jr(i,j) / 10.
           enddo
        enddo

        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta

              ! The condition is on theta only. So there is no need to loop in Psi.
              ! An or statement can be used. Also, what was the point of the effort above
              ! if now everything is overwritten below 60 degrees latitude?
              ! Equatorward Boundary Limit
              if ((IONO_NORTH_Theta(i,j)*cRadToDeg) >= 40.) then
                 !if (iono_north_jr(i,j) <= 1E-09) then 
                 discrete_ae(i,j) = IONO_Min_Ave_E
                 discrete_ef(i,j) = IONO_Min_EFlux
                 diffuse_ae(i,j) = IONO_Min_Ave_E
                 diffuse_ef(i,j) = IONO_Min_EFlux
              end if

              ! Open-Closed Boundary Limit
              ! Why convert to degrees. Radians can be compared just as well.
              if ((IONO_NORTH_Theta(i,j)*cRadToDeg) <= OCFLB(j)*cRadToDeg) then
                 !if (iono_north_jr(i,j) <= 1E-09) then 
                 discrete_ae(i,j) = IONO_Min_Ave_E
                 discrete_ef(i,j) = IONO_Min_EFlux
                 diffuse_ae(i,j) = IONO_Min_Ave_E
                 diffuse_ef(i,j) = IONO_Min_EFlux
              end if

           end do
        end do

        ! What about the min function?
        where (iono_north_ave_e > 100.0) &
             iono_north_ave_e = 100.0

        if (DoTest) then
           write(*,*) 'N_Discr: MinVal Downward FAC', minval(iono_north_jr)
           write(*,*) 'N_Discr: MaxVal NFlux ', maxval(discrete_nf(1:IONO_ntheta-1,:))
           write(*,*) 'N_Discr: MaxVal EFlux (in W/m2)', maxval(discrete_ef)
           write(*,*) 'N_Discr: MaxVal Ave_E (in keV)', maxval(discrete_ae)
        end if
        
        ! Let's add a little conductance on ANY not in the polar cap,
        ! so the code doesn't blow up
        !do j = 1, IONO_nPsi
        !   do i = 1, IONO_nTheta
        !      if (abs(iono_north_jr(i,j)) > 0.0 .and. &
        !           discrete_ef(i,j) < 5.0e-3) then
        !         discrete_ef(i,j) = MulFac_Def/2.5 * &
        !              (iono_north_jr(i,j)*1e6)**2*discrete_k(i,j)
        !      endif
        !   enddo
        !enddo

        where (diffuse_ae < IONO_Min_Ave_E/2) diffuse_ae = IONO_Min_Ave_E/2
        where (discrete_ae < IONO_Min_Ave_E/2) discrete_ae = IONO_Min_Ave_E/2

        ! High and Low Ceilings - discrete
        where (discrete_ae < IONO_Min_Ave_E) &
             discrete_ae = IONO_Min_Ave_E
        where (discrete_ef < IONO_Min_EFlux) &
             discrete_ef = IONO_Min_EFlux
        where (discrete_ef < IONO_Min_EFlux) &
             discrete_ae = IONO_Min_Ave_E
        where (iono_north_eflux < IONO_Min_EFlux) &
             iono_north_ave_e = IONO_Min_Ave_E

        where (diffuse_ae > 50.0) &
             diffuse_ae = IONO_Min_Ave_E
        where (discrete_ae > 50.0) &
             discrete_ae = IONO_Min_Ave_E
        
        !where (iono_north_im_eflux/1000.0 > IONO_Min_EFlux) &
        !     diffuse_ef = iono_north_im_eflux/1000.0
        !where (iono_north_im_eflux/1000.0 > IONO_Min_EFlux) &
        !     diffuse_ae = iono_north_im_avee

        ! -----------------------------------------
        ! Smoothing Filter for Discrete Energy Flux
        ! -----------------------------------------
        !
        ! During our simulations, it appears that a combination of
        ! boundary issues coupled with the discrete flux calculations
        ! leads to localized pockets of high FACs with unphysical
        ! values as the simulation marches in time. This step is
        ! added to make sure such localization is arrested early on.

        ! Smoothing Filter on Boundaries - discrete_ef
        call smooth_lagrange_polar(discrete_ef, IONO_nTheta, IONO_nPsi)

        ! Smoothing Filter on Boundaries - discrete_nf
        call smooth_lagrange_polar(discrete_nf, IONO_nTheta, IONO_nPsi)

        ! -----------------------------------------------------
        ! Calculate net electron fluxes =======================
        ! -----------------------------------------------------

        ! NOTE: This calculation is only for diffuse and discrete
        ! electron precip. For ion diffuse and broadband, the flux
        ! is converted into conductance differently. - AM

        if(DoTest) then ! Debug
           write(*,*)'IONO_nPsi and nTheta = ', IONO_nPsi, IONO_nTheta
           write(*,*)'Let us look at diffuse_ef, discrete_ef, diffuse_nf, discrete_nf'
           do j = 1, IONO_nPsi
              do i = 1, IONO_nTheta
                 write(*,'(2(i4.4,1x),4(E10.5, 1x))') i,j,diffuse_ef(i,j), &
                      discrete_ef(i,j), diffuse_nf(i,j), discrete_nf(i,j)
              end do
           end do
        end if

        ! Energy Flux ================================================
        ! Weight the Average Energy by the Number Flux, which is ef/av
        
        ! These will be parameterized in the future... - AM & DTW
        !iono_north_ave_e = &            ! If only using diffuse precip
        !     diffuse_ef /  & 
        !     (1E03 * 1.6e-19 * diffuse_nf)
        !iono_north_ave_e = discrete_ae  ! If only using discrete precip
        
        iono_north_ave_e =                & ! If using both precip
             (diffuse_ef + discrete_ef) / &
             (1E03 * 1.6e-19 *            & ! Convert this into keV (x 10^3)
             (diffuse_nf + discrete_nf))

        ! Average Energy =========================================
        ! The energy flux should be weighted by the average energy

        ! These will be parameterized in the future... - AM & DTW
        !iono_north_eflux = diffuse_ef      ! If only using diffuse prcip
        !iono_north_eflux = discrete_ef     ! If only using discrete precip
        iono_north_eflux =                 & ! If using both precip
             (1E03 * 1.6e-19 *             &
             (diffuse_nf + discrete_nf)) * &
             iono_north_ave_e

        ! High and Low Ceilings - From iModel 6
        where (iono_north_ave_e < IONO_Min_Ave_E) &
             iono_north_ave_e = IONO_Min_Ave_E
        where (iono_north_ave_e < IONO_Min_Ave_E) &
             iono_north_ave_e = IONO_Min_Ave_E
        where (iono_north_eflux < IONO_Min_EFlux) &
             iono_north_eflux = IONO_Min_EFlux
        where (iono_north_eflux < IONO_Min_EFlux) &
             iono_north_ave_e = IONO_Min_Ave_E
        where (iono_north_eflux < IONO_Min_EFlux) &
             iono_north_ave_e = IONO_Min_Ave_E

        !where ((iono_north_im_eflux/1000.0) > IONO_Min_Eflux) &
        !     diffuse_ef = iono_north_im_eflux/1000.0
        !where ((iono_north_im_eflux/1000.0) > IONO_Min_Eflux) &
        !     diffuse_ae = iono_north_im_avee

        ! The Great Loop of our time...
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              ! Equatorward Boundary Limit
              if ((IONO_NORTH_Theta(i,j)*cRadToDeg) >= 40.) then
                 !if (iono_north_jr(i,j) <= 1E-09) then 
                 iono_north_ave_e(i,j) = IONO_Min_Ave_E
                 iono_north_eflux(i,j) = IONO_Min_EFlux
              end if

              ! Open-Closed Boundary Limit
              if ((IONO_NORTH_Theta(i,j)*cRadToDeg) <= OCFLB(j)*cRadToDeg) then
                 !if (iono_north_jr(i,j) <= 1E-09) then 
                 iono_north_ave_e(i,j) = IONO_Min_Ave_E
                 iono_north_eflux(i,j) = IONO_Min_EFlux
              end if
              
           end do
        end do

        where (iono_north_ave_e > 10.0) &
             iono_north_ave_e = IONO_Min_Ave_E

        ! For Storage onto global variables
        IONO_NORTH_DIFF_Ave_E = diffuse_ae
        IONO_NORTH_DIFF_EFlux = diffuse_ef

        IONO_NORTH_MONO_Ave_E = discrete_ae
        IONO_NORTH_MONO_EFlux = discrete_ef

        ! Smoothing Filter on Boundaries - Total Electron E-Flux
        call smooth_lagrange_polar(iono_north_eflux, IONO_nTheta, IONO_nPsi)

        ! Smoothing Filter on Boundaries - Total Electron Ave-E
        call smooth_lagrange_polar(iono_north_ave_e, IONO_nTheta, IONO_nPsi)
        
        if(DoTest) write(*,*) "MAGNIT: Done with northern aurora."
        
        !! Let's add a little conductance on ANY not in the polar cap,
        !! so the code doesn't blow up
        !do j = 1, IONO_nPsi
        !   do i = 1, IONO_nTheta
        !      if (iono_north_p(i,j) <= 0.0) then
        !         iono_north_ave_e(i,j) = max(iono_north_ave_e(i,j),polarcap_avee)
        !         iono_north_eflux(i,j) = &
        !              max(iono_north_eflux(i,j),polarcap_eflux)
        !      endif
        !   enddo
        !enddo

     else !----------------- SOUTHERN HEMISPHERE --------------------------

        !iono_south_eflux = 1.0e-6
        !iono_south_ave_e = 1.0

        ! The OCFLB Calculation is a little different in the
        ! Southern Hemisphere (because of how the colatitudes
        ! work), so be careful of doing this calculation separately
        ! and not using the solution from the Northern Hemisphere.
        
        do j = 1, IONO_nPsi
           IsPeakFound = .false.
           MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)
        
           AuroraWidth = -1.0
           IsDone = .false.
           i = IONO_nTheta
           OCFLB_S = -1.0
           do while (.not. IsDone)
        
              if (iono_south_p(i,IONO_nPsi-j+1) > 0) then
        
                 if (OCFLB_S(j) == -1) then
                    !OCFLB_temp(j) = cPi - iono_south_theta(i,j)
                    OCFLB_S(j) = iono_south_theta(i,j)
                 endif
                 
                 AuroraWidth = OCFLB_S(j) - iono_south_theta(i,j)
                 !MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)
                 iono_south_eflux(i,j) = MulFac_ef*MaxP
        
                 if (iono_south_p(i,IONO_nPsi-j+1)==MaxP) IsPeakFound = .true.
        
                 if (IsPeakFound .and. AuroraWidth >= MinWidth) then
                    EquatorWardEdge(j) = cPi - iono_south_theta(i,j)
                    IsDone = .true.
                 endif
        
                 if (i == 1) IsDone = .true.
        
              endif
              
              i = i - 1
              
           enddo
           
           if (.not. IsPeakFound) then
              OCFLB_S(j) = cPi - MinWidth
              EquatorwardEdge(j) = cPi - MinWidth*2
           endif

           OCFLB(j) = OCFLB_S(j)
           
        enddo

        !write(*,*) OCFLB
        
        !do j = 1, IONO_nPsi
        !   write(*,*) j, OCFLB(j)
        !enddo
        

        smooth = 0.0
        do j = 1, IONO_nPsi
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + OCFLB(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        OCFLB = smooth/(nHalfSmooth*2+1)
           
        smooth = 0.0
        do j = 1, IONO_nPsi
           do i = j-nHalfSmooth-1, j+nHalfSmooth-1
              smooth(j) = smooth(j) + &
                   EquatorWardEdge(mod(i+IONO_nPsi,IONO_nPsi)+1)
           enddo
        enddo
        EquatorWardEdge = smooth/(nHalfSmooth*2+1)

        !where (iono_south_rho < 1E-21) &
        !     iono_south_rho = 1E-21
        !where (iono_south_p < (IONO_Min_Ave_E/(1E-21 * MulFac_ae))) &
        !     iono_south_p = (IONO_Min_Ave_E/(1E-21 * MulFac_ae))

        !write(*,*) 'Minimum Rho, P = ', 1E-21, (IONO_Min_Ave_E/(1E-21 * MulFac_ae))

        ! Conductance Model No. 6 - Calculation of Energy Flux and Avg. Energy
        do j = 1, IONO_nPsi
           Center = (OCFLB(j)*2.0 + EquatorWardEdge(j))/3.0
           Width = abs(EquatorWardEdge(j) - OCFLB(j))/2

           MaxP = max(maxval(iono_south_p(:,IONO_nPsi-j+1)),1.0e-15)

           !write(*,*) 'OCFLB, EquatorEdge South', OCFLB(j), EquatorWardEdge(j)
           !iono_south_eflux(:,j) = &
           !     MulFac_ef*MaxP * &
           !     exp(-abs(Center-iono_south_theta(:,j))/Width) * &
           !     (0.375*cos(iono_south_psi(:,j)-cPi)+0.625)
         
           iono_south_ave_e(:,j) = &
                iono_south_p(:,IONO_nPsi-j+1) / &
                (iono_south_rho(:,IONO_nPsi-j+1)+1.0e-32) * MulFac_ae

        enddo

        ! ================= SOURCE ==================
        ! Electron Diffuse Precipitation Calculation:
        ! ===========================================
        
        !-----------------------------------------------------------------
        !----------- Using same technique as in Zhang et al. 2015, -------
        !----------- Raeder et al. 2001, Wiltberger et al. 2009    -------
        ! Implemented as a test condition --------------------------------
        !-----------------------------------------------------------------
        
        !
        ! We already have a good estimate of the average energy based on
        ! the MHD temperature. In the following section, we find E-flux
        ! again!
        
        ! Number Flux = N_e * ((kB * T)/(2pi * m_elec))**0.5
        ! N_e = Electron density (amu/m3) = 1.66E-27 kg/m3
        ! T = Ave_E * 1E03 * 11000 (conversion from keV to Kelvins)
        ! m_elec = 9.1E-31 kg
        ! kB = 1.38E-23 J/K
        ! Multiplication Factor = (kB/(2pi * m_elec))**0.5
        !                       = 1553.5632

        ! Energy Flux = Number Flux * Average Energy

        do j = 1, IONO_nPsi
           ! Conversion to amu/cm3
           iono_Ne(:,j) = iono_south_rho(:,IONO_nPsi-j+1) / 1.66e-21
        enddo
        
        !write(*,*) 'MaxVal Number Density (in amu/cm3)', maxval(iono_Ne)
        ! Conversion to Kelvins
        iono_T = iono_south_ave_e * 1e03 * 1.1e04! * 1./6.

        !write(*,*) 'MaxVal Temperature (in K)', maxval(iono_T)
        ! Energy Flux

        iono_south_eflux = iono_Ne * 1.38e-23 * 1e06 * &
             1553.5632  * iono_T**1.5 !* &
             !0.4 ! Beta = Loss Cone Factor (Can be replaced using same method as in RCM)

        if (DoTest) then
           write(*,*) 'S_Diff: MaxVal Ave_E (in keV)', maxval(iono_south_ave_e)
           write(*,*) 'S_Diff: MaxVal EFlux (in W/m2)', maxval(iono_south_eflux)
        end if
        !-----------------------------------------------------------------
        !-----------------------------------------------------------------

        
        ! At this point, we will add EFlux and AveE from IM in addition to
        ! the fluxes calculated from GM. 
        !
        ! Variables:
        ! ==========
        !
        ! IM EFlux              iono_north_im_eflux/1000.0 (values in keV)
        ! IM AveE               iono_north_im_avee
        !
        ! GM EFlux              iono_north_eflux
        ! GM AveE               iono_north_ave_e

        ! Limits on IM & GM Precip
        where (iono_south_im_avee > 20.0) &
          iono_south_im_avee = 20.0
        where (iono_south_im_avee < IONO_Min_Ave_E) &
             iono_south_im_avee = IONO_Min_Ave_E
        where (iono_south_im_eflux < IONO_Min_EFlux) &
             iono_south_im_eflux = IONO_Min_EFlux
        where (iono_south_ave_e < IONO_Min_Ave_E) &
             iono_south_ave_e = IONO_Min_Ave_E
        where (iono_south_eflux < IONO_Min_EFlux) &
             iono_south_eflux = IONO_Min_EFlux
        where (iono_south_eflux < IONO_Min_EFlux) &
             iono_south_ave_e = IONO_Min_Ave_E
        where (iono_south_im_eflux < IONO_Min_EFlux) &
             iono_south_im_avee = IONO_Min_Ave_E
        
        if(iModel.eq.8) then
        !   ! Precipitation from only GM
        !   diffuse_ef = iono_south_eflux
        !   diffuse_ae = iono_south_ave_e
        !else
        !   ! Precipitation from both GM and IM
           !where (iono_south_im_eflux/1000.0 > IONO_Min_EFlux) &
           !     iono_south_eflux = iono_south_im_eflux/1000.0
           !where (iono_south_im_avee > IONO_Min_Ave_E) &
           !     iono_south_ave_e = iono_south_im_avee
           diffuse_ef = iono_south_eflux
           diffuse_ae = iono_south_ave_e

           !diffuse_ef = iono_south_im_eflux
           !diffuse_ae = iono_south_im_avee
           
           !diffuse_ae = &!iono_south_im_avee
           !     (iono_south_eflux + iono_south_im_eflux/1000.0)/ &
           !     (iono_south_eflux/iono_south_ave_e + &
           !     iono_south_im_eflux/(iono_south_im_avee*1000.0))
           
           !diffuse_ef = &!iono_south_im_eflux/1000.0
           !     (iono_south_eflux/iono_south_ave_e + &
           !     iono_south_im_eflux/(iono_south_im_avee*1000.0)) * diffuse_ae
        end if

        !================== Discrete Conductance =========================
        
        discrete_ef = 1e-12
        !discrete_k  = 0.0
        discrete_nf = 1e-8
        discrete_ae = 1e-12

        ! =================== SOURCE ====================
        ! Electron Monoenergetic Precipitation (Discrete)
        ! ===============================================
        !
        ! Following similar procedure as in Knight (1973) and
        ! Fridman and Lemaire (1980). Re-done in Zhang et al.
        ! (2015) and Yu et al. (2016)
        !

        ! DISCRETE PRECIP
        diffuse_nf = iono_Ne * 1e06 * 1553.5632 * iono_T**0.5 ! * 0.15 
        rm = 1.1

        ! Set some floorline values
        ev = 0.0
        numflux_floor = cFACFloor/cElectronCharge
        discrete_nf = numflux_floor
        discrete_ae = IONO_Min_Ave_E
        discrete_ef = IONO_Min_EFlux
        
        ! The Great Loop of our time...
        do j = 1, IONO_nPsi
           ! IONO_nTheta is the equator and the formulas below break down
           do i = 1, IONO_nTheta - 1
              !iono_north_jr(i,j) = iono_north_jr(i,j) * 10.
              !(180. - IONO_SOUTH_Theta(i,j)*cRadToDeg)*cPi/180.
              rm = sqrt(1 - 3*cos(IONO_SOUTH_Theta(i,j)))
              !IONO_SOUTH_Theta(i,j))
              !----------------- DISCRETE NUMBER FLUX -------------------

              ! Discrete N_flux = J/e
              ! Consider only upward FACs
              if (iono_south_jr(i,j) > 0.) then
                 discrete_nf(i,j) = ABS(iono_south_jr(i,j)) / cElectronCharge
              else
                 discrete_nf(i,j) = numflux_floor
              end if

              ! Set a minimum on diffuse number flux
              if (diffuse_nf(i,j) < numflux_floor) then
                 diffuse_nf(i,j) = numflux_floor
              end if

              ! Set a minimum on electron temperature = 1 eV.
              if (iono_T(i,j) < cEVtoK) then
                 iono_T(i,j) = cEVtoK
              end if              

              ! Deal with the following portion carefully
              ! Lots of avenues to get non-physical answers
              if (iono_south_jr(i,j) > 0.) then
                 ! Variable must not be greater than 
                 vari = ABS(iono_south_jr(i,j)) / &
                      (1.6E-19* diffuse_nf(i,j))
                 if ((rm - vari) < 0.) then
                    vari = rm - 5E-05 ! Try
                 end if
                 
                 ! Potential Energy eV
                 eV(i,j) = 1.38E-23 * iono_T(i,j) * (rm - 1) * &
                      LOG((rm - 1)/(rm - vari))
              else
                 eV(i,j) = 0.
              end if

              ! Grand Calculation for the Energy Flux (in W/m2)
              ! See Knight [1973] and Fridman and Lemaire [1980] for details.
              !write(*,*)'!!! South', i,j,ev(i,j),iono_T(i,j),rm

              var_rm = EXP(-eV(i,j)/(1.38e-23 * iono_T(i,j) * &  !
                   (rm - 1)))                                    !
              numerator   = 1 - var_rm                           ! DUMMY VARS
              denominator = 1 + (1 - 1/rm) * var_rm              !
              if (denominator == 0.) then
                 denominator = 1E-15
              end if
              var_tot = eV(i,j) * numerator/denominator          ! 
              
              discrete_ef(i,j) = (var_tot + 2*1.38e-23*iono_T(i,j)) * &
                   discrete_nf(i,j)

              !write(*,*) eV(i,j), vari
              
              ! Average Energy can be estimated as,
              ! Ave_E = E_Flux/N_Flux (in eV)
              discrete_ae(i,j) = discrete_ef(i,j)/&!discrete_nf(i,j)
                   (1.6e-19 * discrete_nf(i,j) )
              
              ! Convert EFlux and Ave_E into units for Robinson relation
              discrete_ef(i,j) = discrete_ef(i,j) !* 10.0 !* 1E03 ! Convert to ergs/cm2
              discrete_ae(i,j) = discrete_ae(i,j) / 1E03 ! Convert to keV

              ! Limit 1 - Equatorward Boundary
              if ((180. - IONO_SOUTH_Theta(i,j)*cRadToDeg) >= 40.) then
                 !if (iono_south_jr(i,j) <= 1E-09) then
                 !write(*,*) IONO_SOUTH_Theta(i,j)*cRadToDeg, 180. - IONO_SOUTH_Theta(i,j)*cRadToDeg
                 !write(*,*) 'Discrete and Diffuse AE and EF (before)', &
                 !     discrete_ae(i,j), discrete_ef(i,j), diffuse_ae(i,j), diffuse_ef(i,j)
                 discrete_ae(i,j) = IONO_Min_Ave_E
                 discrete_ef(i,j) = IONO_Min_EFlux
                 diffuse_ae(i,j) = IONO_Min_Ave_E
                 diffuse_ef(i,j) = IONO_Min_EFlux
                 !write(*,*) 'Discrete and Diffuse AE and EF (before)', &
                 !     discrete_ae(i,j), discrete_ef(i,j), diffuse_ae(i,j), diffuse_ef(i,j)
              end if

              ! Limit 2 - Open-Closed FL Boundary
              if ((180. - IONO_SOUTH_Theta(i,j)*cRadToDeg) <= 180. - OCFLB(j)*cRadToDeg) then
                 !write(*,*) 'CoLat', 180. - IONO_SOUTH_Theta(i,j)*cRadToDeg, 'OCFLB', 180. - OCFLB(j)*cRadToDeg 
                 !write(*,*) 'Discrete and Diffuse AE and EF (before)', &
                 !     discrete_ae(i,j), discrete_ef(i,j), diffuse_ae(i,j), diffuse_ef(i,j)
                 discrete_ae(i,j) = IONO_Min_Ave_E
                 discrete_ef(i,j) = IONO_Min_EFlux
                 diffuse_ae(i,j) = IONO_Min_Ave_E
                 diffuse_ef(i,j) = IONO_Min_EFlux
                 !write(*,*) 'Discrete and Diffuse AE and EF (after)', &
                 !     discrete_ae(i,j), discrete_ef(i,j), diffuse_ae(i,j), diffuse_ef(i,j)
              end if
              
              !iono_north_jr(i,j) = iono_north_jr(i,j) / 10.
           enddo
        enddo

        where (iono_south_ave_e > 100.0) &
          iono_south_ave_e = IONO_Min_Ave_E!100.0

        if (DoTest) then
           !write(*,*) 'S_Discr: MaxVal Downward FAC', minval(iono_south_jr)
           write(*,*) 'S_Discr: MaxVal NFlux ', maxval(discrete_nf)
           write(*,*) 'S_Discr: MaxVal EFlux (in W/m2)', maxval(discrete_ef)
           write(*,*) 'S_Discr: MaxVal Ave_E (in keV)', maxval(discrete_ae)
        end if
        
        ! Let's add a little conductance on ANY not in the polar cap,
        ! so the code doesn't blow up
        !do j = 1, IONO_nPsi
        !   do i = 1, IONO_nTheta
        !      if (abs(iono_south_jr(i,j)) > 0.0 .and. &
        !           discrete_ef(i,j) < 5.0e-3) then
        !         discrete_ef(i,j) = MulFac_Def/2.5 * &
        !              (iono_south_jr(i,j)*1e6)**2*discrete_k(i,j)
        !      endif
        !   enddo
        !enddo

        where (diffuse_ae < IONO_Min_Ave_E/2) diffuse_ae = IONO_Min_Ave_E/2
        where (discrete_ae < IONO_Min_Ave_E/2) discrete_ae = IONO_Min_Ave_E/2

        ! -----------------------------------------
        ! Smoothing Filter for Discrete Energy Flux
        ! -----------------------------------------
        !
        ! During out simulations, it appears that a combination of
        ! boundary issues coupled with the discrete flux calculations
        ! leads to localized pockets of high FACs with unphysical
        ! values as the simulation marches in time. This step is
        ! added to make sure such localization is arrested early on.

        ! Smoothing Filter on Boundaries - discrete_ef
        call smooth_lagrange_polar(discrete_ef, IONO_nTheta, IONO_nPsi)

        ! Smoothing Filter on Boundaries - discrete_nf
        call smooth_lagrange_polar(discrete_nf, IONO_nTheta, IONO_nPsi)

        ! -----------------------------------------------------
        ! Calculate net electron fluxes =======================
        ! -----------------------------------------------------

        ! NOTE: This calculation is only for diffuse and discrete
        ! electron precip. For ion diffuse and broadband, the flux
        ! is converted into conductance differently. - AM

        if(DoTest) then ! Debug
           write(*,*)'IONO_nPsi and nTheta = ', IONO_nPsi, IONO_nTheta
           write(*,*)'Let us look at diffuse_ef, discrete_ef, diffuse_nf, discrete_nf'
           do j = 1, IONO_nPsi
              do i = 1, IONO_nTheta
                 write(*,'(2(i4.4,1x),4(E10.5, 1x))') i,j,diffuse_ef(i,j), &
                      discrete_ef(i,j), diffuse_nf(i,j), discrete_nf(i,j)
              end do
           end do
        end if

        ! Energy Flux ================================================
        ! Weight the Average Energy by the Number Flux, which is ef/av
        
        ! These will be parameterized in the future... - AM & DTW
        !iono_south_ave_e = &            ! If only using diffuse precip
        !     diffuse_ef /  & 
        !     (1E03 * 1.6e-19 * diffuse_nf)
        !iono_south_ave_e = discrete_ae  ! If only using discrete precip
        
        iono_south_ave_e =                & ! If using both precip
             (diffuse_ef + discrete_ef) / &
             (1E03 * 1.6e-19 *            & ! Convert this into keV (x 10^3)
             (diffuse_nf + discrete_nf))

        ! Average Energy =========================================
        ! The energy flux should be weighted by the average energy

        ! These will be parameterized in the future... - AM & DTW
        !iono_south_eflux = diffuse_ef      ! If only using diffuse prcip
        !iono_south_eflux = discrete_ef     ! If only using discrete precip
        iono_south_eflux =                 & ! If using both precip
             (1E03 * 1.6e-19 *             &
             (diffuse_nf + discrete_nf)) * &
             iono_south_ave_e

        where (iono_south_ave_e > 10.) &
          iono_south_ave_e = IONO_Min_Ave_E
        where (iono_south_ave_e < IONO_Min_Ave_E) &
             iono_south_ave_e = IONO_Min_Ave_E
        !where (iono_south_ave_e < IONO_Min_Ave_E) &
        !     iono_south_ave_e = IONO_Min_Ave_E
        !where (iono_south_eflux < IONO_Min_EFlux) &
        !     iono_south_eflux = IONO_Min_EFlux
        !where (iono_south_eflux < IONO_Min_EFlux) &
        !     iono_south_ave_e = IONO_Min_Ave_E
        !where (iono_south_eflux < IONO_Min_EFlux) &
        !     iono_south_ave_e = IONO_Min_Ave_E
        
        ! The Great Loop of our time...
        !do j = 1, IONO_nPsi
        !   do i = 1, IONO_nTheta
        !   end do
        !end do
        
        !      ! Equatorward Boundary Limit
        !      if ((IONO_SOUTH_Theta(i,j)*cRadToDeg) >= 40.) then
        !         !if (iono_south_jr(i,j) <= 1E-09) then 
        !         iono_south_ave_e(i,j) = IONO_Min_Ave_E
        !         iono_south_eflux(i,j) = IONO_Min_EFlux
        !      end if

              !if ((IONO_SOUTH_Theta(i,j)*cRadToDeg) >= 45.) then
              !   !if (iono_north_jr(i,j) <= 1E-09) then 
              !      hal_a1 = hal_a0
              !      ped_a1 = ped_a0
              !   end if

              ! Open-Closed Boundary Limit
              !if ((IONO_SOUTH_Theta(i,j)*cRadToDeg) <= OCFLB(j)*cRadToDeg) then
              !   !if (iono_north_jr(i,j) <= 1E-09) then 
              !   iono_south_ave_e(i,j) = IONO_Min_Ave_E
              !   iono_south_eflux(i,j) = IONO_Min_EFlux
              !end if

        !   end do
        !end do

        ! Equatorward Boundary Limit
        where ((180. - IONO_SOUTH_Theta*cRadToDeg) >= 40.) &
             iono_south_ave_e = IONO_Min_Ave_E
        where ((180. - IONO_SOUTH_Theta*cRadToDeg) >= 40.) &
             iono_south_eflux = IONO_Min_EFlux

        ! For Storage
        IONO_SOUTH_DIFF_Ave_E = diffuse_ae
        IONO_SOUTH_DIFF_EFlux = diffuse_ef

        IONO_SOUTH_MONO_Ave_E = discrete_ae
        IONO_SOUTH_MONO_EFlux = discrete_ef

        
        ! Smoothing Filter on Boundaries - Total Electron E-Flux
        call smooth_lagrange_polar(iono_south_eflux, IONO_nTheta, IONO_nPsi)

        ! Smoothing Filter on Boundaries - Total Electron Ave-E
        call smooth_lagrange_polar(iono_south_ave_e, IONO_nTheta, IONO_nPsi)
        
        if(DoTest) write(*,*) "MAGNIT: Done with southern aurora."

!        do j = 1, IONO_nPsi
!           do i = 1, IONO_nTheta
!              if (iono_south_p(i,j) <= 0.0) then
!                 iono_south_ave_e(i,j) = max(iono_south_ave_e(i,j),polarcap_avee)
!                 iono_south_eflux(i,j) = &
!                      max(iono_south_eflux(i,j),polarcap_eflux)
!              endif
!           enddo
!        enddo

     endif

  endif
  
end subroutine FACs_to_fluxes
!==============================================================================
subroutine ionosphere_conductance(Sigma0, SigmaH, SigmaP, &
     SigmaThTh, SigmaThPs, SigmaPsPs,      &
     dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
     dSigmaPsPs_dTheta,                    &
     dSigmaThTh_dPsi, dSigmaThPs_dPsi,     &
     dSigmaPsPs_dPsi,                      &
     Eflux, Ave_E,                         &
     Theta, Psi, nTheta, nPsi,             &
     dTheta, dPsi,                         &
     iModel, f107)

  ! This subroutine computes the height-integrated field-aligned and
  ! Hall and Pedersen conductances for the ionosphere at each
  ! location of the discretized solution domain.  The gradients of
  ! these quantities are also computed.

  use IE_ModMain
  use ModIonosphere
  use ModProcIE
  use CON_world
  use CON_comp_param
  implicit none

  integer :: nTheta,nPsi
  integer, intent(in) :: iModel
  real, intent(in)    :: f107
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
       Sigma0, SigmaH, SigmaP, &
       SigmaThTh, SigmaThPs, SigmaPsPs, &
       dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
       dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
       Eflux, Ave_E,                                            &
       Theta, Psi, sin_clat, cos_clat, sin_lon, cos_lon, cy,    &
       oSigmaH, oSigmaP, conv_SigmaH, conv_SigmaP, cos_SZA,     &
       tmp_x, tmp_y, tmp_z, tmp_sigp, tmp_sigh

  real, dimension(1:IONO_NTheta) :: dTheta
  real, dimension(1:IONO_NPsi)   :: dPsi

  integer :: i,j
  real :: f107p53, f107p49
  real :: sn, cs, sn2, cs2, cs3, cs4, C
  real :: SigmaH_EUV, SigmaP_EUV, SigmaH_SCAT, SigmaP_SCAT, &
       SigmaH_STAR, SigmaP_STAR, SigmaH_EUV_2, SigmaP_EUV_2
  real :: SigmaH_Particles, SigmaP_Particles, tau
  logical :: old

  real    :: time_delay
  real    :: cos_limit, meeting_value_h, meeting_value_p

  !--------------------------------------------------------------------------

  if (theta(1,1) < 2.0*IONO_Theta_0) then
     north = .true.
     tmp_x = IONO_NORTH_X
     tmp_y = IONO_NORTH_Y
     tmp_z = IONO_NORTH_Z
  else
     north = .false.
     tmp_x = IONO_SOUTH_X
     tmp_y = IONO_SOUTH_Y
     tmp_z = IONO_SOUTH_Z
  endif

  if (use_comp(UA_)) then

     Sigma0 = 1000.00

     ! Use Conductances from the UAM, which are put into the following
     ! variables

     if (north) then
        SigmaH = IONO_NORTH_SigmaH
        SigmaP = IONO_NORTH_SigmaP
     else
        SigmaH = IONO_South_SigmaH
        SigmaP = IONO_South_SigmaP
     endif

  else

     if (iModel > 1) then

        cos_SZA = (tmp_x*cosTHETATilt - tmp_z*sinTHETATilt) &
             / sqrt(tmp_x**2 + tmp_y**2 + tmp_z**2)

        ! We are going to need F10.7 ^ 0.53 and F10.7 ^ 0.49 a lot,
        ! So, let's just store them straight away:

        f107p53 = f107**0.53
        f107p49 = f107**0.49
        cos_limit = cos(70.0*cDegToRad)
        meeting_value_p = f107p49*(0.34*cos_limit+0.93*sqrt(cos_limit))
        meeting_value_h = f107p53*(0.81*cos_limit+0.54*sqrt(cos_limit))

     endif

     if (iModel.eq.0) then

        do j = 1, nPsi
           do i = 1, nTheta
              
              Sigma0(i,j) = 1000.00
              SigmaH(i,j) = 0.00
              SigmaP(i,j) = StarLightPedConductance

           enddo
        enddo

     endif

     if (iModel.eq.10) then

        do j = 1, nPsi
           do i = 1, nTheta
              if (north) then
                 SigmaP(i,j) = PedConductance_North
              else
                 SigmaP(i,j) = PedConductance_South
              endif

              Sigma0(i,j) = 1000.00
              SigmaH(i,j) = 0.00

           enddo
        enddo

     endif

     if (iModel.eq.1) then

        do j = 1, nPsi
           do i = 1, nTheta

              Sigma0(i,j) = 1000.00
              SigmaH(i,j) = PolarCapPedConductance
              SigmaP(i,j) = StarLightPedConductance

           enddo
        enddo

     endif

     if (iModel.eq.2) then

        do j = 1, nPsi
           do i = 1, nTheta

              ! Rasmussen and Schunk model B (JGR, Vol. 92, pp. 4491-4504, 1987).

              Sigma0(i,j) = 1000.00

              if (cos_SZA(i,j) > 0) then
                 SigmaH_EUV=f107p53*(0.81*cos_SZA(i,j)+0.54*sqrt(cos_SZA(i,j)))
                 SigmaP_EUV=f107p49*(0.34*cos_SZA(i,j)+0.93*sqrt(cos_SZA(i,j)))
                 SigmaH_SCAT = 1.00
                 SigmaP_SCAT = 0.50
                 if (cos_SZA(i,j) < cos_limit) then
                    SigmaH_EUV_2 = meeting_value_h *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaP_EUV_2 = meeting_value_p *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaH_EUV = (SigmaH_EUV + SigmaH_EUV_2)/2.0
                    SigmaP_EUV = (SigmaP_EUV + SigmaP_EUV_2)/2.0
                 endif
              else
                 SigmaH_EUV = meeting_value_h *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaP_EUV = meeting_value_p *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaH_SCAT = 1.00*(10.00**cos_SZA(i,j))
                 SigmaP_SCAT = 0.50*(10.00**cos_SZA(i,j))
              end if

              SigmaH_STAR = StarLightPedConductance*2.0
              SigmaP_STAR = StarLightPedConductance

              SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV +                     &
                   SigmaH_SCAT*SigmaH_SCAT +                                 &
                   SigmaH_STAR*SigmaH_STAR)
              SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV +                     &
                   SigmaP_SCAT*SigmaP_SCAT +                                 &
                   SigmaP_STAR*SigmaP_STAR)

           enddo
        enddo

     endif

     if ((iModel.ge.3) .and. (iModel.ne.10)) then

        do j = 1, nPsi
           do i = 1, nTheta

              Sigma0(i,j) = 1000.00

              if (cos_SZA(i,j) > 0) then
                 SigmaH_EUV=f107p53*(0.81*cos_SZA(i,j)+0.54*sqrt(cos_SZA(i,j)))
                 SigmaP_EUV=f107p49*(0.34*cos_SZA(i,j)+0.93*sqrt(cos_SZA(i,j)))
                 SigmaH_SCAT = 1.00
                 SigmaP_SCAT = 0.50
                 if (cos_SZA(i,j) < cos_limit) then
                    SigmaH_EUV_2 = meeting_value_h *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaP_EUV_2 = meeting_value_p *   &
                         exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                    SigmaH_EUV = (SigmaH_EUV + SigmaH_EUV_2)/2.0
                    SigmaP_EUV = (SigmaP_EUV + SigmaP_EUV_2)/2.0
                 endif
              else
                 SigmaH_EUV = meeting_value_h *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaP_EUV = meeting_value_p *   &
                      exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                 SigmaH_SCAT = 1.00*(10.00**cos_SZA(i,j))
                 SigmaP_SCAT = 0.50*(10.00**cos_SZA(i,j))
              end if

              SigmaH_STAR = StarLightPedConductance*2.0
              SigmaP_STAR = StarLightPedConductance

              ! Use Robinson's Formula to convert the Ave_E and E_Flux to 
              ! SigmaP and SigmaH

              SigmaP_Particles = 40.0 * Ave_E(i,j) /                     &
                   (16.0 + Ave_E(i,j)*Ave_E(i,j))  *                     &
                   sqrt(EFlux(i,j)*1000.0)

              SigmaH_Particles = 0.45 * (Ave_E(i,j)**0.85) * SigmaP_Particles

              SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + &
                   SigmaH_SCAT*SigmaH_SCAT + &
                   SigmaH_STAR*SigmaH_STAR + &
                   SigmaH_Particles*SigmaH_Particles)

              SigmaP_EUV = SigmaP_EUV*SigmaP_EUV + SigmaP_SCAT*SigmaP_SCAT + &
                   SigmaP_STAR*SigmaP_STAR

              SigmaP_Particles = SigmaP_Particles*SigmaP_Particles

              SigmaP(i,j) = sqrt(SigmaP_EUV + SigmaP_Particles)

           enddo

        enddo

     endif

     !if (iModel.eq.9) then

     !   do j = 1, nPsi
     !      do i = 2, nTheta-1
     !         tmp_sigp(i,j)=(sigmap(i-1,j)+sigmap(i,j)+sigmap(i+1,j))/3
     !         tmp_sigh(i,j)=(sigmah(i-1,j)+sigmah(i,j)+sigmah(i+1,j))/3
     !      enddo
     !   enddo

     !   sigmap=tmp_sigp
     !   sigmah=tmp_sigh

     !endif
        
     if (north) then

        ! Subsolar is i=nTheta ; j = 0
        ! Just in case we change this in the future, add noon and midnight to
        ! test for conductance

        if (SAVE_NORTH_SigmaH(nTheta,1) +                                     &
             SAVE_NORTH_SigmaH(nTheta,nPsi/2) > 0.0) then

           ! We have "old" conductances

           old = .true.
           SigmaH = SAVE_NORTH_SigmaH
           SigmaP = SAVE_NORTH_SigmaP

        endif

     else

        ! Subsolar is i=0 ; j = 0
        ! Just in case we change this in the future, add noon and midnight to
        ! test for conductance

        if (SAVE_SOUTH_SigmaH(1,1) +                                         &
             SAVE_SOUTH_SigmaH(1,nPsi/2) > 0.0) then

           ! We have "old" conductances

           old = .true.
           SigmaH = SAVE_SOUTH_SigmaH
           SigmaP = SAVE_SOUTH_SigmaP

        endif

     endif

     if (north) then
        SAVE_NORTH_SigmaP = 0.0
        SAVE_NORTH_SigmaH = 0.0
     else
        SAVE_SOUTH_SigmaP = 0.0
        SAVE_SOUTH_SigmaH = 0.0
     endif

  end if

  do j = 1, nPsi
     do i = 1, nTheta

        sn = sin(Theta(i,j))
        cs = cos(Theta(i,j))
        sn2= sn*sn
        cs2 = cs*cs
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)
        C = 4.00*Sigma0(i,j)*cs2 + &
             SigmaP(i,j)*sn2

        SigmaThTh(i,j) = Sigma0(i,j)*SigmaP(i,j)*cs3/C
        SigmaThPs(i,j) = 2.00*Sigma0(i,j)*SigmaH(i,j)* &
             cs*cs4/C
        SigmaPsPs(i,j) = SigmaP(i,j)+ &
             SigmaH(i,j)*SigmaH(i,j)* &
             sn2/C

        dSigmaThTh_dTheta(i,j) = 0.00
        dSigmaThTh_dPsi(i,j) = 0.00
        dSigmaThPs_dTheta(i,j) = 0.00
        dSigmaThPs_dPsi(i,j) = 0.00
        dSigmaPsPs_dTheta(i,j) = 0.00
        dSigmaPsPs_dPsi(i,j) = 0.00

     end do
  end do

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                (dTheta(i))
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/ &
                (dPsi(j))

           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                (dTheta(i))
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/ &
                (dPsi(j))

           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                (dTheta(i))
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/ &
                (dPsi(j))
        end do
     else if (j == 1) then
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                (dTheta(i))
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nPsi-1))/ &
                (dPsi(j))

           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                (dTheta(i))
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nPsi-1))/ &
                (dPsi(j))

           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                (dTheta(i))
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nPsi-1))/ &
                (dPsi(j))
        end do
     else
        do i = 2, nTheta-1
           dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/ &
                (dTheta(i))
           dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,2)-SigmaThTh(i,j-1))/ &
                (dPsi(j))

           dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/ &
                (dTheta(i))
           dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,2)-SigmaThPs(i,j-1))/ &
                (dPsi(j))

           dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/ &
                (dTheta(i))
           dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,2)-SigmaPsPs(i,j-1))/ &
                (dPsi(j))
        end do
     end if
  end do

end subroutine ionosphere_conductance
!==============================================================================
subroutine determine_oval_characteristics(Current_in, Theta_in, Psi_in, &
     Loc_of_Oval, Width_of_Oval, Strength_of_Oval)

  ! This routine calculates everything in radians away from the pole.

  use ModIonosphere
  use IE_ModIo,       ONLY: NameIonoDir
  use ModConductance, ONLY: DoOvalShift
  use IE_ModMain,     ONLY: Time_Array, nSolve
  use ModIoUnit,      ONLY: UnitTMP_
  use ModUtilities,   ONLY: open_file, close_file

  implicit none

  ! inputs:
  real, dimension(1:IONO_nTheta, 1:IONO_nPsi), intent(in) :: &
       Current_in, Theta_in, Psi_in

  ! Outputs:
  real, dimension(1:IONO_nPsi), intent(out) :: &
       Loc_of_Oval, Width_of_Oval, Strength_of_Oval

  real, dimension(1:IONO_nTheta, 1:IONO_nPsi) :: &
       Current, Theta, Psi

  real, dimension(1:8) :: &
       max_fac, max_fac_colat, width, J_Save

  real    :: day_colat, dusk_colat, midnight_colat, dawn_colat
  real    :: day_fac, dusk_fac, midnight_fac, dawn_fac
  real    :: noon_mid, dusk_dawn, day_strength, night_strength
  real    :: mean_colat, dev_colat, sumFAC, Night_Width, Day_Width

  integer :: i, j, n, nloc, dJ, J_Start, J_End

  ! Testing and output variables:
  character(len=100), save    :: NameFile, StringFormat
  logical:: IsFirstWrite = .true.
  logical:: DoTest, DoTestMe

  character(len=*), parameter :: NameSub = 'determine_oval_characteristics'
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  
  ! Reverse the Arrays for Southern Hemisphere:
  if (Theta_in(1,1) < cHalfPi) then
     Current = Current_in
     Theta   = Theta_in
     Psi     = Psi_in
     north   = .true.
  else
     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi
           Current(IONO_nTheta - (i-1), j) = Current_in(i,j)
           Theta(IONO_nTheta - (i-1), j)   = cPi - Theta_in(i,j)
           Psi(IONO_nTheta - (i-1), j)     = Psi_in(i,j)
        enddo
     enddo
     north   = .false.
  endif

  !
  ! Start the Oval Determination
  !

  dJ = IONO_nPsi/9

  do n = 1, 8 

     !
     ! figure out location of auroral oval:
     ! Let's start a little ways away from the pole...
     !

     nloc = 1
     J_Start = max(1,(n-1)*dJ)
     J_End   = min(n*dJ, IONO_nPsi)
     max_fac(n) = abs(Current(3,J_Start))
     max_fac_colat(n) = Theta(3,J_Start)
     J_Save(n) = J_Start 
     do j = J_Start, J_End
        do i = 4, IONO_nTheta
           if (Current(i,j) > max_fac(n)) then
              max_fac(n) = abs(Current(i,j))
              max_fac_colat(n) = Theta(i,j)
              J_Save(n) = j
              nloc = i
           endif
        enddo
     enddo

     !
     ! figure out width
     !

     width(n) = 0.0
     j = J_Save(n)
     do i = nloc, IONO_nTheta
        if (Current(i,j) > max_fac(n)/4.0) then
           width(n) = abs(max_fac_colat(n) - Theta(i,j))
        endif
     enddo

     if (width(n).le.(theta(2,1)-theta(1,1))) width(n) = max_fac_colat(n)/5.0

  enddo

  day_colat = (max_fac_colat(1) + max_fac_colat(8))/2.0
  dusk_colat = (max_fac_colat(2) + max_fac_colat(3))/2.0
  midnight_colat = (max_fac_colat(4) + max_fac_colat(5))/2.0
  dawn_colat = (max_fac_colat(6) + max_fac_colat(7))/2.0

  day_fac = (max_fac(1) + max_fac(8))/2.0
  dusk_fac = (max_fac(2) + max_fac(3))/2.0
  midnight_fac = (max_fac(4) + max_fac(5))/2.0
  dawn_fac = (max_fac(6) + max_fac(7))/2.0

  night_width = 0.0
  sumFAC = 0.0
  mean_colat = 0.0

  do n=1,8
     night_width = night_width + width(n) * max_fac(n)
     mean_colat = mean_colat + max_fac_colat(n) * max_fac(n)
     sumFAC = sumFAC + max_fac(n)
  enddo

  mean_colat = mean_colat/sumFAC
  Night_Width = Night_Width/sumFAC

  if (Night_Width > 6.0*cDegToRad) Night_Width=6.0*cDegToRad
  Day_Width = max(Night_Width/2.0,1.0*cDegToRad)

  if (mean_colat < 15.0*cDegToRad) then
     mean_colat = 15.0*cDegToRad
     dev_colat = 0.0
  else if(.not.DoOvalShift) then
     dev_colat = 0.0
  else
     dev_colat = ((day_colat - mean_colat) * day_fac - &
          (midnight_colat - mean_colat) * midnight_fac) / &
          (day_fac + midnight_fac)

     ! Restrict auroral location a little bit:

     if (abs(dev_colat) > mean_colat/2.0) then
        dev_colat = dev_colat*mean_colat/2.0/abs(dev_colat)
     endif
  endif

  Day_Strength = dawn_fac
  Night_Strength = sumFAC


  !  dev_colat = 1.0*cDegToRad

  !  Day_Width = 2.5*cDegToRad
  !  Night_Width = Day_Width

  !  mean_colat = 2.0*cDegToRad

  do j=1,IONO_nPsi
     Loc_of_Oval(j)   = mean_colat - dev_colat*cos(Psi(1,j))
     Width_of_Oval(j) = Day_Width + &
          (Night_Width - Day_Width)*sin(Psi(1,j)/2.0)
     Strength_of_Oval(j) = Day_Strength + &
          (Night_Strength - Day_Strength)*sin(Psi(1,j)/2.0)
  enddo

  ! For testing purposes, write oval info to file.
  if(.not.DoTestMe .or. .not.north) return
  
  if(IsFirstWrite)then
     ! Open file:
     write(NameFile,'(a,i8.8,a)')trim(NameIonoDir)//'aurora_n',nSolve,'.txt'
     call open_file(FILE=NameFile, STATUS='replace', NameCaller=NameSub)
     ! Write header w/ longitudes:
     write(UnitTmp_, '(a)', advance='NO')'Auroral oval location at Lon='
     do j=1, IONO_nPsi
        write(UnitTmp_,'(1x,f6.2)', advance='NO') Psi(1,j)*cRadToDeg
     end do
     write(UnitTmp_,'(a)')'','YYYY MM DD HH MN SS msc oval_colat'

     ! Set format codes for writing output:
     write(StringFormat,'("(i4,5i3,i4,", i4,"(1x,f6.2))")') IONO_nPsi
     IsFirstWrite=.false.
  else
     ! Open file in append mode:
     call open_file(FILE=NameFile, STATUS='old', POSITION='append', &
          NameCaller=NameSub)
  end if

  ! Write record:
  write(UnitTmp_, StringFormat) Time_Array(1:7), Loc_of_oval(:)*cRadToDeg

  call close_file(NameCaller=NameSub)

end subroutine determine_oval_characteristics

subroutine ADELPHI_Calc(jr,psi,hall,ped)

  use ModIonosphere
  use ModNumConst, ONLY: cPi, cDegToRad

  implicit none

  ! inputs:
  real :: jr, psi

  ! outputs:
  real :: hall, ped

  ! working variables:
  real :: A_h0,B_h0,C_h0,A_h1,B_h1,C_h1
  real :: A_p0,B_p0,C_p0,A_p1,B_p1,C_p1
  real :: phi
  real :: hal_a0,hal_a1,ped_a0,ped_a1

  if (jr > 0) then
     
     A_h0 = 8.70
     B_h0 = 4.60
     C_h0 = 327.10
     A_h1 = 14.80
     B_h1 = -10.40
     C_h1 = 129.90
     A_p0 = 4.20
     B_p0 = 1.10
     C_p0 = 318.60
     A_p1 = 6.80
     B_p1 = -1.50
     C_p1 = 184.90
     
  else
     
     A_h0 = 7.70
     B_h0 = -1.80
     C_h0 = 139.00
     A_h1 = -7.30
     B_h1 = 5.60
     C_h1 = 100.90
     A_p0 = 5.00
     B_p0 = -0.80
     C_p0 = 60.90
     A_p1 = -3.20
     B_p1 = -3.60
     C_p1 = 21.90
     
  endif
    
  phi = mod((psi+cPi),2.0*cPi)
  hal_a0 = A_h0 + B_h0*cos(C_h0*cDegToRad+phi)
  hal_a1 = A_h1 + B_h1*cos(C_h1*cDegToRad+phi)
  ped_a0 = A_p0 + B_p0*cos(C_p0*cDegToRad+phi)
  ped_a1 = A_p1 + B_p1*cos(C_p1*cDegToRad+phi)
  
  ! I try to make sure all coefficients are positive. By Zihan Wang. 09/10/2021
  ! I also add a mask function. By Zihan Wang. 10/08/2021
  if (abs(jr*1.0e6)>0.01) then
     hall=hal_a0+hal_a1*jr*1.0e6
     ped=ped_a0+ped_a1*jr*1.0e6
  else                 
     hall=3
     ped=1.5
  endif
end subroutine ADELPHI_Calc
