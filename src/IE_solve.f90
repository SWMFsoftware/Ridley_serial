!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine IE_solve

  use IE_ModMain
  use IE_ModIO, ONLY: write_prefix
  use ModProcIE
  use ModIonosphere
  use ModNumConst
  use ModMpi
  use ModConductance, ONLY: generate_conductance

  implicit none

  real    :: CurrentSum
  integer :: iBlock
  integer :: nSize, iError

  logical DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'IE_solve'
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  SinThetaTilt = sin(ThetaTilt)
  CosThetaTilt = cos(ThetaTilt)

  if (DoTest .and. iProc==0) &
       write(*,'(a,i4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3)')&
       " "//NameSub//" at ",Time_Array

  do iBlock = 1, 2

     if(DoTest)write(*,*) 'iblock',iblock

     select case(iBlock)
     case(1) ! Northern hemisphere
        if(iProc /= 0) CYCLE

        north = .true.

        CurrentSum = sum(abs(IONO_NORTH_JR))
        if(DoTest)write(*,*)NameSub,': sum(abs(IONO_NORTH_JR))=', CurrentSum
        if(CurrentSum < 1e-6)CYCLE
        if (UseFakeRegion2) then
           call Create_Region2_Currents(1)
           IONO_NORTH_JR = IONO_NORTH_JR + IONO_NORTH_Fake_JR
           if(DoTest)write(*,*)NameSub,': after UseFakeRegion2=',&
                sum(abs(IONO_NORTH_JR))
        end if

        ! Add the IM currents before the conductances are calculated
        ! CANDIDATE FOR REMOVAL: good GM resolution and 2-way GM-IM
        ! coupling makes this unnecessary.
        IONO_NORTH_JR = IONO_NORTH_JR + FractionImJr*iono_north_im_jr

        ! Next two calls replaced with:
        call generate_conductance('north')

        ! Add in ionospheric currents after the conductances
        ! are calculated
        IONO_NORTH_JR = IONO_NORTH_JR - Iono_North_Tgcm_Jr

        call ionosphere_solver(iBlock, &
             IONO_NORTH_JR,     &
             IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs,   &
             IONO_NORTH_SigmaPsPs,                         &
             IONO_NORTH_dSigmaThTh_dTheta,                 &
             IONO_NORTH_dSigmaThPs_dTheta, &
             IONO_NORTH_dSigmaPsPs_dTheta, &
             IONO_NORTH_dSigmaThTh_dPsi,  &
             IONO_NORTH_dSigmaThPs_dPsi, &
             IONO_NORTH_dSigmaPsPs_dPsi, &
             IONO_NORTH_Theta, IONO_NORTH_Psi, &
             dTheta_North, dPsi_North, &
             IONO_NORTH_PHI)

        if(DoTest)then
           call write_prefix;
           write(*,'(a,f8.3,a)') "Northern Cross Polar Cap Potential=", &
                cpcp_north," kV"
        end if

        ! Calculate Currents and Boundary Conditions for North

        call ionosphere_currents(iBlock, &
             IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz,&
             IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
             IONO_NORTH_ETh,IONO_NORTH_EPs, &
             IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
             IONO_NORTH_PHI, &
             IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
             IONO_NORTH_SigmaPsPs, &
             IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
             IONO_NORTH_Theta, IONO_NORTH_Psi, &
             dTheta_North, dPsi_North)

        ! add joule heating for north (JouleHeating = sigmaP * E^2)
        ! Yiqun
        call ionosphere_jouleheating_ionflux(iBlock, &
             IONO_NORTH_ETh, IONO_NORTH_EPs, &
             IONO_NORTH_SigmaP, &
             IONO_NORTH_Joule,  &
             IONO_NORTH_IonNumFlux)

     case(2) ! Southern hemisphere

        if(iProc /= nProc-1) CYCLE

        north = .false.

        CurrentSum = sum(abs(IONO_SOUTH_JR))
        if(DoTest)write(*,*)NameSub,': sum(abs(IONO_SOUTH_JR))=', CurrentSum
        if(CurrentSum < 1e-6)CYCLE

        if (UseFakeRegion2) then
           call Create_Region2_Currents(iBlock)
           IONO_SOUTH_JR = IONO_SOUTH_JR + IONO_SOUTH_Fake_JR
           if(DoTest)write(*,*)NameSub,': after UseFakeRegion2=',&
                sum(abs(IONO_SOUTH_JR))
        end if

        ! Add the IM currents before the conductances are calculated
        IONO_SOUTH_JR = IONO_SOUTH_JR + FractionImJr*iono_south_im_jr

        ! Obtain conductance:
        call generate_conductance('south')

        ! Add in ionospheric currents after the conductances
        ! are calculated
        IONO_SOUTH_JR = IONO_SOUTH_JR - Iono_South_Tgcm_Jr

        call ionosphere_solver(iBlock, &
             IONO_SOUTH_JR, &
             IONO_SOUTH_SigmaThTh, &
             IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_dSigmaThTh_dTheta, &
             IONO_SOUTH_dSigmaThPs_dTheta, &
             IONO_SOUTH_dSigmaPsPs_dTheta, &
             IONO_SOUTH_dSigmaThTh_dPsi, &
             IONO_SOUTH_dSigmaThPs_dPsi, &
             IONO_SOUTH_dSigmaPsPs_dPsi, &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             dTheta_South, dPsi_South, &
             IONO_SOUTH_PHI)

        if(DoTest)then
           call write_prefix;
           write(*,'(a,f8.3,a)') "Southern Cross Polar Cap Potential=",&
                cpcp_south," kV"
        end if

        ! Calculate Currents and Boundary Conditions for South

        call ionosphere_currents(iBlock, &
             IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz,&
             IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
             IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
             IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
             IONO_SOUTH_PHI, &
             IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             dTheta_South, dPsi_South)

        ! add joule heating for north (JouleHeating = sigmaP * E^2)
        ! Yiqun
        call ionosphere_jouleheating_ionflux(iBlock, &
             IONO_SOUTH_ETh, IONO_SOUTH_EPs, &
             IONO_SOUTH_SigmaP, &
             IONO_SOUTH_Joule,  &
             IONO_SOUTH_IonNumFlux)

     end select

  end do

  ! Broadcast the sourthern current and the potential for IM if needed
  if(TypeImCouple /= 'north' .and. nProc > 1)then

     nSize = IONO_nTheta * IONO_nPsi
     call MPI_bcast(IONO_SOUTH_PHI, nSize, MPI_REAL, 1, iComm, iError)
     call MPI_bcast(IONO_SOUTH_JR,  nSize, MPI_REAL, 1, iComm, iError)
     call MPI_bcast(cpcp_south,         1, MPI_REAL, 1, iComm, iError)

  end if

end subroutine IE_solve
!==============================================================================
