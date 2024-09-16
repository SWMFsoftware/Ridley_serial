!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IE_wrapper

  ! Wrapper for Ridley's ionosphere model RIM

  use ModUtilities, ONLY: CON_set_do_test, CON_stop, open_file, close_file
  use ModIoUnit, ONLY: UNITTMP_

  implicit none

  private ! except

  ! CON wrapper
  public:: IE_set_param
  public:: IE_init_session
  public:: IE_run
  public:: IE_save_restart
  public:: IE_finalize

  ! Coupling with GM
  public:: IE_get_for_gm ! Also used by IM/RAM
  public:: IE_put_from_gm

  ! Coupling with IM
  public:: IE_get_for_im
  public:: IE_put_from_im
  public:: IE_put_from_im_complete

  ! Coupling with PS
  public:: IE_get_for_ps

  ! Coupling with PW
  public:: IE_get_for_pw

  ! Coupling with RB
  public:: IE_get_for_rb

  ! Coupling with UA
  public:: IE_get_info_for_ua
  public:: IE_get_for_ua
  public:: IE_put_from_ua

  ! Lookup table for F10.7 flux
  integer:: iTableF107 = -100

contains
  !============================================================================
  subroutine IE_set_param(CompInfo, TypeAction)

    use ModProcIE
    use ModIonosphere
    use IE_ModIo
    use IE_ModMain
!    use ModIonoMagPerturb

    use ModIoUnit
    use CON_comp_info

    ! Arguments
    type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp
    character(len=*), intent(in)      :: TypeAction ! What to do

    character(len=*), parameter :: NameSub = 'IE_set_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use=.true.,                                    &
            NameVersion='Serial Potential Solver (Ridley)',&
            Version=1.1)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

       if( nProc>2 )call CON_stop(NameSub//' IE_ERROR '//&
            'this version can run on 1 or 2 PE-s only!')
       if( NameIonoDir(1:3) /= 'IE/' ) NameIonoDir = 'IE/'//NameIonoDir
    case('READ','CHECK')
       call read_param
    case('STDOUT')
       iUnitOut=STDOUT_
       if(nProc==1)then
          StringPrefix='IE:'
       else
          write(StringPrefix,'(a,i1,a)')'IE',iProc,':'
       end if
    case('FILEOUT')
       call get(CompInfo,iUnitOut=iUnitOut)
       StringPrefix=''
    case('GRID')
       call IE_set_grid
    case default
       call CON_stop(NameSub//' IE_ERROR: invalid TypeAction='//TypeAction)
    end select

  contains
    !==========================================================================

    subroutine read_param

      use ModReadParam
      use ModIE_Interface
      use ModFiles
      use ModConductance, ONLY: DoUseEuvCond, f107_flux, &
           PolarCapPedCond, StarLightCond, SigmaHalConst, SigmaPedConst, &
           imodel_legacy, DoUseAurora, NameAuroraMod
      use ModIeRlm, ONLY: UseOval, UseNewOval, DoOvalShift, &
           UseSubOvalCond, DoFitCircle, FactorHallCMEE, FactorPedCMEE, &
           NameHalFile, NamePedFile, LatNoConductanceSI, UseCMEEFitting
      use ModUtilities,   ONLY: fix_dir_name, check_dir, lower_case

      ! The name of the command
      character (len=100) :: NameCommand

      ! Read parameters
      logical :: UseStrict=.true., IsUninitialized=.true.

      ! Interface variables for legacy #IONOSPHERE command
      logical :: UseFullCurrent ! LEGACY, NOT USED, CANDIDATE FOR REMOVAL.
      integer :: iModLegacy=5
      real :: f10Legacy, StarLightLegacy, PolarCapLegacy

      ! Plot file parameters
      integer :: iFile, iDebugProc
      character (len=50) :: plot_string

      !------------------------------------------------------------------------
      select case(TypeAction)
      case('CHECK')
         IsUninitialized=.false.

         ! We should check and correct parameters here
         if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()

         RETURN
      case('READ')
         if(iProc==0)write(*,*) NameSub,': READ iSession =',i_session_read(),&
              ' iLine=',i_line_read(),' nLine =',n_line_read()

         IsUninitialized=.false.
      end select

      ! Read input data from text via ModReadParam
      do
         if(.not.read_line() ) EXIT
         if(.not.read_command(NameCommand)) CYCLE

         select case(NameCommand)
         case("#STRICT")
            call read_var('UseStrict',UseStrict)

         ! I/O-related params
         case("#IONODIR")
            call read_var("NameIonoDir",NameIonoDir)
            call fix_dir_name(NameIonoDIr)
            if (iProc==0) call check_dir(NameIonoDir)
         case("#SAVEPLOT", "#IE_SAVEPLOT")
            call read_var('nPlotFile',nFile)
            if (nFile > MaxFile)call CON_stop(NameSub//&
                 ' number of ouput files is too large in #IE_SAVEPLOT:'&
                 //' nFile>MaxFile')
            if (nFile>0.and.iProc==0) call check_dir(NameIonoDir)
            do iFile=1,nFile

               call read_var('StringPlot',plot_string)
               call lower_case(plot_string)

               ! Check to see if the ionosphere directory exists...
               if(iProc==0)call check_dir(NameIonoDir)

               ! Plotting frequency
               call read_var('DnSavePlot',dn_output(iFile))
               call read_var('DtSavePlot',dt_output(iFile))

               ! Plot file format
               if(index(plot_string,'idl')>0)then
                  plot_form(iFile)='idl'
               elseif(index(plot_string,'tec')>0)then
                  plot_form(iFile)='tec'
               else
                  call CON_stop(NameSub//&
                       ' IE_ERROR format (idl,tec) missing from plot_string='&
                       //plot_string)
               end if
               if(index(plot_string,'min')>0)then
                  plot_vars(iFile)='minimum'
               elseif(index(plot_string,'max')>0)then
                  plot_vars(iFile)='maximum'
               elseif(index(plot_string,'aur')>0)then
                  plot_vars(iFile)='aur'
               elseif(index(plot_string,'uam')>0)then
                  plot_vars(iFile)='uam'
               elseif(index(plot_string,'xyz')>0)then
                  plot_vars(iFile)='xyz'
               else
                  call CON_stop(NameSub//&
                       ' variable definition missing in #IE_SAVEPLOT'//&
                       ' from plot_string='//plot_string)
               end if
            end do
         case("#SAVEPLOTNAME")
            call read_var('IsPlotName_e',IsPlotName_e)
         case("#SAVELOGNAME")
            call read_var('IsLogName_e',IsLogName_e)
         case("#SAVELOGFILE")
            call read_var('DoSaveLogfile',DoSaveLogfile)
            if(DoSaveLogfile)then
               if(iProc==0)call check_dir(NameIonoDir)
            endif

         ! Conductance related params:
         case('#SOLAREUV')
            call read_var('DoUseEuvCond', DoUseEuvCond)
            if(DoUseEuvCond) call read_var('F10.7 Flux', f107_flux)
         case('#UNIFORMCONDUCTANCE')
            call read_var('SigmaPedConst', SigmaPedConst)
            call read_var('SigmaHalConst', SigmaHalConst)
         case('#BACKGROUNDCOND')
            call read_var('StarLightCond', StarLightCond)
            call read_var('PolarCapPedCond', PolarCapPedCond)
         case('#AURORA')
            call read_var('DoUseAurora', DoUseAurora)
            if(DoUseAurora) call read_var('NameAuroraMod', NameAuroraMod)
         case("#CONDUCTANCEFILES")
            call read_var('NameFileHall',     NameHalFile)
            call read_var('NameFilePedersen', NamePedFile)
         case('#USECMEEFIT')
            call read_var('UseCMEEFitting', UseCMEEFitting)
            if (UseCMEEFitting) then
               call read_var('LatNoConductanceSI', LatNoConductanceSI)
               call read_var('FactorHallCMEE',     FactorHallCMEE)
               call read_var('FactorPedCMEE',      FactorPedCMEE)
            endif
         case("#AURORALOVAL")
            call read_var('UseOval', UseOval)
            if(UseOval)then
               call read_var('UseOvalShift',          DoOvalShift)
               call read_var('UseSubOvalConductance', UseSubOvalCond)
               call read_var('UseAdvancedOval',       UseNewOval)
               if(UseNewOval) call read_var('DoFitCircle', DoFitCircle)
            end if
         case("#RLMCONDUCTANCE")
            call read_var('LatNoConductanceSI', LatNoConductanceSI)
            call read_var('OvalWidthFactor',    OvalWidthFactor)
            call read_var('OvalStrengthFactor', OvalStrengthFactor)
            call read_var('ConductanceFactor',  CondFactor)
            if (trim(NameAuroraMod).eq.'RLM3') then
               write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                    ' command '//trim(NameCommand)// &
                    ' can only be used with conductance model 4 or 5'
               if(UseStrict)call CON_stop('Correct PARAM.in!')
            end if

         ! Physics & solver related params
         case("#IM")
            call read_var('TypeImCouple',TypeImCouple)
            call lower_case(TypeImCouple)
            call read_var('FractionImJr',FractionImJr)
         case("#BOUNDARY")
            call read_var('LatBoundary',LatBoundary)
            LatBoundary = LatBoundary * cDegToRad
         case("#UA")
            call read_var('DoCoupleUaCurrent',DoCoupleUaCurrent)
            if(DoCoupleUaCurrent)then
               call read_var('LatBoundary',LatBoundary)
               LatBoundary = LatBoundary * cDegToRad
            endif
         case("#SPS")
            call read_var('UseSPS',UseSPS)
            IE_NameOfEFieldModel = "SPS"
            UseGridBasedIE = .true.
         case("#SOLVER")
            call read_var('NameSolver',        NameSolver, IsLowerCase=.true.)
         case("#KRYLOV")
            call read_var('UsePreconditioner', UsePreconditioner)
            call read_var('UseInitialGuess',   UseInitialGuess)
            call read_var('Tolerance',         Tolerance)
            call read_var('MaxIteration',      MaxIteration)
         case("#DEBUG")
            call read_var('iDebugLevel',iDebugLevel)
            call read_var('iDebugProc',iDebugProc)
            if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
               iDebugLevel = -1
            endif
         case("#AMIEFILES")
            call read_var('NameAmieFileNorth',AMIEFileNorth)
            call read_var('NameAmieFileSouth',AMIEFileSouth)
            IE_NameOfEFieldModel = "amie"
            UseGridBasedIE = .true.
            UseAMIE = .true.
         case("#BACKGROUND")

            call read_var('NameOfModelDir',IE_NameOfModelDir)
            call read_var('NameOfEFieldModel',IE_NameOfEFieldModel)
            call read_var('NameOfAuroralModel',IE_NameOfAuroralModel)
            call read_var('NameOfSolarModel',IE_NameOfSolarModel)

            if (index(IE_NameOfAuroralModel,'IHP') > 0) &
                 IE_NameOfAuroralModel = 'ihp'
            if (index(IE_NameOfAuroralModel,'PEM') > 0) &
                 IE_NameOfAuroralModel = 'pem'

            if (index(IE_NameOfEFieldModel,'AMIE') > 0) &
                 IE_NameOfEFieldModel = 'amie'

            if (index(IE_NameOfEFieldModel,'weimer01') > 0) &
                 IE_NameOfEFieldModel = 'weimer01'
            if (index(IE_NameOfEFieldModel,'Weimer01') > 0) &
                 IE_NameOfEFieldModel = 'weimer01'
            if (index(IE_NameOfEFieldModel,'WEIMER01') > 0) &
                 IE_NameOfEFieldModel = 'weimer01'

            if (index(IE_NameOfEFieldModel,'weimer') > 0 .and. &
                 index(IE_NameOfEFieldModel,'01') == 0) &
                 IE_NameOfEFieldModel = 'weimer96'
            if (index(IE_NameOfEFieldModel,'Weimer') > 0 .and. &
                 index(IE_NameOfEFieldModel,'01') == 0) &
                 IE_NameOfEFieldModel = 'weimer96'
            if (index(IE_NameOfEFieldModel,'WEIMER') > 0 .and. &
                 index(IE_NameOfEFieldModel,'01') == 0) &
                 IE_NameOfEFieldModel = 'weimer96'

            if (index(IE_NameOfEFieldModel,'weimer96') > 0) &
                 IE_NameOfEFieldModel = 'weimer96'
            if (index(IE_NameOfEFieldModel,'Weimer96') > 0) &
                 IE_NameOfEFieldModel = 'weimer96'
            if (index(IE_NameOfEFieldModel,'WEIMER96') > 0) &
                 IE_NameOfEFieldModel = 'weimer96'

            if (index(IE_NameOfEFieldModel,'SAMIE') > 0) &
                 IE_NameOfEFieldModel = 'samie'

            UseGridBasedIE = .false.

         case("#RESTART")
            call read_var('DoRestart', DoRestart)

         ! The following are LEGACY PARAMS and CANDIDATES FOR REMOVAL.
         case("#IONOSPHERE")
            ! Read LEGACY variables and store locally.
            call read_var('iConductanceModel', iModLegacy)
            call read_var('UseFullCurrent' ,   UseFullCurrent)
            call read_var('UseFakeRegion2' ,   UseFakeRegion2)
            call read_var('F10.7 Flux',        f10Legacy)
            call read_var('StarLightPedConductance', StarLightLegacy)
            call read_var('PolarCapPedConductance',  PolarCapLegacy)
            ! Set options as a function of iModel:
            call imodel_legacy(iModLegacy, f10Legacy, &
                 StarLightLegacy, PolarCapLegacy)
         case("#MAGNETOMETER")
            write(*,*)'IE_WARNING: #MAGNETOMETER COMMAND NOW GM-ONLY.'
         case("#GEOMAGINDICES")
            write(*,*)'IE_WARNING: #GEOMAGINDICES COMMAND NOW GM-ONLY.'

         case default
            if(iProc==0) then
               write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                    ' invalid command '//trim(NameCommand)
               if(UseStrict)call CON_stop('Correct PARAM.in!')
            end if
         end select
      end do

    end subroutine read_param

  end subroutine IE_set_param
  !============================================================================
  subroutine IE_set_grid

    ! Set the grid descriptor for IE
    ! Since IE has a static grid the descriptor has to be set once.
    ! There can be many couplers that attempt to set the descriptor,
    ! so we must check IsInitialized.

    use ModProcIE
    use ModIonosphere
    use IE_ModIo
    use IE_ModMain
    use CON_coupler
    use ModNumConst

    logical :: IsInitialized=.false.
    integer :: iProc_A(2)

    real :: Colat_I(2*IONO_nTheta-1) = 0.
    real :: Longitude_I(IONO_nPsi) = 0.

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'IE_set_grid'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized

    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! IE runs on 1 or 2 PE-s so processor array is (/0,0/) or (/0,1/)
    iProc_A(1)=0
    iProc_A(2)=nProc-1

    ! Coordinates can only be set on the IE processors where nProc > 0
    if (nProc > 0)then
       ! The colatitudes for both hemispheres
       Colat_I(            1:  IONO_nTheta) = IONO_NORTH_Theta(:,1)
       Colat_I(IONO_nTheta:2*IONO_nTheta-1) = IONO_SOUTH_Theta(:,1)
       Longitude_I                          = IONO_NORTH_Psi(1,:)
    end if

    call set_grid_descriptor(                        &
         IE_,                                        &! component index
         nDim=2,                                     &! dimensionality
         nRootBlock_D=[2,1],                         &! north+south hemispheres
         nCell_D =[IONO_nTheta - 1,IONO_nPsi - 1],   &! size of node based grid
         XyzMin_D=[1.0, 1.0],                        &! minimum indexes
         XyzMax_D=[real(2*IONO_nTheta-1.0),          &! maximum indexes
         real(IONO_nPsi)],                           &
         TypeCoord='SMG',                            &! solar magnetic coord.
         Coord1_I=Colat_I,                           &! colatitudes
         Coord2_I=Longitude_I,                       &! longitudes
         Coord3_I=[IONO_Radius + IONO_Height],       &! radial size in meters
         iProc_A = iProc_A)                           ! processor assigment

  end subroutine IE_set_grid
  !============================================================================

  subroutine IE_get_for_gm(Buffer_IIV, iSize, jSize, nVar, NameVar_I, &
       tSimulation)

    ! Put variables listed in NameVar_I into the buffer

    use ModProcIE
    use ModIonosphere

    integer,          intent(in) :: iSize, jSize, nVar
    real,             intent(out):: Buffer_IIV(iSize,jSize,nVar)
    character(len=*), intent(in) :: NameVar_I(nVar)
    real,             intent(in) :: tSimulation

    real:: tSimulationTmp
    integer:: iVar

    character(len=*), parameter:: NameSub = 'IE_get_for_gm'
    !--------------------------------------------------------------------------
    if(iSize /= IONO_nTheta*2-1 .or. jSize /= IONO_nPsi)then
       write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
            ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    ! Make sure that the most recent result is provided
    tSimulationTmp = tSimulation
    call IE_run(tSimulationTmp, tSimulation)

    do iVar = 1, nVar
       select case(NameVar_I(iVar))
       case('potential')
          Buffer_IIV(:,:,iVar) = IONO_Phi
       case('jouleheat')
          Buffer_IIV(:,:,iVar) = IONO_Joule
       case('sigmahall')
          Buffer_IIV(:,:,iVar) = IONO_SigmaH
       case('sigmapedersen')
          Buffer_IIV(:,:,iVar) = IONO_SigmaP
       case default
          call CON_stop(NameSub//': unknown NameVar='//NameVar_I(iVar))
       end select
    end do

  end subroutine IE_get_for_gm
  !============================================================================

  subroutine IE_get_for_pw(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
       tSimulation)

    use ModProcIE
    use ModIonosphere

    integer, intent(in)           :: iSize, jSize, nVar
    real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
    character (len=*),intent(in)  :: NameHem
    character (len=*),intent(in)  :: Name_V(nVar)
    real,             intent(in)  :: tSimulation

    integer :: iVar
    real    :: tSimulationTmp
    character(len=*), parameter :: NameSub='IE_get_for_pw'
    !--------------------------------------------------------------------------
    if(iSize /= 2*IONO_nTheta-1 .or. jSize /= IONO_nPsi)then
       write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
            ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    ! Make sure that the most recent result is provided
    tSimulationTmp = tSimulation
    call IE_run(tSimulationTmp,tSimulation)

    if(iProc /= 0) RETURN

    do iVar = 1, nVar
       select case(Name_V(iVar))
       case('Pot')
          Buffer_IIV(:,:,iVar) = IONO_Phi
       case('Jr')
          Buffer_IIV(:,:,iVar) = IONO_Jr
       case('Ave')
          Buffer_IIV(:,:,iVar) = IONO_Ave_E
       case('Tot')
          Buffer_IIV(:,:,iVar) = IONO_EFlux*1.0e3 ! J/m^2 --> ergs/cm^2
       case default
          call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
       end select
    end do

  end subroutine IE_get_for_pw
  !============================================================================

  subroutine IE_get_for_rb(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
       tSimulation)

    use ModProcIE
    use ModIonosphere

    integer, intent(in)           :: iSize, jSize, nVar
    real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
    character (len=*),intent(in)  :: NameHem
    character (len=*),intent(in)  :: Name_V(nVar)
    real,             intent(in)  :: tSimulation

    integer :: iVar
    real    :: tSimulationTmp
    character(len=*), parameter :: NameSub='IE_get_for_rb'
    !--------------------------------------------------------------------------
    if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
       write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
            ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    ! Make sure that the most recent result is provided
    tSimulationTmp = tSimulation
    call IE_run(tSimulationTmp,tSimulation)

    select case(NameHem)

    case('North')

       if(iProc /= 0) RETURN
       do iVar = 1, nVar
          select case(Name_V(iVar))
          case('Pot')
             Buffer_IIV(:,:,iVar) = IONO_NORTH_Phi
          case default
             call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
          end select
       end do

    case('South')

       if(iProc /= nProc - 1) RETURN
       do iVar = 1, nVar
          select case(Name_V(iVar))
          case('Pot')
             Buffer_IIV(:,:,iVar) = IONO_SOUTH_Phi
          case default
             call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
          end select
       end do

    case default

       call CON_stop(NameSub//' invalid NameHem='//NameHem)

    end select

  end subroutine IE_get_for_rb
  !============================================================================
  subroutine IE_get_info_for_ua(nVar, NameVar_V)
    ! Get number and names of variables for UA to IE coupling.
    ! IE reports what variables it needs from UA.
    ! UA will use this info to create and fill buffers appropriately.

    integer, intent(out) :: nVar
    character(len=*), intent(out), optional :: NameVar_V(:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IE_get_info_for_ua'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Right now, only request Hall & Pedersen conductance.
    ! In the future, neutral wind FAC coupling or other values can/will
    ! be shared from UA.
    nVar = 4
    if(present(NameVar_V)) NameVar_V(1:4) = ['lon','lat','hal','ped']

    if(DoTestMe)then
       write(*,*) NameSub//': nVar=', nVar
       if(present(NameVar_V)) write(*,*) NameSub//': NameVar_V=',NameVar_V
    end if

  end subroutine IE_get_info_for_ua
  !============================================================================

  subroutine IE_get_for_ua(Buffer_IIV, iSize, jSize, nVarIn, NameVar_V, &
       iBlock,tSimulation)

    use ModProcIE
    use ModIonosphere

    integer,          intent(in)  :: iSize, jSize, nVarIn, iBlock
    real,             intent(out) :: Buffer_IIV(iSize,jSize,nVarIn)
    character (len=*),intent(in)  :: NameVar_V(nVarIn)
    real,             intent(in)  :: tSimulation

    character(len=5) :: NameHem
    integer :: iVar
    real    :: tSimulationTmp

    ! Debug variables:
    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'IE_get_for_ua'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)then
       write(*,*)NameSub//': Gathering data for UA at t=', tSimulation
       write(*,*)NameSub//': Current hemisphere block = ', iBlock
       write(*,*)NameSub//': Variable list = '
       do iVar=1, nVarIn
          write(*,'(10x,i2.2,5x,a)')iVar, NameVar_V(iVar)
       end do
       write(*,*)NameSub//': Expected IE grid dimensions (iSize,jSize) =', &
            iSize, jSize
    end if

    if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
       write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
            ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    ! Set hemisphere to transfer:
    if(iBlock==1) then
       NameHem = 'North'
    else if (iBlock==2) then
       NameHem = 'South'
    else
       write(*,*) NameSub//': invalid iBlock = ', iBlock
       call CON_stop(NameSub//': Error coupling with UA')
    end if

    ! Make sure that the most recent result is provided
    tSimulationTmp = tSimulation
    call IE_run(tSimulationTmp,tSimulation)

    select case(NameHem)

    case('North')

       if(iProc /= 0) RETURN

       do iVar=1, nVarIn
          select case(NameVar_V(iVar))

          case('pot')
             Buffer_IIV(:,:,iVar) = IONO_NORTH_Phi
          case('ave')
             Buffer_IIV(:,:,iVar) = IONO_NORTH_Ave_E
          case('tot')
             Buffer_IIV(:,:,iVar) = IONO_NORTH_EFlux
          case default
             call CON_stop(NameSub//' invalid NameVar='//NameVar_V(iVar))
          end select
       end do

    case('South')

       if(iProc /= nProc - 1) RETURN

       do iVar=1, nVarIn
          select case(NameVar_V(iVar))

          case('pot')
             Buffer_IIV(:,:,iVar) = IONO_SOUTH_Phi
          case('ave')
             Buffer_IIV(:,:,iVar) = IONO_SOUTH_Ave_E
          case('tot')
             Buffer_IIV(:,:,iVar) = IONO_SOUTH_EFlux
          case default
             call CON_stop(NameSub//' invalid NameVar='//NameVar_V(iVar))
          end select
       end do
    end select

  end subroutine IE_get_for_ua
  !============================================================================

  subroutine IE_put_from_gm(Buffer_IIV, iSize, jSize, nVar)

    use IE_ModMain, ONLY: IsNewInput
    use ModProcIE
    use ModIonosphere
    use ModConductance, ONLY: GmRhoFloor, GmPFloor

    integer,          intent(in) :: iSize, jSize, nVar
    integer                      :: i, j
    real                         :: Buffer_IIV(iSize, jSize, nVar)

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'IE_put_from_gm'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,' starting with iSize, jSize, nVar=', &
         iSize, jSize, nVar

    IsNewInput = .true.

    ! Here, we change the coupler values to be palatable to the
    ! SPECIFIC code: coord transform, unit conversion, floors, ceilings.

    ! Set minimum acceptable values for density & pressure:
    where (Buffer_IIV(:,:,3) < GmRhoFloor) Buffer_IIV(:,:,3)=GmRhoFloor
    where (Buffer_IIV(:,:,4) < GmPFloor  ) Buffer_IIV(:,:,4)=GmPFloor

    if (iProc == 0) then
       Iono_North_Jr = Buffer_IIV(1:IONO_nTheta,:,1)
       Iono_North_Jr(IONO_nTheta-1:IONO_nTheta,1) = 0.0
       if(nVar>1)then
          Iono_North_invB = Buffer_IIV(1:IONO_nTheta,:,2)
          Iono_North_rho  = Buffer_IIV(1:IONO_nTheta,:,3)
          Iono_North_p    = Buffer_IIV(1:IONO_nTheta,:,4)
          Iono_North_dLat = Buffer_IIV(1:IONO_nTheta,:,5)
          Iono_North_dLon = Buffer_IIV(1:IONO_nTheta,:,6)
          if(DoTest) call write_dataN
       end if
    endif
    if (iProc == nProc-1) then
       Iono_South_Jr = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,1)
       Iono_South_Jr(1:2,1) = 0.0
       if(nVar>1)then
          Iono_South_invB = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,2)
          Iono_South_rho  = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,3)
          Iono_South_p    = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,4)
          Iono_South_dLat = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,5)
          Iono_South_dLon = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,6)
          if(DoTest) call write_dataS
       end if
    endif

    ! Debug statements:
    if(DoTest)then
       ! write statements here.
       do j = 1, IONO_nPsi
          do i = 1, IONO_nTheta
             write(*,'(2(i4.4,1x), 2(E10.5, 1x))') i, j, &
                  Iono_North_p(i,j), Iono_North_rho(i,j)
          end do
       end do

    endif

    if(DoTest)write(*,*)NameSub,' finished'

  contains
    !==========================================================================

    subroutine write_dataN

      ! write values to North plot file

      character(len=80) :: filename
      integer:: i, j
      integer:: nCall=0
      character(len=*), parameter:: NameSub = 'write_dataN'
      !------------------------------------------------------------------------
      nCall=nCall+1
      write(filename,'(a,i5.5,a)')"gm2ie_debugN_",nCall,".dat"
      call open_file(FILE=filename, STATUS='unknown')
      write(UNITTMP_,'(a)') 'TITLE="gm2ie debugN values"'
      write(UNITTMP_,'(a)') &
           'VARIABLES="J", "I", "Theta", "Psi", "JR", "1/B", "rho", "p"'
      write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
           'ZONE T="PE=',iProc,'", I=',IONO_nPsi,&
           ', J=',IONO_nTheta,', K=1, F=POINT'
      do i=1,IONO_nTheta; do j=1, IONO_nPsi
         write(UNITTMP_,'(2i4,6G14.6)') j,i, &
              Iono_North_Theta(i,j),Iono_North_Psi(i,j),Iono_North_Jr(i,j), &
              Iono_North_invB(i,j),Iono_North_rho(i,j),Iono_North_p(i,j)
      end do; end do
      call close_file(NameCaller=NameSub)

    end subroutine write_dataN
    !==========================================================================
    subroutine write_dataS

      ! write values to South plot file

      character(len=80) :: filename
      integer:: i,j
      integer:: nCall=0
      character(len=*), parameter:: NameSub = 'write_dataS'
      !------------------------------------------------------------------------
      nCall=nCall+1
      write(filename,'(a,i5.5,a)')"gm2ie_debugS_",nCall,".dat"
      call open_file(FILE=filename, STATUS='unknown')
      write(UNITTMP_,'(a)') 'TITLE="gm2ie debugS values"'
      write(UNITTMP_,'(a)') &
           'VARIABLES="J", "I", "Theta", "Psi", "JR", "1/B", "rho", "p"'
      write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
           'ZONE T="PE=',iProc,'", I=',IONO_nTheta, &
           ', J=',IONO_nTheta,', K=1, F=POINT'
      do i=1,IONO_nTheta; do j=1,IONO_nPsi
         write(UNITTMP_,'(2i4,6G14.6)') j,i, &
              Iono_South_Theta(i,j),Iono_South_Psi(i,j),Iono_South_Jr(i,j), &
              Iono_South_invB(i,j),Iono_South_rho(i,j),Iono_South_p(i,j)
      end do; end do
      call close_file(NameCaller=NameSub)

    end subroutine write_dataS
    !==========================================================================

  end subroutine IE_put_from_gm
  !============================================================================
  subroutine IE_put_from_ua(Buffer_IIBV, nMLTs, nLats, nVarIn, NameVarUaIn_V)

    use IE_ModMain,     ONLY: IsNewInput, DoCoupleUaCurrent
    use ModConductance, ONLY: StarLightCond
    use ModIonosphere
    use ModConst
    use ModUtilities,   ONLY: check_allocate

    save

    ! Arguments: returning IE variables on a MLT-Lat grid, one
    ! per hemisphere, for nVarIn variables.
    integer,          intent(in) :: nMlts, nLats, nVarIn
    character(len=3), intent(in) :: NameVarUaIn_V(nVarIn)
    real,             intent(in) :: Buffer_IIBV(nMlts, nLats, 2, nVarIn)

    ! UA_Lats and UA_Mlts are the latitudes and magnetic local times of the
    !    UA magnetic grid.
    ! iLat and iMlt are the indices of where to find the points in the
    !    UA magnetic grid.
    ! rLat and rMlt are the multiplication factors to get stuff from the
    !    UA magnetic grid.

    real,    dimension(:,:,:), allocatable :: UA_Lats, UA_Mlts
    integer, dimension(Iono_nTheta,2) :: iLat
    integer, dimension(Iono_nPsi,2)   :: iMlt
    real,    dimension(Iono_nTheta,2) :: rLat
    real,    dimension(Iono_nPsi,2)   :: rMlt

    real, dimension(IONO_nTheta,IONO_nPsi) :: TmpVar_II

    integer :: iError, i, j, ii, jj, iVar, iBlock
    real    :: t, p

    integer, parameter :: Fac_ = 1
    integer, parameter :: Ped_ = 2
    integer, parameter :: Hal_ = 3
    integer, parameter :: Lat_ = 4
    integer, parameter :: Mlt_ = 5

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='IE_put_from_ua'
    !--------------------------------------------------------------------------
    ! Set test variables.
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! Indicate that new information has arrived -> IE will recalculate solution
    IsNewInput=.true.

    ! Set intelligent defaults as needed:
    IONO_NORTH_TGCM_JR = 0.0
    IONO_SOUTH_TGCM_JR = 0.0

    ! Check if arrays are allocated & allocate as necessary
    if (.not.allocated(UA_Lats)) then
       allocate(UA_Lats(nMLTs, nLats,2), UA_Mlts(nMLTs, nLats,2), stat=iError)
       call check_allocate(iError,NameSub//'UA_Lats,UA_Mlts')
    endif

    ! Copy over lat/lon into UA_* vars
    do iVar=1, nVarIn
       select case (NameVarUaIn_V(iVar))
       case('lat')
          UA_Lats(:,:,1) = Buffer_IIBV(:,:,1,iVar)
          UA_Lats(:,:,2) = Buffer_IIBV(:,:,2,iVar)
       case('lon')
          UA_Mlts(:,:,1) = Buffer_IIBV(:,:,1,iVar)
          UA_Mlts(:,:,2) = Buffer_IIBV(:,:,2,iVar)
       case default
          ! Other vars handled below.
       end select
    enddo

    ! Set up index mapping for both hemispheres to
    ! interpolate from UA to IE grid
    BLOCK: do iBlock=1, 2
       ! Mapping for lat/colat:
       COLAT: do i = 1, IONO_nTheta ! convert colat to lat:
          if (iBlock == 1) t = 90.0 - Iono_North_Theta(i,1)*cRadToDeg
          if (iBlock == 2) t = 90.0 - Iono_South_Theta(i,1)*cRadToDeg
          ! In this instance, t is theta -- Aaron Ridley, 1998
          if (t >= maxval(UA_Lats(1,:,iBlock))) then
             ii = nLats-1
             iLat(i,iBlock) = ii
             rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                  (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
          else if (t < minval(UA_Lats(1,:,iBlock))) then
             ii = 1
             iLat(i,iBlock) = ii
             rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                  (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
          else
             ii = 1
             do while (ii < nLats)
                if ((t >= UA_Lats(1,ii,iBlock) .and. &
                     t <  UA_Lats(1,ii+1,iBlock)) .or.  &
                     (t <= UA_Lats(1,ii,iBlock) .and. &
                     t >  UA_Lats(1,ii+1,iBlock))) then
                   iLat(i,iBlock) = ii
                   rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                        (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
                   ii = nLats
                endif
                ii = ii+1
             enddo
          endif
       enddo COLAT

       ! Mappint for MLT/longitude
       LON: do j = 1, IONO_nPsi
          if (iBlock == 1) p = mod(Iono_North_Psi(1,j)*12.0/cPi + 12.0,24.0)
          if (iBlock == 2) p = mod(Iono_South_Psi(1,j)*12.0/cPi + 12.0,24.0)
          ! added a compatible time shift in CON couple

          jj = 1

          do while (jj < nMlts)

             if (p >= UA_Mlts(jj,1,iBlock) .and. &
                  p <  UA_Mlts(jj+1,1,iBlock)) then
                iMlt(j,iBlock) = jj
                rMlt(j,iBlock) = 1.0 - (p - UA_Mlts(jj,1,iBlock)) / &
                     (UA_Mlts(jj+1,1,iBlock) - UA_Mlts(jj,1,iBlock))
                jj = nMlts
             else if (p >= UA_Mlts(nMlts,1,iBlock) .and. &
                  p <  UA_Mlts(1,1,iBlock)) then
                iMlt(j,iBlock) = nMlts
                rMlt(j,iBlock) = 1.0 - (p - UA_Mlts(nMlts,1,iBlock)) / &
                     (UA_Mlts(1,1,iBlock) - UA_Mlts(nMlts,1,iBlock))
                jj = nMlts
             else if (p <= minval(UA_Mlts(:,1,iBlock))) then
                jj = minloc(UA_Mlts(:,1,iBlock),DIM=1)
                if (jj == 1) then
                   iMlt(j,iBlock) = nMlts
                   rMlt(j,iBlock) = 1.0 - (p - (UA_Mlts(nMlts,1,iBlock)-24)) /&
                        (UA_Mlts(jj,1,iBlock) - (UA_Mlts(nMlts,1,iBlock)-24))
                else
                   iMlt(j,iBlock) = jj-1
                   rMlt(j,iBlock) = 1.0 - (p - (UA_Mlts(jj-1,1,iBlock) - 24))  /&
                        (UA_Mlts(jj,1,iBlock) - (UA_Mlts(jj-1,1,iBlock) - 24))
                end if
                jj=nMlts
             else if (p >= maxval(UA_Mlts(:,1,iBlock))) then
                jj = maxloc(UA_Mlts(:,1,iBlock),DIM=1)
                if (jj == nMlts) then
                   iMlt(j,iBlock) = jj
                   rMlt(j,iBlock) = 1.0 - (p - UA_Mlts(jj,1,iBlock)) / &
                        (UA_Mlts(1,1,iBlock) - (UA_Mlts(jj,1,iBlock)-24))
                else
                   iMlt(j,iBlock) = jj
                   rMlt(j,iBlock) = 1.0 - (p - UA_Mlts(jj,1,iBlock)) / &
                        (UA_Mlts(jj+1,1,iBlock) - (UA_Mlts(jj,1,iBlock)-24))
                end if
                jj=nMlts
             end if
             if (iMlt(j,iBlock) == 0) then
                iMlt(j,iBlock) = 1
             end if
             jj = jj+1
          end do
       end do LON

    end do BLOCK

    ! Loop over variables and interpolate from UA to IE grid.
    VAR: do iVar=1, nVarIn
       ! Skip lat/lon
       if ((NameVarUaIn_V(iVar) == 'lat').or.(NameVarUaIn_V(iVar) == 'lon')) &
            CYCLE VAR
       do iBlock=1,2
          TmpVar_II=0.0
          ! Interpolate value to IE grid...
          do i = 1, Iono_nTheta
             ! t is the interpolation coefficient for theta
             ii = iLat(i,iBlock)
             t  = rLat(i,iBlock)

             do j = 1, Iono_nPsi
                ! p is the interpolation coefficient for psi
                jj = iMlt(j,iBlock)
                p  = rMlt(j,iBlock)

                ! Interpolate variables
                TmpVar_II(i,j) =    &
                     (    t)*(    p)*Buffer_IIBV(jj  ,ii  ,iBlock,iVar) + &
                     (1.0-t)*(    p)*Buffer_IIBV(jj  ,ii+1,iBlock,iVar) + &
                     (    t)*(1.0-p)*Buffer_IIBV(jj+1,ii  ,iBlock,iVar) + &
                     (1.0-t)*(1.0-p)*Buffer_IIBV(jj+1,ii+1,iBlock,iVar)
             end do
          end do

          ! Copy result to correct IE variable, taking care
          ! to use the correct hemisphere.  Conducantances are
          ! given floor values to prevent problems.
          select case (NameVarUaIn_V(iVar))
          case('hal') ! Hall conductance
             if(iBlock==1)then
                IONO_NORTH_SigmaH = max(TmpVar_II, 2*StarLightCond)
             else
                IONO_SOUTH_SigmaH = max(TmpVar_II, 2*StarLightCond)
             end if
          case('ped') ! Pedersen conductance
             if(iBlock==1)then
                IONO_NORTH_SigmaP = max(TmpVar_II, StarLightCond)
             else
                IONO_SOUTH_SigmaP = max(TmpVar_II, StarLightCond)
             end if
          case('fac') ! Neutral wind FACs
             if(iBlock==1)then
                IONO_NORTH_TGCM_JR = TmpVar_II
             else
                IONO_SOUTH_TGCM_JR = TmpVar_II
             end if
          case default
             call CON_stop(NameSub//' Unrecognized coupling variable: ', &
                  NameVarUaIn_V(iVar))
          end select
       end do
    end do VAR

  end subroutine IE_put_from_ua
  !============================================================================

  subroutine IE_put_from_im(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

    use IE_ModMain, ONLY: IsNewInput
    use ModIonosphere
    use ModProcIE

    use CON_router,   ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd
    integer:: iLat,iLon
    character(len=*), parameter:: NameSub = 'IE_put_from_im'
    !--------------------------------------------------------------------------
    if(nPoint>1)then
       write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
            nPoint,iPointStart,Weight % Weight_I
       call CON_stop(NameSub//': should be called with 1 point')
    end if
    if(DoAdd)then
       write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
            nPoint,iPointStart,Weight % Weight_I
       write(*,*)NameSub,': WARNING DoAdd is true'
    end if

    iLat = Index % iCB_II(1,iPointStart)
    iLon = Index % iCB_II(2,iPointStart)

    if ( iLat >= 1 .and. iLat <= iono_nTheta .and. &
         iLon >= 1 .and. iLon <= iono_nPsi) then
       if (.not.IsFilledWithIm(iLat,iLon)) then
          iono_north_im_jr(iLat,iLon)    = buff_v(1)
          iono_north_im_eflux(iLat,iLon) = buff_v(2)
          iono_north_im_avee(iLat,iLon)  = buff_v(3)
          IsFilledWithIm(iLat,iLon) = .true.
       endif
    endif

    IsNewInput = .true.

  end subroutine IE_put_from_im
  !============================================================================

  subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

    ! Provide potential and current for IM
    ! The value should be interpolated from nPoints with
    ! indexes stored in Index and weights stored in Weight
    ! The variables should be put into Buff_V

    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi, &
         IONO_NORTH_PHI, IONO_NORTH_JR, IONO_SOUTH_PHI, IONO_SOUTH_JR, &
         IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
         IONO_SOUTH_SigmaH, IONO_SOUTH_SigmaP, &
         cpcp_north, cpcp_south
    use IE_ModMain,    ONLY: TypeImCouple

    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight

    integer :: iBlock, i, j, iSouth, iPoint
    real    :: w

    character(len=*), parameter :: NameSub='IE_get_for_im'
    !--------------------------------------------------------------------------
    Buff_V = 0.0

    do iPoint = iPointStart, iPointStart + nPoint - 1

       i      = Index % iCB_II(1,iPoint)
       j      = Index % iCB_II(2,iPoint)
       iBlock = Index % iCB_II(3,iPoint)
       w      = Weight % Weight_I(iPoint)

       if(iBlock/=1)then
          write(*,*)NameSub,': iPoint,Index % iCB_II=',&
               iPoint,Index%iCB_II(:,iPoint)
          call CON_stop(NameSub//&
               ' SWMF_ERROR iBlock should be 1=North in IE-IM coupling')
       end if

       if(i<1 .or. i>IONO_nTheta .or. j<1 .or. j>IONO_nPsi)then
          write(*,*)'i,j=',i,j
          call CON_stop(NameSub//' SWMF_ERROR index out of range')
       end if

       ! Index for the same latitude on the southern hemisphere
       iSouth = IONO_nTheta + 1 - i

       select case(TypeImCouple)
       case('north')
          Buff_V(1) = Buff_V(1) + w * IONO_NORTH_PHI(i,j)
          Buff_V(2) = Buff_V(2) + w * IONO_NORTH_JR(i,j)
          Buff_V(3) = Buff_V(3) + w * IONO_NORTH_SigmaH(i,j)
          Buff_V(4) = Buff_V(4) + w * IONO_NORTH_SigmaP(i,j)
       case('south')
          Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_PHI(iSouth,j)
          Buff_V(2) = Buff_V(2) + w * IONO_SOUTH_JR(iSouth,j)
          Buff_V(3) = Buff_V(3) + w * IONO_SOUTH_SigmaH(iSouth,j)
          Buff_V(4) = Buff_V(4) + w * IONO_SOUTH_SigmaP(iSouth,j)
       case('cpcpmin')
          if(cpcp_north < cpcp_south)then
             Buff_V(1) = Buff_V(1) + w * IONO_NORTH_PHI(i,j)
             Buff_V(2) = Buff_V(2) + w * IONO_NORTH_JR(i,j)
             Buff_V(3) = Buff_V(3) + w * IONO_NORTH_SigmaH(i,j)
             Buff_V(4) = Buff_V(4) + w * IONO_NORTH_SigmaP(i,j)
          else
             Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_PHI(iSouth,j)
             Buff_V(2) = Buff_V(2) + w * IONO_SOUTH_JR(iSouth,j)
             Buff_V(3) = Buff_V(3) + w * IONO_SOUTH_SigmaH(iSouth,j)
             Buff_V(4) = Buff_V(4) + w * IONO_SOUTH_SigmaP(iSouth,j)
          end if
       case('average')
          Buff_V(1) = Buff_V(1) + w * &
               0.5*(IONO_NORTH_PHI(i,j) + IONO_SOUTH_PHI(iSouth,j))
          Buff_V(2) = Buff_V(2) + w * &
               0.5*(IONO_NORTH_JR(i,j)  + IONO_SOUTH_JR(iSouth,j))
          Buff_V(3) = Buff_V(3) + w * 0.5*( &
               IONO_NORTH_SigmaH(i,j)  + &
               IONO_SOUTH_SigmaH(iSouth,j))
          Buff_V(4) = Buff_V(4) + w * 0.5*( &
               IONO_NORTH_SigmaP(i,j)  + &
               IONO_SOUTH_SigmaP(iSouth,j))
       case default
          call CON_stop(NameSub//' ERROR: Unknown value for TypeImCouple='// &
               TypeImCouple)
       end select
    end do

  contains
    !==========================================================================

    real function minmod(a,b)
      real, intent(in) :: a,b
      !------------------------------------------------------------------------
      minmod = (sign(0.5, a) + sign(0.5, b)) * min(abs(a), abs(b))
    end function minmod
    !==========================================================================

  end subroutine IE_get_for_im
  !============================================================================

  subroutine IE_put_from_im_complete

    use ModProcIE
    use ModIonosphere
    use ModMpi

    integer iError, i
   !---------------------------------------------------------------------------
    iono_north_im_eflux(:,iono_npsi) = iono_north_im_eflux(:,1)
    iono_north_im_avee(:,iono_npsi)  = iono_north_im_avee(:,1)
    iono_north_im_jr(:,iono_npsi)  = iono_north_im_jr(:,1)

    if (nProc > 1) then
       iError = 0
       call MPI_Bcast(iono_north_im_eflux, iono_nTheta*iono_nPsi, &
            MPI_Real, 0, iComm, iError)
       call MPI_Bcast(iono_north_im_avee, iono_nTheta*iono_nPsi, &
            MPI_Real, 0, iComm, iError)
       call MPI_Bcast(iono_north_im_jr, iono_nTheta*iono_nPsi, &
            MPI_Real, 0, iComm, iError)
    endif

    do i = 1, IONO_nTheta
       iono_south_im_jr(i,:) = iono_north_im_jr(Iono_nTheta-i+1,:)
    enddo

    IsFilledWithIm = .false.

  end subroutine IE_put_from_im_complete
  !============================================================================

  subroutine IE_init_session(iSession, tSimulation)

    ! Initialize the Ionosphere Electrostatic (IE) module for session iSession

    use CON_physics,    ONLY: get_time, get_planet, get_axes
    use ModIonosphere,  ONLY: IONO_Bdp, init_mod_ionosphere
    use IE_ModMain,     ONLY: time_accurate, time_simulation, ThetaTilt
    use IE_ModIo,       ONLY: dt_output, t_output_last
    use ModConductance, ONLY: NameAuroraMod
    use ModIeRlm,       ONLY: load_conductances, UseCMEEFitting, &
         NameHalFile, NamePedFile
    use ModProcIE

    integer,  intent(in) :: iSession      ! session number (starting from 1)
    real,     intent(in) :: tSimulation   ! seconds from start time

    ! Initialize the Ionosphere Electrostatic (IE) module for session iSession

    logical :: IsUninitialized=.true.

    logical :: DoTest,DoTestMe
    character(len=*), parameter :: NameSub = 'IE_init_session'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(IsUninitialized)then
       ! Set configurations based on selected auroral models:
       if (trim(NameAuroraMod).eq.'CMEE') then
          ! Switch coefficient input files to CMEE:
          UseCMEEFitting = .true.
          NameHalFile = 'cmee_hal_coeffs.dat'
          NamePedFile = 'cmee_ped_coeffs.dat'
       end if

       call init_mod_ionosphere
       ! Read empirical conductance values from files as necessary
       if((index(NameAuroraMod,'RLM')>0).or.(index(NameAuroraMod,'CMEE')>0)) &
            call load_conductances()
       call ionosphere_fine_grid
       call ionosphere_init

       IsUninitialized = .false.
    end if

    time_simulation = tSimulation

    call get_time(  DoTimeAccurateOut = time_accurate)
    call get_planet(DipoleStrengthOut = IONO_Bdp)
    call get_axes(tSimulation, MagAxisTiltGsmOut = ThetaTilt)

    IONO_Bdp = IONO_Bdp*1.0e9 ! Tesla -> nT

    if(DoTest)write(*,*)NameSub,': IONO_Bdp, ThetaTilt =',IONO_Bdp,ThetaTilt

    ! Reset t_output_last in case the plotting frequency has changed
    if(time_accurate)then
       where(dt_output>0.) &
            t_output_last=int(time_simulation/dt_output)
    end if

  end subroutine IE_init_session
  !============================================================================

  subroutine IE_finalize(tSimulation)

    use ModProcIE
    use IE_ModMain, ONLY: Time_Array, time_simulation, nSolve, DoSaveLogFile
    use IE_ModIo, ONLY: nFile, unitlog
    use CON_physics, ONLY: get_time
    use ModTimeConvert, ONLY: time_real_to_int
    use ModKind, ONLY: Real8_
    use ModIonosphere, ONLY: clean_mod_ionosphere
    use ModConductance, ONLY: DoUseAurora, NameAuroraMod
    use ModIeRLM

    real,     intent(in) :: tSimulation   ! seconds from start time

    integer :: iFile
    real(Real8_) :: tCurrent

    character(len=*), parameter :: NameSub='IE_finalize'
    !--------------------------------------------------------------------------
    call get_time(tCurrentOut = tCurrent)
    call time_real_to_int(tCurrent, Time_Array)
    time_simulation = tSimulation

    if(nSolve>0)then
       do iFile=1,nFile
          if(iProc==0)      call ionosphere_write_output(iFile, 1)
          if(iProc==nProc-1)call ionosphere_write_output(iFile, 2)
       end do
    end if

    if(DoSaveLogfile .and. iProc==0 .and. unitlog>0) &
         call close_file(unitlog, NameCaller=NameSub)

    call clean_mod_ionosphere

    ! If legacy conductance model used, clean associated variables.
    if((DoUseAurora).and.(NameAuroraMod.eq.'FAC2FLUX'))then
       ! Hall conductance coeffs:
       deallocate(hal_a0_up, hal_a1_up, hal_a2_up, &
            hal_a0_do, hal_a1_do, hal_a2_do )
       ! Pedersen conductance coeffs:
       deallocate(ped_a0_up, ped_a1_up, ped_a2_up, &
            ped_a0_do, ped_a1_do, ped_a2_do )
       ! Coefficient grids:
       deallocate(cond_mlts, cond_lats)
    end if

  end subroutine IE_finalize
  !============================================================================

  subroutine IE_save_restart(tSimulation)

    use CON_coupler, ONLY: NameRestartOutDirComp
    use IE_ModIo,   ONLY: NameRestartOutDir

    real, intent(in) :: tSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IE_save_restart'
    !-------------------------------------------------------------------------
    if(NameRestartOutDirComp /= '') NameRestartOutDir = NameRestartOutDirComp

    call ionosphere_write_restart_file

  end subroutine IE_save_restart
  !============================================================================

  subroutine IE_run(tSimulation, tSimulationLimit)

    use ModProcIE
    use IE_ModMain
    use IE_ModIo,       ONLY: DoRestart, iUnitOut, StringPrefix
    use CON_physics,    ONLY: get_time, get_axes, time_real_to_int
    use ModLookupTable, ONLY: i_lookup_table, init_lookup_table, &
         interpolate_lookup_table
    use ModConductance, ONLY: f107_flux
    use ModKind

    real, intent(inout) :: tSimulation   ! current time of component

    real, intent(in) :: tSimulationLimit ! simulation time not to be exceeded

    real(Real8_) :: tStart, tNow
    real         :: tNowReal
    integer      :: nStep

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'IE_run'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,': iProc, tSimulation, tSimulationLimit=',&
         iProc, tSimulation, tSimulationLimit

    ! Store the current time
    time_simulation = tSimulation

    ! Since IE is not a time dependent component, it may advance to the
    ! next coupling time in a time accurate run
    if(time_accurate)tSimulation = tSimulationLimit

    if(DoTest)write(*,*)NameSub,': iProc, IsNewInput, DoRestart=', &
         iProc, IsNewInput, DoRestart

    ! Check if there has been information coming in
    if(.not.IsNewInput)then
       if(DoRestart)then
          ! Read restart files
          call ionosphere_read_restart_file
          if(DoTest)write(*,*) 'read restart done'

          ! Collect information onto proc 0
          call IE_gather
          if(DoTest)write(*,*) 'gather done'

          ! Restart only once
          DoRestart = .false.
       end if

       ! Do not solve if there is no new input
       RETURN
    end if

    ! Make sure that restart files are not used if there is new input
    DoRestart = .false.

    ! Check if we can have a reasonable magnetic field already
    call get_time(nStepOut=nStep)

    if(DoTest)write(*,*)NameSub,': iProc, nStep = ', iProc, nStep

    ! After the solve this input can be considered old
    IsNewInput = .false.

    ! Obtain the position of the magnetix axis
    call get_axes(time_simulation, MagAxisTiltGsmOut=ThetaTilt)

    ! Get the current time
    call get_time(tStartOut=tStart)
    tNow = tStart + time_simulation
    call time_real_to_int(tNow, Time_Array)

    ! If f107_flux is negative, use lookup table
    if(f107_flux < 0)then
       iTableF107 = i_lookup_table('F107')
       if(iTableF107 < 0)then
          ! Load default F107 table from Param/f107.txt
          write(iUnitOut,*) StringPrefix, ' loading Param/f107.txt'
          call init_lookup_table( &
               NameTable = 'F107', &
               NameCommand = 'load', &
               NameFile = 'Param/f107.txt', &
               TypeFile = 'log')
          iTableF107 = i_lookup_table('F107')
       end if
    end if

    ! get F10.7 from lookup table if available
    if(iTableF107 > 0)then
       tNowReal = tNow ! so it compiles with single precision
       call interpolate_lookup_table(iTableF107, tNowReal, f107_flux)
    end if

    if(f107_flux < 0) &
         call CON_stop(NameSub//': provide positive F10.7 value or table')

    nSolve = nSolve + 1

    if(DoTest)write(*,'(a,f8.3)') ' solve with F10.7=', f107_flux

    ! Solve for the ionosphere potential
    call IE_solve

    if(DoTest)write(*,*) 'done with solve'

    call IE_gather
    if(DoTest)write(*,*) 'gather done'

    ! Save solution (plot files) into file if required
    call IE_output

    ! Save logfile if required
    call IE_save_logfile
    if(DoTest)write(*,*) 'done with output'

    if(DoTest)write(*,*) 'done with IE_run'

  end subroutine IE_run
  !============================================================================

  subroutine IE_get_for_ps(Buffer_II, iSize, jSize, tSimulation)

    use ModNumConst,   ONLY: cRadToDeg
    use ModIonosphere, ONLY: IONO_nPsi, IONO_nTheta, &
         IONO_Phi, IONO_NORTH_Theta,IONO_NORTH_Psi

    integer, intent(in) :: iSize, jSize
    real,    intent(out):: Buffer_II(iSize,jSize)
    real,    intent(in) :: tSimulation

    integer :: i, j
    real    :: tSimulationTmp

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub = 'IE_get_for_ps'
    !--------------------------------------------------------------------------
    if(iSize /= IONO_nTheta*2-1 .or. jSize /= IONO_nPsi)then
       write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
            ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
       call CON_stop(NameSub//' ERROR in IE-PS coupling.')
    end if

    ! Make sure that the most recent result is provided
    tSimulationTmp = tSimulation
    call IE_run(tSimulationTmp, tSimulation)

    ! Pass potential to coupler:
    Buffer_II = IONO_Phi

    if(DoTestMe)then
       ! Write info to screen:
       write(*,*) "IE: Preparing potential for PS"
       write(*,'(a,f11.1,a,f11.1,a,f11.1,a)') "IE CPCP == ", &
            maxval(IONO_Phi), ' - ', minval(IONO_Phi), ' = ', &
            maxval(IONO_Phi)    -    minval(IONO_Phi)
       write(*,*) "IE: Size of pot array = ", size(Buffer_II), &
            size(Buffer_II,1), size(Buffer_II,2)
       ! Write potential to file:
       call open_file(FILE='ie_potential.txt', STATUS='replace', &
            NameCaller=NameSub)
       write(UnitTmp_,*)'Colat   Lon   Potential(V)'
       do i=1, iSize
          do j=1, jSize
             write(UnitTmp_, '(f6.2, 1x, f6.2, 1x, f9.1)') &
                  IONO_North_Theta(i,j)*cRadToDeg, &
                  IONO_North_Psi(i,j)*cRadToDeg, Buffer_II(i,j)
          end do
       end do
       call close_file(NameCaller=NameSub)
    end if

  end subroutine IE_get_for_ps
  !============================================================================
end module IE_wrapper
!==============================================================================
