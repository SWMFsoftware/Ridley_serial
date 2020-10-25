!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

Module ModConductance
  implicit none
  save

  integer :: longmx, latmx, ndx
  real    :: steplat

  integer, parameter :: latmdx = 40
  integer, parameter :: lonmdx = 30
  integer, parameter :: mndx   = 10
  integer, parameter :: nConductanceSolutions = 4

  real, dimension(mndx) :: halmin, pedmin, avk50, efx50

  real, dimension(0:lonmdx,0:latmdx,mndx) :: &
       halar, pedar, avkar, efxar

  real, dimension(nConductanceSolutions) :: ConductanceBackground

  real, dimension(4) :: bkgc, bkgcer
  real               :: flxmax, dbymax, dbzmax

  integer, parameter :: pedersen_ = 1
  integer, parameter :: hall_     = 2
  integer, parameter :: eflux_    = 3
  integer, parameter :: avee_     = 4

  ! Auroral oval config for iModel=5:
  logical :: UseOval=.true., DoOvalShift=.true., &
       UseSubOvalCond=.false., DoFitCircle=.true.

  ! Option to turn on new oval:
  logical :: UseNewOval = .false.

  ! Floor values for GM density and pressure, SI units:
  real, parameter :: GmRhoFloor = 1E-21, GmPFloor = 1E-13

contains
  !===========================================================================
  subroutine smooth_lagrange_polar(a_II, iLatSize, jLonSize, tolIn)
    ! Use a simple sliding window technique to smooth a 2D array in polar
    ! coordinates such that periodicity and continuity are enforced across
    ! the LON=360/0 boundary. The max and min latitude points are not
    ! considered.  Only points that are more than *tol* percent greater than
    ! the average of their neighbors are changed.
    !
    ! Arguments
    !
    ! a_II(iLatSize,jLonSize): A two-dimentional array to smooth in place.
    ! iLatSize, jLonSize     : The size of a_II in lat/lon space

    integer, intent(in)        :: iLatSize, jLonSize
    real, intent(inout)        :: a_II(iLatSize, jLonSize)
    real, intent(in), optional :: tolIn

    integer                    :: i, j
    real                       :: aSmooth_II(iLatSize, jLonSize), tol
    !------------------------------------------------------------------------

    ! Set default tolerance to 15% of original value:
    if (present(tolIn)) then
       tol = tolIn
    else
       tol=0.15
    endif

    ! Construct smoothed values
    ! No smoothing at poles:
    aSmooth_II(1,:) = a_II(1,:)
    aSmooth_II(iLatSize,:) = a_II(iLatSize,:)
    ! Smooth all non-pole points:
    colat: do i = 2, iLatSize-1
       ! Perform averaging at longitude boundaries (azimuthal periodicity)
       aSmooth_II(i,1) = (                  & ! At lon=0 degrees
            sum(a_II(i-1:i+1, jLonSize)) +  &
            sum(a_II(i-1:i+1, 2))        +  &
            a_II(i-1,1) + a_II(i+1,1)       &
            ) / 8.0
       
       aSmooth_II(i,jLonSize) = (              & ! At lon=360 degrees
            sum(a_II(i-1:i+1, jLonSize -1)) +  &
            sum(a_II(i-1:i+1, 1))           +  &
            a_II(i-1,jLonSize)              +  &
            a_II(i+1,jLonSize) ) / 8.0

       ! Smooth over all non-longitude boundary points:
       lon: do j = 2, jLonSize-1
          aSmooth_II(i,j) = (             & ! For all other (i,j)
               sum(a_II(i-1:i+1, j-1)) +  & !left of point i,j
               sum(a_II(i-1:i+1, j+1)) +  & !right of point i,j
               a_II(i+1,j) + a_II(i-1,j)  & !above and below
               )/ 8.0
       end do lon
       
    end do colat

    ! Apply smoothing to discrete_nf only where smoothing changes values
    ! by tolerance percentage *tol* or more from their original value:
    where (a_II<=ABS(aSmooth_II - a_II) /tol) a_II = aSmooth_II
    
  end subroutine smooth_lagrange_polar
  !===========================================================================
  
end Module ModConductance
