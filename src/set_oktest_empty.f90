!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine set_oktest(str,DoTest,DoTestMe)

  implicit none

  character (len=*) :: str
  logical :: DoTest, DoTestMe

  !----------------------------------------------------------------------------
  DoTest = .true.
  DoTestMe = .true.

  ! This routine does nothing.

end subroutine set_oktest
!==============================================================================
