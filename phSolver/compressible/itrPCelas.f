c-----------------------------------------------------------------------
c
c    Predict solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrPredictElas (disp)
      
      include "common.h"
      
      real*8   disp(nshg,nelas)
c
      disp = zero
c
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrectElas ( disp, Dy)
      
      include "common.h"
      
      real*8   disp(nshg,nelas), Dy(nshg,nelas)
c      
      disp = disp + Dy
c
      return
      end

