#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

      SUBROUTINE DIAGNOSTICS_INTERP_VERT_OCEAN(
     I        bi,bj, k,
     I        w0,w1,gw_tmp,
     O        tmpFld,
     I        myThid)

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"

C !INPUT PARAMETERS: ===================================================
C  bi,bj                :: tile indices
C  k                    :: vertical level
C  KEscheme             :: spacial discretisation scheme for KE
C  uFld                 :: zonal flow
C  vFld                 :: meridional flow
C  KE                   :: Kinetic Energy
C  myThid               :: thread number
      INTEGER bi,bj,k
      _RL w0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL w1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL gw_tmp(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER KEscheme
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  KE                   :: Kinetic energy
      _RL tmpFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i, j         :: Loop counters
C     k            :: Level index
      INTEGER i,j
CEOP

      IF (K .LT.NR ) THEN
         DO j=1,sNy
           DO i=1,sNx
             tmpFld(i,j) = 0.25 _d 0
     &          *(w0(i,j,k,bi,bj)+w1(i,j,k,bi,bj))
     &          *gw_tmp(i,j,k,bi,bj)
     $          + 0.25 _d 0
     &          *(w0(i,j,k+1,bi,bj)+w1(i,j,k+1,bi,bj))
     &          *gw_tmp(i,j,k+1,bi,bj)
           ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
           DO i=1,sNx
             tmpFld(i,j) = 0.25 _d 0
     &          *(w0(i,j,k,bi,bj)+w1(i,j,k,bi,bj))
     &          *gw_tmp(i,j,k,bi,bj)
           ENDDO
         ENDDO

      ENDIF

      END
