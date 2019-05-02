!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
SUBROUTINE p1p2Main
USE constants
USE atoms
USE coords
USE beta1
USE beta2
USE water
USE input
USE framesNum
USE files
USE p1p2
IMPLICIT NONE
INTEGER :: nw, ibinx, ibiny, ibinz, nAcc, nwj, nnw, n1, n2, n3, n4
REAL*8 :: r, x, y, z
REAL*8 :: rxiOx, ryiOx, rziOx, rxiH1, ryiH1, rziH1, rxiH2, ryiH2, rziH2, drxiOx, dryiOx, drziOx, driOx
REAL*8 :: r1, r2, r3, r4, Cosijk1, Cosijk2, Cosijk3, Cosijk4, Cosijk5, Cosijk6, Qk, ijkCos1
REAL*8 :: rxiDip, ryiDip, rziDip, magDip
REAL*8 :: exiDip, eyiDip, eziDip
REAL*8 :: dotDipi, riOx
REAL*8 :: rResXj, rResYj, rResZj, rResj
REAL*8 :: rxiOxp, ryiOxp, rziOxp, riOxp
REAL*8 :: rDifX, rDifY, rDifZ, rDif
REAL*8 :: rcutP1P2, p1, p2, p21 , FastCut
REAL*8 :: rxOH1, ryOH1, rzOH1, rxOH2, ryOH2, rzOH2, p1OH1, p2OH1, p1OH2, p2OH2
REAL*8 :: exOH1, eyOH1, ezOH1, exOH2, eyOH2, ezOH2
REAL*8 :: rOH1, rOH2, dotOH1, dotOH2
REAL*8 :: cChiw, roh1n, roh2n, rohAn, rohBn
REAL*8 :: rABCrossx, rABCrossy, rABCrossz, rABCrossMag
REAL*8 :: nmCrossx, nmCrossy, nmCrossz, nmCrossMag
REAL*8 :: BrDif, Bnw, Baidbt2
REAL*8 :: EPSLON,SIGMA, SR2, SR6, SR12, V, VIJ, V2,r0
REAL*8 ::  SR20, SR60, SR120, VIJ0
INTEGER :: j, Bj
CHARACTER(LEN=200) :: hbFile
CHARACTER(LEN=10)  :: intString


   EPSLON=0.104260731d0
   SIGMA= 3.75d0
   r0= (( 2**(1.0/6.0)) *SIGMA)
   rcutP1P2 = 18.0d0
   FastCut = 7.0d0
   rDifSave(:) = 0.0d0
   nwSave(1) = 0
   p1Acc(1) = 0.0d0
   p2Acc(1) = 0.0d0
   p21Acc(1) = 0.0d0


! Perturbation treatment of LJ interactions
! According to WCA perturbation approach, short range repulsive part is VIJ0 (Huang+Chandler, JPCB 2002, eq1 ):
             SR20   = (SIGMA**2) / (r0**2)
             SR60= SR20 * SR20 * SR20
             SR120= SR60**2
             VIJ0   = SR120 - SR60



         V2=0.0d0 
   DO j=1,nbeta1 ! loop over beta=1 particles; could be atoms or residues; RIGHT NOW ATOMS
      rResXj = rxbt1(j)
      rResYj = rybt1(j)
      rResZj = rzbt1(j)
         V=0.0d0 
         DO nw=1,nbeta2!,3 ! loop over water oxygen atoms
            rxiOx =  rxbt2(nw)
            ryiOx =  rybt2(nw)
            rziOx =  rzbt2(nw)

            rDifX = rxiOx - rResXj
            rDifY = ryiOx - rResYj
            rDifZ = rziOx - rResZj
            rDif = SQRT(rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ)

            IF (rDif.LT.rcutP1P2.AND.rDif.GT.r0 ) THEN
!               print *,rDif
! Perturbation treatment of LJ interactions
! According to WCA perturbation approach, longer-ranged and often slowly varying attractive  part is VIJ (Huang+Chandler, JPCB 2002, eq1 ):
             SR2   = (SIGMA**2) / (rDif**2)
             SR6= SR2 * SR2 * SR2
             SR12= SR6**2
             VIJ   = SR12 - SR6
             V= V+VIJ           
             !WRITE(*,*) rDif

            ELSE IF (rDif.LT.rcutP1P2.AND.rDif.LE.r0 ) THEN

!             SR2   = (SIGMA**2) / (r0**2)
!             SR6= SR2 * SR2 * SR2
!             SR12= SR6**2
!             VIJ   = SR12 - SR6
             V= V+VIJ0           

            END IF 
         END DO
         V2= V2+ V  

   END DO
             WRITE(*,*)'vdw=', V2*4*EPSLON

RETURN
END SUBROUTINE p1p2Main
!-----------------------------------------------------------------------------------------
