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
INTEGER :: nw, ibinx, ibiny, ibinz, nAcc, nwj
REAL*8 :: r, x, y, z
REAL*8 :: rxiOx, ryiOx, rziOx, rxiH1, ryiH1, rziH1, rxiH2, ryiH2, rziH2
REAL*8 :: rxiDip, ryiDip, rziDip, magDip
REAL*8 :: exiDip, eyiDip, eziDip
REAL*8 :: dotDipi, riOx
REAL*8 :: rResXj, rResYj, rResZj, rResj
REAL*8 :: rxiOxp, ryiOxp, rziOxp, riOxp
REAL*8 :: rDifX, rDifY, rDifZ, rDif
REAL*8 :: rcutP1P2, p1, p2, p21, FastCut
REAL*8 :: rxOH1, ryOH1, rzOH1, rxOH2, ryOH2, rzOH2, p1OH1, p2OH1, p1OH2, p2OH2
REAL*8 :: exOH1, eyOH1, ezOH1, exOH2, eyOH2, ezOH2
REAL*8 :: rOH1, rOH2, dotOH1, dotOH2
REAL*8 :: cChiw, roh1n, roh2n, rohAn, rohBn
REAL*8 :: rABCrossx, rABCrossy, rABCrossz, rABCrossMag
REAL*8 :: nmCrossx, nmCrossy, nmCrossz, nmCrossMag
REAL*8 :: BrDif, Bnw, Baidbt2
INTEGER :: j, Bj
CHARACTER(LEN=200) :: hbFile
CHARACTER(LEN=10)  :: intString

   rcutP1P2 = 3.0d0
   rDifSave(:) = 0.0d0
   nwSave(1) = 0
   p1Acc(1) = 0.0d0
   p2Acc(1) = 0.0d0
   p21Acc(1) = 0.0d0
   ! THE WHOLE PURPOSE OF THIS DOUBLE LOOP IS TO FIND
   ! THE WATERS CLOSEST TO EACH ATOM IN A GIVEN FRAME.
   ! THE NUMBER OF WATERS FOUND SHOULD BE EQUAL TO THE
   ! NUMBER OF PROTEIN ATOMS.
   DO nw=1,nbeta2,3 ! loop over water oxygen atoms
      rxiOx =  rxbt2(nw)
      ryiOx =  rybt2(nw)
      rziOx =  rzbt2(nw)
      DO j=1,nbeta1 ! loop over beta=1 particles; could be atoms or residues; RIGHT NOW ATOMS
         rResXj = rxbt1(j)
         rResYj = rybt1(j)
         rResZj = rzbt1(j)
         rDifX = rxiOx - rResXj
         rDifY = ryiOx - rResYj
         rDifZ = rziOx - rResZj
         rDif = SQRT(rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ)
         IF (j.EQ.1) THEN  ! get first water difference and save
            BrDif = rDif
            Bnw = nw
            Baidbt2 = aidbt2(nw)
            Bj = j
         END IF
         IF (rDif.LT.BrDif.AND.rDif.LT.rcutP1P2) THEN
            BrDif = rDif         ! update difference if this (nw) water is closer to the jth atom.
            Bnw = nw             ! save water oxygen number (here)
            Baidbt2 = aidbt2(nw) ! save water oxygen index (PDB)
            Bj = j
         END IF
      END DO
!!
   FastCut=rcutP1P2+1
    IF (BrDif.LT.FastCut) THEN
!!
!!!!!!!!!!!!!!
      rxiH1 =  rxbt2(nw+1)
      ryiH1 =  rybt2(nw+1)
      rziH1 =  rzbt2(nw+1)
      rxiH2 =  rxbt2(nw+2)
      ryiH2 =  rybt2(nw+2)
      rziH2 =  rzbt2(nw+2)
      Do j=1,nbeta1
         rResXj = rxbt1(j)
         rResYj = rybt1(j)
         rResZj = rzbt1(j)
         rDifX = rxiH1 - rResXj
         rDifY = ryiH1 - rResYj
         rDifZ = rziH1 - rResZj
         rDif = SQRT(rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ)
         IF (rDif.LT.BrDif.AND.rDif.LT.rcutP1P2) THEN
            BrDif = rDif         ! update difference if this (nw) water is
            Bnw = nw             ! save water oxygen number (here)
            Baidbt2 = aidbt2(nw) ! save water oxygen index (PDB)
            Bj = j
         END IF
      END DO
!!
      Do j=1,nbeta1
         rResXj = rxbt1(j)
         rResYj = rybt1(j)
         rResZj = rzbt1(j)
         rDifX = rxiH2 - rResXj
         rDifY = ryiH2 - rResYj
         rDifZ = rziH2 - rResZj
         rDif = SQRT(rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ)
         IF (rDif.LT.BrDif.AND.rDif.LT.rcutP1P2) THEN
            BrDif = rDif         ! update difference if this (nw) water is
            Bnw = nw             ! save water oxygen number (here)
            Baidbt2 = aidbt2(nw) ! save water oxygen index (PDB)
            Bj = j
         END IF
      END DO

!!!!!!!!!!!!!!

           !END DO
           IF (BrDif.LT.rcutP1P2) THEN
              !DO j=1,nbeta1
              rResXj = rxbt1(Bj)
              rResYj = rybt1(Bj)
              rResZj = rzbt1(Bj)
              !IF (rDifSave(j).GT.rcutP1P2) THEN
                 !rDifSave(j) = 0.0d0
                 !nwSave(j) = 0
                 !aidbt2Save(j) = 0
                 !rDifX=0.0d0
                 !rDifY=0.0d0
                 !rDifZ=0.0d0
              !ELSE
                 nwSave(1) = nwSave(1) + 1
                 !nwj = nwSave(j)
                 !rxiOx =  rxbt2(nw)
                 !ryiOx =  rybt2(nw)
                 !rziOx =  rzbt2(nw)
                 riOx  = SQRT(rxiOx*rxiOx + ryiOx*ryiOx + rziOx*rziOx)
                 rDifX = rxiOx - rResXj
                 rDifY = ryiOx - rResYj
                 rDifZ = rziOx - rResZj
                 rDif  = SQRT(rDifX*rDifX + rDifY*rDifY + rDifZ*rDifZ)
                 !rxiH1 =  rxbt2(nw+1)
                 !ryiH1 =  rybt2(nw+1)
                 !rziH1 =  rzbt2(nw+1)
                 !rxiH2 =  rxbt2(nw+2)
                 !ryiH2 =  rybt2(nw+2)
                 !rziH2 =  rzbt2(nw+2)
                 rxOH1 = rxiH1 - rxiOx
                 ryOH1 = ryiH1 - ryiOx
                 rzOH1 = rziH1 - rziOx
                 rOH1  = SQRT(rxOH1*rxOH1 + ryOH1*ryOH1 + rzOH1*rzOH1)
                 exOH1 = rxOH1/rOH1
                 eyOH1 = ryOH1/rOH1
                 ezOH1 = rzOH1/rOH1
                 rxOH2 = rxiH2 - rxiOx
                 ryOH2 = ryiH2 - ryiOx
                 rzOH2 = rziH2 - rziOx
                 rOH2  = SQRT(rxOH2*rxOH2 + ryOH2*ryOH2 + rzOH2*rzOH2)
                 exOH2 = rxOH2/rOH2
                 eyOH2 = ryOH2/rOH2
                 ezOH2 = rzOH2/rOH2
                 !CHECK
                 !WRITE(*,*) ACOS(exOH1*exOH2+eyOH1*eyOH2+ezOH1*ezOH2)*180.0d0/3.14159d0
                 !WRITE(*,*) 'num of water out of :' , aidbt2(nw), nbeta2
                 rxiDip = qow*rxiOx + qhw*rxiH1 + qhw*rxiH2
                 ryiDip = qow*ryiOx + qhw*ryiH1 + qhw*ryiH2
                 rziDip = qow*rziOx + qhw*rziH1 + qhw*rziH2
                 magDip = SQRT(rxiDip*rxiDip + ryiDip*ryiDip + rziDip*rziDip)
                 exiDip = rxiDip/magDip
                 eyiDip = ryiDip/magDip
                 eziDip = rziDip/magDip
                 dotDipi = (exiDip*rDifX + eyiDip*rDifY + eziDip*rDifZ)/rDif
                 !dotOH1  = (exOH1*rDifX + eyOH1*rDifY + ezOH1*rDifZ)/rDif
                 !dotOH2  = (exOH2*rDifX + eyOH2*rDifY + ezOH2*rDifZ)/rDif
                 !!!!
                 roh1n = (exOH1*rDifx + eyOH1*rDify + ezOH1*rDifz)/rDif
                 roh2n = (exOH2*rDifx + eyOH2*rDify + ezOH2*rDifz)/rDif
                 IF (roh1n.GE.0.0d0 .AND. roh2n.LE.0.0d0) THEN
                    rohAn = roh1n
                    rohBn = roh2n
                    rABCrossx = eyOH1*ezOH2 - ezOH1*eyOH2
                    rABCrossy = ezOH1*exOH2 - exOH1*ezOH2
                    rABCrossz = exOH1*eyOH2 - eyOH1*exOH2
                 ELSE
                    rohAn = roh2n
                    rohBn = roh1n
                    rABCrossx = eyOH2*ezOH1 - ezOH2*eyOH1
                    rABCrossy = ezOH2*exOH1 - exOH2*ezOH1
                    rABCrossz = exOH2*eyOH1 - eyOH2*exOH1
                 END IF
                 nmCrossx = (rDifY*eziDip - rDifZ*eyiDip)/rDif
                 nmCrossy = (rDifZ*exiDip - rDifX*eziDip)/rDif
                 nmCrossz = (rDifX*eyiDip - rDifY*exiDip)/rDif
                 nmCrossMag = SQRT(nmCrossx*nmCrossx + nmCrossy*nmCrossy + nmCrossz*nmCrossz)
                 rABCrossMag = SQRT(rABCrossx*rABCrossx + rABCrossy*rABCrossy + rABCrossz*rABCrossz)
                 cChiw = (rABCrossx*nmCrossx + rABCrossy*nmCrossy + rABCrossz*nmCrossz)/nmCrossMag/rABCrossMag
                 !!!!
                 p1 = dotDipi
                 p2 = 0.50d0*(3.0d0*dotDipi*dotDipi - 1.0d0)
                 p21 = 0.50d0*((1.0d0-dotDipi*dotDipi)*(2.0d0*cChiw*cChiw-1.0d0))
                 WRITE(*,*) 'here is the p1 and p2 and p3:' , p1, p2 ,p21
                 p1Acc(1) = p1Acc(1) + p1
                 p2Acc(1) = p2Acc(1) + p2
                 p21Acc(1) = p21Acc(1) + p21
              END IF
           !END DO
!!
     END IF
!!
    END DO
RETURN
END SUBROUTINE p1p2Main
!-----------------------------------------------------------------------------------------
