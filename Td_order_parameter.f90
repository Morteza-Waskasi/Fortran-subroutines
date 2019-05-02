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
INTEGER :: j, Bj
CHARACTER(LEN=200) :: hbFile
CHARACTER(LEN=10)  :: intString

   rcutP1P2 = 8.50d0
   FastCut = 9.50d0
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
         !IF (rDif.LT.BrDif.AND.rDif.LT.rcutP1P2) THEN
         IF (rDif.LT.BrDif.AND.rDif.LT.FastCut) THEN
            BrDif = rDif         ! update difference if this (nw) water is closer to the jth atom.
            Bnw = nw             ! save water oxygen number (here)
            Baidbt2 = aidbt2(nw) ! save water oxygen index (PDB)
            Bj = j
         END IF
      END DO
      !!
      !FastCut=rcutP1P2+1
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
!!!!!!!!!!!!!
                         r1=30d0
                         r2=30d0
                         r3=30d0
                         r4=30d0
                         DO nnw=1,nbeta2,3 ! loop over water oxygen atoms
                         IF (nnw.NE.nw) THEN
                             drxiOx =  rxbt2(nnw) - rxiOx
                             dryiOx =  rybt2(nnw) - ryiOx
                             drziOx =  rzbt2(nnw) - rziOx
                             driOx = SQRT(drxiOx*drxiOx + dryiOx*dryiOx + drziOx*drziOx)
                             IF (driOx.LT.r1) THEN
                                 IF (driOx.LT.r2) THEN
                                     IF (driOx.LT.r3) THEN
                                         IF (driOx.LT.r4) THEN
                                             n1 = n2
                                             r1 = r2
                                             n2 = n3
                                             r2 = r3
                                             n3 = n4
                                             r3 = r4
                                             n4 = nnw
                                             r4 = driOx
                                         ELSE
                                             n1 = n2
                                             r1 = r2
                                             n2 = n3
                                             r2 = r3
                                             n3 = nnw
                                             r3 = driOx
                                         END IF
                                     ELSE
                                         n1 = n2
                                         r1 = r2
                                         n2 = nnw
                                         r2 = driOx
                                     END IF
                                 ELSE
                                     n1 = nnw
                                     r1 = driOx
                                 END IF
                             END IF
                         END IF
                         END DO
                         Cosijk1 = (rxbt2(n1) - rxiOx)*(rxbt2(n2) - rxiOx)/(r1*r2) + & 
                          (rybt2(n1) - ryiOx)*(rybt2(n2) - ryiOx)/(r1*r2) + &
                          (rzbt2(n1) - rziOx)*(rzbt2(n2) - rziOx)/(r1*r2)
                         Cosijk2 = (rxbt2(n1) - rxiOx)*(rxbt2(n3) - rxiOx)/(r1*r3) + &
                          (rybt2(n1) - ryiOx)*(rybt2(n3) - ryiOx)/(r1*r3) + &
                          (rzbt2(n1) - rziOx)*(rzbt2(n3) - rziOx)/(r1*r3)
                         Cosijk3 = (rxbt2(n1) - rxiOx)*(rxbt2(n4) - rxiOx)/(r1*r4) + &
                          (rybt2(n1) - ryiOx)*(rybt2(n4) - ryiOx)/(r1*r4) + &
                          (rzbt2(n1) - rziOx)*(rzbt2(n4) - rziOx)/(r1*r4)
                         Cosijk4 = (rxbt2(n2) - rxiOx)*(rxbt2(n3) - rxiOx)/(r2*r3) + &
                          (rybt2(n2) - ryiOx)*(rybt2(n3) - ryiOx)/(r2*r3) + &
                          (rzbt2(n2) - rziOx)*(rzbt2(n3) - rziOx)/(r2*r3)
                         Cosijk5 = (rxbt2(n2) - rxiOx)*(rxbt2(n4) - rxiOx)/(r2*r4) + &
                          (rybt2(n2) - ryiOx)*(rybt2(n4) - ryiOx)/(r2*r4) + &
                          (rzbt2(n2) - rziOx)*(rzbt2(n4) - rziOx)/(r2*r4)
                         Cosijk6 = (rxbt2(n3) - rxiOx)*(rxbt2(n4) - rxiOx)/(r3*r4) + &
                          (rybt2(n3) - ryiOx)*(rybt2(n4) - ryiOx)/(r3*r4) + &
                          (rzbt2(n3) - rziOx)*(rzbt2(n4) - rziOx)/(r3*r4)
                         Qk = 1 - (3.0/8)*((cosijk1+1.0/3)*(cosijk1+1.0/3) + (cosijk2+1.0/3)*(cosijk2+1.0/3) + &
                          (cosijk3+1.0/3)*(cosijk3+1.0/3) + (cosijk4+1.0/3)*(cosijk4+1.0/3) + &
                          (cosijk5+1.0/3)*(cosijk5+1.0/3) + (cosijk6+1.0/3)*(cosijk6+1.0/3))
                         WRITE(*,*) Qk
!!!!!!!!!!!!!
                     END IF
                   !END DO
      !!
      END IF
      !!
   END DO
RETURN
END SUBROUTINE p1p2Main
!-----------------------------------------------------------------------------------------
