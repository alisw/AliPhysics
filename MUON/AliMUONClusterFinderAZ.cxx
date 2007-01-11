/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// -------------------------------
// Class AliMUONClusterFinderAZ
// -------------------------------
// Clusterizer class based on the Expectation-Maximization algorithm
// Author: Alexander Zinchenko, JINR Dubna

#include <stdlib.h>
#include <Riostream.h>
#include <TH2.h>
#include <TMinuit.h>
#include <TMatrixD.h>
#include <TRandom.h>

#include "AliMUONClusterFinderAZ.h"
#include "AliMUONClusterDrawAZ.h"
#include "AliMUONVGeometryDESegmentation.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONClusterInput.h"
#include "AliMUONPixel.h"
#include "AliMUONMathieson.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONClusterFinderAZ)
/// \endcond
 
 const Double_t AliMUONClusterFinderAZ::fgkCouplMin = 1.e-3; // threshold on coupling 
 const Double_t AliMUONClusterFinderAZ::fgkZeroSuppression = 6; // average zero suppression value
 const Double_t AliMUONClusterFinderAZ::fgkSaturation = 3000; // average saturation level
 AliMUONClusterFinderAZ* AliMUONClusterFinderAZ::fgClusterFinder = 0x0;
 TMinuit* AliMUONClusterFinderAZ::fgMinuit = 0x0;
//FILE *lun1 = fopen("nxny.dat","w");

//_____________________________________________________________________________
AliMUONClusterFinderAZ::AliMUONClusterFinderAZ(Bool_t draw)
  : AliMUONClusterFinderVS(),
    fZpad(0),
    fNpar(0),
    fQtot(0),
    fReco(1),
    fCathBeg(0),
    fDraw(0x0),
    fPixArray(0x0),
    fnCoupled(0),
    fDebug(0)
{
/// Constructor
  fnPads[0]=fnPads[1]=0;
  
  for (Int_t i=0; i<7; i++)
    for (Int_t j=0; j<fgkDim; j++)
      fXyq[i][j]= 9999.;

  for (Int_t i=0; i<4; i++)
    for (Int_t j=0; j<fgkDim; j++) 
      fPadIJ[i][j]=-1;

  for (Int_t i=0; i<2; i++)
    for (Int_t j=0; j<fgkDim; j++) 
      fUsed[i][j] = 0;

  fSegmentation[1] = fSegmentation[0] = 0x0; 

  fPadBeg[0] = fPadBeg[1] = 0;

  if (!fgMinuit) fgMinuit = new TMinuit(8);
  if (!fgClusterFinder) fgClusterFinder = this;
  fPixArray = new TObjArray(20); 

  if (draw) {
    fDebug = 1;
    fReco = 0;
    fDraw = new AliMUONClusterDrawAZ(this);
  }
  cout << " *** Running AZ cluster finder *** " << endl;
}

//_____________________________________________________________________________
AliMUONClusterFinderAZ::~AliMUONClusterFinderAZ()
{
/// Destructor
  delete fgMinuit; fgMinuit = 0; delete fPixArray; fPixArray = 0;
  delete fDraw;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FindRawClusters()
{
/// To provide the same interface as in AliMUONClusterFinderVS

  ResetRawClusters(); 
  EventLoop (fEvtNumber, fInput->Chamber());
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::EventLoop(Int_t nev, Int_t ch)
{
/// Loop over digits
  
  if (fDraw && !fDraw->FindEvCh(nev, ch)) return;

  fSegmentation[0] = (AliMUONVGeometryDESegmentation*) fInput->
    Segmentation2(0)->GetDESegmentation(fInput->DetElemId());
  fSegmentation[1] = (AliMUONVGeometryDESegmentation*) fInput->
    Segmentation2(1)->GetDESegmentation(fInput->DetElemId());
    
  Int_t ndigits[2] = {9,9}, nShown[2] = {0};
  if (fReco != 2) { // skip initialization for the combined cluster / track
    fCathBeg = fPadBeg[0] = fPadBeg[1] = 0;
    for (Int_t i = 0; i < 2; i++) {
      for (Int_t j = 0; j < fgkDim; j++) { fUsed[i][j] = kFALSE; }
    }
  }

next:
  if (fReco == 2 && (nShown[0] || nShown[1])) return; // only one precluster for the combined finder
  if (ndigits[0] == nShown[0] && ndigits[1] == nShown[1]) return;

  Float_t xpad, ypad, zpad, zpad0;
  Bool_t first = kTRUE;
  if (fDebug) cout << " *** Event # " << nev << " chamber: " << ch << endl;
  fnPads[0] = fnPads[1] = 0;
  for (Int_t i = 0; i < fgkDim; i++) fPadIJ[1][i] = 0; 

  for (Int_t iii = fCathBeg; iii < 2; iii++) { 
    Int_t cath = TMath::Odd(iii);
    ndigits[cath] = fInput->NDigits(cath); 
    if (!ndigits[0] && !ndigits[1]) return;
    if (ndigits[cath] == 0) continue;
    if (fDebug) cout << " ndigits: " << ndigits[cath] << " " << cath << endl;

    AliMUONDigit  *mdig;
    Int_t digit;

    Bool_t eEOC = kTRUE; // end-of-cluster
    for (digit = fPadBeg[cath]; digit < ndigits[cath]; digit++) {
      mdig = AliMUONClusterInput::Instance()->Digit(cath,digit); 
      if (first) {
	// Find first unused pad
	if (fUsed[cath][digit]) continue;
	//if (!fSegmentation[cath]->GetPadC(fInput->DetElemId(),mdig->PadX(),mdig->PadY(),xpad,ypad,zpad0)) { 
	if (!fSegmentation[cath]->HasPad(mdig->PadX(), mdig->PadY())) {
	  // Handle "non-existing" pads
	  fUsed[cath][digit] = kTRUE; 
	  continue; 
	} 
	fSegmentation[cath]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad0); 
      } else {
	if (fUsed[cath][digit]) continue;
	//if (!fSegmentation[cath]->GetPadC(fInput->DetElemId(),mdig->PadX(),mdig->PadY(),xpad,ypad,zpad)) {
	if (!fSegmentation[cath]->HasPad(mdig->PadX(), mdig->PadY())) {
	  // Handle "non-existing" pads
	  fUsed[cath][digit] = kTRUE; 
	  continue; 
	} 
	fSegmentation[cath]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad); 
	//if (TMath::Abs(zpad-zpad0) > 0.1) continue; // different slats
	// Find a pad overlapping with the cluster
	if (!Overlap(cath,mdig)) continue;
      }
      // Add pad - recursive call
      AddPad(cath,digit);
      //AZ !!!!!! Temporary fix of St1 overlap regions !!!!!!!! 
      /*
      if (cath && ch < 2) {
	Int_t npads = fnPads[0] + fnPads[1] - 1;
	Int_t cath1 = fPadIJ[0][npads];
	Int_t idig = TMath::Nint (fXyq[5][npads]);
	mdig = AliMUONClusterInput::Instance()->Digit(cath1,idig);
	//fSegmentation[cath1]->GetPadC(fInput->DetElemId(),mdig->PadX(),mdig->PadY(),xpad,ypad,zpad);
	fSegmentation[cath1]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad);
	if (TMath::Abs(zpad-zpad0) > 0.1) zpad0 = zpad;
      }	
      */
      eEOC = kFALSE;
      if (digit >= 0) break;
    }
    if (first && eEOC) {
      // No more unused pads 
      if (cath == 0) continue; // on cathode #0 - check #1
      else return; // No more clusters
    }
    if (eEOC) break; // cluster found
    first = kFALSE;
    if (fDebug) cout << " nPads: " << fnPads[cath] << " " << nShown[cath]+fnPads[cath] << " " << cath << endl;
  } // for (Int_t iii = 0;

  fZpad = zpad0;
  if (fDraw) fDraw->DrawCluster(); 

  // Use MLEM for cluster finder
  Int_t nMax = 1, localMax[100], maxPos[100];
  Double_t maxVal[100];
  
  if (CheckPrecluster(nShown)) {
    BuildPixArray();
    //*
    if (fnPads[0]+fnPads[1] > 50) nMax = FindLocalMaxima(fPixArray, localMax, maxVal);
    if (nMax > 1) TMath::Sort(nMax, maxVal, maxPos, kTRUE); // in decreasing order
    Int_t iSimple = 0, nInX = -1, nInY;
    PadsInXandY(nInX, nInY);
    if (fDebug) cout << "Pads in X and Y: " << nInX << " " << nInY << endl;
    if (nMax == 1 && nInX < 4 && nInY < 4) iSimple = 1; //1; // simple cluster
    //*/
    /* For test
    Int_t iSimple = 0, nInX = -1, nInY;
    PadsInXandY(nInX, nInY);
    if (fDebug) cout << "Pads in X and Y: " << nInX << " " << nInY << endl;
    if (nMax == 1 && nInX < 4 && nInY < 4) iSimple = 1; //1; // simple cluster
    if (!iSimple) nMax = FindLocalMaxima(fPixArray, localMax, maxVal);
    nMax = 1;
    if (nMax > 1) TMath::Sort(nMax, maxVal, maxPos, kTRUE); // in decreasing order
    */
    for (Int_t i=0; i<nMax; i++) {
      if (nMax > 1) FindCluster(localMax, maxPos[i]);
      MainLoop(iSimple);
      if (i < nMax-1) {
	for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
	  if (fPadIJ[1][j] == 0) continue; // pad charge was not modified
	  fPadIJ[1][j] = 0;
	  fXyq[2][j] = fXyq[6][j]; // use backup charge value
	}
      }
    } // for (Int_t i=0; i<nMax;
    if (nMax > 1) ((TH2D*) gROOT->FindObject("anode"))->Delete();
    TH2D *mlem = (TH2D*) gROOT->FindObject("mlem");
    if (mlem) mlem->Delete();
  }
  if (!fDraw || fDraw->Next()) goto next;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddPad(Int_t cath, Int_t digit)
{
/// Add pad to the cluster

  AliMUONDigit *mdig = fInput->Digit(cath,digit); 

  Float_t charge = mdig->Signal();
  // get the center of the pad
  Float_t xpad, ypad, zpad0;
  //if (!fSegmentation[cath]->GetPadC(fInput->DetElemId(),mdig->PadX(),mdig->PadY(),xpad,ypad,zpad0)) {	  // Handle "non-existing" pads
  if (!fSegmentation[cath]->HasPad(mdig->PadX(), mdig->PadY())) {
    fUsed[cath][digit] = kTRUE; 
    return; 
  } 
  fSegmentation[cath]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad0);
  Int_t isec = fSegmentation[cath]->Sector(mdig->PadX(), mdig->PadY());
  Int_t nPads = fnPads[0] + fnPads[1];
  fXyq[0][nPads] = xpad;
  fXyq[1][nPads] = ypad;
  fXyq[2][nPads] = charge;
  fXyq[3][nPads] = fSegmentation[cath]->Dpx(isec)/2;
  fXyq[4][nPads] = fSegmentation[cath]->Dpy(isec)/2;
  fXyq[5][nPads] = digit;
  fXyq[6][nPads] = 0;
  fPadIJ[0][nPads] = cath;
  fPadIJ[1][nPads] = 0;
  fPadIJ[2][nPads] = mdig->PadX();
  fPadIJ[3][nPads] = mdig->PadY();
  fUsed[cath][digit] = kTRUE;
  if (fDebug) printf(" bbb %d %d %f %f %f %f %f %f %3d %3d \n", nPads, cath, 
                     xpad, ypad, zpad0, fXyq[3][nPads]*2, fXyq[4][nPads]*2, 
                     charge, mdig->PadX(), mdig->PadY());
  fnPads[cath]++;

  // Check neighbours
  Int_t nn, ix, iy, xList[10], yList[10];
  AliMUONDigit  *mdig1;

  Int_t ndigits = fInput->NDigits(cath); 
  fSegmentation[cath]->Neighbours(mdig->PadX(), mdig->PadY(), &nn, xList, yList); 
  for (Int_t in = 0; in < nn; in++) {
    ix = xList[in];
    iy = yList[in];
    for (Int_t digit1 = 0; digit1 < ndigits; digit1++) {
      if (digit1 == digit) continue;
      mdig1 = fInput->Digit(cath,digit1); 
      if (!fUsed[cath][digit1] && mdig1->PadX() == ix && mdig1->PadY() == iy) {
	fUsed[cath][digit1] = kTRUE;
	// Add pad - recursive call
	AddPad(cath,digit1);
      }
    } //for (Int_t digit1 = 0;
  } // for (Int_t in = 0;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::Overlap(Int_t cath, AliMUONDigit *mdig)
{
/// Check if the pad from one cathode overlaps with a pad 
/// in the precluster on the other cathode

  Float_t xpad, ypad, zpad;
  fSegmentation[cath]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad);
  Int_t isec = fSegmentation[cath]->Sector(mdig->PadX(), mdig->PadY());

  Float_t xy1[4], xy12[4];
  xy1[0] = xpad - fSegmentation[cath]->Dpx(isec)/2;
  xy1[1] = xy1[0] + fSegmentation[cath]->Dpx(isec);
  xy1[2] = ypad - fSegmentation[cath]->Dpy(isec)/2;
  xy1[3] = xy1[2] + fSegmentation[cath]->Dpy(isec);
  //cout << " ok " << fnPads[0]+fnPads[1] << xy1[0] << xy1[1] << xy1[2] << xy1[3] << endl;

  Int_t cath1 = TMath::Even(cath);
  for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
    if (fPadIJ[0][i] != cath1) continue;
    if (Overlap(xy1, i, xy12, 0)) return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::Overlap(Float_t *xy1, Int_t iPad, Float_t *xy12, Int_t iSkip)
{
/// Check if the pads xy1 and iPad overlap and return overlap area

  Float_t xy2[4];
  xy2[0] = fXyq[0][iPad] - fXyq[3][iPad];
  xy2[1] = fXyq[0][iPad] + fXyq[3][iPad];
  if (xy1[0] > xy2[1]-1.e-4 || xy1[1] < xy2[0]+1.e-4) return kFALSE;
  xy2[2] = fXyq[1][iPad] - fXyq[4][iPad];
  xy2[3] = fXyq[1][iPad] + fXyq[4][iPad];
  if (xy1[2] > xy2[3]-1.e-4 || xy1[3] < xy2[2]+1.e-4) return kFALSE;
  if (!iSkip) return kTRUE; // just check overlap (w/out computing the area)
  xy12[0] = TMath::Max (xy1[0],xy2[0]);
  xy12[1] = TMath::Min (xy1[1],xy2[1]);
  xy12[2] = TMath::Max (xy1[2],xy2[2]);
  xy12[3] = TMath::Min (xy1[3],xy2[3]);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::CheckPrecluster(Int_t *nShown)
{
/// Check precluster in order to attempt to simplify it (mostly for
/// two-cathode preclusters)

  Int_t i1, i2, cath=0, digit=0;
  Float_t xy1[4], xy12[4];
  
  Int_t npad = fnPads[0] + fnPads[1];
  if (npad == 1) { 
    // Disregard one-pad clusters (leftovers from splitting)
    nShown[0] += fnPads[0]; 
    nShown[1] += fnPads[1]; 
    return kFALSE;
  }

  // If pads have the same size take average of pads on both cathodes 
  //Int_t sameSize = (fnPads[0] && fnPads[1]) ? 1 : 0;
  Int_t sameSize = 0; //AZ - 17-01-06
  
  if (sameSize) {
    Double_t xSize = -1, ySize = 0;
    for (Int_t i=0; i<npad; i++) {
      if (fXyq[2][i] < 0) continue;
      if (xSize < 0) { xSize = fXyq[3][i]; ySize = fXyq[4][i]; }
      if (TMath::Abs(xSize-fXyq[3][i]) > 1.e-4 ||  TMath::Abs(ySize-fXyq[4][i]) > 1.e-4) { sameSize = 0; break; }
    }
  } // if (sameSize)
  if (sameSize && fnPads[0] == 1 && fnPads[1] == 1) sameSize = 0; //AZ
  // Handle shift by half a pad in Station 1
  if (sameSize) {
    Int_t cath0 = fPadIJ[0][0];
    for (Int_t i = 1; i < npad; i++) {
      if (fPadIJ[0][i] == cath0) continue;
      Double_t dx = TMath::Abs ((fXyq[0][i] - fXyq[0][0]) / fXyq[3][i] / 2);
      Int_t idx = (Int_t) TMath::Abs ((fXyq[0][i] - fXyq[0][0]) / fXyq[3][i] / 2);
      if (TMath::Abs (dx - idx) > 0.001) sameSize = 0;
      break;
    }
  } // if (sameSize)
   
  if (sameSize && (fnPads[0] >= 2 || fnPads[1] >= 2)) {
    nShown[0] += fnPads[0];
    nShown[1] += fnPads[1];
    fnPads[0] = fnPads[1] = 0;
    Int_t div;
    for (Int_t i=0; i<npad; i++) {
      if (fXyq[2][i] < 0) continue; // used pad
      fXyq[2][fnPads[0]] = fXyq[2][i];
      div = 1;
      cath = fPadIJ[0][i];
      for (Int_t j=i+1; j<npad; j++) {
	if (fPadIJ[0][j] == fPadIJ[0][i]) continue; // same cathode
	if (TMath::Abs(fXyq[0][j]-fXyq[0][i]) > 1.e-4) continue;
	if (TMath::Abs(fXyq[1][j]-fXyq[1][i]) > 1.e-4) continue;
	fXyq[2][fnPads[0]] += fXyq[2][j];
	div = 2;
	fXyq[2][j] = -2;
	if (cath) fXyq[5][fnPads[0]] = fXyq[5][j]; // save digit number for cath 0
	break;
      }
      // Flag that the digit from the other cathode
      if (cath && div == 1) fXyq[5][fnPads[0]] = -fXyq[5][i] - 1; 
      // If low pad charge take the other equal to 0 
      //if (div == 1 && fXyq[2][fnPads[0]] < fgkZeroSuppression + 1.5*3) div = 2;
      fXyq[2][fnPads[0]] /= div;
      fXyq[0][fnPads[0]] = fXyq[0][i];
      fXyq[1][fnPads[0]] = fXyq[1][i];
      fPadIJ[2][fnPads[0]] = fPadIJ[2][i];
      fPadIJ[3][fnPads[0]] = fPadIJ[3][i];
      fPadIJ[0][fnPads[0]++] = 0;
    }
  } // if (sameSize)

  // Check if one-cathode precluster
  i1 = fnPads[0]!=0 ? 0 : 1;
  i2 = fnPads[1]!=0 ? 1 : 0;

  if (i1 != i2) { // two-cathode 

    Int_t *flags = new Int_t[npad];
    for (Int_t i=0; i<npad; i++) { flags[i] = 0; }

    // Check pad overlaps
    for (Int_t i=0; i<npad; i++) {
      if (fPadIJ[0][i] != i1) continue;
      xy1[0] = fXyq[0][i] - fXyq[3][i];
      xy1[1] = fXyq[0][i] + fXyq[3][i];
      xy1[2] = fXyq[1][i] - fXyq[4][i];
      xy1[3] = fXyq[1][i] + fXyq[4][i];
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[0][j] != i2) continue;
	if (!Overlap(xy1, j, xy12, 0)) continue;
	flags[i] = flags[j] = 1; // mark overlapped pads
      } // for (Int_t j=0;
    } // for (Int_t i=0;

    // Check if all pads overlap
    Int_t nFlags=0;
    for (Int_t i=0; i<npad; i++) {
      if (flags[i]) continue;
      nFlags ++;
      if (fDebug) cout << i << " " << fPadIJ[0][i] << " " << fXyq[0][i] << " " << fXyq[1][i] << endl;
    }
    if (fDebug && nFlags) cout << " nFlags = " << nFlags << endl;
    //if (nFlags > 2 || (Float_t)nFlags / npad > 0.2) { // why 2 ??? - empirical choice
    if (nFlags > 0) {
      for (Int_t i=0; i<npad; i++) {
	if (flags[i]) continue;
	digit = TMath::Nint (fXyq[5][i]);
	cath = fPadIJ[0][i];
	// Check for edge effect (missing pads on the other cathode)
	Int_t cath1 = TMath::Even(cath), ix, iy;
	ix = iy = 0;
        //if (!fSegmentation[cath1]->GetPadI(fInput->DetElemId(),fXyq[0][i],fXyq[1][i],fZpad,ix,iy)) continue;
	if (!fSegmentation[cath1]->HasPad(fXyq[0][i], fXyq[1][i], fZpad)) continue;
	if (nFlags == 1 && fXyq[2][i] < fgkZeroSuppression * 3) continue;
	fUsed[cath][digit] = kFALSE; // release pad
	fXyq[2][i] = -2;
	fnPads[cath]--;
      }
      if (fDraw) fDraw->UpdateCluster(npad);
    } // if (nFlags > 2)

    // Check correlations of cathode charges
    if (fnPads[0] && fnPads[1]) { // two-cathode
      Double_t sum[2]={0};
      Int_t over[2] = {1, 1};
      for (Int_t i=0; i<npad; i++) {
	cath = fPadIJ[0][i];
	if (fXyq[2][i] > 0) sum[cath] += fXyq[2][i];
	if (fXyq[2][i] > fgkSaturation-1) over[cath] = 0;
      }
      if (fDebug) cout << " Total charge: " << sum[0] << " " << sum[1] << endl;
      if ((over[0] || over[1]) && TMath::Abs(sum[0]-sum[1])/(sum[0]+sum[1])*2 > 1) { // 3 times difference
	if (fDebug) cout << " Release " << endl;
	// Big difference
	cath = sum[0] > sum[1] ? 0 : 1;
	Int_t imax = 0, imin = 0;
	Double_t cmax = -1, cmin = 9999, dxMin = 0, dyMin = 0;
	Double_t *dist = new Double_t[npad];
	for (Int_t i = 0; i < npad; i++) {
	  if (fPadIJ[0][i] != cath || fXyq[2][i] < 0) continue;
	  if (fXyq[2][i] < cmin) {
	    cmin = fXyq[2][i];
	    imin = i;
	  }
	  if (fXyq[2][i] < cmax) continue;
	  cmax = fXyq[2][i];
	  imax = i;
	}
	// Arrange pads according to their distance to the max, 
	// normalized to the pad size
	for (Int_t i = 0; i < npad; i++) {
	  dist[i] = 0;
	  if (fPadIJ[0][i] != cath || fXyq[2][i] < 0) continue;
	  if (i == imax) continue; 
	  Double_t dx = (fXyq[0][i] - fXyq[0][imax]) / fXyq[3][imax] / 2;
	  Double_t dy = (fXyq[1][i] - fXyq[1][imax]) / fXyq[4][imax] / 2;
	  dist[i] = TMath::Sqrt (dx * dx + dy * dy);
	  if (i == imin) {
	    cmin = dist[i] + 0.001; // distance to the pad with minimum charge 
	    dxMin = dx;
	    dyMin = dy;
	  }
	}
	TMath::Sort(npad, dist, flags, kFALSE); // in increasing order
	Int_t indx;
	Double_t xmax = -1;
	for (Int_t i = 0; i < npad; i++) {
	  indx = flags[i];
	  if (fPadIJ[0][indx] != cath || fXyq[2][indx] < 0) continue;
	  if (dist[indx] > cmin) {
	    // Farther than the minimum pad
	    Double_t dx = (fXyq[0][indx] - fXyq[0][imax]) / fXyq[3][imax] / 2;
	    Double_t dy = (fXyq[1][indx] - fXyq[1][imax]) / fXyq[4][imax] / 2;
	    dx *= dxMin;
	    dy *= dyMin;
	    if (dx >= 0 && dy >= 0) continue;
	    if (TMath::Abs(dx) > TMath::Abs(dy) && dx >= 0) continue;
	    if (TMath::Abs(dy) > TMath::Abs(dx) && dy >= 0) continue;
	  }
	  if (fXyq[2][indx] <= cmax || TMath::Abs(dist[indx]-xmax) < 1.e-3) {
	    // Release pads
	    if (TMath::Abs(dist[indx]-xmax) < 1.e-3) 
                cmax = TMath::Max((Double_t)(fXyq[2][indx]),cmax);
	    else cmax = fXyq[2][indx];
	    xmax = dist[indx];
	    digit = TMath::Nint (fXyq[5][indx]);
	    fUsed[cath][digit] = kFALSE; 
	    fXyq[2][indx] = -2;
	    fnPads[cath]--;
	  } 
	} // for (Int_t i = 0; i < npad;

	// Check pad overlaps once more
	for (Int_t j = 0; j < npad; j++) flags[j] = 0; 
	for (Int_t k = 0; k < npad; k++) {
	  if (fXyq[2][k] < 0 || fPadIJ[0][k] != i1) continue;
	  xy1[0] = fXyq[0][k] - fXyq[3][k];
	  xy1[1] = fXyq[0][k] + fXyq[3][k];
	  xy1[2] = fXyq[1][k] - fXyq[4][k];
	  xy1[3] = fXyq[1][k] + fXyq[4][k];
	  for (Int_t j = 0; j < npad; j++) {
	    if (fXyq[2][j] < 0) continue;
	    if (fPadIJ[0][j] != i2) continue;
	    if (!Overlap(xy1, j, xy12, 0)) continue;
	    flags[k] = flags[j] = 1; // mark overlapped pads
	  } // for (Int_t j = 0;
	} // for (Int_t k = 0;
	nFlags = 0;
	for (Int_t j = 0; j < npad; j++) {
	  if (fXyq[2][j] < 0 || flags[j]) continue;
	  nFlags ++;
	}
	if (nFlags == fnPads[0] + fnPads[1]) {
	  // No overlap
	  for (Int_t j = 0; j < npad; j++) {
	    if (fXyq[2][j] < 0 || fPadIJ[0][j] != cath) continue;
	    fXyq[2][j] = -2;
	    fnPads[cath]--;
	  }
	}
	delete [] dist; dist = 0;
	if (fDraw) fDraw->UpdateCluster(npad);
      } // TMath::Abs(sum[0]-sum[1])...
    } // if (fnPads[0] && fnPads[1])
    delete [] flags; flags = 0;
  } // if (i1 != i2) 

  if (!sameSize) { nShown[0] += fnPads[0]; nShown[1] += fnPads[1]; }

  // Move released pads to the right
  Int_t beg = 0, end = npad-1, padij;
  Double_t xyq;
  while (beg < end) {
    if (fXyq[2][beg] > 0) { beg++; continue; }
    for (Int_t j=end; j>beg; j--) {
      if (fXyq[2][j] < 0) continue;
      end = j - 1;
      for (Int_t j1=0; j1<4; j1++) {
	padij = fPadIJ[j1][beg]; 
	fPadIJ[j1][beg] = fPadIJ[j1][j];
	fPadIJ[j1][j] = padij;
      }
      for (Int_t j1=0; j1<6; j1++) {
	xyq = fXyq[j1][beg]; 
	fXyq[j1][beg] = fXyq[j1][j];
	fXyq[j1][j] = xyq;
      }
      break;
    } // for (Int_t j=end;
    beg++;
  } // while
  npad = fnPads[0] + fnPads[1];
  if (npad > 500) { 
    AliWarning(Form(" *** Too large cluster. Give up. %d ", npad));
    return kFALSE; 
  }
  // Back up charge value
  for (Int_t j = 0; j < npad; j++) fXyq[6][j] = fXyq[2][j];

  return kTRUE;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::BuildPixArray()
{
/// Build pixel array for MLEM method
  
  Int_t nPix=0, i1, i2;
  Float_t xy1[4], xy12[4];
  AliMUONPixel *pixPtr=0;

  Int_t npad = fnPads[0] + fnPads[1];

  // One cathode is empty
  i1 = fnPads[0]!=0 ? 0 : 1;
  i2 = fnPads[1]!=0 ? 1 : 0;

  // Build array of pixels on anode plane
  if (i1 == i2) { // one-cathode precluster
    for (Int_t j=0; j<npad; j++) {
      pixPtr = new AliMUONPixel();
      for (Int_t i=0; i<2; i++) {
	pixPtr->SetCoord(i, fXyq[i][j]); // pixel coordinates
	pixPtr->SetSize(i, fXyq[i+3][j]); // pixel size
      }
      pixPtr->SetCharge(fXyq[2][j]); // charge
      fPixArray->Add((TObject*)pixPtr);
      nPix++;
    }
  } else { // two-cathode precluster    
    i1 = fPadIJ[0][0];
    i2 = TMath::Even (i1);
    for (Int_t i = 0; i < npad; i++) {
      if (fPadIJ[0][i] != i1) continue;
      xy1[0] = fXyq[0][i] - fXyq[3][i];
      xy1[1] = fXyq[0][i] + fXyq[3][i];
      xy1[2] = fXyq[1][i] - fXyq[4][i];
      xy1[3] = fXyq[1][i] + fXyq[4][i];
      for (Int_t j = 1; j < npad; j++) {
	if (fPadIJ[0][j] != i2) continue;
	if (!Overlap(xy1, j, xy12, 1)) continue;
	pixPtr = new AliMUONPixel();
	for (Int_t k=0; k<2; k++) {
	  pixPtr->SetCoord(k, (xy12[2*k]+xy12[2*k+1])/2); // pixel coordinates
	  pixPtr->SetSize(k, xy12[2*k+1]-pixPtr->Coord(k)); // size
	}
	pixPtr->SetCharge(TMath::Min (fXyq[2][i],fXyq[2][j])); //charge
	fPixArray->Add((TObject*)pixPtr);
	//cout << nPix << " " << pixPtr->Coord(0) << " " << pixPtr->Size(0) << " " << pixPtr->Coord(1) << " " << pixPtr->Size(1) << " " << pixPtr->Charge() << endl;
	nPix++;
      } // for (Int_t j=0;
    } // for (Int_t i=0;
  } // else

  Float_t xPadMin = 999, yPadMin = 999;
  for (Int_t i = 0; i < npad; i++) {
    xPadMin = TMath::Min (xPadMin, fXyq[3][i]);
    yPadMin = TMath::Min (yPadMin, fXyq[4][i]);
  }
  if (fDebug) cout << xPadMin << " " << yPadMin << endl;

  Float_t wxmin = 999, wymin = 999;
  for (Int_t i = 0; i < nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    wxmin = TMath::Min ((Double_t)wxmin, pixPtr->Size(0));
    wymin = TMath::Min ((Double_t)wymin, pixPtr->Size(1));
  }
  if (fDebug) cout << wxmin << " " << wymin << endl;
  wxmin = TMath::Abs (wxmin - xPadMin/2) > 0.001 ? xPadMin : xPadMin / 2;
  wymin = TMath::Abs (wymin - yPadMin/2) > 0.001 ? yPadMin : yPadMin / 2;
  //wxmin = xPadMin; wymin = yPadMin;

  // Check if small pixel X-size
  AdjustPixel(wxmin, 0);
  // Check if small pixel Y-size
  AdjustPixel(wymin, 1);
  // Check if large pixel size
  AdjustPixel(wxmin, wymin);

  // Remove discarded pixels
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    //pixPtr->Print();
    if (pixPtr->Charge() < 1) { fPixArray->RemoveAt(i); delete pixPtr; }// discarded pixel
  }
  fPixArray->Compress();
  nPix = fPixArray->GetEntriesFast();

  if (nPix > npad) {
    if (fDebug) cout << nPix << endl;
    // Too many pixels - sort and remove pixels with the lowest signal
    fPixArray->Sort();
    for (Int_t i=npad; i<nPix; i++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
      //pixPtr->Print();
      fPixArray->RemoveAt(i);
      delete pixPtr;
    }
    nPix = npad;
  } // if (nPix > npad)

  // Set pixel charges to the same value (for MLEM)
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    //pixPtr->SetCharge(10);
    if (fDebug) cout << i+1 << " " << pixPtr->Coord(0) << " " << pixPtr->Coord(1) << " " << pixPtr->Size(0) << " " << pixPtr->Size(1) << endl;
  }
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AdjustPixel(Float_t width, Int_t ixy)
{
/// Check if some pixels have small size (adjust if necessary)

  AliMUONPixel *pixPtr, *pixPtr1 = 0;
  Int_t ixy1 = TMath::Even(ixy);
  Int_t nPix = fPixArray->GetEntriesFast();

  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    if (pixPtr->Size(ixy)-width < -1.e-4) {
      // try to merge 
      if (fDebug) cout << i << " Small X or Y: " << ixy << " " << pixPtr->Size(ixy) << " " << width << " " << pixPtr->Coord(0) << " " << pixPtr->Coord(1) << endl;
      for (Int_t j=i+1; j<nPix; j++) {
	pixPtr1 = (AliMUONPixel*) fPixArray->UncheckedAt(j);
	if (pixPtr1->Charge() < 1) continue; // discarded pixel
	if (TMath::Abs(pixPtr1->Size(ixy)-width) < 1.e-4) continue; // right size 
	if (TMath::Abs(pixPtr1->Coord(ixy1)-pixPtr->Coord(ixy1)) > 1.e-4) continue; // different rows/columns
	if (TMath::Abs(pixPtr1->Coord(ixy)-pixPtr->Coord(ixy)) < 2*width) {
	  // merge
	  Double_t tmp = pixPtr->Coord(ixy) + pixPtr1->Size(ixy) * 
	    TMath::Sign (1., pixPtr1->Coord(ixy) - pixPtr->Coord(ixy));
	  pixPtr->SetCoord(ixy, tmp);
	  pixPtr->SetSize(ixy, width);
	  pixPtr->SetCharge(TMath::Min (pixPtr->Charge(),pixPtr1->Charge()));
	  pixPtr1->SetCharge(0);
	  pixPtr1 = 0;
	  break;
	}
      } // for (Int_t j=i+1;
      //if (!pixPtr1) { cout << " I am here!" << endl; pixPtr->SetSize(ixy, width); } // ???
      //else if (pixPtr1->Charge() > 0.5 || i == nPix-1) {
      if (pixPtr1 || i == nPix-1) {
	// edge pixel - just increase its size
	if (fDebug) cout << " Edge ..." << endl; 
	for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
          //if (fPadIJ[0][j] != ixy1) continue;
	  //???-check if (TMath::Abs(pixPtr->Coord(ixy1)-fXyq[ixy1][j]) > 1.e-4) continue;
	  if (pixPtr->Coord(ixy) < fXyq[ixy][j]) 
	    //pixPtr->Shift(ixy, -pixPtr->Size(ixy));
	    pixPtr->Shift(ixy, pixPtr->Size(ixy)-width);
	  //else pixPtr->Shift(ixy, pixPtr->Size(ixy));
	  else pixPtr->Shift(ixy, -pixPtr->Size(ixy)+width);
  	  pixPtr->SetSize(ixy, width);
	  break;
	}
      }
    } // if (pixPtr->Size(ixy)-width < -1.e-4)
  } // for (Int_t i=0; i<nPix;
  return;
}
  
//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AdjustPixel(Float_t wxmin, Float_t wymin)
{
/// Check if some pixels have large size (adjust if necessary)

  Int_t n1[2], n2[2], iOK = 1, nPix = fPixArray->GetEntriesFast();
  AliMUONPixel *pixPtr, pix;
  Double_t xy0[2] = {9999, 9999}, wxy[2], dist[2] = {0};

  // Check if large pixel size
  for (Int_t i = 0; i < nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    if (pixPtr->Size(0) - wxmin < 1.e-4) {
      if (xy0[0] > 9998) xy0[0] = pixPtr->Coord(0); // position of a "normal" pixel
      if (pixPtr->Size(1) - wymin < 1.e-4) { 
	if (xy0[1] > 9998) xy0[1] = pixPtr->Coord(1); // position of a "normal" pixel
	continue;
      } else iOK = 0; // large pixel
    } else {
      iOK = 0; // large pixel
      if (xy0[1] > 9998 && pixPtr->Size(1) - wymin < 1.e-4) xy0[1] = pixPtr->Coord(1); // "normal" pixel
    }      
    if (xy0[0] < 9998 && xy0[1] < 9998) break;
  }
  if (iOK) return;

  wxy[0] = wxmin;
  wxy[1] = wymin;
  //cout << xy0[0] << " " << xy0[1] << endl;
  for (Int_t i = 0; i < nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    n1[0] = n1[1] = 999;
    n2[0] = n2[1] = 1;
    for (Int_t j = 0; j < 2; j++) {
      if (pixPtr->Size(j) - wxy[j] < 1.e-4) continue;
      dist[j] = (pixPtr->Coord(j) - xy0[j]) / wxy[j] / 2; // normalized distance to "normal" pixel
      n2[j] = TMath::Nint (pixPtr->Size(j) / wxy[j]);
      n1[j] = n2[j] == 1 ? TMath::Nint(dist[j]) : (Int_t)dist[j];
    }
    if (n1[0] > 998 && n1[1] > 998) continue;
    if (fDebug) cout << " Different " << pixPtr->Size(0) << " " << wxy[0] << " "
		     << pixPtr->Size(1) << " " << wxy[1] <<endl;
    
    if (n2[0] > 2 || n2[1] > 2) { 
      //cout << n2[0] << " " << n2[1] << endl; 
      if (n2[0] > 2 && n1[0] < 999) n1[0]--;
      if (n2[1] > 2 && n1[1] < 999) n1[1]--;
    }
    //cout << n1[0] << " " << n2[0] << " " << n1[1] << " " << n2[1] << endl;
    pix = *pixPtr;
    pix.SetSize(0, wxy[0]); pix.SetSize(1, wxy[1]);
    //pixPtr->Print();
    for (Int_t ii = 0; ii < n2[0]; ii++) {
      if (n1[0] < 999) pix.SetCoord(0, xy0[0] + (n1[0] + TMath::Sign(1.,dist[0]) * ii) * 2 * wxy[0]);
      for (Int_t jj = 0; jj < n2[1]; jj++) {
	if (n1[1] < 999) pix.SetCoord(1, xy0[1] + (n1[1] + TMath::Sign(1.,dist[1]) * jj) * 2 * wxy[1]);
	fPixArray->Add(new AliMUONPixel(pix));
	//pix.Print();
      }
    }
    pixPtr->SetCharge(0);
  } // for (Int_t i = 0; i < nPix;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::MainLoop(Int_t iSimple)
{
/// Repeat MLEM algorithm until pixel size becomes sufficiently small
  
  TH2D *mlem;

  Int_t ix, iy;
  //Int_t nn, xList[10], yList[10];
  Int_t nPix = fPixArray->GetEntriesFast();
  AliMUONPixel *pixPtr = 0;
  Double_t *coef = 0, *probi = 0; 
  AddVirtualPad(); // add virtual pads if necessary
  Int_t npadTot = fnPads[0] + fnPads[1], npadOK = 0;
  for (Int_t i = 0; i < npadTot; i++) if (fPadIJ[1][i] == 0) npadOK++;
  if (fDraw) fDraw->ResetMuon();

  while (1) {

    mlem = (TH2D*) gROOT->FindObject("mlem");
    if (mlem) mlem->Delete();
    // Calculate coefficients
    if (fDebug) cout << " nPix, npadTot, npadOK " << nPix << " " << npadTot << " " << npadOK << endl;

    // Calculate coefficients and pixel visibilities
    coef = new Double_t [npadTot*nPix];
    probi = new Double_t [nPix];
    for (Int_t ipix=0; ipix<nPix; ipix++) probi[ipix] = 0;
    Int_t indx = 0, indx1 = 0, cath = 0;

    for (Int_t j=0; j<npadTot; j++) {
      indx = j*nPix;
      if (fPadIJ[1][j] == 0) {
	cath = fPadIJ[0][j];
	ix = fPadIJ[2][j];
	iy = fPadIJ[3][j];
	fSegmentation[cath]->SetPad(ix, iy);
	/*
	  fSegmentation[cath]->Neighbours(fInput->DetElemId(),ix,iy,&nn,xList,yList); 
	  if (nn != 4) {
	  cout << nn << ": ";
	  for (Int_t i=0; i<nn; i++) {cout << xList[i] << " " << yList[i] << ", ";}
	  cout << endl;
	  }
	*/
      }

      for (Int_t ipix=0; ipix<nPix; ipix++) {
	indx1 = indx + ipix;
	if (fPadIJ[1][j] < 0) { coef[indx1] = 0; continue; }
	pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
	fSegmentation[cath]->SetHit(pixPtr->Coord(0), pixPtr->Coord(1), fZpad);
	coef[indx1] = fInput->Mathieson()->IntXY(fInput->DetElemId(),fInput->Segmentation2(cath));
	probi[ipix] += coef[indx1];
      } // for (Int_t ipix=0;
    } // for (Int_t j=0;
    for (Int_t ipix=0; ipix<nPix; ipix++) if (probi[ipix] < 0.01) pixPtr->SetCharge(0); // "invisible" pixel

    // MLEM algorithm
    Mlem(coef, probi, 15);

    Double_t xylim[4] = {999, 999, 999, 999};
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      //cout << ipix+1; pixPtr->Print();
      for (Int_t i=0; i<4; i++) 
	xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
    }
    for (Int_t i=0; i<4; i++) {
      xylim[i] -= pixPtr->Size(i/2); if (fDebug) cout << (i%2 ? -1 : 1)*xylim[i] << " "; }
    if (fDebug) cout << endl;

    // Adjust histogram to approximately the same limits as for the pads
    // (for good presentation)
    if (fDraw) fDraw->AdjustHist(xylim, pixPtr);

    Int_t nx = TMath::Nint ((-xylim[1]-xylim[0])/pixPtr->Size(0)/2);
    Int_t ny = TMath::Nint ((-xylim[3]-xylim[2])/pixPtr->Size(1)/2);
    
    mlem = new TH2D("mlem","mlem",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      mlem->Fill(pixPtr->Coord(0),pixPtr->Coord(1),pixPtr->Charge());
    }
    if (fDraw) fDraw->DrawHist("c2", mlem);

    // Check if the total charge of pixels is too low
    Double_t qTot = 0;
    for (Int_t i=0; i<nPix; i++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
      qTot += pixPtr->Charge();
    }
    if (qTot < 1.e-4 || npadOK < 3 && qTot < 7) {
      delete [] coef; delete [] probi; coef = 0; probi = 0;
      fPixArray->Delete(); 
      for (Int_t i=0; i<npadTot; i++) if (fPadIJ[1][i] == 0) fPadIJ[1][i] = -1;
      return kFALSE; 
    }

    // Plot data - expectation
    /*
    Double_t x, y, cont;
    for (Int_t j=0; j<npadTot; j++) {
      Double_t sum1 = 0;
      for (Int_t i=0; i<nPix; i++) {
        // Caculate expectation
	pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
	sum1 += pixPtr->Charge()*coef[j*nPix+i];
      }
      sum1 = TMath::Min (sum1,fgkSaturation);
      x = fXyq[0][j];
      y = fXyq[1][j];
      cath = fPadIJ[0][j];
      Int_t ihist = cath*2;
      ix = fHist[ihist]->GetXaxis()->FindBin(x);
      iy = fHist[ihist]->GetYaxis()->FindBin(y);
      cont = fHist[ihist]->GetCellContent(ix,iy);
      if (cont == 0 && fHist[ihist+1]) {
        ihist += 1;
	ix = fHist[ihist]->GetXaxis()->FindBin(x);
	iy = fHist[ihist]->GetYaxis()->FindBin(y);
      }
      fHist[ihist]->SetBinContent(ix,iy,fXyq[2][j]-sum1);
    }
    ((TCanvas*)gROOT->FindObject("c1"))->cd(1);
    //gPad->SetTheta(55);
    //gPad->SetPhi(30);
    //mlem->Draw("lego1");
    gPad->Modified();
    ((TCanvas*)gROOT->FindObject("c1"))->cd(2);
    gPad->Modified();
    */

    if (iSimple) {
      // Simple cluster - skip further passes thru EM-procedure
      Simple();
      delete [] coef; delete [] probi; coef = 0; probi = 0;
      fPixArray->Delete(); 
      return kTRUE;
    }

    // Calculate position of the center-of-gravity around the maximum pixel
    Double_t xyCOG[2];
    FindCOG(mlem, xyCOG);

    if (TMath::Min(pixPtr->Size(0),pixPtr->Size(1)) < 0.07 && pixPtr->Size(0) > pixPtr->Size(1)) break;
    //if (TMath::Min(pixPtr->Size(0),pixPtr->Size(1)) < 0.007 && pixPtr->Size(0) > pixPtr->Size(1)) break;
    //if (TMath::Min(pixPtr->Size(0),pixPtr->Size(1)) >= 0.07 || pixPtr->Size(0) < pixPtr->Size(1)) {
    // Sort pixels according to the charge
    fPixArray->Sort();
    /*
    for (Int_t i=0; i<nPix; i++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
      cout << i+1; pixPtr->Print();
    }
    */
    Double_t pixMin = 0.01*((AliMUONPixel*)fPixArray->UncheckedAt(0))->Charge();
    pixMin = TMath::Min (pixMin,50.);

    // Decrease pixel size and shift pixels to make them centered at 
    // the maximum one
    indx = (pixPtr->Size(0)>pixPtr->Size(1)) ? 0 : 1;
    Double_t width = 0, shift[2]={0};
    ix = 1;
    for (Int_t i=0; i<4; i++) xylim[i] = 999;
    Int_t nPix1 = nPix; nPix = 0;
    for (Int_t ipix=0; ipix<nPix1; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      if (nPix >= npadOK) { // too many pixels already
	fPixArray->RemoveAt(ipix); 
	delete pixPtr; 
	continue;
      }
      if (pixPtr->Charge() < pixMin) { // low charge
	fPixArray->RemoveAt(ipix); 
	delete pixPtr; 
	continue;
      }
      for (Int_t i=0; i<2; i++) {
	if (!i) {
	  pixPtr->SetCharge(10);
	  pixPtr->SetSize(indx, pixPtr->Size(indx)/2);
	  width = -pixPtr->Size(indx);
	  pixPtr->Shift(indx, width);
	  // Shift pixel position
	  if (ix) {
	    ix = 0;
	    for (Int_t j=0; j<2; j++) {
	      shift[j] = pixPtr->Coord(j) - xyCOG[j];
	      shift[j] -= ((Int_t)(shift[j]/pixPtr->Size(j)/2))*pixPtr->Size(j)*2;
	    }
	    //cout << ipix << " " << i << " " << shift[0] << " " << shift[1] << endl;
	  } // if (ix)
	  pixPtr->Shift(0, -shift[0]);
	  pixPtr->Shift(1, -shift[1]);
	} else {
	  pixPtr = new AliMUONPixel(*pixPtr);
	  pixPtr->Shift(indx, -2*width);
	  fPixArray->Add((TObject*)pixPtr);
	} // else
	//pixPtr->Print();
	for (Int_t i=0; i<4; i++) 
	  xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
      } // for (Int_t i=0; i<2;
      nPix += 2;
    } // for (Int_t ipix=0;

    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();

    // Remove excessive pixels
    if (nPix > npadOK) {
      for (Int_t ipix=npadOK; ipix<nPix; ipix++) { 
	pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
	fPixArray->RemoveAt(ipix); 
	delete pixPtr;
      }
    } else {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(0);
      // add pixels if the maximum is at the limit of pixel area
      // start from Y-direction
      Int_t j = 0;
      for (Int_t i=3; i>-1; i--) {
	if (nPix < npadOK && 
	    TMath::Abs((i%2 ? -1 : 1)*xylim[i]-xyCOG[i/2]) < pixPtr->Size(i/2)) {
	  pixPtr = new AliMUONPixel(*pixPtr);
	  pixPtr->SetCoord(i/2, xyCOG[i/2]+(i%2 ? 2:-2)*pixPtr->Size(i/2));
	  j = TMath::Even (i/2);
	  pixPtr->SetCoord(j, xyCOG[j]);
	  fPixArray->Add((TObject*)pixPtr);
	  nPix++;
	}
      }
    } // else    

    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();
    delete [] coef; delete [] probi; coef = 0; probi = 0;
  } // while (1)

  // remove pixels with low signal or low visibility
  // Cuts are empirical !!!
  Double_t thresh = TMath::Max (mlem->GetMaximum()/100.,1.);
  thresh = TMath::Min (thresh,50.);
  Double_t cmax = -1, charge = 0;
  for (Int_t i=0; i<nPix; i++) cmax = TMath::Max (cmax,probi[i]); 
  //cout << thresh << " " << cmax << " " << cmax*0.9 << endl;
  // Mark pixels which should be removed
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    charge = pixPtr->Charge();
    if (charge < thresh) pixPtr->SetCharge(-charge);
    //else if (cmax > 1.91) {
    //  if (probi[i] < 1.9) pixPtr->SetCharge(-charge);
    //}
    //AZ else if (probi[i] < cmax*0.9) pixPtr->SetCharge(-charge);
    //18-01-06 else if (probi[i] < cmax*0.8) pixPtr->SetCharge(-charge);
    //cout << i << " " << pixPtr->Coord(0) << " " << pixPtr->Coord(1) << " " << charge << " " << probi[i] << endl;
  }
  // Move charge of removed pixels to their nearest neighbour (to keep total charge the same)
  Int_t near = 0;
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    charge = pixPtr->Charge();
    if (charge > 0) continue;
    near = FindNearest(pixPtr);
    pixPtr->SetCharge(0);
    probi[i] = 0; // make it "invisible"
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(near);
    pixPtr->SetCharge(pixPtr->Charge() + (-charge));
  }
  Mlem(coef,probi,2);
  // Update histogram
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    ix = mlem->GetXaxis()->FindBin(pixPtr->Coord(0));
    iy = mlem->GetYaxis()->FindBin(pixPtr->Coord(1));
    mlem->SetBinContent(ix, iy, pixPtr->Charge());
  }
  if (fDraw) fDraw->DrawHist("c2", mlem);

  // Try to split into clusters
  Bool_t ok = kTRUE;
  if (mlem->GetSum() < 1) ok = kFALSE;
  else Split(mlem, coef);
  delete [] coef; delete [] probi; coef = 0; probi = 0;
  fPixArray->Delete(); 
  return ok;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Mlem(Double_t *coef, Double_t *probi, Int_t nIter)
{
/// Use MLEM to find pixel charges
  
  Int_t nPix = fPixArray->GetEntriesFast();
  Int_t npad = fnPads[0] + fnPads[1];
  Double_t *probi1 = new Double_t [nPix];
  Double_t probMax = 0;
  Int_t indx, indx1;
  AliMUONPixel *pixPtr;

  for (Int_t ipix=0; ipix<nPix; ipix++) if (probi[ipix] > probMax) probMax = probi[ipix];
  for (Int_t iter=0; iter<nIter; iter++) {
    // Do iterations
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      // Correct each pixel
      if (probi[ipix] < 0.01) continue; // skip "invisible" pixel
      Double_t sum = 0;
      //probi1[ipix] = probi[ipix];
      probi1[ipix] = probMax;
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[1][j] < 0) continue; 
	Double_t sum1 = 0;
	indx1 = j*nPix;
	indx = indx1 + ipix;
	for (Int_t i=0; i<nPix; i++) {
	  // Caculate expectation
	  pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
	  sum1 += pixPtr->Charge()*coef[indx1+i];
	} // for (Int_t i=0;
	if (fXyq[2][j] > fgkSaturation-1 && sum1 > fXyq[2][j]) { probi1[ipix] -= coef[indx]; continue; } // correct for pad charge overflows
	//cout << sum1 << " " << fXyq[2][j] << " " << coef[j*nPix+ipix] << endl;
	if (coef[indx] > 1.e-6) sum += fXyq[2][j]*coef[indx]/sum1;
      } // for (Int_t j=0;
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      if (probi1[ipix] > 1.e-6) pixPtr->SetCharge(pixPtr->Charge()*sum/probi1[ipix]);
    } // for (Int_t ipix=0;
  } // for (Int_t iter=0;
  delete [] probi1;
  return;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FindCOG(TH2D *mlem, Double_t *xyc)
{
/// Calculate position of the center-of-gravity around the maximum pixel

  Int_t ixmax, iymax, ix, nsumx=0, nsumy=0, nsum=0;
  Int_t i1 = -9, j1 = -9;
  mlem->GetMaximumBin(ixmax,iymax,ix);
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t thresh = mlem->GetMaximum()/10;
  Double_t x, y, cont, xq=0, yq=0, qq=0;

  for (Int_t i=TMath::Max(1,iymax-1); i<=TMath::Min(ny,iymax+1); i++) {
    y = mlem->GetYaxis()->GetBinCenter(i);
    for (Int_t j=TMath::Max(1,ixmax-1); j<=TMath::Min(nx,ixmax+1); j++) {
      cont = mlem->GetCellContent(j,i);
      if (cont < thresh) continue;
      if (i != i1) {i1 = i; nsumy++;}
      if (j != j1) {j1 = j; nsumx++;}
      x = mlem->GetXaxis()->GetBinCenter(j);
      xq += x*cont;
      yq += y*cont;
      qq += cont;
      nsum++;
    }
  }

  Double_t cmax = 0;
  Int_t i2 = 0, j2 = 0;
  x = y = 0;
  if (nsumy == 1) {
    // one bin in Y - add one more (with the largest signal)
    for (Int_t i=TMath::Max(1,iymax-1); i<=TMath::Min(ny,iymax+1); i++) {
      if (i == iymax) continue;
      for (Int_t j=TMath::Max(1,ixmax-1); j<=TMath::Min(nx,ixmax+1); j++) {
	cont = mlem->GetCellContent(j,i);
	if (cont > cmax) {
	  cmax = cont;
	  x = mlem->GetXaxis()->GetBinCenter(j);
	  y = mlem->GetYaxis()->GetBinCenter(i);
	  i2 = i;
	  j2 = j;
	}
      }
    }
    xq += x*cmax;
    yq += y*cmax;
    qq += cmax;
    if (i2 != i1) nsumy++;
    if (j2 != j1) nsumx++;
    nsum++;
  } // if (nsumy == 1)

  if (nsumx == 1) {
    // one bin in X - add one more (with the largest signal)
    cmax = x = y = 0;
    for (Int_t j=TMath::Max(1,ixmax-1); j<=TMath::Min(nx,ixmax+1); j++) {
      if (j == ixmax) continue;
      for (Int_t i=TMath::Max(1,iymax-1); i<=TMath::Min(ny,iymax+1); i++) {
	cont = mlem->GetCellContent(j,i);
	if (cont > cmax) {
	  cmax = cont;
	  x = mlem->GetXaxis()->GetBinCenter(j);
	  y = mlem->GetYaxis()->GetBinCenter(i);
	  i2 = i;
	  j2 = j;
	}
      }
    }
    xq += x*cmax;
    yq += y*cmax;
    qq += cmax;
    if (i2 != i1) nsumy++;
    if (j2 != j1) nsumx++;
    nsum++;
  } // if (nsumx == 1)

  xyc[0] = xq/qq; xyc[1] = yq/qq;
  if (fDebug) cout << xyc[0] << " " << xyc[1] << " " << qq << " " << nsum << " " << nsumx << " " << nsumy << endl;
  return;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::FindNearest(AliMUONPixel *pixPtr0)
{
/// Find the pixel nearest to the given one
/// (algorithm may be not very efficient)

  Int_t nPix = fPixArray->GetEntriesFast(), imin = 0;
  Double_t rmin = 99999, dx = 0, dy = 0, r = 0;
  Double_t xc = pixPtr0->Coord(0), yc = pixPtr0->Coord(1);
  AliMUONPixel *pixPtr;

  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 0.5) continue;
    dx = (xc - pixPtr->Coord(0)) / pixPtr->Size(0);
    dy = (yc - pixPtr->Coord(1)) / pixPtr->Size(1);
    r = dx *dx + dy * dy;
    if (r < rmin) { rmin = r; imin = i; }
  }
  return imin;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Split(TH2D *mlem, Double_t *coef)
{
/// The main steering function to work with clusters of pixels in anode
/// plane (find clusters, decouple them from each other, merge them (if
/// necessary), pick up coupled pads, call the fitting function)
  
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Int_t nPix = fPixArray->GetEntriesFast();

  Bool_t *used = new Bool_t[ny*nx];
  Double_t cont;
  Int_t nclust = 0, indx, indx1;

  for (Int_t i=0; i<ny*nx; i++) used[i] = kFALSE; 

  TObjArray *clusters[200]={0};
  TObjArray *pix;

  // Find clusters of histogram bins (easier to work in 2-D space)
  for (Int_t i=1; i<=ny; i++) {
    for (Int_t j=1; j<=nx; j++) {
      indx = (i-1)*nx + j - 1;
      if (used[indx]) continue;
      cont = mlem->GetCellContent(j,i);
      if (cont < 0.5) continue;
      pix = new TObjArray(20);
      used[indx] = 1;
      pix->Add(BinToPix(mlem,j,i));
      AddBin(mlem, i, j, 0, used, pix); // recursive call
      if (nclust >= 200) AliFatal(" Too many clusters !!!");
      clusters[nclust++] = pix;
    } // for (Int_t j=1; j<=nx; j++) {
  } // for (Int_t i=1; i<=ny;
  if (fDebug) cout << nclust << endl;
  delete [] used; used = 0;
  
  // Compute couplings between clusters and clusters to pads
  Int_t npad = fnPads[0] + fnPads[1];

  // Write out some information for algorithm development
  Int_t cath=0, npadx[2]={0}, npady[2]={0};
  Double_t xlow[2]={9999,9999}, xhig[2]={-9999,-9999};
  Double_t ylow[2]={9999,9999}, yhig[2]={-9999,-9999};
  for (Int_t j=0; j<npad; j++) {
    if (fXyq[3][j] < 0) continue; // exclude virtual pads
    cath = fPadIJ[0][j];
    if (fXyq[0][j] < xlow[cath]-0.001) { 
      if (fXyq[0][j]+fXyq[3][j] <= xlow[cath] && npadx[cath]) npadx[cath]++;
      xlow[cath] = fXyq[0][j];
    }
    if (fXyq[0][j] > xhig[cath]+0.001) { 
      if (fXyq[0][j]-fXyq[3][j] >= xhig[cath]) npadx[cath]++; 
      xhig[cath] = fXyq[0][j]; 
    }
    if (fXyq[1][j] < ylow[cath]-0.001) { 
      if (fXyq[1][j]+fXyq[4][j] <= ylow[cath] && npady[cath]) npady[cath]++;
      ylow[cath] = fXyq[1][j];
    }
    if (fXyq[1][j] > yhig[cath]+0.001) { 
      if (fXyq[1][j]-fXyq[4][j] >= yhig[cath]) npady[cath]++;
      yhig[cath] = fXyq[1][j]; 
    }
  }
  //if (lun1) fprintf(lun1," %4d %2d %3d %3d %3d %3d \n",gAlice->GetHeader()->GetEvent(),AliMUONClusterInput::Instance()->Chamber(), npadx[0], npadx[1], npady[0], npady[1]);

  // Exclude pads with overflows
  for (Int_t j=0; j<npad; j++) {
    if (fXyq[2][j] > fgkSaturation-1) fPadIJ[1][j] = -5;
    else fPadIJ[1][j] = 0;
  }

  // Compute couplings of clusters to pads
  TMatrixD *aijclupad = new TMatrixD(nclust,npad);
  *aijclupad = 0;
  Int_t npxclu;
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    pix = clusters[iclust];
    npxclu = pix->GetEntriesFast();
    for (Int_t i=0; i<npxclu; i++) {
      indx = fPixArray->IndexOf(pix->UncheckedAt(i));
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[1][j] < 0 && fPadIJ[1][j] != -5) continue;
	if (coef[j*nPix+indx] < fgkCouplMin) continue;
	(*aijclupad)(iclust,j) += coef[j*nPix+indx];
      }
    }
  }
  // Compute couplings between clusters
  TMatrixD *aijcluclu = new TMatrixD(nclust,nclust);
  *aijcluclu = 0;
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    for (Int_t j=0; j<npad; j++) {
      // Exclude overflows
      if (fPadIJ[1][j] < 0) continue;
      if ((*aijclupad)(iclust,j) < fgkCouplMin) continue;
      for (Int_t iclust1=iclust+1; iclust1<nclust; iclust1++) {
	if ((*aijclupad)(iclust1,j) < fgkCouplMin) continue;
	(*aijcluclu)(iclust,iclust1) += 
	  TMath::Sqrt ((*aijclupad)(iclust,j)*(*aijclupad)(iclust1,j));
      }
    }
  }
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    for (Int_t iclust1=iclust+1; iclust1<nclust; iclust1++) {
      (*aijcluclu)(iclust1,iclust) = (*aijcluclu)(iclust,iclust1);
    }
  }

  if (fDebug && nclust > 1) aijcluclu->Print();

  // Find groups of coupled clusters
  used = new Bool_t[nclust];
  for (Int_t i=0; i<nclust; i++) used[i] = kFALSE;
  Int_t *clustNumb = new Int_t[nclust];
  Int_t nCoupled, nForFit, minGroup[3], clustFit[3], nfit = 0;
  Double_t parOk[8];

  for (Int_t igroup=0; igroup<nclust; igroup++) {
    if (used[igroup]) continue;
    used[igroup] = kTRUE;
    clustNumb[0] = igroup;
    nCoupled = 1;
    // Find group of coupled clusters
    AddCluster(igroup, nclust, aijcluclu, used, clustNumb, nCoupled); // recursive
    if (fDebug) {
      cout << " nCoupled: " << nCoupled << endl;
      for (Int_t i=0; i<nCoupled; i++) cout << clustNumb[i] << " "; cout << endl; 
    }
    fnCoupled = nCoupled;

    while (nCoupled > 0) {

      if (nCoupled < 4) {
	nForFit = nCoupled;
	for (Int_t i=0; i<nCoupled; i++) clustFit[i] = clustNumb[i];
      } else {
	// Too many coupled clusters to fit - try to decouple them
	// Find the lowest coupling of 1, 2, min(3,nLinks/2) pixels with 
	// all the others in the group 
	for (Int_t j=0; j<3; j++) minGroup[j] = -1;
	Double_t coupl = MinGroupCoupl(nCoupled, clustNumb, aijcluclu, minGroup);

	// Flag clusters for fit
	nForFit = 0;
	while (minGroup[nForFit] >= 0 && nForFit < 3) {
	  if (fDebug) cout << clustNumb[minGroup[nForFit]] << " ";
	  clustFit[nForFit] = clustNumb[minGroup[nForFit]];
	  clustNumb[minGroup[nForFit]] -= 999;
	  nForFit++;
	}
	if (fDebug) cout << nForFit << " " << coupl << endl;
      } // else

      // Select pads for fit. 
      if (SelectPad(nCoupled, nForFit, clustNumb, clustFit, aijclupad) < 3 && nCoupled > 1) {
	// Deselect pads
	for (Int_t j=0; j<npad; j++) {
	  if (TMath::Abs(fPadIJ[1][j]) == 1) fPadIJ[1][j] = 0;
	  if (TMath::Abs(fPadIJ[1][j]) == -9) fPadIJ[1][j] = -5;
	}
	// Merge the failed cluster candidates (with too few pads to fit) with 
	// the one with the strongest coupling
	Merge(nForFit, nCoupled, clustNumb, clustFit, clusters, aijcluclu, aijclupad);
      } else {
	// Do the fit
	nfit = Fit(0, nForFit, clustFit, clusters, parOk);
      }

      // Subtract the fitted charges from pads with strong coupling and/or
      // return pads for further use
      UpdatePads(nfit, parOk);

      // Mark used pads
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[1][j] == 1) fPadIJ[1][j] = -1;
	if (fPadIJ[1][j] == -9) fPadIJ[1][j] = -5;
      }

      // Sort the clusters (move to the right the used ones)
      Int_t beg = 0, end = nCoupled - 1;
      while (beg < end) {
        if (clustNumb[beg] >= 0) { beg++; continue; }
        for (Int_t j=end; j>beg; j--) {
          if (clustNumb[j] < 0) continue;
          end = j - 1;
          indx = clustNumb[beg];
          clustNumb[beg] = clustNumb[j];
          clustNumb[j] = indx;
          break;
        }
	beg++;
      }

      nCoupled -= nForFit;
      if (nCoupled > 3) {
	// Remove couplings of used clusters
	for (Int_t iclust=nCoupled; iclust<nCoupled+nForFit; iclust++) {
	  indx = clustNumb[iclust] + 999;
	  for (Int_t iclust1=0; iclust1<nCoupled; iclust1++) {
	    indx1 = clustNumb[iclust1];
	    (*aijcluclu)(indx,indx1) = (*aijcluclu)(indx1,indx) = 0;
	  }
	}

	// Update the remaining clusters couplings (exclude couplings from 
	// the used pads)
	for (Int_t j=0; j<npad; j++) {
	  if (fPadIJ[1][j] != -1) continue;
	  for (Int_t iclust=0; iclust<nCoupled; iclust++) {
	    indx = clustNumb[iclust];
	    if ((*aijclupad)(indx,j) < fgkCouplMin) continue;
	    for (Int_t iclust1=iclust+1; iclust1<nCoupled; iclust1++) {
	      indx1 = clustNumb[iclust1];
	      if ((*aijclupad)(indx1,j) < fgkCouplMin) continue;
	      // Check this
	      (*aijcluclu)(indx,indx1) -= 
		TMath::Sqrt ((*aijclupad)(indx,j)*(*aijclupad)(indx1,j));
	      (*aijcluclu)(indx1,indx) = (*aijcluclu)(indx,indx1);
	    }
	  }
	  fPadIJ[1][j] = -8;
	} // for (Int_t j=0; j<npad;
      } // if (nCoupled > 3)
    } // while (nCoupled > 0)
  } // for (Int_t igroup=0; igroup<nclust;

  aijcluclu->Delete(); aijclupad->Delete();
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    pix = clusters[iclust]; 
    pix->Clear();
    delete pix; pix = 0;
  }
  delete [] clustNumb; clustNumb = 0; delete [] used; used = 0;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddBin(TH2D *mlem, Int_t ic, Int_t jc, Int_t mode, Bool_t *used, TObjArray *pix)
{
/// Add a bin to the cluster

  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t cont1, cont = mlem->GetCellContent(jc,ic);
  AliMUONPixel *pixPtr = 0;

  for (Int_t i=TMath::Max(ic-1,1); i<=TMath::Min(ic+1,ny); i++) {
    for (Int_t j=TMath::Max(jc-1,1); j<=TMath::Min(jc+1,nx); j++) {
      if (i != ic && j != jc) continue;
      if (used[(i-1)*nx+j-1]) continue;
      cont1 = mlem->GetCellContent(j,i);
      if (mode && cont1 > cont) continue;
      used[(i-1)*nx+j-1] = kTRUE;
      if (cont1 < 0.5) continue;
      if (pix) pix->Add(BinToPix(mlem,j,i)); 
      else {
	pixPtr = new AliMUONPixel (mlem->GetXaxis()->GetBinCenter(j), 
				   mlem->GetYaxis()->GetBinCenter(i), 0, 0, cont1);
	fPixArray->Add((TObject*)pixPtr);
      }
      AddBin(mlem, i, j, mode, used, pix); // recursive call
    }
  }
}

//_____________________________________________________________________________
TObject* AliMUONClusterFinderAZ::BinToPix(TH2D *mlem, Int_t jc, Int_t ic)
{
/// Translate histogram bin to pixel 
  
  Double_t yc = mlem->GetYaxis()->GetBinCenter(ic);
  Double_t xc = mlem->GetXaxis()->GetBinCenter(jc);
  
  Int_t nPix = fPixArray->GetEntriesFast();
  AliMUONPixel *pixPtr = NULL;

  // Compare pixel and bin positions
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 0.5) continue;
    if (TMath::Abs(pixPtr->Coord(0)-xc)<1.e-4 && TMath::Abs(pixPtr->Coord(1)-yc)<1.e-4) return (TObject*) pixPtr;
  }
  AliError(Form(" Something wrong ??? %f %f ", xc, yc));
  return NULL;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddCluster(Int_t ic, Int_t nclust, TMatrixD *aijcluclu, Bool_t *used, Int_t *clustNumb, Int_t &nCoupled)
{
/// Add a cluster to the group of coupled clusters

  for (Int_t i=0; i<nclust; i++) {
    if (used[i]) continue;
    if ((*aijcluclu)(i,ic) < fgkCouplMin) continue;
    used[i] = kTRUE;
    clustNumb[nCoupled++] = i;
    AddCluster(i, nclust, aijcluclu, used, clustNumb, nCoupled);
  }
}

//_____________________________________________________________________________
Double_t AliMUONClusterFinderAZ::MinGroupCoupl(Int_t nCoupled, Int_t *clustNumb, TMatrixD *aijcluclu, Int_t *minGroup)
{
/// Find group of clusters with minimum coupling to all the others

  Int_t i123max = TMath::Min(3,nCoupled/2); 
  Int_t indx, indx1, indx2, indx3, nTot = 0;
  Double_t *coupl1 = 0, *coupl2 = 0, *coupl3 = 0;

  for (Int_t i123=1; i123<=i123max; i123++) {

    if (i123 == 1) {
      coupl1 = new Double_t [nCoupled];
      for (Int_t i=0; i<nCoupled; i++) coupl1[i] = 0;
    }
    else if (i123 == 2) {
      nTot = nCoupled*nCoupled;
      coupl2 = new Double_t [nTot];
      for (Int_t i=0; i<nTot; i++) coupl2[i] = 9999;
    } else {
      nTot = nTot*nCoupled;
      coupl3 = new Double_t [nTot];
      for (Int_t i=0; i<nTot; i++) coupl3[i] = 9999;
    } // else

    for (Int_t i=0; i<nCoupled; i++) {
      indx1 = clustNumb[i];
      for (Int_t j=i+1; j<nCoupled; j++) {
	indx2 = clustNumb[j];
	if (i123 == 1) {
	  coupl1[i] += (*aijcluclu)(indx1,indx2);
	  coupl1[j] += (*aijcluclu)(indx1,indx2);
	} 
	else if (i123 == 2) {
	  indx = i*nCoupled + j;
	  coupl2[indx] = coupl1[i] + coupl1[j];
	  coupl2[indx] -= 2 * ((*aijcluclu)(indx1,indx2));
	} else {
	  for (Int_t k=j+1; k<nCoupled; k++) {
	    indx3 = clustNumb[k];
	    indx = i*nCoupled*nCoupled + j*nCoupled + k;
	    coupl3[indx] = coupl2[i*nCoupled+j] + coupl1[k];
	    coupl3[indx] -= 2 * ((*aijcluclu)(indx1,indx3)+(*aijcluclu)(indx2,indx3));
	  }
	} // else
      } // for (Int_t j=i+1;
    } // for (Int_t i=0;
  } // for (Int_t i123=1;

  // Find minimum coupling
  Double_t couplMin = 9999;
  Int_t locMin = 0;

  for (Int_t i123=1; i123<=i123max; i123++) {
    if (i123 == 1) {
      locMin = TMath::LocMin(nCoupled, coupl1);
      couplMin = coupl1[locMin];
      minGroup[0] = locMin;
      delete [] coupl1; coupl1 = 0;
    } 
    else if (i123 == 2) {
      locMin = TMath::LocMin(nCoupled*nCoupled, coupl2);
      if (coupl2[locMin] < couplMin) {
	couplMin = coupl2[locMin];
	minGroup[0] = locMin/nCoupled;
	minGroup[1] = locMin%nCoupled;
      }
      delete [] coupl2; coupl2 = 0;
    } else {
      locMin = TMath::LocMin(nTot, coupl3);
      if (coupl3[locMin] < couplMin) {
	couplMin = coupl3[locMin];
	minGroup[0] = locMin/nCoupled/nCoupled;
	minGroup[1] = locMin%(nCoupled*nCoupled)/nCoupled;
	minGroup[2] = locMin%nCoupled;
      }
      delete [] coupl3; coupl3 = 0;
    } // else
  } // for (Int_t i123=1;
  return couplMin;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::SelectPad(Int_t nCoupled, Int_t nForFit, Int_t *clustNumb, Int_t *clustFit, TMatrixD *aijclupad)
{
/// Select pads for fit. If too many coupled clusters, find pads giving 
/// the strongest coupling with the rest of clusters and exclude them from the fit.

  Int_t npad = fnPads[0] + fnPads[1];
  Double_t *padpix = 0;

  if (nCoupled > 3) {
    padpix = new Double_t[npad];
    for (Int_t i=0; i<npad; i++) padpix[i] = 0; 
  }

  Int_t nOK = 0, indx, indx1;
  for (Int_t iclust=0; iclust<nForFit; iclust++) {
    indx = clustFit[iclust];
    for (Int_t j=0; j<npad; j++) {
      if ((*aijclupad)(indx,j) < fgkCouplMin) continue;
      if (fPadIJ[1][j] == -5) fPadIJ[1][j] = -9; // flag overflow
      if (fPadIJ[1][j] < 0) continue; // exclude overflows and used pads
      if (!fPadIJ[1][j]) { fPadIJ[1][j] = 1; nOK++; } // pad to be used in fit
      if (nCoupled > 3) {
	// Check other clusters
	for (Int_t iclust1=0; iclust1<nCoupled; iclust1++) {
	  indx1 = clustNumb[iclust1];
	  if (indx1 < 0) continue;
	  if ((*aijclupad)(indx1,j) < fgkCouplMin) continue;
	  padpix[j] += (*aijclupad)(indx1,j);
	}
      } // if (nCoupled > 3)
    } // for (Int_t j=0; j<npad;
  } // for (Int_t iclust=0; iclust<nForFit
  if (nCoupled < 4) return nOK;

  Double_t aaa = 0;
  for (Int_t j=0; j<npad; j++) {
    if (padpix[j] < fgkCouplMin) continue;
    if (fDebug) cout << j << " " << padpix[j] << " " << fXyq[0][j] << " " << fXyq[1][j] << endl;
    aaa += padpix[j];
    fPadIJ[1][j] = -1; // exclude pads with strong coupling to the other clusters
    nOK--;
  }
  delete [] padpix; padpix = 0;
  return nOK;
}
  
//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Merge(Int_t nForFit, Int_t nCoupled, Int_t *clustNumb, Int_t *clustFit, TObjArray **clusters, TMatrixD *aijcluclu, TMatrixD *aijclupad)
{
/// Merge the group of clusters with the one having the strongest coupling with them

  Int_t indx, indx1, npxclu, npxclu1, imax=0;
  TObjArray *pix, *pix1;
  Double_t couplMax;

  for (Int_t icl=0; icl<nForFit; icl++) {
    indx = clustFit[icl];
    pix = clusters[indx];
    npxclu = pix->GetEntriesFast();
    couplMax = -1;
    for (Int_t icl1=0; icl1<nCoupled; icl1++) {
      indx1 = clustNumb[icl1];
      if (indx1 < 0) continue;
      if ((*aijcluclu)(indx,indx1) > couplMax) {
	couplMax = (*aijcluclu)(indx,indx1);
	imax = indx1;
      }
    } // for (Int_t icl1=0;
    /*if (couplMax < fgkCouplMin) {
      cout << " Oops " << couplMax << endl;
      aijcluclu->Print();
      cout << icl << " " << indx << " " << npxclu << " " << nLinks << endl;
      ::exit(0);
      }*/
    // Add to it
    pix1 = clusters[imax];
    npxclu1 = pix1->GetEntriesFast();
    // Add pixels 
    for (Int_t i=0; i<npxclu; i++) { pix1->Add(pix->UncheckedAt(i)); pix->RemoveAt(i); }
    if (fDebug) cout << " New number of pixels: " << npxclu1 << " " << pix1->GetEntriesFast() << endl;
    //Add cluster-to-cluster couplings
    //aijcluclu->Print();
    for (Int_t icl1=0; icl1<nCoupled; icl1++) {
      indx1 = clustNumb[icl1];
      if (indx1 < 0 || indx1 == imax) continue;
      (*aijcluclu)(indx1,imax) += (*aijcluclu)(indx,indx1);
      (*aijcluclu)(imax,indx1) = (*aijcluclu)(indx1,imax);
    }
    (*aijcluclu)(indx,imax) = (*aijcluclu)(imax,indx) = 0;
    //aijcluclu->Print();
    //Add cluster-to-pad couplings
    for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
      if (fPadIJ[1][j] < 0 && fPadIJ[1][j] != -5) continue; // exclude used pads
      (*aijclupad)(imax,j) += (*aijclupad)(indx,j);
      (*aijclupad)(indx,j) = 0;
    }
  } // for (Int_t icl=0; icl<nForFit;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::Fit(Int_t iSimple, Int_t nfit, Int_t *clustFit, TObjArray **clusters, Double_t *parOk)
{
/// Find selected clusters to selected pad charges
  
  TH2D *mlem = (TH2D*) gROOT->FindObject("mlem");
  Double_t xmin = mlem->GetXaxis()->GetXmin() - mlem->GetXaxis()->GetBinWidth(1);
  Double_t xmax = mlem->GetXaxis()->GetXmax() + mlem->GetXaxis()->GetBinWidth(1);
  Double_t ymin = mlem->GetYaxis()->GetXmin() - mlem->GetYaxis()->GetBinWidth(1);
  Double_t ymax = mlem->GetYaxis()->GetXmax() + mlem->GetYaxis()->GetBinWidth(1);
  Double_t step[3]={0.01,0.002,0.02}, xPad = 0, yPad = 99999;

  // Number of pads to use and number of virtual pads
  Int_t npads = 0, nVirtual = 0, nfit0 = nfit;
  for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
    if (fXyq[3][i] < 0) nVirtual++;
    if (fPadIJ[1][i] != 1) continue;
    if (fXyq[3][i] > 0) {
      npads++;
      if (yPad > 9999) { 
	xPad = fXyq[0][i]; 
	yPad = fXyq[1][i]; 
      } else {
	if (fXyq[4][i] < fXyq[3][i]) yPad = fXyq[1][i]; 
	else xPad = fXyq[0][i]; 
      }
    }
  }
  if (fDebug) {
    for (Int_t i=0; i<nfit; i++) {cout << i+1 << " " << clustFit[i] << " ";}
    cout << nfit << endl;
    cout << " Number of pads to fit: " << npads << endl;
  }
  fNpar = 0;
  fQtot = 0;
  if (npads < 2) return 0; 
  
  Int_t digit = 0;
  AliMUONDigit *mdig = 0;
  Int_t tracks[3] = {-1, -1, -1};
  for (Int_t cath=0; cath<2; cath++) {  
    for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
      if (fPadIJ[0][i] != cath) continue;
      if (fPadIJ[1][i] != 1) continue;
      if (fXyq[3][i] < 0) continue; // exclude virtual pads
      digit = TMath::Nint (fXyq[5][i]);
      if (digit >= 0) mdig = fInput->Digit(cath,digit);
      else mdig = fInput->Digit(TMath::Even(cath),-digit-1);
      //if (!mdig) mdig = fInput->Digit(TMath::Even(cath),digit);
      if (!mdig) continue; // protection for cluster display
      if (mdig->Hit() >= 0) {
	if (tracks[0] < 0) {
	  tracks[0] = mdig->Hit();
	  tracks[1] = mdig->Track(0);
	} else if (mdig->Track(0) < tracks[1]) {
	  tracks[0] = mdig->Hit();
	  tracks[1] = mdig->Track(0);
	}
      }
      if (mdig->Track(1) >= 0 && mdig->Track(1) != tracks[1]) {
	if (tracks[2] < 0) tracks[2] = mdig->Track(1);
	else tracks[2] = TMath::Min (tracks[2], mdig->Track(1));
      }
      //if (!mdig) break;
      //cout << mdig->Hit() << " " << mdig->Track(0) << " " << mdig->Track(1) <<endl;
    } // for (Int_t i=0;
  } // for (Int_t cath=0;
  //cout << tracks[0] << " " << tracks[1] << " " << tracks[2] <<endl;
  
  // Get number of pads in X and Y 
  Int_t nInX = 0, nInY;
  PadsInXandY(nInX, nInY);
  //cout << " nInX and Y: " << nInX << " " << nInY << endl;

  Int_t nfitMax = 3; 
  nfitMax = TMath::Min (nfitMax, (npads + 1) / 3);
  if (nfitMax > 1) {
    if (nInX < 3 && nInY < 3 || nInX == 3 && nInY < 3 || nInX < 3 && nInY == 3) nfitMax = 1; // not enough pads in each direction
  }
  if (nfit > nfitMax) nfit = nfitMax;

  // Take cluster maxima as fitting seeds
  TObjArray *pix;
  AliMUONPixel *pixPtr;
  Int_t npxclu;
  Double_t cont, cmax = 0, xseed = 0, yseed = 0, errOk[8], qq = 0;
  Double_t xyseed[3][2], qseed[3], xyCand[3][2] = {{0},{0}}, sigCand[3][2] = {{0},{0}};

  for (Int_t ifit=1; ifit<=nfit0; ifit++) {
    cmax = 0;
    pix = clusters[clustFit[ifit-1]];
    npxclu = pix->GetEntriesFast();
    //qq = 0;
    for (Int_t clu=0; clu<npxclu; clu++) {
      pixPtr = (AliMUONPixel*) pix->UncheckedAt(clu);
      cont = pixPtr->Charge();
      fQtot += cont;
      if (cont > cmax) { 
	cmax = cont; 
	xseed = pixPtr->Coord(0);
	yseed = pixPtr->Coord(1);
      }
      qq += cont;
      /*
      xyCand[ifit-1][0] += pixPtr->Coord(0) * cont;
      xyCand[ifit-1][1] += pixPtr->Coord(1) * cont;
      sigCand[ifit-1][0] += pixPtr->Coord(0) * pixPtr->Coord(0) * cont;
      sigCand[ifit-1][1] += pixPtr->Coord(1) * pixPtr->Coord(1) * cont;
      */
      xyCand[0][0] += pixPtr->Coord(0) * cont;
      xyCand[0][1] += pixPtr->Coord(1) * cont;
      sigCand[0][0] += pixPtr->Coord(0) * pixPtr->Coord(0) * cont;
      sigCand[0][1] += pixPtr->Coord(1) * pixPtr->Coord(1) * cont;
    }
    xyseed[ifit-1][0] = xseed;
    xyseed[ifit-1][1] = yseed;
    qseed[ifit-1] = cmax;
    /*
    xyCand[ifit-1][0] /= qq; // <x>
    xyCand[ifit-1][1] /= qq; // <y>
    sigCand[ifit-1][0] = sigCand[ifit-1][0]/qq - xyCand[ifit-1][0]*xyCand[ifit-1][0]; // <x^2> - <x>^2
    sigCand[ifit-1][0] = sigCand[ifit-1][0] > 0 ? TMath::Sqrt (sigCand[ifit-1][0]) : 0;
    sigCand[ifit-1][1] = sigCand[ifit-1][1]/qq - xyCand[ifit-1][1]*xyCand[ifit-1][1]; // <y^2> - <y>^2
    sigCand[ifit-1][1] = sigCand[ifit-1][1] > 0 ? TMath::Sqrt (sigCand[ifit-1][1]) : 0;
    cout << xyCand[ifit-1][0] << " " << xyCand[ifit-1][1] << " " << sigCand[ifit-1][0] << " " << sigCand[ifit-1][1] << endl;
    */
  } // for (Int_t ifit=1;

  xyCand[0][0] /= qq; // <x>
  xyCand[0][1] /= qq; // <y>
  sigCand[0][0] = sigCand[0][0]/qq - xyCand[0][0]*xyCand[0][0]; // <x^2> - <x>^2
  sigCand[0][0] = sigCand[0][0] > 0 ? TMath::Sqrt (sigCand[0][0]) : 0;
  sigCand[0][1] = sigCand[0][1]/qq - xyCand[0][1]*xyCand[0][1]; // <y^2> - <y>^2
  sigCand[0][1] = sigCand[0][1] > 0 ? TMath::Sqrt (sigCand[0][1]) : 0;
  if (fDebug) cout << xyCand[0][0] << " " << xyCand[0][1] << " " << sigCand[0][0] << " " << sigCand[0][1] << endl;

  Int_t nDof, maxSeed[3], nMax = 0;
  Double_t fmin, chi2o = 9999, chi2n;

  TMath::Sort(nfit0, qseed, maxSeed, kTRUE); // in decreasing order
  /*
  Int_t itmp[100], localMax[100];
  Double_t maxVal[100];
  if (!iSimple && nfit < nfitMax) {
    // Try to split pixel cluster according to local maxima
    Int_t nfit1 = nfit;
    for (Int_t iclus = 0; iclus < nfit1; iclus++) {
      nMax = FindLocalMaxima (clusters[clustFit[maxSeed[iclus]]], localMax, maxVal);
      TH2D *hist = (TH2D*) gROOT->FindObject("anode1");
      if (nMax == 1) { hist->Delete(); continue; }
      // Add extra fitting seeds from local maxima
      Int_t ixseed = hist->GetXaxis()->FindBin(xyseed[maxSeed[iclus]][0]);
      Int_t iyseed = hist->GetYaxis()->FindBin(xyseed[maxSeed[iclus]][1]);
      Int_t nx = hist->GetNbinsX();
      TMath::Sort(nMax, maxVal, itmp, kTRUE); // in decreasing order
      for (Int_t j = 0; j < nMax; j++) {
	Int_t iyc = localMax[itmp[j]] / nx + 1;
	Int_t ixc = localMax[itmp[j]] % nx + 1;
	if (ixc == ixseed && iyc == iyseed) continue; // local max already taken for seeding
	xyseed[nfit][0] = hist->GetXaxis()->GetBinCenter(ixc);
	xyseed[nfit][1] = hist->GetYaxis()->GetBinCenter(iyc);
	qseed[nfit] = maxVal[itmp[j]];
	maxSeed[nfit] = nfit++;
	if (nfit >= nfitMax) break;
      }
      hist->Delete();
      if (nfit >= nfitMax) break;
    } // for (Int_t iclus = 0;
    //nfit0 = nfit;
    //TMath::Sort(nfit0, qseed, maxSeed, kTRUE); // in decreasing order
  } //if (!iSimple && nfit < nfitMax)
  */

  Double_t *gin = 0, func0, func1, param[8], step0[8];
  Double_t param0[2][8]={{0},{0}}, deriv[2][8]={{0},{0}}; 
  Double_t shift[8], stepMax, derMax, parmin[8], parmax[8], func2[2], shift0;
  Double_t delta[8], scMax, dder[8], estim, shiftSave = 0;
  Int_t min, max, nCall = 0, memory[8] = {0}, nLoop, idMax = 0, iestMax = 0, nFail;
  Double_t rad, dist[3] = {0};

  // Try to fit with one-track hypothesis, then 2-track. If chi2/dof is 
  // lower, try 3-track (if number of pads is sufficient).
  for (Int_t iseed=0; iseed<nfit; iseed++) {

    if (iseed) { for (Int_t j=0; j<fNpar; j++) param[j] = parOk[j]; } // for bounded params
    for (Int_t j=0; j<3; j++) step0[fNpar+j] = shift[fNpar+j] = step[j];
    if (nfit == 1) param[fNpar] = xyCand[0][0]; // take COG
    else param[fNpar] = xyseed[maxSeed[iseed]][0];
    parmin[fNpar] = xmin; 
    parmax[fNpar++] = xmax; 
    if (nfit == 1) param[fNpar] = xyCand[0][1]; // take COG
    else param[fNpar] = xyseed[maxSeed[iseed]][1];
    parmin[fNpar] = ymin; 
    parmax[fNpar++] = ymax; 
    if (fNpar > 2) {
      param[fNpar] = fNpar == 4 ? 0.5 : 0.3;
      parmin[fNpar] = 0; 
      parmax[fNpar++] = 1; 
    }
    if (iseed) { for (Int_t j=0; j<fNpar; j++) param0[1][j] = 0; }

    // Try new algorithm
    min = nLoop = 1; stepMax = func2[1] = derMax = 999999; nFail = 0;

    while (1) {
      max = !min;
      Fcn1(fNpar, gin, func0, param, 1); nCall++;
      //cout << " Func: " << func0 << endl;

      func2[max] = func0;
      for (Int_t j=0; j<fNpar; j++) {
	param0[max][j] = param[j];
	delta[j] = step0[j];
	param[j] += delta[j] / 10;
	if (j > 0) param[j-1] -= delta[j-1] / 10;
	Fcn1(fNpar, gin, func1, param, 1); nCall++;
	deriv[max][j] = (func1 - func0) / delta[j] * 10; // first derivative
	//cout << j << " " << deriv[max][j] << endl;
	dder[j] = param0[0][j] != param0[1][j] ? (deriv[0][j] - deriv[1][j]) / 
	  (param0[0][j] - param0[1][j]) : 0; // second derivative
      }
      param[fNpar-1] -= delta[fNpar-1] / 10;
      if (nCall > 2000) break;

      min = func2[0] < func2[1] ? 0 : 1;
      nFail = min == max ? 0 : nFail + 1;

      stepMax = derMax = estim = 0;
      for (Int_t j=0; j<fNpar; j++) { 
	// Estimated distance to minimum
	shift0 = shift[j];
	if (nLoop == 1) shift[j] = TMath::Sign (step0[j], -deriv[max][j]); // first step
	else if (TMath::Abs(deriv[0][j]) < 1.e-3 && TMath::Abs(deriv[1][j]) < 1.e-3) shift[j] = 0;
	else if (deriv[min][j]*deriv[!min][j] > 0 && TMath::Abs(deriv[min][j]) > TMath::Abs(deriv[!min][j])
		 //|| TMath::Abs(deriv[0][j]-deriv[1][j]) < 1.e-3) {
	  || TMath::Abs(deriv[0][j]-deriv[1][j]) < 1.e-3 || TMath::Abs(dder[j]) < 1.e-6) {
	  shift[j] = -TMath::Sign (shift[j], (func2[0]-func2[1]) * (param0[0][j]-param0[1][j]));
	  if (min == max) { 
	    if (memory[j] > 1) { shift[j] *= 2; } //cout << " Memory " << memory[j] << " " << shift[j] << endl; }
	    memory[j]++;
	  }
	} else {
	  shift[j] = dder[j] != 0 ? -deriv[min][j] / dder[j] : 0;
	  memory[j] = 0;
	}
	if (TMath::Abs(shift[j])/step0[j] > estim) { 
	  estim = TMath::Abs(shift[j])/step0[j];
	  iestMax = j;
	}

	// Too big step
	if (TMath::Abs(shift[j])/step0[j] > 10) shift[j] = TMath::Sign(10.,shift[j]) * step0[j]; // 

	// Failed to improve minimum
	if (min != max) {
	  memory[j] = 0;
	  param[j] = param0[min][j];
	  if (TMath::Abs(shift[j]+shift0) > 0.1*step0[j]) shift[j] = (shift[j] + shift0) / 2;
	  else shift[j] /= -2;
	} 

	// Too big step
	if (TMath::Abs(shift[j]*deriv[min][j]) > func2[min]) 
	  shift[j] = TMath::Sign (func2[min]/deriv[min][j], shift[j]);

	// Introduce step relaxation factor
	if (memory[j] < 3) {
	  scMax = 1 + 4 / TMath::Max(nLoop/2.,1.);
	  if (TMath::Abs(shift0) > 0 && TMath::Abs(shift[j]/shift0) > scMax) 
	    shift[j] = TMath::Sign (shift0*scMax, shift[j]);
	}
	param[j] += shift[j]; 
	//AZ Check parameter limits 27-12-2004
	if (param[j] < parmin[j]) { 
	  shift[j] = parmin[j] - param[j]; 
	  param[j] = parmin[j]; 
	} else if (param[j] > parmax[j]) {
          shift[j] = parmax[j] - param[j];
          param[j] = parmax[j];
	}
	//cout << " xxx " << j << " " << shift[j] << " " << param[j] << endl;
	stepMax = TMath::Max (stepMax, TMath::Abs(shift[j]/step0[j]));
	if (TMath::Abs(deriv[min][j]) > derMax) {
	  idMax = j;
	  derMax = TMath::Abs (deriv[min][j]);
	}
      } // for (Int_t j=0; j<fNpar;
      //cout << max << " " << func2[min] << " " << derMax << " " << stepMax << " " << estim << " " << iestMax << " " << nCall << endl;
      if (estim < 1 && derMax < 2 || nLoop > 150) break; // minimum was found

      nLoop++;
      // Check for small step
      if (shift[idMax] == 0) { shift[idMax] = step0[idMax]/10; param[idMax] += shift[idMax]; continue; }
      if (!memory[idMax] && derMax > 0.5 && nLoop > 10) {
	//cout << " ok " << deriv[min][idMax] << " " << deriv[!min][idMax] << " " << dder[idMax]*shift[idMax] << " " << shift[idMax] << endl;
	if (dder[idMax] != 0 && TMath::Abs(deriv[min][idMax]/dder[idMax]/shift[idMax]) > 10) {
	  if (min == max) dder[idMax] = -dder[idMax];
	  shift[idMax] = -deriv[min][idMax] / dder[idMax] / 10; 
	  param[idMax] += shift[idMax];
	  stepMax = TMath::Max (stepMax, TMath::Abs(shift[idMax])/step0[idMax]);
	  //cout << shift[idMax] << " " << param[idMax] << endl;
	  if (min == max) shiftSave = shift[idMax];
	}
	if (nFail > 10) {
	  param[idMax] -= shift[idMax];
	  shift[idMax] = 4 * shiftSave * (gRandom->Rndm(0) - 0.5);
	  param[idMax] += shift[idMax];
	  //cout << shift[idMax] << endl;
	}
      }      
    } // while (1)
    fmin = func2[min];

    nDof = npads - fNpar + nVirtual;
    if (!nDof) nDof++;
    chi2n = fmin / nDof;
    if (fDebug) cout << " Chi2 " << chi2n << " " << fNpar << endl;

    if (chi2n*1.2+1.e-6 > chi2o ) { fNpar -= 3; break; }

    // Save parameters and errors

    if (nInX == 1) {
      // One pad per direction 
      for (Int_t i=0; i<fNpar; i++) if (i == 0 || i == 2 || i == 5) param0[min][i] = xPad;
    }
    if (nInY == 1) {
      // One pad per direction 
      for (Int_t i=0; i<fNpar; i++) if (i == 1 || i == 3 || i == 6) param0[min][i] = yPad;
    }

    /*
    if (iseed > 0) {
      // Find distance to the nearest neighbour
      dist[0] = dist[1] = TMath::Sqrt ((param0[min][0]-param0[min][2])*
				       (param0[min][0]-param0[min][2])
				      +(param0[min][1]-param0[min][3])*
				       (param0[min][1]-param0[min][3]));
      if (iseed > 1) {
	dist[2] = TMath::Sqrt ((param0[min][0]-param0[min][5])*
			       (param0[min][0]-param0[min][5])
			      +(param0[min][1]-param0[min][6])*
			       (param0[min][1]-param0[min][6]));
	rad = TMath::Sqrt ((param0[min][2]-param0[min][5])*
			   (param0[min][2]-param0[min][5])
			  +(param0[min][3]-param0[min][6])*
			   (param0[min][3]-param0[min][6]));
	if (dist[2] < dist[0]) dist[0] = dist[2];
	if (rad < dist[1]) dist[1] = rad;
	if (rad < dist[2]) dist[2] = rad;
      }
      cout << dist[0] << " " << dist[1] << " " << dist[2] << endl;
      if (dist[TMath::LocMin(iseed+1,dist)] < 1.) { fNpar -= 3; break; }
    }
    */

    for (Int_t i=0; i<fNpar; i++) {
      parOk[i] = param0[min][i];
      //errOk[i] = fmin;
      errOk[i] = chi2n;
      // Bounded params
      parOk[i] = TMath::Max (parOk[i], parmin[i]);
      parOk[i] = TMath::Min (parOk[i], parmax[i]);
    }

    chi2o = chi2n;
    if (fmin < 0.1) break; // !!!???
  } // for (Int_t iseed=0; 

  if (fDebug) {
    for (Int_t i=0; i<fNpar; i++) {
      if (i == 4 || i == 7) {
	if (i == 7 || i == 4 && fNpar < 7) cout << parOk[i] << endl;
	else cout << parOk[i] * (1-parOk[7]) << endl;
	continue;
      }
      cout << parOk[i] << " " << errOk[i] << endl;
    }
  }
  nfit = (fNpar + 1) / 3;
  dist[0] = dist[1] = dist[2] = 0;

  if (nfit > 1) {
    // Find distance to the nearest neighbour
    dist[0] = dist[1] = TMath::Sqrt ((parOk[0]-parOk[2])*
				     (parOk[0]-parOk[2])
				    +(parOk[1]-parOk[3])*
				     (parOk[1]-parOk[3]));
    if (nfit > 2) {
      dist[2] = TMath::Sqrt ((parOk[0]-parOk[5])*
			     (parOk[0]-parOk[5])
			    +(parOk[1]-parOk[6])*
			     (parOk[1]-parOk[6]));
      rad = TMath::Sqrt ((parOk[2]-parOk[5])*
			 (parOk[2]-parOk[5])
			+(parOk[3]-parOk[6])*
			 (parOk[3]-parOk[6]));
      if (dist[2] < dist[0]) dist[0] = dist[2];
      if (rad < dist[1]) dist[1] = rad;
      if (rad < dist[2]) dist[2] = rad;
    }
  }

  Int_t indx;
  fnPads[1] -= nVirtual;
  if (!fDraw) {
    Double_t coef = 0;
    if (iSimple) fnCoupled = 0;
    //for (Int_t j=0; j<nfit; j++) {
    for (Int_t j=nfit-1; j>=0; j--) {
      indx = j<2 ? j*2 : j*2+1;  
      if (nfit == 1) coef = 1;
      else coef = j==nfit-1 ? parOk[indx+2] : 1-coef;
      coef = TMath::Max (coef, 0.);
      if (nfit == 3 && j < 2) coef = j==1 ? coef*parOk[indx+2] : coef - parOk[7];
      coef = TMath::Max (coef, 0.);
      AddRawCluster (parOk[indx], parOk[indx+1], coef*fQtot, errOk[indx], nfit0+10*nfit+100*nMax+10000*fnCoupled, tracks,
		     //sigCand[maxSeed[j]][0], sigCand[maxSeed[j]][1]);
		     //sigCand[0][0], sigCand[0][1], dist[j]);
		     sigCand[0][0], sigCand[0][1], dist[TMath::LocMin(nfit,dist)]);
    }
  } else fDraw->FillMuon(nfit, parOk, errOk);
  return nfit;
}  

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Fcn1(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
/// Fit for one track
/// AZ for Muinuit AliMUONClusterFinderAZ& c = *(AliMUONClusterFinderAZ::fgClusterFinder);    

  AliMUONClusterFinderAZ& c = *this; //AZ
  
  Int_t cath, ix, iy, indx, npads=0;
  Double_t charge, delta, coef=0, chi2=0, qTot = 0;
  for (Int_t j=0; j<c.fnPads[0]+c.fnPads[1]; j++) {
    if (c.fPadIJ[1][j] != 1) continue;
    cath = c.fPadIJ[0][j];
    if (c.fXyq[3][j] > 0) npads++; // exclude virtual pads
    qTot += c.fXyq[2][j];
    ix = c.fPadIJ[2][j];
    iy = c.fPadIJ[3][j];
    c.fSegmentation[cath]->SetPad(ix, iy);
    charge = 0;
    for (Int_t i=c.fNpar/3; i>=0; i--) { // sum over tracks
      indx = i<2 ? 2*i : 2*i+1;
      c.fSegmentation[cath]->SetHit(par[indx], par[indx+1], c.fZpad);
      if (c.fNpar == 2) coef = 1;
      else coef = i==c.fNpar/3 ? par[indx+2] : 1-coef;
      coef = TMath::Max (coef, 0.);
      if (c.fNpar == 8 && i < 2) coef = i==1 ? coef*par[indx+2] : coef - par[7];
      coef = TMath::Max (coef, 0.);
      charge += c.fInput->Mathieson()->IntXY(fInput->DetElemId(), c.fInput->Segmentation2(cath))*coef;
    }
    charge *= c.fQtot;
    delta = charge - c.fXyq[2][j];
    delta *= delta;
    delta /= c.fXyq[2][j];
    //if (cath) delta /= 5; // just for test
    chi2 += delta;
  } // for (Int_t j=0;
  f = chi2; 
  Double_t qAver = qTot/npads; //(c.fnPads[0]+c.fnPads[1]);
  f = chi2/qAver;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::UpdatePads(Int_t /*nfit*/, Double_t *par)
{
/// Subtract the fitted charges from pads with strong coupling

  Int_t cath, ix, iy, indx;
  Double_t charge, coef=0;
  for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
    if (fPadIJ[1][j] != -1) continue;
    if (fNpar != 0) {
      cath = fPadIJ[0][j];
      ix = fPadIJ[2][j];
      iy = fPadIJ[3][j];
      fSegmentation[cath]->SetPad(ix, iy);
      charge = 0;
      for (Int_t i=fNpar/3; i>=0; i--) { // sum over tracks
	indx = i<2 ? 2*i : 2*i+1;
	fSegmentation[cath]->SetHit(par[indx], par[indx+1], fZpad);
	if (fNpar == 2) coef = 1;
	else coef = i==fNpar/3 ? par[indx+2] : 1-coef;
	coef = TMath::Max (coef, 0.);
	if (fNpar == 8 && i < 2) coef = i==1 ? coef*par[indx+2] : coef - par[7];
	coef = TMath::Max (coef, 0.);
	charge += fInput->Mathieson()->IntXY(fInput->DetElemId(),fInput->Segmentation2(cath))*coef;
      }
      charge *= fQtot;
      fXyq[2][j] -= charge;
    } // if (fNpar != 0)
    if (fXyq[2][j] > fgkZeroSuppression) fPadIJ[1][j] = 0; // return pad for further using
  } // for (Int_t j=0;
}  

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::TestTrack(Int_t /*t*/) const 
{
/// Test if track was user selected

  return kTRUE;
  /*
    if (fTrack[0]==-1 || fTrack[1]==-1) {
	return kTRUE;
    } else if (t==fTrack[0] || t==fTrack[1]) {
	return kTRUE;
    } else {
	return kFALSE;
    }
  */
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddRawCluster(Double_t x, Double_t y, Double_t qTot, Double_t fmin, Int_t nfit, Int_t *tracks, Double_t /*sigx*/, Double_t /*sigy*/, Double_t /*dist*/)
{
/// Add a raw cluster copy to the list

  if (qTot <= 0.501) return; 
  AliMUONRawCluster cnew;

  Int_t cath, npads[2] = {0}, nover[2] = {0};
  for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
    cath = fPadIJ[0][j];
    // There was an overflow
    if (fPadIJ[1][j] == -9) nover[cath]++;
    if (fPadIJ[1][j] != 1 && fPadIJ[1][j] != -9) continue;
    cnew.SetMultiplicity(cath,cnew.GetMultiplicity(cath)+1);
    if (fXyq[2][j] > cnew.GetPeakSignal(cath)) cnew.SetPeakSignal(cath,fXyq[2][j]);
    //cnew.SetCharge(cath,cnew.GetCharge(cath) + TMath::Nint (fXyq[2][j]));
    cnew.SetContrib(npads[cath],cath,fXyq[2][j]);
    cnew.SetIndex(npads[cath],cath,TMath::Nint (fXyq[5][j]));
    cnew.SetDetElemId(fInput->DetElemId());
    npads[cath]++;
  }

  cnew.SetClusterType(nover[0] + nover[1] * 100);
  for (Int_t j=0; j<3; j++) cnew.SetTrack(j,tracks[j]);

  Double_t xg, yg, zg;
  for (cath=0; cath<2; cath++) {
    // Perform local-to-global transformation
    fInput->Segmentation2(cath)->GetTransformer()->Local2Global(fInput->DetElemId(), x, y, fZpad, xg, yg, zg);
    cnew.SetX(cath, xg);
    cnew.SetY(cath, yg);
    cnew.SetZ(cath, zg);
    cnew.SetCharge(cath, TMath::Nint(qTot));
    //cnew.SetPeakSignal(cath,20);
    //cnew.SetMultiplicity(cath, 5);
    cnew.SetNcluster(cath, nfit);
    cnew.SetChi2(cath, fmin); //0.;1
  } 
  // Evaluate measurement errors
  //AZ Errors(&cnew);

  cnew.SetGhost(nfit); //cnew.SetX(1,sigx); cnew.SetY(1,sigy); cnew.SetZ(1,dist);
  //cnew.fClusterType=cnew.PhysicsContribution();
  new((*fRawClusters)[fNRawClusters++]) AliMUONRawCluster(cnew); 
  if (fDebug) cout << fNRawClusters << " " << fInput->Chamber() << endl;
  //fNPeaks++;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::FindLocalMaxima(TObjArray *pixArray, Int_t *localMax, Double_t *maxVal)
{
/// Find local maxima in pixel space for large preclusters in order to
/// try to split them into smaller pieces (to speed up the MLEM procedure)
/// or to find additional fitting seeds if clusters were not completely resolved  

  TH2D *hist = NULL;
  //if (pixArray == fPixArray) hist = (TH2D*) gROOT->FindObject("anode");
  //else { hist = (TH2D*) gROOT->FindObject("anode1"); cout << hist << endl; }
  //if (hist) hist->Delete();

  Double_t xylim[4] = {999, 999, 999, 999};
  Int_t nPix = pixArray->GetEntriesFast();
  AliMUONPixel *pixPtr = 0;
  for (Int_t ipix=0; ipix<nPix; ipix++) {
    pixPtr = (AliMUONPixel*) pixArray->UncheckedAt(ipix);
    for (Int_t i=0; i<4; i++) 
         xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
  }
  for (Int_t i=0; i<4; i++) xylim[i] -= pixPtr->Size(i/2); 

  Int_t nx = TMath::Nint ((-xylim[1]-xylim[0])/pixPtr->Size(0)/2);
  Int_t ny = TMath::Nint ((-xylim[3]-xylim[2])/pixPtr->Size(1)/2);
  if (pixArray == fPixArray) hist = new TH2D("anode","anode",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  else hist = new TH2D("anode1","anode1",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  for (Int_t ipix=0; ipix<nPix; ipix++) {
    pixPtr = (AliMUONPixel*) pixArray->UncheckedAt(ipix);
    hist->Fill(pixPtr->Coord(0), pixPtr->Coord(1), pixPtr->Charge());
  }
  if (fDraw && pixArray == fPixArray) fDraw->DrawHist("c2", hist);

  Int_t nMax = 0, indx;
  Int_t *isLocalMax = new Int_t[ny*nx];
  for (Int_t i=0; i<ny*nx; i++) isLocalMax[i] = 0;

  for (Int_t i=1; i<=ny; i++) {
    indx = (i-1) * nx;
    for (Int_t j=1; j<=nx; j++) {
      if (hist->GetCellContent(j,i) < 0.5) continue;
      //if (isLocalMax[indx+j-1] < 0) continue;
      if (isLocalMax[indx+j-1] != 0) continue;
      FlagLocalMax(hist, i, j, isLocalMax);
    }
  }

  for (Int_t i=1; i<=ny; i++) {
    indx = (i-1) * nx;
    for (Int_t j=1; j<=nx; j++) {
      if (isLocalMax[indx+j-1] > 0) { 
	localMax[nMax] = indx + j - 1; 
	maxVal[nMax++] = hist->GetCellContent(j,i);
	if (nMax > 99) AliFatal(" Too many local maxima !!!");
      }
    }
  }
  if (fDebug) cout << " Local max: " << nMax << endl;
  delete [] isLocalMax; isLocalMax = 0;
  return nMax;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax)
{
/// Flag pixels (whether or not local maxima)

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Int_t cont = TMath::Nint (hist->GetCellContent(j,i));
  Int_t cont1 = 0, indx = (i-1)*nx+j-1, indx1 = 0, indx2 = 0;

  for (Int_t i1=i-1; i1<i+2; i1++) {
    if (i1 < 1 || i1 > ny) continue;
    indx1 = (i1 - 1) * nx;
    for (Int_t j1=j-1; j1<j+2; j1++) {
      if (j1 < 1 || j1 > nx) continue;
      if (i == i1 && j == j1) continue;
      indx2 = indx1 + j1 - 1;
      cont1 = TMath::Nint (hist->GetCellContent(j1,i1));
      if (cont < cont1) { isLocalMax[indx] = -1; return; }
      else if (cont > cont1) isLocalMax[indx2] = -1;
      else { // the same charge
	isLocalMax[indx] = 1; 
	if (isLocalMax[indx2] == 0) {
	  FlagLocalMax(hist, i1, j1, isLocalMax);
	  if (isLocalMax[indx2] < 0) { isLocalMax[indx] = -1; return; }
	  else isLocalMax[indx2] = -1;
	}
      } 
    }
  }
  isLocalMax[indx] = 1; // local maximum
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FindCluster(Int_t *localMax, Int_t iMax)
{
/// Find pixel cluster around local maximum \a iMax and pick up pads
/// overlapping with it

  TH2D *hist = (TH2D*) gROOT->FindObject("anode");
  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Int_t ic = localMax[iMax] / nx + 1;
  Int_t jc = localMax[iMax] % nx + 1;
  Bool_t *used = new Bool_t[ny*nx];
  for (Int_t i=0; i<ny*nx; i++) used[i] = kFALSE;

  // Drop all pixels from the array - pick up only the ones from the cluster
  fPixArray->Delete();

  Double_t wx = hist->GetXaxis()->GetBinWidth(1)/2; 
  Double_t wy = hist->GetYaxis()->GetBinWidth(1)/2;  
  Double_t yc = hist->GetYaxis()->GetBinCenter(ic);
  Double_t xc = hist->GetXaxis()->GetBinCenter(jc);
  Double_t cont = hist->GetCellContent(jc,ic);
  AliMUONPixel *pixPtr = new AliMUONPixel (xc, yc, wx, wy, cont);
  fPixArray->Add((TObject*)pixPtr);
  used[(ic-1)*nx+jc-1] = kTRUE;
  AddBin(hist, ic, jc, 1, used, (TObjArray*)0); // recursive call

  Int_t nPix = fPixArray->GetEntriesFast(), npad = fnPads[0] + fnPads[1];
  for (Int_t i=0; i<nPix; i++) {
    ((AliMUONPixel*)fPixArray->UncheckedAt(i))->SetSize(0,wx); 
    ((AliMUONPixel*)fPixArray->UncheckedAt(i))->SetSize(1,wy); 
  }
  if (fDebug) cout << iMax << " " << nPix << endl;

  Float_t xy[4], xy12[4];
  // Pick up pads which overlap with found pixels
  for (Int_t i=0; i<npad; i++) fPadIJ[1][i] = -1;
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    for (Int_t j=0; j<4; j++) 
      xy[j] = pixPtr->Coord(j/2) + (j%2 ? 1 : -1)*pixPtr->Size(j/2);
    for (Int_t j=0; j<npad; j++) 
      if (Overlap(xy, j, xy12, 0)) fPadIJ[1][j] = 0; // flag for use
  }

  delete [] used; used = 0;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddVirtualPad()
{
/// Add virtual pad (with small charge) to improve fit for some
/// clusters (when pad with max charge is at the extreme of the cluster)

  // Get number of pads in X and Y-directions
  Int_t nInX = -1, nInY;
  PadsInXandY(nInX, nInY);
  //return;

  // Add virtual pad only if number of pads per direction == 2
  if (nInX != 2 && nInY != 2) return;

  // Find pads with max charge
  Int_t maxpad[2][2] = {{-1, -1}, {-1, -1}}, cath;
  Double_t sigmax[2] = {0}, aamax[2] = {0};
  for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
    if (fPadIJ[1][j] != 0) continue;
    cath = fPadIJ[0][j];
    if (fXyq[2][j] > sigmax[cath]) {
      maxpad[cath][1] = maxpad[cath][0];
      aamax[cath] = sigmax[cath];
      sigmax[cath] = fXyq[2][j];
      maxpad[cath][0] = j;
    }
  }
  if (maxpad[0][0] >= 0 && maxpad[0][1] < 0 || maxpad[1][0] >= 0 && maxpad[1][1] < 0) {
    for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
      if (fPadIJ[1][j] != 0) continue;
      cath = fPadIJ[0][j];
      if (j == maxpad[cath][0] || j == maxpad[cath][1]) continue;
      if (fXyq[2][j] > aamax[cath]) {
	aamax[cath] = fXyq[2][j];
	maxpad[cath][1] = j;
      }
    }
  }
  // Check for mirrors (side X on cathode 0) 
  Bool_t mirror = kFALSE;
  if (maxpad[0][0] >= 0 && maxpad[1][0] >= 0) {
    mirror = fXyq[3][maxpad[0][0]] < fXyq[4][maxpad[0][0]]; 
    if (!mirror && TMath::Abs(fXyq[3][maxpad[0][0]]-fXyq[3][maxpad[1][0]]) < 0.001) {
      // Special case when pads on both cathodes have the same size
      Int_t yud[2] = {0};
      for (Int_t j = 0; j < fnPads[0]+fnPads[1]; j++) {
	cath = fPadIJ[0][j];
	if (j == maxpad[cath][0]) continue;
	if (fPadIJ[2][j] != fPadIJ[2][maxpad[cath][0]]) continue;
	if (fPadIJ[3][j] + 1 == fPadIJ[3][maxpad[cath][0]] || 
	    fPadIJ[3][j] - 1 == fPadIJ[3][maxpad[cath][0]]) yud[cath]++;
      }
      if (!yud[0]) mirror = kTRUE; // take the other cathode
    } // if (!mirror &&...
  } // if (maxpad[0][0] >= 0 && maxpad[1][0] >= 0)

  // Find neughbours of pads with max charges
  Int_t nn, xList[10], yList[10], ix0, iy0, ix, iy, neighb;
  for (cath=0; cath<2; cath++) {
    if (!cath && maxpad[0][0] < 0) continue; // one-sided cluster - cathode 1
    if (cath && maxpad[1][0] < 0) break; // one-sided cluster - cathode 0
    if (maxpad[1][0] >= 0) {
      if (!mirror) {
	if (!cath && nInY != 2) continue;
	if (cath && nInX != 2 && (maxpad[0][0] >= 0 || nInY != 2)) continue;
      } else {
	if (!cath && nInX != 2) continue;
	if (cath && nInY != 2 && (maxpad[0][0] >= 0 || nInX != 2)) continue;
      }
    }

    Int_t iAddX = 0, iAddY = 0, ix1 = 0, iy1 = 0, iPad = 0;
    if (maxpad[0][0] < 0) iPad = 1;

    for (iPad=0; iPad<2; iPad++) {
      if (maxpad[cath][iPad] < 0) continue;
      if (iPad && !iAddX && !iAddY) break;
      if (iPad && fXyq[2][maxpad[cath][1]] / sigmax[cath] < 0.5) break;

      Int_t neighbx = 0, neighby = 0;
      ix0 = fPadIJ[2][maxpad[cath][iPad]];
      iy0 = fPadIJ[3][maxpad[cath][iPad]];
      fSegmentation[cath]->Neighbours(ix0, iy0, &nn, xList, yList); 
      Float_t zpad; 
      for (Int_t j=0; j<nn; j++) {
	if (TMath::Abs(xList[j]-ix0) == 1 || xList[j]*ix0 == -1) neighbx++;
	if (TMath::Abs(yList[j]-iy0) == 1 || yList[j]*iy0 == -1) neighby++;
      }
      if (!mirror) {
	if (cath) neighb = neighbx; 
	else neighb = neighby;
	if (maxpad[0][0] < 0) neighb += neighby;
	else if (maxpad[1][0] < 0) neighb += neighbx;
      } else {
	if (!cath) neighb = neighbx; 
	else neighb = neighby;
	if (maxpad[0][0] < 0) neighb += neighbx;
	else if (maxpad[1][0] < 0) neighb += neighby;
      }

      for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
	if (fPadIJ[0][j] != cath) continue;
	ix = fPadIJ[2][j];
	iy = fPadIJ[3][j];
	if (iy == iy0 && ix == ix0) continue; 
	for (Int_t k=0; k<nn; k++) {
	  if (xList[k] != ix || yList[k] != iy) continue;
	  if (!mirror) {
	    if ((!cath || maxpad[0][0] < 0) && 
		(TMath::Abs(iy-iy0) == 1 || iy*iy0 == -1)) {
	      if (!iPad && TMath::Abs(ix-ix0) == 1 || ix*ix0 == -1) ix1 = xList[k]; //19-12-05 
	      xList[k] = yList[k] = 0; 
	      neighb--;
	      break;
	    }
	    if ((cath || maxpad[1][0] < 0) && 
		(TMath::Abs(ix-ix0) == 1 || ix*ix0 == -1)) {
	      if (!iPad) ix1 = xList[k]; //19-12-05
	      xList[k] = yList[k] = 0; 
	      neighb--;
	    }
	  } else {
	    if ((!cath || maxpad[0][0] < 0) && 
		(TMath::Abs(ix-ix0) == 1 || ix*ix0 == -1)) {
	      if (!iPad) ix1 = xList[k]; //19-12-05
	      xList[k] = yList[k] = 0; 
	      neighb--;
	      break;
	    }
	    if ((cath || maxpad[1][0] < 0) && 
		(TMath::Abs(iy-iy0) == 1 || iy*iy0 == -1)) {
	      xList[k] = yList[k] = 0; 
	      neighb--;
	    }
	  }
	  break;
	} // for (Int_t k=0; k<nn;
	if (!neighb) break;
      } // for (Int_t j=0; j<fnPads[0]+fnPads[1];
      if (!neighb) continue;
      
      // Add virtual pad
      Int_t npads, isec;
      isec = 0;
      for (Int_t j=0; j<nn; j++) {
	if (xList[j] == 0 && yList[j] == 0) continue;
	npads = fnPads[0] + fnPads[1];
	fPadIJ[0][npads] = cath;
	fPadIJ[1][npads] = 0;
	ix = xList[j]; 
	iy = yList[j];
	if (TMath::Abs(ix-ix0) == 1 || ix*ix0 == -1) {
	  if (iy != iy0) continue; // new segmentation - check
	  if (nInX != 2) continue; // new
	  if (!mirror) {
	    if (!cath && maxpad[1][0] >= 0) continue;
	  } else {
	    if (cath && maxpad[0][0] >= 0) continue;
	  }
	  if (iPad && !iAddX) continue;
          fSegmentation[cath]->GetPadC(ix, iy, fXyq[0][npads], fXyq[1][npads], zpad);
	  if (fXyq[0][npads] > 1.e+5) continue; // temporary fix
	  if (ix == ix1) continue; //19-12-05
	  if (ix1 == ix0) continue;
	  if (maxpad[1][0] < 0 || mirror && maxpad[0][0] >= 0) {
	    if (!iPad) fXyq[2][npads] = TMath::Min (sigmax[0]/100, 5.);
	    else fXyq[2][npads] = TMath::Min (aamax[0]/100, 5.);
	  }
	  else {
	    if (!iPad) fXyq[2][npads] = TMath::Min (sigmax[1]/100, 5.);
	    else fXyq[2][npads] = TMath::Min (aamax[1]/100, 5.);
	  }
	  fXyq[2][npads] = TMath::Max (fXyq[2][npads], (float)1);
	  fXyq[3][npads] = -2; // flag
	  fPadIJ[2][npads] = ix;
	  fPadIJ[3][npads] = iy;
	  fnPads[1]++;
	  iAddX = npads;
	  if (fDebug) printf(" ***** Add virtual pad in X ***** %f %f %f %3d %3d \n", fXyq[2][npads], 
			     fXyq[0][npads], fXyq[1][npads], ix, iy);
	  ix1 = ix0;
	  continue;
	} 
	if (nInY != 2) continue;
	if (!mirror && cath && maxpad[0][0] >= 0) continue;
	if (mirror && !cath && maxpad[1][0] >= 0) continue;
	if (TMath::Abs(iy-iy0) == 1 || TMath::Abs(iy*iy0) == 1) {
	  if (ix != ix0) continue; // new segmentation - check
	  if (iPad && !iAddY) continue;
          fSegmentation[cath]->GetPadC(ix, iy, fXyq[0][npads], fXyq[1][npads], zpad);
	  if (iy1 == iy0) continue;
	  //if (iPad && iy1 == iy0) continue;
	  if (maxpad[0][0] < 0 || mirror && maxpad[1][0] >= 0) {
	    if (!iPad) fXyq[2][npads] = TMath::Min (sigmax[1]/15, fgkZeroSuppression);
	    else fXyq[2][npads] = TMath::Min (aamax[1]/15, fgkZeroSuppression);
	  }
	  else {
	    if (!iPad) fXyq[2][npads] = TMath::Min (sigmax[0]/15, fgkZeroSuppression);
	    else fXyq[2][npads] = TMath::Min (aamax[0]/15, fgkZeroSuppression);
	  }
	  fXyq[2][npads] = TMath::Max (fXyq[2][npads], (float)1);
	  fXyq[3][npads] = -2; // flag
	  fPadIJ[2][npads] = ix;
	  fPadIJ[3][npads] = iy;
	  fnPads[1]++;
	  iAddY = npads;
	  if (fDebug) printf(" ***** Add virtual pad in Y ***** %f %f %f %3d %3d \n", fXyq[2][npads], 
			     fXyq[0][npads], fXyq[1][npads], ix, iy);
	  iy1 = iy0;
	}
      } // for (Int_t j=0; j<nn;
    } // for (Int_t iPad=0;
  } // for (cath=0; cath<2;
  return;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::PadsInXandY(Int_t &nInX, Int_t &nInY)
{
/// Find number of pads in X and Y-directions (excluding virtual ones and
/// overflows)

  static Int_t nXsaved = 0, nYsaved = 0;
  nXsaved = nYsaved = 0;
  //if (nInX >= 0) {nInX = nXsaved; nInY = nYsaved; return; }
  Float_t *xPad0 = NULL, *yPad0 = NULL, *xPad1 = NULL, *yPad1 = NULL;
  Float_t wMinX[2] = {99, 99}, wMinY[2] = {99, 99};
  Int_t *nPad0 = NULL, *nPad1 = NULL;
  Int_t nPads = fnPads[0] + fnPads[1];
  if (fnPads[0]) {
    xPad0 = new Float_t[nPads];
    yPad0 = new Float_t[nPads];
    nPad0 = new Int_t[nPads];
  }
  if (fnPads[1]) {
    xPad1 = new Float_t[nPads];
    yPad1 = new Float_t[nPads];
    nPad1 = new Int_t[nPads];
  }
  Int_t n0 = 0, n1 = 0, cath, npadx[2] = {1, 1}, npady[2] = {1, 1};
  for (Int_t j = 0; j < nPads; j++) {
    if (nInX < 0 && fPadIJ[1][j] != 0) continue; // before fit
    else if (nInX == 0 && fPadIJ[1][j] != 1) continue; // fit - exclude overflows
    else if (nInX > 0 && fPadIJ[1][j] != 1 && fPadIJ[1][j] != -9) continue; // exclude non-marked
    if (nInX <= 0 && fXyq[2][j] > fgkSaturation-1) continue; // skip overflows
    cath = fPadIJ[0][j];
    if (fXyq[3][j] > 0) { // exclude virtual pads
      wMinX[cath] = TMath::Min (wMinX[cath], fXyq[3][j]);
      wMinY[cath] = TMath::Min (wMinY[cath], fXyq[4][j]);
      //20-12-05 }
      if (cath) { xPad1[n1] = fXyq[0][j]; yPad1[n1++] = fXyq[1][j]; }
      else { xPad0[n0] = fXyq[0][j]; yPad0[n0++] = fXyq[1][j]; }
    }
  }

  // Sort
  if (n0) { 
    TMath::Sort (n0, xPad0, nPad0); // in X
    for (Int_t i = 1; i < n0; i++) 
      if (xPad0[nPad0[i]] - xPad0[nPad0[i-1]] < -0.01) npadx[0]++;
    TMath::Sort (n0, yPad0, nPad0); // in Y
    for (Int_t i = 1; i < n0; i++) 
      if (yPad0[nPad0[i]] - yPad0[nPad0[i-1]] < -0.01) npady[0]++;
  }
  
  if (n1) { 
    TMath::Sort (n1, xPad1, nPad1); // in X
    for (Int_t i = 1; i < n1; i++) 
      if (xPad1[nPad1[i]] - xPad1[nPad1[i-1]] < -0.01) npadx[1]++;
    TMath::Sort (n1, yPad1, nPad1); // in Y
    for (Int_t i = 1; i < n1; i++) 
      if (yPad1[nPad1[i]] - yPad1[nPad1[i-1]] < -0.01) npady[1]++;
  }
  if (fnPads[0]) { delete [] xPad0; delete [] yPad0; delete [] nPad0; }
  if (fnPads[1]) { delete [] xPad1; delete [] yPad1; delete [] nPad1; }
  if (TMath::Abs (wMinY[0] - wMinY[1]) < 1.e-3) nInY = TMath::Max (npady[0], npady[1]);
  else nInY = wMinY[0] < wMinY[1] ? npady[0] : npady[1];
  if (TMath::Abs (wMinX[0] - wMinX[1]) < 1.e-3) nInX = TMath::Max (npadx[0], npadx[1]);
  else nInX = wMinX[0] < wMinX[1] ? npadx[0] : npadx[1];
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Simple()
{
/// Process simple cluster (small number of pads) without EM-procedure

  Int_t nForFit = 1, clustFit[1] = {0}, nfit;
  Double_t parOk[3] = {0.}; 
  TObjArray *clusters[1]; 
  clusters[0] = fPixArray;
  for (Int_t i = 0; i < fnPads[0]+fnPads[1]; i++) {
    if (fXyq[2][i] > fgkSaturation-1) fPadIJ[1][i] = -9;
    else fPadIJ[1][i] = 1;
  }
  nfit = Fit(1, nForFit, clustFit, clusters, parOk);
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Errors(AliMUONRawCluster *clus)
{
/// Correct reconstructed coordinates for some clusters and evaluate errors

  Double_t qTot = clus->GetCharge(0), fmin = clus->GetChi2(0);
  Double_t xreco = clus->GetX(0), yreco = clus->GetY(0), zreco = clus->GetZ(0);
  Double_t sigmax[2] = {0};

  Int_t nInX = 1, nInY, maxdig[2] ={-1, -1}, digit, cath1, isec;
  PadsInXandY(nInX, nInY);

  // Find pad with maximum signal
  for (Int_t cath = 0; cath < 2; cath++) {
    for (Int_t j = 0; j < clus->GetMultiplicity(cath); j++) {
      cath1 = cath;
      digit = clus->GetIndex(j, cath);
      if (digit < 0) { cath1 = TMath::Even(cath); digit = -digit - 1; } // from the other cathode

      if (clus->GetContrib(j,cath) > sigmax[cath1]) {
	sigmax[cath1] = clus->GetContrib(j,cath);
	maxdig[cath1] = digit;
      }
    }
  }

  // Size of pad with maximum signal and reco coordinate distance from the pad center
  AliMUONDigit *mdig = 0;
  Double_t wx[2], wy[2], dxc[2], dyc[2];
  Float_t xpad, ypad, zpad;
  Int_t ix, iy;
  for (Int_t cath = 0; cath < 2; cath++) {
    if (maxdig[cath] < 0) continue;
    mdig = fInput->Digit(cath,maxdig[cath]);
    isec = fSegmentation[cath]->Sector(mdig->PadX(), mdig->PadY());
    wx[cath] = fSegmentation[cath]->Dpx(isec);
    wy[cath] = fSegmentation[cath]->Dpy(isec);
    fSegmentation[cath]->GetPadI(xreco, yreco, zreco, ix, iy);
    isec = fSegmentation[cath]->Sector(ix, iy);
    if (isec > 0) {
      fSegmentation[cath]->GetPadC(ix, iy, xpad, ypad, zpad);
      dxc[cath] = xreco - xpad;
      dyc[cath] = yreco - ypad;
    }
  }

  // Check if pad with max charge at the edge (number of neughbours)
  Int_t nn, xList[10], yList[10], neighbx[2][2] = {{0,0}, {0,0}}, neighby[2][2]= {{0,0}, {0,0}};
  for (Int_t cath = 0; cath < 2; cath++) {
    if (maxdig[cath] < 0) continue;
    mdig = fInput->Digit(cath,maxdig[cath]);
    fSegmentation[cath]->Neighbours(mdig->PadX(), mdig->PadY(), &nn, xList, yList); 
    isec = fSegmentation[cath]->Sector(mdig->PadX(), mdig->PadY());
    /*??
    Float_t sprX = fResponse->SigmaIntegration() * fResponse->ChargeSpreadX();
    Float_t sprY = fResponse->SigmaIntegration() * fResponse->ChargeSpreadY();
    //fSegmentation[cath]->FirstPad(fInput->DetElemId(),muons[ihit][1], muons[ihit][2], muons[ihit][3], sprX, sprY);  
    //fSegmentation[cath]->FirstPad(fInput->DetElemId(),xreco, yreco, zreco, sprX, sprY);  
    fSegmentation[cath]->FirstPad(xreco, yreco, zreco, sprX, sprY);  
    Int_t border = 0;
    //if (fSegmentation[cath]->Sector(fInput->DetElemId(),fSegmentation[cath]->Ix(),fSegmentation[cath]->Iy()) <= 0) {
    if (fSegmentation[cath]->Sector(fSegmentation[cath]->Ix(), fSegmentation[cath]->Iy()) <= 0) {
      //fSegmentation[cath]->NextPad(fInput->DetElemId());
      fSegmentation[cath]->NextPad();
      border = 1;
    } 
    */
    for (Int_t j=0; j<nn; j++) {
      //if (border && yList[j] < fSegmentation[cath]->Iy()) continue;
      fSegmentation[cath]->GetPadC(xList[j], yList[j], xpad, ypad, zpad);
      //cout << ch << " " << xList[j] << " " << yList[j] << " " << border << " " << x << " " << y << " " << xpad << " " << ypad << endl;
      if (TMath::Abs(xpad) < 1 && TMath::Abs(ypad) < 1) continue;
      if (xList[j] == mdig->PadX()-1 || mdig->PadX() == 1 && 
	  xList[j] == -1) neighbx[cath][0] = 1;
      else if (xList[j] == mdig->PadX()+1 || mdig->PadX() == -1 && 
	       xList[j] == 1) neighbx[cath][1] = 1;
      if (yList[j] == mdig->PadY()-1 || mdig->PadY() == 1 &&
	  yList[j] == -1) neighby[cath][0] = 1;
      else if (yList[j] == mdig->PadY()+1 || mdig->PadY() == -1 &&
	       yList[j] == 1) neighby[cath][1] = 1;
    } // for (Int_t j=0; j<nn;
    if (neighbx[cath][0] && neighbx[cath][1]) neighbx[cath][0] = 0;
    else if (neighbx[cath][1]) neighbx[cath][0] = -1;
    else neighbx[cath][0] = 1;
    if (neighby[cath][0] && neighby[cath][1]) neighby[cath][0] = 0;
    else if (neighby[cath][1]) neighby[cath][0] = -1;
    else neighby[cath][0] = 1;
  }

  Int_t iOver = clus->GetClusterType();
  // One-sided cluster
  if (!clus->GetMultiplicity(0)) {
    neighby[0][0] = neighby[1][0];
    wy[0] = wy[1];
    if (iOver < 99) iOver += 100 * iOver;
    dyc[0] = dyc[1];
  } else if (!clus->GetMultiplicity(1)) {
    neighbx[1][0] = neighbx[0][0];
    wx[1] = wx[0];
    if (iOver < 99) iOver += 100 * iOver;
    dxc[1] = dxc[0];
  }

  // Apply corrections and evaluate errors
  Double_t errY, errX;
  Errors(nInY, nInX, neighby[0][0],neighbx[1][0], fmin, wy[0]*10, wx[1]*10, iOver, 
	 dyc[0], dxc[1], qTot, yreco, xreco, errY, errX);
  errY = TMath::Max (errY, 0.01);
  //errY = 0.01;
  //errX = TMath::Max (errX, 0.144);
  clus->SetX(0, xreco); clus->SetY(0, yreco);
  clus->SetErrX(errX); clus->SetErrY(errY);
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Errors(Int_t ny, Int_t nx, Int_t iby, Int_t ibx, Double_t fmin,
				    Double_t wy, Double_t wx, Int_t iover, 
				    Double_t dyc, Double_t /*dxc*/, Double_t qtot, 
				    Double_t &yrec, Double_t &xrec, Double_t &erry, Double_t &errx)
{
/// Correct reconstructed coordinates for some clusters and evaluate errors

    erry = 0.01;
    errx = 0.144;
    Int_t iovery = iover % 100;
    Double_t corr = 0;

/* ---> Ny = 1 */
    if (ny == 1) {
      if (iby != 0) {
	// edge effect 
	yrec += iby * (0.1823+0.2008)/2;
	erry = 0.04587;
      } else {
	// Find "effective pad width" 
	Double_t width = 0.218 / (1.31e-4 * TMath::Exp (2.688 * TMath::Log(qtot)) + 1) * 2;
	width = TMath::Min (width, 0.4);
	erry = width / TMath::Sqrt(12.);
	erry = TMath::Max (erry, 0.01293);
      }
      goto x; //return;
    }

/* ---> "Bad" fit */
    if (fmin > 0.4) {
      erry = 0.1556;
      if (ny == 5) erry = 0.06481;
      goto x; //return;
    }

/* ---> By != 0 */
    if (iby != 0) {
      if (ny > 2) {
	erry = 0.00417; //0.01010
      } else {
        // ny = 2 
	if (dyc * iby > -0.05) {
	  Double_t dyc2 = dyc * dyc;
	  if (iby < 0) {
	    corr = 0.019 - 0.602 * dyc + 8.739 * dyc2 - 44.209 * dyc2 * dyc;
	    corr = TMath::Min (corr, TMath::Abs(-0.25-dyc));
	    yrec -= corr;
	    //dyc -= corr;
	    erry = 0.00814;
	  } else {
	    corr = 0.006 + 0.300 * dyc + 6.147 * dyc2 + 42.039 * dyc2 * dyc;
	    corr = TMath::Min (corr, 0.25-dyc);
	    yrec += corr;
	    //dyc += corr;
	    erry = 0.01582;
	  }
	} else {
	  erry = (0.00303 + 0.00296) / 2;
	}
      }
      goto x; //return;
    }

/* ---> Overflows */
    if (iovery != 0) {
      if (qtot < 3000) {
	erry = 0.0671;
      } else {
	if (iovery > 1) {
	  erry = 0.09214;
	} else if (TMath::Abs(wy - 5) < 0.1) {
	  erry = 0.061; //0.06622
	} else {
	  erry = 0.00812; // 0.01073 
	}
      }
      goto x; //return;
    }

/* ---> "Good" but very high signal */
    if (qtot > 4000) {
      if (TMath::Abs(wy - 4) < 0.1) {
	erry = 0.00117;
      } else if (fmin < 0.03 && qtot < 6000) {
	erry = 0.01003;
      } else {
	erry = 0.1931;
      }
      goto x; //return;
    }

/* ---> "Good" clusters */
    if (ny > 3) {
      if (TMath::Abs(wy - 5) < 0.1) {
	erry = 0.0011; //0.00304 
      } else if (qtot < 400.) {
	erry = 0.0165;
      } else {
	erry = 0.00135; // 0.00358 
      }
    } else if (ny == 3) {
      if (TMath::Abs(wy - 4) < 0.1) {
	erry = 35.407 / (1 + TMath::Exp(5.511*TMath::Log(qtot/265.51))) + 11.564;
	//erry = 83.512 / (1 + TMath::Exp(3.344*TMath::Log(qtot/211.58))) + 12.260;
      } else {
	erry = 147.03 / (1 + TMath::Exp(1.713*TMath::Log(qtot/73.151))) + 9.575;
	//erry = 91.743 / (1 + TMath::Exp(2.332*TMath::Log(qtot/151.67))) + 11.453;
      }
      erry *= 1.e-4;
    } else {
      // ny = 2 
      if (TMath::Abs(wy - 4) < 0.1) {
	erry = 60.800 / (1 + TMath::Exp(3.305*TMath::Log(qtot/104.53))) + 11.702;
	//erry = 73.128 / (1 + TMath::Exp(5.676*TMath::Log(qtot/120.93))) + 17.839;
      } else {
	erry = 117.98 / (1 + TMath::Exp(2.005*TMath::Log(qtot/37.649))) + 21.431;
	//erry = 99.066 / (1 + TMath::Exp(4.900*TMath::Log(qtot/107.57))) + 25.315;
      }
      erry *= 1.e-4;
    }
    //return;

 x:
/* ---> X-coordinate */
/* ---> Y-side */    
    if (wx > 11) { 
      errx = 0.0036;
      xrec -= 0.1385;
      return;
    }
/* ---> Nx = 1 */
    if (nx == 1) {
      if (TMath::Abs(wx - 6) < 0.1) {
	if (qtot < 40) errx = 0.1693;
	else errx = 0.06241;
      } else if (TMath::Abs(wx - 7.5) < 0.1) {
	if (qtot < 40) errx = 0.2173;
	else errx = 0.07703;
      } else if (TMath::Abs(wx - 10) < 0.1) {
	if (ibx == 0) {
	  if (qtot < 40) errx = 0.2316;
	  else errx = 0.1426;
	} else {
	  xrec += (0.2115 + 0.1942) / 2 * ibx;
	  errx = 0.1921;
	}
      }
      return;
    }
/* ---> "Bad" fit */
    if (fmin > 0.5) {
      errx = 0.1591;
      return;
    }
/* ---> Bx != 0 */
    if (ibx != 0) {
      if (ibx > 0) { errx = 0.06761; xrec -= 0.03832; }
      else { errx = 0.06653; xrec += 0.02581; }
      return;
    }
/* ---> Overflows */
    if (iover != 0) {
      if (TMath::Abs(wx - 6) < 0.1) errx = 0.06979;
      else if (TMath::Abs(wx - 7.5) < 0.1) errx = 0.1089;
      else if (TMath::Abs(wx - 10) < 0.1) errx = 0.09847;
      return;
    }
/* ---> Good */
    if (TMath::Abs(wx - 6) < 0.1) errx = 0.06022;
    else if (TMath::Abs(wx - 7.5) < 0.1) errx = 0.07247;
    else if (TMath::Abs(wx - 10) < 0.1) errx = 0.07359;
}

