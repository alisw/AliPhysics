/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
//____________________________________________________________________
//
// AliITSTrackleterSPDEff - find SPD chips efficiencies by using tracklets.
//
// This class has been developed from AliITSMultReconstructor (see
// it for more details). It is the class for the Trackleter used to estimate
// SPD plane efficiency.
// The trackleter prediction is built using the vertex and 1 cluster.
//
//
//  Author :  Giuseppe Eugenio Bruno, based on the skeleton of Reconstruct method  provided by Tiziano Virgili
//  email:    giuseppe.bruno@ba.infn.it
//
//____________________________________________________________________

/* $Id$ */

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TClonesArray.h>

#include "AliTracker.h"
#include "AliITSTrackleterSPDEff.h"
#include "AliITSgeomTGeo.h"
#include "AliLog.h"
#include "AliITSPlaneEffSPD.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliRunLoader.h"
#include "AliITSReconstructor.h"
#include "AliITSRecPoint.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
//____________________________________________________________________
ClassImp(AliITSTrackleterSPDEff)


//____________________________________________________________________
AliITSTrackleterSPDEff::AliITSTrackleterSPDEff():
AliTracker(),
//
fClustersLay1(0),
fClustersLay2(0),
fTracklets(0),
fAssociationFlag(0),
fNClustersLay1(0),
fNClustersLay2(0),
fNTracklets(0),
fOnlyOneTrackletPerC2(0),
fPhiWindowL2(0),
fZetaWindowL2(0),
fPhiOverlapCut(0),
fZetaOverlapCut(0),
fHistOn(0),
fhClustersDPhiAcc(0),
fhClustersDThetaAcc(0),
fhClustersDZetaAcc(0),
fhClustersDPhiAll(0),
fhClustersDThetaAll(0),
fhClustersDZetaAll(0),
fhDPhiVsDThetaAll(0),
fhDPhiVsDThetaAcc(0),
fhDPhiVsDZetaAll(0),
fhDPhiVsDZetaAcc(0),
fhetaTracklets(0),
fhphiTracklets(0),
fhetaClustersLay1(0),
fhphiClustersLay1(0),
//
fAssociationFlag1(0),
fChipPredOnLay2(0),
fChipPredOnLay1(0),
fNTracklets1(0),
fPhiWindowL1(0),
fZetaWindowL1(0),
fOnlyOneTrackletPerC1(0),
fUpdateOncePerEventPlaneEff(0),
fMinContVtx(0),
fChipUpdatedInEvent(0),
fPlaneEffSPD(0),
fPlaneEffBkg(0),
fReflectClusterAroundZAxisForLayer0(kFALSE),
fReflectClusterAroundZAxisForLayer1(kFALSE),
fLightBkgStudyInParallel(kFALSE),
fMC(0),
fUseOnlyPrimaryForPred(0),
fUseOnlySecondaryForPred(0), 
fUseOnlySameParticle(0),
fUseOnlyDifferentParticle(0),
fUseOnlyStableParticle(1),
fPredictionPrimary(0),
fPredictionSecondary(0),
fClusterPrimary(0),
fClusterSecondary(0),
fSuccessPP(0),
fSuccessTT(0),
fSuccessS(0),
fSuccessP(0),
fFailureS(0),
fFailureP(0),
fRecons(0),
fNonRecons(0),
fhClustersDPhiInterpAcc(0),
fhClustersDThetaInterpAcc(0),
fhClustersDZetaInterpAcc(0),
fhClustersDPhiInterpAll(0),
fhClustersDThetaInterpAll(0),
fhClustersDZetaInterpAll(0),
fhDPhiVsDThetaInterpAll(0),
fhDPhiVsDThetaInterpAcc(0),
fhDPhiVsDZetaInterpAll(0),
fhDPhiVsDZetaInterpAcc(0),
fhetaClustersLay2(0),
fhphiClustersLay2(0),
fhClustersInChip(0),
fhClustersInModuleLay1(0),
fhClustersInModuleLay2(0)
{
   // default constructor
// from AliITSMultReconstructor
  Init();
}
//______________________________________________________________________
void AliITSTrackleterSPDEff::Init() {
  SetPhiWindowL2();
  SetZetaWindowL2();
  SetOnlyOneTrackletPerC2();
  fClustersLay1       = new Float_t*[300000];
  fClustersLay2       = new Float_t*[300000];
  fTracklets          = new Float_t*[300000];
  fAssociationFlag    = new Bool_t[300000];
//
  SetPhiWindowL1();
  SetZetaWindowL1();
  SetOnlyOneTrackletPerC1();

  fAssociationFlag1   = new Bool_t[300000];
  fChipPredOnLay2     = new UInt_t[300000];
  fChipPredOnLay1     = new UInt_t[300000];
  fChipUpdatedInEvent = new Bool_t[1200];

  for(Int_t i=0; i<300000; i++) {
    // from AliITSMultReconstructor
    fClustersLay1[i]       = new Float_t[6];
    fClustersLay2[i]       = new Float_t[6];
    fTracklets[i]          = new Float_t[5];
    fAssociationFlag[i]    = kFALSE;
    //
    fAssociationFlag1[i]   = kFALSE;
  }
  for(Int_t i=0;i<1200; i++) fChipUpdatedInEvent[i] = kFALSE;

  if (GetHistOn()) BookHistos();

  fPlaneEffSPD = new AliITSPlaneEffSPD();
  SetLightBkgStudyInParallel();
}
//______________________________________________________________________
AliITSTrackleterSPDEff::AliITSTrackleterSPDEff(const AliITSTrackleterSPDEff &mr) :  
AliTracker(mr),
// from AliITSMultReconstructor
fClustersLay1(mr.fClustersLay1),
fClustersLay2(mr.fClustersLay2),
fTracklets(mr.fTracklets),
fAssociationFlag(mr.fAssociationFlag),
fNClustersLay1(mr.fNClustersLay1),
fNClustersLay2(mr.fNClustersLay2),
fNTracklets(mr.fNTracklets),
fOnlyOneTrackletPerC2(mr.fOnlyOneTrackletPerC2),
fPhiWindowL2(mr.fPhiWindowL2),
fZetaWindowL2(mr.fZetaWindowL2),
fPhiOverlapCut(mr.fPhiOverlapCut),
fZetaOverlapCut(mr.fZetaOverlapCut),
fHistOn(mr.fHistOn),
fhClustersDPhiAcc(mr.fhClustersDPhiAcc),
fhClustersDThetaAcc(mr.fhClustersDThetaAcc),
fhClustersDZetaAcc(mr.fhClustersDZetaAcc),
fhClustersDPhiAll(mr.fhClustersDPhiAll),
fhClustersDThetaAll(mr.fhClustersDThetaAll),
fhClustersDZetaAll(mr.fhClustersDZetaAll),
fhDPhiVsDThetaAll(mr.fhDPhiVsDThetaAll),
fhDPhiVsDThetaAcc(mr.fhDPhiVsDThetaAcc),
fhDPhiVsDZetaAll(mr.fhDPhiVsDZetaAll),
fhDPhiVsDZetaAcc(mr.fhDPhiVsDZetaAcc),
fhetaTracklets(mr.fhetaTracklets),
fhphiTracklets(mr.fhphiTracklets),
fhetaClustersLay1(mr.fhetaClustersLay1),
fhphiClustersLay1(mr.fhphiClustersLay1),
//
fAssociationFlag1(mr.fAssociationFlag1),
fChipPredOnLay2(mr.fChipPredOnLay2),
fChipPredOnLay1(mr.fChipPredOnLay1),
fNTracklets1(mr.fNTracklets1),
fPhiWindowL1(mr.fPhiWindowL1),
fZetaWindowL1(mr.fZetaWindowL1),
fOnlyOneTrackletPerC1(mr.fOnlyOneTrackletPerC1),
fUpdateOncePerEventPlaneEff(mr.fUpdateOncePerEventPlaneEff),
fMinContVtx(mr.fMinContVtx),
fChipUpdatedInEvent(mr.fChipUpdatedInEvent),
fPlaneEffSPD(mr.fPlaneEffSPD),
fPlaneEffBkg(mr.fPlaneEffBkg),
fReflectClusterAroundZAxisForLayer0(mr.fReflectClusterAroundZAxisForLayer0),
fReflectClusterAroundZAxisForLayer1(mr.fReflectClusterAroundZAxisForLayer1),
fLightBkgStudyInParallel(mr.fLightBkgStudyInParallel),
fMC(mr.fMC),
fUseOnlyPrimaryForPred(mr.fUseOnlyPrimaryForPred),
fUseOnlySecondaryForPred(mr.fUseOnlySecondaryForPred),
fUseOnlySameParticle(mr.fUseOnlySameParticle),
fUseOnlyDifferentParticle(mr.fUseOnlyDifferentParticle),
fUseOnlyStableParticle(mr.fUseOnlyStableParticle),
fPredictionPrimary(mr.fPredictionPrimary),
fPredictionSecondary(mr.fPredictionSecondary),
fClusterPrimary(mr.fClusterPrimary),
fClusterSecondary(mr.fClusterSecondary),
fSuccessPP(mr.fSuccessPP),
fSuccessTT(mr.fSuccessTT),
fSuccessS(mr.fSuccessS),
fSuccessP(mr.fSuccessP),
fFailureS(mr.fFailureS),
fFailureP(mr.fFailureP),
fRecons(mr.fRecons),
fNonRecons(mr.fNonRecons),
fhClustersDPhiInterpAcc(mr.fhClustersDPhiInterpAcc),
fhClustersDThetaInterpAcc(mr.fhClustersDThetaInterpAcc),
fhClustersDZetaInterpAcc(mr.fhClustersDZetaInterpAcc),
fhClustersDPhiInterpAll(mr.fhClustersDPhiInterpAll),
fhClustersDThetaInterpAll(mr.fhClustersDThetaInterpAll),
fhClustersDZetaInterpAll(mr.fhClustersDZetaInterpAll),
fhDPhiVsDThetaInterpAll(mr.fhDPhiVsDThetaInterpAll),
fhDPhiVsDThetaInterpAcc(mr.fhDPhiVsDThetaInterpAcc),
fhDPhiVsDZetaInterpAll(mr.fhDPhiVsDZetaInterpAll),
fhDPhiVsDZetaInterpAcc(mr.fhDPhiVsDZetaInterpAcc),
fhetaClustersLay2(mr.fhetaClustersLay2),
fhphiClustersLay2(mr.fhphiClustersLay2),
fhClustersInChip(mr.fhClustersInChip),
fhClustersInModuleLay1(mr.fhClustersInModuleLay1),
fhClustersInModuleLay2(mr.fhClustersInModuleLay2)
{
  // Copy constructor
}

//______________________________________________________________________
AliITSTrackleterSPDEff& AliITSTrackleterSPDEff::operator=(const AliITSTrackleterSPDEff& mr){
  // Assignment operator
  this->~AliITSTrackleterSPDEff();
  new(this) AliITSTrackleterSPDEff(mr);
  return *this;
}
//______________________________________________________________________
AliITSTrackleterSPDEff::~AliITSTrackleterSPDEff(){
  // Destructor
// from AliITSMultReconstructor
 // delete arrays
  for(Int_t i=0; i<300000; i++) {
    delete [] fClustersLay1[i];
    delete [] fClustersLay2[i];
    delete [] fTracklets[i];
  }
  delete [] fClustersLay1;
  delete [] fClustersLay2;
  delete [] fTracklets;
  delete [] fAssociationFlag;
//
  // delete histograms
  DeleteHistos();

  delete [] fAssociationFlag1;

  delete [] fChipPredOnLay2;
  delete [] fChipPredOnLay1;

  delete [] fChipUpdatedInEvent;

  delete [] fPredictionPrimary;  
  delete [] fPredictionSecondary; 
  delete [] fClusterPrimary;  
  delete [] fClusterSecondary; 
  delete [] fSuccessPP;
  delete [] fSuccessTT;
  delete [] fSuccessS;
  delete [] fSuccessP;
  delete [] fFailureS;
  delete [] fFailureP;
  delete [] fRecons;
  delete [] fNonRecons;

  // delete PlaneEff
  delete fPlaneEffSPD;
  fPlaneEffSPD=0;
  if(fPlaneEffBkg) {
    delete fPlaneEffBkg;
    fPlaneEffBkg=0;

  }
}
//____________________________________________________________________
void
AliITSTrackleterSPDEff::Reconstruct(AliStack *pStack, TTree *tRef, Bool_t lbkg) {
  //
  // - you have to take care of the following, before of using Reconstruct
  //   1) call LoadClusters(TTree* cl) that finds the position of the clusters (in global coord)
  //   and  convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). 
  //   2) call SetVertex(vtxPos, vtxErr) which set the position of the vertex
  // - Find the extrapolation/interpolation point.
  // - Find the chip corresponding to that
  // - Check if there is a cluster near that point  
  //
  // reset counters
  if(lbkg && !GetLightBkgStudyInParallel()) {
    AliError("You asked for lightBackground in the Reconstruction without proper call to SetLightBkgStudyInParallel(1)"); 
    return;
  }
  AliITSPlaneEffSPD *pe;
  if(lbkg) {
    pe=fPlaneEffBkg;
  } else {
    pe=fPlaneEffSPD;
  }
  fNTracklets = 0; 
  // retrieve the vertex position
  Float_t vtx[3];
  vtx[0]=(Float_t)GetX();
  vtx[1]=(Float_t)GetY();
  vtx[2]=(Float_t)GetZ();
  // to study residual background (i.e. contribution from TT' to measured efficiency) 
  if(fReflectClusterAroundZAxisForLayer0 && !lbkg) ReflectClusterAroundZAxisForLayer(0);
  if(fReflectClusterAroundZAxisForLayer1 && !lbkg) ReflectClusterAroundZAxisForLayer(1);
  //
  if(fMC && !pStack && !lbkg) {AliError("You asked for MC infos but AliStack not properly loaded"); return;}
  if(fMC && !tRef   && !lbkg) {AliError("You asked for MC infos but TrackRef Tree not properly loaded"); return;}
  Bool_t found;
  Int_t nfTraPred1=0;  Int_t ntTraPred1=0;
  Int_t nfTraPred2=0;  Int_t ntTraPred2=0;
  Int_t nfClu1=0;      Int_t ntClu1=0; 
  Int_t nfClu2=0;      Int_t ntClu2=0;
  
  // Set fChipUpdatedInEvent=kFALSE for all the chips (none of the chip efficiency already updated 
  // for this new event)
  for(Int_t i=0;i<1200;i++) fChipUpdatedInEvent[i] = kFALSE;

  // find the tracklets
  AliDebug(1,"Looking for tracklets... ");  
  AliDebug(1,Form("Reconstruct: vtx[0] = %f, vtx[1] = %f, vtx[2] = %f",vtx[0],vtx[1],vtx[2]));

  //###########################################################
  // Loop on layer 1 : finding theta, phi and z 
  UInt_t key;
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    
    Float_t x = fClustersLay1[iC1][0] - vtx[0];
    Float_t y = fClustersLay1[iC1][1] - vtx[1];
    Float_t z = fClustersLay1[iC1][2] - vtx[2];

    Float_t r    = TMath::Sqrt(x*x + y*y +z*z); 
    
    fClustersLay1[iC1][0] = TMath::ACos(z/r);                   // Store Theta
    fClustersLay1[iC1][1] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi
    fClustersLay1[iC1][2] = z;                  // Store z

    // find the Radius and the chip corresponding to the extrapolation point

    found=FindChip(key, 1, vtx, fClustersLay1[iC1][0],fClustersLay1[iC1][1]);
    if (!found) {
      AliDebug(1,Form("Reconstruct: cannot find chip prediction on outer layer for cluster %d on the inner layer",iC1)); 
      key=999999;               // also some other actions should be taken if not Found 
    }
    nfTraPred2+=(Int_t)found; // this for debugging purpose
    ntTraPred2++;             // to check efficiency of the method FindChip
    fChipPredOnLay2[iC1] = key;
    fAssociationFlag1[iC1] = kFALSE;
 
    if (fHistOn && !lbkg) {
      Float_t eta=fClustersLay1[iC1][0];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaClustersLay1->Fill(eta);
      fhphiClustersLay1->Fill(fClustersLay1[iC1][1]);
      fhClustersInChip->Fill(fhClustersInChip->GetBinCenter(key+1)); // if found=kFALSE -> overflow
    }      
  }
  // Loop on layer 2 : finding theta, phi and r   
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {    
    Float_t x = fClustersLay2[iC2][0] - vtx[0];
    Float_t y = fClustersLay2[iC2][1] - vtx[1];
    Float_t z = fClustersLay2[iC2][2] - vtx[2];
   
    Float_t r    = TMath::Sqrt(x*x + y*y +z*z);
    
    fClustersLay2[iC2][0] = TMath::ACos(z/r);                   // Store Theta
    fClustersLay2[iC2][1] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi (done properly in the range [0,2pi])
    fClustersLay2[iC2][2] = z;                  // Store z

    // find the Radius and the chip corresponding to the extrapolation point

    found=FindChip(key, 0, vtx, fClustersLay2[iC2][0],fClustersLay2[iC2][1]);
    if (!found) {
      AliDebug(1,Form("Reconstruct: cannot find chip prediction on inner layer for cluster %d on the outer layer",iC2)); 
      key=999999;
    }
    nfTraPred1+=(Int_t)found; // this for debugging purpose
    ntTraPred1++;             // to check efficiency of the method FindChip
    fChipPredOnLay1[iC2] = key;
    fAssociationFlag[iC2] = kFALSE;
 
    if (fHistOn && !lbkg) {
      Float_t eta=fClustersLay2[iC2][0];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaClustersLay2->Fill(eta);
      fhphiClustersLay2->Fill(fClustersLay2[iC2][1]);
      fhClustersInChip->Fill(fhClustersInChip->GetBinCenter(key+1)); // if found=kFALSE -> overflow
    }
  }  
  
  //###########################################################

 // First part : Extrapolation to Layer 2 

  // Loop on layer 1 
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    

    // here the control to check whether the efficiency of the chip traversed by this tracklet
    // prediction has already been updated in this event using another tracklet prediction
    if(fUpdateOncePerEventPlaneEff && fChipPredOnLay2[iC1]<1200 && fChipUpdatedInEvent[fChipPredOnLay2[iC1]]) continue;
  
    // reset of variables for multiple candidates
    Int_t  iC2WithBestDist = 0;     // reset 
    Float_t distmin        = 100.;  // just to put a huge number! 
    Float_t dPhimin        = 0.;  // Used for histograms only! 
    Float_t dThetamin      = 0.;  // Used for histograms only! 
    Float_t dZetamin       = 0.;  // Used for histograms only! 

    // in any case, if MC has been required, store statistics of primaries and secondaries
    Bool_t primary=kFALSE; Bool_t secondary=kFALSE; // it is better to have both since chip might not be found
    if (fMC && !lbkg) {
       Int_t lab1=(Int_t)fClustersLay1[iC1][3];
       Int_t lab2=(Int_t)fClustersLay1[iC1][4];
       Int_t lab3=(Int_t)fClustersLay1[iC1][5];
       // do it always as a function of the chip number used to built the prediction
       found=FindChip(key,0,vtx,fClustersLay1[iC1][0],fClustersLay1[iC1][1],fClustersLay1[iC1][2]);
       if (!found) {AliDebug(1,
         Form("Reconstruct MC: cannot find chip on inner layer for cluster %d",iC1)); }
       else {
         if((lab1 != -2  &&  PrimaryTrackChecker(lab1,pStack) ) ||
            (lab2 != -2  &&  PrimaryTrackChecker(lab2,pStack) ) ||
            (lab3 != -2  &&  PrimaryTrackChecker(lab3,pStack))) 
         { // this cluster is from a primary particle
           fClusterPrimary[key]++;
           primary=kTRUE;
           if(fUseOnlySecondaryForPred) continue; // skip this tracklet built with a primary track
         } else { // this cluster is from a secondary particle
            fClusterSecondary[key]++;
            secondary=kTRUE;
            if(fUseOnlyPrimaryForPred) continue; // skip this tracklet built with a secondary track
         }
       }
       // do it as a function of the chip number where you exspect the cluster (i.e. tracklet prediction)
       // (in case the prediction is reliable)
       if( fChipPredOnLay2[iC1]<1200) {
         if((lab1 != -2  &&  PrimaryTrackChecker(lab1,pStack) ) ||
            (lab2 != -2  &&  PrimaryTrackChecker(lab2,pStack) ) ||
            (lab3 != -2  &&  PrimaryTrackChecker(lab3,pStack))) fPredictionPrimary[fChipPredOnLay2[iC1]]++;
         else fPredictionSecondary[fChipPredOnLay2[iC1]]++;
         if((lab1 != -2  &&  IsReconstructableAt(1,iC1,lab1,vtx,pStack,tRef)) ||
            (lab2 != -2  &&  IsReconstructableAt(1,iC1,lab2,vtx,pStack,tRef)) ||
            (lab3 != -2  &&  IsReconstructableAt(1,iC1,lab3,vtx,pStack,tRef))) fRecons[fChipPredOnLay2[iC1]]++;
         else fNonRecons[fChipPredOnLay2[iC1]]++;
       }
    }
    
    // Loop on layer 2 
    for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {      
      
      // The following excludes double associations
      if (!fAssociationFlag[iC2]) {
	
	// find the difference in angles
	Float_t dTheta = fClustersLay2[iC2][0] - fClustersLay1[iC1][0];
	Float_t dPhi   = TMath::Abs(fClustersLay2[iC2][1] - fClustersLay1[iC1][1]);
        // take into account boundary condition
        if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;	

	// find the difference in z (between linear projection from layer 1
	// and the actual point: Dzeta= z1/r1*r2 -z2) 	
	Float_t r2    = fClustersLay2[iC2][2]/TMath::Cos(fClustersLay2[iC2][0]);
        Float_t dZeta = TMath::Cos(fClustersLay1[iC1][0])*r2 - fClustersLay2[iC2][2];

 	if (fHistOn && !lbkg) {
	  fhClustersDPhiAll->Fill(dPhi);    
	  fhClustersDThetaAll->Fill(dTheta);    
	  fhClustersDZetaAll->Fill(dZeta);    
	  fhDPhiVsDThetaAll->Fill(dTheta, dPhi);
	  fhDPhiVsDZetaAll->Fill(dZeta, dPhi);
	}

	// make "elliptical" cut in Phi and Zeta! 
	Float_t d = TMath::Sqrt(dPhi*dPhi/fPhiWindowL2/fPhiWindowL2 + 
                                dZeta*dZeta/fZetaWindowL2/fZetaWindowL2);

	if (d>1) continue;      
	
	//look for the minimum distance: the minimum is in iC2WithBestDist
       	if (TMath::Sqrt(dZeta*dZeta+(r2*dPhi*r2*dPhi)) < distmin ) {
	  distmin=TMath::Sqrt(dZeta*dZeta + (r2*dPhi*r2*dPhi));
	  dPhimin = dPhi;
	  dThetamin = dTheta;
	  dZetamin = dZeta; 
	  iC2WithBestDist = iC2;
	}
      } 
    } // end of loop over clusters in layer 2 
    
    if (distmin<100) { // This means that a cluster in layer 2 was found that matches with iC1

      if (fHistOn && !lbkg) {
	fhClustersDPhiAcc->Fill(dPhimin);
	fhClustersDThetaAcc->Fill(dThetamin);    
	fhClustersDZetaAcc->Fill(dZetamin);    
	fhDPhiVsDThetaAcc->Fill(dThetamin, dPhimin);
	fhDPhiVsDZetaAcc->Fill(dZetamin, dPhimin);
      }
      
      if (fOnlyOneTrackletPerC2) fAssociationFlag[iC2WithBestDist] = kTRUE; 
       // flag the association
      
      // store the tracklet
      
      // use the theta from the clusters in the first layer
      fTracklets[fNTracklets][0] = fClustersLay1[iC1][0];
      // use the phi from the clusters in the first layer
      fTracklets[fNTracklets][1] = fClustersLay1[iC1][1];
      // Store the difference between phi1 and phi2
      fTracklets[fNTracklets][2] = fClustersLay1[iC1][1] - fClustersLay2[iC2WithBestDist][1];

      // find labels
      Int_t label1 = 0;
      Int_t label2 = 0;
      while (label2 < 3)
      {
        if ((Int_t) fClustersLay1[iC1][3+label1] != -2 && (Int_t) fClustersLay1[iC1][3+label1] == (Int_t) fClustersLay2[iC2WithBestDist][3+label2])
          break;
        label1++;
        if (label1 == 3)
        {
          label1 = 0;
          label2++;
        }
      }

      if (label2 < 3)
      {
        fTracklets[fNTracklets][3] = fClustersLay1[iC1][3+label1];
      }
      else
      {
        fTracklets[fNTracklets][3] = -2;
      }

      if (fHistOn && !lbkg) {
	Float_t eta=fTracklets[fNTracklets][0];
	eta= TMath::Tan(eta/2.);
	eta=-TMath::Log(eta);
	fhetaTracklets->Fill(eta);    
	fhphiTracklets->Fill(fTracklets[fNTracklets][1]);    
      }

// Check that this cluster is still in the same chip (here you pass also Zvtx for better computation)
      found=FindChip(key,1,vtx,fClustersLay2[iC2WithBestDist][0],fClustersLay2[iC2WithBestDist][1],fClustersLay2[iC2WithBestDist][2]);
      if(!found){
        AliDebug(1,
         Form("Reconstruct: cannot find chip on outer layer for cluster %d",iC2WithBestDist));
        key=999999;
      }
      nfClu2+=(Int_t)found; // this for debugging purpose
      ntClu2++;             // to check efficiency of the method FindChip
      if(key<1200) { // the Chip has been found
        if(fMC && !lbkg) { // this part only for MC
          // Int_t labc1=(Int_t)fClustersLay2[iC2WithBestDist][3];
          // Int_t labc2=(Int_t)fClustersLay2[iC2WithBestDist][4];
          // Int_t labc3=(Int_t)fClustersLay2[iC2WithBestDist][5];
          if (label2 < 3) {
            fSuccessTT[key]++;
            if(primary) fSuccessPP[key]++;
          }
          if (fUseOnlyDifferentParticle && label2 < 3) continue; // same label (reject it)
          if (fUseOnlySameParticle && label2 == 3) continue;      // different label (reject it)
        }

        if (key==fChipPredOnLay2[iC1]) { // this control seems too loose: has to be checked !
          				 // OK, success
                pe->UpDatePlaneEff(kTRUE,key); // success
                fChipUpdatedInEvent[key]=kTRUE; 
                if(fMC && !lbkg) {
                  if(primary)   fSuccessP[key]++;
                  if(secondary) fSuccessS[key]++;
                }
        }
        else {
                pe->UpDatePlaneEff(kTRUE,key); // this should not be a failure
                fChipUpdatedInEvent[key]=kTRUE;          // (might be in the tracking tollerance)
                if(fMC && !lbkg) {
                  if(primary)   fSuccessP[key]++;
                  if(secondary) fSuccessS[key]++;
                }
        }
      }

      fNTracklets++;

    } // if any cluster found --> increment statistics by 1 failure (provided you have chip prediction)
    else if (fChipPredOnLay2[iC1]<1200) {
      pe->UpDatePlaneEff(kFALSE,fChipPredOnLay2[iC1]);
      fChipUpdatedInEvent[fChipPredOnLay2[iC1]]=kTRUE;
      if(fMC && !lbkg) {
        if(primary)   fFailureP[fChipPredOnLay2[iC1]]++;
        if(secondary) fFailureS[fChipPredOnLay2[iC1]]++;
      }
    }
  } // end of loop over clusters in layer 1

    fNTracklets1=fNTracklets;

//###################################################################

  // Second part : Interpolation to Layer 1 

  // Loop on layer 2 
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {    

    // here the control to check whether the efficiency of the chip traversed by this tracklet
    // prediction has already been updated in this event using another tracklet prediction
    if(fUpdateOncePerEventPlaneEff && fChipPredOnLay1[iC2]<1200 && fChipUpdatedInEvent[fChipPredOnLay1[iC2]]) continue;

    // reset of variables for multiple candidates
    Int_t  iC1WithBestDist = 0;     // reset 
    Float_t distmin        = 100.;  // just to put a huge number! 
    Float_t dPhimin        = 0.;  // Used for histograms only! 
    Float_t dThetamin      = 0.;  // Used for histograms only! 
    Float_t dZetamin       = 0.;  // Used for histograms only! 

    // in any case, if MC has been required, store statistics of primaries and secondaries
    Bool_t primary=kFALSE; Bool_t secondary=kFALSE;
    if (fMC && !lbkg) {
       Int_t lab1=(Int_t)fClustersLay2[iC2][3];
       Int_t lab2=(Int_t)fClustersLay2[iC2][4];
       Int_t lab3=(Int_t)fClustersLay2[iC2][5];
       // do it always as a function of the chip number used to built the prediction
       found=FindChip(key,1,vtx,fClustersLay2[iC2][0],fClustersLay2[iC2][1],fClustersLay2[iC2][2]);
       if (!found) {AliDebug(1,
         Form("Reconstruct MC: cannot find chip on outer layer for cluster %d",iC2)); }
       else {
         if((lab1 != -2  &&  PrimaryTrackChecker(lab1,pStack) ) ||
            (lab2 != -2  &&  PrimaryTrackChecker(lab2,pStack) ) ||
            (lab3 != -2  &&  PrimaryTrackChecker(lab3,pStack))) 
         {  // this cluster is from a primary particle
            fClusterPrimary[key]++;
            primary=kTRUE;
            if(fUseOnlySecondaryForPred) continue; //  skip this tracklet built with a primary track
         } else { // this cluster is from a secondary particle
           fClusterSecondary[key]++;
           secondary=kTRUE;
           if(fUseOnlyPrimaryForPred) continue; //  skip this tracklet built with a secondary track
         }
       }
       // do it as a function of the chip number where you exspect the cluster (i.e. tracklet prediction)
       // (in case the prediction is reliable)
       if( fChipPredOnLay1[iC2]<1200) {
         if((lab1 != -2  &&  PrimaryTrackChecker(lab1,pStack) ) ||
            (lab2 != -2  &&  PrimaryTrackChecker(lab2,pStack) ) ||
            (lab3 != -2  &&  PrimaryTrackChecker(lab3,pStack)))   fPredictionPrimary[fChipPredOnLay1[iC2]]++;
         else fPredictionSecondary[fChipPredOnLay1[iC2]]++;
         if((lab1 != -2  &&  IsReconstructableAt(0,iC2,lab1,vtx,pStack,tRef)) ||
            (lab2 != -2  &&  IsReconstructableAt(0,iC2,lab2,vtx,pStack,tRef)) ||
            (lab3 != -2  &&  IsReconstructableAt(0,iC2,lab3,vtx,pStack,tRef))) fRecons[fChipPredOnLay1[iC2]]++;
         else fNonRecons[fChipPredOnLay1[iC2]]++;
       }
    }
    
    // Loop on layer 1 
    for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {
      
      // The following excludes double associations
      if (!fAssociationFlag1[iC1]) {
	
	// find the difference in angles
	Float_t dTheta = fClustersLay2[iC2][0] - fClustersLay1[iC1][0];
	Float_t dPhi   = TMath::Abs(fClustersLay2[iC2][1] - fClustersLay1[iC1][1]);
        // take into account boundary condition
        if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;	


	// find the difference in z (between linear projection from layer 2
	// and the actual point: Dzeta= z2/r2*r1 -z1) 	
	Float_t r1    = fClustersLay1[iC1][2]/TMath::Cos(fClustersLay1[iC1][0]);
        Float_t dZeta = TMath::Cos(fClustersLay2[iC2][0])*r1 - fClustersLay1[iC1][2];


 	if (fHistOn && !lbkg) {
	  fhClustersDPhiInterpAll->Fill(dPhi);    
	  fhClustersDThetaInterpAll->Fill(dTheta);    
	  fhClustersDZetaInterpAll->Fill(dZeta);    
	  fhDPhiVsDThetaInterpAll->Fill(dTheta, dPhi);
	  fhDPhiVsDZetaInterpAll->Fill(dZeta, dPhi);
	}
	// make "elliptical" cut in Phi and Zeta! 
	Float_t d = TMath::Sqrt(dPhi*dPhi/fPhiWindowL1/fPhiWindowL1 + 
                                dZeta*dZeta/fZetaWindowL1/fZetaWindowL1);

	if (d>1) continue;      
	
	//look for the minimum distance: the minimum is in iC1WithBestDist
       	if (TMath::Sqrt(dZeta*dZeta+(r1*dPhi*r1*dPhi)) < distmin ) {
	  distmin=TMath::Sqrt(dZeta*dZeta + (r1*dPhi*r1*dPhi));
	  dPhimin = dPhi;
	  dThetamin = dTheta;
	  dZetamin = dZeta; 
	  iC1WithBestDist = iC1;
	}
      } 
    } // end of loop over clusters in layer 1 
    
    if (distmin<100) { // This means that a cluster in layer 1 was found that matches with iC2

      if (fHistOn && !lbkg) {
	fhClustersDPhiInterpAcc->Fill(dPhimin);
	fhClustersDThetaInterpAcc->Fill(dThetamin);    
	fhClustersDZetaInterpAcc->Fill(dZetamin);    
	fhDPhiVsDThetaInterpAcc->Fill(dThetamin, dPhimin);
	fhDPhiVsDZetaInterpAcc->Fill(dZetamin, dPhimin);
      }
      
      if (fOnlyOneTrackletPerC1) fAssociationFlag1[iC1WithBestDist] = kTRUE; // flag the association
       // flag the association
      
      // store the tracklet
      
      // use the theta from the clusters in the first layer
      fTracklets[fNTracklets][0] = fClustersLay2[iC2][0];
      // use the phi from the clusters in the first layer
      fTracklets[fNTracklets][1] = fClustersLay2[iC2][1];
      // Store the difference between phi1 and phi2
      fTracklets[fNTracklets][2] = fClustersLay2[iC2][1] - fClustersLay1[iC1WithBestDist][1];

      // find labels
      Int_t label1 = 0;
      Int_t label2 = 0;
      while (label2 < 3)
      {
        if ((Int_t) fClustersLay2[iC2][3+label1] != -2 && (Int_t) fClustersLay2[iC2][3+label1] == (Int_t) fClustersLay1[iC1WithBestDist][3+label2])
          break;
        label1++;
        if (label1 == 3)
        {
          label1 = 0;
          label2++;
        }
      }

      if (label2 < 3)
      {
        fTracklets[fNTracklets][3] = fClustersLay2[iC2][3+label1];
      }
      else
      {
        fTracklets[fNTracklets][3] = -2;
      }

// Check that this cluster is still in the same chip (here you pass also Zvtx for better computation)
      found=FindChip(key,0,vtx,fClustersLay1[iC1WithBestDist][0],fClustersLay1[iC1WithBestDist][1],fClustersLay1[iC1WithBestDist][2]);
      if(!found){
        AliDebug(1,
         Form("Reconstruct: cannot find chip on inner layer for cluster %d",iC1WithBestDist));
        key=999999;
      }
      nfClu1+=(Int_t)found; // this for debugging purpose
      ntClu1++;             // to check efficiency of the method FindChip
      if(key<1200) {
        if(fMC && !lbkg) { // this part only for MC
          // Int_t labc1=(Int_t)fClustersLay1[iC1WithBestDist][3];
          // Int_t labc2=(Int_t)fClustersLay1[iC1WithBestDist][4];
          // Int_t labc3=(Int_t)fClustersLay1[iC1WithBestDist][5];
          if (label2 < 3) { // same label 
            fSuccessTT[key]++;
            if(primary) fSuccessPP[key]++;
          }
          if (fUseOnlyDifferentParticle && label2 < 3) continue; // same label (reject it)
          if (fUseOnlySameParticle && label2 == 3) continue;      // different label (reject it)
        }

        if (key==fChipPredOnLay1[iC2]) { // this control seems too loose: has to be checked !
          				 // OK, success
                pe->UpDatePlaneEff(kTRUE,key); // success
                fChipUpdatedInEvent[key]=kTRUE;
                if(fMC && !lbkg) {
                  if(primary)   fSuccessP[key]++;
                  if(secondary) fSuccessS[key]++;
                }
        } else {
                pe->UpDatePlaneEff(kTRUE,key); // this should not be a failure
                fChipUpdatedInEvent[key]=kTRUE;          // (might be in the tracking tollerance)
                if(fMC && !lbkg) {
                  if(primary)   fSuccessP[key]++;
                  if(secondary) fSuccessS[key]++;
                }
        }
      }

    fNTracklets++;

    } // if no cluster found --> increment statistics by 1 failure (provided you have chip prediction)
    else if (fChipPredOnLay1[iC2]<1200) {
      pe->UpDatePlaneEff(kFALSE,fChipPredOnLay1[iC2]);
      fChipUpdatedInEvent[fChipPredOnLay1[iC2]]=kTRUE;
      if(fMC && !lbkg) {
        if(primary)   fFailureP[fChipPredOnLay1[iC2]]++;
        if(secondary) fFailureS[fChipPredOnLay1[iC2]]++;
      }
    }
  } // end of loop over clusters in layer 2
  
  AliDebug(1,Form("%d tracklets found", fNTracklets));
  AliDebug(1,Form(("Eff. of method FindChip for Track pred. on lay 1 = %d / %d"),nfTraPred1,ntTraPred1));
  AliDebug(1,Form(("Eff. of method FindChip for Track pred. on lay 2 = %d / %d"),nfTraPred2,ntTraPred2));
  AliDebug(1,Form(("Eff. of method FindChip for Cluster on lay 1 = %d / %d"),nfClu1,ntClu1));
  AliDebug(1,Form(("Eff. of method FindChip for Cluster on lay 2 = %d / %d"),nfClu2,ntClu2));
}
//____________________________________________________________________
Bool_t AliITSTrackleterSPDEff::FindChip(UInt_t &key, Int_t layer,const  Float_t* vtx, 
                                  Float_t thetaVtx, Float_t phiVtx, Float_t zVtx) {
//
// Input: a) layer number in the range [0,1]
//        b) vtx[3]: actual vertex 
//        c) zVtx     \ z of the cluster (-999 for tracklet) computed with respect to vtx
//        d) thetaVtx  > theta and phi of the cluster/tracklet computed with respect to vtx
//        e) phiVtx   /
// Output: Unique key to locate a chip
// return: kTRUE if succesfull

    if(layer<0 || layer >1) {AliWarning("Wrong layer: should be 0 or 1!"); return kFALSE;}
    Double_t r=GetRLayer(layer);
    //AliInfo(Form("Radius on layer %d  is %f cm",layer,r));

  // set phiVtx in the range [0,2pi]
  if(!SetAngleRange02Pi(phiVtx)) return kFALSE ;
  
  Double_t zAbs,phiAbs; // those are the polar coordinate, in the Absolute ALICE Reference 
                        // of the intersection of the tracklet with the pixel layer.  
  if (TMath::Abs(zVtx)<100) zAbs=zVtx + vtx[2]; // this is fine only for the cluster, not for the track prediction
  else {
    if(TMath::Abs(thetaVtx)<1E-6) return kFALSE;
    zAbs=r/TMath::Tan(thetaVtx) + vtx[2]; // this is the only way to do for the tracklet prediction
  }
  AliDebug(1,Form("FindChip: vtx[0] = %f, vtx[1] = %f, vtx[2] = %f",vtx[0],vtx[1],vtx[2]));
  Double_t vtxy[2]={vtx[0],vtx[1]};
  if (vtxy[0]*vtxy[1]+vtxy[1]*vtxy[1]>0) { // this method holds only for displaced vertices 
    // this method gives you two interceptions
    if (!FindIntersectionPolar(vtxy,(Double_t)phiVtx,r,phiAbs)) return kFALSE;
    // set phiAbs in the range [0,2pi]
    if(!SetAngleRange02Pi(phiAbs)) return kFALSE; 
    // since Vtx is very close to the ALICE origin, then phiVtx and phiAbs are very close; 
    // therofore you can select the right intersection (among phiAbs1 and phiAbs2) by 
    // taking the closest one to phiVtx
    AliDebug(1,Form("PhiVtx= %f, PhiAbs= %f",phiVtx,phiAbs));
  } else phiAbs=phiVtx;
  Int_t idet=FindDetectorIndex(layer,phiAbs,zAbs); // this is the detector number 

  // now you need to locate the chip within the idet detector, 
  // starting from the local coordinates in such a detector

  Float_t locx; // local Cartesian coordinate (to be determined) corresponding to 
  Float_t locz; // the Global Cilindrica coordinate (r,phiAbs,zAbs) . 
  if(!FromGloCilToLocCart(layer,idet,r,phiAbs,zAbs, locx, locz)) return kFALSE; 

  key=fPlaneEffSPD->GetKeyFromDetLocCoord(layer,idet,locx,locz); 
  return kTRUE;
}
//______________________________________________________________________________
Double_t AliITSTrackleterSPDEff::GetRLayer(Int_t layer) {
//
//  Return the average radius of a layer from Geometry
//
    if(layer<0 || layer >1) {AliWarning("Wrong layer: should be 0 or 1!"); return -999.;}
    Int_t i=layer+1; // in AliITSgeomTGeo you count from 1 to 6 !

    Double_t xyz[3], &x=xyz[0], &y=xyz[1];
    AliITSgeomTGeo::GetOrigTranslation(i,1,1,xyz);
    Double_t r=TMath::Sqrt(x*x + y*y);

    AliITSgeomTGeo::GetOrigTranslation(i,1,2,xyz);
    r += TMath::Sqrt(x*x + y*y);
    AliITSgeomTGeo::GetOrigTranslation(i,2,1,xyz);
    r += TMath::Sqrt(x*x + y*y);
    AliITSgeomTGeo::GetOrigTranslation(i,2,2,xyz);
    r += TMath::Sqrt(x*x + y*y);
    r*=0.25;
    return r;
}
//______________________________________________________________________________
Bool_t AliITSTrackleterSPDEff::FromGloCilToLocCart(Int_t ilayer,Int_t idet, Double_t r, Double_t phi, Double_t z, 
                           Float_t &xloc, Float_t &zloc) {
  // this method transform Global Cilindrical coordinates into local (i.e. module) 
  // cartesian coordinates
  //
  //Compute Cartesian Global Coordinate
  Double_t xyzGlob[3],xyzLoc[3];
  xyzGlob[2]=z;
  xyzGlob[0]=r*TMath::Cos(phi);
  xyzGlob[1]=r*TMath::Sin(phi);

  xloc=0.;
  zloc=0.;

  if(idet<0)  return kFALSE;

  Int_t ndet=AliITSgeomTGeo::GetNDetectors(ilayer+1); // layers from 1 to 6
  Int_t lad = Int_t(idet/ndet) + 1;
  Int_t det = idet - (lad-1)*ndet + 1;

  AliITSgeomTGeo::GlobalToLocal(ilayer+1,lad,det,xyzGlob,xyzLoc);

  xloc = (Float_t)xyzLoc[0];
  zloc = (Float_t)xyzLoc[2];

return kTRUE;
}
//______________________________________________________________________________
Int_t AliITSTrackleterSPDEff::FindDetectorIndex(Int_t layer, Double_t phi, Double_t z) {
  //--------------------------------------------------------------------
  // This function finds the detector crossed by the track
  // Input: layer in range [0,1]
  //        phi   in ALICE absolute reference system
  //         z     "  "       "         "        "
  //--------------------------------------------------------------------
    if(layer<0 || layer >1) {AliWarning("Wrong layer: should be 0 or 1!"); return -1;}
    Int_t i=layer+1; // in AliITSgeomTGeo you count from 1 to 6 !
    Int_t nladders=AliITSgeomTGeo::GetNLadders(i);
    Int_t ndetectors=AliITSgeomTGeo::GetNDetectors(i);

    Double_t xyz[3], &x=xyz[0], &y=xyz[1], &z2=xyz[2];
    AliITSgeomTGeo::GetOrigTranslation(i,1,1,xyz);
    Double_t phiOffset=TMath::ATan2(y,x);
    Double_t zOffset=z2;

  Double_t dphi;
  if (zOffset<0)            // old geometry
    dphi = -(phi-phiOffset);
  else                       // new geometry
    dphi = phi-phiOffset;

  if      (dphi <  0) dphi += 2*TMath::Pi();
  else if (dphi >= 2*TMath::Pi()) dphi -= 2*TMath::Pi();
  Int_t np=Int_t(dphi*nladders*0.5/TMath::Pi()+0.5);
  if (np>=nladders) np-=nladders;
  if (np<0)          np+=nladders;

  Double_t dz=zOffset-z;
  Double_t nnz = dz*(ndetectors-1)*0.5/zOffset+0.5;
  Int_t nz = (nnz<0 ? -1 : (Int_t)nnz);
  if (nz>=ndetectors) {AliDebug(1,Form("too  long: nz =%d",nz)); return -1;}
  if (nz<0)           {AliDebug(1,Form("too short: nz =%d",nz)); return -1;}

  return np*ndetectors + nz;
}
//____________________________________________________________
Bool_t AliITSTrackleterSPDEff::FindIntersectionPolar(Double_t vtx[2],Double_t phiVtx, Double_t R,Double_t &phi) {
// this method find the intersection in xy between a tracklet (straight line) and 
// a circonference (r=R), using polar coordinates. 
/*
Input: - vtx[2]: actual vertex w.r.t. ALICE reference system
       - phiVtx: phi angle of the line (tracklet) computed w.r.t. vtx
       - R: radius of the circle
Output: - phi : phi angle of the unique interception in the ALICE Global ref. system 

Correct method below: you have the equation of a circle (in polar coordinate) w.r.t. Actual vtx:
r^2-2*r*r0*cos(phi-phi0) + r0^2 = R^2 , where (r0,phi0) is the centre of the circle
In the same system, the equation of a semi-line is: phi=phiVtx;
Hence you get one interception only: P=(r,phiVtx)
Finally you want P in the ABSOLUTE ALICE reference system.
*/
Double_t rO=TMath::Sqrt(vtx[0]*vtx[0]+vtx[1]*vtx[1]); // polar coordinates of the ALICE origin
Double_t phiO=TMath::ATan2(-vtx[1],-vtx[0]);          // in the system with vtx[2] as origin
Double_t bB=-2.*rO*TMath::Cos(phiVtx-phiO);
Double_t cC=rO*rO-R*R;
Double_t dDelta=bB*bB-4*cC;
if(dDelta<0) return kFALSE;
Double_t r1,r2;
r1=(-bB-TMath::Sqrt(dDelta))/2;
r2=(-bB+TMath::Sqrt(dDelta))/2;
if(r1*r2>0) { printf("allora non hai capito nulla \n"); return kFALSE;}
Double_t r=TMath::Max(r1,r2); // take the positive
Double_t pvtx[2]; // Cartesian coordinates of the interception w.r.t. vtx
Double_t pP[2]; // Cartesian coordinates of the interception w.r.t. ALICE origin
pvtx[0]=r*TMath::Cos(phiVtx);
pvtx[1]=r*TMath::Sin(phiVtx);
pP[0]=vtx[0]+pvtx[0];
pP[1]=vtx[1]+pvtx[1];
phi=TMath::ATan2(pP[1],pP[0]);
return kTRUE;
}
//___________________________________________________________
Bool_t AliITSTrackleterSPDEff::SetAngleRange02Pi(Double_t &angle) const {
//
//  simple method to reduce all angles (in rad)
//  in range [0,2pi[
//
//
while(angle >=2*TMath::Pi() || angle<0) {
  if(angle >= 2*TMath::Pi()) angle-=2*TMath::Pi();
  if(angle < 0) angle+=2*TMath::Pi();
}
return kTRUE;
}
//___________________________________________________________
Bool_t AliITSTrackleterSPDEff::PrimaryTrackChecker(Int_t ipart,AliStack* stack) {
//
//  This method check if a particle is primary; i.e.  
//  it comes from the main vertex and it is a "stable" particle, according to 
//  AliStack::IsPhysicalPrimary() (note that there also Sigma0 are considered as 
//  a stable particle: it has no effect on this analysis). 
//  This method can be called only for MC events, where Kinematics is available.
//  if fUseOnlyStableParticle is kTRUE (via SetUseOnlyStableParticle) then it 
//  returns kTRUE if also AliITSTrackleterSPDEff::DecayingTrackChecker() return 0.
//  The latter (see below) try to verify if a primary particle is also "detectable".
//
if(!fMC) {AliError("This method works only if SetMC() has been called"); return kFALSE;}
if(!stack) {AliError("null pointer to MC stack"); return kFALSE;}
if(ipart >= stack->GetNtrack()) {AliError("this track label is not in MC stack"); return kFALSE;}
// return stack->IsPhysicalPrimary(ipart); // looking at AliStack.cxx this does not seem to be complete (e.g. Pi0 Dalitz)
 if(!stack->IsPhysicalPrimary(ipart)) return kFALSE; 
 // like below: as in the correction for Multiplicity (i.e. by hand in macro)
 TParticle* part = stack->Particle(ipart);
 TParticle* part0 = stack->Particle(0); // first primary
 TParticle* partl = stack->Particle(stack->GetNprimary()-1); //last primary
 if (part0->Vx()-partl->Vx()>0) AliDebug(1,Form("Difference in vtx position between 1th and last primaries %f %f %f",
            part0->Vx()-partl->Vx(),part0->Vy()-partl->Vy(), part0->Vz()-partl->Vz() ));

 if (!part || strcmp(part->GetName(),"XXX")==0) {AliWarning("String , not particle ??") ;return kFALSE; }
 TParticlePDG* pdgPart = part->GetPDG();
 if (TMath::Abs(pdgPart->Charge()) < 3) {AliWarning("This seems a quark"); return kFALSE;}
 
  Double_t distx = part->Vx() - part0->Vx();
  Double_t disty = part->Vy() - part0->Vy();
  Double_t distz = part->Vz() - part0->Vz();
  Double_t distR=TMath::Sqrt(distx*distx + disty*disty + distz*distz);

  if (distR > 0.05) {AliDebug(1,Form("True vertex should be %f %f, this particle from %f %f ",
                                 part0->Vx(),part0->Vy(),part->Vx(),part->Vy()));
                      return kFALSE; }// primary if within 500 microns from true Vertex

 if(fUseOnlyStableParticle && DecayingTrackChecker(ipart,stack)>0) return kFALSE; 
 return kTRUE;
}
//_____________________________________________________________________________________________
Int_t AliITSTrackleterSPDEff::DecayingTrackChecker(Int_t ipart,AliStack* stack) {
//
// This private method can be applied on MC particles (if stack is available),  
// provided they have been identified as "primary" from PrimaryTrackChecker() (see above).
//   
// It define "detectable" a primary particle according to the following criteria:
//
// - if no decay products can be found in the stack (note that this does not 
//     means it is stable, since a particle is stored in stack if it has at least 1 hit in a 
//     sensitive detector)
// - if it has at least one decay daughter produced outside or just on the outer pixel layer 
// - if the last decay particle is an electron (or a muon) which is not produced in-between 
//     the two pixel layers (this is likely to be a kink).
if(!fMC) {AliError("This method works only if SetMC() has been called"); return 0;}
if(!stack) {AliError("null pointer to MC stack"); return 0;}
if(ipart >= stack->GetNtrack()) {AliError("this track label is not in MC stack"); return 0;}

TParticle* part = stack->Particle(ipart);
//TParticle* part0 = stack->Particle(0); // first primary

  Int_t nret=0;
  TParticle* dau = 0;
  Int_t nDau = 0;
  Int_t pdgDau;
  Int_t firstDau = part->GetFirstDaughter(); // if no daugther stored then no way to understand i
                                             // its real fate ! But you have to take it !
  if (firstDau > 0) { // if it has daugther(s) try to infer if it is "detectable" as a tracklet
    Int_t lastDau = part->GetLastDaughter();
    nDau = lastDau - firstDau + 1;
    Double_t distMax=0.;
    Int_t jmax=0;
    for(Int_t j=firstDau; j<=lastDau; j++)  {
      dau = stack->Particle(j);
      Double_t distx = dau->Vx();
      Double_t disty = dau->Vy();
      //Double_t distz = dau->Vz();
      Double_t distR = TMath::Sqrt(distx*distx+disty*disty);
      if(distR<distMax) continue; // considere only the daughter produced at largest radius
      distMax=distR;
      jmax=j;
    }
    dau = stack->Particle(jmax);
    pdgDau=dau->GetPdgCode();
    if (pdgDau == 11 || pdgDau == 13 ) {
       if(distMax < GetRLayer(1)-0.25 && distMax > GetRLayer(0)+0.27) nret=1; // can be a kink (reject it)
       else nret =0; // delta-ray emission in material  (keep it)
    }
    else {// not ele or muon
      if (distMax < GetRLayer(1)-0.25 )  nret= 1;}  // decay before the second pixel layer (reject it)
    }
return nret;
}
//_________________________________________________________________
void AliITSTrackleterSPDEff::InitPredictionMC() {
//
// this method allocate memory for the MC related informations
// all the counters are set to 0
//
//
if(!fMC) {AliError("This method works only if SetMC() has been called"); return;}
fPredictionPrimary   = new Int_t[1200];
fPredictionSecondary = new Int_t[1200];
fClusterPrimary      = new Int_t[1200];
fClusterSecondary    = new Int_t[1200];
fSuccessPP           = new Int_t[1200];
fSuccessTT           = new Int_t[1200];
fSuccessS            = new Int_t[1200];
fSuccessP            = new Int_t[1200];
fFailureS            = new Int_t[1200];
fFailureP            = new Int_t[1200];
fRecons              = new Int_t[1200];
fNonRecons           = new Int_t[1200];
for(Int_t i=0; i<1200; i++) {
 fPredictionPrimary[i]=0;
 fPredictionSecondary[i]=0; 
 fPredictionSecondary[i]=0;
 fClusterSecondary[i]=0;
 fSuccessPP[i]=0;
 fSuccessTT[i]=0;
 fSuccessS[i]=0;
 fSuccessP[i]=0;
 fFailureS[i]=0;
 fFailureP[i]=0;
 fRecons[i]=0;
 fNonRecons[i]=0;
}
return;
}
//_________________________________________________________________
void AliITSTrackleterSPDEff::DeletePredictionMC() {
//
// this method deallocate memory for the MC related informations
// all the counters are set to 0
//
//
if(fMC) {AliInfo("This method works only if fMC=kTRUE"); return;}
if(fPredictionPrimary) {
  delete fPredictionPrimary; fPredictionPrimary=0;
}
if(fPredictionSecondary) {
  delete fPredictionSecondary; fPredictionSecondary=0;
}
if(fClusterPrimary) {
  delete fClusterPrimary; fClusterPrimary=0;
}
if(fClusterSecondary) {
  delete fClusterSecondary; fClusterSecondary=0;
}
if(fSuccessPP) {
  delete fSuccessPP; fSuccessPP=0;
}
if(fSuccessTT) {
  delete fSuccessTT; fSuccessTT=0;
}
if(fSuccessS) {
  delete fSuccessS; fSuccessS=0;
}
if(fSuccessP) {
  delete fSuccessP; fSuccessP=0;
}
if(fFailureS) {
  delete fFailureS; fFailureS=0;
}
if(fFailureP) {
  delete fFailureP; fFailureP=0;
}
if(fRecons) {
  delete fRecons; fRecons=0;
}
if(fNonRecons) {
  delete fNonRecons; fNonRecons=0;
}
return;
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetPredictionPrimary(const UInt_t key) const {
//
// This method return the Data menmber fPredictionPrimary [1200].
// You can call it only for MC events.
// fPredictionPrimary[key] contains the number of tracklet predictions on the
// given chip key built using  a cluster on the other layer produced (at least)
// from a primary particle.
// Key refers to the chip crossed by the prediction 
//
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fPredictionPrimary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetPredictionSecondary(const UInt_t key) const {
//
// This method return the Data menmber fPredictionSecondary [1200].
// You can call it only for MC events.
// fPredictionSecondary[key] contains the number of tracklet predictions on the
// given chip key built using  a cluster on the other layer produced (only)
// from a secondary particle
// Key refers to the chip crossed by the prediction 
//
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fPredictionSecondary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetClusterPrimary(const UInt_t key) const {
//
// This method return the Data menmber fClusterPrimary [1200].
// You can call it only for MC events.
// fClusterPrimary[key] contains the number of tracklet predictions 
// built using  a cluster on that layer produced (only)
// from a primary particle
// Key refers to the chip used to build the prediction
//
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fClusterPrimary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetClusterSecondary(const UInt_t key) const {
//
// This method return the Data menmber fClusterSecondary [1200].
// You can call it only for MC events.
// fClusterSecondary[key] contains the number of tracklet predictions
// built using  a cluster on that layer produced (only)
// from a secondary particle
// Key refers to the chip used to build the prediction
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fClusterSecondary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetSuccessPP(const UInt_t key) const {
//
// This method return the Data menmber fSuccessPP [1200].
// You can call it only for MC events.
// fSuccessPP[key] contains the number of successes (i.e. a tracklet prediction matching
// with a cluster on the other layer) built by using the same primary particle
// the unique chip key refers to the chip which get updated its efficiency
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fSuccessPP[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetSuccessTT(const UInt_t key) const {
//
// This method return the Data menmber fSuccessTT [1200].
// You can call it only for MC events.
// fSuccessTT[key] contains the number of successes (i.e. a tracklet prediction matching
// with a cluster on the other layer) built by using the same  particle (whatever)
// the unique chip key refers to the chip which get updated its efficiency
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fSuccessTT[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetSuccessS(const UInt_t key) const {
//
// This method return the Data menmber fSuccessS [1200].
// You can call it only for MC events.
// fSuccessS[key] contains the number of successes (i.e. a tracklet prediction matching
// with a cluster on the other layer) built by using a secondary particle
// the unique chip key refers to the chip which get updated its efficiency
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fSuccessS[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetSuccessP(const UInt_t key) const {
//
// This method return the Data menmber fSuccessP [1200].
// You can call it only for MC events.
// fSuccessP[key] contains the number of successes (i.e. a tracklet prediction matching
// with a cluster on the other layer) built by using a primary particle
// the unique chip key refers to the chip which get updated its efficiency
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fSuccessP[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetFailureS(const UInt_t key) const {
//
// This method return the Data menmber fFailureS [1200].
// You can call it only for MC events.
// fFailureS[key] contains the number of failures (i.e. a tracklet prediction not matching
// with a cluster on the other layer) built by using a secondary particle
// the unique chip key refers to the chip which get updated its efficiency
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fFailureS[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetFailureP(const UInt_t key) const {
//
// This method return the Data menmber fFailureP [1200].
// You can call it only for MC events.
// fFailureP[key] contains the number of failures (i.e. a tracklet prediction not matching
// with a cluster on the other layer) built by using a primary particle
// the unique chip key refers to the chip which get updated its efficiency
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fFailureP[(Int_t)key];
}
//_____________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetRecons(const UInt_t key) const {
//
// This method return the Data menmber fRecons [1200].
// You can call it only for MC events.
// fRecons[key] contains the number of reconstractable tracklets (i.e. a tracklet prediction which
// has an hit in the detector)
// the unique chip key refers to the chip where fall the prediction
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fRecons[(Int_t)key];
}
//_____________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetNonRecons(const UInt_t key) const {
//
// This method return the Data menmber fNonRecons [1200].
// You can call it only for MC events.
// fRecons[key] contains the number of unreconstractable tracklets (i.e. a tracklet prediction which
// has not any hit in the detector)
// the unique chip key refers to the chip where fall the prediction
//
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fNonRecons[(Int_t)key];
}
//______________________________________________________________________
void AliITSTrackleterSPDEff::PrintAscii(ostream *os)const{
    // Print out some class data values in Ascii Form to output stream
    // Inputs:
    //   ostream *os   Output stream where Ascii data is to be writen
    // Outputs:
    //   none.
    // Return:
    //   none.
    *os << fPhiWindowL1 <<" "<< fZetaWindowL1 << " " << fPhiWindowL2 <<" "<< fZetaWindowL2 
        << " " << fOnlyOneTrackletPerC1 << " " << fOnlyOneTrackletPerC2 
        << " " << fUpdateOncePerEventPlaneEff << " " << fMinContVtx 
        << " " << fReflectClusterAroundZAxisForLayer0
        << " " << fReflectClusterAroundZAxisForLayer1;
    *os << " " << fMC;
    if(!fMC) {AliInfo("Writing only cuts, no MC info"); return;}
    *os << " " << fUseOnlyPrimaryForPred << " " << fUseOnlySecondaryForPred
        << " " << fUseOnlySameParticle   << " " << fUseOnlyDifferentParticle
        << " " << fUseOnlyStableParticle ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetPredictionPrimary(i)  ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetPredictionSecondary(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetClusterPrimary(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetClusterSecondary(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetSuccessPP(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetSuccessTT(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetSuccessS(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetSuccessP(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetFailureS(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetFailureP(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetRecons(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetNonRecons(i) ;
    return;
}
//______________________________________________________________________
void AliITSTrackleterSPDEff::ReadAscii(istream *is){
    // Read in some class data values in Ascii Form to output stream
    // Inputs:
    //   istream *is   Input stream where Ascii data is to be read in from
    // Outputs:
    //   none.
    // Return:
    //   none.

    Bool_t tmp= fMC;
    *is >> fPhiWindowL1 >> fZetaWindowL1 >> fPhiWindowL2 >> fZetaWindowL2 
        >> fOnlyOneTrackletPerC1 >> fOnlyOneTrackletPerC2  
        >> fUpdateOncePerEventPlaneEff >> fMinContVtx 
        >> fReflectClusterAroundZAxisForLayer0
        >> fReflectClusterAroundZAxisForLayer1;
    //if(!fMC) {AliInfo("Reading only cuts, no MC info available");return;}
    *is >> fMC;
    if(!fMC) {AliInfo("Reading only cuts, no MC info"); if(tmp) SetMC(kFALSE); }
    else {
      if(!tmp) {AliInfo("Calling SetMC() to read this file wtih MC info"); SetMC();}
      *is >> fUseOnlyPrimaryForPred >> fUseOnlySecondaryForPred
          >> fUseOnlySameParticle   >> fUseOnlyDifferentParticle
          >> fUseOnlyStableParticle;
      for(Int_t i=0;i<1200;i++) *is >> fPredictionPrimary[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fPredictionSecondary[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fClusterPrimary[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fClusterSecondary[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fSuccessPP[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fSuccessTT[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fSuccessS[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fSuccessP[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fFailureS[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fFailureP[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fRecons[i] ;
      for(Int_t i=0;i<1200;i++) *is >> fNonRecons[i] ;
    } 
    return;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,const AliITSTrackleterSPDEff &s){
    // Standard output streaming function
    // Inputs:
    //   ostream            &os  output steam
    //   AliITSTrackleterSPDEff &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer

    s.PrintAscii(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTrackleterSPDEff &s){
    // Standard inputput streaming function
    // Inputs:
    //   istream            &is  input steam
    //   AliITSTrackleterSPDEff &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer

    //printf("prova %d \n", (Int_t)s.GetMC());
    s.ReadAscii(&is);
    return is;
}
//______________________________________________________________________
void AliITSTrackleterSPDEff::SavePredictionMC(TString filename) const {
//
// This Method write into an either asci or root file 
// the used cuts and the statistics  of the MC related quantities
// The method SetMC() has to be called before 
// Input TString filename: name of file for output (it deletes already existing 
// file)
// Output: none
//
//
 //if(!fMC) {CallWarningMC(); return;}
 if (!filename.Contains(".root")) {
   ofstream out(filename.Data(),ios::out | ios::binary);
   out << *this;
   out.close();
   return;
 }
 else {
    TFile* mcfile = TFile::Open(filename, "RECREATE");
    TH1F* cuts = new TH1F("cuts", "list of cuts", 11, 0, 11); // TH1I containing cuts 
    cuts->SetBinContent(1,fPhiWindowL1);
    cuts->SetBinContent(2,fZetaWindowL1);
    cuts->SetBinContent(3,fPhiWindowL2);
    cuts->SetBinContent(4,fZetaWindowL2);
    cuts->SetBinContent(5,fOnlyOneTrackletPerC1);
    cuts->SetBinContent(6,fOnlyOneTrackletPerC2);
    cuts->SetBinContent(7,fUpdateOncePerEventPlaneEff);
    cuts->SetBinContent(8,fMinContVtx);
    cuts->SetBinContent(9,fReflectClusterAroundZAxisForLayer0);
    cuts->SetBinContent(10,fReflectClusterAroundZAxisForLayer1);
    cuts->SetBinContent(11,fMC);
    cuts->Write();
    delete cuts;
    if(!fMC) {AliInfo("Writing only cuts, no MC info");}
    else {
      TH1C* mc0 = new TH1C("mc0", "mc cuts", 5, 0, 5);
      mc0->SetBinContent(1,fUseOnlyPrimaryForPred);
      mc0->SetBinContent(2,fUseOnlySecondaryForPred);
      mc0->SetBinContent(3,fUseOnlySameParticle);
      mc0->SetBinContent(4,fUseOnlyDifferentParticle);
      mc0->SetBinContent(5,fUseOnlyStableParticle);
      mc0->Write();
      delete mc0;
      TH1I *mc1;
      mc1 = new TH1I("mc1", "mc info PredictionPrimary", 1200, 0, 1200); 
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetPredictionPrimary(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc2", "mc info PredictionSecondary", 1200, 0, 1200); 
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetPredictionSecondary(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc3", "mc info ClusterPrimary", 1200, 0, 1200); 
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetClusterPrimary(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc4", "mc info ClusterSecondary", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetClusterSecondary(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc5", "mc info SuccessPP", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetSuccessPP(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc6", "mc info SuccessTT", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetSuccessTT(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc7", "mc info SuccessS", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetSuccessS(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc8", "mc info SuccessP", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetSuccessP(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc9", "mc info FailureS", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetFailureS(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc10", "mc info FailureP", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetFailureP(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc11", "mc info Recons", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetRecons(i)) ;
      mc1->Write();
      mc1 = new TH1I("mc12", "mc info NonRecons", 1200, 0, 1200);
      for(Int_t i=0;i<1200;i++)  mc1->SetBinContent(i+1,GetNonRecons(i)) ;
      mc1->Write();
      delete mc1;
   }
   mcfile->Close();
 }
return;
}
//____________________________________________________________________
void AliITSTrackleterSPDEff::ReadPredictionMC(TString filename) {
//
// This Method read from an asci file (do not know why binary does not work)
// the cuts to be used and the statistics  of the MC related quantities
// Input TString filename: name of input file for output 
// The method SetMC() has to be called before
// Output: none
//
//
 //if(!fMC) {CallWarningMC(); return;}
 if( gSystem->AccessPathName( filename.Data() ) ) {
      AliError( Form( "file (%s) not found", filename.Data() ) );
      return;
   }

 if (!filename.Contains(".root")) {
   ifstream in(filename.Data(),ios::in | ios::binary);
   in >> *this;
   in.close();
   return;
 }
 else {
    Bool_t tmp= fMC;
    TFile *mcfile = TFile::Open(filename);
    TH1F *cuts = (TH1F*)mcfile->Get("cuts"); 
    fPhiWindowL1=(Float_t)cuts->GetBinContent(1);
    fZetaWindowL1=(Float_t)cuts->GetBinContent(2);
    fPhiWindowL2=(Float_t)cuts->GetBinContent(3);
    fZetaWindowL2=(Float_t)cuts->GetBinContent(4);
    fOnlyOneTrackletPerC1=(Bool_t)cuts->GetBinContent(5);
    fOnlyOneTrackletPerC2=(Bool_t)cuts->GetBinContent(6);
    fUpdateOncePerEventPlaneEff=(Bool_t)cuts->GetBinContent(7);
    fMinContVtx=(Int_t)cuts->GetBinContent(8);
    fReflectClusterAroundZAxisForLayer0=(Bool_t)cuts->GetBinContent(9);
    fReflectClusterAroundZAxisForLayer1=(Bool_t)cuts->GetBinContent(10);
    fMC=(Bool_t)cuts->GetBinContent(11);
    if(!fMC) {AliInfo("Reading only cuts, no MC info"); if(tmp) SetMC(kFALSE); }
    else { // only if file with MC predictions 
      if(!tmp) {AliInfo("Calling SetMC() to read this file wtih MC info"); SetMC();}
      TH1C *mc0 = (TH1C*)mcfile->Get("mc0");
      fUseOnlyPrimaryForPred=(Bool_t)mc0->GetBinContent(1);
      fUseOnlySecondaryForPred=(Bool_t)mc0->GetBinContent(2);
      fUseOnlySameParticle=(Bool_t)mc0->GetBinContent(3);
      fUseOnlyDifferentParticle=(Bool_t)mc0->GetBinContent(4);
      fUseOnlyStableParticle=(Bool_t)mc0->GetBinContent(5);
      TH1I *mc1;
      mc1 =(TH1I*)mcfile->Get("mc1");
      for(Int_t i=0;i<1200;i++)  fPredictionPrimary[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc2");
      for(Int_t i=0;i<1200;i++)  fPredictionSecondary[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc3");
      for(Int_t i=0;i<1200;i++)  fClusterPrimary[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc4");
      for(Int_t i=0;i<1200;i++)  fClusterSecondary[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc5");
      for(Int_t i=0;i<1200;i++)  fSuccessPP[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc6");
      for(Int_t i=0;i<1200;i++)  fSuccessTT[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc7");
      for(Int_t i=0;i<1200;i++)  fSuccessS[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc8");
      for(Int_t i=0;i<1200;i++)  fSuccessP[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc9");
      for(Int_t i=0;i<1200;i++)  fFailureS[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc10");
      for(Int_t i=0;i<1200;i++)  fFailureP[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc11");
      for(Int_t i=0;i<1200;i++)  fRecons[i]=(Int_t)mc1->GetBinContent(i+1) ;
      mc1 =(TH1I*)mcfile->Get("mc12");
      for(Int_t i=0;i<1200;i++)  fNonRecons[i]=(Int_t)mc1->GetBinContent(i+1) ;
    }
   mcfile->Close();
 }
 return;
}
//____________________________________________________________________
Bool_t AliITSTrackleterSPDEff::SaveHists() {
  // This (private) method save the histograms on the output file
  // (only if fHistOn is TRUE).
  // Also the histograms from the base class are saved through the 
  // AliITSMultReconstructor::SaveHists() call

  if (!GetHistOn()) return kFALSE;

//  AliITSMultReconstructor::SaveHists(); // this save the histograms of the base class
  fhClustersDPhiAll->Write();
  fhClustersDThetaAll->Write();
  fhClustersDZetaAll->Write();
  fhDPhiVsDThetaAll->Write();
  fhDPhiVsDZetaAll->Write();

  fhClustersDPhiAcc->Write();
  fhClustersDThetaAcc->Write();
  fhClustersDZetaAcc->Write();
  fhDPhiVsDThetaAcc->Write();
  fhDPhiVsDZetaAcc->Write();

  fhetaTracklets->Write();
  fhphiTracklets->Write();
  fhetaClustersLay1->Write();
  fhphiClustersLay1->Write();

  fhClustersDPhiInterpAll->Write();
  fhClustersDThetaInterpAll->Write();
  fhClustersDZetaInterpAll->Write();
  fhDPhiVsDThetaInterpAll->Write();
  fhDPhiVsDZetaInterpAll->Write();

  fhClustersDPhiInterpAcc->Write();
  fhClustersDThetaInterpAcc->Write();
  fhClustersDZetaInterpAcc->Write();
  fhDPhiVsDThetaInterpAcc->Write();
  fhDPhiVsDZetaInterpAcc->Write();

  fhetaClustersLay2->Write();
  fhphiClustersLay2->Write();
  fhClustersInChip->Write();
  for (Int_t nhist=0;nhist<80;nhist++){
    fhClustersInModuleLay1[nhist]->Write(); 
  }
  for (Int_t nhist=0;nhist<160;nhist++){
    fhClustersInModuleLay2[nhist]->Write(); 
  }
  return kTRUE;
}
//__________________________________________________________
Bool_t AliITSTrackleterSPDEff::WriteHistosToFile(TString filename, Option_t* option) {
  //
  // Saves the histograms into a tree and saves the trees into a file
  // Also the histograms from the base class are saved 
  //
  if (!GetHistOn()) return kFALSE;
  if (!strcmp(filename.Data(),"")) {
     AliWarning("WriteHistosToFile: null output filename!");
     return kFALSE;
  }
  TFile *hFile=new TFile(filename.Data(),option,
                         "The File containing the histos for SPD efficiency studies with tracklets");
  if(!SaveHists()) return kFALSE; 
  hFile->Write();
  hFile->Close();
  return kTRUE;
}
//____________________________________________________________
void AliITSTrackleterSPDEff::BookHistos() {
//
// This method books addtitional histograms 
// w.r.t. those of the base class.
// In particular, the differences of cluster coordinate between the two SPD
// layers are computed in the interpolation phase
//
  if (! GetHistOn()) { AliInfo("Call SetHistOn(kTRUE) first"); return;}
//
  fhClustersDPhiAcc   = new TH1F("dphiacc",  "dphi",  100,0.,0.1);
  fhClustersDPhiAcc->SetDirectory(0);
  fhClustersDThetaAcc = new TH1F("dthetaacc","dtheta",100,-0.1,0.1);
  fhClustersDThetaAcc->SetDirectory(0);
  fhClustersDZetaAcc = new TH1F("dzetaacc","dzeta",100,-1.,1.);
  fhClustersDZetaAcc->SetDirectory(0);

  fhDPhiVsDZetaAcc = new TH2F("dphiVsDzetaacc","",100,-1.,1.,100,0.,0.1);
  fhDPhiVsDZetaAcc->SetDirectory(0);
  fhDPhiVsDThetaAcc = new TH2F("dphiVsDthetaAcc","",100,-0.1,0.1,100,0.,0.1);
  fhDPhiVsDThetaAcc->SetDirectory(0);

  fhClustersDPhiAll   = new TH1F("dphiall",  "dphi",  100,0.0,0.5);
  fhClustersDPhiAll->SetDirectory(0);
  fhClustersDThetaAll = new TH1F("dthetaall","dtheta",100,-0.5,0.5);
  fhClustersDThetaAll->SetDirectory(0);
  fhClustersDZetaAll = new TH1F("dzetaall","dzeta",100,-5.,5.);
  fhClustersDZetaAll->SetDirectory(0);

  fhDPhiVsDZetaAll = new TH2F("dphiVsDzetaall","",100,-5.,5.,100,0.,0.5);
  fhDPhiVsDZetaAll->SetDirectory(0);
  fhDPhiVsDThetaAll = new TH2F("dphiVsDthetaAll","",100,-0.5,0.5,100,0.,0.5);
  fhDPhiVsDThetaAll->SetDirectory(0);

  fhetaTracklets  = new TH1F("etaTracklets",  "eta",  100,-2.,2.);
  fhetaTracklets->SetDirectory(0);
  fhphiTracklets  = new TH1F("phiTracklets",  "phi",  100, 0., 2*TMath::Pi());
  fhphiTracklets->SetDirectory(0);
  fhetaClustersLay1  = new TH1F("etaClustersLay1",  "etaCl1",  100,-2.,2.);
  fhetaClustersLay1->SetDirectory(0);
  fhphiClustersLay1  = new TH1F("phiClustersLay1", "phiCl1", 100, 0., 2*TMath::Pi());
  fhphiClustersLay1->SetDirectory(0);
//
  fhClustersDPhiInterpAcc   = new TH1F("dphiaccInterp",  "dphi Interpolation phase",  100,0.,0.1);
  fhClustersDPhiInterpAcc->SetDirectory(0);
  fhClustersDThetaInterpAcc = new TH1F("dthetaaccInterp","dtheta Interpolation phase",100,-0.1,0.1);
  fhClustersDThetaInterpAcc->SetDirectory(0);
  fhClustersDZetaInterpAcc = new TH1F("dzetaaccInterp","dzeta Interpolation phase",100,-1.,1.);
  fhClustersDZetaInterpAcc->SetDirectory(0);

  fhDPhiVsDZetaInterpAcc = new TH2F("dphiVsDzetaaccInterp","dphiVsDzeta Interpolation phase",100,-1.,1.,100,0.,0.1);
  fhDPhiVsDZetaInterpAcc->SetDirectory(0);
  fhDPhiVsDThetaInterpAcc = new TH2F("dphiVsDthetaAccInterp","dphiVsDtheta Interpolation phase",100,-0.1,0.1,100,0.,0.1);
  fhDPhiVsDThetaInterpAcc->SetDirectory(0);

  fhClustersDPhiInterpAll   = new TH1F("dphiallInterp",  "dphi Interpolation phase",  100,0.0,0.5);
  fhClustersDPhiInterpAll->SetDirectory(0);
  fhClustersDThetaInterpAll = new TH1F("dthetaallInterp","dtheta Interpolation phase",100,-0.5,0.5);
  fhClustersDThetaInterpAll->SetDirectory(0);
  fhClustersDZetaInterpAll = new TH1F("dzetaallInterp","dzeta Interpolation phase",100,-5.,5.);
  fhClustersDZetaInterpAll->SetDirectory(0);

  fhDPhiVsDZetaInterpAll = new TH2F("dphiVsDzetaallInterp","dphiVsDzeta Interpolation phase",100,-5.,5.,100,0.,0.5);
  fhDPhiVsDZetaInterpAll->SetDirectory(0);
  fhDPhiVsDThetaInterpAll = new TH2F("dphiVsDthetaAllInterp","dphiVsDtheta Interpolation phase",100,-0.5,0.5,100,0.,0.5);
  fhDPhiVsDThetaInterpAll->SetDirectory(0);

  fhetaClustersLay2  = new TH1F("etaClustersLay2",  "etaCl2",  100,-2.,2.);
  fhetaClustersLay2->SetDirectory(0);
  fhphiClustersLay2  = new TH1F("phiClustersLay2", "phiCl2", 100, 0., 2*TMath::Pi());
  fhphiClustersLay2->SetDirectory(0);
  fhClustersInChip = new TH1F("fhClustersInChip", "ClustersPerChip", 1200, -0.5, 1199.5);
  fhClustersInChip->SetDirectory(0);
// each chip is divided 8(z) x 4(y), i.e. in 32 squares, each containing 4 columns and 64 rows.
  Float_t bz[160]; const Float_t kconv = 1.0E-04; // converts microns to cm.
  for(Int_t i=0;i<160;i++) bz[i] = 425.0; // most are 425 microns except below
  bz[ 31] = bz[ 32] = 625.0; // first chip boundry
  bz[ 63] = bz[ 64] = 625.0; // first chip boundry
  bz[ 95] = bz[ 96] = 625.0; // first chip boundry
  bz[127] = bz[128] = 625.0; // first chip boundry
  Double_t xbins[41]; // each bin in x (Z loc coordinate) includes 4 columns
  //xbins[0]=0;
  Float_t xmn,xmx,zmn,zmx;
  if(!fPlaneEffSPD->GetBlockBoundaries(0,xmn,xmx,zmn,zmx)) AliWarning("Could not book histo properly");
  xbins[0]=(Double_t)zmn;
  for(Int_t i=0;i<40;i++) {
   xbins[i+1]=xbins[i] + (bz[4*i]+bz[4*i+1]+bz[4*i+2]+bz[4*i+3])*kconv; 
  }
  TString histname="ClustersLay1_mod_",aux;
  fhClustersInModuleLay1 =new TH2F*[80];
  for (Int_t nhist=0;nhist<80;nhist++){
    aux=histname;
    aux+=nhist;
    //  
    fhClustersInModuleLay1[nhist]=new TH2F("histname","histname",40,xbins,4,(Double_t)xmn,(Double_t)xmx); 
    fhClustersInModuleLay1[nhist]->SetName(aux.Data());
    fhClustersInModuleLay1[nhist]->SetTitle(aux.Data());
    fhClustersInModuleLay1[nhist]->SetDirectory(0);
  }
  histname="ClustersLay2_mod_";
  fhClustersInModuleLay2 =new TH2F*[160];
  for (Int_t nhist=0;nhist<160;nhist++){
    aux=histname;
    aux+=nhist;
    fhClustersInModuleLay2[nhist]=new TH2F("histname","histname",40,xbins,4,(Double_t)xmn,(Double_t)xmx);
    fhClustersInModuleLay2[nhist]->SetName(aux.Data());
    fhClustersInModuleLay2[nhist]->SetTitle(aux.Data());
    fhClustersInModuleLay2[nhist]->SetDirectory(0);
  }
//
  return;
}
//____________________________________________________________
void AliITSTrackleterSPDEff::DeleteHistos() {
//
// Private method to delete Histograms from memory 
// it is called. e.g., by the destructor.
//
// form AliITSMultReconstructor
  if(fhClustersDPhiAcc) {delete fhClustersDPhiAcc; fhClustersDPhiAcc=0;}
  if(fhClustersDThetaAcc) {delete fhClustersDThetaAcc; fhClustersDThetaAcc=0;}
  if(fhClustersDZetaAcc) {delete fhClustersDZetaAcc; fhClustersDZetaAcc=0;}
  if(fhClustersDPhiAll) {delete fhClustersDPhiAll; fhClustersDPhiAll=0;}
  if(fhClustersDThetaAll) {delete fhClustersDThetaAll; fhClustersDThetaAll=0;}
  if(fhClustersDZetaAll) {delete fhClustersDZetaAll; fhClustersDZetaAll=0;}
  if(fhDPhiVsDThetaAll) {delete fhDPhiVsDThetaAll; fhDPhiVsDThetaAll=0;}
  if(fhDPhiVsDThetaAcc) {delete fhDPhiVsDThetaAcc; fhDPhiVsDThetaAcc=0;}
  if(fhDPhiVsDZetaAll) {delete fhDPhiVsDZetaAll; fhDPhiVsDZetaAll=0;}
  if(fhDPhiVsDZetaAcc) {delete fhDPhiVsDZetaAcc; fhDPhiVsDZetaAcc=0;}
  if(fhetaTracklets) {delete fhetaTracklets; fhetaTracklets=0;}
  if(fhphiTracklets) {delete fhphiTracklets; fhphiTracklets=0;}
  if(fhetaClustersLay1) {delete fhetaClustersLay1; fhetaClustersLay1=0;}
  if(fhphiClustersLay1) {delete fhphiClustersLay1; fhphiClustersLay1=0;}
//
    if(fhClustersDPhiInterpAcc) {delete fhClustersDPhiInterpAcc; fhClustersDPhiInterpAcc=0;}
    if(fhClustersDThetaInterpAcc) {delete fhClustersDThetaInterpAcc; fhClustersDThetaInterpAcc=0;}
    if(fhClustersDZetaInterpAcc) {delete fhClustersDZetaInterpAcc; fhClustersDZetaInterpAcc=0;}
    if(fhClustersDPhiInterpAll) {delete fhClustersDPhiInterpAll; fhClustersDPhiInterpAll=0;}
    if(fhClustersDThetaInterpAll) {delete fhClustersDThetaInterpAll; fhClustersDThetaInterpAll=0;}
    if(fhClustersDZetaInterpAll) {delete fhClustersDZetaInterpAll; fhClustersDZetaInterpAll=0;}
    if(fhDPhiVsDThetaInterpAll) {delete fhDPhiVsDThetaInterpAll; fhDPhiVsDThetaInterpAll=0;}
    if(fhDPhiVsDThetaInterpAcc) {delete fhDPhiVsDThetaInterpAcc; fhDPhiVsDThetaInterpAcc=0;}
    if(fhDPhiVsDZetaInterpAll) {delete fhDPhiVsDZetaInterpAll; fhDPhiVsDZetaInterpAll=0;}
    if(fhDPhiVsDZetaInterpAcc) {delete fhDPhiVsDZetaInterpAcc; fhDPhiVsDZetaInterpAcc=0;}
    if(fhetaClustersLay2) {delete fhetaClustersLay2; fhetaClustersLay2=0;}
    if(fhphiClustersLay2) {delete fhphiClustersLay2; fhphiClustersLay2=0;}
    if(fhClustersInChip) {delete fhClustersInChip; fhClustersInChip=0;}
    if(fhClustersInModuleLay1) {
      for (Int_t i=0; i<80; i++ ) delete fhClustersInModuleLay1[i];
      delete [] fhClustersInModuleLay1; fhClustersInModuleLay1=0;
    }
    if(fhClustersInModuleLay2) {
      for (Int_t i=0; i<160; i++ ) delete fhClustersInModuleLay2[i];
      delete [] fhClustersInModuleLay2; fhClustersInModuleLay2=0;
    }
}
//_______________________________________________________________
Bool_t AliITSTrackleterSPDEff::IsReconstructableAt(Int_t layer,Int_t iC,Int_t ipart,
                                                   const Float_t* vtx, const AliStack *stack, TTree *ref) {
// This (private) method can be used only for MC events, where both AliStack and the TrackReference
// are available. 
// It is used to asses whether a tracklet prediction is reconstructable or not at the other layer
// Input: 
//      - Int_t layer (either 0 or 1): layer which you want to check if the tracklete can be 
//                                     reconstructed at
//      - Int_t iC : cluster index used to build the tracklet prediction 
//                   if layer=0 ==> iC=iC2 ; elseif layer=1 ==> iC=iC1
//      - Float_t* vtx: actual event vertex
//      - stack: pointer to Stack
//      - ref:   pointer to TTRee of TrackReference
Bool_t ret=kFALSE; // returned value
Float_t trefLayExtr[3]; // equivalent to fClustersLay1/fClustersLay2 but for the track reference
if(!fMC) {AliError("This method works only if SetMC() has been called"); return ret;}
if(!stack) {AliError("null pointer to MC stack"); return ret;}
if(!ref)  {AliError("null pointer to TrackReference Tree"); return ret;}
if(ipart >= stack->GetNtrack()) {AliError("this track label is not in MC stack"); return ret;}
if(layer<0 || layer>1) {AliError("You can extrapolate either at lay 0 or at lay 1"); return ret;}

AliTrackReference *tref=0x0;
Int_t imatch=-100; // index of the track in TrackReference which matches with ipart
Int_t nentries = (Int_t)ref->GetEntries();
TClonesArray *tcaRef = new TClonesArray("AliTrackReference");
TBranch *br = ref->GetBranch("TrackReferences");
br->SetAddress(&tcaRef);
for(Int_t itrack=0;itrack<nentries;itrack++) { // loop over all Tracks in TrackReference to match the ipart one
  br->GetEntry(itrack);
  Int_t nref=tcaRef->GetEntriesFast();
  if(nref>0) { //it is enough to look at the first one
    tref=(AliTrackReference*)tcaRef->At(0); // it is enough to look at the first one
    if(tref->GetTrack()==ipart) {imatch=itrack; break;}
  }
}
if(imatch<0) {AliWarning(Form("Could not find AliTrackReference for particle %d",ipart)); return kFALSE;}
br->GetEntry(imatch); // redundant, nevertheless ...
Int_t nref=tcaRef->GetEntriesFast();
for(Int_t iref=0;iref<nref;iref++) { // loop over all the refs of the matching track
  tref=(AliTrackReference*)tcaRef->At(iref);
  if(tref->R()>10) continue; // not SPD ref
  if(layer==0 && tref->R()>5) continue; // ref on SPD outer layer
  if(layer==1 && tref->R()<5) continue; // ref on SPD inner layer

// compute the proper quantities for this tref, as was done for fClustersLay1/2
  Float_t x = tref->X() - vtx[0];
  Float_t y = tref->Y() - vtx[1];
  Float_t z = tref->Z() - vtx[2];

  Float_t r    = TMath::Sqrt(x*x + y*y +z*z);

  trefLayExtr[0] = TMath::ACos(z/r);                   // Store Theta
  trefLayExtr[1] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi
  trefLayExtr[2] = z;                                    // Store z

  if(layer==1) { // try to see if it is reconstructable at the outer layer
// find the difference in angles
    Float_t dPhi   = TMath::Abs(trefLayExtr[1] - fClustersLay1[iC][1]);
    // take into account boundary condition
    if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;

    // find the difference in z (between linear projection from layer 1
    // and the actual point: Dzeta= z1/r1*r2 -z2)
    Float_t r2    = trefLayExtr[2]/TMath::Cos(trefLayExtr[0]);
    Float_t dZeta = TMath::Cos(fClustersLay1[iC][0])*r2 - trefLayExtr[2];

    // make "elliptical" cut in Phi and Zeta!
    Float_t d = TMath::Sqrt(dPhi*dPhi/fPhiWindowL2/fPhiWindowL2 +
                              dZeta*dZeta/fZetaWindowL2/fZetaWindowL2);
    if (d<1) {ret=kTRUE; break;}
  }
  if(layer==0) { // try to see if it is reconstructable at the inner layer

    // find the difference in angles
    Float_t dPhi   = TMath::Abs(fClustersLay2[iC][1] - trefLayExtr[1]);
    // take into account boundary condition
    if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;

    // find the difference in z (between linear projection from layer 2
    // and the actual point: Dzeta= z2/r2*r1 -z1)
    Float_t r1    = trefLayExtr[2]/TMath::Cos(trefLayExtr[0]);
    Float_t dZeta = TMath::Cos(fClustersLay2[iC][0])*r1 - trefLayExtr[2];

    // make "elliptical" cut in Phi and Zeta!
    Float_t d = TMath::Sqrt(dPhi*dPhi/fPhiWindowL1/fPhiWindowL1 +
                            dZeta*dZeta/fZetaWindowL1/fZetaWindowL1);
    if (d<1) {ret=kTRUE; break;};
  }
}
delete tcaRef;
return ret;
}
//_________________________________________________________________________
void AliITSTrackleterSPDEff::ReflectClusterAroundZAxisForLayer(Int_t ilayer){
//
// this method apply a rotation by 180 degree around the Z (beam) axis to all 
// the RecPoints in a given layer to be used to build tracklets.
// **************** VERY IMPORTANT:: ***************
// It must be called just after LoadClusterArrays, since afterwards the datamember
// fClustersLay1[iC1][0] and fClustersLay1[iC1][1] are redefined using polar coordinate 
// instead of Cartesian
//
if(ilayer<0 || ilayer>1) {AliInfo("Input argument (ilayer) should be either 0 or 1: nothing done"); return ;}
AliDebug(3,Form("Applying a rotation by 180 degree around z axiz to all clusters on layer %d",ilayer));
if(ilayer==0) {
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {
    fClustersLay1[iC1][0]*=-1;
    fClustersLay1[iC1][1]*=-1;
  }
}
if(ilayer==1) {
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {
    fClustersLay2[iC2][0]*=-1;
    fClustersLay2[iC2][1]*=-1;
  }
}
return;
}
//____________________________________________________________________________
Int_t AliITSTrackleterSPDEff::Clusters2Tracks(AliESDEvent *esd){
// This method is used to find the tracklets. 
// It is called from AliReconstruction
// The vertex is supposed to be associated to the Tracker (i.e. to this) already
// The cluster is supposed to be associated to the Tracker already
// In case Monte Carlo is required, the appropriate linking to Stack and TrackRef is attempted 
//
  Int_t rc=1;
  // apply cuts on the vertex quality
  const AliESDVertex *vertex = esd->GetVertex();
  if(vertex->GetNContributors()<fMinContVtx) return 0;
  //
  AliRunLoader* runLoader = AliRunLoader::Instance();
  if (!runLoader) {
    Error("Clusters2Tracks", "no run loader found");
    return rc;
  }
  AliStack *pStack=0x0; TTree *tRefTree=0x0;
  if(GetMC()) {
    runLoader->LoadKinematics("read");
    runLoader->LoadTrackRefs("read");
    pStack= runLoader->Stack();
    tRefTree= runLoader->TreeTR();
  }
  Reconstruct(pStack,tRefTree);

  if (GetLightBkgStudyInParallel()) {
    AliStack *dummy1=0x0; TTree *dummy2=0x0;
    ReflectClusterAroundZAxisForLayer(1);
    Reconstruct(dummy1,dummy2,kTRUE);
  }
  return 0;
}
//____________________________________________________________________________
Int_t AliITSTrackleterSPDEff::PostProcess(AliESDEvent *){
// 
// It is called from AliReconstruction
// 
// 
// 
//
  Int_t rc=0;
  if(GetMC()) SavePredictionMC("TrackletsMCpred.root");
  if(GetHistOn()) rc=(Int_t)WriteHistosToFile();
  if(GetLightBkgStudyInParallel()) {
    TString name="AliITSPlaneEffSPDtrackletBkg.root";
    TFile* pefile = TFile::Open(name, "RECREATE");
    rc*=fPlaneEffBkg->Write();
    pefile->Close();
  }
  return rc;
}
//____________________________________________________________________
void
AliITSTrackleterSPDEff::LoadClusterArrays(TTree* itsClusterTree) {
  // This method
  // - gets the clusters from the cluster tree
  // - convert them into global coordinates
  // - store them in the internal arrays
  // - count the number of cluster-fired chips

  //AliDebug(1,"Loading clusters and cluster-fired chips ...");

  fNClustersLay1 = 0;
  fNClustersLay2 = 0;

  TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
  TBranch* itsClusterBranch=itsClusterTree->GetBranch("ITSRecPoints");

  itsClusterBranch->SetAddress(&itsClusters);

  Int_t nItsSubs = (Int_t)itsClusterTree->GetEntries();
  Float_t cluGlo[3]={0.,0.,0.};

  // loop over the its subdetectors
  for (Int_t iIts=0; iIts < nItsSubs; iIts++) {

    if (!itsClusterTree->GetEvent(iIts))
      continue;

    Int_t nClusters = itsClusters->GetEntriesFast();

    // number of clusters in each chip of the current module
    Int_t layer = 0;

    // loop over clusters
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);

      layer = cluster->GetLayer();
      if (layer>1) continue;

      cluster->GetGlobalXYZ(cluGlo);
      Float_t x = cluGlo[0];
      Float_t y = cluGlo[1];
      Float_t z = cluGlo[2];

      if (layer==0) {
        fClustersLay1[fNClustersLay1][0] = x;
        fClustersLay1[fNClustersLay1][1] = y;
        fClustersLay1[fNClustersLay1][2] = z;

        for (Int_t i=0; i<3; i++)
                fClustersLay1[fNClustersLay1][3+i] = cluster->GetLabel(i);
        fNClustersLay1++;
        if(fHistOn) { 
          Int_t det=cluster->GetDetectorIndex();
          if(det<0 || det>79) {AliError("Cluster with det. index out of boundaries"); return;}
          fhClustersInModuleLay1[det]->Fill((Double_t)cluster->GetDetLocalZ(),(Double_t)cluster->GetDetLocalX());
        }
      }
      if (layer==1) {
        fClustersLay2[fNClustersLay2][0] = x;
        fClustersLay2[fNClustersLay2][1] = y;
        fClustersLay2[fNClustersLay2][2] = z;

        for (Int_t i=0; i<3; i++)
                fClustersLay2[fNClustersLay2][3+i] = cluster->GetLabel(i);
        fNClustersLay2++;
        if(fHistOn) {
          Int_t det=cluster->GetDetectorIndex();
          if(det<0 || det>159) {AliError("Cluster with det. index out of boundaries"); return;}
          fhClustersInModuleLay2[det]->Fill((Double_t)cluster->GetDetLocalZ(),(Double_t)cluster->GetDetLocalX());
        }
      }

    }// end of cluster loop

  } // end of its "subdetector" loop
  if (itsClusters) {
    itsClusters->Delete();
    delete itsClusters;
    itsClusters = 0;
  }
  AliDebug(1,Form("(clusters in layer 1 : %d,  layer 2: %d)",fNClustersLay1,fNClustersLay2));
}
//_________________________________________________________________________
void
AliITSTrackleterSPDEff::SetLightBkgStudyInParallel(Bool_t b) {
//     This method:
//  - set Bool_t fLightBackgroundStudyInParallel = b 
//    a) if you set this kTRUE, then the estimation of the 
//      SPD efficiency is done as usual for data, but in 
//      parallel a light (i.e. without control histograms, etc.) 
//      evaluation of combinatorial background is performed
//      with the usual ReflectClusterAroundZAxisForLayer method.
//    b) if you set this kFALSE, then you would not have a second 
//      container for PlaneEfficiency statistics to be used for background 
//      (fPlaneEffBkg=0). If you want to have a full evaluation of the 
//      background (with all control histograms and additional data 
//      members referring to the background) then you have to call the 
//      method SetReflectClusterAroundZAxisForLayer(kTRUE) esplicitily
  fLightBkgStudyInParallel=b; 
  if(fLightBkgStudyInParallel) {
    if(!fPlaneEffBkg) fPlaneEffBkg = new AliITSPlaneEffSPD();   
  }
  else {
    delete fPlaneEffBkg;
    fPlaneEffBkg=0;
  }
}
//______________________________________________________________
void AliITSTrackleterSPDEff::SetReflectClusterAroundZAxisForLayer(Int_t ilayer,Bool_t b){  
//
// method to study residual background:
// Input b= KTRUE --> reflect the clusters 
//      ilayer (either 0 or 1) --> which SPD layers should be reflected
//
    if(b) {AliInfo(Form("All clusters on layer %d will be rotated by 180 deg around z",ilayer));
           SetLightBkgStudyInParallel(kFALSE);}
    if(ilayer==0) fReflectClusterAroundZAxisForLayer0=b;                   // a rotation by 180degree around the Z axis  
    else if(ilayer==1) fReflectClusterAroundZAxisForLayer1=b;              // (x->-x; y->-y) to all RecPoints on a 
    else AliInfo("Nothing done: input argument (ilayer) either 0 or 1");   // given layer is applied. In such a way 
  }                    
