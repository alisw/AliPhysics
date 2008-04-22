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

/* $Id$ */

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

#include <TFile.h>
#include <TParticle.h>
#include <TSystem.h>
#include <Riostream.h>

#include "AliITSMultReconstructor.h"
#include "AliITSTrackleterSPDEff.h"
#include "AliITSgeomTGeo.h"
#include "AliLog.h"
#include "AliITSPlaneEffSPD.h"
#include "AliStack.h"

//____________________________________________________________________
ClassImp(AliITSTrackleterSPDEff)


//____________________________________________________________________
AliITSTrackleterSPDEff::AliITSTrackleterSPDEff():
AliITSMultReconstructor(),
fAssociationFlag1(0),
fChipPredOnLay2(0),
fChipPredOnLay1(0),
fNTracklets1(0),
fPhiWindowL1(0),
fZetaWindowL1(0),
fOnlyOneTrackletPerC1(0),
fPlaneEffSPD(0),
fMC(0),
fUseOnlyPrimaryForPred(0),
fUseOnlySecondaryForPred(0), 
fUseOnlySameParticle(0),
fUseOnlyDifferentParticle(0),
fUseOnlyStableParticle(0),
fPredictionPrimary(0),
fPredictionSecondary(0),
fClusterPrimary(0),
fClusterSecondary(0),
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
fhphiClustersLay2(0)
{

  // Method to check the SPD chips efficiencies by using tracklets 


  SetPhiWindowL1();
  SetZetaWindowL1();
  SetOnlyOneTrackletPerC1();

  fAssociationFlag1   = new Bool_t[300000];
  fChipPredOnLay2     = new UInt_t[300000];
  fChipPredOnLay1     = new UInt_t[300000];

  for(Int_t i=0; i<300000; i++) {
    fAssociationFlag1[i]   = kFALSE;
  }

  if (GetHistOn()) BookHistos();

  fPlaneEffSPD = new AliITSPlaneEffSPD();
}
//______________________________________________________________________
AliITSTrackleterSPDEff::AliITSTrackleterSPDEff(const AliITSTrackleterSPDEff &mr) : AliITSMultReconstructor(mr),
fAssociationFlag1(mr.fAssociationFlag1),
fChipPredOnLay2(mr.fChipPredOnLay2),
fChipPredOnLay1(mr.fChipPredOnLay1),
fNTracklets1(mr.fNTracklets1),
fPhiWindowL1(mr.fPhiWindowL1),
fZetaWindowL1(mr.fZetaWindowL1),
fOnlyOneTrackletPerC1(mr.fOnlyOneTrackletPerC1),
fPlaneEffSPD(mr.fPlaneEffSPD),
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
fhphiClustersLay2(mr.fhphiClustersLay2)
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

  // delete histograms
  DeleteHistos();

  delete [] fAssociationFlag1;

  delete [] fChipPredOnLay2;
  delete [] fChipPredOnLay1;

  delete [] fPredictionPrimary;  
  delete [] fPredictionSecondary; 
  delete [] fClusterPrimary;  
  delete [] fClusterSecondary; 

  // delete PlaneEff
  delete fPlaneEffSPD;
}
//____________________________________________________________________
void
AliITSTrackleterSPDEff::Reconstruct(TTree* clusterTree, Float_t* vtx, Float_t* /* vtxRes*/,AliStack *pStack) {
  //
  // - calls LoadClusterArray that finds the position of the clusters
  //   (in global coord) 
  // - convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). Find the extrapolation/interpolation point.
  // - Find the chip corresponding to that
  // - Check if there is a cluster near that point  
  //

  // reset counters
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0; 
  // loading the clusters 
  LoadClusterArrays(clusterTree);
  if(fMC && !pStack) {AliError("You asked for MC infos but AliStack not properly loaded"); return;}
  Bool_t found;
  Int_t nfTraPred1=0;  Int_t ntTraPred1=0;
  Int_t nfTraPred2=0;  Int_t ntTraPred2=0;
  Int_t nfClu1=0;      Int_t ntClu1=0; 
  Int_t nfClu2=0;      Int_t ntClu2=0;
  

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
 
    if (fHistOn) {
      Float_t eta=fClustersLay1[iC1][0];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaClustersLay1->Fill(eta);    
      fhphiClustersLay1->Fill(fClustersLay1[iC1][1]);
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
      AliWarning(Form("Reconstruct: cannot find chip prediction on inner layer for cluster %d on the outer layer",iC2)); 
      key=999999;
    }
    nfTraPred1+=(Int_t)found; // this for debugging purpose
    ntTraPred1++;             // to check efficiency of the method FindChip
    fChipPredOnLay1[iC2] = key;
    fAssociationFlag[iC2] = kFALSE;
 
    if (fHistOn) {
      Float_t eta=fClustersLay2[iC2][0];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaClustersLay2->Fill(eta);
      fhphiClustersLay2->Fill(fClustersLay2[iC2][1]);
    }
  }  
  
  //###########################################################

 // First part : Extrapolation to Layer 2 

  // Loop on layer 1 
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    

    // reset of variables for multiple candidates
    Int_t  iC2WithBestDist = 0;     // reset 
    Float_t distmin        = 100.;  // just to put a huge number! 
    Float_t dPhimin        = 0.;  // Used for histograms only! 
    Float_t dThetamin      = 0.;  // Used for histograms only! 
    Float_t dZetamin       = 0.;  // Used for histograms only! 

    // in any case, if MC has been required, store statistics of primaries and secondaries
    if (fMC) {
       Int_t lab1=(Int_t)fClustersLay1[iC1][3];
       Int_t lab2=(Int_t)fClustersLay1[iC1][4];
       Int_t lab3=(Int_t)fClustersLay1[iC1][5];
       // do it always as a function of the chip number used to built the prediction
       found=FindChip(key,0,vtx,fClustersLay1[iC1][0],fClustersLay1[iC1][1],fClustersLay1[iC1][2]);
       if (!found) {AliWarning(
         Form("Reconstruct MC: cannot find chip on inner layer for cluster %d",iC1)); }
       else {
         if((lab1 != -2  &&  PrimaryTrackChecker(lab1,pStack) ) ||
            (lab2 != -2  &&  PrimaryTrackChecker(lab2,pStack) ) ||
            (lab3 != -2  &&  PrimaryTrackChecker(lab3,pStack))) 
         { // this cluster is from a primary particle
           fClusterPrimary[key]++;
           if(fUseOnlySecondaryForPred) continue; // skip this tracklet built with a primary track
         } else { // this cluster is from a secondary particle
            fClusterSecondary[key]++;
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

 	if (fHistOn) {
	  fhClustersDPhiAll->Fill(dPhi);    
	  fhClustersDThetaAll->Fill(dTheta);    
	  fhClustersDZetaAll->Fill(dZeta);    
	  fhDPhiVsDThetaAll->Fill(dTheta, dPhi);
	  fhDPhiVsDZetaAll->Fill(dZeta, dPhi);
	}

	// make "elliptical" cut in Phi and Zeta! 
	Float_t d = TMath::Sqrt(dPhi*dPhi/fPhiWindow/fPhiWindow + 
                                dZeta*dZeta/fZetaWindow/fZetaWindow);

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

      if (fHistOn) {
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

      if (fHistOn) {
	Float_t eta=fTracklets[fNTracklets][0];
	eta= TMath::Tan(eta/2.);
	eta=-TMath::Log(eta);
	fhetaTracklets->Fill(eta);    
	fhphiTracklets->Fill(fTracklets[fNTracklets][1]);    
      }

// Check that this cluster is still in the same chip (here you pass also Zvtx for better computation)
      found=FindChip(key,1,vtx,fClustersLay2[iC2WithBestDist][0],fClustersLay2[iC2WithBestDist][1],fClustersLay2[iC2WithBestDist][2]);
      if(!found){
        AliWarning(
         Form("Reconstruct: cannot find chip on outer layer for cluster %d",iC2WithBestDist));
        key=999999;
      }
      nfClu2+=(Int_t)found; // this for debugging purpose
      ntClu2++;             // to check efficiency of the method FindChip
      if(key<1200) { // the Chip has been found
        if(fMC) { // this part only for MC
          // Int_t labc1=(Int_t)fClustersLay2[iC2WithBestDist][3];
          // Int_t labc2=(Int_t)fClustersLay2[iC2WithBestDist][4];
          // Int_t labc3=(Int_t)fClustersLay2[iC2WithBestDist][5];
          if (fUseOnlyDifferentParticle && label2 < 3) continue; // same label (reject it)
          if (fUseOnlySameParticle && label2 == 3) continue;      // different label (reject it)
        }

        if (key==fChipPredOnLay2[iC1]) { // this control seems too loose: has to be checked !
          // OK, success
                fPlaneEffSPD->UpDatePlaneEff(kTRUE,key); // success
        }
        else {
                fPlaneEffSPD->UpDatePlaneEff(kTRUE,key); // this should not be a failure
                                                         // (might be in the tracking tollerance)
        }
      }

      fNTracklets++;

    } // if any cluster found --> increment statistics by 1 failure (provided you have chip prediction)
    else if (fChipPredOnLay2[iC1]<1200) fPlaneEffSPD->UpDatePlaneEff(kFALSE,fChipPredOnLay2[iC1]);

  } // end of loop over clusters in layer 1

    fNTracklets1=fNTracklets;

//###################################################################

  // Second part : Interpolation to Layer 1 

  // Loop on layer 2 
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {    

    // reset of variables for multiple candidates
    Int_t  iC1WithBestDist = 0;     // reset 
    Float_t distmin        = 100.;  // just to put a huge number! 
    Float_t dPhimin        = 0.;  // Used for histograms only! 
    Float_t dThetamin      = 0.;  // Used for histograms only! 
    Float_t dZetamin       = 0.;  // Used for histograms only! 

    // in any case, if MC has been required, store statistics of primaries and secondaries
    if (fMC) {
       Int_t lab1=(Int_t)fClustersLay2[iC2][3];
       Int_t lab2=(Int_t)fClustersLay2[iC2][4];
       Int_t lab3=(Int_t)fClustersLay2[iC2][5];
       // do it always as a function of the chip number used to built the prediction
       found=FindChip(key,1,vtx,fClustersLay2[iC2][0],fClustersLay2[iC2][1],fClustersLay2[iC2][2]);
       if (!found) {AliWarning(
         Form("Reconstruct MC: cannot find chip on outer layer for cluster %d",iC2)); }
       else {
         if((lab1 != -2  &&  PrimaryTrackChecker(lab1,pStack) ) ||
            (lab2 != -2  &&  PrimaryTrackChecker(lab2,pStack) ) ||
            (lab3 != -2  &&  PrimaryTrackChecker(lab3,pStack))) 
         {  // this cluster is from a primary particle
            fClusterPrimary[key]++;
            if(fUseOnlySecondaryForPred) continue; //  skip this tracklet built with a primary track
         } else { // this cluster is from a secondary particle
           fClusterSecondary[key]++;
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


 	if (fHistOn) {
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
    
    if (distmin<100) { // This means that a cluster in layer 1 was found that mathes with iC2

      if (fHistOn) {
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
        AliWarning(
         Form("Reconstruct: cannot find chip on inner layer for cluster %d",iC1WithBestDist));
        key=999999;
      }
      nfClu1+=(Int_t)found; // this for debugging purpose
      ntClu1++;             // to check efficiency of the method FindChip
      if(key<1200) {
        if(fMC) { // this part only for MC
          // Int_t labc1=(Int_t)fClustersLay1[iC1WithBestDist][3];
          // Int_t labc2=(Int_t)fClustersLay1[iC1WithBestDist][4];
          // Int_t labc3=(Int_t)fClustersLay1[iC1WithBestDist][5];
          if (fUseOnlyDifferentParticle && label2 < 3) continue; // same label (reject it)
          if (fUseOnlySameParticle && label2 == 3) continue;      // different label (reject it)
        }

        if (key==fChipPredOnLay1[iC2]) { // this control seems too loose: has to be checked !
          // OK, success
                fPlaneEffSPD->UpDatePlaneEff(kTRUE,key); // success
        } else {
                fPlaneEffSPD->UpDatePlaneEff(kTRUE,key); // this should not be a failure
                                                         // (might be in the tracking tollerance)
        }
      }

    fNTracklets++;

    } // if no cluster found --> increment statistics by 1 failure (provided you have chip prediction)
    else if (fChipPredOnLay1[iC2]<1200) fPlaneEffSPD->UpDatePlaneEff(kFALSE,fChipPredOnLay1[iC2]);

  } // end of loop over clusters in layer 2
  
  AliDebug(1,Form("%d tracklets found", fNTracklets));
  AliDebug(1,Form(("Eff. of method FindChip for Track pred. on lay 1 = %d / %d"),nfTraPred1,ntTraPred1));
  AliDebug(1,Form(("Eff. of method FindChip for Track pred. on lay 2 = %d / %d"),nfTraPred2,ntTraPred2));
  AliDebug(1,Form(("Eff. of method FindChip for Cluster on lay 1 = %d / %d"),nfClu1,ntClu1));
  AliDebug(1,Form(("Eff. of method FindChip for Cluster on lay 2 = %d / %d"),nfClu2,ntClu2));
}
//____________________________________________________________________
Bool_t AliITSTrackleterSPDEff::FindChip(UInt_t &key, Int_t layer,  Float_t* vtx, 
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
  else zAbs=r/TMath::Tan(thetaVtx) + vtx[2]; // this is the only way to do for the tracklet prediction
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
  //This function finds the detector crossed by the track
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
Finally you want P in the ABSOLUTE ALICE system.
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
Bool_t AliITSTrackleterSPDEff::SetAngleRange02Pi(Double_t &angle) {
while(angle >=2*TMath::Pi() || angle<0) {
  if(angle >= 2*TMath::Pi()) angle-=2*TMath::Pi();
  if(angle < 0) angle+=2*TMath::Pi();
}
return kTRUE;
}
//___________________________________________________________
Bool_t AliITSTrackleterSPDEff::PrimaryTrackChecker(Int_t ipart,AliStack* stack) {
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

 if(fUseOnlyStableParticle && DecayingTrackChecker(ipart,stack)<2) return kFALSE; 
 return kTRUE;
}
//_____________________________________________________________________________________________
Int_t AliITSTrackleterSPDEff::DecayingTrackChecker(Int_t ipart,AliStack* stack) {
if(!fMC) {AliError("This method works only if SetMC() has been called"); return 0;}
if(!stack) {AliError("null pointer to MC stack"); return 0;}
if(ipart >= stack->GetNtrack()) {AliError("this track label is not in MC stack"); return 0;}

TParticle* part = stack->Particle(ipart);
//TParticle* part0 = stack->Particle(0); // first primary

  Int_t nret=0;
  TParticle* dau = 0;
  Int_t nDau = 0;
  Int_t firstDau = part->GetFirstDaughter();
  if (firstDau > 0) {
    Int_t lastDau = part->GetLastDaughter();
    nDau = lastDau - firstDau + 1;
    //printf("number of daugthers %d \n",nDau);
    if (nDau > 0) {
      //for(Int_t j=firstDau; j<=lastDau; j++) 
      for(Int_t j=firstDau; j<=firstDau; j++) 
                                              { // only first one 
        dau = stack->Particle(j);
        Double_t distx = dau->Vx()-part->Vx();
        Double_t disty = dau->Vy()-part->Vy();
        Double_t distz = dau->Vz()-part->Vz();
        Double_t distR = TMath::Sqrt(distx*distx+disty*disty+distz*distz);
        if (distR > GetRLayer(0)+0.5)  nret=1;  // decay after first pixel layer
        if (distR > GetRLayer(1)+0.5)  nret=2;  // decay after second pixel layer
      }
    }
  } else nret = 3; // stable particle
return nret; 
}
//_________________________________________________________________
void AliITSTrackleterSPDEff::InitPredictionMC() {
if(!fMC) {AliError("This method works only if SetMC() has been called"); return;}
fPredictionPrimary   = new Int_t[1200];
fPredictionSecondary = new Int_t[1200];
fClusterPrimary      = new Int_t[1200];
fClusterSecondary    = new Int_t[1200];
for(Int_t i=0; i<1200; i++) {
 fPredictionPrimary[i]=0;
 fPredictionSecondary[i]=0; 
 fPredictionSecondary[i]=0;
 fClusterSecondary[i]=0;
}
return;
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetPredictionPrimary(const UInt_t key) const {
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fPredictionPrimary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetPredictionSecondary(const UInt_t key) const {
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fPredictionSecondary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetClusterPrimary(const UInt_t key) const {
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fClusterPrimary[(Int_t)key];
}
//______________________________________________________________________
Int_t AliITSTrackleterSPDEff::GetClusterSecondary(const UInt_t key) const {
if (!fMC) {CallWarningMC(); return 0;}
if (key>=1200) {AliWarning("You asked for a non existing chip"); return -999;}
return fClusterSecondary[(Int_t)key];
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
    *os << fPhiWindowL1 <<" "<< fZetaWindowL1 << " " << fPhiWindow <<" "<< fZetaWindow ;
    *os << " " << fMC;
    if(!fMC) {AliInfo("Writing only cuts, no MC info"); return;}
    *os << " " << fUseOnlyPrimaryForPred << " " << fUseOnlySecondaryForPred
        << " " << fUseOnlySameParticle   << " " << fUseOnlyDifferentParticle
        << " " << fUseOnlyStableParticle ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetPredictionPrimary(i)  ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetPredictionSecondary(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetClusterPrimary(i) ;
    for(Int_t i=0;i<1200;i++) *os <<" "<< GetClusterSecondary(i) ;
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

    *is >> fPhiWindowL1 >> fZetaWindowL1 >> fPhiWindow >> fZetaWindow;
    *is >> fMC;
    if(!fMC) {AliInfo("Reading only cuts, no MC info available");return;}
    *is >> fUseOnlyPrimaryForPred >> fUseOnlySecondaryForPred
        >> fUseOnlySameParticle   >> fUseOnlyDifferentParticle
        >> fUseOnlyStableParticle;
    for(Int_t i=0;i<1200;i++) *is >> fPredictionPrimary[i] ;
    for(Int_t i=0;i<1200;i++) *is >> fPredictionSecondary[i] ;
    for(Int_t i=0;i<1200;i++) *is >> fClusterPrimary[i] ;
    for(Int_t i=0;i<1200;i++) *is >> fClusterSecondary[i] ;
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
 if(!fMC) {CallWarningMC(); return;}
 ofstream out(filename.Data(),ios::out | ios::binary);
 out << *this;
 out.close();
return;
}
//____________________________________________________________________
void AliITSTrackleterSPDEff::ReadPredictionMC(TString filename) {
 if(!fMC) {CallWarningMC(); return;}
 if( gSystem->AccessPathName( filename.Data() ) ) {
      AliError( Form( "file (%s) not found", filename.Data() ) );
      return;
   }

 ifstream in(filename.Data(),ios::in | ios::binary);
 in >> *this;
 in.close();
 return;
}
//____________________________________________________________________
Bool_t AliITSTrackleterSPDEff::SaveHists() {
  // This method save the histograms on the output file
  // (only if fHistOn is TRUE).

  if (!GetHistOn()) return kFALSE;

  AliITSMultReconstructor::SaveHists(); // this save the histograms of the base class

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
  return kTRUE;
}
//__________________________________________________________
Bool_t AliITSTrackleterSPDEff::WriteHistosToFile(TString filename, Option_t* option) {
  //
  // Saves the histograms into a tree and saves the trees into a file
  //
  if (!GetHistOn()) return kFALSE;
  if (filename.Data()=="") {
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
  if (! GetHistOn()) { AliInfo("Call SetHistOn(kTRUE) first"); return;}
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
  return;
}
//____________________________________________________________
void AliITSTrackleterSPDEff::DeleteHistos() {
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
}
//_______________________________________________________________
