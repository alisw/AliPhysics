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

//_________________________________________________________________________
// 
//        Implementation of the ITS-SPD trackleter class
//
// It retrieves clusters in the pixels (theta and phi) and finds tracklets.
// These can be used to extract charged particle multiplicity from the ITS.
//
// A tracklet consists of two ITS clusters, one in the first pixel layer and 
// one in the second. The clusters are associated if the differences in 
// Phi (azimuth) and Theta (polar angle) are within fiducial windows.
// In case of multiple candidates the candidate with minimum
// distance is selected. 
//
// Two methods return the number of tracklets and the number of unassociated 
// clusters (i.e. not used in any tracklet) in the first SPD layer
// (GetNTracklets and GetNSingleClusters)
//
// The cuts on phi and theta depend on the interacting system (p-p or Pb-Pb)
// and can be set via AliITSRecoParam class
// (SetPhiWindow and SetThetaWindow)  
// 
// Origin: Tiziano Virgili 
//
// Current support and development: 
//         Domenico Elia, Maria Nicassio (INFN Bari) 
//         Domenico.Elia@ba.infn.it, Maria.Nicassio@ba.infn.it
//
// Most recent updates:
//     - multiple association forbidden (fOnlyOneTrackletPerC2 = kTRUE)    
//     - phi definition changed to ALICE convention (0,2*TMath::pi()) 
//     - cluster coordinates taken with GetGlobalXYZ()
//     - fGeometry removed
//     - number of fired chips on the two layers
//     - option to cut duplicates in the overlaps
//     - options and fiducial cuts via AliITSRecoParam
//     - move from DeltaZeta to DeltaTheta cut
//     - update to the new algorithm by Mariella and Jan Fiete
//     - store also DeltaTheta in the ESD 
//     - less new and delete calls when creating the needed arrays
//
//     - RS: to decrease the number of new/deletes the clusters data are stored 
//           not in float[6] attached to float**, but in 1-D array.
//     - RS: Clusters are sorted in Z in roder to have the same numbering as in the ITS reco
//     - RS: Clusters used by ESDtrack are flagged, this information is passed to AliMulitiplicity object 
//           when storing the tracklets and single cluster info
//_________________________________________________________________________

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBits.h>
#include <TArrayI.h>

#include "AliITSMultReconstructor.h"
#include "AliITSReconstructor.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeom.h"
#include "AliITSgeomTGeo.h"
#include "AliITSDetTypeRec.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#include "AliLog.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"

//____________________________________________________________________
ClassImp(AliITSMultReconstructor)


//____________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor():
fDetTypeRec(0),fESDEvent(0),fTreeRP(0),fUsedClusLay1(0),fUsedClusLay2(0),
fClustersLay1(0),
fClustersLay2(0),
fDetectorIndexClustersLay1(0),
fDetectorIndexClustersLay2(0),
fOverlapFlagClustersLay1(0),
fOverlapFlagClustersLay2(0),
fTracklets(0),
fSClusters(0),
fNClustersLay1(0),
fNClustersLay2(0),
fNTracklets(0),
fNSingleCluster(0),
fPhiWindow(0),
fThetaWindow(0),
fPhiShift(0),
fRemoveClustersFromOverlaps(0),
fPhiOverlapCut(0),
fZetaOverlapCut(0),
fHistOn(0),
fhClustersDPhiAcc(0),
fhClustersDThetaAcc(0),
fhClustersDPhiAll(0),
fhClustersDThetaAll(0),
fhDPhiVsDThetaAll(0),
fhDPhiVsDThetaAcc(0),
fhetaTracklets(0),
fhphiTracklets(0),
fhetaClustersLay1(0),
fhphiClustersLay1(0){

  fNFiredChips[0] = 0;
  fNFiredChips[1] = 0;
  // Method to reconstruct the charged particles multiplicity with the 
  // SPD (tracklets).

  SetHistOn();

  if(AliITSReconstructor::GetRecoParam()) { 
    SetPhiWindow(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiWindow());
    SetThetaWindow(AliITSReconstructor::GetRecoParam()->GetTrackleterThetaWindow());
    SetPhiShift(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiShift());
    SetRemoveClustersFromOverlaps(AliITSReconstructor::GetRecoParam()->GetTrackleterRemoveClustersFromOverlaps());
    SetPhiOverlapCut(AliITSReconstructor::GetRecoParam()->GetTrackleterPhiOverlapCut());
    SetZetaOverlapCut(AliITSReconstructor::GetRecoParam()->GetTrackleterZetaOverlapCut());
  } else {
    SetPhiWindow();
    SetThetaWindow();
    SetPhiShift();
    SetRemoveClustersFromOverlaps();
    SetPhiOverlapCut();
    SetZetaOverlapCut();
  } 
  
  fClustersLay1              = 0;
  fClustersLay2              = 0;
  fDetectorIndexClustersLay1 = 0;
  fDetectorIndexClustersLay2 = 0;
  fOverlapFlagClustersLay1   = 0;
  fOverlapFlagClustersLay2   = 0;
  fTracklets                 = 0;
  fSClusters                 = 0;

  // definition of histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fhClustersDPhiAcc   = new TH1F("dphiacc",  "dphi",  100,-0.1,0.1);
  fhClustersDThetaAcc = new TH1F("dthetaacc","dtheta",100,-0.1,0.1);

  fhDPhiVsDThetaAcc = new TH2F("dphiVsDthetaAcc","",100,-0.1,0.1,100,-0.1,0.1);

  fhClustersDPhiAll   = new TH1F("dphiall",  "dphi",  100,0.0,0.5);
  fhClustersDThetaAll = new TH1F("dthetaall","dtheta",100,0.0,0.5);

  fhDPhiVsDThetaAll = new TH2F("dphiVsDthetaAll","",100,0.,0.5,100,0.,0.5);

  fhetaTracklets  = new TH1F("etaTracklets",  "eta",  100,-2.,2.);
  fhphiTracklets  = new TH1F("phiTracklets",  "phi",  100, 0., 2*TMath::Pi());
  fhetaClustersLay1  = new TH1F("etaClustersLay1",  "etaCl1",  100,-2.,2.);
  fhphiClustersLay1  = new TH1F("phiClustersLay1", "phiCl1", 100, 0., 2*TMath::Pi());
  
  TH1::AddDirectory(oldStatus);
}

//______________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor(const AliITSMultReconstructor &mr) : 
AliTrackleter(mr),
fDetTypeRec(0),fESDEvent(0),fTreeRP(0),fUsedClusLay1(0),fUsedClusLay2(0),
fClustersLay1(0),
fClustersLay2(0),
fDetectorIndexClustersLay1(0),
fDetectorIndexClustersLay2(0),
fOverlapFlagClustersLay1(0),
fOverlapFlagClustersLay2(0),
fTracklets(0),
fSClusters(0),
fNClustersLay1(0),
fNClustersLay2(0),
fNTracklets(0),
fNSingleCluster(0),
fPhiWindow(0),
fThetaWindow(0),
fPhiShift(0),
fRemoveClustersFromOverlaps(0),
fPhiOverlapCut(0),
fZetaOverlapCut(0),
fHistOn(0),
fhClustersDPhiAcc(0),
fhClustersDThetaAcc(0),
fhClustersDPhiAll(0),
fhClustersDThetaAll(0),
fhDPhiVsDThetaAll(0),
fhDPhiVsDThetaAcc(0),
fhetaTracklets(0),
fhphiTracklets(0),
fhetaClustersLay1(0),
fhphiClustersLay1(0)
 {
  // Copy constructor :!!! RS ATTENTION: old c-tor reassigned the pointers instead of creating a new copy -> would crash on delete
   AliError("May not use");
}

//______________________________________________________________________
AliITSMultReconstructor& AliITSMultReconstructor::operator=(const AliITSMultReconstructor& mr){
  // Assignment operator
  if (this != &mr) {
    this->~AliITSMultReconstructor();
    new(this) AliITSMultReconstructor(mr);
  }
  return *this;
}

//______________________________________________________________________
AliITSMultReconstructor::~AliITSMultReconstructor(){
  // Destructor

  // delete histograms
  delete fhClustersDPhiAcc;
  delete fhClustersDThetaAcc;
  delete fhClustersDPhiAll;
  delete fhClustersDThetaAll;
  delete fhDPhiVsDThetaAll;
  delete fhDPhiVsDThetaAcc;
  delete fhetaTracklets;
  delete fhphiTracklets;
  delete fhetaClustersLay1;
  delete fhphiClustersLay1;
  delete[] fUsedClusLay1;
  delete[] fUsedClusLay2;
  // delete arrays    
  for(Int_t i=0; i<fNTracklets; i++)
    delete [] fTracklets[i];
    
  for(Int_t i=0; i<fNSingleCluster; i++)
    delete [] fSClusters[i];
    
  delete [] fClustersLay1;
  delete [] fClustersLay2;
  delete [] fDetectorIndexClustersLay1;
  delete [] fDetectorIndexClustersLay2;
  delete [] fOverlapFlagClustersLay1;
  delete [] fOverlapFlagClustersLay2;
  delete [] fTracklets;
  delete [] fSClusters;
}

//____________________________________________________________________
void AliITSMultReconstructor::Reconstruct(AliESDEvent* esd, TTree* treeRP) 
{
  //
  // - calls LoadClusterArrays that finds the position of the clusters
  //   (in global coord) 
  // - convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). 
  // - makes an array of tracklets 
  //   
  // After this method has been called, the clusters of the two layers
  // and the tracklets can be retrieved by calling the Get'er methods.

  // reset counters
  if (fMult) delete fMult; fMult = 0;
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0;
  //
  if (!treeRP) { AliError(" Invalid ITS cluster tree !\n"); return; }
  //
  fESDEvent = esd;
  fTreeRP = treeRP;
  //
  // >>>> RS: this part is equivalent to former AliITSVertexer::FindMultiplicity
  //
  // see if there is a SPD vertex 
  Bool_t isVtxOK=kTRUE, isCosmics=kFALSE;
  AliESDVertex* vtx = (AliESDVertex*)fESDEvent->GetPrimaryVertexSPD();
  if (!vtx && vtx->GetNContributors()<0) isVtxOK = kFALSE;
  if (vtx && strstr(vtx->GetTitle(),"cosmics")) {
    isVtxOK = kFALSE;
    isCosmics = kTRUE;
  }
  //
  if (!isVtxOK) {
    if (!isCosmics) {
      AliDebug(1,"Tracklets multiplicity not determined because the primary vertex was not found");
      AliDebug(1,"Just counting the number of cluster-fired chips on the SPD layers");
    }
    vtx = 0;
  }
  float vtxf[3] = {vtx->GetX(),vtx->GetY(),vtx->GetZ()};
  FindTracklets(vtxf);
  //
  CreateMultiplicityObject();
}

//____________________________________________________________________
void AliITSMultReconstructor::Reconstruct(TTree* clusterTree, Float_t* vtx, Float_t* /* vtxRes*/) {
  //
  // RS NOTE - this is old reconstructor invocation, to be used from VertexFinder
  //
  // - calls LoadClusterArray that finds the position of the clusters
  //   (in global coord) 
  // - convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). 
  // - makes an array of tracklets 
  //   
  // After this method has been called, the clusters of the two layers
  // and the tracklets can be retrieved by calling the Get'er methods.
  if (fMult) delete fMult; fMult = 0;
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0;
  //
  if (!clusterTree) { AliError(" Invalid ITS cluster tree !\n"); return; }
  //
  fESDEvent = 0;
  fTreeRP = clusterTree;
  //
  FindTracklets(vtx);
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::FindTracklets(const Float_t *vtx) 
{
  // Find tracklets converging to vertex
  //
  LoadClusterArrays(fTreeRP);
  // flag clusters used by ESD tracks
  ProcessESDTracks();

  if (!vtx) return;

  const Double_t pi = TMath::Pi();
  
  // dPhi shift is field dependent
  // get average magnetic field
  Float_t bz = 0;
  AliMagF* field = 0;
  if (TGeoGlobalMagField::Instance()) field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
  if (!field)
  {
    AliError("Could not retrieve magnetic field. Assuming no field. Delta Phi shift will be deactivated in AliITSMultReconstructor.")
  }
  else
    bz = TMath::Abs(field->SolenoidField());
  
  const Double_t dPhiShift = fPhiShift / 5 * bz; 
  AliDebug(1, Form("Using phi shift of %f", dPhiShift));
  
  const Double_t dPhiWindow2 = fPhiWindow * fPhiWindow;
  const Double_t dThetaWindow2 = fThetaWindow * fThetaWindow;
  
  Int_t* partners = new Int_t[fNClustersLay2];
  Float_t* minDists = new Float_t[fNClustersLay2];
  Int_t* associatedLay1 = new Int_t[fNClustersLay1];
  TArrayI** blacklist = new TArrayI*[fNClustersLay1];

  for (Int_t i=0; i<fNClustersLay2; i++) {
    partners[i] = -1;
    minDists[i] = 2;
  }
  for (Int_t i=0; i<fNClustersLay1; i++)
    associatedLay1[i] = 0; 
  for (Int_t i=0; i<fNClustersLay1; i++)
    blacklist[i] = 0;

  // find the tracklets
  AliDebug(1,"Looking for tracklets... ");  
  
  //###########################################################
  // Loop on layer 1 : finding theta, phi and z 
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    
    float *clPar = GetClusterLayer1(iC1);
    Float_t x = clPar[kClTh] - vtx[0];
    Float_t y = clPar[kClPh] - vtx[1];
    Float_t z = clPar[kClZ]  - vtx[2];

    Float_t r    = TMath::Sqrt(x*x + y*y + z*z);
    
    clPar[kClTh] = TMath::ACos(z/r);                   // Store Theta
    clPar[kClPh] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi
    
    if (fHistOn) {
      Float_t eta = clPar[kClTh];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaClustersLay1->Fill(eta);    
      fhphiClustersLay1->Fill(clPar[kClPh]);
    }      
  }
  
  // Loop on layer 2 : finding theta, phi and r   
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {    
    float *clPar = GetClusterLayer2(iC2);
    Float_t x = clPar[kClTh] - vtx[0];
    Float_t y = clPar[kClPh] - vtx[1];
    Float_t z = clPar[kClZ]  - vtx[2];
   
    Float_t r    = TMath::Sqrt(x*x + y*y + z*z);

    clPar[kClTh] = TMath::ACos(z/r);                   // Store Theta
    clPar[kClPh] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi    
  }  
  
  //###########################################################
  Int_t found = 1;
  while (found > 0) {
    found = 0;

    // Step1: find all tracklets allowing double assocation
    // Loop on layer 1 
    for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    

      // already used or in the overlap ?
      if (associatedLay1[iC1] != 0 || fOverlapFlagClustersLay1[iC1]) continue;

      found++;

      // reset of variables for multiple candidates
      Int_t  iC2WithBestDist = -1;   // reset
      Double_t minDist       =  2;   // reset
      float* clPar1 = GetClusterLayer1(iC1);

      // Loop on layer 2 
      for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {      

        // in the overlap ?
        if (fOverlapFlagClustersLay2[iC2]) continue;
	float* clPar2 = GetClusterLayer2(iC2);

        if (blacklist[iC1]) {
          Bool_t blacklisted = kFALSE;
          for (Int_t i=blacklist[iC1]->GetSize(); i--;) {
            if (blacklist[iC1]->At(i) == iC2) {
              blacklisted = kTRUE;
              break;
            }
          }
          if (blacklisted) continue;
        }

	// find the difference in angles
	Double_t dTheta = TMath::Abs(clPar2[kClTh] - clPar1[kClTh]);
	Double_t dPhi   = TMath::Abs(clPar2[kClPh] - clPar1[kClPh]);
        // take into account boundary condition
        if (dPhi>pi) dPhi=2.*pi-dPhi;
        
 	if (fHistOn) {
	  fhClustersDPhiAll->Fill(dPhi);
	  fhClustersDThetaAll->Fill(dTheta);    
	  fhDPhiVsDThetaAll->Fill(dTheta, dPhi);
	}
        
        dPhi -= dPhiShift;
                
	// make "elliptical" cut in Phi and Theta! 
	Float_t d = dPhi*dPhi/dPhiWindow2 + dTheta*dTheta/dThetaWindow2;

	// look for the minimum distance: the minimum is in iC2WithBestDist
       	if (d<1 && d<minDist) {
	  minDist=d;
	  iC2WithBestDist = iC2;
	}
      } // end of loop over clusters in layer 2 
    
      if (minDist<1) { // This means that a cluster in layer 2 was found that matches with iC1

        if (minDists[iC2WithBestDist] > minDist) {
          Int_t oldPartner = partners[iC2WithBestDist];
          partners[iC2WithBestDist] = iC1;
          minDists[iC2WithBestDist] = minDist;

          // mark as assigned
          associatedLay1[iC1] = 1;
          
          if (oldPartner != -1) {
            // redo partner search for cluster in L0 (oldPartner), putting this one (iC2WithBestDist) on its blacklist
            if (blacklist[oldPartner] == 0) {
              blacklist[oldPartner] = new TArrayI(1);
            } else blacklist[oldPartner]->Set(blacklist[oldPartner]->GetSize()+1);

            blacklist[oldPartner]->AddAt(iC2WithBestDist, blacklist[oldPartner]->GetSize()-1);

            // mark as free
            associatedLay1[oldPartner] = 0;
          }
        } else {
          // try again to find a cluster without considering iC2WithBestDist 
          if (blacklist[iC1] == 0) {
            blacklist[iC1] = new TArrayI(1);
          }
          else 
            blacklist[iC1]->Set(blacklist[iC1]->GetSize()+1);
   
          blacklist[iC1]->AddAt(iC2WithBestDist, blacklist[iC1]->GetSize()-1);
        } 

      } else // cluster has no partner; remove
        associatedLay1[iC1] = 2;   
    } // end of loop over clusters in layer 1
  }  
 
  // Step2: store tracklets; remove used clusters
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {

    if (partners[iC2] == -1) continue;

    if (fRemoveClustersFromOverlaps) FlagClustersInOverlapRegions (partners[iC2],iC2);


    if (fOverlapFlagClustersLay1[partners[iC2]] || fOverlapFlagClustersLay2[iC2]) continue;

    float* clPar2 = GetClusterLayer2(iC2);
    float* clPar1 = GetClusterLayer1(partners[iC2]);

    Float_t* tracklet = fTracklets[fNTracklets] = new Float_t[kTrNPar]; // RS Add also the cluster id's
  
    // use the theta from the clusters in the first layer
    tracklet[kTrTheta] = clPar1[kClTh];
    // use the phi from the clusters in the first layer
    tracklet[kTrPhi] = clPar1[kClPh];
    // store the difference between phi1 and phi2
    tracklet[kTrDPhi] = clPar1[kClPh] - clPar2[kClPh];

    // define dphi in the range [0,pi] with proper sign (track charge correlated)
    if (tracklet[kTrDPhi] > TMath::Pi())   tracklet[kTrDPhi] = tracklet[kTrDPhi]-2.*TMath::Pi();
    if (tracklet[kTrDPhi] < -TMath::Pi())  tracklet[kTrDPhi] = tracklet[kTrDPhi]+2.*TMath::Pi();

    // store the difference between theta1 and theta2
    tracklet[kTrDTheta] = clPar1[kClTh] - clPar2[kClTh];

    if (fHistOn) {
      fhClustersDPhiAcc->Fill(tracklet[kTrDPhi]); 
      fhClustersDThetaAcc->Fill(tracklet[kTrDTheta]);    
      fhDPhiVsDThetaAcc->Fill(tracklet[kTrDTheta],tracklet[kTrDPhi]);
    }

    // find label
    // if equal label in both clusters found this label is assigned
    // if no equal label can be found the first labels of the L1 AND L2 cluster are assigned
    Int_t label1 = 0;
    Int_t label2 = 0;
    while (label2 < 3) {
      if ((Int_t) clPar1[kClMC0+label1] != -2 && (Int_t) clPar1[kClMC0+label1] == (Int_t) clPar2[kClMC0+label2])
        break;
      label1++;
      if (label1 == 3) {
        label1 = 0;
        label2++;
      }
    }
    if (label2 < 3) {
      AliDebug(AliLog::kDebug, Form("Found label %d == %d for tracklet candidate %d\n", (Int_t) clPar1[kClMC0+label1], (Int_t) clPar1[kClMC0+label2], fNTracklets));
      tracklet[kTrLab1] = clPar1[kClMC0+label1];
      tracklet[kTrLab2] = clPar2[kClMC0+label2];
    } else {
      AliDebug(AliLog::kDebug, Form("Did not find label %d %d %d %d %d %d for tracklet candidate %d\n", (Int_t) clPar1[kClMC0], (Int_t) clPar1[kClMC1], (Int_t) clPar1[kClMC2], (Int_t) clPar2[kClMC0], (Int_t) clPar2[kClMC1], (Int_t) clPar2[kClMC2], fNTracklets));
      tracklet[kTrLab1] = clPar1[kClMC0];
      tracklet[kTrLab2] = clPar2[kClMC0];
    }

    if (fHistOn) {
      Float_t eta = tracklet[kTrTheta];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaTracklets->Fill(eta);
      fhphiTracklets->Fill(tracklet[kTrPhi]);
    }
    //
    tracklet[kClID1] = partners[iC2];
    tracklet[kClID2] = iC2;
    //
    AliDebug(1,Form(" Adding tracklet candidate %d ", fNTracklets));
    AliDebug(1,Form(" Cl. %d of Layer 1 and %d of Layer 2", partners[iC2], iC2));
    fNTracklets++;

    associatedLay1[partners[iC2]] = 1;
  }
  
  // Delete the following else if you do not want to save Clusters! 
  // store the cluster
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {

    float* clPar1 = GetClusterLayer1(iC1);

    if (associatedLay1[iC1]==2||associatedLay1[iC1]==0) { 
      fSClusters[fNSingleCluster] = new Float_t[kClNPar];
      fSClusters[fNSingleCluster][kSCTh] = clPar1[kClTh];
      fSClusters[fNSingleCluster][kSCPh] = clPar1[kClPh];
      fSClusters[fNSingleCluster][kSCID] = iC1;
      AliDebug(1,Form(" Adding a single cluster %d (cluster %d  of layer 1)",
                fNSingleCluster, iC1));
      fNSingleCluster++;
    }
  }

  delete[] partners;
  delete[] minDists;

  for (Int_t i=0; i<fNClustersLay1; i++)
    if (blacklist[i])
      delete blacklist[i];
  delete[] blacklist;

  AliDebug(1,Form("%d tracklets found", fNTracklets));
}

//____________________________________________________________________
void AliITSMultReconstructor::CreateMultiplicityObject()
{
  // create AliMultiplicity object and store it in the ESD event
  //
  TBits fastOrFiredMap,firedChipMap;
  if (fDetTypeRec) {
   fastOrFiredMap  = fDetTypeRec->GetFastOrFiredMap();
   firedChipMap    = fDetTypeRec->GetFiredChipMap(fTreeRP);
  }
  //
  fMult = new AliMultiplicity(fNTracklets,fNSingleCluster,fNFiredChips[0],fNFiredChips[1],fastOrFiredMap);
  fMult->SetFiredChipMap(firedChipMap);
  AliITSRecPointContainer* rcont = AliITSRecPointContainer::Instance();
  fMult->SetITSClusters(0,rcont->GetNClustersInLayer(1,fTreeRP));
  for(Int_t kk=2;kk<=6;kk++) fMult->SetITSClusters(kk-1,rcont->GetNClustersInLayerFast(kk));
  //
  for (int i=fNTracklets;i--;)  {
    float* tlInfo = fTracklets[i];
    fMult->SetTrackletData(i,tlInfo, fUsedClusLay1[int(tlInfo[kClID1])]|fUsedClusLay2[int(tlInfo[kClID2])]);
  }
  //  
  for (int i=fNSingleCluster;i--;) {
    float* clInfo = fSClusters[i];
    fMult->SetSingleClusterData(i,clInfo,fUsedClusLay1[int(clInfo[kSCID])]);
  }
  fMult->CompactBits();
  //
}


//____________________________________________________________________
void AliITSMultReconstructor::LoadClusterArrays(TTree* itsClusterTree) 
{
  // This method
  // - gets the clusters from the cluster tree 
  // - convert them into global coordinates 
  // - store them in the internal arrays
  // - count the number of cluster-fired chips
  //
  // RS: This method was strongly modified wrt original by Jan Fiete. In order to have the same numbering
  // of clusters as in the ITS reco I had to introduce sorting in Z
  // Also note that now the clusters data are stored not in float[6] attached to float**, but in 1-D array
  
  AliDebug(1,"Loading clusters and cluster-fired chips ...");
  
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNFiredChips[0] = 0;
  fNFiredChips[1] = 0;
  
  AliITSsegmentationSPD seg;

  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  TClonesArray* itsClusters=rpcont->FetchClusters(0,itsClusterTree);
  if(!rpcont->IsSPDActive()){
    AliWarning("No SPD rec points found, multiplicity not calculated");
    return;
  } 
  //
  // count clusters
  // loop over the SPD subdetectors
  TObjArray clArr(100);
  for (int il=0;il<2;il++) {
    int nclLayer = 0;
    int detMin = AliITSgeomTGeo::GetModuleIndex(il+1,1,1);
    int detMax = AliITSgeomTGeo::GetModuleIndex(il+2,1,1);
    for (int idt=detMin;idt<detMax;idt++) {
      itsClusters=rpcont->UncheckedGetClusters(idt);
      int nClusters = itsClusters->GetEntriesFast();
      if (!nClusters) continue;
      Int_t nClustersInChip[5] = {0,0,0,0,0};
      while(nClusters--) {
	AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
	if (!cluster) continue;
	clArr[nclLayer++] = cluster;
	nClustersInChip[ seg.GetChipFromLocal(0,cluster->GetDetLocalZ()) ]++; 
      }
      for(Int_t ifChip=5;ifChip--;) if (nClustersInChip[ifChip]) fNFiredChips[il]++;
    }
    // sort the clusters in Z (to have the same numbering as in ITS reco
    Float_t *z    = new Float_t[nclLayer];
    Int_t * index = new Int_t[nclLayer];
    for (int ic=0;ic<nclLayer;ic++) z[ic] = ((AliITSRecPoint*)clArr[ic])->GetZ();
    TMath::Sort(nclLayer,z,index,kFALSE);
    Float_t*   clustersLay              = new Float_t[nclLayer*kClNPar];
    Int_t*     detectorIndexClustersLay = new Int_t[nclLayer];
    Bool_t*    overlapFlagClustersLay   = new Bool_t[nclLayer];
    Char_t*    usedClusLay              = new Char_t[nclLayer];
    //
    for (int ic=0;ic<nclLayer;ic++) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)clArr[index[ic]];
      float* clPar = &clustersLay[ic*kClNPar];
      //      
      cluster->GetGlobalXYZ( clPar );
      detectorIndexClustersLay[ic] = cluster->GetDetectorIndex(); 
      overlapFlagClustersLay[ic]   = kFALSE;
      usedClusLay[ic]              = 0;
      for (Int_t i=3;i--;) clPar[kClMC0+i] = cluster->GetLabel(i);
    }
    clArr.Clear();
    delete[] z;
    delete[] index;
    //
    if (il==0) {
      fClustersLay1              = clustersLay;
      fOverlapFlagClustersLay1   = overlapFlagClustersLay;
      fDetectorIndexClustersLay1 = detectorIndexClustersLay;
      fUsedClusLay1              = usedClusLay;
      fNClustersLay1             = nclLayer;
    }
    else {
      fClustersLay2              = clustersLay;
      fOverlapFlagClustersLay2   = overlapFlagClustersLay;
      fDetectorIndexClustersLay2 = detectorIndexClustersLay;
      fUsedClusLay2              = usedClusLay;
      fNClustersLay2             = nclLayer;
    }
  }
  //
  // no double association allowed
  int nmaxT                  = TMath::Min(fNClustersLay1, fNClustersLay2);
  fTracklets                 = new Float_t*[nmaxT];
  fSClusters                 = new Float_t*[fNClustersLay1]; 
  for (Int_t i=nmaxT;i--;) fTracklets[i] = 0;
  //
  AliDebug(1,Form("(clusters in layer 1 : %d,  layer 2: %d)",fNClustersLay1,fNClustersLay2));
  AliDebug(1,Form("(cluster-fired chips in layer 1 : %d,  layer 2: %d)",fNFiredChips[0],fNFiredChips[1]));
}
//____________________________________________________________________
void
AliITSMultReconstructor::LoadClusterFiredChips(TTree* itsClusterTree) {
  // This method
  // - gets the clusters from the cluster tree 
  // - counts the number of (cluster)fired chips
  
  AliDebug(1,"Loading cluster-fired chips ...");
  
  fNFiredChips[0] = 0;
  fNFiredChips[1] = 0;
  
  AliITSsegmentationSPD seg;
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  TClonesArray* itsClusters=rpcont->FetchClusters(0,itsClusterTree);
  if(!rpcont->IsSPDActive()){
    AliWarning("No SPD rec points found, multiplicity not calculated");
    return;
  } 

  // loop over the its subdetectors
  Int_t nSPDmodules=AliITSgeomTGeo::GetModuleIndex(3,1,1);
  for (Int_t iIts=0; iIts < nSPDmodules; iIts++) {
    itsClusters=rpcont->UncheckedGetClusters(iIts);
    Int_t nClusters = itsClusters->GetEntriesFast();

    // number of clusters in each chip of the current module
    Int_t nClustersInChip[5] = {0,0,0,0,0};
    Int_t layer = 0;
    
    // loop over clusters
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
      
      layer = cluster->GetLayer();
      if (layer>1) continue;            

      // find the chip for the current cluster
      Float_t locz = cluster->GetDetLocalZ();
      Int_t iChip = seg.GetChipFromLocal(0,locz);
      nClustersInChip[iChip]++; 
      
    }// end of cluster loop

    // get number of fired chips in the current module
    for(Int_t ifChip=0; ifChip<5; ifChip++) {
      if(nClustersInChip[ifChip] >= 1)  fNFiredChips[layer]++;
    }

  } // end of its "subdetector" loop  
  

  AliDebug(1,Form("(cluster-fired chips in layer 1 : %d,  layer 2: %d)",fNFiredChips[0],fNFiredChips[1]));
}
//____________________________________________________________________
void
AliITSMultReconstructor::SaveHists() {
  // This method save the histograms on the output file
  // (only if fHistOn is TRUE). 
  
  if (!fHistOn)
    return;

  fhClustersDPhiAll->Write();
  fhClustersDThetaAll->Write();
  fhDPhiVsDThetaAll->Write();

  fhClustersDPhiAcc->Write();
  fhClustersDThetaAcc->Write();
  fhDPhiVsDThetaAcc->Write();

  fhetaTracklets->Write();
  fhphiTracklets->Write();
  fhetaClustersLay1->Write();
  fhphiClustersLay1->Write();
}

//____________________________________________________________________
void 
AliITSMultReconstructor::FlagClustersInOverlapRegions (Int_t iC1, Int_t iC2WithBestDist) {

  Float_t distClSameMod=0.;
  Float_t distClSameModMin=0.;
  Int_t   iClOverlap =0;
  Float_t meanRadiusLay1 = 3.99335; // average radius inner layer
  Float_t meanRadiusLay2 = 7.37935; // average radius outer layer;

  Float_t zproj1=0.;
  Float_t zproj2=0.;
  Float_t deZproj=0.;
  Float_t* clPar1  = GetClusterLayer1(iC1);
  Float_t* clPar2B = GetClusterLayer2(iC2WithBestDist);
  // Loop on inner layer clusters
  for (Int_t iiC1=0; iiC1<fNClustersLay1; iiC1++) {
    if (!fOverlapFlagClustersLay1[iiC1]) {
      // only for adjacent modules
      if ((TMath::Abs(fDetectorIndexClustersLay1[iC1]-fDetectorIndexClustersLay1[iiC1])==4)||
         (TMath::Abs(fDetectorIndexClustersLay1[iC1]-fDetectorIndexClustersLay1[iiC1])==76)) {
	Float_t *clPar11 = GetClusterLayer1(iiC1);
        Float_t dePhi=TMath::Abs(clPar11[kClPh]-clPar1[kClPh]);
        if (dePhi>TMath::Pi()) dePhi=2.*TMath::Pi()-dePhi;

        zproj1=meanRadiusLay1/TMath::Tan(clPar1[kClTh]);
        zproj2=meanRadiusLay1/TMath::Tan(clPar11[kClTh]);

        deZproj=TMath::Abs(zproj1-zproj2);

        distClSameMod = TMath::Sqrt(TMath::Power(deZproj/fZetaOverlapCut,2)+TMath::Power(dePhi/fPhiOverlapCut,2));
        if (distClSameMod<=1.) fOverlapFlagClustersLay1[iiC1]=kTRUE;

//        if (distClSameMod<=1.) {
//          if (distClSameModMin==0. || distClSameMod<distClSameModMin) {
//            distClSameModMin=distClSameMod;
//            iClOverlap=iiC1;
//          } 
//        }


      } // end adjacent modules
    } 
  } // end Loop on inner layer clusters

//  if (distClSameModMin!=0.) fOverlapFlagClustersLay1[iClOverlap]=kTRUE;

  distClSameMod=0.;
  distClSameModMin=0.;
  iClOverlap =0;
  // Loop on outer layer clusters
  for (Int_t iiC2=0; iiC2<fNClustersLay2; iiC2++) {
    if (!fOverlapFlagClustersLay2[iiC2]) {
      // only for adjacent modules
      Float_t *clPar2 = GetClusterLayer2(iiC2);
      if ((TMath::Abs(fDetectorIndexClustersLay2[iC2WithBestDist]-fDetectorIndexClustersLay2[iiC2])==4) ||
         (TMath::Abs(fDetectorIndexClustersLay2[iC2WithBestDist]-fDetectorIndexClustersLay2[iiC2])==156)) {
        Float_t dePhi=TMath::Abs(clPar2[kClPh]-clPar2B[kClPh]);
        if (dePhi>TMath::Pi()) dePhi=2.*TMath::Pi()-dePhi;

        zproj1=meanRadiusLay2/TMath::Tan(clPar2B[kClTh]);
        zproj2=meanRadiusLay2/TMath::Tan(clPar2[kClTh]);

        deZproj=TMath::Abs(zproj1-zproj2);
        distClSameMod = TMath::Sqrt(TMath::Power(deZproj/fZetaOverlapCut,2)+TMath::Power(dePhi/fPhiOverlapCut,2));
        if (distClSameMod<=1.) fOverlapFlagClustersLay2[iiC2]=kTRUE;

//        if (distClSameMod<=1.) {
//          if (distClSameModMin==0. || distClSameMod<distClSameModMin) {
//            distClSameModMin=distClSameMod;
//            iClOverlap=iiC2;
//          }
//        }

      } // end adjacent modules
    }
  } // end Loop on outer layer clusters

//  if (distClSameModMin!=0.) fOverlapFlagClustersLay2[iClOverlap]=kTRUE;

}

//____________________________________________________________________
void AliITSMultReconstructor::ProcessESDTracks()
{
  // Flag the clusters used by ESD tracks
  // Flag primary tracks to be used for multiplicity counting 
  //
  AliESDVertex* vtx = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();
  if (!vtx) vtx = (AliESDVertex*)fESDEvent->GetPrimaryVertexSPD();
  if (!vtx) {
    AliDebug(1,"No primary vertex: cannot flag primary tracks");
    return;
  }
  Int_t ntracks = fESDEvent->GetNumberOfTracks();
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack* track = fESDEvent->GetTrack(itr);
    if (!track->IsOn(AliESDtrack::kITSin)) continue; // use only tracks propagated in ITS to vtx
    FlagTrackClusters(track);
    FlagIfPrimary(track,vtx);
  }
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::FlagTrackClusters(const AliESDtrack* track)
{
  // RS: flag the SPD clusters of the track if it is useful for the multiplicity estimation
  //
  Int_t idx[12];
  if ( track->GetITSclusters(idx)<3 ) return; // at least 3 clusters must be used in the fit
  //
  char mark = track->IsOn(AliESDtrack::kITSpureSA) ? kITSSAPBit : kITSTPCBit;
  char *uClus[2] = {fUsedClusLay1,fUsedClusLay2};
  for (int i=AliESDfriendTrack::kMaxITScluster;i--;) {
    // note: i>=6 is for extra clusters
    if (idx[i]<0) continue;
    int layID= (idx[i] & 0xf0000000) >> 28; 
    if (layID>1) continue; // SPD only
    int clID = (idx[i] & 0x0fffffff);
    uClus[layID][clID] |= mark;
  }
  //
}

//____________________________________________________________________
void AliITSMultReconstructor::FlagIfPrimary(AliESDtrack* track, const AliVertex* vtx)
{
  // RS: check if the track is primary and set the flag
  const double kPDCASPD1 = 0.1;
  const double kPDCASPD0 = 0.3;
  //
  double cut = (track->HasPointOnITSLayer(0)||track->HasPointOnITSLayer(1)) ? kPDCASPD1 : kPDCASPD0;
  // in principle, the track must already have been propagated to vertex
  /*
  Double_t dzRec[2]={0,0}, covdzRec[3];
  track->PropagateToDCA(vtx, fESDEvent->GetMagneticField(), 3.0, dzRec, covdzRec);
  */
  double dist = track->GetD(vtx->GetX(),vtx->GetY(),fESDEvent->GetMagneticField());
  if (TMath::Abs(dist*track->P())<cut) track->SetStatus(AliESDtrack::kMultPrimary);
}
