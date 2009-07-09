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
//_________________________________________________________________________

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include "TArrayI.h"

#include "AliITSMultReconstructor.h"
#include "AliITSReconstructor.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliLog.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"

//____________________________________________________________________
ClassImp(AliITSMultReconstructor)


//____________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor():
TObject(),
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
AliITSMultReconstructor::AliITSMultReconstructor(const AliITSMultReconstructor &mr) : TObject(mr),
fClustersLay1(mr.fClustersLay1),
fClustersLay2(mr.fClustersLay2),
fDetectorIndexClustersLay1(mr.fDetectorIndexClustersLay1),
fDetectorIndexClustersLay2(mr.fDetectorIndexClustersLay2),
fOverlapFlagClustersLay1(mr.fOverlapFlagClustersLay1),
fOverlapFlagClustersLay2(mr.fOverlapFlagClustersLay2),
fTracklets(mr.fTracklets),
fSClusters(mr.fSClusters),
fNClustersLay1(mr.fNClustersLay1),
fNClustersLay2(mr.fNClustersLay2),
fNTracklets(mr.fNTracklets),
fNSingleCluster(mr.fNSingleCluster),
fPhiWindow(mr.fPhiWindow),
fThetaWindow(mr.fThetaWindow),
fPhiShift(mr.fPhiShift),
fRemoveClustersFromOverlaps(mr.fRemoveClustersFromOverlaps),
fPhiOverlapCut(mr.fPhiOverlapCut),
fZetaOverlapCut(mr.fZetaOverlapCut),
fHistOn(mr.fHistOn),
fhClustersDPhiAcc(mr.fhClustersDPhiAcc),
fhClustersDThetaAcc(mr.fhClustersDThetaAcc),
fhClustersDPhiAll(mr.fhClustersDPhiAll),
fhClustersDThetaAll(mr.fhClustersDThetaAll),
fhDPhiVsDThetaAll(mr.fhDPhiVsDThetaAll),
fhDPhiVsDThetaAcc(mr.fhDPhiVsDThetaAcc),
fhetaTracklets(mr.fhetaTracklets),
fhphiTracklets(mr.fhphiTracklets),
fhetaClustersLay1(mr.fhetaClustersLay1),
fhphiClustersLay1(mr.fhphiClustersLay1) {
  // Copy constructor

}

//______________________________________________________________________
AliITSMultReconstructor& AliITSMultReconstructor::operator=(const AliITSMultReconstructor& mr){
  // Assignment operator
  this->~AliITSMultReconstructor();
  new(this) AliITSMultReconstructor(mr);
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

  // delete arrays
  for(Int_t i=0; i<fNClustersLay1; i++)
    delete [] fClustersLay1[i];
    
  for(Int_t i=0; i<fNClustersLay2; i++)
    delete [] fClustersLay2[i];
    
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
void AliITSMultReconstructor::Reconstruct(TTree* clusterTree, Float_t* vtx, Float_t* /* vtxRes*/) {
  //
  // - calls LoadClusterArray that finds the position of the clusters
  //   (in global coord) 
  // - convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). 
  // - makes an array of tracklets 
  //   
  // After this method has been called, the clusters of the two layers
  // and the tracklets can be retrieved by calling the Get'er methods.

  // reset counters
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNTracklets = 0; 
  fNSingleCluster = 0;

  // loading the clusters 
  LoadClusterArrays(clusterTree);

  const Double_t pi = TMath::Pi();
  
  // dPhi shift is field dependent
  // get average magnetic field
  Float_t bz = 0;
  AliMagF* field = 0;
  if (TGeoGlobalMagField::Instance())
    field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
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
    Float_t x = fClustersLay1[iC1][0] - vtx[0];
    Float_t y = fClustersLay1[iC1][1] - vtx[1];
    Float_t z = fClustersLay1[iC1][2] - vtx[2];

    Float_t r    = TMath::Sqrt(x*x + y*y + z*z);
    
    fClustersLay1[iC1][0] = TMath::ACos(z/r);                   // Store Theta
    fClustersLay1[iC1][1] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi
    
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
   
    Float_t r    = TMath::Sqrt(x*x + y*y + z*z);
    
    fClustersLay2[iC2][0] = TMath::ACos(z/r);                   // Store Theta
    fClustersLay2[iC2][1] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi
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
     
      // Loop on layer 2 
      for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {      

        // in the overlap ?
        if (fOverlapFlagClustersLay2[iC2]) continue;

        if (blacklist[iC1]) {
          Bool_t blacklisted = kFALSE;
          for (Int_t i=0; i<blacklist[iC1]->GetSize(); i++) {
            if (blacklist[iC1]->At(i) == iC2) {
              blacklisted = kTRUE;
              break;
            }
          }
          if (blacklisted) continue;
        }

	// find the difference in angles
	Double_t dTheta = TMath::Abs(fClustersLay2[iC2][0] - fClustersLay1[iC1][0]);
	Double_t dPhi  = TMath::Abs(fClustersLay2[iC2][1] - fClustersLay1[iC1][1]);
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

    if (fOverlapFlagClustersLay1[partners[iC2]] || fOverlapFlagClustersLay2[iC2]) continue;
    if (fRemoveClustersFromOverlaps) FlagClustersInOverlapRegions (partners[iC2],iC2);

    fTracklets[fNTracklets] = new Float_t[6];
  
    // use the theta from the clusters in the first layer
    fTracklets[fNTracklets][0] = fClustersLay1[partners[iC2]][0];
    // use the phi from the clusters in the first layer
    fTracklets[fNTracklets][1] = fClustersLay1[partners[iC2]][1];
    // store the difference between phi1 and phi2
    fTracklets[fNTracklets][2] = fClustersLay1[partners[iC2]][1] - fClustersLay2[iC2][1];

    // define dphi in the range [0,pi] with proper sign (track charge correlated)
    if (fTracklets[fNTracklets][2] > TMath::Pi())
      fTracklets[fNTracklets][2] = fTracklets[fNTracklets][2]-2.*TMath::Pi();
    if (fTracklets[fNTracklets][2] < -TMath::Pi())
      fTracklets[fNTracklets][2] = fTracklets[fNTracklets][2]+2.*TMath::Pi();

    // store the difference between theta1 and theta2
    fTracklets[fNTracklets][3] = fClustersLay1[partners[iC2]][0] - fClustersLay2[iC2][0];

    if (fHistOn) {
      fhClustersDPhiAcc->Fill(fTracklets[fNTracklets][2]); 
      fhClustersDThetaAcc->Fill(fTracklets[fNTracklets][3]);    
      fhDPhiVsDThetaAcc->Fill(fTracklets[fNTracklets][3],fTracklets[fNTracklets][2]);
    }

    // find label
    // if equal label in both clusters found this label is assigned
    // if no equal label can be found the first labels of the L1 AND L2 cluster are assigned
    Int_t label1 = 0;
    Int_t label2 = 0;
    while (label2 < 3) {
      if ((Int_t) fClustersLay1[partners[iC2]][3+label1] != -2 && (Int_t) fClustersLay1[partners[iC2]][3+label1] == (Int_t) fClustersLay2[iC2][3+label2])
        break;
      label1++;
      if (label1 == 3) {
        label1 = 0;
        label2++;
      }
    }
    if (label2 < 3) {
      AliDebug(AliLog::kDebug, Form("Found label %d == %d for tracklet candidate %d\n", (Int_t) fClustersLay1[partners[iC2]][3+label1], (Int_t) fClustersLay2[iC2][3+label2], fNTracklets));
      fTracklets[fNTracklets][4] = fClustersLay1[partners[iC2]][3+label1];
      fTracklets[fNTracklets][5] = fClustersLay2[iC2][3+label2];
    } else {
      AliDebug(AliLog::kDebug, Form("Did not find label %d %d %d %d %d %d for tracklet candidate %d\n", (Int_t) fClustersLay1[partners[iC2]][3], (Int_t) fClustersLay1[partners[iC2]][4], (Int_t) fClustersLay1[partners[iC2]][5], (Int_t) fClustersLay2[iC2][3], (Int_t) fClustersLay2[iC2][4], (Int_t) fClustersLay2[iC2][5], fNTracklets));
      fTracklets[fNTracklets][4] = fClustersLay1[partners[iC2]][3];
      fTracklets[fNTracklets][5] = fClustersLay2[iC2][3];
    }

    if (fHistOn) {
      Float_t eta=fTracklets[fNTracklets][0];
      eta= TMath::Tan(eta/2.);
      eta=-TMath::Log(eta);
      fhetaTracklets->Fill(eta);
      fhphiTracklets->Fill(fTracklets[fNTracklets][1]);
    }

    AliDebug(1,Form(" Adding tracklet candidate %d ", fNTracklets));
    AliDebug(1,Form(" Cl. %d of Layer 1 and %d of Layer 2", partners[iC2], iC2));
    fNTracklets++;

    associatedLay1[partners[iC2]] = 1;
  }
  
  // Delete the following else if you do not want to save Clusters! 
  // store the cluster
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {
    if (associatedLay1[iC1]==2||associatedLay1[iC1]==0) { 
      fSClusters[fNSingleCluster] = new Float_t[2];
      fSClusters[fNSingleCluster][0] = fClustersLay1[iC1][0];
      fSClusters[fNSingleCluster][1] = fClustersLay1[iC1][1];
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
void
AliITSMultReconstructor::LoadClusterArrays(TTree* itsClusterTree) {
  // This method
  // - gets the clusters from the cluster tree 
  // - convert them into global coordinates 
  // - store them in the internal arrays
  // - count the number of cluster-fired chips
  
  AliDebug(1,"Loading clusters and cluster-fired chips ...");
  
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNFiredChips[0] = 0;
  fNFiredChips[1] = 0;
  
  AliITSsegmentationSPD *seg = new AliITSsegmentationSPD();

  TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
  TBranch* itsClusterBranch=itsClusterTree->GetBranch("ITSRecPoints");

  itsClusterBranch->SetAddress(&itsClusters);

  Int_t nItsSubs = (Int_t)itsClusterTree->GetEntries();  
  Float_t cluGlo[3]={0.,0.,0.};
 

  // count clusters
  // loop over the its subdetectors
  for (Int_t iIts=0; iIts < nItsSubs; iIts++) {
    if (!itsClusterTree->GetEvent(iIts)) 
      continue;
    
    Int_t nClusters = itsClusters->GetEntriesFast();
    // loop over clusters
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
      
      Int_t layer = cluster->GetLayer();
      if (layer == 0)
        fNClustersLay1++;
      else if (layer == 1)
        fNClustersLay2++;
    }
  }
  
  // create arrays
  fClustersLay1              = new Float_t*[fNClustersLay1];
  fDetectorIndexClustersLay1 = new Int_t[fNClustersLay1];
  fOverlapFlagClustersLay1   = new Bool_t[fNClustersLay1];
  
  fClustersLay2              = new Float_t*[fNClustersLay2];
  fDetectorIndexClustersLay2 = new Int_t[fNClustersLay2];
  fOverlapFlagClustersLay2   = new Bool_t[fNClustersLay2];
  
  // no double association allowed
  fTracklets                 = new Float_t*[TMath::Min(fNClustersLay1, fNClustersLay2)];
  fSClusters                 = new Float_t*[fNClustersLay1];
  
  for (Int_t i=0; i<fNClustersLay1; i++) {
    fClustersLay1[i]       = new Float_t[6];
    fOverlapFlagClustersLay1[i]   = kFALSE;
    fSClusters[i] = 0;
  } 
 
  for (Int_t i=0; i<fNClustersLay2; i++) {
    fClustersLay2[i]       = new Float_t[6];
    fOverlapFlagClustersLay2[i]   = kFALSE;
  } 
 
  for (Int_t i=0; i<TMath::Min(fNClustersLay1, fNClustersLay2); i++)
    fTracklets[i] = 0;
    
  // fill clusters
  // loop over the its subdetectors
  fNClustersLay1 = 0; // reset to 0
  fNClustersLay2 = 0;
  for (Int_t iIts=0; iIts < nItsSubs; iIts++) {
    
    if (!itsClusterTree->GetEvent(iIts)) 
      continue;
    
    Int_t nClusters = itsClusters->GetEntriesFast();

    // number of clusters in each chip of the current module
    Int_t nClustersInChip[5] = {0,0,0,0,0};
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

      // find the chip for the current cluster
      Float_t locz = cluster->GetDetLocalZ();
      Int_t iChip = seg->GetChipFromLocal(0,locz);
      nClustersInChip[iChip]++; 
      
      if (layer==0) {
	fClustersLay1[fNClustersLay1][0] = x;
	fClustersLay1[fNClustersLay1][1] = y;
	fClustersLay1[fNClustersLay1][2] = z;
 
        fDetectorIndexClustersLay1[fNClustersLay1]=iIts;  

	for (Int_t i=0; i<3; i++)
		fClustersLay1[fNClustersLay1][3+i] = cluster->GetLabel(i);
	fNClustersLay1++;
      }
      if (layer==1) {
	fClustersLay2[fNClustersLay2][0] = x;
	fClustersLay2[fNClustersLay2][1] = y;
	fClustersLay2[fNClustersLay2][2] = z;

        fDetectorIndexClustersLay2[fNClustersLay2]=iIts;
 
	for (Int_t i=0; i<3; i++)
		fClustersLay2[fNClustersLay2][3+i] = cluster->GetLabel(i);
	fNClustersLay2++;
      }
      
    }// end of cluster loop

    // get number of fired chips in the current module
    if(layer<2)
    for(Int_t ifChip=0; ifChip<5; ifChip++) {
      if(nClustersInChip[ifChip] >= 1)  fNFiredChips[layer]++;
    }

  } // end of its "subdetector" loop  

  if (itsClusters) {
    itsClusters->Delete();
    delete itsClusters;
    delete seg;
    itsClusters = 0;
  }
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
  
  AliITSsegmentationSPD *seg = new AliITSsegmentationSPD();

  TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
  TBranch* itsClusterBranch=itsClusterTree->GetBranch("ITSRecPoints");

  itsClusterBranch->SetAddress(&itsClusters);

  Int_t nItsSubs = (Int_t)itsClusterTree->GetEntries();  
 
  // loop over the its subdetectors
  for (Int_t iIts=0; iIts < nItsSubs; iIts++) {
    
    if (!itsClusterTree->GetEvent(iIts)) 
      continue;
    
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
      Int_t iChip = seg->GetChipFromLocal(0,locz);
      nClustersInChip[iChip]++; 
      
    }// end of cluster loop

    // get number of fired chips in the current module
    if(layer<2)
    for(Int_t ifChip=0; ifChip<5; ifChip++) {
      if(nClustersInChip[ifChip] >= 1)  fNFiredChips[layer]++;
    }

  } // end of its "subdetector" loop  
  
  if (itsClusters) {
    itsClusters->Delete();
    delete itsClusters;
    delete seg;
    itsClusters = 0;
  }
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

  // Loop on inner layer clusters
  for (Int_t iiC1=0; iiC1<fNClustersLay1; iiC1++) {
    if (!fOverlapFlagClustersLay1[iiC1]) {
      // only for adjacent modules
      if ((TMath::Abs(fDetectorIndexClustersLay1[iC1]-fDetectorIndexClustersLay1[iiC1])==4)||
         (TMath::Abs(fDetectorIndexClustersLay1[iC1]-fDetectorIndexClustersLay1[iiC1])==76)) {
        Float_t dePhi=TMath::Abs(fClustersLay1[iiC1][1]-fClustersLay1[iC1][1]);
        if (dePhi>TMath::Pi()) dePhi=2.*TMath::Pi()-dePhi;

        zproj1=meanRadiusLay1/TMath::Tan(fClustersLay1[iC1][0]);
        zproj2=meanRadiusLay1/TMath::Tan(fClustersLay1[iiC1][0]);

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
      if ((TMath::Abs(fDetectorIndexClustersLay2[iC2WithBestDist]-fDetectorIndexClustersLay2[iiC2])==4) ||
         (TMath::Abs(fDetectorIndexClustersLay2[iC2WithBestDist]-fDetectorIndexClustersLay2[iiC2])==156)) {
        Float_t dePhi=TMath::Abs(fClustersLay2[iiC2][1]-fClustersLay2[iC2WithBestDist][1]);
        if (dePhi>TMath::Pi()) dePhi=2.*TMath::Pi()-dePhi;

        zproj1=meanRadiusLay2/TMath::Tan(fClustersLay2[iC2WithBestDist][0]);
        zproj2=meanRadiusLay2/TMath::Tan(fClustersLay2[iiC2][0]);

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
