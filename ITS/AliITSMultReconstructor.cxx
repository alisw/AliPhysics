//____________________________________________________________________
// 
// AliITSMultReconstructor - find clusters in the pixels (theta and
// phi) and tracklets.
// 
// These can be used to extract charged particles multiplcicity from the ITS.
//
// A tracklet consist of two ITS clusters, one in the first pixel
// layer and one in the second. The clusters are associates if the 
// differencies in Phi (azimuth) and Zeta (longitudinal) are inside 
// a fiducial volume. In case of multiple candidates it is selected the
// candidate with minimum distance in Phi. 
// The parameter AssociationChoice allows to control if two clusters 
// in layer 2 can be associated to the same cluster in layer 1 or not.
//
// -----------------------------------------------------------------
// 
// TODO: 
// 
// - Introduce a rough pt estimation from the difference in phi ? 
// - Allow for a more refined selection criterium in case of multiple 
//   candidates (for instance by introducing weights for the difference 
//   in Phi and Zeta). 
//
//____________________________________________________________________

#include "AliITSMultReconstructor.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"


#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliLog.h"

//____________________________________________________________________
ClassImp(AliITSMultReconstructor)

//____________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor() {

  fGeometry =0;

  SetHistOn();
  SetPhiWindow();
  SetZetaWindow();
  SetOnlyOneTrackletPerC2();

  fClustersLay1       = new Float_t*[300000];
  fClustersLay2       = new Float_t*[300000];
  fTracklets          = new Float_t*[300000];
  fAssociationFlag    = new Bool_t[300000];

  for(Int_t i=0; i<300000; i++) {
    fClustersLay1[i]       = new Float_t[3];
    fClustersLay2[i]       = new Float_t[3];
    fTracklets[i]          = new Float_t[3];
    fAssociationFlag[i]    = kFALSE;
  }

  // definition of histograms
  fhClustersDPhi   = new TH1F("dphi",  "dphi",  200,-0.1,0.1);
  fhClustersDPhi->SetDirectory(0);
  fhClustersDTheta = new TH1F("dtheta","dtheta",200,-0.1,0.1);
  fhClustersDTheta->SetDirectory(0);
  fhClustersDZeta = new TH1F("dzeta","dzeta",200,-0.2,0.2);
  fhClustersDZeta->SetDirectory(0);

  fhDPhiVsDThetaAll = new TH2F("dphiVsDthetaAll","",200,-0.1,0.1,200,-0.1,0.1);
  fhDPhiVsDThetaAll->SetDirectory(0);
  fhDPhiVsDThetaAcc = new TH2F("dphiVsDthetaAcc","",200,-0.1,0.1,200,-0.1,0.1);
  fhDPhiVsDThetaAcc->SetDirectory(0);

}


//____________________________________________________________________
void
AliITSMultReconstructor::Reconstruct(TTree* clusterTree, Float_t* vtx, Float_t* /* vtxRes*/) {
  //
  // - calls LoadClusterArray that finds the position of the clusters
  //   (in global coord) 
  // - convert the cluster coordinates to theta, phi (seen from the
  //   interaction vertex). The third coordinate is used for ....
  // - makes an array of tracklets 
  //   
  // After this method has been called, the clusters of the two layers
  // and the tracklets can be retrieved by calling the Get'er methods.

  // reset counters
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  fNTracklets = 0; 

  // loading the clusters 
  LoadClusterArrays(clusterTree);
  
  // find the tracklets
  AliDebug(1,"Looking for tracklets... ");  

  //###########################################################
  // Loop on layer 1 : finding theta, phi and z 
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    
    Float_t x = fClustersLay1[iC1][0] - vtx[0];
    Float_t y = fClustersLay1[iC1][1] - vtx[1];
    Float_t z = fClustersLay1[iC1][2] - vtx[2];
    
    Float_t r    = TMath::Sqrt(TMath::Power(x,2) +
			       TMath::Power(y,2) +
			       TMath::Power(z,2));
    
    fClustersLay1[iC1][0] = TMath::ACos(z/r);  // Store Theta
    fClustersLay1[iC1][1] = TMath::ATan(y/x);  // Store Phi
    fClustersLay1[iC1][2] = z/r;               // Store scaled z 
  }
  
  // Loop on layer 2 : finding theta, phi and r   
  for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {    
    Float_t x = fClustersLay2[iC2][0] - vtx[0];
    Float_t y = fClustersLay2[iC2][1] - vtx[1];
    Float_t z = fClustersLay2[iC2][2] - vtx[2];
    
    Float_t r    = TMath::Sqrt(TMath::Power(x,2) +
			       TMath::Power(y,2) +
			       TMath::Power(z,2));
    
    fClustersLay2[iC2][0] = TMath::ACos(z/r);  // Store Theta
    fClustersLay2[iC2][1] = TMath::ATan(y/x);  // Store Phi
    fClustersLay2[iC2][2] = z;                 // Store z

    // this only needs to be initialized for the fNClustersLay2 first associations
    fAssociationFlag[iC2] = kFALSE;
  }  
  
  //###########################################################
  // Loop on layer 1 
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    

    // reset of variables for multiple candidates
    Int_t   iC2WithBestPhi = 0;     // reset 
    Float_t dPhimin        = 100.;  // just to put a huge number! 
    
    // Loop on layer 2 
    for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {      
      
      // The following excludes double associations
      if (!fAssociationFlag[iC2]) {
	
	// find the difference in angles
	Float_t dTheta = fClustersLay2[iC2][0] - fClustersLay1[iC1][0];
	Float_t dPhi   = fClustersLay2[iC2][1] - fClustersLay1[iC1][1];
	
	// find the difference in z (between linear projection from layer 1
	// and the actual point: Dzeta= z1/r1*r2 -z2) 	
	Float_t r2     = fClustersLay2[iC2][2]/TMath::Cos(fClustersLay2[iC2][0]);
	Float_t dZeta  = fClustersLay2[iC1][2]*r2 - fClustersLay2[iC2][2]; 
	
	if (fHistOn) {
	  fhClustersDPhi->Fill(dPhi);    
	  fhClustersDTheta->Fill(dTheta);    
	  fhClustersDZeta->Fill(dZeta);    
	  fhDPhiVsDThetaAll->Fill(dTheta, dPhi);
	}
	// make "elliptical" cut in Phi and Zeta! 
	Float_t d = TMath::Sqrt(TMath::Power(dPhi/fPhiWindow,2) + TMath::Power(dZeta/fZetaWindow,2));
	if (d>1) continue;      
	
	//look for the minimum distance in Phi: the minimum is in iC2WithBestPhi
	if (TMath::Abs(dPhi) < dPhimin) {
	  dPhimin = TMath::Abs(dPhi);
	  iC2WithBestPhi = iC2;
	}
      } 
    } // end of loop over clusters in layer 2 
    
    if (dPhimin<100) { // This means that a cluster in layer 2 was found that mathes with iC1
      
      if (fOnlyOneTrackletPerC2) fAssociationFlag[iC2WithBestPhi] = kTRUE; // flag the association
      
      // store the tracklet
      
      // use the average theta from the clusters in the two layers
      fTracklets[fNTracklets][0] = 0.5*(fClustersLay1[iC1][0]+fClustersLay2[iC2WithBestPhi][0]);
      // use the phi from the clusters in the first layer 
      fTracklets[fNTracklets][1] = fClustersLay1[iC1][1];
      // Store the difference between phi1 and phi2
      fTracklets[fNTracklets][2] = fClustersLay1[iC1][1] - fClustersLay2[iC2WithBestPhi][1];         
      fNTracklets++;
      
      AliDebug(1,Form(" Adding tracklet candidate %d (cluster %d  of layer 1 and %d  of layer 2)", fNTracklets, iC1));
    }
  } // end of loop over clusters in layer 1
  
  AliDebug(1,Form("%d tracklets found", fNTracklets));
}

//____________________________________________________________________
void
AliITSMultReconstructor::LoadClusterArrays(TTree* itsClusterTree) {
  // This method
  // - gets the clusters from the cluster tree 
  // - convert them into global coordinates 
  // - store them in the internal arrays
  
  AliDebug(1,"Loading clusters ...");
  
  fNClustersLay1 = 0;
  fNClustersLay2 = 0;
  
  TClonesArray* itsClusters = new TClonesArray("AliITSclusterV2");
  TBranch* itsClusterBranch=itsClusterTree->GetBranch("Clusters");
  itsClusterBranch->SetAddress(&itsClusters);
  
  Int_t nItsSubs = (Int_t)itsClusterTree->GetEntries();  
  
  // loop over the its subdetectors
  for (Int_t iIts=0; iIts < nItsSubs; iIts++) {
    
    if (!itsClusterTree->GetEvent(iIts)) 
      continue;
    
    Int_t nClusters = itsClusters->GetEntriesFast();
    
    // stuff needed to get the global coordinates
    Double_t rot[9];   fGeometry->GetRotMatrix(iIts,rot);
    Int_t lay,lad,det; fGeometry->GetModuleId(iIts,lay,lad,det);
    Float_t tx,ty,tz;  fGeometry->GetTrans(lay,lad,det,tx,ty,tz);
    
    // Below:
    // "alpha" is the angle from the global X-axis to the
    //         local GEANT X'-axis  ( rot[0]=cos(alpha) and rot[1]=sin(alpha) )
    // "phi" is the angle from the global X-axis to the
    //       local cluster X"-axis
    
    Double_t alpha   = TMath::ATan2(rot[1],rot[0])+TMath::Pi();
    Double_t itsPhi = TMath::Pi()/2+alpha;
    
    if (lay==1) itsPhi+=TMath::Pi();
    Double_t cp=TMath::Cos(itsPhi), sp=TMath::Sin(itsPhi);
    Double_t r=tx*cp+ty*sp;
    
    // loop over clusters
    while(nClusters--) {
      AliITSclusterV2* cluster = (AliITSclusterV2*)itsClusters->UncheckedAt(nClusters);	
      
      if (cluster->GetLayer()>1) 
	continue;            
      
      Float_t x = r*cp - cluster->GetY()*sp;
      Float_t y = r*sp + cluster->GetY()*cp;
      Float_t z = cluster->GetZ();      
      
      if (cluster->GetLayer()==0) {
	fClustersLay1[fNClustersLay1][0] = x;
	fClustersLay1[fNClustersLay1][1] = y;
	fClustersLay1[fNClustersLay1][2] = z;
	fNClustersLay1++;
      }
      if (cluster->GetLayer()==1) {	
	fClustersLay2[fNClustersLay2][0] = x;
	fClustersLay2[fNClustersLay2][1] = y;
	fClustersLay2[fNClustersLay2][2] = z;
	fNClustersLay2++;
      }
      
    }// end of cluster loop
  } // end of its "subdetector" loop  
  
  AliDebug(1,Form("(clusters in layer 1 : %d,  layer 2: %d)",fNClustersLay1,fNClustersLay2));
}
//____________________________________________________________________
void
AliITSMultReconstructor::SaveHists() {
  
  if (!fHistOn)
    return;

  fhClustersDPhi->Write();
  fhClustersDTheta->Write();
  fhClustersDZeta->Write();
  fhDPhiVsDThetaAll->Write();
  fhDPhiVsDThetaAcc->Write();
}
