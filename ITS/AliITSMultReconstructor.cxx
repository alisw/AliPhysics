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
// NOTE: The cuts on phi and zeta depends on the interacting system (p-p  
//  or Pb-Pb). Please, check the file AliITSMultReconstructor.h and be 
//  sure that SetPhiWindow and SetZetaWindow are defined accordingly.
// 
//  
//  
//
//____________________________________________________________________

#include "AliITSMultReconstructor.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliITSRecPoint.h"
#include "AliITSgeom.h"
#include "AliLog.h"

//____________________________________________________________________
ClassImp(AliITSMultReconstructor)


//____________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor() {
  // Method to reconstruct the charged particles multiplicity with the 
  // SPD (tracklets).

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
  fhClustersDPhiAcc   = new TH1F("dphiacc",  "dphi",  100,-0.1,0.1);
  fhClustersDPhiAcc->SetDirectory(0);
  fhClustersDThetaAcc = new TH1F("dthetaacc","dtheta",100,-0.1,0.1);
  fhClustersDThetaAcc->SetDirectory(0);
  fhClustersDZetaAcc = new TH1F("dzetaacc","dzeta",100,-1.,1.);
  fhClustersDZetaAcc->SetDirectory(0);

  fhDPhiVsDZetaAcc = new TH2F("dphiVsDzetaacc","",100,-1.,1.,100,-0.1,0.1);
  fhDPhiVsDZetaAcc->SetDirectory(0);
  fhDPhiVsDThetaAcc = new TH2F("dphiVsDthetaAcc","",100,-0.1,0.1,100,-0.1,0.1);
  fhDPhiVsDThetaAcc->SetDirectory(0);

  fhClustersDPhiAll   = new TH1F("dphiall",  "dphi",  100,-0.5,0.5);
  fhClustersDPhiAll->SetDirectory(0);
  fhClustersDThetaAll = new TH1F("dthetaall","dtheta",100,-0.5,0.5);
  fhClustersDThetaAll->SetDirectory(0);
  fhClustersDZetaAll = new TH1F("dzetaall","dzeta",100,-5.,5.);
  fhClustersDZetaAll->SetDirectory(0);

  fhDPhiVsDZetaAll = new TH2F("dphiVsDzetaall","",100,-5.,5.,100,-0.5,0.5);
  fhDPhiVsDZetaAll->SetDirectory(0);
  fhDPhiVsDThetaAll = new TH2F("dphiVsDthetaAll","",100,-0.5,0.5,100,-0.5,0.5);
  fhDPhiVsDThetaAll->SetDirectory(0);

  fhetaTracklets  = new TH1F("etaTracklets",  "eta",  100,-2.,2.);
  fhphiTracklets  = new TH1F("phiTracklets",  "phi",  100,-3.14159,3.14159);
  fhetaClustersLay1  = new TH1F("etaClustersLay1",  "etaCl1",  100,-2.,2.);
  fhphiClustersLay1  = new TH1F("phiClustersLay1", "phiCl1", 100,-3.141,3.141);

}

//______________________________________________________________________
AliITSMultReconstructor::AliITSMultReconstructor(const AliITSMultReconstructor &mr) : TObject(mr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSMultReconstructor","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSMultReconstructor& AliITSMultReconstructor::operator=(const AliITSMultReconstructor& /* mr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}

//______________________________________________________________________
AliITSMultReconstructor::~AliITSMultReconstructor(){
  // Destructor

  // delete histograms
  delete fhClustersDPhiAcc;
  delete fhClustersDThetaAcc;
  delete fhClustersDZetaAcc;
  delete fhClustersDPhiAll;
  delete fhClustersDThetaAll;
  delete fhClustersDZetaAll;
  delete fhDPhiVsDThetaAll;
  delete fhDPhiVsDThetaAcc;
  delete fhDPhiVsDZetaAll;
  delete fhDPhiVsDZetaAcc;
  delete fhetaTracklets;
  delete fhphiTracklets;
  delete fhetaClustersLay1;
  delete fhphiClustersLay1;

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
    fClustersLay1[iC1][1] = TMath::ATan2(x,y);  // Store Phi
    fClustersLay1[iC1][2] = z/r;               // Store scaled z 
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
   
    Float_t r    = TMath::Sqrt(TMath::Power(x,2) +
			       TMath::Power(y,2) +
			       TMath::Power(z,2));
    
    fClustersLay2[iC2][0] = TMath::ACos(z/r);  // Store Theta
    fClustersLay2[iC2][1] = TMath::ATan2(x,y);  // Store Phi
    fClustersLay2[iC2][2] = z;                 // Store z

 // this only needs to be initialized for the fNClustersLay2 first associations
    fAssociationFlag[iC2] = kFALSE;
  }  
  
  //###########################################################
  // Loop on layer 1 
  for (Int_t iC1=0; iC1<fNClustersLay1; iC1++) {    

    // reset of variables for multiple candidates
    Int_t  iC2WithBestDist = 0;     // reset 
    Float_t distmin        = 100.;  // just to put a huge number! 
    Float_t dPhimin        = 0.;  // Used for histograms only! 
    Float_t dThetamin      = 0.;  // Used for histograms only! 
    Float_t dZetamin       = 0.;  // Used for histograms only! 
    
    // Loop on layer 2 
    for (Int_t iC2=0; iC2<fNClustersLay2; iC2++) {      
      
      // The following excludes double associations
      if (!fAssociationFlag[iC2]) {
	
	// find the difference in angles
	Float_t dTheta = fClustersLay2[iC2][0] - fClustersLay1[iC1][0];
	Float_t dPhi   = fClustersLay2[iC2][1] - fClustersLay1[iC1][1];
	
	// find the difference in z (between linear projection from layer 1
	// and the actual point: Dzeta= z1/r1*r2 -z2) 	
	Float_t r2   = fClustersLay2[iC2][2]/TMath::Cos(fClustersLay2[iC2][0]);
        Float_t dZeta  = fClustersLay1[iC1][2]*r2 - fClustersLay2[iC2][2]; 

 	if (fHistOn) {
	  fhClustersDPhiAll->Fill(dPhi);    
	  fhClustersDThetaAll->Fill(dTheta);    
	  fhClustersDZetaAll->Fill(dZeta);    
	  fhDPhiVsDThetaAll->Fill(dTheta, dPhi);
	  fhDPhiVsDZetaAll->Fill(dZeta, dPhi);
	}
	// make "elliptical" cut in Phi and Zeta! 
	Float_t d = TMath::Sqrt(TMath::Power(dPhi/fPhiWindow,2) + TMath::Power(dZeta/fZetaWindow,2));

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
    
    if (distmin<100) { // This means that a cluster in layer 2 was found that mathes with iC1

      if (fHistOn) {
	fhClustersDPhiAcc->Fill(dPhimin);    
	fhClustersDThetaAcc->Fill(dThetamin);    
	fhClustersDZetaAcc->Fill(dZetamin);    
	fhDPhiVsDThetaAcc->Fill(dThetamin, dPhimin);
	fhDPhiVsDZetaAcc->Fill(dZetamin, dPhimin);
      }
      
      if (fOnlyOneTrackletPerC2) fAssociationFlag[iC2WithBestDist] = kTRUE; // flag the association
      
      // store the tracklet
      
      // use the theta from the clusters in the first layer 
      fTracklets[fNTracklets][0] = fClustersLay1[iC1][0];
      // use the phi from the clusters in the first layer 
      fTracklets[fNTracklets][1] = fClustersLay1[iC1][1];
      // Store the difference between phi1 and phi2
      fTracklets[fNTracklets][2] = fClustersLay1[iC1][1] - fClustersLay2[iC2WithBestDist][1];       
  
      if (fHistOn) {
	Float_t eta=fTracklets[fNTracklets][0];
	eta= TMath::Tan(eta/2.);
	eta=-TMath::Log(eta);
	fhetaTracklets->Fill(eta);    
	fhphiTracklets->Fill(fTracklets[fNTracklets][1]);    
      }
      
      AliDebug(1,Form(" Adding tracklet candidate %d ", fNTracklets));
      AliDebug(1,Form(" Cl. %d of Layer 1 and %d of Layer 2", iC1, 
		      iC2WithBestDist));
      fNTracklets++;
    }

    // Delete the following else if you do not want to save Clusters! 

    else { // This means that the cluster has not been associated 

      // store the cluster
       
      fTracklets[fNTracklets][0] = fClustersLay1[iC1][0];
      fTracklets[fNTracklets][1] = fClustersLay1[iC1][1];
      // Store a flag. This will indicate that the "tracklet" 
      // was indeed a single cluster! 
      fTracklets[fNTracklets][2] = -999999.;       
      AliDebug(1,Form(" Adding a single cluster %d (cluster %d  of layer 1)", 
		      fNTracklets, iC1));
      fNTracklets++;
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
  
  TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
  TBranch* itsClusterBranch=itsClusterTree->GetBranch("ITSRecPoints");

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
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);	
      
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
  // This method save the histograms on the output file
  // (only if fHistOn is TRUE). 
  
  if (!fHistOn)
    return;

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
}
