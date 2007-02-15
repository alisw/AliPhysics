#include "TFile.h"
#include "TTree.h"

#include "AliAODEvent.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"

#include <iostream>

void CreateAODfromESD(const char *inFileName = "AliESD.root", const char *outFileName = "AliAOD.root") {

  // create an AliAOD object 
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();

  // open the file
  TFile *outFile = TFile::Open(outFileName, "RECREATE");

  // create the tree
  TTree *aodTree = new TTree("AOD", "AliAOD tree");
  aodTree->Branch(aod->GetList());

  // connect to ESD
  TFile *inFile = TFile::Open(inFileName, "READ");
  TTree *t = (TTree*) inFile->Get("esdTree");
  TBranch *b = t->GetBranch("ESD");
  AliESD *esd = 0;
  b->SetAddress(&esd);

  Int_t nEvents = b->GetEntries();

  // loop over events and fill them
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    b->GetEntry(iEvent);

    // Multiplicity information needed by the header (to be revised!)
    Int_t nTracks   = esd->GetNumberOfTracks();
    Int_t nPosTracks = 0;
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) 
      if (esd->GetTrack(iTrack)->GetSign()> 0) nPosTracks++;

    // create the header
    aod->AddHeader(new AliAODHeader(esd ->GetEventNumber(), 
				    esd->GetRunNumber(),
				    nTracks,
				    nPosTracks,
				    nTracks-nPosTracks,
				    esd->GetMagneticField(),
				    -999., // centrality; to be filled, still
				    esd->GetTriggerMask(),
				    esd->GetTriggerCluster(),
				    esd->GetEventType()));

    Int_t nV0s      = esd->GetNumberOfV0s();
    Int_t nCascades = esd->GetNumberOfCascades();
    Int_t nKinks    = esd->GetNumberOfKinks();
    Int_t nVertices = nV0s + nCascades + nKinks;
    
    aod->ResetStd(nTracks, nVertices);
    

    // Array to take into account the tracks already added to the AOD
    Bool_t * usedTrack = NULL;
    if (nTracks>0) {
      usedTrack = new Bool_t[nTracks];
      for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) usedTrack[iTrack]=kFALSE;
    }
    // Array to take into account the V0s already added to the AOD
    Bool_t * usedV0 = NULL;
    if (nV0s>0) {
      usedV0 = new Bool_t[nV0s];
      for (Int_t iV0=0; iV0<nV0s; ++iV0) usedV0[iV0]=kFALSE;
    }

    // Access to the AOD container of vertices
    TClonesArray &vertices = *(aod->GetVertices());
    Int_t jVertices=0;

    // Access to the AOD container of tracks
    TClonesArray &tracks = *(aod->GetTracks());
    Int_t jTracks=0; 
  
    // Add primary vertex. The primary tracks will be defined
    // after the loops on the composite objects (V0, cascades, kinks)
    const AliESDVertex *vtx = esd->GetPrimaryVertex();
      
    Double_t pos[3];
    vtx->GetXYZ(pos); // position
    Double_t cov[6];
    vtx->GetCovMatrix(cov); //covariance matrix

    AliAODVertex * primary = new(vertices[jVertices++])
      AliAODVertex(pos, cov, vtx->GetChi2(), NULL, AliAODVertex::kPrimary);
         
    // Create vertices starting from the most complex objects
      
    // Cascades
    for (Int_t nCascade = 0; nCascade < nCascades; ++nCascade) {
      AliESDcascade *cascade = esd->GetCascade(nCascade);
      
      Double_t posXi[3];
      cascade->GetXYZ(posXi[0], posXi[1], posXi[2]);
      Double_t covXi[6];
      cascade->GetPosCovXi(covXi);
     
      // Add the cascade vertex
      AliAODVertex * vcascade = new(vertices[jVertices++]) AliAODVertex(posXi,
									covXi,
									cascade->GetChi2Xi(),
									primary,
									AliAODVertex::kCascade);

      primary->AddDaughter(vcascade);

      // Add the V0 from the cascade. The ESD class have to be optimized...
      // Now we have to search for the corresponding Vo in the list of V0s
      // using the indeces of the positive and negative tracks

      Int_t posFromV0 = cascade->GetPindex();
      Int_t negFromV0 = cascade->GetNindex();


      AliESDv0 * v0 = 0x0;
      Int_t indV0 = -1;

      for (Int_t iV0=0; iV0<nV0s; ++iV0) {

	v0 = esd->GetV0(iV0);
	Int_t pos = v0->GetPindex();
	Int_t neg = v0->GetNindex();

	if (pos==posFromV0 && neg==negFromV0) {
	  indV0 = iV0;
	  break;
	}
      }

      AliAODVertex * vV0FromCascade = 0x0;

      if (indV0>-1) {
	
	usedV0[indV0] = kTRUE;
	
	Double_t posV0[3];
	v0->GetXYZ(posV0[0], posV0[1], posV0[2]);
	Double_t covV0[6];
	v0->GetPosCov(covV0);
      
	vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(posV0,
								 covV0,
								 v0->GetChi2V0(),
								 vcascade,
								 AliAODVertex::kV0);
      } else {

	Double_t posV0[3];
	cascade->GetXYZ(posV0[0], posV0[1], posV0[2]);
	Double_t covV0[6];
	cascade->GetPosCov(covV0);
      
	vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(posV0,
								 covV0,
								 v0->GetChi2V0(),
								 vcascade,
								 AliAODVertex::kV0);
	vcascade->AddDaughter(vV0FromCascade);
      }

      // Add the positive tracks from the V0

      usedTrack[posFromV0] = kTRUE;
      {

	AliESDtrack *esdTrack = esd->GetTrack(posFromV0);
      
	Double_t p[3];
	esdTrack->GetPxPyPz(p);
	
	Double_t x[3];
	esdTrack->GetXYZ(x);
	
	Double_t cov[21];
	esdTrack->GetCovarianceXYZPxPyPz(cov);

	Double_t pid[10];
	esdTrack->GetESDpid(pid);
	
	vV0FromCascade->AddDaughter(
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(), 
					   p, kTRUE,
					   x,
					   kFALSE,
					   cov, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0FromCascade,
					   AliAODTrack::kSecondary)
		);
      }

      // Add the negative tracks from the V0

      usedTrack[negFromV0] = kTRUE;
      {

	AliESDtrack *esdTrack = esd->GetTrack(negFromV0);
      
	Double_t p[3];
	esdTrack->GetPxPyPz(p);
	
	Double_t x[3];
	esdTrack->GetXYZ(x);
	
	Double_t cov[21];
	esdTrack->GetCovarianceXYZPxPyPz(cov);
	
	Double_t pid[10];
	esdTrack->GetESDpid(pid);

	vV0FromCascade->AddDaughter(
                new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   x,
					   kFALSE,
					   cov, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0FromCascade,
					   AliAODTrack::kSecondary)
		);
      }

      // Add the bachelor track from the cascade

      Int_t bachelor = cascade->GetBindex();

      usedTrack[bachelor] = kTRUE;
      {

	AliESDtrack *esdTrack = esd->GetTrack(bachelor);
      
	Double_t p[3];
	esdTrack->GetPxPyPz(p);
	
	Double_t x[3];
	esdTrack->GetXYZ(x);
	
	Double_t cov[21];
	esdTrack->GetCovarianceXYZPxPyPz(cov);
	
	Double_t pid[10];
	esdTrack->GetESDpid(pid);

	vcascade->AddDaughter(
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   x,
					   kFALSE,
					   cov, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vcascade,
					   AliAODTrack::kSecondary)
		);
      }

      // Add the primary track of the cascade (if any)

    }
    
    // V0s
        
    for (Int_t nV0 = 0; nV0 < nV0s; ++nV0) {

      if (usedV0[nV0]) continue; // skip if aready added to the AOD

      AliESDv0 *v0 = esd->GetV0(nV0);
      
      Double_t posV0[3];
      v0->GetXYZ(posV0[0], posV0[1], posV0[2]);
      Double_t covV0[6];
      v0->GetPosCov(covV0);

      AliAODVertex * vV0 = 
	new(vertices[jVertices++]) AliAODVertex(posV0,
						covV0,
						v0->GetChi2V0(),
						primary,
						AliAODVertex::kV0);
      primary->AddDaughter(vV0);

      Int_t posFromV0 = v0->GetPindex();
      Int_t negFromV0 = v0->GetNindex();

      // Add the positive tracks from the V0

      usedTrack[posFromV0] = kTRUE;
      {

	AliESDtrack *esdTrack = esd->GetTrack(posFromV0);
      
	Double_t p[3];
	esdTrack->GetPxPyPz(p);
	
	Double_t x[3];
	esdTrack->GetXYZ(x);
	
	Double_t cov[21];
	esdTrack->GetCovarianceXYZPxPyPz(cov);

	Double_t pid[10];
	esdTrack->GetESDpid(pid);
	
	vV0->AddDaughter(
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(), 
					   p, kTRUE,
					   x,
					   kFALSE,
					   cov, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0,
					   AliAODTrack::kSecondary)
		);
      }

      // Add the negative tracks from the V0

      usedTrack[negFromV0] = kTRUE;
      {

	AliESDtrack *esdTrack = esd->GetTrack(negFromV0);
      
	Double_t p[3];
	esdTrack->GetPxPyPz(p);
	
	Double_t x[3];
	esdTrack->GetXYZ(x);
	
	Double_t cov[21];
	esdTrack->GetCovarianceXYZPxPyPz(cov);
	
	Double_t pid[10];
	esdTrack->GetESDpid(pid);

	vV0->AddDaughter(
                new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   x,
					   kFALSE,
					   cov, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0,
					   AliAODTrack::kSecondary)
		);
      }


    }
    
    // Kinks
    for (Int_t nKink = 0; nKink < nKinks; ++nKink) {

      AliESDkink *kink = esd->GetKink(nKink);
      
      new(vertices[jVertices++]) AliAODVertex(kink->GetPosition(),
					      NULL,
					      0.,
					      NULL,
					      AliAODVertex::kKink); // create a vtx      
    }
    
    // create tracks
      
    for (Int_t nTrack = 0; nTrack < nTracks; ++nTrack) {
	

      if (usedTrack[nTrack]) continue;

      AliESDtrack *esdTrack = esd->GetTrack(nTrack);
      
      Double_t p[3];
      esdTrack->GetPxPyPz(p);
      
      Double_t x[3];
      esdTrack->GetXYZ(x);
      
      Double_t cov[21];
      esdTrack->GetCovarianceXYZPxPyPz(cov);
	
      Double_t pid[10];
      esdTrack->GetESDpid(pid);
      
      primary->AddDaughter(
      new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					 esdTrack->GetLabel(),
					 p,
					 kTRUE,
					 x,
					 kFALSE,
					 cov, 
					 (Short_t)esdTrack->GetSign(),
					 esdTrack->GetITSClusterMap(), 
					 pid,
					 primary,
					 AliAODTrack::kPrimary)
      );
      
    }
      
    // fill the tree for this event
    aodTree->Fill();
  }
  
  aodTree->GetUserInfo()->Add(aod);

  // close ESD file
  inFile->Close();
  
  // write the tree to the specified file
  outFile = aodTree->GetCurrentFile();
  outFile->cd();
  aodTree->Write();
  outFile->Close();

}
