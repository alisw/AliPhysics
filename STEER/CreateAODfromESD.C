#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODCluster.h"

#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDCaloCluster.h"

#endif

void CreateAODfromESD(const char *inFileName = "AliESDs.root",
		      const char *outFileName = "AliAOD.root") {

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
    aod->AddHeader(new AliAODHeader(esd->GetRunNumber(),
				    esd->GetBunchCrossNumber(),
				    esd->GetOrbitNumber(),
				    nTracks,
				    nPosTracks,
				    nTracks-nPosTracks,
				    esd->GetMagneticField(),
				    -999., // fill muon magnetic field
				    -999., // centrality; to be filled, still
				    esd->GetZDCN1Energy(),
				    esd->GetZDCP1Energy(),
				    esd->GetZDCN2Energy(),
				    esd->GetZDCP2Energy(),
				    esd->GetZDCEMEnergy(),
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
    // Array to take into account the kinks already added to the AOD
    Bool_t * usedKink = NULL;
    if (nKinks>0) {
      usedKink = new Bool_t[nKinks];
      for (Int_t iKink=0; iKink<nKinks; ++iKink) usedKink[iKink]=kFALSE;
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
    Double_t covVtx[6]; // We have to give changing names to the variables (like cov?, x?, and p?) because CINT doesn't recognize blocks correctly.
    vtx->GetCovMatrix(covVtx); //covariance matrix

    AliAODVertex * primary = new(vertices[jVertices++])
      AliAODVertex(pos, covVtx, vtx->GetChi2(), NULL, AliAODVertex::kPrimary);
         
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
	Int_t posV0 = v0->GetPindex();
	Int_t negV0 = v0->GetNindex();

	if (posV0==posFromV0 && negV0==negFromV0) {
	  indV0 = iV0;
	  break;
	}
      }

      AliAODVertex * vV0FromCascade = 0x0;

      if (indV0>-1 && !usedV0[indV0] ) {
	
	// the V0 exists in the array of V0s and is not used

	usedV0[indV0] = kTRUE;
	
	Double_t posV0[3];
	v0->GetXYZ(posV0[0], posV0[1], posV0[2]);
	Double_t covV0_1[6];
	v0->GetPosCov(covV0_1);
	
	vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(posV0,
								 covV0_1,
								 v0->GetChi2V0(),
								 vcascade,
								 AliAODVertex::kV0);
      } else {

	// the V0 doesn't exist in the array of V0s or was used
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " The V0 " << indV0 
	     << " doesn't exist in the array of V0s or was used!" << endl;

	Double_t posV0_2[3];
	cascade->GetXYZ(posV0_2[0], posV0_2[1], posV0_2[2]);
	Double_t covV0_2[6];
	cascade->GetPosCov(covV0_2);
      
	vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(posV0_2,
								 covV0_2,
								 v0->GetChi2V0(),
								 vcascade,
								 AliAODVertex::kV0);
	vcascade->AddDaughter(vV0FromCascade);
      }

      // Add the positive tracks from the V0

      if (! usedTrack[posFromV0]) {

	usedTrack[posFromV0] = kTRUE;

	AliESDtrack *esdTrack = esd->GetTrack(posFromV0);
	
	Double_t p1[3];
	esdTrack->GetPxPyPz(p1);
	
	Double_t x1[3];
	esdTrack->GetXYZ(x1);
	
	Double_t covV0PosTr[21];
	esdTrack->GetCovarianceXYZPxPyPz(covV0PosTr);
	
	Double_t pid1[10];
	esdTrack->GetESDpid(pid1);
	
	vV0FromCascade->AddDaughter(
				    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(), 
					   p1, 
					   kTRUE,
					   x1,
					   kFALSE,
					   covV0PosTr, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid1,
					   vV0FromCascade,
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
      }
      else {
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " track " << posFromV0 << " has already been used!" << endl;
      }

      // Add the negative tracks from the V0

      if (!usedTrack[negFromV0]) {
	
	usedTrack[negFromV0] = kTRUE;
	
	AliESDtrack *esdTrack = esd->GetTrack(negFromV0);
	
	Double_t p2[3];
	esdTrack->GetPxPyPz(p2);
	
	Double_t x2[3];
	esdTrack->GetXYZ(x2);
	
	Double_t covV0NegTr[21];
	esdTrack->GetCovarianceXYZPxPyPz(covV0NegTr);
	
	Double_t pid2[10];
	esdTrack->GetESDpid(pid2);
	
	vV0FromCascade->AddDaughter(
                new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p2,
					   kTRUE,
					   x2,
					   kFALSE,
					   covV0NegTr, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid2,
					   vV0FromCascade,
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
      }
      else {
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " track " << negFromV0 << " has already been used!" << endl;
      }

      // Add the bachelor track from the cascade

      Int_t bachelor = cascade->GetBindex();
      
      if(!usedTrack[bachelor]) {
      
	usedTrack[bachelor] = kTRUE;
	
	AliESDtrack *esdTrack = esd->GetTrack(bachelor);
	
	Double_t p3[3];
	esdTrack->GetPxPyPz(p3);
	
	Double_t x3[3];
	esdTrack->GetXYZ(x3);
	
	Double_t covXiTr[21];
	esdTrack->GetCovarianceXYZPxPyPz(covXiTr);
	
	Double_t pid3[10];
	esdTrack->GetESDpid(pid3);

	vcascade->AddDaughter(
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p3,
					   kTRUE,
					   x3,
					   kFALSE,
					   covXiTr, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid3,
					   vcascade,
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
      }
      else {
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " track " << bachelor << " has already been used!" << endl;
      }

      // Add the primary track of the cascade (if any)

    } // end of the loop on cascades
    
    // V0s
        
    for (Int_t nV0 = 0; nV0 < nV0s; ++nV0) {

      if (usedV0[nV0]) continue; // skip if aready added to the AOD

      AliESDv0 *v0 = esd->GetV0(nV0);
      
      Double_t posV0_3[3];
      v0->GetXYZ(posV0_3[0], posV0_3[1], posV0_3[2]);
      Double_t covV0_3[6];
      v0->GetPosCov(covV0_3);

      AliAODVertex * vV0 = 
	new(vertices[jVertices++]) AliAODVertex(posV0_3,
						covV0_3,
						v0->GetChi2V0(),
						primary,
						AliAODVertex::kV0);
      primary->AddDaughter(vV0);

      Int_t posFromV0 = v0->GetPindex();
      Int_t negFromV0 = v0->GetNindex();
      
      // Add the positive tracks from the V0

      if (!usedTrack[posFromV0]) {
	
	usedTrack[posFromV0] = kTRUE;

	AliESDtrack *esdTrack = esd->GetTrack(posFromV0);
      
	Double_t p4[3];
	esdTrack->GetPxPyPz(p4);
	
	Double_t x4[3];
	esdTrack->GetXYZ(x4);
	
	Double_t covV0PosTr_2[21];
	esdTrack->GetCovarianceXYZPxPyPz(covV0PosTr_2);
	
	Double_t pid4[10];
	esdTrack->GetESDpid(pid4);
	
	vV0->AddDaughter(
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(), 
					   p4, 
					   kTRUE,
					   x4,
					   kFALSE,
					   covV0PosTr_2, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid4,
					   vV0,
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
      }
      else {
	cerr << "Error: event " << iEvent << " V0 " << nV0
	     << " track " << posFromV0 << " has already been used!" << endl;
      }

      // Add the negative tracks from the V0

      if (!usedTrack[negFromV0]) {

	usedTrack[negFromV0] = kTRUE;

	AliESDtrack *esdTrack = esd->GetTrack(negFromV0);
      
	Double_t p5[3];
	esdTrack->GetPxPyPz(p5);
	
	Double_t x5[3];
	esdTrack->GetXYZ(x5);
	
	Double_t covV0NegTr_2[21];
	esdTrack->GetCovarianceXYZPxPyPz(covV0NegTr_2);
	
	Double_t pid5[10];
	esdTrack->GetESDpid(pid5);

	vV0->AddDaughter(
                new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p5,
					   kTRUE,
					   x5,
					   kFALSE,
					   covV0NegTr_2, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid5,
					   vV0,
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
      }
      else {
	cerr << "Error: event " << iEvent << " V0 " << nV0
	     << " track " << negFromV0 << " has already been used!" << endl;
      }

    } // end of the loop on V0s
    
    // Kinks: it is a big mess the access to the information in the kinks
    // The loop is on the tracks in order to find the mother and daugther of each kink


    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {


      AliESDtrack * track = esd->GetTrack(iTrack);


      Int_t ikink = track->GetKinkIndex(0);

      if (ikink) {
	// Negative kink index: mother, positive: daughter

	// Search for the second track of the kink

	for (Int_t jTrack = iTrack+1; jTrack<nTracks; ++jTrack) {

	  AliESDtrack * track1 = esd->GetTrack(jTrack);

	  Int_t jkink = track1->GetKinkIndex(0);

	  if ( TMath::Abs(ikink)==TMath::Abs(jkink) ) {

	    // The two tracks are from the same kink
	  
	    if (usedKink[TMath::Abs(ikink)-1]) continue; // skip used kinks

	    Int_t imother = -1;
	    Int_t idaughter = -1;

	    if (ikink<0 && jkink>0) {

	      imother = iTrack;
	      idaughter = jTrack;
	    }
	    else if (ikink>0 && jkink<0) {

	      imother = jTrack;
	      idaughter = iTrack;
	    }
	    else {
	      cerr << "Error: Wrong combination of kink indexes: "
	      << ikink << " " << jkink << endl;
	      continue;
	    }

	    // Add the mother track

	    AliAODTrack * mother = NULL;

	    if (!usedTrack[imother]) {
	
	      usedTrack[imother] = kTRUE;
	
	      AliESDtrack *esdTrack = esd->GetTrack(imother);
	
	      Double_t p6[3];
	      esdTrack->GetPxPyPz(p6);
	      
	      Double_t x6[3];
	      esdTrack->GetXYZ(x6);
	      
	      Double_t covKinkMother[21];
	      esdTrack->GetCovarianceXYZPxPyPz(covKinkMother);
	      
	      Double_t pid6[10];
	      esdTrack->GetESDpid(pid6);

	      mother = 
		new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p6,
					   kTRUE,
					   x6,
					   kFALSE,
					   covKinkMother, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid6,
					   primary,
					   kTRUE, // check if this is right
					   AliAODTrack::kPrimary);
	      primary->AddDaughter(mother);
	    }
	    else {
	      cerr << "Error: event " << iEvent << " kink " << TMath::Abs(ikink)-1
	      << " track " << imother << " has already been used!" << endl;
	    }

	    // Add the kink vertex
	    AliESDkink * kink = esd->GetKink(TMath::Abs(ikink)-1);

	    AliAODVertex * vkink = 
	    new(vertices[jVertices++]) AliAODVertex(kink->GetPosition(),
						    NULL,
						    0.,
						    mother,
						    AliAODVertex::kKink);
	    // Add the daughter track

	    AliAODTrack * daughter = NULL;

	    if (!usedTrack[idaughter]) {
	
	      usedTrack[idaughter] = kTRUE;
	
	      AliESDtrack *esdTrack = esd->GetTrack(idaughter);
	
	      Double_t p7[3];
	      esdTrack->GetPxPyPz(p7);
	      
	      Double_t x7[3];
	      esdTrack->GetXYZ(x7);
	      
	      Double_t covKinkDaughter[21];
	      esdTrack->GetCovarianceXYZPxPyPz(covKinkDaughter);
	      
	      Double_t pid7[10];
	      esdTrack->GetESDpid(pid7);

	      daughter = 
		new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p7,
					   kTRUE,
					   x7,
					   kFALSE,
					   covKinkDaughter, 
					   (Short_t)esdTrack->GetSign(),
					   esdTrack->GetITSClusterMap(), 
					   pid7,
					   vkink,
					   kTRUE, // check if this is right
					   AliAODTrack::kPrimary);
	      vkink->AddDaughter(daughter);
	    }
	    else {
	      cerr << "Error: event " << iEvent << " kink " << TMath::Abs(ikink)-1
	      << " track " << idaughter << " has already been used!" << endl;
	    }


	  }
	}

      }      

    }

    
    // Tracks (primary and orphan)
      
    for (Int_t nTrack = 0; nTrack < nTracks; ++nTrack) {
	

      if (usedTrack[nTrack]) continue;

      AliESDtrack *esdTrack = esd->GetTrack(nTrack);
      
      Double_t p8[3];
      esdTrack->GetPxPyPz(p8);
      
      Double_t x8[3];
      esdTrack->GetXYZ(x8);
      
      Double_t covTr[21];
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
 	
      Double_t pid8[10];
      esdTrack->GetESDpid(pid8);

      Float_t impactXY, impactZ;

      esdTrack->GetImpactParameters(impactXY,impactZ);

      if (impactXY<3) {
	// track inside the beam pipe
      
	primary->AddDaughter(
	    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					 esdTrack->GetLabel(),
					 p8,
					 kTRUE,
					 x8,
					 kFALSE,
					 covTr, 
					 (Short_t)esdTrack->GetSign(),
					 esdTrack->GetITSClusterMap(), 
					 pid8,
					 primary,
					 kTRUE, // check if this is right
					 AliAODTrack::kPrimary)
	    );
      }
      else {
	// outside the beam pipe: orphan track
	    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					 esdTrack->GetLabel(),
					 p8,
					 kTRUE,
					 x8,
					 kFALSE,
					 covTr, 
					 (Short_t)esdTrack->GetSign(),
					 esdTrack->GetITSClusterMap(), 
					 pid8,
					 NULL,
					 kFALSE, // check if this is right
					 AliAODTrack::kOrphan);
      }	
    } // end of loop on tracks


    // Access to the AOD container of clusters
    TClonesArray &clusters = *(aod->GetClusters());
    Int_t jClusters=0;

    // Calo Clusters
    Int_t nClusters    = esd->GetNumberOfCaloClusters();

    for (Int_t iClust=0; iClust<nClusters; ++iClust) {

      AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);

      Int_t id = cluster->GetID();
      Int_t label = -1;
      Float_t energy = cluster->GetClusterEnergy();
      Float_t x9[3];
      cluster->GetGlobalPosition(x9);
      Float_t * covMatrix = NULL;
      Float_t * pid9 = NULL; 
      AliAODVertex *prodVertex = primary;
      AliAODTrack *primTrack = NULL;
      Char_t ttype=AliAODCluster::kUndef;

      if (cluster->IsPHOS()) ttype=AliAODCluster::kPHOSNeutral;
      else if (cluster->IsEMCAL()) {

	if (cluster->GetClusterType() == AliESDCaloCluster::kPseudoCluster)
	  ttype = AliAODCluster::kEMCALPseudoCluster;
	else
	  ttype = AliAODCluster::kEMCALClusterv1;

      }
      
      new(clusters[jClusters++]) AliAODCluster(id,
					       label,
					       energy,
					       x9,
					       covMatrix,
					       pid9,
					       prodVertex,
					       primTrack,
					       ttype);

    } // end of loop on calo clusters


    delete [] usedTrack;
    delete [] usedV0;
    delete [] usedKink;

    // fill the tree for this event
    aodTree->Fill();
  } // end of event loop
  
  aodTree->GetUserInfo()->Add(aod);

  // close ESD file
  inFile->Close();
  
  // write the tree to the specified file
  outFile = aodTree->GetCurrentFile();
  outFile->cd();
  aodTree->Write();
  outFile->Close();

}
