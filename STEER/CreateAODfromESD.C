#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TArrayS.h>
#include <TArrayD.h>

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODCaloCluster.h"
#include "AliAODPmdCluster.h"
#include "AliAODTracklets.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDkink.h"
#include "AliESDcascade.h"
#include "AliESDCaloCluster.h"
#include "AliESDPmdTrack.h"
#include "AliMultiplicity.h"

#endif

void CreateAODfromESD(const char *inFileName = "AliESDs.root",
		      const char *outFileName = "AliAOD.root") {

  // open input file
  TFile *inFile = TFile::Open(inFileName, "READ");

  // create an AliAOD object 
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();

  // open output file
  TFile *outFile = TFile::Open(outFileName, "RECREATE");
  outFile->cd();

  // create the tree
  TTree *aodTree = new TTree("aodTree", "AliAOD tree");
  aodTree->Branch(aod->GetList());

  // connect to ESD
  TTree *t = (TTree*) inFile->Get("esdTree");
  AliESDEvent *esd = new AliESDEvent();
  esd->ReadFromTree(t);

  Int_t nEvents = t->GetEntries();

  // set arrays and pointers
  Float_t posF[3];
  Double_t pos[3];
  Double_t p[3];
  Double_t p_pos[3];
  Double_t p_neg[3];
  Double_t covVtx[6];
  Double_t covTr[21];
  Double_t pid[10];

  // loop over events and fill them
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {
    //cout << "event: " << iEvent << endl;
    t->GetEntry(iEvent);

    // Multiplicity information needed by the header (to be revised!)
    Int_t nTracks   = esd->GetNumberOfTracks();
    Int_t nPosTracks = 0;
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) 
      if (esd->GetTrack(iTrack)->Charge()> 0) nPosTracks++;

    // Access the header
    AliAODHeader *header = aod->GetHeader();

    // fill the header
    header->SetRunNumber       (esd->GetRunNumber()       );
    header->SetBunchCrossNumber(esd->GetBunchCrossNumber());
    header->SetOrbitNumber     (esd->GetOrbitNumber()     );
    header->SetPeriodNumber    (esd->GetPeriodNumber()    );
    header->SetTriggerMask     (esd->GetTriggerMask()     ); 
    header->SetTriggerCluster  (esd->GetTriggerCluster()  );
    header->SetEventType       (esd->GetEventType()       );
    header->SetMagneticField   (esd->GetMagneticField()   );
    header->SetZDCN1Energy     (esd->GetZDCN1Energy()     );
    header->SetZDCP1Energy     (esd->GetZDCP1Energy()     );
    header->SetZDCN2Energy     (esd->GetZDCN2Energy()     );
    header->SetZDCP2Energy     (esd->GetZDCP2Energy()     );
    header->SetZDCEMEnergy     (esd->GetZDCEMEnergy()     );
    header->SetRefMultiplicity   (nTracks);
    header->SetRefMultiplicityPos(nPosTracks);
    header->SetRefMultiplicityNeg(nTracks - nPosTracks);
    header->SetMuonMagFieldScale(-999.); // FIXME
    header->SetCentrality(-999.);        // FIXME

    Int_t nV0s      = esd->GetNumberOfV0s();
    Int_t nCascades = esd->GetNumberOfCascades();
    Int_t nKinks    = esd->GetNumberOfKinks();
    Int_t nVertices = nV0s + nCascades + nKinks + 1 /* = prim. vtx*/;
    Int_t nJets     = 0;
    Int_t nCaloClus = esd->GetNumberOfCaloClusters();
    Int_t nFmdClus  = 0;
    Int_t nPmdClus  = esd->GetNumberOfPmdTracks();
   
    aod->ResetStd(nTracks, nVertices, nV0s+nCascades, nJets, nCaloClus, nFmdClus, nPmdClus);
    
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
   
    // Access to the AOD container of V0s
    TClonesArray &V0s = *(aod->GetV0s());
    Int_t jV0s=0;
    
    // Add primary vertex. The primary tracks will be defined
    // after the loops on the composite objects (V0, cascades, kinks)
    const AliESDVertex *vtx = esd->GetPrimaryVertex();
      
    vtx->GetXYZ(pos); // position
    vtx->GetCovMatrix(covVtx); //covariance matrix

    AliAODVertex * primary = new(vertices[jVertices++])
      AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);
         

    AliAODTrack *aodTrack = 0x0;
    
    // Create vertices starting from the most complex objects

    // Cascades
    for (Int_t nCascade = 0; nCascade < nCascades; ++nCascade) {
      AliESDcascade *cascade = esd->GetCascade(nCascade);
      
      cascade->GetXYZ(pos[0], pos[1], pos[2]);
      cascade->GetPosCovXi(covVtx);
     
      // Add the cascade vertex
      AliAODVertex * vcascade = new(vertices[jVertices++]) AliAODVertex(pos,
									covVtx,
									cascade->GetChi2Xi(), // = chi2/NDF since NDF = 2*2-3
									primary,
									nCascade,
									AliAODVertex::kCascade);

      primary->AddDaughter(vcascade); // the cascade 'particle' (represented by a vertex) is added as a daughter to the primary vertex

      // Add the V0 from the cascade. The ESD class have to be optimized...
      // Now we have to search for the corresponding V0 in the list of V0s
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

      if (indV0>-1 && !usedV0[indV0]) {
	
	// the V0 exists in the array of V0s and is not used

	usedV0[indV0] = kTRUE;
	
	v0->GetXYZ(pos[0], pos[1], pos[2]);
	v0->GetPosCov(covVtx);
	
	vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
								 covVtx,
								 v0->GetChi2V0(), // = chi2/NDF since NDF = 2*2-3
								 vcascade,
								 indV0,
								 AliAODVertex::kV0);
      } else {

	// the V0 doesn't exist in the array of V0s or was used
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " The V0 " << indV0 
	     << " doesn't exist in the array of V0s or was used!" << endl;

	cascade->GetXYZ(pos[0], pos[1], pos[2]);
	cascade->GetPosCov(covVtx);
      
	vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
								 covVtx,
								 v0->GetChi2V0(), // = chi2/NDF since NDF = 2*2-3
								 vcascade,
								 indV0,
								 AliAODVertex::kV0);
	vcascade->AddDaughter(vV0FromCascade);

      }

      // Add the positive tracks from the V0

      if (! usedTrack[posFromV0]) {

	usedTrack[posFromV0] = kTRUE;

	AliESDtrack *esdTrack = esd->GetTrack(posFromV0);
	esdTrack->GetPxPyPz(p_pos);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);
	
	vV0FromCascade->AddDaughter(aodTrack =
				    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(), 
					   p_pos, 
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0FromCascade,
					   kTRUE,  // check if this is right
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
	aodTrack->ConvertAliPIDtoAODPID();
      }
      else {
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " track " << posFromV0 << " has already been used!" << endl;
      }

      // Add the negative tracks from the V0

      if (!usedTrack[negFromV0]) {
	
	usedTrack[negFromV0] = kTRUE;
	
	AliESDtrack *esdTrack = esd->GetTrack(negFromV0);
	esdTrack->GetPxPyPz(p_neg);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);
	
	vV0FromCascade->AddDaughter(aodTrack =
                new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p_neg,
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0FromCascade,
					   kTRUE,  // check if this is right
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
	aodTrack->ConvertAliPIDtoAODPID();
      }
      else {
	cerr << "Error: event " << iEvent << " cascade " << nCascade
	     << " track " << negFromV0 << " has already been used!" << endl;
      }

      // add it to the V0 array as well
      Double_t d0[2] = { -999., -99.};
      // counting is probably wrong
      new(V0s[jV0s++]) AliAODv0(vV0FromCascade, -999., -99., p_pos, p_neg, d0); // to be refined

      // Add the bachelor track from the cascade

      Int_t bachelor = cascade->GetBindex();
      
      if(!usedTrack[bachelor]) {
      
	usedTrack[bachelor] = kTRUE;
	
	AliESDtrack *esdTrack = esd->GetTrack(bachelor);
	esdTrack->GetPxPyPz(p);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);

	vcascade->AddDaughter(aodTrack =
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vcascade,
					   kTRUE,  // check if this is right
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
 	aodTrack->ConvertAliPIDtoAODPID();
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
     
      v0->GetXYZ(pos[0], pos[1], pos[2]);
      v0->GetPosCov(covVtx);

      AliAODVertex * vV0 = 
	new(vertices[jVertices++]) AliAODVertex(pos,
						covVtx,
						v0->GetChi2V0(), // = chi2/NDF since NDF = 2*2-3
						primary,
						nV0,
						AliAODVertex::kV0);
      primary->AddDaughter(vV0);

      Int_t posFromV0 = v0->GetPindex();
      Int_t negFromV0 = v0->GetNindex();
      
      // Add the positive tracks from the V0

      if (!usedTrack[posFromV0]) {
	
	usedTrack[posFromV0] = kTRUE;

	AliESDtrack *esdTrack = esd->GetTrack(posFromV0);
	esdTrack->GetPxPyPz(p_pos);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);
	
	vV0->AddDaughter(aodTrack =
        	new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(), 
					   p_pos, 
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0,
					   kTRUE,  // check if this is right
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
	aodTrack->ConvertAliPIDtoAODPID();
      }
      else {
	cerr << "Error: event " << iEvent << " V0 " << nV0
	     << " track " << posFromV0 << " has already been used!" << endl;
      }

      // Add the negative tracks from the V0

      if (!usedTrack[negFromV0]) {

	usedTrack[negFromV0] = kTRUE;

	AliESDtrack *esdTrack = esd->GetTrack(negFromV0);
	esdTrack->GetPxPyPz(p_neg);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);

	vV0->AddDaughter(aodTrack =
                new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p_neg,
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vV0,
					   kTRUE,  // check if this is right
					   kFALSE, // check if this is right
					   AliAODTrack::kSecondary)
		);
	aodTrack->ConvertAliPIDtoAODPID();
      }
      else {
	cerr << "Error: event " << iEvent << " V0 " << nV0
	     << " track " << negFromV0 << " has already been used!" << endl;
      }

      // add it to the V0 array as well
      Double_t d0[2] = { 999., 99.};
      new(V0s[jV0s++]) AliAODv0(vV0, 999., 99., p_pos, p_neg, d0); // to be refined
    } // end of the loop on V0s
    
    // Kinks: it is a big mess the access to the information in the kinks
    // The loop is on the tracks in order to find the mother and daugther of each kink


    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {


      AliESDtrack * esdTrack = esd->GetTrack(iTrack);

      Int_t ikink = esdTrack->GetKinkIndex(0);

      if (ikink) {
	// Negative kink index: mother, positive: daughter

	// Search for the second track of the kink

	for (Int_t jTrack = iTrack+1; jTrack<nTracks; ++jTrack) {

	  AliESDtrack * esdTrack1 = esd->GetTrack(jTrack);

	  Int_t jkink = esdTrack1->GetKinkIndex(0);

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
	      esdTrack->GetPxPyPz(p);
	      esdTrack->GetXYZ(pos);
	      esdTrack->GetCovarianceXYZPxPyPz(covTr);
	      esdTrack->GetESDpid(pid);

	      mother = 
		new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   primary,
					   kTRUE, // check if this is right
					   kTRUE, // check if this is right
					   AliAODTrack::kPrimary);
	      primary->AddDaughter(mother);
	      mother->ConvertAliPIDtoAODPID();
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
						    esdTrack->GetID(), // This is the track ID of the mother's track!
						    AliAODVertex::kKink);
	    // Add the daughter track

	    AliAODTrack * daughter = NULL;

	    if (!usedTrack[idaughter]) {
	
	      usedTrack[idaughter] = kTRUE;
	
	      AliESDtrack *esdTrack = esd->GetTrack(idaughter);
	      esdTrack->GetPxPyPz(p);
	      esdTrack->GetXYZ(pos);
	      esdTrack->GetCovarianceXYZPxPyPz(covTr);
	      esdTrack->GetESDpid(pid);

	      daughter = 
		new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					   esdTrack->GetLabel(),
					   p,
					   kTRUE,
					   pos,
					   kFALSE,
					   covTr, 
					   (Short_t)esdTrack->Charge(),
					   esdTrack->GetITSClusterMap(), 
					   pid,
					   vkink,
					   kTRUE, // check if this is right
					   kTRUE, // check if this is right
					   AliAODTrack::kPrimary);
	      vkink->AddDaughter(daughter);
	      daughter->ConvertAliPIDtoAODPID();
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
      esdTrack->GetPxPyPz(p);
      esdTrack->GetXYZ(pos);
      esdTrack->GetCovarianceXYZPxPyPz(covTr);
      esdTrack->GetESDpid(pid);

      Float_t impactXY, impactZ;

      esdTrack->GetImpactParameters(impactXY,impactZ);

      if (impactXY<3.) {
	// track inside the beam pipe
      
	primary->AddDaughter(aodTrack =
	    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
					 esdTrack->GetLabel(),
					 p,
					 kTRUE,
					 pos,
					 kFALSE,
					 covTr, 
					 (Short_t)esdTrack->Charge(),
					 esdTrack->GetITSClusterMap(), 
					 pid,
					 primary,
					 kTRUE, // check if this is right
					 kTRUE, // check if this is right
					 AliAODTrack::kPrimary)
	    );
	aodTrack->ConvertAliPIDtoAODPID();
      }
      else {
	// outside the beam pipe: orphan track
	// Don't write them anymore!
	continue;
      }	
    } // end of loop on tracks
    
    // muon tracks
    Int_t nMuTracks = esd->GetNumberOfMuonTracks();
    for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
      
      AliESDMuonTrack *esdMuTrack = esd->GetMuonTrack(nMuTrack);     
      p[0] = esdMuTrack->Px(); 
      p[1] = esdMuTrack->Py(); 
      p[2] = esdMuTrack->Pz();
      pos[0] = primary->GetX(); 
      pos[1] = primary->GetY(); 
      pos[2] = primary->GetZ();
      
      // has to be changed once the muon pid is provided by the ESD
      for (Int_t i = 0; i < 10; pid[i++] = 0.); pid[AliAODTrack::kMuon]=1.;
      
      primary->AddDaughter(aodTrack =
	  new(tracks[jTracks++]) AliAODTrack(0, // no ID provided
					     0, // no label provided
					     p,
					     kTRUE,
					     pos,
					     kFALSE,
					     NULL, // no covariance matrix provided
					     esdMuTrack->Charge(),
					     0, // ITSClusterMap is set below
					     pid,
					     primary,
 					     kFALSE,  // muon tracks are not used to fit the primary vtx
					     kFALSE,  // not used for vertex fit
					     AliAODTrack::kPrimary)
	  );

      aodTrack->SetHitsPatternInTrigCh(esdMuTrack->GetHitsPatternInTrigCh());
      Int_t track2Trigger = esdMuTrack->GetMatchTrigger();
      aodTrack->SetMatchTrigger(track2Trigger);
      if (track2Trigger) 
	aodTrack->SetChi2MatchTrigger(esdMuTrack->GetChi2MatchTrigger());
      else 
	aodTrack->SetChi2MatchTrigger(0.);
    }
   
    // Access to the AOD container of PMD clusters
    TClonesArray &pmdClusters = *(aod->GetPmdClusters());
    Int_t jPmdClusters=0;
  
    for (Int_t iPmd = 0; iPmd < nPmdClus; ++iPmd) {
      // file pmd clusters, to be revised!
      AliESDPmdTrack *pmdTrack = esd->GetPmdTrack(iPmd);
      Int_t nLabel = 0;
      Int_t *label = 0x0;
      Double_t pos[3] = { pmdTrack->GetClusterX(), pmdTrack->GetClusterY(), pmdTrack->GetClusterZ() };
      Double_t pid[9] = { 0., 0., 0., 0., 0., 0., 0., 0., 0. }; // to be revised!
      // type not set!
      // assoc cluster not set
      new(pmdClusters[jPmdClusters++]) AliAODPmdCluster(iPmd, nLabel, label, pmdTrack->GetClusterADC(), pos, pid);
    }

    // Access to the AOD container of clusters
    TClonesArray &caloClusters = *(aod->GetCaloClusters());
    Int_t jClusters=0;

    // Calo Clusters
    TArrayS EMCCellNumber(15000);
    TArrayD EMCCellAmplitude(15000);
    Int_t nEMCCells = 0;
    const Float_t fEMCAmpScale = 1./500;
 
    for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {

      AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);

      Int_t id = cluster->GetID();
      Int_t nLabel = 0;
      Int_t *label = 0x0;
      Float_t energy = cluster->E();
      cluster->GetPosition(posF);
      Char_t ttype=AliAODCluster::kUndef;

      if (cluster->GetClusterType() == AliESDCaloCluster::kPHOSCluster) {
	ttype=AliAODCluster::kPHOSNeutral;
      } 
      else if (cluster->GetClusterType() == AliESDCaloCluster::kEMCALClusterv1) {
	ttype = AliAODCluster::kEMCALClusterv1;
      }
      else if (cluster->GetClusterType() == AliESDCaloCluster::kEMCALPseudoCluster) {
	// Collect raw tower info
	for (Int_t iDig = 0; iDig < cluster->GetNumberOfDigits(); iDig++) {
	  EMCCellNumber[nEMCCells] = cluster->GetDigitIndex()->At(iDig);
	  EMCCellAmplitude[nEMCCells] = fEMCAmpScale*cluster->GetDigitAmplitude()->At(iDig);
	  nEMCCells++;
	}
	// don't write cluster data (it's just a pseudo cluster, holding the tower information)
	continue; 
      }
      
      AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) AliAODCaloCluster(id,
											nLabel,
											label,
											energy,
											pos,
											NULL,
											ttype);
      
      caloCluster->SetCaloCluster(); // to be refined!

    } // end of loop on calo clusters

    // fill EMC cell info
    AliAODCaloCells &EMCCells = *(aod->GetCaloCells());
    EMCCells.CreateContainer(nEMCCells);
    EMCCells.SetType(AliAODCaloCells::kEMCAL);
    for (Int_t iCell = 0; iCell < nEMCCells; iCell++) {      
      EMCCells.SetCell(iCell,EMCCellNumber[iCell],EMCCellAmplitude[iCell]);
    }
    EMCCells.Sort();

    // tracklets    
    AliAODTracklets &SPDTracklets = *(aod->GetTracklets());
    const AliMultiplicity *mult = esd->GetMultiplicity();
    if (mult) {
      if (mult->GetNumberOfTracklets()>0) {
	SPDTracklets.CreateContainer(mult->GetNumberOfTracklets());

	for (Int_t n=0; n<mult->GetNumberOfTracklets(); n++) {
	  SPDTracklets.SetTracklet(n, mult->GetTheta(n), mult->GetPhi(n), mult->GetDeltaPhi(n), mult->GetLabel(n));
	}
      }
    } else {
      Printf("ERROR: AliMultiplicity could not be retrieved from ESD");
    }

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
