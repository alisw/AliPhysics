//
// A macro to test the Cooked Matrix tracker out of the reco framework
//

#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <Riostream.h>
   #include <TSystem.h>
   #include <TFile.h>
   #include <TTree.h>
   #include <TGeoGlobalMagField.h>
   #include <TStopwatch.h>

   #include "AliRunLoader.h"
   #include "AliGeomManager.h"
   #include "AliLoader.h"
   #include "AliHeader.h"
   #include "AliGenEventHeader.h"
   #include "AliMagF.h"
   #include "AliESDEvent.h"
   #include "AliITSUTrackerCooked.h"
   #include "AliITSUReconstructor.h"
#endif

extern TSystem *gSystem;

const AliESDVertex *SetMCvertex(const AliRunLoader *rl, AliTracker *tr);

void testTrackerCooked() {
    gSystem->Load("libITSUpgradeBase");
    gSystem->Load("libITSUpgradeSim");
    gSystem->Load("libITSUpgradeRec");

    TGeoGlobalMagField::Instance()->
          SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
    AliGeomManager::LoadGeometry("geometry.root");

    AliITSUReconstructor rec;
    rec.Init();
    AliITSUTrackerCooked *tracker=new AliITSUTrackerCooked(&rec);
    tracker->SetSAonly(kTRUE);

    TFile *esdFile=TFile::Open("AliESDs.root","recreate");
    TTree *esdTree = new TTree("esdTree", "Tree with ESD objects");
    AliESDEvent *esd=new AliESDEvent();
    esd->CreateStdContent();
    esd->WriteToTree(esdTree);
    
    AliRunLoader *rl = AliRunLoader::Open("galice.root","something");
    rl->LoadHeader();
    AliLoader *itsl = rl->GetLoader("ITSLoader");
    itsl->LoadRecPoints();

    TStopwatch timer;    
    Int_t nEvents=rl->GetNumberOfEvents();
    for (Int_t i=0; i<nEvents; i++) {
        cout<<"\nEvent number "<<i<<endl;
        rl->GetEvent(i);

        const AliESDVertex *vtx=SetMCvertex(rl,tracker);
        esd->SetPrimaryVertexSPD(vtx);              

        TTree *cTree=itsl->TreeR();
        tracker->LoadClusters(cTree);
        tracker->Clusters2Tracks(esd);
        tracker->PropagateBack(esd);
        tracker->RefitInward(esd);
        tracker->UnloadClusters();

        Int_t n=esd->GetNumberOfTracks();
        for (Int_t t=0; t<n; t++) {
	  AliESDtrack *track=esd->GetTrack(t);
          if (!track->RelateToVertex(vtx, tracker->GetBz(), 33)) continue;
          //Double_t chi2=track->GetITSchi2();
          //Double_t ncl=track->GetITSclusters(0);
	  //cout<<ncl<<' '<<chi2<<endl;
	}

        esdTree->Fill();
        esd->Reset();
        delete vtx;
    }
    timer.Stop(); timer.Print();
 
    delete tracker;
    esdFile->cd();
    esdTree->Write();
    delete esd;
    delete esdFile;
}

const AliESDVertex *SetMCvertex(const AliRunLoader *rl, AliTracker *tracker) {
    AliGenEventHeader *h=rl->GetHeader()->GenEventHeader();
    TArrayF vtx(3);
    h->PrimaryVertex(vtx);
    cout<<"Vertex "<<vtx[0]<<' '<<vtx[1]<<' '<<vtx[2]<<endl;
    Double_t xyz[]={vtx[0],vtx[1],vtx[2]};
    Double_t ers[]={2.e-4,2.e-4,2.e-2};
    tracker->SetVertex(xyz,ers);
    AliESDVertex *vertex=new AliESDVertex(xyz,ers);
    vertex->SetNContributors(33); // Set some dummy number of contributors
    return vertex;
}

