#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>

#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliESDEvent.h"

#include "AliITSLoader.h"
#include "AliITSMultReconstructor.h"
#include "AliGeomManager.h"                                                                                 

#endif

  void testITSMultReco(Char_t* dir = ".") {

  Char_t fileName[256];

  // defining pointers
  AliRunLoader* runLoader;
  TTree* esdTree = 0;
  AliESDEvent* esd = new AliESDEvent();

  // get runloader

  if (gAlice) {
    delete AliRunLoader::Instance();
    delete gAlice;
    gAlice=0;
  }

  sprintf(fileName,"%s/galice.root",dir);
  runLoader = AliRunLoader::Open(fileName);
/*  if (runLoader == 0x0) {
    cout << "Can not open session"<<endl;
    return;
  }*/
  
  // get geometry (here geometry.root is used, change it if needed)
  if (!gGeoManager) {
    sprintf(fileName,"%s/geometry.root",dir);
    AliGeomManager::LoadGeometry(fileName);
  }

  // open the ESD file and get the tree

  sprintf(fileName,"%s/AliESDs.root",dir);
  TFile esdFile(fileName, "READ");
  esdTree = (TTree*)esdFile.Get("esdTree");
  esd->ReadFromTree(esdTree);

  // setup ITS stuff

  AliITSLoader* itsLoader = (AliITSLoader*)runLoader->GetLoader("ITSLoader");
  if (!itsLoader) {
    cout << " Can't get the ITS loader!" << endl;
    return ;
  }
  itsLoader->LoadRecPoints("read");

  AliITSMultReconstructor* multReco = new AliITSMultReconstructor();
//  multReco->SetGeometry(itsGeo);

  // getting number of events

  Int_t nEvents = (Int_t)runLoader->GetNumberOfEvents();
  Int_t nESDEvents = esdTree->GetEntries();

  if (nEvents!=nESDEvents) {
    cout << " Different number of events from runloader and esdtree!!!" 
	 << nEvents << " / " << nESDEvents << endl;
    return;
  }

  // loop over number of events
  cout << nEvents << " event(s) found in the file set" << endl;
  for (Int_t iEv=0; iEv<nEvents; ++iEv) {
    
    cout << "-------------------------" << endl << " event# " << iEv << endl;
    
    runLoader->GetEvent(iEv);
    esdTree->GetEvent(iEv);

    // get the ESD vertex

    const AliESDVertex* vtxESD = esd->GetVertex();
    Double_t vtx[3];
    vtxESD->GetXYZ(vtx);   
    Float_t esdVtx[3];
    esdVtx[0] = vtx[0];
    esdVtx[1] = vtx[1];
    esdVtx[2] = vtx[2];
//    cout<<"vertex Z->"<<esdVtx[2]<<endl;    

    // get ITS clusters 

    TTree* itsClusterTree = itsLoader->TreeR();
    if (!itsClusterTree) {
      cerr<< " Can't get the ITS cluster tree !\n";
      return;
    }

    multReco->SetHistOn(kTRUE);
    multReco->Reconstruct(itsClusterTree, esdVtx, esdVtx);

    cout <<"Number of tracklets: "<<multReco->GetNTracklets()<<endl;     
    for (Int_t itr=0; itr<multReco->GetNTracklets(); itr++) {
      
      cout << "  tracklet "    << itr 
	   << " , theta = "    << multReco->GetTracklet(itr)[0]
	   << " , phi = "      << multReco->GetTracklet(itr)[1] 
           << " , DeltaPhi = " << multReco->GetTracklet(itr)[2]<< endl; 

    }
    cout<<""<<endl;
    cout <<"Number of single clusters (inner layer): "<<multReco->GetNSingleClusters()<<endl;     
    for (Int_t iscl=0; iscl<multReco->GetNSingleClusters(); iscl++) {
      
      cout << "  cluster "  << iscl
	   << " , theta = " << multReco->GetCluster(iscl)[0]
	   << " , phi = "   << multReco->GetCluster(iscl)[1] << endl; 
    }

  }
 
  TFile* fout = new TFile("out.root","RECREATE");  

  multReco->SaveHists();
  fout->Write();
  fout->Close();


}
