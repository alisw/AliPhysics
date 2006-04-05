#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "AliRunLoader.h"
#include "AliESD.h"
#include "AliRun.h"

#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITSMultReconstructor.h"

#endif

void testITSMultReco(Char_t* dir = ".") {

  Char_t str[256];

  // ########################################################
  // defining pointers
  AliRunLoader* runLoader;
  TFile* esdFile = 0;
  TTree* esdTree = 0;
  AliESD* esd = 0;

  // #########################################################
  // setup galice and runloader

  if (gAlice) {
    delete gAlice->GetRunLoader();
    delete gAlice;
    gAlice=0;
  }

  sprintf(str,"%s/galice.root",dir);
  runLoader = AliRunLoader::Open(str);
  if (runLoader == 0x0) {
    cout << "Can not open session"<<endl;
    return;
  }
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();
  runLoader->LoadKinematics();
  runLoader->LoadHeader();

  // #########################################################
  // open esd file and get the tree

  // close it first to avoid memory leak
  if (esdFile)
    if (esdFile->IsOpen())
      esdFile->Close();

  sprintf(str,"%s/AliESDs.root",dir);
  esdFile = TFile::Open(str);
  esdTree = (TTree*)esdFile->Get("esdTree");
  TBranch * esdBranch = esdTree->GetBranch("ESD");
  esdBranch->SetAddress(&esd);


  // #########################################################
  // setup its stuff

  AliITS* its=(AliITS*)runLoader->GetAliRun()->GetDetector("ITS");
  if (!its) {
    cout << " Can't get the ITS!" << endl;
    return ;
  }
  AliITSgeom* itsGeo=its->GetITSgeom();
  if (!itsGeo) {
    cout << " Can't get the ITS geometry!" << endl;
    return ;
  }
  AliITSLoader* itsLoader = (AliITSLoader*)runLoader->GetLoader("ITSLoader");
  if (!itsLoader) {
    cout << " Can't get the ITS loader!" << endl;
    return ;
  }
  itsLoader->LoadRecPoints("read");

  // #########################################################
  AliITSMultReconstructor* multReco = new AliITSMultReconstructor();
  multReco->SetGeometry(itsGeo);

  // #########################################################
  // getting number of events

  Int_t nEvents = (Int_t)runLoader->GetNumberOfEvents();
  Int_t nESDEvents = esdBranch->GetEntries();

  if (nEvents!=nESDEvents) {
    cout << " Different number of events from runloader and esdtree!!!" 
	 << nEvents << " / " << nESDEvents << endl;
    return;
  }

  // ########################################################
  // loop over number of events
  cout << nEvents << " event(s) found in the file set" << endl;
  for(Int_t i=0; i<nEvents; i++) {
    
    cout << "-------------------------" << endl << " event# " << i << endl;
    
    runLoader->GetEvent(i);
    esdBranch->GetEntry(i);

    // ########################################################
    // get the EDS vertex
    const AliESDVertex* vtxESD = esd->GetVertex();
    Double_t vtx[3];
    vtxESD->GetXYZ(vtx);   
    Float_t esdVtx[3];
    esdVtx[0] = vtx[0];
    esdVtx[1] = vtx[1];
    esdVtx[2] = vtx[2];
    
    ///#########################################################
    // get ITS clusters 
    TTree* itsClusterTree = itsLoader->TreeR();
    if (!itsClusterTree) {
      cerr<< " Can't get the ITS cluster tree !\n";
      return;
    }
    multReco->SetHistOn(kTRUE);
    multReco->Reconstruct(itsClusterTree, esdVtx, esdVtx);

    
    for (Int_t t=0; t<multReco->GetNTracklets(); t++) {
      
      cout << "  tracklet " << t 
	   << " , theta = " << multReco->GetTracklet(t)[0]
	   << " , phi = " << multReco->GetTracklet(t)[1] << endl; 
    }

  }
 
  TFile* fout = new TFile("out.root","RECREATE");  

  multReco->SaveHists();
  fout->Write();
  fout->Close();


}
