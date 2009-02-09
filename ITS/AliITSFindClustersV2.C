#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliITSInitGeometry.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliCDBManager.h"
#include "AliITSDetTypeRec.h"
#include "AliGeomManager.h"

#endif

/*
$Id$ 
*/
/***************************************************************
 *  This macro performs the ITS local reconstruction           *
 *  It is intended for special purposes and tests.             *
 *  The reccomended way to reconstruct ITs is via the          *
 *  class STEER/AliReconstruction                              *
 *  Present version: M.Masera - Previous version: J. Belikov   *
 ***************************************************************/

void Reconstruct(AliRunLoader* runLoader,Option_t *opt);
void Reconstruct(AliRunLoader* runLoader, AliRawReader *rawreader,Option_t *opt);

Int_t AliITSFindClustersV2(char *inputRawData = NULL,TString filename="galice.root",Option_t *opt="All"){
  // if kine tree is available MC labels are used for rec points
  // set opt equal to "SPD" or to "SDD" or to "SSD" do build
  // rec points for individual subdetectors 
  if (gAlice) {
    delete AliRunLoader::Instance();
    delete gAlice;
    gAlice = NULL;
  }

  // Get geometry
  AliGeomManager::LoadGeometry("geometry.root");

  //Get Run loader and ITS loader - set geometry
  AliRunLoader* rl = AliRunLoader::Open(filename.Data());

  AliITSInitGeometry initgeom;
  AliITSgeom *geom = initgeom.CreateAliITSgeom();
  printf("Geometry name: %s \n",(initgeom.GetGeometryName()).Data());
  AliITSLoader* loader = static_cast<AliITSLoader*>(rl->GetLoader("ITSLoader"));
  if (!loader) {
    Error("Init", "ITS loader not found");
    return -1;
  }
  loader->SetITSgeom(geom);

  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number 0\n");
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(0);
  }
  else {
    printf("Using deafult storage \n");
  }

  AliRawReader *rawreader = NULL;
  TString fileRaw(inputRawData);
  if(!fileRaw.IsNull()){
    if (fileRaw.EndsWith(".root")) {
      cout<<"Raw data format - ROOT file \n"; 
      rawreader = new AliRawReaderRoot(fileRaw); // root format
    }
    else {
      cout<<"Raw data format - DATE file \n";
      rawreader = new AliRawReaderDate(fileRaw);  // DATE format
    }
    //    if (!fEquipIdMap.IsNull() && fRawReader)fRawReader->LoadEquipmentIdsMap(fEquipIdMap);
    Reconstruct(rl,rawreader,opt);
  }
  else {
    cout<< "Starting from DIGITS \n";
    Reconstruct(rl,opt);
  }

  return 0;
}

//________________________________________________________________________
void Reconstruct(AliRunLoader* runLoader,Option_t *opt){
// reconstruct clusters starting from DIGITS
// MC truth if available is used to label clusters according to the particles
// originating them

  AliITSLoader* loader = 
                static_cast<AliITSLoader*>(runLoader->GetLoader("ITSLoader"));
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetITSgeom(loader->GetITSgeom());
  rec->SetDefaults();

  runLoader->LoadKinematics();
  TTree *trK=(TTree*)runLoader->TreeK();
  if(trK){
    cout<<"kine tree found - MC labels will be used in RP \n";
    if(runLoader->LoadgAlice())gAlice = runLoader->GetAliRun();
  }
  else{
    cout<<"kine tree not found - MC labels will not b used\n";
  }

  Int_t nEvents = runLoader->GetNumberOfEvents();
  // loop on the events

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    runLoader->GetEvent(iEvent);
    cout<<">>>>>>>   Processing event number "<<iEvent+1<<endl;
    loader->LoadRecPoints("update");
    loader->CleanRecPoints();
    loader->MakeRecPointsContainer();
    loader->LoadDigits("read");
    TTree *tR = loader->TreeR();
    TTree *tD = loader->TreeD();
    if(!tR){
      cout<<"Tree R pointer not found - Abort \n";
      break;
    }
    rec->SetTreeAddressD(tD);
    rec->MakeBranch(tR,"R");
    rec->SetTreeAddressR(tR);
    rec->DigitsToRecPoints(tD,tR,0,opt,kTRUE);    
    loader->WriteRecPoints("OVERWRITE");
    loader->UnloadRecPoints();
    loader->UnloadDigits();
    runLoader->UnloadKinematics();
  }

}

void Reconstruct(AliRunLoader* runLoader, AliRawReader *rawreader,Option_t *opt){
// reconstruct clusters starting from raw data (root or DATE format)
// MC truth if available is used to label clusters according to the particles
// originating them

  AliITSLoader* loader = static_cast<AliITSLoader*>(runLoader->GetLoader("ITSLoader"));
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetITSgeom(loader->GetITSgeom());
  rec->SetDefaults();
  // direct clusterfinding starting from raw data is implemented only 
  // in AliITSClusterFinderV2
  rec->SetDefaultClusterFindersV2(kTRUE);

  runLoader->LoadKinematics();
  TTree *trK=(TTree*)runLoader->TreeK();
  if(trK){
    cout<<"kine tree found - MC labels will be used in RP \n";
    if(runLoader->LoadgAlice())gAlice = runLoader->GetAliRun();
  }
  Int_t nEvents = runLoader->GetNumberOfEvents();
  rawreader->RewindEvents();
  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    rawreader->NextEvent();
    runLoader->GetEvent(iEvent);
    cout<<">>>>>>>   Processing event number: "<<iEvent<<endl;
    loader->LoadRecPoints("update");
    loader->CleanRecPoints();
    loader->MakeRecPointsContainer();
    TTree *tR = loader->TreeR();
    if(!tR){
      cout<<"Tree R pointer not found - Abort \n";
      break;
    }
    rec->DigitsToRecPoints(rawreader,tR,opt);
    rec->ResetRecPoints();
    rec->ResetClusters();    
    loader->WriteRecPoints("OVERWRITE");
    loader->UnloadRecPoints();
    loader->UnloadDigits();
    runLoader->UnloadKinematics();
  }

}
