#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TObjArray.h>
#include <TGeoManager.h>
#include <TProfile.h>

#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"

#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITSRecPoint.h"

#include "AliITSTrackleterSPDEff.h"
#include "AliITSPlaneEffSPD.h"

#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliStack.h>

#endif

void EvaluateSPDEffWithTracklets(Char_t* dir=".", Bool_t mc=kTRUE, Bool_t bckg=kFALSE,
                                 TString cdburi="") {

  Char_t str[256];

  AliITSTrackleterSPDEff* trackleterSPDEff = new AliITSTrackleterSPDEff();
// outer layer (estrapolation)
  trackleterSPDEff->SetPhiWindow(0.07);
  trackleterSPDEff->SetZetaWindow(0.4);
// inner layer (interpolation)
  trackleterSPDEff->SetPhiWindowL1(0.10);
  trackleterSPDEff->SetZetaWindowL1(0.6);
//
  trackleterSPDEff->SetUpdateOncePerEventPlaneEff();
// Study the residual background: reflect outer RecPoints
  if(bckg) trackleterSPDEff->SetReflectClusterAroundZAxisForLayer(1,kTRUE);
//
// this special setting for MC
  if(mc) trackleterSPDEff->SetMC();
  if(trackleterSPDEff->GetMC()) trackleterSPDEff->SetUseOnlyStableParticle();
//
// this for having histograms (both from base class and the new ones)
  trackleterSPDEff->SetHistOn();
//
//
  const Int_t minCont=3;
  const Bool_t VtxMC=kFALSE;
//
  const Bool_t misalign=kTRUE;
//
  // Defining pointers
  AliRunLoader* runLoader;

    if (gAlice) {
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice=0;
    }

    sprintf(str,"%s/galice.root",dir);
    runLoader = AliRunLoader::Open(str);
    runLoader->LoadgAlice();
    gAlice = runLoader->GetAliRun();

    runLoader->LoadKinematics("read");
    runLoader->LoadTrackRefs("read");
    Int_t retval = runLoader->LoadHeader();
    if (retval){
      cerr<<"LoadHeader returned error"<<endl;
      return;
    }

    // open the new ESD file
    sprintf(str,"%s/AliESDs.root",dir);

    TFile inFile(str, "READ");
    TTree *esdTree = (TTree*)inFile.Get("esdTree");
    AliESDEvent *esd = new AliESDEvent();
    esd->ReadFromTree(esdTree);

    // Set OfflineConditionsDataBase if needed
    AliCDBManager* man = AliCDBManager::Instance();
    if (cdburi.Length() > 0) {
      printf("Default CDB storage is set to: %s\n", cdburi.Data());
      man->SetDefaultStorage(cdburi);
      man->SetRun(0);
    }
    if (!man->IsDefaultStorageSet()) {
      printf("Setting a local default CDB storage and run number 0\n");
      man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      man->SetRun(0);
    }
    // retrives geometry
    if (!gGeoManager) {
      sprintf(str,"%s/geometry.root",dir);
      AliGeomManager::LoadGeometry(str);
    }
    // apply misalignement
    if (misalign) AliGeomManager::ApplyAlignObjsFromCDB("ITS"); 

    AliITSLoader* ITSloader =  (AliITSLoader*) runLoader->GetLoader("ITSLoader");
    if (!ITSloader){
      cerr<<"ITS loader not found"<<endl;
      return;
    }
    ITSloader->LoadRecPoints("read");

    // getting number of events
    Int_t nEvents = (Int_t)runLoader->GetNumberOfEvents();

    // loop over events
    for (Int_t iev=0; iev<nEvents; iev++) {

      runLoader->GetEvent(iev);
      // read events
      esdTree->GetEvent(iev);

      // get the ESD vertex
      const AliESDVertex* vtxESD = esd->GetVertex();
      Double_t esdvtx[3];
      vtxESD->GetXYZ(esdvtx);
      Int_t ncont=vtxESD->GetNContributors();
      if(ncont <= minCont) continue;
      Float_t ESDvtx[3];
      ESDvtx[0]=esdvtx[0];
      ESDvtx[1]=esdvtx[1];
      ESDvtx[2]=esdvtx[2];

      // get the MC vertex
      TArrayF vertex(3);
      runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);

     // Read the generated particles
     AliStack *pStack=0x0; TTree *tRefTree=0x0;
     if (trackleterSPDEff->GetMC()) {
       pStack= runLoader->Stack();
       tRefTree= runLoader->TreeTR();
     }

     TTree *itsClusterTree = ITSloader->TreeR();

     if(!VtxMC) {
      if (ESDvtx[2]!=0.) {
        if(trackleterSPDEff->GetMC()) trackleterSPDEff->Reconstruct(itsClusterTree, ESDvtx, ESDvtx, pStack,t RefTree);
        else trackleterSPDEff->Reconstruct(itsClusterTree, ESDvtx, ESDvtx); }
     }
     else {
       Float_t vtx[3]={0.,0.,vertex[2]};
       if(trackleterSPDEff->GetMC()) trackleterSPDEff->Reconstruct(itsClusterTree, vtx, vtx, pStack,tRefTree );
     }

   } // end loop over events

   runLoader->UnloadAll();
   delete runLoader;

if(trackleterSPDEff->GetMC()) trackleterSPDEff->SavePredictionMC("TrackletsMCpred.root");
if(!bckg && !trackleterSPDEff->WriteHistosToFile()) printf("cannot write histos to file \n");
//trackleterSPDEff->GetPlaneEff()->WriteIntoCDB();
const char* name="AliITSPlaneEffSPDtracklet.root";
TFile* pefile = TFile::Open(name, "RECREATE");
Int_t nb=trackleterSPDEff->GetPlaneEff()->Write();
if(nb>0) printf("Writing PlaneEfficiency to file %s\n",name);
pefile->Close();
return;
}
