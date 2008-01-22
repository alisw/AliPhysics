// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveVSDCreator.h"

#include <TEveTreeTools.h>

#include <AliStack.h>
#include <AliITSLoader.h>
#include <AliTPCTrackHitsV2.h>
#include <AliPDG.h>
#include <AliHit.h>
#include <AliSimDigits.h>
#include <AliKalmanTrack.h>
#include <AliESDEvent.h>
#include <AliESDv0.h>
#include <AliTPCclusterMI.h>
#include <AliTPCClustersRow.h>
#include <AliITS.h>
#include <AliITSclusterV2.h>
#include <AliTrackReference.h>
#include <AliESDkink.h>
#include <AliESDtrack.h>

#include <AliRun.h>
#include <AliTPCParam.h>

#include <TSystem.h>
#include <TFile.h>
#include <TError.h>

//______________________________________________________________________________
// AliEveVSDCreator
//

ClassImp(AliEveVSDCreator)

AliEveVSDCreator::AliEveVSDCreator(const Text_t* name, const Text_t* title) :
  TEveVSD(name, title),

  mKineType (kKT_Standard),
  mDataDir  ("."),
  mEvent    (0),

  mTPCHitRes (2),
  mTRDHitRes (2),

  mDebugLevel (0),
  mGenInfoMap (),

  pRunLoader (0)
{
  // Particles not in ROOT's PDG database occuring in ALICE
  AliPDG::AddParticlesToPdgDataBase();
  {
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    // const Int_t kspe=50000000;
    const Int_t kion=10000000;

    const Double_t kAu2Gev=0.9314943228;
    const Double_t khSlash = 1.0545726663e-27;
    const Double_t kErg2Gev = 1/1.6021773349e-3;
    const Double_t khShGev = khSlash*kErg2Gev;
    const Double_t kYear2Sec = 3600*24*365.25;

    pdgDB->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,
		       0,1,"Ion",kion+10020);
    pdgDB->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE,
		       khShGev/(12.33*kYear2Sec),1,"Ion",kion+10030);
    pdgDB->AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE,
		       khShGev/(12.33*kYear2Sec),2,"Ion",kion+20040);
    pdgDB->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,
		       0,2,"Ion",kion+20030);
  }

  // AliKalmanTrack::SetConvConst(1);
}

/******************************************************************************/

void AliEveVSDCreator::CreateVSD(const Text_t* data_dir, Int_t event,
			   const Text_t* vsd_file)
{
  static const TEveException eH("AliEveVSDCreator::CreateVSD ");

  mDataDir = data_dir;
  mEvent   = event;

  string galice_file (Form("%s/galice.root", mDataDir.Data()));

  if(mDebugLevel > 0)
    printf("%s opening %s \n", eH.Data(), galice_file.c_str());

  if(gSystem->AccessPathName(galice_file.c_str(), kReadPermission)) {
    throw(eH + "Can not read file '" + galice_file + "'.");
  }
  pRunLoader = AliRunLoader::Open(galice_file.c_str());
  if(pRunLoader == 0)
    throw(eH + "AliRunLoader::Open failed.");

  pRunLoader->LoadgAlice();
  Int_t status = pRunLoader->GetEvent(mEvent);
  if(status)
    throw(eH + Form("GetEvent(%d) failed, exit code %s.", mEvent, status));

  if(mDebugLevel > 0)
    printf("%s open seems ok. Now loading sim data.\n", eH.Data());

  pRunLoader->LoadHeader();
  pRunLoader->LoadKinematics();
  pRunLoader->LoadTrackRefs();
  pRunLoader->LoadHits();

  // GledNS::PushFD();

  if(mDebugLevel > 0)
    printf("%s opening output TEveVSD.\n", eH.Data());

  TFile* file = TFile::Open(vsd_file, "RECREATE", "ALICE VisualizationDataSummary");
  fDirectory = new TDirectoryFile("Event0", "");

  if(mDebugLevel > 0)
    printf("%s creating trees now ...\n", eH.Data());

  CreateTrees();

  if(mDebugLevel > 0)
    printf("%s trees created, closing files.\n", eH.Data());

  file->Write();
  file->Close();
  delete file;
  fDirectory =0;

  //GledNS::PopFD();

  // clean after the TEveVSD data was sucessfuly written
  fTreeK      = 0;
  fTreeH      = 0;
  //fTreeTR     = 0;
  fTreeC      = 0;
  fTreeV0     = 0;
  fTreeKK     = 0;
  fTreeR      = 0;
  fTreeGI     = 0;

  pRunLoader->UnloadAll();
  delete pRunLoader;
  if(gAlice) {
    delete gAlice; gAlice = 0;
  }
  pRunLoader = 0;

  if(mDebugLevel > 0)
    printf("%s all done.\n", eH.Data());
}

void AliEveVSDCreator::CreateTrees()
{
  static const TEveException eH("AliEveVSDCreator::CreateTrees ");

  if(fDirectory == 0)
    throw(eH + "output directory not set.");

  try {
    if(mDebugLevel > 1)
      printf("%sConvertKinematics.\n", eH.Data());
    ConvertKinematics();
  } catch(TEveException& exc) { Warning(eH, exc); }

  try {
    if(mDebugLevel > 1)
      printf("%sConvertHits.\n", eH.Data());
    ConvertHits();
  } catch(TEveException& exc) { Warning(eH, exc); }

  try {
    if(mDebugLevel > 1)
      printf("%sConvertClusters.\n", eH.Data());
    ConvertClusters();
  } catch(TEveException& exc) { Warning(eH, exc); }

  try {
    if(mDebugLevel > 1)
      printf("%sConvertRecTracks.\n", eH.Data());
    ConvertRecTracks();
  } catch(TEveException& exc) {
    Warning(exc, "skipping AliEveV0 extraction.");
    goto end_esd_processing;
  }

  try {
    if(mDebugLevel > 1)
      printf("%sConvertV0.\n", eH.Data());
    ConvertV0();
  } catch(TEveException& exc) { Warning(eH, exc); }

  try {
    if(mDebugLevel > 1)
      printf("%sConvertKinks.\n", eH.Data());
    ConvertKinks();
  } catch(TEveException& exc) { Warning(eH, exc); }

end_esd_processing:

  try {
    if(mDebugLevel > 1)
      printf("%sConvertGenInfo.\n", eH.Data());
    ConvertGenInfo();
  } catch(TEveException& exc) { Warning(eH, exc); }

  return;
}

/******************************************************************************/
// Kinematics
/******************************************************************************/

void AliEveVSDCreator::ConvertKinematics()
{
  static const TEveException eH("AliEveVSDCreator::ConvertKinematics ");

  if(fTreeK != 0)
    throw (eH + "kinematics already converted");

  AliStack* stack = pRunLoader->Stack();
  if(stack == 0)
    throw(eH + "stack is null.");

  fDirectory->cd();
  fTreeK = new TTree("Kinematics", "TParticles sorted by Label");

  Int_t nentries = stack->GetNtrack();
  std::vector<TEveMCTrack>  vmc(nentries);
  for (Int_t idx=0; idx<nentries; idx++) {
    TParticle* tp = stack->Particle(idx);
    vmc[idx]        = *tp;
    vmc[idx].fLabel = idx;
  }

  // read track refrences
  // functionality now in AliEveKineTools.
  /*
  TTree* fTreeTR =  pRunLoader->TreeTR();

  if(fTreeTR == 0) {
    Warning(eH, "no TrackRefs; some data will not be available.");
  } else {
    TClonesArray* RunArrayTR = 0;
    fTreeTR->SetBranchAddress("AliRun", &RunArrayTR);

    Int_t nPrimaries = (Int_t) fTreeTR->GetEntries();
    for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
      // printf("T0 fTreeTR->GetEntry(%d) \n",iPrimPart);
      fTreeTR->GetEntry(iPrimPart);
      // printf("END fTreeTR->GetEntry(%d) \n",iPrimPart);

      for (Int_t iTrackRef = 0; iTrackRef < RunArrayTR->GetEntriesFast(); iTrackRef++) {
	AliTrackReference *trackRef = (AliTrackReference*)RunArrayTR->At(iTrackRef);
	Int_t track = trackRef->GetTrack();
	if(track < nentries && track > 0){
	  TEveMCTrack& mct = vmc[track];
	  if(trackRef->TestBit(kNotDeleted)) {
	    mct.decayed   = true;
	    mct.t_decay   = trackRef->GetTime();
	    mct.V_decay.x = trackRef->X();
	    mct.V_decay.y = trackRef->Y();
	    mct.V_decay.z = trackRef->Z();
	    mct.P_decay.x = trackRef->Px();
	    mct.P_decay.y = trackRef->Py();
	    mct.P_decay.z = trackRef->Pz();
	    if(TMath::Abs(mct.GetPdgCode()) == 11)
	      mct.decayed = false; // a bug in TreeTR
	  }
	}
      }
    }
  }
  */

  fTreeK->Branch("K", "TEveMCTrack",  &fpK, fBuffSize);

  printf("sizeofvmc = %d\n", vmc.size());
  for(std::vector<TEveMCTrack>::iterator k=vmc.begin(); k!=vmc.end(); ++k) {
    TEveMCTrack& mct = *k;
    fK = mct;

    TParticle* m  = &mct;
    Int_t      mi = mct.fLabel;
    int cnt = 0;
    while(m->GetMother(0) != -1) {
      if(cnt > 100) {
	printf("cnt %d mi=%d, mo=%d\n", cnt, mi, m->GetMother(0));
      }
      mi = m->GetMother(0);
      m = &vmc[mi];
      ++cnt;
    }
    fK.fEvaLabel = mi;

    fTreeK->Fill();
  }

  fTreeK->BuildIndex("label");
}

/******************************************************************************/
// Hits
/******************************************************************************/

namespace {

  struct Detector
  {
    const char*   name;
    const char*   hitbranch;
    unsigned char detidx;
  };

  Detector detects[] = {
    { "ITS",  "AliITShit",         0 },
    { "TPC",  "AliTPCTrackHitsV2", 1 },
    { "TRD",  "AliTRDhit",         2 },
    { "TOF",  "AliTOFhit",         3 }
    // { "HMPID", "AliHMPIDhit",        4 },
  };

}

/******************************************************************************/

void AliEveVSDCreator::ConvertHits()
{
  static const TEveException eH("AliEveVSDCreator::ConvertHits ");

  if(fTreeH != 0)
    throw(eH + "hits already converted.");

  fDirectory->cd();
  fTreeH =  new TTree("Hits", "Combined detector hits.");
  fTreeH->Branch("H", "TEveHit", &fpH, fBuffSize);

  std::map<Int_t, Int_t> hmap;
  // parameters for ITS, TPC hits filtering
  Float_t x,y,z, x1,y1,z1;
  Float_t tpc_sqr_res = mTPCHitRes*mTPCHitRes;
  Float_t trd_sqr_res = mTRDHitRes*mTRDHitRes;

  int l=0;
  // load hits from the rest of detectors
  while(detects[l].name != 0) {
    Detector& det = detects[l++];

    switch(det.detidx) {
    case 1: {
      Int_t count = 0;
      TTree* treeh = pRunLoader->GetTreeH(det.name, false);
      if(treeh == 0) {
	Warning(eH, Form("no hits for %s.", det.name));
	continue;
      }
      AliTPCTrackHitsV2 hv2, *_hv2=&hv2;
      treeh->SetBranchAddress("TPC2", &_hv2);
      Int_t np = treeh->GetEntries();
      for(Int_t i=0; i<np; i++){
	treeh->GetEntry(i);
	Int_t eva_idx = np -i -1;
	if (hv2.First() == 0) continue;
        x = y = z = 0;
	do {
	  AliHit* ah = hv2.GetHit();
	  x1 = ah->X(); y1 = ah->Y(); z1 = ah->Z();
	  if ((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) > tpc_sqr_res)
	  {
	    fH.fDetId    = det.detidx;
	    fH.fSubdetId = 0;
	    fH.fLabel    = ah->Track();
	    fH.fEvaLabel = eva_idx;
	    fH.fV.fX = x1; fH.fV.fY = y1; fH.fV.fZ = z1;
	    fTreeH->Fill();
	    hmap[fH.fLabel]++;
	    x = x1; y = y1; z = z1;
	    count++;
	  }
	} while (hv2.Next());
      }
      // printf("%d entries in TPChits \n",count);
      break;
    }
    default: {
      TTree* treeh = pRunLoader->GetTreeH(det.name, false);
      if(treeh == 0) {
	Warning(eH, Form("no hits for %s.", det.name));
	continue;
      }
      TClonesArray *arr = new TClonesArray(det.hitbranch);
      treeh->SetBranchAddress(det.name, &arr);
      Int_t np = treeh->GetEntries();
      // in TreeH files hits are grouped in clones arrays
      // each eva particle has its own clone array
      for (Int_t i=0; i<np; i++) {
	treeh->GetEntry(i);
	Int_t eva_idx = np -i -1;
	Int_t nh=arr->GetEntriesFast();
	x = y = z = 0;
	// printf("%d entry %d hits for primary %d \n", i, nh, eva_idx);
	for (Int_t j=0; j<nh; j++) {
	  AliHit* ali_hit = (AliHit*)arr->UncheckedAt(j);
	  fH.fDetId    = det.detidx;
	  fH.fSubdetId = 0;
	  fH.fLabel     = ali_hit->GetTrack();
	  fH.fEvaLabel = eva_idx;
	  fH.fV.Set(ali_hit->X(), ali_hit->Y(), ali_hit->Z());
	  if(det.detidx == 2) {
	    x1=ali_hit->X();y1=ali_hit->Y();z1=ali_hit->Z();
	    if((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1) < trd_sqr_res) continue;
	    x=x1; y=y1; z=z1;
	  }
	  hmap[fH.fLabel]++;
	  fTreeH->Fill();
	}
      }
      delete arr;
      break;
    } // end default
    } // end switch
  } // end while


  //set geninfo
  for(std::map<Int_t, Int_t>::iterator j=hmap.begin(); j!=hmap.end(); ++j) {
    GetGeninfo(j->first)->fNHits += j->second;
  }
}

/******************************************************************************/
// Clusters
/******************************************************************************/

void AliEveVSDCreator::ConvertClusters()
{
  static const TEveException eH("AliEveVSDCreator::ConvertClusters ");

  if(fTreeC != 0)
    throw(eH + "clusters already converted.");

  fDirectory->cd();
  fTreeC =  new TTree("Clusters", "rec clusters");
  fTreeC->Branch("C", "TEveCluster", &fpC, fBuffSize);

  try {
    ConvertITSClusters();
  } catch(TEveException& exc) { Warning(eH, exc); }

  try {
    ConvertTPCClusters();
  } catch(TEveException& exc) { Warning(eH, exc); }
}

/******************************************************************************/

void AliEveVSDCreator::ConvertTPCClusters()
{
  static const TEveException eH("AliEveVSDCreator::ConvertTPCClusters ");

  auto_ptr<TFile> f
    ( TFile::Open(Form("%s/TPC.RecPoints.root", mDataDir.Data())) );
  if(!f.get())
    throw(eH + "can not open 'TPC.RecPoints.root' file.");

  auto_ptr<TDirectory> d
    ( (TDirectory*) f->Get(Form("AliEveEventManager%d", mEvent)) );
  if(!d.get())
    throw(eH + Form("event directory '%d' not found.", 0));

  auto_ptr<TTree> tree( (TTree*) d->Get("TreeR") );
  if(!tree.get())
    throw(eH + "'TreeR' not found.");

  auto_ptr<AliTPCParam> par( GetTpcParam(eH) );

  AliTPCClustersRow  clrow, *_clrow=&clrow;
  AliTPCclusterMI   *cl;
  _clrow->SetClass("AliTPCclusterMI");
  tree->SetBranchAddress("Segment", &_clrow);

  // count clusters
  Int_t nClusters = 0;
  Int_t n_ent = tree->GetEntries();
  for (Int_t n=0; n<n_ent; n++) {
    tree->GetEntry(n);
    nClusters += _clrow->GetArray()->GetEntriesFast();
  }

  // calculate xyz for a cluster and add it to container
  Double_t x,y,z;
  Float_t cs, sn, tmp;
  std::map<Int_t, Int_t> cmap;

  for (Int_t n=0; n<tree->GetEntries(); n++) {
    tree->GetEntry(n);
    Int_t ncl = _clrow->GetArray()->GetEntriesFast();
    if(ncl > 0) {
      Int_t sec,row;
      par->AdjustSectorRow(_clrow->GetID(),sec,row);
      while (ncl--) {
	if(_clrow->GetArray()) {
	  // cl = new AliTPCclusterMI(*(AliTPCclusterMI*)_clrow->GetArray()->UncheckedAt(ncl));
	  cl = (AliTPCclusterMI*)_clrow->GetArray()->UncheckedAt(ncl);
          if(cl->GetLabel(0) >= 0)
	  {
	    x = par->GetPadRowRadii(sec,row); y = cl->GetY(); z = cl->GetZ();
	    par->AdjustCosSin(sec,cs,sn);
	    tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp;

	    fC.fDetId    = 1;
	    fC.fSubdetId = 0;
	    fC.fLabel[0] = cl->GetLabel(0);
	    fC.fLabel[1] = cl->GetLabel(1);
	    fC.fLabel[2] = cl->GetLabel(2);
	    fC.fV.Set(x, y, z);

	    fTreeC->Fill();
	    {
	      int i = 0;
	      while(i < 3 && fC.fLabel[i])
		cmap[fC.fLabel[i++]]++;
	    }
	  }
	}
      }
    }
  }
  //set geninfo
  for(std::map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j) {
    GetGeninfo(j->first)->fNClus += j->second;
  }
}

/******************************************************************************/

void AliEveVSDCreator::ConvertITSClusters()
{
  static const TEveException eH("AliEveVSDCreator::ConvertITSClusters ");

  auto_ptr<TFile> f
    ( TFile::Open(Form("%s/ITS.RecPoints.root", mDataDir.Data())) );
  if(!f.get())
    throw(eH + "can not open 'ITS.RecPoints.root' file.");

  auto_ptr<TDirectory> d
    ( (TDirectory*) f->Get(Form("AliEveEventManager%d", mEvent)) );
  if(!d.get())
    throw(eH + Form("event directory '%d' not found.", 0));

  auto_ptr<TTree> tree( (TTree*) d->Get("TreeR") );
  if(!tree.get())
    throw(eH + "'TreeR' not found.");

  AliITSLoader* ITSld =  (AliITSLoader*) pRunLoader->GetLoader("ITSLoader");
  //AliITS* pITS = ITSld->GetITS();
  AliITSgeom* geom = ITSld->GetITSgeom();
  //AliITSgeom* geom = new AliITSgeom();
  //geom->ReadNewFile("/home/aljam/ITSgeometry.det");

  //printf("alice ITS geom %p \n",geom );

  if(!geom)
    throw(eH + "can not find ITS geometry");

  TClonesArray *arr = new TClonesArray("AliITSclusterV2");
  tree->SetBranchAddress("Clusters", &arr);
  Int_t nmods = tree->GetEntries();
  Float_t gc[3];
  std::map<Int_t, Int_t> cmap;

  for (Int_t mod=0; mod<nmods; mod++) {
    tree->GetEntry(mod);
    Int_t nc=arr->GetEntriesFast();
    for (Int_t j=0; j<nc; j++) {
      AliITSclusterV2* recp = (AliITSclusterV2*)arr->UncheckedAt(j);

      Double_t rot[9];
      geom->GetRotMatrix(mod,rot);
      Int_t lay,lad,det;
      geom->GetModuleId(mod,lay,lad,det);
      Float_t tx,ty,tz;
      geom->GetTrans(lay,lad,det,tx,ty,tz);

      Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
      Double_t phi1=TMath::Pi()/2+alpha;
      if(lay==1) phi1+=TMath::Pi();

      Float_t cp=TMath::Cos(phi1), sp=TMath::Sin(phi1);
      Float_t  r=tx*cp+ty*sp;
      gc[0] = r*cp - recp->GetY()*sp;
      gc[1] = r*sp + recp->GetY()*cp;
      gc[2] = recp->GetZ();

      fC.fDetId    = 0;
      fC.fSubdetId = 0;
      fC.fLabel[0] = recp->GetLabel(0);
      fC.fLabel[1] = recp->GetLabel(1);
      fC.fLabel[2] = recp->GetLabel(2);
      fC.fV.fX     = r*cp - recp->GetY()*sp;
      fC.fV.fY     = r*sp + recp->GetY()*cp;
      fC.fV.fZ     = recp->GetZ();
      fTreeC->Fill();
      { int i = 0;
	while(i < 3 && fC.fLabel[i])
	  cmap[fC.fLabel[i++]]++;
      }
    }

    for(std::map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j) {
      GetGeninfo(j->first)->fNClus += j->second;
    }
  }
  delete arr;
}

/******************************************************************************/
// ESD
/******************************************************************************/

void AliEveVSDCreator::ConvertRecTracks()
{
  static const TEveException eH("AliEveVSDCreator::ConvertRecTracks ");

  if(fTreeR != 0)
    throw(eH + "tracks already converted.");

  fDirectory->cd();
  fTreeR =  new TTree("RecTracks", "rec tracks");

  fTreeR->Branch("R", "TEveRecTrack", &fpR, 512*1024,1);

  TFile f(Form("%s/AliESDs.root", mDataDir.Data()));
  if(!f.IsOpen())
    throw(eH + "no AliESDs.root file.");

  TTree* tree = (TTree*) f.Get("esdTree");
  if (tree == 0)
    throw(eH + "no esdTree.");


  AliESDEvent *fEvent= new AliESDEvent();
  fEvent->ReadFromTree(tree);
  tree->GetEntry(mEvent);
  if(fEvent->GetAliESDOld())fEvent->CopyFromOldESD();


  // reconstructed tracks
  AliESDtrack* esd_t;
  Double_t     dbuf[3];
  for (Int_t n=0; n<fEvent->GetNumberOfTracks(); n++) {
    esd_t = fEvent->GetTrack(n);

    fR.fLabel  = esd_t->GetLabel();
    fR.fStatus = (Int_t) esd_t->GetStatus();
    fR.fSign   = (Int_t) esd_t->GetSign();
    esd_t->GetXYZ(dbuf);    fR.fV.Set(dbuf);
    esd_t->GetPxPyPz(dbuf); fR.fP.Set(dbuf);
    Double_t ep = esd_t->GetP();
    fR.fBeta = ep/TMath::Sqrt(ep*ep + TMath::C()*TMath::C()*esd_t->GetMass()*esd_t->GetMass());
    fTreeR->Fill();
  }
  fTreeR->BuildIndex("label");
  delete fEvent;
}

/******************************************************************************/

void AliEveVSDCreator::ConvertV0()
{
  static const TEveException eH("AliEveVSDCreator::ConvertV0 ");

  if(fTreeV0 != 0)
    throw(eH + "AliEveV0 already converted.");

  fDirectory->cd();
  fTreeV0 =  new TTree("AliEveV0", "AliEveV0 points");

  fTreeV0->Branch("AliEveV0", "TEveRecV0", &fpV0, 512*1024,1);

  TFile f(Form("%s/AliESDs.root", mDataDir.Data()));
  if(!f.IsOpen()){
    throw(eH + "no AliESDs.root file.");
  }

  TTree* tree = (TTree*) f.Get("esdTree");
  if (tree == 0)
    throw(eH + "no esdTree.");

  AliESDEvent *fEvent= new AliESDEvent();
  fEvent->ReadFromTree(tree);
  tree->GetEntry(mEvent);
  if(fEvent->GetAliESDOld())fEvent->CopyFromOldESD();

  for (Int_t n =0; n< fEvent->GetNumberOfV0s(); n++)
  {
    AliESDv0    *av     = fEvent->GetV0(n);
    AliESDtrack *trackN = fEvent->GetTrack(av->GetNindex()); // negative daughter
    AliESDtrack *trackP = fEvent->GetTrack(av->GetPindex()); // positive daughter

    Double_t pos[3];

    fV0.fStatus = av->GetStatus();
    // Point of closest approach
    av->GetXYZ(pos[0],pos[1],pos[2]);
    fV0.fVCa.fX = pos[0];
    fV0.fVCa.fY = pos[1];
    fV0.fVCa.fZ = pos[2];
    // set birth vertex of neutral particle
    av->GetXYZ(pos[0], pos[1], pos[2]);
    fV0.fV0Birth.Set(pos);

    // momentum and position of negative particle
    av->GetParamN()->GetPxPyPz(pos);
    fV0.fPNeg.Set(pos);
    av->GetParamN()->GetXYZ(pos);
    fV0.fVNeg.Set(pos);

    // momentum and position of positive particle
    av->GetParamP()->GetPxPyPz(pos);
    fV0.fPPos.Set(pos);
    av->GetParamP()->GetXYZ(pos);
    fV0.fVPos.Set(pos);

    fV0.fLabel = 0; // !!!! mother label unknown
    fV0.fPdg   = av->GetPdgCode();

    // daughter indices
    fV0.fDLabel[0] = TMath::Abs(trackN->GetLabel());
    fV0.fDLabel[1] = TMath::Abs(trackP->GetLabel());

    // printf("AliEveV0 convert labels(%d,%d) index(%d,%d)\n",
    //	   fV0.d_label[0],  fV0.d_label[1],
    //	   av->GetNIndex(), av->GetPIndex());

    fTreeV0->Fill();
  }
  // if(fEvent->GetNumberOfV0s()) fTreeV0->BuildIndex("label");
  delete fEvent;
}

/******************************************************************************/

void AliEveVSDCreator::ConvertKinks()
{
  static const TEveException eH("AliEveVSDCreator::ConvertKinks ");

  if(fTreeKK != 0)
    throw(eH + "Kinks already converted.");

  fDirectory->cd();
  fTreeKK =  new TTree("Kinks", "ESD Kinks");

  fTreeKK->Branch("KK", "TEveRecKink", &fpKK, fBuffSize);

  TFile f(Form("%s/AliESDs.root", mDataDir.Data()));
  if(!f.IsOpen()){
    throw(eH + "no AliESDs.root file.");
  }

  TTree* tree = (TTree*) f.Get("esdTree");
  if (tree == 0)
    throw(eH + "no esdTree.");


  AliESDEvent *fEvent= new AliESDEvent();
  fEvent->ReadFromTree(tree);
  tree->GetEntry(mEvent);
  if(fEvent->GetAliESDOld())fEvent->CopyFromOldESD();


  //  printf("CONVERT KINK Read %d entries in tree kinks \n",  fEvent->GetNumberOfKinks());
  for (Int_t n =0; n< fEvent->GetNumberOfKinks(); n++) {
    AliESDkink* kk = fEvent->GetKink(n);

    Double_t pos[3];

    fKK.fLabel  = kk->GetLabel(0);
    fKK.fStatus = Int_t(kk->GetStatus(1) << 8 + kk->GetStatus(2));

    // reconstructed kink position
    fKK.fLabelSec = kk->GetLabel(1);
    fKK.fVKink.Set(kk->GetPosition());

    const AliExternalTrackParam& tp_mother = kk->RefParamMother();
    // momentum and position of mother
    tp_mother.GetPxPyPz(pos);
    fKK.fP.Set(pos);
    tp_mother.GetXYZ(pos);
    fKK.fV.Set(pos);
    const Double_t* par =  tp_mother.GetParameter();
    // printf("KINK Pt %f, %f \n",1/tp_mother.Pt(),par[4] );
    fKK.fSign = (par[4] < 0) ? -1 : 1;

    const AliExternalTrackParam& tp_daughter = kk->RefParamDaughter();
    // momentum and position of daughter
    tp_daughter.GetPxPyPz(pos);
    fKK.fPSec.Set(pos);
    tp_daughter.GetXYZ(pos);
    fKK.fVEnd.Set(pos);

    fTreeKK->Fill();
  }
  if(fEvent->GetNumberOfKinks()) fTreeKK->BuildIndex("label");
  delete fEvent;
}
/******************************************************************************/
// TEveMCRecCrossRef
/******************************************************************************/

void AliEveVSDCreator::ConvertGenInfo()
{
  static const TEveException eH("AliEveVSDCreator::ConvertGenInfo ");

  if(fTreeGI != 0)
    throw(eH + "GI already converted.");

  fDirectory->cd();
  fTreeGI = new TTree("TEveMCRecCrossRef", "Objects prepared for cross querry");

  TEveMCRecCrossRef::Class()->IgnoreTObjectStreamer(true);
  fTreeGI->Branch("GI", "TEveMCRecCrossRef",  &fpGI, fBuffSize);
  fTreeGI->Branch("K.", "TEveMCTrack",  &fpK);
  fTreeGI->Branch("R.", "TEveRecTrack", &fpR);

  for (std::map<Int_t, TEveMCRecCrossRef*>::iterator j=mGenInfoMap.begin(); j!=mGenInfoMap.end(); ++j) {
    fGI        = *(j->second);
    fGI.fLabel = j->first;
    fTreeK->GetEntry(j->first);

    if (fTreeR) {
      Int_t re = fTreeR->GetEntryNumberWithIndex(j->first);
      if(re != -1)
	fGI.fIsRec = true;
    }
    //    Int_t has_v0 =  fTreeV0->GetEntryNumberWithIndex(j->first);
    //if (has_v0 != -1)
    //  fGI.has_AliEveV0 = true;
    if (fTreeKK) {
      Int_t has_kk =  fTreeKK->GetEntryNumberWithIndex(j->first);
      if (has_kk != -1)
	fGI.fHasKink = true;
    }
    fTreeGI->Fill();
  }
  mGenInfoMap.clear();
}

/******************************************************************************/
/******************************************************************************/
// Protected methods
/******************************************************************************/
/******************************************************************************/

AliTPCParam* AliEveVSDCreator::GetTpcParam(const TEveException& eh)
{
  auto_ptr<TFile> fp( TFile::Open(Form("%s/galice.root", mDataDir.Data())) );
  if(!fp.get())
    throw(eh + "can not open 'galice.root' file.");
  AliTPCParam* par = (AliTPCParam *) fp->Get("75x40_100x60_150x60");
  if(!par)
    throw(eh + "TPC data not found.");
  return par;
}



TEveMCRecCrossRef* AliEveVSDCreator::GetGeninfo(Int_t label)
{
  // printf("get_geninfo %d\n", label);
  TEveMCRecCrossRef* gi;
  std::map<Int_t, TEveMCRecCrossRef*>::iterator i = mGenInfoMap.find(label);
  if(i == mGenInfoMap.end()) {
    gi =  new TEveMCRecCrossRef();
    mGenInfoMap[label] = gi;
  } else {
    gi = i->second;
  }
  return gi;
}
