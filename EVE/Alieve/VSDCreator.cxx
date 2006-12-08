// $Header$

#include "VSDCreator.h"

#include <Reve/TTreeTools.h>

#include <AliStack.h>
#include <AliITSLoader.h>
#include <AliTPCTrackHitsV2.h>
#include <AliPDG.h>
#include <AliHit.h>
#include <AliSimDigits.h>
#include <AliKalmanTrack.h>
#include <AliESD.h>
#include <AliESDv0.h>
#include <AliTPCclusterMI.h>
#include <AliTPCClustersRow.h>
#include <AliITS.h>
#include <AliITSclusterV2.h>
#include <AliTrackReference.h>
#include <AliESDkink.h>

#include <AliRun.h>
#include <AliTPCParam.h>

#include <TSystem.h>
#include <TFile.h>

using namespace Reve;
using namespace Alieve;

using namespace std;

//______________________________________________________________________
// VSDCreator
//

ClassImp(VSDCreator)

VSDCreator::VSDCreator(const Text_t* name, const Text_t* title) :
  VSD(name, title),

  mKineType (KT_Standard),
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

/**************************************************************************/

void VSDCreator::CreateVSD(const Text_t* data_dir, Int_t event,
			   const Text_t* vsd_file)
{
  static const Exc_t eH("VSDCreator::CreateVSD ");

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
    printf("%s opening output VSD.\n", eH.Data());

  TFile* file = TFile::Open(vsd_file, "RECREATE", "ALICE VisualizationDataSummary");
  mDirectory = new TDirectory("Event0", "");

  if(mDebugLevel > 0)
    printf("%s creating trees now ...\n", eH.Data());

  CreateTrees();

  if(mDebugLevel > 0)
    printf("%s trees created, closing files.\n", eH.Data());

  file->Write();
  file->Close();
  delete file; 
  mDirectory =0;

  //GledNS::PopFD();

  // clean after the VSD data was sucessfuly written
  mTreeK      = 0;
  mTreeH      = 0;
  //mTreeTR     = 0;
  mTreeC      = 0;
  mTreeV0     = 0;
  mTreeKK     = 0;
  mTreeR      = 0;
  mTreeGI     = 0;

  pRunLoader->UnloadAll();
  delete pRunLoader;
  if(gAlice) {
    delete gAlice; gAlice = 0;
  }
  pRunLoader = 0;

  if(mDebugLevel > 0)
    printf("%s all done.\n", eH.Data());
}

void VSDCreator::CreateTrees()
{
  static const Exc_t eH("VSDCreator::CreateTrees ");

  if(mDirectory == 0)
    throw(eH + "output directory not set.");

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertKinematics.\n", eH.Data());
    ConvertKinematics();
  } catch(Exc_t& exc) { WarnCaller(exc); }

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertHits.\n", eH.Data());
    ConvertHits();
  } catch(Exc_t& exc) { WarnCaller(exc); }

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertClusters.\n", eH.Data());
    ConvertClusters();
  } catch(Exc_t& exc) { WarnCaller(exc); }

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertRecTracks.\n", eH.Data());
    ConvertRecTracks();
  } catch(Exc_t& exc) {
    WarnCaller(exc + " Skipping V0 extraction.");
    goto end_esd_processing;
  }

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertV0.\n", eH.Data());
    ConvertV0();
  } catch(Exc_t& exc) { WarnCaller(exc); }

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertKinks.\n", eH.Data());
    ConvertKinks();
  } catch(Exc_t& exc) { WarnCaller(exc); }

end_esd_processing:

  try {
    if(mDebugLevel > 1)
      printf("%s ConvertGenInfo.\n", eH.Data());
    ConvertGenInfo();
  } catch(Exc_t& exc) { WarnCaller(exc); }

  return;
}

/**************************************************************************/
// Kinematics
/**************************************************************************/

void VSDCreator::ConvertKinematics()
{
  static const Exc_t eH("VSDCreator::ConvertKinematics ");

  if(mTreeK != 0) 
    throw (eH + "kinematics already converted");

  AliStack* stack = pRunLoader->Stack();
  if(stack == 0)
    throw(eH + "stack is null.");

  mDirectory->cd();
  mTreeK = new TTree("Kinematics", "TParticles sorted by Label");
 
  Int_t nentries = stack->GetNtrack();
  vector<MCTrack>  vmc(nentries);
  for (Int_t idx=0; idx<nentries; idx++) {
    TParticle* tp = stack->Particle(idx);
    vmc[idx]       = *tp;
    vmc[idx].label = idx;
  }

  // read track refrences 
  TTree* mTreeTR =  pRunLoader->TreeTR();

  if(mTreeTR == 0) {
    WarnCaller(eH + "no TrackRefs; some data will not be available.");
  } else {
    TClonesArray* RunArrayTR = 0;
    mTreeTR->SetBranchAddress("AliRun", &RunArrayTR);

    Int_t nPrimaries = (Int_t) mTreeTR->GetEntries();
    for (Int_t iPrimPart = 0; iPrimPart<nPrimaries; iPrimPart++) {
      // printf("T0 mTreeTR->GetEntry(%d) \n",iPrimPart);
      mTreeTR->GetEntry(iPrimPart);
      // printf("END mTreeTR->GetEntry(%d) \n",iPrimPart);
    
      for (Int_t iTrackRef = 0; iTrackRef < RunArrayTR->GetEntriesFast(); iTrackRef++) {
	AliTrackReference *trackRef = (AliTrackReference*)RunArrayTR->At(iTrackRef); 
	Int_t track = trackRef->GetTrack();
	if(track < nentries && track > 0){ 
	  MCTrack& mct = vmc[track];	
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

  mTreeK->Branch("K", "Reve::MCTrack",  &mpK, fBuffSize);

  printf("sizeofvmc = %d\n", vmc.size());
  for(vector<MCTrack>::iterator k=vmc.begin(); k!=vmc.end(); ++k) {
    MCTrack& mct = *k;
    mK = mct;

    TParticle* m  = &mct;
    Int_t      mi = mct.label;
    int cnt = 0;
    while(m->GetMother(0) != -1) {
      if(cnt > 100) {
	printf("cnt %d mi=%d, mo=%d\n", cnt, mi, m->GetMother(0));
      }
      mi = m->GetMother(0);
      m = &vmc[mi];
      ++cnt;
    }
    mK.eva_label = mi;

    mTreeK->Fill();
  }

  mTreeK->BuildIndex("label");
}

/**************************************************************************/
// Hits
/**************************************************************************/

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

/**************************************************************************/

void VSDCreator::ConvertHits()
{
  static const Exc_t eH("VSDCreator::ConvertHits ");

  if(mTreeH != 0)
    throw(eH + "hits already converted.");

  mDirectory->cd();
  mTreeH =  new TTree("Hits", "Combined detector hits.");
  mTreeH->Branch("H", "Reve::Hit", &mpH, fBuffSize);
 
  map<Int_t, Int_t> hmap;
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
	WarnCaller(eH + "no hits for "+ det.name +".");
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
	  x1=ah->X();y1=ah->Y();z1=ah->Z();
	  if((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1) > tpc_sqr_res) {
	    mH.det_id    = det.detidx;
	    mH.subdet_id = 0;
	    mH.label     = ah->Track();
	    mH.eva_label = eva_idx;
	    mH.V.x = x1; mH.V.y = y1; mH.V.z = z1;
	    mTreeH->Fill();
	    hmap[mH.label]++;
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
	WarnCaller(eH + "no hits for "+ det.name +".");
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
	  mH.det_id    = det.detidx;
	  mH.subdet_id = 0;
	  mH.label     = ali_hit->GetTrack();
	  mH.eva_label = eva_idx;
	  mH.V.Set(ali_hit->X(), ali_hit->Y(), ali_hit->Z());
	  if(det.detidx == 2) {
	    x1=ali_hit->X();y1=ali_hit->Y();z1=ali_hit->Z();
	    if((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1) < trd_sqr_res) continue;
	    x=x1; y=y1; z=z1;
	  } 
	  hmap[mH.label]++;
	  mTreeH->Fill(); 
	}
      }
      delete arr;
      break;
    } // end default 
    } // end switch
  } // end while
  

  //set geninfo
  for(map<Int_t, Int_t>::iterator j=hmap.begin(); j!=hmap.end(); ++j) {
    GetGeninfo(j->first)->n_hits += j->second;
  }
  
}

/**************************************************************************/
// Clusters
/**************************************************************************/

void VSDCreator::ConvertClusters()
{
  static const Exc_t eH("VSDCreator::ConvertClusters ");

  if(mTreeC != 0)
    throw(eH + "clusters already converted.");

  mDirectory->cd();
  mTreeC =  new TTree("Clusters", "rec clusters");
  mTreeC->Branch("C", "Reve::Cluster", &mpC, fBuffSize);

  try {
    ConvertITSClusters();
  } catch(Exc_t& exc) { WarnCaller(exc); }

  try {
    ConvertTPCClusters();
  } catch(Exc_t& exc) { WarnCaller(exc); }
}

/**************************************************************************/

void VSDCreator::ConvertTPCClusters()
{
  static const Exc_t eH("VSDCreator::ConvertTPCClusters ");

  auto_ptr<TFile> f 
    ( TFile::Open(Form("%s/TPC.RecPoints.root", mDataDir.Data())) );
  if(!f.get())
    throw(eH + "can not open 'TPC.RecPoints.root' file.");
    
  auto_ptr<TDirectory> d
    ( (TDirectory*) f->Get(Form("Event%d", mEvent)) );
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
  map<Int_t, Int_t> cmap;

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
          if(cl->GetLabel(0) >= 0){
	    x = par->GetPadRowRadii(sec,row); y = cl->GetY(); z = cl->GetZ();
	    par->AdjustCosSin(sec,cs,sn);
	    tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp; 

	    mC.det_id    = 1;
	    mC.subdet_id = 0;
	    mC.label[0]  = cl->GetLabel(0);
	    mC.label[1]  = cl->GetLabel(1);
	    mC.label[2]  = cl->GetLabel(2);
	    mC.V.Set(x, y, z);

	    mTreeC->Fill();
	    { int i = 0;
	      while(i < 3 && mC.label[i])
		cmap[mC.label[i++]]++;
	    }
	  }
	}
      }
    }
  }
  //set geninfo
  for(map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j) {
    GetGeninfo(j->first)->n_clus += j->second;
  }
}

/**************************************************************************/

void VSDCreator::ConvertITSClusters()
{
  static const Exc_t eH("VSDCreator::ConvertITSClusters ");

  auto_ptr<TFile> f 
    ( TFile::Open(Form("%s/ITS.RecPoints.root", mDataDir.Data())) );
  if(!f.get())
    throw(eH + "can not open 'ITS.RecPoints.root' file.");
    
  auto_ptr<TDirectory> d
    ( (TDirectory*) f->Get(Form("Event%d", mEvent)) );
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
  map<Int_t, Int_t> cmap;

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
      gc[0]= r*cp - recp->GetY()*sp;
      gc[1]= r*sp + recp->GetY()*cp;
      gc[2]= recp->GetZ();

      mC.det_id    = 0;
      mC.subdet_id = 0;
      mC.label[0]  = recp->GetLabel(0);
      mC.label[1]  = recp->GetLabel(1);
      mC.label[2]  = recp->GetLabel(2);
      mC.V.x       = r*cp - recp->GetY()*sp;
      mC.V.y       = r*sp + recp->GetY()*cp;
      mC.V.z       = recp->GetZ();
      mTreeC->Fill();
      { int i = 0;
	while(i < 3 && mC.label[i])
	  cmap[mC.label[i++]]++;
      }
    } 

    for(map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j) {
      GetGeninfo(j->first)->n_clus += j->second;
    }
  }
  delete arr;
}

/**************************************************************************/
// ESD
/**************************************************************************/

void VSDCreator::ConvertRecTracks()
{
  static const Exc_t eH("VSDCreator::ConvertRecTracks ");

  if(mTreeR != 0)
    throw(eH + "tracks already converted.");

  mDirectory->cd();
  mTreeR =  new TTree("RecTracks", "rec tracks");

  mTreeR->Branch("R", "Reve::RecTrack", &mpR, 512*1024,1);
 
  TFile f(Form("%s/AliESDs.root", mDataDir.Data()));
  if(!f.IsOpen())
    throw(eH + "no AliESDs.root file.");

  TTree* tree = (TTree*) f.Get("esdTree");
  if (tree == 0) 
    throw(eH + "no esdTree.");

 
  AliESD *fEvent=0;  
  tree->SetBranchAddress("ESD", &fEvent);
  tree->GetEntry(mEvent); 

 
  // reconstructed tracks
  AliESDtrack* esd_t;
  Double_t     dbuf[3];
  for (Int_t n=0; n<fEvent->GetNumberOfTracks(); n++) {
    esd_t = fEvent->GetTrack(n);

    mR.label  = esd_t->GetLabel();
    mR.status = (Int_t) esd_t->GetStatus();
    mR.sign   = (Int_t) esd_t->GetSign();
    esd_t->GetXYZ(dbuf);    mR.V.Set(dbuf);
    esd_t->GetPxPyPz(dbuf); mR.P.Set(dbuf);
    Double_t ep = esd_t->GetP();
    mR.beta = ep/TMath::Sqrt(ep*ep + TMath::C()*TMath::C()*esd_t->GetMass()*esd_t->GetMass());
    mTreeR->Fill();
  }
  mTreeR->BuildIndex("label");
}

/**************************************************************************/

void VSDCreator::ConvertV0()
{
  static const Exc_t eH("VSDCreator::ConvertV0 ");

  if(mTreeV0 != 0)
    throw(eH + "V0 already converted.");

  mDirectory->cd();
  mTreeV0 =  new TTree("V0", "V0 points");

  mTreeV0->Branch("V0", "Reve::RecV0", &mpV0, 512*1024,1);

  TFile f(Form("%s/AliESDs.root", mDataDir.Data()));
  if(!f.IsOpen()){
    throw(eH + "no AliESDs.root file.");
  }

  TTree* tree = (TTree*) f.Get("esdTree");
  if (tree == 0) 
    throw(eH + "no esdTree.");

  AliESD *fEvent=0;  
  tree->SetBranchAddress("ESD", &fEvent);
  tree->GetEntry(mEvent); 

  for (Int_t n =0; n< fEvent->GetNumberOfV0s(); n++)
  {
    AliESDv0    *av     = fEvent->GetV0(n);
    AliESDtrack *trackN = fEvent->GetTrack(av->GetNindex()); // negative daughter
    AliESDtrack *trackP = fEvent->GetTrack(av->GetPindex()); // positive daughter

    Double_t pos[3];

    mV0.status = av->GetStatus();
    // Point of closest approach
    av->GetXYZ(pos[0],pos[1],pos[2]);
    mV0.V_ca.x = pos[0]; 
    mV0.V_ca.y = pos[1];
    mV0.V_ca.z = pos[2];
    // set birth vertex of neutral particle     
    av->GetXYZ(pos[0], pos[1], pos[2]);
    mV0.V0_birth.Set(pos);

    // momentum and position of negative particle
    av->GetParamN()->GetPxPyPz(pos);
    mV0.P_neg.Set(pos);
    av->GetParamN()->GetXYZ(pos);
    mV0.V_neg.Set(pos);

    // momentum and position of positive particle
    av->GetParamP()->GetPxPyPz(pos);
    mV0.P_pos.Set(pos);
    av->GetParamP()->GetXYZ(pos);
    mV0.V_pos.Set(pos);

    mV0.label = 0; // !!!! mother label unknown
    mV0.pdg   = av->GetPdgCode();

    // daughter indices
    mV0.d_label[0] = TMath::Abs(trackN->GetLabel());
    mV0.d_label[1] = TMath::Abs(trackP->GetLabel());

    // printf("V0 convert labels(%d,%d) index(%d,%d)\n", 
    //	   mV0.d_label[0],  mV0.d_label[1],
    //	   av->GetNIndex(), av->GetPIndex());

    mTreeV0->Fill();
  }
  // if(fEvent->GetNumberOfV0s()) mTreeV0->BuildIndex("label");
}

/**************************************************************************/

void VSDCreator::ConvertKinks()
{
  static const Exc_t eH("VSDCreator::ConvertKinks ");

  if(mTreeKK != 0)
    throw(eH + "Kinks already converted.");

  mDirectory->cd();
  mTreeKK =  new TTree("Kinks", "ESD Kinks");

  mTreeKK->Branch("KK", "Reve::RecKink", &mpKK, fBuffSize);

  TFile f(Form("%s/AliESDs.root", mDataDir.Data()));
  if(!f.IsOpen()){
    throw(eH + "no AliESDs.root file.");
  }

  TTree* tree = (TTree*) f.Get("esdTree");
  if (tree == 0) 
    throw(eH + "no esdTree.");

  AliESD *fEvent=0;  
  tree->SetBranchAddress("ESD", &fEvent);
  tree->GetEntry(mEvent); 

  //  printf("CONVERT KINK Read %d entries in tree kinks \n",  fEvent->GetNumberOfKinks());
  for (Int_t n =0; n< fEvent->GetNumberOfKinks(); n++) {
    AliESDkink* kk = fEvent->GetKink(n);

    Double_t pos[3];


    mKK.label  = kk->GetLabel(0);
    mKK.status = Int_t(kk->GetStatus(1) << 8 + kk->GetStatus(2));

    // reconstructed kink position
    mKK.label_sec = kk->GetLabel(1);
    mKK.V_kink.Set(kk->GetPosition());

    const AliExternalTrackParam& tp_mother = kk->RefParamMother();
    // momentum and position of mother 
    tp_mother.GetPxPyPz(pos);
    mKK.P.Set(pos);
    tp_mother.GetXYZ(pos);
    mKK.V.Set(pos);
    const Double_t* par =  tp_mother.GetParameter();
    // printf("KINK Pt %f, %f \n",1/tp_mother.Pt(),par[4] );
    mKK.sign = (par[4] < 0) ? -1 : 1;
   
    const AliExternalTrackParam& tp_daughter = kk->RefParamDaughter();
    // momentum and position of daughter 
    tp_daughter.GetPxPyPz(pos);
    mKK.P_sec.Set(pos);
    tp_daughter.GetXYZ(pos);
    mKK.V_end.Set(pos);

    mTreeKK->Fill();
  }
  if(fEvent->GetNumberOfKinks()) mTreeKK->BuildIndex("label");
}
/**************************************************************************/
// GenInfo
/**************************************************************************/

void VSDCreator::ConvertGenInfo()
{
  static const Exc_t eH("VSDCreator::ConvertGenInfo ");

  if(mTreeGI != 0)
    throw(eH + "GI already converted.");

  mDirectory->cd();
  mTreeGI = new TTree("GenInfo", "Objects prepared for cross querry");

  GenInfo::Class()->IgnoreTObjectStreamer(true);
  mTreeGI->Branch("GI", "Reve::GenInfo",  &mpGI, fBuffSize);
  mTreeGI->Branch("K.", "Reve::MCTrack",  &mpK);
  mTreeGI->Branch("R.", "Reve::RecTrack", &mpR);

  for(map<Int_t, GenInfo*>::iterator j=mGenInfoMap.begin(); j!=mGenInfoMap.end(); ++j) {
    mGI = *(j->second);
    mGI.label = j->first;
    mTreeK->GetEntry(j->first);

    if(mTreeR) {
      Int_t re = mTreeR->GetEntryNumberWithIndex(j->first);
      if(re != -1) 
	mGI.is_rec = true;
    }
    //    Int_t has_v0 =  mTreeV0->GetEntryNumberWithIndex(j->first);
    //if (has_v0 != -1)
    //  mGI.has_V0 = true;
    if (mTreeKK) {
      Int_t has_kk =  mTreeKK->GetEntryNumberWithIndex(j->first);
      if (has_kk != -1)
	mGI.has_kink = true;
    }
    mTreeGI->Fill();
  }
  mGenInfoMap.clear();
}

/**************************************************************************/
/**************************************************************************/
// Protected methods
/**************************************************************************/
/**************************************************************************/

AliTPCParam* VSDCreator::GetTpcParam(const Exc_t& eh)
{
  auto_ptr<TFile> fp( TFile::Open(Form("%s/galice.root", mDataDir.Data())) );
  if(!fp.get())
    throw(eh + "can not open 'galice.root' file.");
  AliTPCParam* par = (AliTPCParam *) fp->Get("75x40_100x60_150x60");
  if(!par)
    throw(eh + "TPC data not found.");
  return par;
}



GenInfo* VSDCreator::GetGeninfo(Int_t label)
{
  // printf("get_geninfo %d\n", label);
  GenInfo* gi;
  map<Int_t, GenInfo*>::iterator i = mGenInfoMap.find(label);
  if(i == mGenInfoMap.end()) {
    gi =  new GenInfo();
    mGenInfoMap[label] = gi;
  } else {
    gi = i->second;
  }
  return gi;
}
