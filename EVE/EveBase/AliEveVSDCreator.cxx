// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <memory>
#include "AliEveVSDCreator.h"

#include "AliEveEventManager.h"

#include <AliStack.h>
#include <AliITSLoader.h>
#include <AliTPCTrackHitsV2.h>
#include <AliPDG.h>
#include <AliHit.h>
#include <AliESDEvent.h>
#include <AliESDv0.h>
#include <AliTPCclusterMI.h>
#include <AliTPCClustersRow.h>
#include <AliITSclusterV2.h>
#include <AliESDkink.h>
#include <AliESDtrack.h>

#include <AliRun.h>
#include <AliRunLoader.h>
#include <AliTPCParam.h>

#include <TSystem.h>
#include <TFile.h>

//______________________________________________________________________________
//
// Create VSD file from ALICE data.

ClassImp(AliEveVSDCreator)

AliEveVSDCreator::AliEveVSDCreator(const Text_t* name, const Text_t* title) :
  TEveVSD(name, title),

  fTPCHitRes  (2),
  fTRDHitRes  (2),

  fDebugLevel (0),
  fRunLoader  (0),
  fGenInfoMap ()
{
  // Constructor.

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

void AliEveVSDCreator::CreateVSD(const Text_t* vsdFile)
{
  // Create the VSD for current event in AliEveEventManager.
  // Result is stored in vsdFile.
  //
  // Needs to be extended to support conversion of multiple events.

  static const TEveException kEH("AliEveVSDCreator::CreateVSD ");

  fRunLoader = AliEveEventManager::AssertRunLoader();

  if(fDebugLevel > 0)
    printf("%s open seems ok. Now loading sim data.\n", kEH.Data());

  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();
  fRunLoader->LoadTrackRefs();
  fRunLoader->LoadHits();

  // GledNS::PushFD();

  if(fDebugLevel > 0)
    printf("%s opening output TEveVSD.\n", kEH.Data());

  TFile* file = TFile::Open(vsdFile, "RECREATE", "ALICE Visualization Summary Data");
  fDirectory = new TDirectoryFile("Event0", "");

  if(fDebugLevel > 0)
    printf("%s creating trees now ...\n", kEH.Data());

  CreateTrees();

  if(fDebugLevel > 0)
    printf("%s trees created, closing files.\n", kEH.Data());

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

  fRunLoader = 0;

  if(fDebugLevel > 0)
    printf("%s all done.\n", kEH.Data());
}

void AliEveVSDCreator::CreateTrees()
{
  // Create and fill the output trees by calling all the
  // ConvertXyzz() functions.
  // Exceptions from individual functions are displayed as warnings.

  static const TEveException kEH("AliEveVSDCreator::CreateTrees ");

  if (fDirectory == 0)
    throw(kEH + "output directory not set.");

  try {
    if (fDebugLevel > 1)
      printf("%sConvertKinematics.\n", kEH.Data());
    ConvertKinematics();
  } catch(TEveException& exc) { Warning(kEH, exc); }

  Warning(kEH, "Explicitly abandoning further conversion.");
  return;

  try {
    if (fDebugLevel > 1)
      printf("%sConvertHits.\n", kEH.Data());
    ConvertHits();
  } catch(TEveException& exc) { Warning(kEH, exc); }

  try {
    if (fDebugLevel > 1)
      printf("%sConvertClusters.\n", kEH.Data());
    ConvertClusters();
  } catch(TEveException& exc) { Warning(kEH, exc); }

  try {
    if (fDebugLevel > 1)
      printf("%sConvertRecTracks.\n", kEH.Data());
    ConvertRecTracks();
  } catch(TEveException& exc) {
    Warning(exc, "skipping AliEveV0 extraction.");
    goto end_esd_processing;
  }

  try {
    if (fDebugLevel > 1)
      printf("%sConvertV0.\n", kEH.Data());
    ConvertV0();
  } catch(TEveException& exc) { Warning(kEH, exc); }

  try {
    if (fDebugLevel > 1)
      printf("%sConvertKinks.\n", kEH.Data());
    ConvertKinks();
  } catch(TEveException& exc) { Warning(kEH, exc); }

end_esd_processing:

  try {
    if (fDebugLevel > 1)
      printf("%sConvertGenInfo.\n", kEH.Data());
    ConvertGenInfo();
  } catch(TEveException& exc) { Warning(kEH, exc); }

  return;
}

/******************************************************************************/
// Kinematics
/******************************************************************************/

void AliEveVSDCreator::ConvertKinematics()
{
  // Convert kinematics.
  // Track references are not stored, they should be.

  static const TEveException kEH("AliEveVSDCreator::ConvertKinematics ");

  if(fTreeK != 0)
    throw (kEH + "kinematics already converted");

  AliStack* stack = fRunLoader->Stack();
  if(stack == 0)
    throw(kEH + "stack is null.");

  fDirectory->cd();
  fTreeK = new TTree("Kinematics", "TParticles sorted by Label");

  Int_t nentries = stack->GetNtrack();
  std::vector<TEveMCTrack>  vmc(nentries);
  for (Int_t idx=0; idx<nentries; idx++) {
    TParticle*   tp = stack->Particle(idx);
    vmc[idx]        = *tp;
    vmc[idx].fLabel = idx;
  }

  // read track refrences
  // functionality now in AliEveKineTools.
  /*
  TTree* fTreeTR =  fRunLoader->TreeTR();

  if(fTreeTR == 0) {
    Warning(kEH, "no TrackRefs; some data will not be available.");
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

  printf("sizeofvmc = %d\n", (Int_t) vmc.size());
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

  fTreeK->BuildIndex("fLabel");
}

/******************************************************************************/
// Hits
/******************************************************************************/

namespace {

  struct Detector_t
  {
    const char*   fName;      // Detector name.
    const char*   fHitbranch; // Name of branch containing hits.
    unsigned char fDetidx;    // Index identifying the detector internally.
  };

  Detector_t fgDetectors[] = {
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
  // Convert MC hits.
  // TPC hits are handled specially as they are compressed - only mayor
  // hits are stored

  static const TEveException kEH("AliEveVSDCreator::ConvertHits ");

  if (fTreeH != 0)
    throw(kEH + "hits already converted.");

  fDirectory->cd();
  fTreeH =  new TTree("Hits", "Combined detector hits.");
  fTreeH->Branch("H", "TEveHit", &fpH, fBuffSize);

  std::map<Int_t, Int_t> hmap;
  // parameters for ITS, TPC hits filtering
  Float_t x,y,z, x1,y1,z1;
  Float_t tpcSqrRes = fTPCHitRes*fTPCHitRes;
  Float_t trdSqrRes = fTRDHitRes*fTRDHitRes;

  int l=0;
  // load hits from the rest of detectors
  while (fgDetectors[l].fName != 0)
  {
    Detector_t& det = fgDetectors[l++];

    switch(det.fDetidx)
    {
      case 1:
      {
	Int_t count = 0;
	TTree* treeh = fRunLoader->GetTreeH(det.fName, false);
	if(treeh == 0) {
	  Warning(kEH, Form("no hits for %s.", det.fName));
	  continue;
	}
	AliTPCTrackHitsV2 hv2, *hv2p = &hv2;
	treeh->SetBranchAddress("TPC2", &hv2p);
	Int_t np = treeh->GetEntries();
	for (Int_t i = 0; i < np; ++i)
	{
	  treeh->GetEntry(i);
	  Int_t evaIdx = np -i -1;
	  if (hv2.First() == 0) continue;
	  x = y = z = 0;
	  do {
	    AliHit* ah = hv2.GetHit();
	    x1 = ah->X(); y1 = ah->Y(); z1 = ah->Z();
	    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) > tpcSqrRes)
	    {
	      fH.fDetId    = det.fDetidx;
	      fH.fSubdetId = 0;
	      fH.fLabel    = ah->Track();
	      fH.fEvaLabel = evaIdx;
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
      default:
      {
	TTree* treeh = fRunLoader->GetTreeH(det.fName, false);
	if (treeh == 0) {
	  Warning(kEH, Form("no hits for %s.", det.fName));
	  continue;
	}
	TClonesArray *arr = new TClonesArray(det.fHitbranch);
	treeh->SetBranchAddress(det.fName, &arr);
	Int_t np = treeh->GetEntries();
	// in TreeH files hits are grouped in clones arrays
	// each eva particle has its own clone array
	for (Int_t i = 0; i < np; ++i)
	{
	  treeh->GetEntry(i);
	  Int_t evaIdx = np -i -1;
	  Int_t nh=arr->GetEntriesFast();
	  x = y = z = 0;
	  // printf("%d entry %d hits for primary %d \n", i, nh, evaIdx);
	  for (Int_t j = 0; j < nh; ++j)
	  {
	    AliHit* aliHit = (AliHit*)arr->UncheckedAt(j);
	    fH.fDetId    = det.fDetidx;
	    fH.fSubdetId = 0;
	    fH.fLabel     = aliHit->GetTrack();
	    fH.fEvaLabel = evaIdx;
	    fH.fV.Set(aliHit->X(), aliHit->Y(), aliHit->Z());
	    if (det.fDetidx == 2)
	    {
	      x1=aliHit->X();y1=aliHit->Y();z1=aliHit->Z();
	      if ((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1) < trdSqrRes) continue;
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
  for(std::map<Int_t, Int_t>::iterator j=hmap.begin(); j!=hmap.end(); ++j)
  {
    GetGeninfo(j->first)->fNHits += j->second;
  }
}

/******************************************************************************/
// Clusters
/******************************************************************************/

void AliEveVSDCreator::ConvertClusters()
{
  // Convert clusters.
  //
  // Only supported for ITS and TPC at the moment, see dedicated
  // functions ConvertITSClusters() and ConvertTPCClusters().
  //
  // It should be possible now to do this in a general manner (with
  // the alignment framework).

  static const TEveException kEH("AliEveVSDCreator::ConvertClusters ");

  if(fTreeC != 0)
    throw(kEH + "clusters already converted.");

  fDirectory->cd();
  fTreeC =  new TTree("Clusters", "rec clusters");
  fTreeC->Branch("C", "TEveCluster", &fpC, fBuffSize);

  try {
    ConvertITSClusters();
  } catch(TEveException& exc) { Warning(kEH, exc); }

  try {
    ConvertTPCClusters();
  } catch(TEveException& exc) { Warning(kEH, exc); }
}

/******************************************************************************/

void AliEveVSDCreator::ConvertTPCClusters()
{
  // Convert TPC clusters and transform them to global coordinates.

  static const TEveException kEH("AliEveVSDCreator::ConvertTPCClusters ");

  fRunLoader->LoadRecPoints("TPC");
  TTree* tree = fRunLoader->GetTreeR("TPC", false);
  if (!tree)
    throw(kEH + "'TreeR' not found.");

  AliTPCClustersRow  clrow, *clrowp = &clrow;
  AliTPCclusterMI   *cl;
  clrow.SetClass("AliTPCclusterMI");
  tree->SetBranchAddress("Segment", &clrowp);

  // count clusters
  Int_t nClusters = 0;
  Int_t nEnt = tree->GetEntries();
  for (Int_t n = 0; n < nEnt; ++n)
  {
    tree->GetEntry(n);
    nClusters += clrow.GetArray()->GetEntriesFast();
  }

  std::map<Int_t, Int_t> cmap;

  for (Int_t n = 0; n < tree->GetEntries(); ++n)
  {
    tree->GetEntry(n);
    Int_t ncl = clrow.GetArray()->GetEntriesFast();
    if (ncl > 0)
    {
      while (ncl--)
      {
	if (clrow.GetArray())
	{
	  // cl = new AliTPCclusterMI(*(AliTPCclusterMI*)clrow.GetArray()->UncheckedAt(ncl));
	  cl = (AliTPCclusterMI*)clrow.GetArray()->UncheckedAt(ncl);
          if (cl->GetLabel(0) >= 0)
	  {
            Float_t g[3]; //global coordinates
            cl->GetGlobalXYZ(g);

	    fC.fDetId    = 1;
	    fC.fSubdetId = clrow.GetID();
	    fC.fLabel[0] = cl->GetLabel(0);
	    fC.fLabel[1] = cl->GetLabel(1);
	    fC.fLabel[2] = cl->GetLabel(2);
	    fC.fV.Set(g);

	    fTreeC->Fill();
	    {
	      int i = 0;
	      while (i < 3 && fC.fLabel[i])
		cmap[fC.fLabel[i++]]++;
	    }
	  }
	}
      }
    }
  }
  //set geninfo
  for (std::map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j)
  {
    GetGeninfo(j->first)->fNClus += j->second;
  }
}

/******************************************************************************/

void AliEveVSDCreator::ConvertITSClusters()
{
  // Convert ITS clusters and transform them to global coordinates.

  static const TEveException kEH("AliEveVSDCreator::ConvertITSClusters ");

  fRunLoader->LoadRecPoints("ITS");
  TTree* tree = fRunLoader->GetTreeR("ITS", false);
  if (!tree)
    throw(kEH + "'TreeR' not found.");

  // 
  AliITSLoader *itsLd = (AliITSLoader*) fRunLoader->GetLoader("ITSLoader");
  AliITSgeom   *geom  = itsLd->GetITSgeom();

  //printf("alice ITS geom %p \n",geom );

  if (!geom)
    throw(kEH + "can not find ITS geometry");

  TClonesArray *arr = new TClonesArray("AliITSclusterV2");
  tree->SetBranchAddress("Clusters", &arr);
  Int_t nmods = tree->GetEntries();
  Float_t gc[3];
  std::map<Int_t, Int_t> cmap;

  for (Int_t mod = 0; mod < nmods; ++mod)
  {
    tree->GetEntry(mod);
    Int_t nc=arr->GetEntriesFast();
    for (Int_t j = 0; j < nc; ++j)
    {
      AliITSclusterV2* recp = (AliITSclusterV2*) arr->UncheckedAt(j);

      Double_t rot[9];
      geom->GetRotMatrix(mod,rot);
      Int_t lay,lad,det;
      geom->GetModuleId(mod,lay,lad,det);
      Float_t tx,ty,tz;
      geom->GetTrans(lay,lad,det,tx,ty,tz);

      Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
      Double_t phi1=TMath::Pi()/2+alpha;
      if (lay == 1) phi1+=TMath::Pi();

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

    for (std::map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j)
    {
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
  // Convert reconstructed tracks.

  static const TEveException kEH("AliEveVSDCreator::ConvertRecTracks ");

  if (fTreeR != 0)
    throw(kEH + "tracks already converted.");

  AliESDEvent* esdEvent = AliEveEventManager::AssertESD();

  fDirectory->cd();
  fTreeR =  new TTree("RecTracks", "rec tracks");

  fTreeR->Branch("R", "TEveRecTrack", &fpR, 512*1024,1);

  // reconstructed tracks
  AliESDtrack* esdTrack;
  Double_t     dbuf[3];
  for (Int_t n = 0; n < esdEvent->GetNumberOfTracks(); ++n)
  {
    esdTrack = esdEvent->GetTrack(n);

    fR.fLabel  = esdTrack->GetLabel();
    fR.fStatus = (Int_t) esdTrack->GetStatus();
    fR.fSign   = (Int_t) esdTrack->GetSign();
    esdTrack->GetXYZ(dbuf);    fR.fV.Set(dbuf);
    esdTrack->GetPxPyPz(dbuf); fR.fP.Set(dbuf);
    Double_t ep = esdTrack->GetP();
    fR.fBeta = ep/TMath::Sqrt(ep*ep + TMath::C()*TMath::C()*esdTrack->GetMass()*esdTrack->GetMass());
    fTreeR->Fill();
  }
  fTreeR->BuildIndex("label");
  delete esdEvent;
}

/******************************************************************************/

void AliEveVSDCreator::ConvertV0()
{
  // Convert reconstructed V0s.

  static const TEveException kEH("AliEveVSDCreator::ConvertV0 ");

  if(fTreeV0 != 0)
    throw(kEH + "AliEveV0 already converted.");

  AliESDEvent* esdEvent = AliEveEventManager::AssertESD();

  fDirectory->cd();
  fTreeV0 =  new TTree("AliEveV0", "AliEveV0 points");

  fTreeV0->Branch("AliEveV0", "TEveRecV0", &fpV0, 512*1024,1);

  for (Int_t n = 0; n < esdEvent->GetNumberOfV0s(); ++n)
  {
    AliESDv0    *av     = esdEvent->GetV0(n);
    AliESDtrack *trackN = esdEvent->GetTrack(av->GetNindex()); // negative daughter
    AliESDtrack *trackP = esdEvent->GetTrack(av->GetPindex()); // positive daughter

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
  // if (esdEvent->GetNumberOfV0s()) fTreeV0->BuildIndex("label");
  delete esdEvent;
}

/******************************************************************************/

void AliEveVSDCreator::ConvertKinks()
{
  // Convert reconstructed kinks.

  static const TEveException kEH("AliEveVSDCreator::ConvertKinks ");

  if (fTreeKK != 0)
    throw(kEH + "Kinks already converted.");

  AliESDEvent* esdEvent = AliEveEventManager::AssertESD();

  fDirectory->cd();
  fTreeKK =  new TTree("Kinks", "ESD Kinks");

  fTreeKK->Branch("KK", "TEveRecKink", &fpKK, fBuffSize);

  //  printf("CONVERT KINK Read %d entries in tree kinks \n",  esdEvent->GetNumberOfKinks());
  for (Int_t n = 0; n < esdEvent->GetNumberOfKinks(); ++n)
  {
    AliESDkink* kk = esdEvent->GetKink(n);

    Double_t pos[3];

    fKK.fLabel  = kk->GetLabel(0);
    fKK.fStatus = 0; // status is Char_t[12] ... have no idea how/what to extract.

    // reconstructed kink position
    fKK.fLabelSec = kk->GetLabel(1);
    fKK.fVKink.Set(kk->GetPosition());

    const AliExternalTrackParam& tpMother = kk->RefParamMother();
    // momentum and position of mother
    tpMother.GetPxPyPz(pos);
    fKK.fP.Set(pos);
    tpMother.GetXYZ(pos);
    fKK.fV.Set(pos);
    const Double_t* par =  tpMother.GetParameter();
    // printf("KINK Pt %f, %f \n",1/tpMother.Pt(),par[4] );
    fKK.fSign = (par[4] < 0) ? -1 : 1;

    const AliExternalTrackParam& tpDaughter = kk->RefParamDaughter();
    // momentum and position of daughter
    tpDaughter.GetPxPyPz(pos);
    fKK.fPSec.Set(pos);
    tpDaughter.GetXYZ(pos);
    fKK.fVEnd.Set(pos);

    fTreeKK->Fill();
  }
  if (esdEvent->GetNumberOfKinks()) fTreeKK->BuildIndex("label");
  delete esdEvent;
}
/******************************************************************************/
// TEveMCRecCrossRef
/******************************************************************************/

void AliEveVSDCreator::ConvertGenInfo()
{
  // Build simulation-reconstruction cross-reference table.
  // In a rather poor state at the moment.

  static const TEveException kEH("AliEveVSDCreator::ConvertGenInfo ");

  if(fTreeGI != 0)
    throw(kEH + "GI already converted.");

  fDirectory->cd();
  fTreeGI = new TTree("TEveMCRecCrossRef", "Objects prepared for cross querry");

  TEveMCRecCrossRef::Class()->IgnoreTObjectStreamer(true);
  fTreeGI->Branch("GI", "TEveMCRecCrossRef",  &fpGI, fBuffSize);
  fTreeGI->Branch("K.", "TEveMCTrack",  &fpK);
  fTreeGI->Branch("R.", "TEveRecTrack", &fpR);

  for (std::map<Int_t, TEveMCRecCrossRef*>::iterator j=fGenInfoMap.begin(); j!=fGenInfoMap.end(); ++j) {
    fGI        = *(j->second);
    fGI.fLabel = j->first;
    fTreeK->GetEntry(j->first);

    if (fTreeR) {
      Int_t re = fTreeR->GetEntryNumberWithIndex(j->first);
      if(re != -1)
	fGI.fIsRec = true;
    }
    //    Int_t hasV0 =  fTreeV0->GetEntryNumberWithIndex(j->first);
    //if (hasV0 != -1)
    //  fGI.has_AliEveV0 = true;
    if (fTreeKK) {
      Int_t hasKk =  fTreeKK->GetEntryNumberWithIndex(j->first);
      if (hasKk != -1)
	fGI.fHasKink = true;
    }
    fTreeGI->Fill();
  }
  fGenInfoMap.clear();
}

/******************************************************************************/
/******************************************************************************/
// Protected methods
/******************************************************************************/
/******************************************************************************/

TEveMCRecCrossRef* AliEveVSDCreator::GetGeninfo(Int_t label)
{
  // Return the cross-reference structure for given label.
  // If the object does not exist it is created.

  // printf("get_geninfo %d\n", label);
  TEveMCRecCrossRef* gi;
  std::map<Int_t, TEveMCRecCrossRef*>::iterator i = fGenInfoMap.find(label);
  if (i == fGenInfoMap.end()) {
    gi =  new TEveMCRecCrossRef();
    fGenInfoMap[label] = gi;
  } else {
    gi = i->second;
  }
  return gi;
}
