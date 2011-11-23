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

#include <TDatabasePDG.h>
#include <TSystem.h>
#include <TFile.h>

//______________________________________________________________________________
//
// Create VSD file from ALICE data.

ClassImp(AliEveVSDCreator)

AliEveVSDCreator::AliEveVSDCreator(const Text_t* name, const Text_t* title) :
  TEveVSD(name, title),

  fTPCHitRes  (2),
  fTRDHitRes  (1),

  fDebugLevel (0),
  fRunLoader  (0),
  fGenInfoMap ()
{
  // Constructor.

  // Particles not in ROOT's PDG database occuring in ALICE
  AliPDG::AddParticlesToPdgDataBase();
}

/******************************************************************************/

void AliEveVSDCreator::CreateVSD(const Text_t* vsdFile)
{
  // Create the VSD for current event in AliEveEventManager.
  // Result is stored in directory "Event%04d".

  static const TEveException kEH("AliEveVSDCreator::CreateVSD ");

  AliEveEventManager::AssertGeometry();

  fRunLoader = AliEveEventManager::AssertRunLoader();

  if (fDebugLevel > 0)
    printf("%s open seems ok. Now loading sim data.\n", kEH.Data());

  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();
  fRunLoader->LoadTrackRefs();
  fRunLoader->LoadHits();

  if (fDebugLevel > 0)
    printf("%s opening output TEveVSD.\n", kEH.Data());

  TFile* file = TFile::Open(vsdFile, "UPDATE", "ALICE Visualization Summary Data");

  Int_t curEvent = AliEveEventManager::CurrentEventId();

  TString eventDir; eventDir.Form("Event%04d", curEvent);

  if (file->Get(eventDir))
  {
    Warning(kEH, "Event-dir '%s' existed -- replacing.", eventDir.Data());
    file->Delete(eventDir + ";*");
  }

  fDirectory = new TDirectoryFile(eventDir, Form("Data for event %d", curEvent));

  if (fDebugLevel > 0)
    printf("%s creating trees now ...\n", kEH.Data());

  CreateTrees();

  if (fDebugLevel > 0)
    printf("%s trees created, closing files.\n", kEH.Data());

  file->Write();
  file->Close();
  delete file;
  fDirectory = 0;

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

  if (fDebugLevel > 0)
    printf("%s all done.\n", kEH.Data());
}

void AliEveVSDCreator::CreateTrees()
{
  // Create and fill the output trees by calling all the
  // ConvertXyzz() functions.
  // Exceptions from individual functions are displayed as warnings.

  static const TEveException kEH("AliEveVSDCreator::CreateTrees ");

  if (fDirectory == 0)
    throw kEH + "output directory not set.";

  try
  {
    if (fDebugLevel > 1)
      printf("%sConvertKinematics.\n", kEH.Data());
    ConvertKinematics();
  }
  catch(TEveException& exc)
  {
    Warning(kEH, "%s", exc.Data());
  }

  try
  {
    if (fDebugLevel > 1)
      printf("%sConvertHits.\n", kEH.Data());
    ConvertHits();
  }
  catch(TEveException& exc)
  {
    Warning(kEH, "%s", exc.Data());
  }

  try
  {
    if (fDebugLevel > 1)
      printf("%sConvertClusters.\n", kEH.Data());
    ConvertClusters();
  }
  catch(TEveException& exc)
  {
    Warning(kEH, "%s", exc.Data());
  }

  try 
  {
    if (fDebugLevel > 1)
      printf("%sConvertRecTracks.\n", kEH.Data());
    ConvertRecTracks();
  }
  catch(TEveException& exc)
  {
    Warning(exc, "skipping AliEveV0 extraction.");
    goto end_esd_processing;
  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Warning(kEH, "Explicitly abandoning further conversion.");
  return;

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*
  try
  {
    if (fDebugLevel > 1)
      printf("%sConvertV0.\n", kEH.Data());
    ConvertV0();
  }
  catch(TEveException& exc)
  {
    Warning(kEH, exc);
  }

  try
  {
    if (fDebugLevel > 1)
      printf("%sConvertKinks.\n", kEH.Data());
    ConvertKinks();
  }
  catch(TEveException& exc)
  {
    Warning(kEH, exc);
  }
*/
end_esd_processing:

  try
  {
    if (fDebugLevel > 1)
      printf("%sConvertGenInfo.\n", kEH.Data());
    ConvertGenInfo();
  }
  catch(TEveException& exc)
  {
    Warning(kEH, "%s", exc.Data());
  }

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

  if (fTreeK != 0)
    throw kEH + "kinematics already converted";

  AliStack* stack = fRunLoader->Stack();
  if (stack == 0)
    throw kEH + "stack is null.";

  fDirectory->cd();
  fTreeK = new TTree("Kinematics", "Particles created during simulation.");

  Int_t nentries = stack->GetNtrack();
  std::vector<TEveMCTrack>  vmc(nentries);
  for (Int_t idx = 0; idx < nentries; ++idx)
  {
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
  for(std::vector<TEveMCTrack>::iterator k = vmc.begin(); k != vmc.end(); ++k)
  {
    TEveMCTrack& mct = *k;
    fK = mct;

    TParticle* m  = &mct;
    Int_t      mi = mct.fLabel;
    int cnt = 0;
    while (m->GetMother(0) != -1)
    {
      if (cnt > 100)
	printf("cnt %d mi=%d, mo=%d\n", cnt, mi, m->GetMother(0));

      mi = m->GetMother(0);
      m = &vmc[mi];
      ++cnt;
    }
    fK.fEvaLabel = mi;

    fTreeK->Fill();
  }

  fTreeK->BuildIndex("fLabel");
}


//==============================================================================
// Hits
//==============================================================================

void AliEveVSDCreator::ConvertHits()
{
  // Convert MC hits.
  // TPC hits are handled specially as they are compressed - only mayor
  // hits are stored

  static const TEveException kEH("AliEveVSDCreator::ConvertHits ");

  if (fTreeH != 0)
    throw kEH + "hits already converted.";

  fDirectory->cd();
  fTreeH = new TTree("Hits", "Simulated energy depositions in detector.");
  fTreeH->Branch("H", "TEveHit", &fpH, fBuffSize);

  try
  {
    ConvertAnyHits("ITS", "AliITShit", 0, 0);
    ConvertTPCHits(                    1, fTPCHitRes*fTPCHitRes);
    ConvertAnyHits("TRD", "AliTRDhit", 2, fTRDHitRes*fTRDHitRes);
    ConvertAnyHits("TOF", "AliTOFhit", 3, 0);

  }
  catch(TEveException& exc)
  {
    Warning(kEH, "%s", exc.Data());
  }
}

//------------------------------------------------------------------------------

void AliEveVSDCreator::ConvertAnyHits(const TString& detector,
                                      const TString& clones_class,
                                      Int_t det_id, Float_t minDistSqr)
{
  // Convert hits for detector.

  static const TEveException kEH("AliEveVSDCreator::ConvertAnyHits ");

  Float_t x,y,z, x1,y1,z1;

  TTree* treeh = fRunLoader->GetTreeH(detector, kFALSE);
  if (treeh == 0)
  {
    Warning(kEH, "no hits for %s.", detector.Data());
    return;
  }

  TClonesArray *arr = new TClonesArray(clones_class);
  treeh->SetBranchAddress(detector, &arr);
  Int_t np = treeh->GetEntries();

  // In TreeH files hits are grouped in clones arrays
  // each eva particle has its own clone array.

  std::map<Int_t, Int_t> hmap;

  for (Int_t i = 0; i < np; ++i)
  {
    treeh->GetEntry(i);
    Int_t evaIdx = np - i - 1;
    Int_t nh = arr->GetEntriesFast();
    x = y = z = 0;
    // printf("%d entry %d hits for primary %d \n", i, nh, evaIdx);
    for (Int_t j = 0; j < nh; ++j)
    {
      AliHit* aliHit = (AliHit*)arr->UncheckedAt(j);

      x1 = aliHit->X(); y1 = aliHit->Y(); z1 = aliHit->Z();
      if (minDistSqr == 0 || (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) > minDistSqr)
      {
        fH.fDetId      = det_id;
        fH.fSubdetId   = 0;
        fH.fLabel      = aliHit->GetTrack();
        fH.fEvaLabel   = evaIdx;
        fH.fV.Set(aliHit->X(), aliHit->Y(), aliHit->Z());

        fTreeH->Fill();

        hmap[fH.fLabel]++;

        x=x1; y=y1; z=z1;
      }
    }
  }
  delete arr;

  // Set geninfo.
  for (std::map<Int_t, Int_t>::iterator j = hmap.begin(); j != hmap.end(); ++j)
  {
    GetGeninfo(j->first)->fNHits += j->second;
  }
}

//------------------------------------------------------------------------------

void AliEveVSDCreator::ConvertTPCHits(Int_t det_id, Float_t minDistSqr)
{
  // Convert hits for TPC.

  static const TEveException kEH("AliEveVSDCreator::ConvertTPCHits ");

  Float_t x,y,z, x1,y1,z1;

  TTree* treeh = fRunLoader->GetTreeH("TPC", false);
  if (treeh == 0)
  {
    Warning(kEH, "no hits for %s.", "TPC");
    return;
  }

  AliTPCTrackHitsV2 hv2, *hv2p = &hv2;
  treeh->SetBranchAddress("TPC2", &hv2p);
  Int_t np = treeh->GetEntries();

  std::map<Int_t, Int_t> hmap;

  for (Int_t i = 0; i < np; ++i)
  {
    treeh->GetEntry(i);
    Int_t evaIdx = np -i -1;
    if (hv2.First() == 0) continue;
    x = y = z = 0;
    do
    {
      AliHit* ah = hv2.GetHit();
      x1 = ah->X(); y1 = ah->Y(); z1 = ah->Z();
      if (minDistSqr == 0 || (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) > minDistSqr)
      {
        fH.fDetId    = det_id;
        fH.fSubdetId = 0;
        fH.fLabel    = ah->Track();
        fH.fEvaLabel = evaIdx;
        fH.fV.fX = x1; fH.fV.fY = y1; fH.fV.fZ = z1;

        fTreeH->Fill();

        hmap[fH.fLabel]++;

        x = x1; y = y1; z = z1;
      }
    } while (hv2.Next());
  }

  // Set geninfo.
  for (std::map<Int_t, Int_t>::iterator j = hmap.begin(); j != hmap.end(); ++j)
  {
    GetGeninfo(j->first)->fNHits += j->second;
  }
}


//==============================================================================
// Clusters
//==============================================================================

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

  if (fTreeC != 0)
    throw kEH + "clusters already converted.";

  fDirectory->cd();
  fTreeC = new TTree("Clusters", "Reconstructed points of particle passage.");
  fTreeC->Branch("C", "TEveCluster", &fpC, fBuffSize);

  try
  {
    ConvertAnyClusters("ITS", "ITSRecPoints", 0);
    ConvertTPCClusters(                       1);
    ConvertAnyClusters("TRD", "TRDcluster",   2);
    ConvertAnyClusters("TOF", "TOF",          3);

  }
  catch(TEveException& exc)
  {
    Warning(kEH, "%s", exc.Data());
  }
}

//------------------------------------------------------------------------------

void AliEveVSDCreator::ConvertAnyClusters(const TString& detector,
                                          const TString& branch_name,
                                          Int_t det_id)
{
  // Convert clusters for detector and transform them to global coordinates.

  static const TEveException kEH("AliEveVSDCreator::ConvertITSClusters ");

  fRunLoader->LoadRecPoints(detector);
  TTree* tree = fRunLoader->GetTreeR(detector, false);
  if (!tree)
    throw kEH + "'TreeR' not found.";

  TClonesArray *cl = 0;
  TBranch *branch  = tree->GetBranch(branch_name);
  branch->SetAddress(&cl);

  Int_t nmods = branch->GetEntries();

  std::map<Int_t, Int_t> cmap;

  for (Int_t mod = 0; mod < nmods; ++mod)
  {
    branch->GetEntry(mod);
    Int_t nc = cl->GetEntriesFast();
    for (Int_t j = 0; j < nc; ++j)
    {
      AliCluster *c = (AliCluster*) cl->UncheckedAt(j);

      fC.fDetId    = det_id;
      fC.fSubdetId = mod;
      fC.fLabel[0] = c->GetLabel(0);
      fC.fLabel[1] = c->GetLabel(1);
      fC.fLabel[2] = c->GetLabel(2);

      c->GetGlobalXYZ(fC.fV.Arr());

      fTreeC->Fill();

      {
        int i = 0;
	while (i < 3 && fC.fLabel[i])
	  cmap[fC.fLabel[i++]]++;
      }
    }
  }
  delete cl;

  for (std::map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j)
  {
    GetGeninfo(j->first)->fNClus += j->second;
  }
}

//------------------------------------------------------------------------------

void AliEveVSDCreator::ConvertTPCClusters(Int_t det_id)
{
  // Convert TPC clusters and transform them to global coordinates.

  static const TEveException kEH("AliEveVSDCreator::ConvertTPCClusters ");

  const Int_t kMaxCl = 100*160;

  fRunLoader->LoadRecPoints("TPC");
  TTree* tree = fRunLoader->GetTreeR("TPC", false);
  if (!tree)
    throw kEH + "'TreeR' not found.";

  AliTPCClustersRow *clrow = new AliTPCClustersRow;
  clrow->SetClass("AliTPCclusterMI");
  clrow->SetArray(kMaxCl);
  tree->SetBranchAddress("Segment", &clrow);

  Int_t nEnt = tree->GetEntries();

  std::map<Int_t, Int_t> cmap;

  for (Int_t j = 0; j < nEnt; ++j)
  {
    if (!tree->GetEvent(j)) continue;

    TClonesArray *cl = clrow->GetArray();
    Int_t ncl = cl->GetEntriesFast();

    while (ncl--)
    {
      AliCluster *c = (AliCluster*) cl->UncheckedAt(ncl);

      fC.fDetId    = det_id;
      fC.fSubdetId = j;
      fC.fLabel[0] = c->GetLabel(0);
      fC.fLabel[1] = c->GetLabel(1);
      fC.fLabel[2] = c->GetLabel(2);

      c->GetGlobalXYZ(fC.fV.Arr());

      fTreeC->Fill();

      {
        int i = 0;
	while (i < 3 && fC.fLabel[i])
	  cmap[fC.fLabel[i++]]++;
      }
    }

    cl->Clear();
  }
  delete clrow;

  //set geninfo
  for (std::map<Int_t, Int_t>::iterator j=cmap.begin(); j!=cmap.end(); ++j)
  {
    GetGeninfo(j->first)->fNClus += j->second;
  }
}


/******************************************************************************/
// ESD
/******************************************************************************/

void AliEveVSDCreator::ConvertRecTracks()
{
  // Convert reconstructed tracks.

  static const TEveException kEH("AliEveVSDCreator::ConvertRecTracks ");

  if (fTreeR != 0)
    throw kEH + "tracks already converted.";

  AliESDEvent* esdEvent = AliEveEventManager::AssertESD();

  fDirectory->cd();
  fTreeR =  new TTree("RecTracks", "Reconstructed particle trajectories.");
  fTreeR->Branch("R", "TEveRecTrack", &fpR, 512*1024);

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

  fTreeR->BuildIndex("fLabel");
}

/******************************************************************************/

void AliEveVSDCreator::ConvertV0()
{
  // Convert reconstructed V0s.

  static const TEveException kEH("AliEveVSDCreator::ConvertV0 ");

  if (fTreeV0 != 0)
    throw kEH + "AliEveV0 already converted.";

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
    throw kEH + "Kinks already converted.";

  throw kEH + "Currently non-supported - TEveRecKink being updated.";

  /*
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
  */
}
/******************************************************************************/
// TEveMCRecCrossRef
/******************************************************************************/

void AliEveVSDCreator::ConvertGenInfo()
{
  // Build simulation-reconstruction cross-reference table.
  // In a rather poor state at the moment.

  static const TEveException kEH("AliEveVSDCreator::ConvertGenInfo ");

  if (fTreeGI != 0)
    throw kEH + "GI already converted.";

  fDirectory->cd();
  fTreeGI = new TTree("TEveMCRecCrossRef", "Objects prepared for cross querry");

  TEveMCRecCrossRef::Class()->IgnoreTObjectStreamer(true);
  fTreeGI->Branch("GI", "TEveMCRecCrossRef",  &fpGI, fBuffSize);
  fTreeGI->Branch("K.", "TEveMCTrack",  &fpK);
  fTreeGI->Branch("R.", "TEveRecTrack", &fpR);

  for (std::map<Int_t, TEveMCRecCrossRef*>::iterator j=fGenInfoMap.begin(); j!=fGenInfoMap.end(); ++j)
  {
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
    if (fTreeKK)
    {
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
  if (i == fGenInfoMap.end())
  {
    gi =  new TEveMCRecCrossRef();
    fGenInfoMap[label] = gi;
  }
  else
  {
    gi = i->second;
  }
  return gi;
}
