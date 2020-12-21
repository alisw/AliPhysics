/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTRDefficiency.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TStyle.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TH3.h>
#include <THStack.h>
#include "TTreeStream.h"
#include "TPDGCode.h"

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliExternalTrackParam.h"
#include "AliTracker.h"
#include "AliAnalysisManager.h"

#include "AliTRDgeometry.h"
#include "AliTRDtrackV1.h"
#include "AliTRDCalPID.h"
#include "AliTRDefficiency.h"
#include "AliTRDtrackInfo.h"

ClassImp(AliTRDefficiency)

//____________________________________________________________________
AliTRDefficiency::AliTRDefficiency()
  :AliTRDrecoTask()
  ,fMissed(NULL)
  ,fProj(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDefficiency", "TRD barrel tracking efficiency checker");
}

//____________________________________________________________________
AliTRDefficiency::AliTRDefficiency(char* name)
  :AliTRDrecoTask(name, "TRD barrel tracking efficiency checker")
  ,fMissed(NULL)
  ,fProj(NULL)
{
  //
  // Default constructor
  //
}

//____________________________________________________________________
AliTRDefficiency::~AliTRDefficiency()
{
  // Destructor
  if(fMissed){
    fMissed->Delete();
    delete fMissed;
  }
}

// //____________________________________________________________________
// void  AliTRDefficiency::UserCreateOutputObjects()
// {
//   //
//   // Create output objects
//   //
// 
//   const Int_t nbins = AliTRDCalPID::kNMom;
//   Float_t xbins[nbins+1] = {.5, .7, .9, 1.3, 1.7, 2.4, 3.5, 4.5, 5.5, 7., 9., 11.};
// 
//   TH1 *h = NULL;
//   fContainer = new TObjArray(); fContainer->SetOwner();
//   for(Int_t is=0; is<AliPID::kSPECIES; is++){
//     fContainer->Add(h = new TProfile(Form("h%s", AliTRDCalPID::GetPartSymb(is)), AliPID::ParticleShortName(is), nbins, xbins));
//     h->SetLineColor(AliTRDCalPID::GetPartColor(is));
//     h->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
//     h->SetMarkerStyle(24);
//   }
//   fContainer->Add(h = new TProfile("h", "", nbins, xbins));
//   h->SetMarkerStyle(7);
//   PostData(1, fContainer);
// }

//____________________________________________________________________
TH1* AliTRDefficiency::PlotBasicEff(const AliTRDtrackV1 *track)
{
// plot TRD efficiency based on ESD info

  if(!fkESD){
    AliDebug(4, "No ESD info.");
    return NULL;
  }

  THnSparse *H(NULL);
  if(!fContainer || !(H = ((THnSparse*)fContainer->FindObject("hEFF")))){
    AliWarning("No output container defined.");
    return NULL;
  }
  if(track) fkTrack = track;

  Double_t val[8]; memset(val, 0, 8*sizeof(Double_t));
  ULong_t status(fkESD->GetStatus());
//     printf("  %3d ITS[%c] TPC[%c] TRDin[%c] TRDout[%c] TRDStop[%c]\n", fkESD->GetId(),
//       (status&AliESDtrack::kITSout)?'y':'n',
//       (status&AliESDtrack::kTPCout)?'y':'n',
//       (status&AliESDtrack::kTRDin)?'y':'n',
//       (status&AliESDtrack::kTRDout)?'y':'n',
//       (status&AliESDtrack::kTRDStop)?'y':'n');


  val[0] =((status&AliESDtrack::kTRDin)?1:0) +
          ((status&AliESDtrack::kTRDStop)?1:0) +
          ((status&AliESDtrack::kTRDout)?2:0);
  val[1] = fkESD->Phi();
  val[2] = fkESD->Eta();
  val[3] = GetPtBin(fPt>0.?fPt:fkESD->Pt());
  val[4] = -2.; val[5] = fkMC?fkMC->GetPt():-1.;
  if(fkMC && fkMC->GetNTrackRefs()>=2){ // TODO define trackable tracks more exactly
    val[4] = -1.;
    if(fkTrack && fkTrack->GetNumberOfTracklets()>0){ // TODO better define 
      if(fkMC->GetLabel() == fkMC->GetTRDlabel()) val[4] = 0.;        // good tracks
      else if(fkMC->GetLabel() == -fkMC->GetTRDlabel()) val[4] = 1.;  // acceptable tracks
      else val[4] = 2.;                                               // wrongly matched tracks
    }
    val[5] = GetPtBin(fkMC->GetTrackRef()->Pt());
  }

  val[6] = fkTrack?fkTrack->GetNumberOfTracklets():0;
  // down scale PID resolution
  Int_t spc(fSpecies); if(spc==-6) spc=0; if(spc==3) spc=2; if(spc==-3) spc=-2; if(spc>3) spc=3; if(spc<-3) spc=-3;
  val[7] = spc;
  if(fkMC){
    Int_t exactPID = 0;
    switch(fkMC->GetPDG()){
    case kElectron: exactPID = -1;break;
    case kPositron: exactPID = 1;break;
    case kPiPlus: exactPID = 2;break;
    case kPiMinus: exactPID = -2;break;
    case kProton: exactPID = 3;break;
    case kProtonBar: exactPID = -3;break;
    }
    val[7] = exactPID;
  }
//     printf("%3d label[%2d %+2d] PDG[%+4d] TR[%d] Trklts[%d] species[%d %d] pt%d[%6.4f] MC_pt%d[%6.4f] in[%c] out[%c] stop[%c] flag[%d]\n",
//     fkESD->GetId(), fkMC->GetLabel(), fkMC->GetTRDlabel(), fkMC->GetPDG(), fkMC->GetNTrackRefs(), Int_t(val[6]),
//     fSpecies, spc, Int_t(val[3]), fPt>0.?fPt:fkESD->Pt(), Int_t(val[5]), fkMC->GetTrackRef()?fkMC->GetTrackRef()->Pt():fkMC->GetPt(),
//     (status&AliESDtrack::kTRDin)?'x':'o', (status&AliESDtrack::kTRDout)?'x':'o', (status&AliESDtrack::kTRDStop)?'x':'o', Int_t(val[4]));

/*  if(val[0]<0.5 && TMath::Abs(val[4]+1)<0.1){

    if(fkTrack) printf("  Seed[%s] tracklets[%d]\n", fkTrack->IsTPCseeded()?"TPC":"ITS", Int_t(val[5]));
  }*/
  if(fkTrack) AliDebug(2, Form("%3d[%s] tracklets[%d] label[%2d %+2d] species[%d %d] in[%c] out[%c] stop[%c]",
    fkESD->GetId(), fkTrack->IsTPCseeded()?"TPC":"ITS", Int_t(val[5]), fkMC->GetLabel(), fkMC->GetTRDlabel(), fSpecies, spc,
    (status&AliESDtrack::kTRDin)?'x':'o', (status&AliESDtrack::kTRDout)?'x':'o', (status&AliESDtrack::kTRDStop)?'x':'o'));
  H->Fill(val);
  return NULL;
}

// //____________________________________________________________________
// TH1* AliTRDefficiency::PlotMC(const AliTRDtrackV1 *track)
// {
// // plot TRD efficiency based on MC info
// 
//   if(!HasMC()) return NULL;
//   if(!fkESD){
//     AliDebug(4, "No ESD info.");
//     return NULL;
//   }
//   if(!fkMC){
//     AliDebug(4, "No MC info.");
//     return NULL;
//   }
// 
//   THnSparse *H(NULL);
//   if(!fContainer || !(H = ((THnSparse*)fContainer->FindObject("hMC")))){
//     AliWarning("No output container defined.");
//     return NULL;
//   }
//   if(track) fkTrack = track;
//   Double_t val[11]; memset(val, 0, 11*sizeof(Double_t));
//   ULong_t status(fkESD->GetStatus());
//   val[0] =((status&AliESDtrack::kTRDin)?1:0) +
//           ((status&AliESDtrack::kTRDStop)?1:0) +
//           ((status&AliESDtrack::kTRDout)?2:0);
//   val[1] = fkESD->Phi();
//   val[2] = fkESD->Eta();
//   val[3] = DebugLevel()>=1?GetPtBin(fkESD->Pt()):GetPtBinSignificant(fkESD->Pt());
//   if(fkTrack){ // read track status in debug mode with friends
//     val[4] = fkTrack->GetStatusTRD(-1);
//     for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++) val[5+ily]=fkTrack->GetStatusTRD(ily);
//   }
//   H->Fill(val);
// 
// }

//____________________________________________________________________
void AliTRDefficiency::LocalUserExec(Option_t *)
{
  //
  // Do it obsolete
  //

  Int_t labelsacc[10000];
  memset(labelsacc, 0, sizeof(Int_t) * 10000);
	
  fTracks = dynamic_cast<TObjArray *>(GetInputData(1));
  if(!fTracks) return;
  if(!fTracks->GetEntriesFast()) return;
  else AliDebug(2, Form("Tracks[%d] for %s", fTracks->GetEntriesFast(), GetName()));
  if(!fMissed){ 
    fMissed = new TClonesArray("AliTRDtrackInfo", 10);
    fMissed->SetOwner();
  }

  Float_t mom;
  Int_t selection[10000], nselect = 0;
  ULong_t status; Int_t pidx;
  Int_t nTRD = 0, nTPC = 0, nMiss = 0;
  AliTRDtrackInfo     *track = NULL;
  AliTrackReference     *ref = NULL;
  AliExternalTrackParam *esd = NULL;
  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);

		if(!track->HasESDtrack()) continue;
    status = track->GetStatus();

    // missing TPC propagation - interesting for SA
    if(!(status&AliESDtrack::kTPCout)) continue;

    // missing MC info.
    if(HasMCdata() && track->GetNTrackRefs() <= 1) continue;
   
    nTPC++;
    selection[nselect++]=itrk;
    ref  = track->GetTrackRef(0);
    esd  = track->GetESDinfo()->GetOuterParam();
    mom  = ref ? ref->P(): esd->P();
    pidx = AliTRDCalPID::GetPartIndex(track->GetPDG());
    pidx = TMath::Max(pidx, 0);
    AliDebug(4, Form("PID: %d", pidx));

    //Int_t n = track->GetNumberOfClusters(); 
    // where are this tracklets ???
    //if(ncls0 > ncls1) printf("%3d ESD[%3d] TRD[%3d|%3d]\n", itrk, ncls0, ncls1, n);
    if(track->GetNumberOfClustersRefit()){ 
      ((TProfile*)fContainer->At(pidx))->Fill(mom, 1.);
			labelsacc[nTRD] = track->GetLabel();
      nTRD++;
      continue;
    }



    Float_t xmed, xleng;
    Int_t iref = 1; Bool_t found = kFALSE;
    while((ref = track->GetTrackRef(iref))){
      xmed = .5*(ref->LocalX() + track->GetTrackRef(iref-1)->LocalX());
      xleng= (ref->LocalX() - track->GetTrackRef(iref-1)->LocalX());
      if(TMath::Abs(xmed - 298.5) < .5 &&
        TMath::Abs(xleng - 3.7) < .1){ 
        found = kTRUE;
        break;
      }
      iref++;
    }
    if(!found){ 
      nTPC--;
      // track missing first layer. Maybe interesting for SA.
      continue;
    }
    nselect--;
    new ((*fMissed)[nMiss]) AliTRDtrackInfo(*track);
    nMiss++;
  }
  AliDebug(2, Form("%3d Tracks: ESD[%3d] TPC[%3d] TRD[%3d | %5.2f%%] Off[%d]", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), fTracks->GetEntriesFast(), nTPC, nTRD, nTPC ? 1.E2*nTRD/float(nTPC) : 0., fMissed->GetEntriesFast()));


  // Find double tracks
  Float_t threshold = 10.;
  AliTrackReference *refMiss = NULL;
  AliExternalTrackParam *op = NULL;
  AliTRDtrackInfo       *tt = NULL;
  for(Int_t imiss=0; imiss<nMiss; imiss++){
    //printf("Searching missing %d ...\n", imiss);

    // get outer param of missed
    tt = (AliTRDtrackInfo*)fMissed->UncheckedAt(imiss);
    op = tt->GetESDinfo()->GetOuterParam();
    Double_t alpha = op->GetAlpha(), cosa = TMath::Cos(alpha), sina = TMath::Sin(alpha);

    Double_t xyz[3], x0, y0, z0, x, y, z, dx, dy, dz, d;

    Bool_t bFOUND = kFALSE;
    for(Int_t iselect=0; iselect<nselect; iselect++){
      track = (AliTRDtrackInfo*)fTracks->UncheckedAt(selection[iselect]);

      // check first MC ... if available
      d = 0;
      for(Int_t iref=0; iref<track->GetNTrackRefs(); iref++){
        if(!(ref = track->GetTrackRef(iref))) continue;
        if((refMiss = tt->GetTrackRef(iref))){
          dy = ref->LocalY() - refMiss->LocalY();
          dz = ref->Z() - refMiss->Z();
        } else {
          // compare missOP with refTrackRef in LTC
          x0 = ref->LocalX();
          op->GetYAt(x0, AliTracker::GetBz(), y0);
          op->GetZAt(x0, AliTracker::GetBz(), z0);
          dy = y0 - ref->LocalY();
          dz = z0 - ref->Z();
        }
        d += (dy*dy + dz*dz);
      }
      //printf("\td[%d] = %f N[%d]\n", selection[iselect], d, track->GetNTrackRefs());
      if((track->GetNTrackRefs())){ 
        d /= track->GetNTrackRefs();
        if(d < threshold){
          //printf("\t\tFound %2d in ref[%3d] : d[%f]\n", imiss, selection[iselect], d/track->GetNTrackRefs());
          bFOUND = kTRUE; break;
        }
      }

      // process outer param ... always available
      // compare missOP with OP in GTC
      esd = track->GetESDinfo()->GetOuterParam();
      esd->GetXYZ(xyz);
      x0 = esd->GetX();
      op->GetYAt(x0, AliTracker::GetBz(), y0);
      op->GetZAt(x0, AliTracker::GetBz(), z0);
      x = x0*cosa - y0*sina;
      y = x0*sina + y0*cosa;
      z = z0;
      dx=xyz[0]-x;
      dy=xyz[1]-y;
      dz=xyz[2]-z;
      d = dx*dx+dy*dy+dz*dz;
      //printf("\td[%d] = %f op\n", selection[iselect], d);
      if(d < threshold){
        //printf("\t\tFound %2d in op[%3d]  : d[%f] dx[%5.2f] dy[%5.2f] dz[%5.2f]\n", imiss, selection[iselect], d, dx, dy, dz);
        bFOUND = kTRUE; break;
      }
    }
    if(bFOUND) nTPC--;
    else{ 
      ref = tt->GetTrackRef(0);
      mom = ref ? ref->P(): op->P();
      pidx = AliTRDCalPID::GetPartIndex(tt->GetPDG());
      pidx = TMath::Max(pidx, 0);
      ((TProfile*)fContainer->At(pidx))->Fill(mom, 0.);
      AliDebug(2, Form("  NOT bFOUND Id[%d] Mom[%f]\n", tt->GetTrackId(), mom));
    }
  }

  AliDebug(2, Form("%3d Tracks: ESD[%3d] TPC[%3d] TRD[%3d | %5.2f%%] Off[%d]", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), fTracks->GetEntriesFast(), nTPC, nTRD, nTPC ? 1.E2*nTRD/float(nTPC) : 0., fMissed->GetEntriesFast()));

  //fMissed->Delete();
	// check for double countings
	Int_t indices[10000]; memset(indices, 0, sizeof(Int_t) * 10000);
	TMath::Sort(nTRD, labelsacc, indices);
	if(DebugLevel() > 2){
	for(Int_t itk = 0; itk < nTRD - 1; itk++)
		if(labelsacc[indices[itk]] ==labelsacc[indices[itk + 1]]) printf("Double counted MC track: %d\n", labelsacc[indices[itk]]);
	}
}

//____________________________________________________________________
Bool_t AliTRDefficiency::GetRefFigure(Int_t ifig)
{
// Steer reference figures

  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }
  gPad->SetLogx();

  TLegend *leg(NULL);
  Bool_t bFIRST(kTRUE);
  TProfile *h(NULL);
  switch(ifig){
  case 0:
    h = (TProfile*)fContainer->At(AliPID::kSPECIES);
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      h->Add((TProfile*)fContainer->At(is));
    }
    h->SetMarkerStyle(24);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetTitle("TRD Efficiency integrated");
    h->SetXTitle("p [GeV/c]");
    h->GetXaxis()->SetMoreLogLabels();
    h->SetYTitle("Efficiency");
    h->GetYaxis()->CenterTitle();
    h->Draw("e1");
    break;
  case 1:
    bFIRST = kTRUE;
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      if(!(h = (TProfile*)fContainer->At(is))) continue;
      h->SetMarkerStyle(24);
      if(bFIRST){
        h->Draw("e1");
        h->SetXTitle("p [GeV/c]");
        h->GetXaxis()->SetMoreLogLabels();
        h->SetYTitle("Efficiency");
        h->GetYaxis()->CenterTitle();
        h->GetYaxis()->SetRangeUser(0.8, 1.05);
        leg=new TLegend(.7, .2, .98, .6);
        leg->SetHeader("Species");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h, h->GetTitle(), "pl");
      } else {
        leg->AddEntry(h, h->GetTitle(), "pl");
        h->Draw("same e1");
      }
      bFIRST = kFALSE;
    }
    if(leg) leg->Draw();
    break;
  }
  return kTRUE;
}

//________________________________________________________
TObjArray* AliTRDefficiency::Histos()
{
  //
  // Define histograms
  //

  if(fContainer) return fContainer;

  fContainer  = new TObjArray(1); fContainer->SetOwner(kTRUE);
  THnSparse *H(NULL);
  TString st;

  //++++++++++++++++++++++
  // cluster to detector
  if(!(H = (THnSparseI*)gROOT->FindObject("hEFF"))){
    const Int_t mdim(8);
    Int_t nlabel(3);
    const Char_t *eTitle[mdim] = {"status", "#phi [rad]", "eta", "p_{t} [bin]", "label", "p_{t}^{MC} [bin]", "N_{trklt}", "chg*spec*rc"};
    const Int_t eNbins[mdim]   = {  5,         144,        45,     fNpt-1,      nlabel,     fNpt-1, AliTRDgeometry::kNlayer-2, 5};
    const Double_t eMin[mdim]  = { -0.5,    -TMath::Pi(),  -.9,      -0.5,      -1.5,         -0.5,         1.5,             -2.5},
                   eMax[mdim]  = {  4.5,     TMath::Pi(),   .9,    fNpt-1.5,     1.5,       fNpt-1.5,       5.5,              2.5};
    st = "basic efficiency;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?mdim:(HasMCdata()?6:4);
    for(Int_t idim(0); idim<ndim; idim++){ st += eTitle[idim]; st+=";";}
    H = new THnSparseI("hEFF", st.Data(), ndim, eNbins, eMin, eMax);
  } else H->Reset();
  fContainer->AddAt(H, 0);

  return fContainer;
}

//____________________________________________________________________
Bool_t AliTRDefficiency::PostProcess()
{
// Fit, Project, Combine, Extract values from the containers filled during execution

  if (!fContainer) {
    AliError("ERROR: list not available");
    return kFALSE;
  }
  if(!fProj){
    AliInfo("Building array of projections ...");
    fProj = new TObjArray(2000); fProj->SetOwner(kTRUE);
  }
  // set pt/p segmentation. guess from data
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject("hEFF"))){
    AliError("Missing/Wrong data @ hEFF.");
    return kFALSE;
  }
  fNpt=H->GetAxis(3)->GetNbins()+1;
  if(!MakeMomSegmentation()) return kFALSE;

  if(!MakeProjectionBasicEff()) return kFALSE;
  return kTRUE;
}

//____________________________________________________________________
Bool_t AliTRDefficiency::MakeProjectionBasicEff()
{
// Make basic efficiency plots

  if(!fContainer || !fProj){
    AliError("Missing data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject("hEFF"))){
    AliError("Missing/Wrong data @ hEFF.");
    return kFALSE;
  }
  Int_t ndim(H->GetNdimensions()); //Bool_t debug(ndim>Int_t(kNdimCl));
  TAxis *aa[11], *al(NULL), *apt(NULL); memset(aa, 0, sizeof(TAxis*) * 11);
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(H->GetNdimensions() > 4) al = H->GetAxis(4);
  if(H->GetNdimensions() > 5) apt = H->GetAxis(5);
  Int_t nlab=al?5:1;

  // define rebinning strategy
  //const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  const Int_t nprojs(50);
  AliTRDrecoProjection hp[nprojs];  TObjArray php(nprojs);
  const Char_t *stat[] = {"!TRDin", "TRDin", "TRDin&TRDStop", "TRDin&TRDout", "TRDin&TRDout&TRDStop"};
  const Char_t *lab[] = {"MC-Miss", "MC-Trkble", "MC-Good", "MC-Accept", "MC-Wrong"};
  Int_t ih(0);
  for(Int_t ilab(0); ilab<nlab; ilab++){
    for(Int_t istat(0); istat<5; istat++){
      hp[ih].Build(Form("HEff0%d%d", ilab, istat),
                    Form("Efficiency ::  Stat[#bf{%s}] %s", stat[istat], nlab>1?lab[ilab]:""), 2, 1, 3, aa);
      //hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php.AddLast(&hp[ih++]);
      if(!apt) continue;
      hp[ih].Build(Form("HEff1%d%d", ilab, istat),
                    Form("Efficiency [MC] ::  Stat[#bf{%s}] %s", stat[istat], lab[ilab]), 2, 1, 5, aa);
      //hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php.AddLast(&hp[ih++]);
    }
  }
  AliInfo(Form("Build %3d 3D projections.", ih));

  AliTRDrecoProjection *pr0(NULL), *pr1(NULL);
  Int_t istatus, ilab(0), coord[11]; memset(coord, 0, sizeof(Int_t) * 11); Double_t v = 0.;
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord); if(v<1.) continue;
    istatus = coord[0]-1; if(istatus>4) continue;
    ilab=al?coord[4]:0;
    Int_t isel = (apt?2:1)*(ilab*5+istatus);
    if(isel>=ih){
      AliError(Form("Wrong selection %d [%3d] {stat[%d] lab[%d]}", isel, ih, istatus, ilab));
      return kFALSE;
    }
    if(!(pr0=(AliTRDrecoProjection*)php.At(isel))) {
      AliError(Form("Missing projection %d", isel));
      return kFALSE;
    }
    if(strcmp(pr0->H()->GetName(), Form("HEff0%d%d", ilab, istatus))!=0){
      AliError(Form("Projection mismatch :: request[HEff0%d%d] found[%s]", ilab, istatus, pr0->H()->GetName()));
      return kFALSE;
    }
    //for(Int_t jh(0); jh<1/*np[isel]*/; jh++) ((AliTRDrecoProjection*)php.At(isel+jh))->Increment(coord, v);
    AliDebug(2, Form("Found %s for selection stat[%d] lab[%d]", pr0->H()->GetName(), istatus, ilab));
    ((AliTRDrecoProjection*)php.At(isel))->Increment(coord, v);
    if(apt) ((AliTRDrecoProjection*)php.At(isel+1))->Increment(coord, v);
  }
  if(HasDump3D()){
    TDirectory *cwd = gDirectory;
    TFile::Open(Form("Dump%s_%s.root", GetName(), H->GetName()), "RECREATE");
    for(Int_t ip(0); ip<php.GetEntriesFast(); ip++){
      if(!(pr0 = (AliTRDrecoProjection*)php.At(ip))) continue;
      if(!pr0->H()) continue;
      TH3 *h3=(TH3*)pr0->H()->Clone();
      h3->Write();
    }
    gFile->Close();
    cwd->cd();
  }

  Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].H()) continue;
    if(hp[ih].H()->Integral()<5.) continue;
    for(Int_t ipt(0); ipt<=fNpt; ipt++) fProj->AddAt(hp[ih].Projection2Dbin(ipt), jh++);
  }

  // PROCESS MC LABEL, CLONE 
  if(nlab>1){
    // MISSED TRACKS
    for(Int_t imc(0); imc<(apt?2:1); imc++){
      AliTRDrecoProjection prMis, prMisOK, prS, prSs, prTRDf;
      for(Int_t istat(0); istat<5; istat++) {
        if(!(pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HEff%d0%d", imc, istat)))) continue;
        if(!prMis.H()){
          prMis = (*pr0);
          prMis.SetNameTitle(Form("HEff%dM", imc), "Missed tracks");
          prMis.H()->SetNameTitle(Form("HEff%dM", imc), Form("Efficiency [%s] :: Missed tracks", imc?"MC":"REC"));
        } else prMis+=(*pr0);
        if(istat==0) {             // MC miss & !TRDin
          prMisOK = (*pr0);
          prMisOK.SetNameTitle(Form("HEff%dMok", imc), "Missed tracks/Not propagated");
          prMisOK.H()->SetNameTitle(Form("HEff%dMok", imc), Form("Efficiency [%s] :: Missed tracks/Not propagated", imc?"MC":"REC"));
        } else if(istat==1) {             // MC miss & TRDin
          prS = (*pr0);
          prS.SetNameTitle(Form("HEff%dMS", imc), "Missed seeded tracks");
          prS.H()->SetNameTitle(Form("HEff%dMS", imc), Form("Efficiency [%s] :: Missed tracks/Seeder propagated", imc?"MC":"REC"));
        } else if(istat==2) {             // MC miss & TRDstop
          prSs = (*pr0);
          prSs.SetNameTitle(Form("HEff%dMSs", imc), "Stop seeded tracks");
          prSs.H()->SetNameTitle(Form("HEff%dMSs", imc), Form("Efficiency [%s] :: Missed tracks/Stop seeded tracks", imc?"MC":"REC"));
        } else if(istat==3 || istat==4) {             // MC miss & (TRDout || TRDout+TRDstop)
          prTRDf = (*pr0);
          prTRDf.SetNameTitle(Form("HEff%dMTtrd", imc), "Missed tracks/TRD fooled");
          prTRDf.H()->SetNameTitle(Form("HEff%dMTtrd", imc), Form("Efficiency [%s] :: Missed tracks/TRD fooled", imc?"MC":"REC"));
        } else {
          if(Int_t(pr0->H()->Integral()) > 0) AliWarning(Form("Detected %3d MISSED entries for %s[%s]", Int_t(pr0->H()->Integral()), stat[istat], imc?"MC":"REC"));
        }
      }
      for(Int_t ipt(-1); ipt<=fNpt; ipt++){
        fProj->AddAt(prMis.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prMisOK.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prS.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prSs.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prTRDf.Projection2Dbin(ipt, imc), jh++);
      }
    
      // TRACKABLE TRACKS
      AliTRDrecoProjection prAll, prTRDns, prTRDfail, prTRD, prTRDbad;//, prTRDstop;
      for(ilab=1; ilab<=4; ilab++){
        for(Int_t istat(0); istat<5; istat++) {
          if(!(pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HEff%d%d%d", imc, ilab, istat)))){
            AliError(Form("Missing prj. HEff%d%d%d", imc, ilab, istat));
            continue;
          }
          if(!prAll.H()){
            prAll = (*pr0);
            prAll.SetNameTitle(Form("HEff%dAT", imc), "Trackable tracks");
            prAll.H()->SetNameTitle(Form("HEff%dAT", imc), Form("Efficiency [%s] :: Trackable tracks", imc?"MC":"REC"));
          } else prAll+=(*pr0);
          if(ilab==1 && istat==0){      // MC trackable & !TRDin
            prTRDns = (*pr0);
            prTRDns.SetNameTitle(Form("HEff%dAnStrd", imc), "Trackable tracks/TRD not seeded");
            prTRDns.H()->SetNameTitle(Form("HEff%dAnStrd", imc), Form("Efficiency [%s] :: Trackable tracks/TRD not seeded", imc?"MC":"REC"));
          } else if(ilab==1 && istat==1){      // MC trackable & TRDin
            prTRDfail = (*pr0);
            prTRDfail.SetNameTitle(Form("HEff%dAFtrd", imc), "Trackable tracks/TRD failed");
            prTRDfail.H()->SetNameTitle(Form("HEff%dAFtrd", imc), Form("Efficiency [%s] :: Trackable tracks/TRD failed", imc?"MC":"REC"));
          } else if((ilab==2 && istat>=3) ||   // MC good & (TRDout || TRDout+TRDstop)
                    (ilab==3 && istat>=3) ||   // MC accept & (TRDout || TRDout+TRDstop)
                    (ilab==4 && istat==2)) {   // MC bad & TRDin+TRDstop
            if(!prTRD.H()){
              prTRD = (*pr0);
              prTRD.SetNameTitle(Form("HEff%dAtrd", imc), "Trackable tracks/TRD OK");
              prTRD.H()->SetNameTitle(Form("HEff%dAtrd", imc), Form("Efficiency [%s] :: Trackable tracks/TRD OK", imc?"MC":"REC"));
            } else prTRD+=(*pr0);
          } else if(ilab==4 && istat>=3){      // MC trackable & (TRDout (+TRDstop))
            if(!prTRDbad.H()){
              prTRDbad = (*pr0);
              prTRDbad.SetNameTitle(Form("HEff%dABtrd", imc), "Trackable tracks/TRD bad");
              prTRDbad.H()->SetNameTitle(Form("HEff%dABtrd", imc), Form("Efficiency [%s] :: Trackable tracks/TRD bad", imc?"MC":"REC"));
            } else  prTRDbad+=(*pr0);
          } else {
            if(Int_t(pr0->H()->Integral()) > 0) AliWarning(Form("Detected %3d TRACKABLE entries for %s -> %s[%s]", Int_t(pr0->H()->Integral()), lab[ilab], stat[istat], imc?"MC":"REC"));
          }
        }
      }
      for(Int_t ipt(-1); ipt<=fNpt; ipt++){
        fProj->AddAt(prAll.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prTRDns.Projection2Dbin(ipt), jh++);
        fProj->AddAt(prTRDfail.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prTRD.Projection2Dbin(ipt, imc), jh++);
        fProj->AddAt(prTRDbad.Projection2Dbin(ipt, imc), jh++);
      }
    } // END pt source LOOP
  } // END PROCESS EFFICIENCY MC

  // PROCESS EFFICIENCY NOMC
  const char suffix[] = {'A', 'S', 'T'};
  const char *sname[] = {"All", "Stopped", "Tracked [Stop]"};
  for(Int_t istat(0); istat<5; istat++) {
    if(!(pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HEff00%d", istat)))){
      AliError(Form("Missing prj. HEff00%d", istat));
      continue;
    }
    // sum over MC labels if available
    for(ilab=1; ilab<nlab; ilab++){
      if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HEff0%d%d", ilab, istat)))){
        AliError(Form("Missing prj. HEff0%d%d", ilab, istat));
        continue;
      }
      (*pr0)+=(*pr1);
    }
    pr0->H()->SetNameTitle(Form("HEff%d", istat), Form("Efficiency :: Stat[#bf{%s}]", stat[istat]));
    for(Int_t ipt(0); ipt<=fNpt; ipt++) fProj->AddAt(pr0->Projection2Dbin(ipt), jh++);

    if(istat>1 && (pr1 = (AliTRDrecoProjection*)php.FindObject("HEff001"))) (*pr1)+=(*pr0);
    //if(istat>2 && (pr1 = (AliTRDrecoProjection*)php.FindObject("HEff02"))) (*pr1)+=(*pr0);
    if(istat>3 && (pr1 = (AliTRDrecoProjection*)php.FindObject("HEff003"))) (*pr1)+=(*pr0);
  }
  // Project 2D tracks
  for(Int_t istat(0); istat<3; istat++){
    if(!(pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HEff00%d", istat+1)))){
      AliError(Form("Missing prj. HEff00%d", istat+1));
      continue;
    }
    pr0->H()->SetNameTitle(Form("HEff%c", suffix[istat]), Form("Efficiency :: %s Tracks", sname[istat]));
    for(Int_t ipt(-1); ipt<=fNpt; ipt++) fProj->AddAt(pr0->Projection2Dbin(ipt), jh++);
  }

  // Efficiency
/*  TH2F *h2T(NULL), *h2P(NULL);
  for(Int_t ipt(-1); ipt<=fNpt; ipt++){
    if(!(h2T = (TH2F*)fProj->FindObject(Form("HEffT%d_2D", ipt)))) continue;
    if(!(h2P = (TH2F*)fProj->FindObject(Form("HEffP%d_2D", ipt)))) continue;
    h2P->Divide(h2T); PutTrendValue(ipt<0?"pt":(Form("pt%d", ipt)), GetMeanStat(h2P, 0.01, ">"));
  }*/
  AliInfo(Form("Done %3d 2D projections.", jh));
  return kTRUE;
}

//____________________________________________________________________
void AliTRDefficiency::MakeSummary()
{
//  Build summary plots
  if(!fProj){
    AliError("Missing results");
    return;
  }
  TVirtualPad *p(NULL); TCanvas *cOut(NULL);
  TH2 *h2(NULL);
  gStyle->SetPalette(1);

  Int_t nbins(DebugLevel()==0?3:20);
  //calculate true pt bin
  Float_t ptBins[23]; ptBins[0] = 0.;
  if(nbins==3){ // significant bins
    ptBins[1] = 0.5;
    ptBins[2] = 0.8;
    ptBins[3] = 1.5;
    ptBins[4] = 5.;
    ptBins[5] = 10.;
  } else if(nbins==20){ // logarithmic bins
    ptBins[1] = 0.5;
    Float_t dpt(0.002);
    for(Int_t ipt(1); ipt<21; ipt++) ptBins[ipt+1] = ptBins[ipt]+(TMath::Exp(ipt*ipt*dpt)-1.);
    ptBins[22] = 10.;
  } else {
    AliError(Form("Unknown no.[%d] of bins in the p_t spectrum", nbins));
    return;// kFALSE;
  }

  const Char_t cid[]={'T','P'};
  cOut = new TCanvas(Form("%s_Eff", GetName()), "TRD Efficiency", 1536, 1536); cOut->Divide(2,2,1.e-5,1.e-5);
  // tracking eff :: eta-phi dependence
  for(Int_t it(0); it<2; it++){
    if(!(h2 = (TH2*)fProj->FindObject(Form("HEff%cEnN", cid[it])))) {
      AliError(Form("Missing \"HEff%cEnN\".", cid[it]));
      continue;
    }
    h2->SetContour(10); h2->Scale(1.e2); SetRangeZ(h2, 80, 100, 5);
    h2->GetZaxis()->SetTitle("Efficiency [%]"); h2->GetZaxis()->CenterTitle();
    p=cOut->cd(it+1);p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
    h2->Draw("colz");
    MakeDetectorPlot();
  }
  if(!(h2 = (TH2*)fProj->FindObject("HEff0En"))) {
    AliError("Missing \"HEff0En\".");
    return;
  }
  p=cOut->cd(3);p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
  h2->Draw("colz"); MakeDetectorPlot();
  // tracking eff :: pt dependence
  TH1 *h[2] = {0};
  if(!(h[0] = (TH1*)fProj->FindObject("HEffP_zN"))){
    AliError("Missing \"HEffP_zN\".");
    return;
  }
  if(!(h[1] = (TH1*)fProj->FindObject("HEffT_zN"))){
    AliError("Missing \"HEffT_zN\".");
    return;
  }
  TH1 *h1[3] = {0};
  Color_t color[] = {kGreen, kBlue, kRed};
  for(Int_t il=0;il<3;il++){
    h1[il]=new TH1F(Form("h1Eff%d", il), "", nbins+2, ptBins);
    h1[il]->SetFillColor(color[il]);
    h1[il]->SetFillStyle(il==2?3002:1001);
    h1[il]->SetLineColor(color[il]);
    h1[il]->SetMarkerStyle(4);
    h1[il]->SetMarkerColor(color[il]);
    h1[il]->SetLineWidth(1);
  }
  for(Int_t ip(0);ip<=(nbins+1);ip++){
    h1[0]->SetBinContent(ip+1, 1.e2*h[0]->GetBinContent(ip)); // propagated
    h1[1]->SetBinContent(ip+1, 1.e2*(h[1]->GetBinContent(ip) - h[0]->GetBinContent(ip))); // stopped
    h1[2]->SetBinContent(ip+1, 1.e2*(1 - h[1]->GetBinContent(ip))); // missed
  }
  const Char_t *labEff[] = {"Propagated", "Stopped", "Missed"};
  THStack *hs = new THStack("hEff","Tracking Efficiency;p_{t} [GeV/c];Efficiency [%]");
  TLegend *leg = new TLegend(0.671371,0.1313559,0.9576613,0.2923729,NULL,"brNDC");
  leg->SetBorderSize(0); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
  for(Int_t ic(0); ic<3;ic++){ hs->Add(h1[ic]);leg->AddEntry(h1[ic], labEff[ic], "f");}
  p=cOut->cd(4); p->SetLeftMargin(0.08266129); p->SetRightMargin(0.0141129);p->SetTopMargin(0.006355932);p->SetLogx();
  hs->Draw(); //hs->Draw("same nostack,e1p");
  leg->Draw();
  hs->GetXaxis()->SetRangeUser(0.6,10.);
  hs->GetXaxis()->SetMoreLogLabels();
  hs->GetXaxis()->SetTitleOffset(1.2);
  hs->GetYaxis()->SetNdivisions(513);
  hs->SetMinimum(75.);
  hs->GetYaxis()->CenterTitle();
  cOut->SaveAs(Form("%s.gif", cOut->GetName()));

  cOut = new TCanvas(Form("%s_MC", GetName()), "TRD Label", 1536, 1536); cOut->Divide(2,2,1.e-5,1.e-5);
  for(Int_t ipad(0); ipad<3; ipad++){
    p=cOut->cd(ipad+1);p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
    if(!(h2 = (TH2*)fProj->FindObject(Form("HEffLb%dEnN", ipad)))) continue;
    h2->SetContour(10);
    h2->Scale(1.e2); SetRangeZ(h2, ipad==1?50:0., ipad==1?90.:50., ipad==1?0.01:0.01);
    h2->GetZaxis()->SetTitle("Efficiency [%]"); h2->GetZaxis()->CenterTitle();
    h2->Draw("colz");
    MakeDetectorPlot();
  }
  color[0] = kRed; color[1] = kGreen; color[2] = kBlue;
  for(Int_t il=0;il<3;il++){
    if(!(h[il] = (TH1D*)fProj->FindObject(Form("HEffLb%d_zN", il)))) continue;
    h1[il]=new TH1F(Form("h1Lab%d", il), "", nbins+2, ptBins);
    for(Int_t ip(0);ip<=(nbins+1);ip++) h1[il]->SetBinContent(ip+1, 1.e2*h[il]->GetBinContent(ip));
    h1[il]->SetFillColor(color[il]);
    h1[il]->SetFillStyle(il==2?3002:1001);
    h1[il]->SetLineColor(color[il]);
    h1[il]->SetLineWidth(1);
  }
  const Char_t *labMC[] = {"TRD != ESD [bad]", "TRD == ESD [good]", "TRD == -ESD [accept]"};
  leg = new TLegend(0.671371,0.1313559,0.9576613,0.2923729,NULL,"brNDC");
  leg->SetBorderSize(0); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
  hs = new THStack("hLab","TRD Label;p_{t} [GeV/c];Efficiency [%]");
  hs->Add(h1[1]);leg->AddEntry(h1[1], labMC[1], "f"); // good
  hs->Add(h1[2]);leg->AddEntry(h1[2], labMC[2], "f"); // accept
  hs->Add(h1[0]);leg->AddEntry(h1[0], labMC[0], "f"); // bad
  p=cOut->cd(4); p->SetLeftMargin(0.08266129); p->SetRightMargin(0.0141129);p->SetTopMargin(0.006355932); p->SetLogx();
  hs->Draw(/*"nostack,e1p"*/); leg->Draw();
  cOut->Modified();cOut->Update();
  hs->GetXaxis()->SetRangeUser(0.6,10.);
  hs->GetXaxis()->SetMoreLogLabels();
  hs->GetXaxis()->SetTitleOffset(1.2);
  hs->GetYaxis()->SetNdivisions(513);
  hs->SetMinimum(50.);
  hs->GetYaxis()->CenterTitle();
  cOut->SaveAs(Form("%s.gif", cOut->GetName()));
}
