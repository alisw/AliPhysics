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

/* $Id: AliTRDefficiencyMC.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TPad.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include "TTreeStream.h"

#include "AliMagF.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliMathBase.h"
#include "AliTrackReference.h"
#include "AliAnalysisManager.h"

#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "Cal/AliTRDCalPID.h"
#include "info/AliTRDtrackInfo.h"
#include "AliTRDinfoGen.h"
#include "AliTRDefficiencyMC.h"

ClassImp(AliTRDefficiencyMC)
Float_t AliTRDefficiencyMC::fgPCut   = 0.2; //[GeV/c]
Float_t AliTRDefficiencyMC::fgPhiCut = 50.; //[deg]
Float_t AliTRDefficiencyMC::fgThtCut = 50.; //[deg]
//_____________________________________________________________________________
AliTRDefficiencyMC::AliTRDefficiencyMC()
  :AliTRDrecoTask()
{
  //
  // Default constructor
  //
}

AliTRDefficiencyMC::AliTRDefficiencyMC(char* name)
  :AliTRDrecoTask(name, "Combined Tracking Efficiency")
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
void AliTRDefficiencyMC::UserCreateOutputObjects(){
  //
  // Create output objects
  //

  fContainer = Histos();
}

//_____________________________________________________________________________
void AliTRDefficiencyMC::UserExec(Option_t *){
  //
  // Execute the task:
  //
  // Loop over TrackInfos
  // 1st: check if there is a trackTRD
  // 2nd: put conditions on the track:
  //      - check if we did not register it before
  //      - check if we have Track References for the track 
  // 3rd: Register track:
  //      - accepted if both conditions are fulfilled
  //      - contamination if at least one condition is not fulfilled
  // 4th: check Monte-Carlo Track wheter findable or not if there is no TRD track in this track info
  // 5th: register MC track as rejected if findable and not jet registered
  // Match the registers accepted and rejected and clean register rejected
  // Fill the histograms
  //
  const Int_t kArraySize = 10000;     // Store indices of track references in an array
  Int_t indexAccept[kArraySize],
        indexReject[kArraySize],
        indexContam[kArraySize];
  memset(indexAccept, 0, sizeof(Int_t) * kArraySize);
  memset(indexReject, 0, sizeof(Int_t) * kArraySize);
  memset(indexContam, 0, sizeof(Int_t) * kArraySize);
  Int_t naccept(0), 
        nreject(0), 
        nfindnt(0), 
        nkink(0), 
        ncontam(0);
  Bool_t isContamination = kFALSE;
  
  fTracks = dynamic_cast<TObjArray *>(GetInputData(1));
  if(!fTracks) return;
  Int_t nTrackInfos = fTracks->GetEntriesFast();
  AliDebug(2, Form("   CANDIDATE TRACKS[%d]", nTrackInfos));

  AliTRDtrackV1 *trackTRD(NULL);
  AliTRDtrackInfo *trkInf(NULL);
  for(Int_t itinf = 0; itinf < nTrackInfos; itinf++){
    trkInf = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(itinf));
    if(!trkInf) continue;

    if(trkInf->GetTrack() || trkInf->GetNumberOfClustersRefit()){
      isContamination = (IsRegistered(trkInf,indexAccept,naccept)>=0);
      if(!trkInf->GetNTrackRefs()){
        // We reject the track since the Monte Carlo Information is missing
        AliDebug(2, Form("MC(Track Reference) missing @ label[%d]", trkInf->GetLabel()));
        isContamination = kTRUE;
        // Debugging
        if(trackTRD && DebugLevel()>5) FillStreamTrackWOMC(trkInf);
      }
      if(isContamination){
        // reject kink (we count these only once)
        if(trkInf->GetKinkIndex()){ 
          AliDebug(4, Form("  track @ idx[%d] MC[%d] is kink.", itinf, trkInf->GetLabel()));
          nkink++;
          continue;
        }
        // Register track as contamination
        AliDebug(4, Form("  track @ idx[%d] MC[%d] is contamination.", itinf, trkInf->GetLabel()));
        indexContam[ncontam++]=itinf;
        continue;
      }
      // Accept track
      AliDebug(4, Form("  track @ idx[%d] is ACCEPTED.", itinf));
      // Register track as accepted
      indexAccept[naccept++] = itinf;
    }else{
      Int_t code(0);
      if((code=IsFindableNot(trkInf))){
        AliDebug(4, Form("  track @ idx[%d] MC[%d] not findable [%d].", itinf, trkInf->GetLabel(), code));
        nfindnt++;
      } else {
        // register track as rejected if not already registered there
        // Attention:
        // These track infos are not!!! registered as contamination
        if(IsRegistered(trkInf, indexReject, nreject)<0){ 
          AliDebug(4, Form("  track @ idx[%d] MC[%d] is missed.", itinf, trkInf->GetLabel()));
          indexReject[nreject++] = itinf;
        }
      }
    }
  }
  AliDebug(2, Form("TRACKS STATISTICS naccept[%d] ncontam[%d] nkink[%d] nmissed[%d] nfindnt[%d] ALL[%d] LOST[%d]", naccept, ncontam, nkink, nreject, nfindnt, naccept+ncontam+nkink, nTrackInfos-(naccept+nreject+ncontam+nkink+nfindnt)));

  // we have to check if the rejected tracks are registered as found
  // a possible source for this:
  // first the track is not found by the barrel tracking but it is later found
  // by the stand alone tracking, then two track info objects with the same 
  // label would be created
  // Attention:
  // these tracks are not! registered as contamination
  Int_t tmprejected[kArraySize]; Int_t nrej = nreject;
  memcpy(tmprejected, indexReject, sizeof(Int_t) * nreject);
  nreject = 0;
  for(Int_t irej = 0; irej < nrej; irej++){
    trkInf = dynamic_cast<AliTRDtrackInfo *>(fTracks->At(tmprejected[irej]));
    Int_t idx(-1);
    if((idx = IsRegistered(trkInf,indexAccept,naccept))<0){
      indexReject[nreject++] = tmprejected[irej];
    }else{
      //printf("tracks @ accept[%d] missed[%d] are the same.\n", indexAccept[idx], tmprejected[irej]);
    }
  }

  // Fill Histograms
  FillHistograms(naccept, &indexAccept[0], kAccept);
  FillHistograms(nreject, &indexReject[0], kMiss);
  FillHistograms(ncontam, &indexContam[0], kFake);

  Int_t nall(naccept + nreject);
  AliInfo(Form("%3d Tracks: MC[%3d] TRD[%3d | %5.2f%%] Fake[%3d | %5.2f%%]", 
    (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), 
    nall, 
    naccept, 
    (nall ? 1.E2*Float_t(naccept)/Float_t(nall) : 0.), 
    ncontam, 
    (nall ? 1.E2*Float_t(ncontam)/Float_t(nall) : 0.)));

  PostData(1, fContainer);
}


//_____________________________________________________________________________
Bool_t AliTRDefficiencyMC::PostProcess()
{
  //
  // Post Process 
  //
  // Change the histogram style
  // For species histograms apply the colors connected with the given particle species
  fNRefFigures = 8;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDefficiencyMC::GetRefFigure(Int_t ifig){
  //
  // Plot the histograms
  //
  if(ifig >= fNRefFigures) return kFALSE;
  if(!gPad) return kFALSE;
  gPad->SetLogx(kTRUE);
  if(ifig < 2){
    (dynamic_cast<TH1 *>(fContainer->At(ifig)))->Draw("e1");
    return kTRUE;
  }
  TH1 *h(NULL); 
  TLegend *leg=new TLegend(.65, .12, .85, .3);
  leg->SetHeader("Charge");
  leg->SetBorderSize(1);leg->SetFillColor(kWhite);
  switch(ifig){
  case 2:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+1));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  case 3:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+2));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+3));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  case 4:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+4));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+5));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  case 5:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+6));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+7));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  case 6:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+8));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+9));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  case 7:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+10));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+11));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  case 8:
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+12));
    h->Draw("e1"); leg->AddEntry(h, "  -", "pl");
    h=dynamic_cast<TH1 *>(fContainer->At(kEfficiencySpeciesHistogram+13));
    h->Draw("e1same"); leg->AddEntry(h, "  +", "pl");
    leg->Draw();
    break;
  }
  return kTRUE;
}

//_____________________________________________________________________________
TObjArray *AliTRDefficiencyMC::Histos(){
  //
  // Create the histograms
  //

  if(fContainer) return fContainer;
  const Int_t nbins = AliTRDCalPID::kNMom;
  Float_t xbins[nbins+1] = {fgPCut, .7, .9, 1.3, 1.7, 2.4, 3.5, 4.5, 5.5, 7., 9., 11.};
  const Int_t marker[2][AliPID::kSPECIES+1] = {
    {20, 21, 22, 23, 29, 2},
    {24, 25, 26, 27, 30, 5}
  };

  fContainer = new TObjArray();fContainer->Expand(14);

  TH1 *h(NULL);
  fContainer->AddAt(h=new TProfile("hEff", "Tracking Efficiency ALL", nbins, xbins), kEfficiencyHistogram);
  h->SetMarkerStyle(22);
  h->SetMarkerColor(kBlue);
  h->GetXaxis()->SetTitle("p [GeV/c]");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetTitle("Efficiency");
  h->GetYaxis()->SetRangeUser(0.2, 1.1);
  fContainer->AddAt(h=new TProfile("hFake", "Fake Tracks", nbins, xbins), kContaminationHistogram);
  h->SetMarkerStyle(22);
  h->SetMarkerColor(kBlue);
  h->GetXaxis()->SetTitle("p [GeV/c]");
  h->GetXaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetTitle("Contamination");

  Char_t sign[]={'+', '-'};
  for(Int_t isign = 0; isign < 2; isign++){
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      fContainer->AddAt(h=new TProfile(
        Form("hEff_%s%c", AliPID::ParticleShortName(ispec), sign[isign]), 
        Form("Tracking Efficiency for %s", AliPID::ParticleName(ispec)), nbins, xbins), 
        kEfficiencySpeciesHistogram+ispec*2+isign);
      h->SetMarkerStyle(marker[isign][ispec]);
      h->SetLineColor(AliTRDCalPID::GetPartColor(ispec));
      h->SetMarkerColor(kBlack);
      h->GetXaxis()->SetTitle("p [GeV/c]");
      h->GetXaxis()->SetMoreLogLabels();
      h->GetYaxis()->SetTitle("Efficiency");
      h->GetYaxis()->SetRangeUser(0.2, 1.1);
    }

    fContainer->AddAt(h=new TProfile(Form("hEff_PID%c", sign[isign]), "Tracking Efficiency no PID", nbins, xbins), kEfficiencySpeciesHistogram+AliPID::kSPECIES*2+isign);
    h->SetMarkerStyle(marker[isign][AliPID::kSPECIES]);
    h->SetMarkerColor(kBlack);h->SetLineColor(kBlack);
    h->GetXaxis()->SetTitle("p [GeV/c]");
    h->GetXaxis()->SetMoreLogLabels();
    h->GetYaxis()->SetTitle("Efficiency");
    h->GetYaxis()->SetRangeUser(0.2, 1.1);
  }
  return fContainer;
}

//_____________________________________________________________________________
Int_t AliTRDefficiencyMC::IsFindableNot(AliTRDtrackInfo * const trkInf){
  //
  // Apply Cuts on the Monte Carlo track references
  // return whether track is findable or not
  //
  

  const Float_t chmbHght = AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght();
  const Float_t eps(1.E-3);
  Int_t ntr(trkInf->GetNTrackRefs());

  AliDebug(10, Form("  CANDIDATE TrackRefs[%d]", ntr));
  // Check if track is findable
  Double_t mom(0.), phi(0.), tht(0.);
  Float_t xmin = 10000.0, xmax = 0.0; 
  Float_t ymin = 0.0, ymax = 0.0;
  Float_t zmin = 0.0, zmax = 0.0;
  Float_t lastx = 0.0, x = 0.0;
  Int_t nLayers(0), ntrTRD(0);
  Int_t sector[20];
  AliTrackReference *trackRef(NULL);
  for(Int_t itr(0); itr<ntr; itr++){
    if(!(trackRef = trkInf->GetTrackRef(itr))) continue;
    x = trackRef->LocalX(); 
    // Be Sure that we are inside TRD
    if(x < AliTRDinfoGen::GetEndTPC() || x > AliTRDinfoGen::GetEndTRD()) continue;	
    sector[ntrTRD] = Int_t(trackRef->Alpha()/AliTRDgeometry::GetAlpha());
    AliDebug(10, Form("    [%2d] x[%7.2f] y[%7.2f] z[%7.2f] Sec[%2d]", itr, trackRef->LocalX(), trackRef->LocalY(), trackRef->Z(), sector[ntrTRD]));
    if(x < xmin){
      xmin = trackRef->LocalX();
      ymin = trackRef->LocalY();
      zmin = trackRef->Z();
      mom  = trackRef->P();
    } else if(x > xmax){
      xmax = trackRef->LocalX();
      ymax = trackRef->LocalY();
      zmax = trackRef->Z();
    }
    if(itr > 0){
      Float_t dist = TMath::Abs(x - lastx);
      if(TMath::Abs(dist - chmbHght) < eps && sector[ntrTRD]==sector[0]){ 
        AliDebug(10, Form("    dx = %7.2f", dist));
        nLayers++;
      }
    }
    lastx = x;
    ntrTRD++; if(ntrTRD>=20) break;
  }
  Double_t dx(xmax - xmin);
  if(TMath::Abs(dx)<eps) return kNoChmb;

  phi = (ymax -ymin)/dx;
  tht = (zmax -zmin)/dx;
  phi=TMath::ATan(phi)*TMath::RadToDeg();
  tht=TMath::ATan(tht)*TMath::RadToDeg();
  Bool_t primary = trkInf->IsPrimary();
  const AliTRDtrackInfo::AliESDinfo *esd(trkInf->GetESDinfo());
  AliDebug(10, Form("    p=%6.3f[GeV/c] phi=%6.2f[deg] theta=%6.2f[deg] nLy[%d]", 
      mom, phi, tht, nLayers));
  if(DebugLevel()){
    (*DebugStream()) << "IsFindable"
      << "P="       << mom
      << "Phi="     << phi
      << "Tht="     << tht
      << "Ntr="     << ntrTRD
      << "NLy="     << nLayers
      << "Primary=" << primary
      << "\n";
  }

  // Apply cuts
  if(!nLayers) return kNoChmb;
  if(xmax < xmin) return kCurved;
  if(mom < fgPCut) return kPCut;


  if(TMath::Abs(phi) > fgPhiCut) return kPhiCut;
  if(TMath::Abs(tht) > fgThtCut) return kThtCut;

  if(nLayers < 4){
    if(!esd)return kLayer;
    if(!(esd->GetStatus() & AliESDtrack::kTPCout)) return kLayer;
  }

  //if(!trkInf->IsPrimary()) {failCode=kPrimary; return kFALSE;}

  return kFindable;
}

//_____________________________________________________________________________
void AliTRDefficiencyMC::FillHistograms(Int_t nTracks, Int_t *indices, ETRDefficiencyMCstatus mode){
  //
  // Fill Histograms in three different modes:
  // 1st tracks which are found and accepted
  // 2nd tracks which are not found and not already counted
  // 3rd contaminating tracks: either double counts (kinks!) or tracks with no MC hit inside TRD
  //
  
  TDatabasePDG *dbPDG(TDatabasePDG::Instance());
  Double_t trkmom(0.);   // the track momentum
  Int_t trkpdg(-1);      // particle PDG code
  AliTRDtrackInfo *trkInf(NULL);
  for(Int_t itk = 0; itk < nTracks; itk++){
    trkInf = dynamic_cast<AliTRDtrackInfo *>(fTracks->At(indices[itk]));
    if(trkInf->GetNTrackRefs()){
      // use Monte-Carlo Information for Momentum and PID
      trkmom = trkInf->GetTrackRef(0)->P();
      trkpdg = trkInf->GetPDG();
    }else{
      // Use TPC Momentum
      trkmom = trkInf->GetTrack()->P();
    }

    const Char_t *cmode(NULL);
    switch(mode){
      case kAccept:
        (dynamic_cast<TProfile *>(fContainer->At(kEfficiencyHistogram)))->Fill(trkmom, 1);
        (dynamic_cast<TProfile *>(fContainer->At(kContaminationHistogram)))->Fill(trkmom, 0);
        cmode="ACCEPT";
        break;
      case kMiss:
        (dynamic_cast<TProfile *>(fContainer->At(kEfficiencyHistogram)))->Fill(trkmom, 0);
        (dynamic_cast<TProfile *>(fContainer->At(kContaminationHistogram)))->Fill(trkmom, 0);
        cmode="MISS";
        break;
      case kFake:
        (dynamic_cast<TProfile *>(fContainer->At(kContaminationHistogram)))->Fill(trkmom, 1);
        cmode="FAKE";
        break;
    }
    AliDebug(3, Form(" track[%d] MC[%d] Mode[%s]", indices[itk], trkInf->GetLabel(), cmode));

    // Fill species histogram
    Int_t idxSpec = AliTRDpidUtil::Pdg2Pid(TMath::Abs(trkpdg));
    Int_t sign = dbPDG->GetParticle(trkpdg)->Charge() > 0. ? 1 : 0;
    //printf("[%d]%s pdg[%d] sign[%d]\n", idxSpec, AliPID::ParticleName(idxSpec), trkpdg, sign);
    if(idxSpec < 0) idxSpec = AliPID::kSPECIES;
    (dynamic_cast<TProfile *>(fContainer->At(kEfficiencySpeciesHistogram + idxSpec*2+sign)))->Fill(trkmom, mode==kAccept?1:0);
  }
}

//_____________________________________________________________________________
void AliTRDefficiencyMC::FillStreamTrackWOMC(AliTRDtrackInfo * const trkInf){
  // fill debug stream
  // we want to know:
  //  1. The event number
  //  2. The track label
  //  3. The TRD track label
  //  4. The frequency of the TRD Label
  //  5. Momentum from TPC (NO STAND ALONE TRACK)
  //  6. TPC Phi angle
  //  7. the TRD track
  //  8. Monte Carlo PID
  //  9. We check the Labels of the TRD track according to them we search the maching Monte-Carlo track.
  //     From the matching Monte-Carlo track we store trackRefs, phi and momentum
  // 10. We may also want to keep the kink index
  Double_t mom = trkInf->GetESDinfo()->GetOuterParam()->P();
  Int_t event = (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry();
  Int_t label = trkInf->GetLabel();
  Int_t kinkIndex = trkInf->GetKinkIndex();
  Int_t pdg = trkInf->GetPDG();
  Double_t phiTPC = trkInf->GetESDinfo()->GetOuterParam()->Phi();
  Int_t labelsTRD[180];	// Container for the cluster labels
  Int_t sortlabels[360];	// Cluster Labels sorted according their occurancy
  AliTRDseedV1 *tracklet(NULL);
  AliTRDcluster *c(NULL);
  Int_t nclusters(0);
  AliTRDtrackV1 *trackTRD = trkInf->GetTrack();
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    tracklet = trackTRD->GetTracklet(il);
    if(!tracklet) continue;
    tracklet->ResetClusterIter();
    c = NULL;
    while((c = tracklet->NextCluster())) labelsTRD[nclusters++] = c->GetLabel(0);
  }
  // Determine Label and Frequency
  AliMathBase::Freq(nclusters, const_cast<const Int_t *>(&labelsTRD[0]), &sortlabels[0], kTRUE);
  Int_t labelTRD = sortlabels[0];
  Int_t freqTRD = sortlabels[1];
  // find the track info object matching to the TRD track
  AliTRDtrackInfo *realtrack = 0;
  TObjArrayIter rtiter(fTracks);
  while((realtrack = (AliTRDtrackInfo *)rtiter())){
    if(realtrack->GetLabel() != labelTRD) continue;
    break;
  }
  TClonesArray trackRefs("AliTrackReference");
  Int_t realPdg = -1;
  Double_t realP = 0.;
  Double_t realPhi = 0.;
  if(realtrack){
    // pack the track references into the trackRefsContainer
    for(Int_t iref = 0; iref < realtrack->GetNTrackRefs(); iref++){
    new(trackRefs[iref])AliTrackReference(*(realtrack->GetTrackRef(iref)));
    }
    realPdg = realtrack->GetPDG();
    if(realtrack->GetNTrackRefs()){
      realP = realtrack->GetTrackRef(0)->P();
      realPhi = realtrack->GetTrackRef(0)->Phi();
    }
  }
  (*DebugStream()) << "EffMCfake"
    << "Event="	<< event
    << "Label=" << label
    << "labelTRD=" << labelTRD
    << "FreqTRDlabel=" << freqTRD
    << "TPCp="	<< mom
    << "phiTPC=" << phiTPC
    << "trackTRD=" << trackTRD
    << "PDG="	<< pdg
    << "TrackRefs=" << &trackRefs
    << "RealPDG=" << realPdg
    << "RealP="	<< realP
    << "RealPhi" << realPhi
    << "KinkIndex=" << kinkIndex
    << "\n";
}

//_____________________________________________________________________________
Int_t AliTRDefficiencyMC::IsRegistered(AliTRDtrackInfo * const trkInf, Int_t *indices, Int_t nTracks){
  //
  // Checks if track is registered in a given mode
  //

  Int_t label(trkInf->GetLabel());
  for(Int_t il(nTracks); il--;){
    if((dynamic_cast<AliTRDtrackInfo *>(fTracks->At(indices[il])))->GetLabel() == label) return il;
  }
  return -1;
}

