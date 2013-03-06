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
//
// static dNdPt helper functions
//
// basic functionality to select events and tracks 
// for dNdPt analysis
//
// Origin: Jan Fiete Grosse-Oetringhaus
// Modified and Extended: Jacek Otwinowski 19/11/2009
// last change: 2013-02-05 by M.Knichel
// 

#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <AliHeader.h>
#include <AliStack.h>
#include <AliLog.h>
#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliESDVertex.h>
#include <AliVertexerTracks.h>
#include <AliMathBase.h>
#include <AliESDtrackCuts.h>
#include <AliTracker.h>
#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include "AlidNdPtHelper.h"

//____________________________________________________________________
ClassImp(AlidNdPtHelper)

//____________________________________________________________________
const AliESDVertex* AlidNdPtHelper::GetVertex(AliESDEvent* const aEsd, const AlidNdPtEventCuts *const evtCuts, const AlidNdPtAcceptanceCuts *const accCuts, const AliESDtrackCuts *const trackCuts, AnalysisMode analysisMode, Bool_t debug, Bool_t bRedoTPC, Bool_t bUseMeanVertex)
{
  // Get the vertex from the ESD and returns it if the vertex is valid
  //
  // Second argument decides which vertex is used (this selects
  // also the quality criteria that are applied)

  if(!aEsd) 
  { 
    ::Error("AlidNdPtHelper::GetVertex()","esd event is NULL");
    return NULL;  
  }
 
  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetVertex()","cuts not available");
    return NULL;  
  }

  const AliESDVertex* vertex = 0;
  AliESDVertex *initVertex = 0;
  if (analysisMode == kSPD || 
      analysisMode == kTPCSPDvtx || analysisMode == kTPCSPDvtxUpdate || analysisMode == kTPCITSHybrid)
  {
    vertex = aEsd->GetPrimaryVertexSPD();
    if (debug)
      Printf("AlidNdPtHelper::GetVertex: Returning SPD vertex");
  }  
  else if (analysisMode == kTPCITS  || analysisMode == kTPCTrackSPDvtx || analysisMode == kTPCTrackSPDvtxUpdate || 
           analysisMode == kTPCITSHybridTrackSPDvtx || analysisMode == kTPCITSHybridTrackSPDvtxDCArPt || 
	   analysisMode == kITSStandAloneTrackSPDvtx ||  analysisMode == kITSStandAloneTPCTrackSPDvtx)
  {
    vertex = aEsd->GetPrimaryVertexTracks();
    if(!vertex) return NULL;
    if(vertex->GetNContributors()<1) {
      // SPD vertex
      vertex = aEsd->GetPrimaryVertexSPD();
    }
  }
  else if (analysisMode == kTPC)
  {
    if(bRedoTPC) {

      Double_t kBz = aEsd->GetMagneticField();
      AliVertexerTracks vertexer(kBz);

      if(bUseMeanVertex) {
	 Double_t pos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
	 Double_t err[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
	 initVertex = new AliESDVertex(pos,err);
	 vertexer.SetVtxStart(initVertex);
	 vertexer.SetConstraintOn();
      }

      Double_t maxDCAr = accCuts->GetMaxDCAr();
      Double_t maxDCAz = accCuts->GetMaxDCAz();
      Int_t minTPCClust = trackCuts->GetMinNClusterTPC();

      //vertexer.SetTPCMode(Double_t dcacut=0.1, Double_t dcacutIter0=1.0, Double_t maxd0z0=5.0, Int_t minCls=10, Int_t mintrks=1, Double_t nsigma=3., Double_t mindetfitter=0.1, Double_t maxtgl=1.5, Double_t fidR=3., Double_t fidZ=30., Int_t finderAlgo=1, Int_t finderAlgoIter0=4);
      vertexer.SetTPCMode(0.1,1.0,5.0,minTPCClust,1,3.,0.1,2.0,maxDCAr,maxDCAz,1,4);

      // TPC track preselection
      Int_t ntracks = aEsd->GetNumberOfTracks();
      TObjArray array(ntracks);
      UShort_t *id = new UShort_t[ntracks];


      Int_t count=0;
      for (Int_t i=0;i <ntracks; i++) {
        AliESDtrack *t = aEsd->GetTrack(i);
        if (!t) continue;
        if (t->Charge() == 0) continue;
        if (!t->GetTPCInnerParam()) continue;
        if (t->GetTPCNcls()<vertexer.GetMinClusters()) continue;
        AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
	if(tpcTrack) { 
	  array.AddLast(tpcTrack);
	  id[count] = (UShort_t)t->GetID();
	  count++;
	}
      } 
      AliESDVertex *vTPC = vertexer.VertexForSelectedTracks(&array,id, kTRUE, kTRUE, bUseMeanVertex);
      if(!vTPC) { 
        delete [] id; id=NULL;
        return 0;
      }
      
      // set recreated TPC vertex
      aEsd->SetPrimaryVertexTPC(vTPC);

      for (Int_t i=0; i<aEsd->GetNumberOfTracks(); i++) {
	AliESDtrack *t = aEsd->GetTrack(i);
	if(!t) continue;

        Double_t x[3]; t->GetXYZ(x);
        Double_t b[3]; AliTracker::GetBxByBz(x,b);
	t->RelateToVertexTPCBxByBz(vTPC, b, kVeryBig);
      }
      
      delete vTPC;
      array.Delete();
      delete [] id; id=NULL;

    }
    vertex = aEsd->GetPrimaryVertexTPC();
    if (debug)
     Printf("AlidNdPtHelper::GetVertex: Returning vertex from tracks");
   }
   else
     Printf("AlidNdPtHelper::GetVertex: ERROR: Invalid second argument %d", analysisMode);

    if (!vertex) {
     if (debug)
      Printf("AlidNdPtHelper::GetVertex: No vertex found in ESD");
      return 0;
    }

  if (debug)
  {
    Printf("AlidNdPtHelper::GetVertex: Returning valid vertex: %s", vertex->GetTitle());
    vertex->Print();
  }
  
  if(initVertex) delete initVertex; initVertex=NULL;
  return vertex;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::TestRecVertex(const AliESDVertex* vertex, const AliESDVertex* vertexSPD, AnalysisMode analysisMode, Bool_t debug)
{
  // Checks if a vertex meets the needed quality criteria
  if(!vertex) return kFALSE;
  if(!vertex->GetStatus()) return kFALSE;

  Float_t requiredZResolution = -1;
  if (analysisMode == kSPD || analysisMode == kTPCITS || 
      analysisMode == kTPCSPDvtx || analysisMode == kTPCSPDvtxUpdate || analysisMode == kTPCITSHybrid ||
      analysisMode == kTPCTrackSPDvtx || analysisMode == kTPCTrackSPDvtxUpdate || analysisMode == kTPCITSHybridTrackSPDvtx || analysisMode == kTPCITSHybridTrackSPDvtxDCArPt
	   || analysisMode == kITSStandAloneTrackSPDvtx ||  analysisMode == kITSStandAloneTPCTrackSPDvtx)
  {
    requiredZResolution = 1000;
  }
  else if (analysisMode == kTPC)
    requiredZResolution = 10.;

  // check resolution
  Double_t zRes = vertex->GetZRes();

  if (zRes > requiredZResolution) {
    if (debug)
      Printf("AlidNdPtHelper::TestVertex: Resolution too poor %f (required: %f", zRes, requiredZResolution);
    return kFALSE;
  }

  // always check for SPD vertex
  if(!vertexSPD) return kFALSE;
  if(!vertexSPD->GetStatus()) return kFALSE;
  if (vertexSPD->IsFromVertexerZ())
  {
    if (vertexSPD->GetDispersion() > 0.02) 
    {
      if (debug)
        Printf("AliPWG0Helper::TestVertex: Delta Phi too large in Vertexer Z: %f (required: %f", vertex->GetDispersion(), 0.02);
      return kFALSE;
    }
  }


  return kTRUE;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsGoodImpPar(const AliESDtrack *const track)
{
//
// check whether particle has good DCAr(Pt) impact
// parameter. Only for TPC+ITS tracks (7*sigma cut)
// Origin: Andrea Dainese
//

Float_t d0z0[2],covd0z0[3];
track->GetImpactParameters(d0z0,covd0z0);
Float_t sigma= 0.0050+0.0060/TMath::Power(track->Pt(),0.9);
Float_t d0max = 7.*sigma;
if(TMath::Abs(d0z0[0]) < d0max) return kTRUE;

return kFALSE;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsPrimaryParticle(AliStack* const stack, Int_t idx, ParticleMode particleMode)
{
//
// check primary particles 
// depending on the particle mode
//
  if(!stack) return kFALSE;

  TParticle* particle = stack->Particle(idx);
  if (!particle) return  kFALSE;

  // only charged particles
  Double_t charge = particle->GetPDG()->Charge()/3.;
  if (TMath::Abs(charge) < 0.001) return kFALSE;

  Int_t pdg = TMath::Abs(particle->GetPdgCode());

  // physical primary
  Bool_t prim = stack->IsPhysicalPrimary(idx);

  if(particleMode==kMCPion) {
    if(prim && pdg==kPiPlus) return kTRUE;
    else return kFALSE;
  } 

  if (particleMode==kMCKaon) {
    if(prim && pdg==kKPlus) return kTRUE;
    else return kFALSE;
  }
    
  if (particleMode==kMCProton) {
    if(prim && pdg==kProton) return kTRUE;
    else return kFALSE;
  }

  if(particleMode==kMCRest) {
    if(prim && pdg!=kPiPlus && pdg!=kKPlus && pdg!=kProton) return kTRUE;
    else return kFALSE;
  }

return prim;
}

//____________________________________________________________________
Bool_t AlidNdPtHelper::IsCosmicTrack(AliESDtrack *const track1, AliESDtrack *const track2)
{
//
// check cosmic tracks
//
  if(!track1) return kFALSE;
  if(!track2) return kFALSE;

  //
  // cosmic tracks in TPC
  //
  //if( TMath::Abs( track1->Theta() - track2->Theta() ) < 0.004  && 
  //  ((TMath::Abs(track1->Phi()-track2->Phi()) - TMath::Pi() )<0.004) && 
  //if( track1->Pt() > 4.0 && (TMath::Abs(track1->Phi()-track2->Phi())-TMath::Pi())<0.1 
  //    && (TMath::Abs(track1->Eta()+track2->Eta())-0.1) < 0.0 && (track1->Charge()+track2->Charge()) == 0)

  //Float_t scaleF= 6.0;
  if ( track1->Pt() > 4 && track2->Pt() > 4 && 
       //(TMath::Abs(track1->GetSnp()-track2->GetSnp())-2.) < scaleF * TMath::Sqrt(track1->GetSigmaSnp2()+track2->GetSigmaSnp2()) &&
       //TMath::Abs(track1->GetTgl()-track2->GetTgl())   < scaleF * TMath::Sqrt(track1->GetSigmaTgl2()+track2->GetSigmaTgl2()) &&
       //TMath::Abs(track1->OneOverPt()-track2->OneOverPt()) < scaleF * TMath::Sqrt(track1->GetSigma1Pt2()+track2->GetSigma1Pt2()) && 
       (track1->Charge()+track2->Charge()) == 0 && 
       track1->Eta()*track2->Eta()<0.0 && TMath::Abs(track1->Eta()+track2->Eta())<0.03 &&
       TMath::Abs(TMath::Abs(track1->Phi()-track2->Phi())-TMath::Pi())<0.1
     )  
  {
    printf("COSMIC  candidate \n");

    printf("track1->Pt() %f, track1->Theta() %f, track1->Eta() %f, track1->Phi() %f, track1->Charge() %d  \n", track1->Pt(), track1->Theta(), track1->Eta(), track1->Phi(), track1->Charge());
    printf("track2->Pt() %f, track2->Theta() %f, track2->Eta() %f, track2->Phi() %f, track2->Charge() %d  \n", track2->Pt(), track2->Theta(), track2->Eta(), track2->Phi(), track2->Charge());
    printf("dtheta %f, deta %f, dphi %f, dq %d  \n", track1->Theta()-track2->Theta(),  track1->Eta()-track2->Eta(), track1->Phi()-track2->Phi(), track1->Charge()+track2->Charge()); 
    printf("dsphi %f, errsphi %f, dtanl %f, errtanl %f  \n", TMath::Abs(track1->GetSnp()-track2->GetSnp()), TMath::Sqrt(track1->GetSigmaSnp2()+track2->GetSigmaSnp2()), TMath::Abs(track1->GetTgl()-track2->GetTgl()), TMath::Sqrt(track1->GetSigmaTgl2()+track2->GetSigmaTgl2())); 
    return kTRUE;
  }
     
return kFALSE; 
}

//____________________________________________________________________
void AlidNdPtHelper::PrintConf(AnalysisMode analysisMode, AliTriggerAnalysis::Trigger trigger)
{
  //
  // Prints the given configuration
  //

  TString str(">>>> Running with ");

  switch (analysisMode)
  {
    case kInvalid: str += "invalid setting"; break;
    case kSPD : str += "SPD-only"; break;
    case kTPC : str += "TPC-only"; break;
    case kTPCITS : str += "Global tracking"; break;
    case kTPCSPDvtx : str += "TPC tracking + SPD event vertex"; break;
    case kTPCSPDvtxUpdate : str += "TPC tracks updated with SPD event vertex point"; break;
    case kTPCTrackSPDvtx : str += "TPC tracking + Tracks event vertex or SPD event vertex"; break;
    case kTPCTrackSPDvtxUpdate : str += "TPC tracks updated with Track or SPD event vertex point"; break;
    case kTPCITSHybrid : str += "TPC tracking + ITS refit + >1 SPD cluster"; break;
    case kTPCITSHybridTrackSPDvtx : str += "TPC tracking + ITS refit + >1 SPD cluster + Tracks event vertex or SPD event vertex"; break;
    case kTPCITSHybridTrackSPDvtxDCArPt : str += "TPC tracking + ITS refit + >1 SPD cluster + Tracks event vertex or SPD event vertex + DCAr(pt)"; break;
    case kITSStandAloneTrackSPDvtx : str += "ITS standalone + Tracks event vertex or SPD event vertex + DCAr(pt)"; break;
    case kITSStandAloneTPCTrackSPDvtx : str += "kITSStandAloneTPCTrackSPDvtx + TPC stand alone track  + Tracks event vertex or SPD event vertex + DCAr(pt)"; break;
    case kMCRec : str += "TPC tracking + Replace rec. with MC values"; break;
  }
  str += " and trigger ";

  str += AliTriggerAnalysis::GetTriggerName(trigger);

  str += " <<<<";

  Printf("%s", str.Data());
}

//____________________________________________________________________
Int_t AlidNdPtHelper::ConvertPdgToPid(const TParticle *const particle) {
//
// Convert Abs(pdg) to pid 
// (0 - e, 1 - muons, 2 - pions, 3 - kaons, 4 - protons, 5 -all rest)
//
Int_t pid=-1;

  if (TMath::Abs(particle->GetPdgCode()) == kElectron)         { pid = 0; }
  else if (TMath::Abs(particle->GetPdgCode()) == kMuonMinus) { pid = 1; }
  else if (TMath::Abs(particle->GetPdgCode()) == kPiPlus)    { pid = 2; }
  else if (TMath::Abs(particle->GetPdgCode()) == kKPlus)     { pid = 3; }
  else if (TMath::Abs(particle->GetPdgCode()) == kProton)    { pid = 4; }
  else                                                       { pid = 5; }

return pid;
}

//_____________________________________________________________________________
TH1F* AlidNdPtHelper::CreateResHisto(TH2F* const hRes2, TH1F **phMean, Int_t integ,  Bool_t drawBinFits, Int_t minHistEntries)
{
//
// Create mean and resolution 
// histograms
//
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  //Bool_t overflowBinFits = kFALSE;
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  TF1 * fitFunc = new TF1("G","[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  //Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nPads = nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  //Int_t dBin = ((overflowBinFits) ? 1 : 0);
  Int_t dBin = 0;
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    Int_t bin0=TMath::Max(bin-integ,0);
    Int_t bin1=TMath::Min(bin+integ,nBins);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin0, bin1);
    //    
    if (hBin->GetEntries() > minHistEntries) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	snprintf(name,256, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	snprintf(name,256, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	snprintf(name,256, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  *phMean = hMean;
  return hRes;
}

//_____________________________________________________________________________
TH1F* AlidNdPtHelper::MakeResol(TH2F * his, Int_t integ, Bool_t type, Bool_t drawBins, Int_t minHistEntries){
// Create resolution histograms
  
     TH1F *hisr=0, *hism=0;
     if (!gPad) new TCanvas;
         hisr = CreateResHisto(his,&hism,integ,drawBins,minHistEntries);
         if (type) return hism;
         else return hisr;

return hisr;	 
}

//_____________________________________________________________________________
TObjArray* AlidNdPtHelper::GetAllChargedTracks(AliESDEvent *esdEvent, AnalysisMode analysisMode)
{
  //
  // all charged TPC particles 
  //
  TObjArray *allTracks = new TObjArray();
  if(!allTracks) return allTracks;

  AliESDtrack *track=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    if(analysisMode == AlidNdPtHelper::kTPC) { 
      //
      // track must be deleted by user 
      // esd track parameters are replaced by TPCinner
      //
      track = AliESDtrackCuts::GetTPCOnlyTrack(esdEvent,iTrack);
      if(!track) continue;
    } 
    else if (analysisMode == AlidNdPtHelper::kTPCSPDvtx || analysisMode == AlidNdPtHelper::kTPCSPDvtxUpdate)
    {
      //
      // track must be deleted by the user 
      // esd track parameters are replaced by TPCinner
      //
      track = AlidNdPtHelper::GetTPCOnlyTrackSPDvtx(esdEvent,iTrack,kFALSE);
      if(!track) continue;
    }
    else if (analysisMode == AlidNdPtHelper::kTPCTrackSPDvtx || analysisMode == AlidNdPtHelper::kTPCTrackSPDvtxUpdate)
    {
      //
      // track must be deleted by the user 
      // esd track parameters are replaced by TPCinner
      //
      track = AlidNdPtHelper::GetTPCOnlyTrackTrackSPDvtx(esdEvent,iTrack,kFALSE);
      if(!track) continue;
    }
    else if(analysisMode == AlidNdPtHelper::kTPCITSHybrid )
    {
      track = AlidNdPtHelper::GetTrackSPDvtx(esdEvent,iTrack,kFALSE);
    }
    else if(analysisMode == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || analysisMode == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt || analysisMode == AlidNdPtHelper::kITSStandAloneTrackSPDvtx || analysisMode ==AlidNdPtHelper::kITSStandAloneTPCTrackSPDvtx)
    {
      track = AlidNdPtHelper::GetTrackTrackSPDvtx(esdEvent,iTrack,kFALSE);
    }
    else 
    {
      track = esdEvent->GetTrack(iTrack);
    }

    if(!track) continue;

    if(track->Charge()==0) { 
      if(analysisMode == AlidNdPtHelper::kTPC || 
         analysisMode == AlidNdPtHelper::kTPCSPDvtx || analysisMode == AlidNdPtHelper::kTPCTrackSPDvtx  ||
         analysisMode == AlidNdPtHelper::kTPCSPDvtxUpdate || analysisMode == AlidNdPtHelper::kTPCTrackSPDvtxUpdate) 
      {
        delete track; continue; 
      } else {
        continue;
      } 
    }

    allTracks->Add(track);
  }

  if(analysisMode == AlidNdPtHelper::kTPC || 
     analysisMode == AlidNdPtHelper::kTPCSPDvtx || analysisMode == AlidNdPtHelper::kTPCTrackSPDvtx || 
     analysisMode == AlidNdPtHelper::kTPCSPDvtxUpdate || analysisMode == AlidNdPtHelper::kTPCTrackSPDvtxUpdate) {
     
     allTracks->SetOwner(kTRUE);
  }

return allTracks;
}

//_____________________________________________________________________________
AliESDtrack *AlidNdPtHelper::GetTPCOnlyTrackSPDvtx(const AliESDEvent* esdEvent, Int_t iTrack, Bool_t bUpdate)
{
//
// Create ESD tracks from TPCinner parameters.
// Propagte to DCA to SPD vertex.
// Update using SPD vertex point (parameter)
//
// It is user responsibility to delete these tracks
//

  if (!esdEvent) return NULL;
  if (!esdEvent->GetPrimaryVertexSPD() ) { return NULL; }
  if (!esdEvent->GetPrimaryVertexSPD()->GetStatus() ) { return  NULL; }
   
  // 
  AliESDtrack* track = esdEvent->GetTrack(iTrack);
  if (!track)
    return NULL;

  Bool_t isOK = kFALSE;
  Double_t x[3]; track->GetXYZ(x);
  Double_t b[3]; AliTracker::GetBxByBz(x,b);

  // create new ESD track
  AliESDtrack *tpcTrack = new AliESDtrack();
 
  // relate TPC-only tracks (TPCinner) to SPD vertex
  AliExternalTrackParam cParam;
  if(bUpdate) {  
    isOK = track->RelateToVertexTPCBxByBz(esdEvent->GetPrimaryVertexSPD(),b,kVeryBig,&cParam);
    track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());

    // reject fake tracks
    if(track->Pt() > 10000.)  {
      ::Error("Exclude no physical tracks","pt>10000. GeV");
      delete tpcTrack; 
      return NULL;
    }
  }
  else {
    isOK = track->RelateToVertexTPCBxByBz(esdEvent->GetPrimaryVertexSPD(), b, kVeryBig);
  }

  // only true if we have a tpc track
  if (!track->FillTPCOnlyTrack(*tpcTrack))
  {
    delete tpcTrack;
    return NULL;
  }
  
  if(!isOK) return NULL;

return tpcTrack;
} 

//_____________________________________________________________________________
AliESDtrack *AlidNdPtHelper::GetTPCOnlyTrackTrackSPDvtx(const AliESDEvent* esdEvent, Int_t iTrack, Bool_t bUpdate)
{
//
// Create ESD tracks from TPCinner parameters.
// Propagte to DCA to Track or SPD vertex.
// Update using SPD vertex point (parameter)
//
// It is user responsibility to delete these tracks
//
  if (!esdEvent) return NULL;
  const AliESDVertex *vertex = esdEvent->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = esdEvent->GetPrimaryVertexSPD();
  }
  if(!vertex) return NULL;
 
  // 
  AliESDtrack* track = esdEvent->GetTrack(iTrack);
  if (!track)
    return NULL;

  Bool_t isOK = kFALSE;
  Double_t x[3]; track->GetXYZ(x);
  Double_t b[3]; AliTracker::GetBxByBz(x,b);

  // create new ESD track
  AliESDtrack *tpcTrack = new AliESDtrack();
 
  // relate TPC-only tracks (TPCinner) to SPD vertex
  AliExternalTrackParam cParam;
  if(bUpdate) {  
    isOK = track->RelateToVertexTPCBxByBz(vertex,b,kVeryBig,&cParam);
    track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());

    // reject fake tracks
    if(track->Pt() > 10000.)  {
      ::Error("Exclude no physical tracks","pt>10000. GeV");
      delete tpcTrack; 
      return NULL;
    }
  }
  else {
    isOK = track->RelateToVertexTPCBxByBz(vertex, b, kVeryBig);
  }

  // only true if we have a tpc track
  if (!track->FillTPCOnlyTrack(*tpcTrack))
  {
    delete tpcTrack;
    return NULL;
  }
  
  if(!isOK) return NULL;

return tpcTrack;
} 

//_____________________________________________________________________________
AliESDtrack *AlidNdPtHelper::GetTrackSPDvtx(const AliESDEvent* esdEvent, Int_t iTrack, Bool_t bUpdate)
{
//
// Propagte track to DCA to SPD vertex.
// Update using SPD vertex point (parameter)
//
  if (!esdEvent) return NULL;
  if (!esdEvent->GetPrimaryVertexSPD() ) { return NULL; }
  if (!esdEvent->GetPrimaryVertexSPD()->GetStatus() ) { return  NULL; }
   
  // 
  AliESDtrack* track = esdEvent->GetTrack(iTrack);
  if (!track)
    return NULL;

  Bool_t isOK = kFALSE;
  Double_t x[3]; track->GetXYZ(x);
  Double_t b[3]; AliTracker::GetBxByBz(x,b);

  // relate tracks to SPD vertex
  AliExternalTrackParam cParam;
  if(bUpdate) {  
    isOK = track->RelateToVertexBxByBz(esdEvent->GetPrimaryVertexSPD(),b,kVeryBig,&cParam);
    track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());

    // reject fake tracks
    if(track->Pt() > 10000.)  {
      ::Error("Exclude no physical tracks","pt>10000. GeV");
      return NULL;
    }
  }
  else {
    isOK = track->RelateToVertexBxByBz(esdEvent->GetPrimaryVertexSPD(), b, kVeryBig);
  }
 
  if(!isOK) return NULL;

return track;
} 

//_____________________________________________________________________________
AliESDtrack *AlidNdPtHelper::GetTrackTrackSPDvtx(const AliESDEvent* esdEvent, Int_t iTrack, Bool_t bUpdate)
{
//
// Propagte track to DCA to Track or SPD vertex.
// Update using SPD vertex point (parameter)
//
  if (!esdEvent) return NULL;

  const AliESDVertex *vertex = esdEvent->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = esdEvent->GetPrimaryVertexSPD();
  }
  if(!vertex) return NULL;

  // 
  AliESDtrack* track = esdEvent->GetTrack(iTrack);
  if (!track)
    return NULL;

  Bool_t isOK = kFALSE;
  Double_t x[3]; track->GetXYZ(x);
  Double_t b[3]; AliTracker::GetBxByBz(x,b);

  // relate tracks to SPD vertex
  AliExternalTrackParam cParam;
  if(bUpdate) {  
    isOK = track->RelateToVertexBxByBz(vertex,b,kVeryBig,&cParam);
    track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());

    // reject fake tracks
    if(track->Pt() > 10000.)  {
      ::Error("Exclude no physical tracks","pt>10000. GeV");
      return NULL;
    }
  }
  else {
    isOK = track->RelateToVertexBxByBz(vertex, b, kVeryBig);
  }
 
  if(!isOK) return NULL;

return track;
} 

//_____________________________________________________________________________
Bool_t AlidNdPtHelper::SelectEvent(const AliESDEvent* const esdEvent, AliESDtrackCuts* const esdTrackCuts) {
// select events with at least
// one reconstructed primary track in acceptance
// pT>0.5 GeV/c, |eta|<0.8 for cross section studies

if(!esdEvent) return kFALSE;
if(!esdTrackCuts) return kFALSE;

  AliESDtrack *track=0;
  Int_t count = 0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
  { 
    track = esdEvent->GetTrack(iTrack);
    if(!track) continue;
    if(track->Charge()==0) continue;
    if(!esdTrackCuts->AcceptTrack(track)) continue;
    if(track->Pt() < 0.5) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;

    count++;
  }

  if(count > 0) return kTRUE;
  else return kFALSE;

return kFALSE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtHelper::SelectMCEvent(AliMCEvent* const mcEvent) {
//
// select events with at least
// one prompt (MC primary) track in acceptance
// pT>0.5 GeV/c, |eta|<0.8 for cross section studies
//

if(!mcEvent) return kFALSE;
AliStack* stack = mcEvent->Stack(); 
if(!stack) return kFALSE;

  Int_t count = 0;
  for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;

    // only charged particles
    if(!particle->GetPDG()) continue;
    Double_t charge = particle->GetPDG()->Charge()/3.;
    if(charge == 0) continue;

    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    if(particle->Pt() < 0.5) continue;
    if(TMath::Abs(particle->Eta()) > 0.8) continue;

    count++;
  }

  if(count > 0) return kTRUE;
  else return kFALSE;

return kFALSE;
}

//_____________________________________________________________________________
Int_t AlidNdPtHelper::GetTPCMBTrackMult(const AliESDEvent *const esdEvent,const AlidNdPtEventCuts *const evtCuts, const AlidNdPtAcceptanceCuts *const accCuts,const  AliESDtrackCuts *const trackCuts)
{
  //
  // get MB event track multiplicity
  //
  if(!esdEvent) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBTrackMult()","esd event is NULL");
    return 0;  
  }
 
  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBTrackMult()","cuts not available");
    return 0;  
  }

  //
  Double_t pos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
  Double_t err[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
  AliESDVertex vtx0(pos,err);

  //
  Float_t maxDCAr = accCuts->GetMaxDCAr();
  Float_t maxDCAz = accCuts->GetMaxDCAz();
  Float_t minTPCClust = trackCuts->GetMinNClusterTPC();
  //
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  Double_t dca[2],cov[3];
  Int_t mult=0;
  for (Int_t i=0;i <ntracks; i++){
    AliESDtrack *t = esdEvent->GetTrack(i);
    if (!t) continue;
    if (t->Charge() == 0) continue;
    if (!t->GetTPCInnerParam()) continue;
    if (t->GetTPCNcls()<minTPCClust) continue;
    //
    Double_t x[3]; t->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    const Double_t kMaxStep = 1;   //max step over the material

    AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
    if(!tpcTrack) return 0;

    if (!tpcTrack->PropagateToDCABxByBz(&vtx0,b,kMaxStep,dca,cov)) 
    {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }
    //
    if (TMath::Abs(dca[0])>maxDCAr || TMath::Abs(dca[1])>maxDCAz) {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }

    mult++;    

    if(tpcTrack) delete tpcTrack; 
  }

return mult;  
}

//_____________________________________________________________________________
Int_t AlidNdPtHelper::GetTPCMBPrimTrackMult(const AliESDEvent *const esdEvent, AliStack *const  stack, const AlidNdPtEventCuts *const evtCuts, const AlidNdPtAcceptanceCuts *const accCuts, const AliESDtrackCuts *const trackCuts)
{
  //
  // get MB primary event track multiplicity
  //
  if(!esdEvent) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBPrimTrackMult()","esd event is NULL");
    return 0;  
  }

  if(!stack) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBPrimTrackMult()","esd event is NULL");
    return 0;  
  }
 
  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetTPCMBPrimTrackMult()","cuts not available");
    return 0;  
  }

  //
  Double_t pos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
  Double_t err[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
  AliESDVertex vtx0(pos,err);

  //
  Float_t maxDCAr = accCuts->GetMaxDCAr();
  Float_t maxDCAz = accCuts->GetMaxDCAz();
  Float_t minTPCClust = trackCuts->GetMinNClusterTPC();

  //
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  Double_t dca[2],cov[3];
  Int_t mult=0;
  for (Int_t i=0;i <ntracks; i++){
    AliESDtrack *t = esdEvent->GetTrack(i);
    if (!t) continue;
    if (t->Charge() == 0) continue;
    if (!t->GetTPCInnerParam()) continue;
    if (t->GetTPCNcls()<minTPCClust) continue;
    //
    Double_t x[3]; t->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    const Double_t kMaxStep = 1;   //max step over the material

    AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
    if(!tpcTrack) return 0;

    if (!tpcTrack->PropagateToDCABxByBz(&vtx0,b,kMaxStep,dca,cov)) 
    {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }
    //
    if (TMath::Abs(dca[0])>maxDCAr || TMath::Abs(dca[1])>maxDCAz) {
      if(tpcTrack) delete tpcTrack; 
      continue;
    }

    Int_t label = TMath::Abs(t->GetLabel());
    TParticle *part = stack->Particle(label);
    if(!part) { 
      if(tpcTrack) delete tpcTrack; 
      continue;
    }
    if(!stack->IsPhysicalPrimary(label)) 
    { 
      if(tpcTrack) delete tpcTrack; 
      continue;
    }

    mult++;    

    if(tpcTrack) delete tpcTrack; 
  }

return mult;  
}


//_____________________________________________________________________________
Double_t AlidNdPtHelper::GetStrangenessCorrFactor(const Double_t pt)
{
// data driven correction factor for secondaries
// underestimated secondaries with strangeness in Pythia (A. Dainese)

//
// pt=0.17; fact=1
// pt=0.4; fact=1.07
// pt=0.6; fact=1.25
// pt>=1.2; fact=1.5
//

if (pt <= 0.17) return 1.0;
if (pt <= 0.4) return GetLinearInterpolationValue(0.17,1.0,0.4,1.07, pt);
if (pt <= 0.6) return GetLinearInterpolationValue(0.4,1.07,0.6,1.25, pt);
if (pt <= 1.2) return GetLinearInterpolationValue(0.6,1.25,1.2,1.5,  pt);
return 1.5;

}

//_____________________________________________________________________________
Double_t AlidNdPtHelper::GetStrangenessCorrFactorPbPb(const Double_t pt)
{
// data driven correction factor for secondaries (PbPb)

if (pt <= 0.25) return 1.0;
if (pt <= 0.5) return GetLinearInterpolationValue(0.25,1.0,0.5,1.4, pt);
if (pt <= 1.0) return GetLinearInterpolationValue(0.5,1.4,1.0,1.47, pt);
if (pt <= 2.0) return GetLinearInterpolationValue(1.0,1.47,2.0,1.56,  pt);
if (pt <= 5.0) return GetLinearInterpolationValue(2.0,1.56,5.0,1.67,  pt);
return 1.67;

}

//___________________________________________________________________________
Double_t AlidNdPtHelper::GetLinearInterpolationValue(const Double_t x1,const  Double_t y1,const  Double_t x2,const  Double_t y2, const Double_t pt)
{
//
// linear interpolation
//
  return ((y2-y1)/(x2-x1))*pt+(y2-(((y2-y1)/(x2-x1))*x2)); 
}

//_____________________________________________________________________________
Int_t AlidNdPtHelper::GetMCTrueTrackMult(AliMCEvent *const mcEvent, AlidNdPtEventCuts *const evtCuts, AlidNdPtAcceptanceCuts *const accCuts)
{
  //
  // calculate mc event true track multiplicity
  //
  if(!mcEvent) return 0;

  AliStack* stack = 0;
  Int_t mult = 0;

  // MC particle stack
  stack = mcEvent->Stack();
  if (!stack) return 0;

  //
  //printf("minZv %f, maxZv %f \n", evtCuts->GetMinZv(), evtCuts->GetMaxZv());
  //

  Bool_t isEventOK = evtCuts->AcceptMCEvent(mcEvent);
  if(!isEventOK) return 0; 

  Int_t nPart  = stack->GetNtrack();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
     TParticle* particle = stack->Particle(iMc);
     if (!particle)
     continue;

     // only charged particles
     if(!particle->GetPDG()) continue;
     Double_t charge = particle->GetPDG()->Charge()/3.;
     if (TMath::Abs(charge) < 0.001)
     continue;
      
     // physical primary
     Bool_t prim = stack->IsPhysicalPrimary(iMc);
     if(!prim) continue;

     // checked accepted without pt cut
     //if(accCuts->AcceptTrack(particle)) 
     if( particle->Eta() > accCuts->GetMinEta() && particle->Eta() < accCuts->GetMaxEta() ) 
     {
       mult++;
     }
  }

return mult;  
}

//_______________________________________________________________________
void  AlidNdPtHelper::PrintMCInfo(AliStack *const pStack,Int_t label)
{
// print information about particles in the stack

  if(!pStack)return;
  label = TMath::Abs(label);
  TParticle *part = pStack->Particle(label);
  Printf("########################");
  Printf("%s:%d %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,label,part->GetUniqueID(),part->GetPdgCode(),part->P());
  part->Print();
  TParticle* mother = part;
  Int_t imo = part->GetFirstMother();
  Int_t nprim = pStack->GetNprimary();

  while((imo >= nprim)) {
      mother =  pStack->Particle(imo);
      Printf("Mother %s:%d Label %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,imo,mother->GetUniqueID(),mother->GetPdgCode(),mother->P());
      mother->Print();
      imo =  mother->GetFirstMother();
 }

 Printf("########################");
}


//_____________________________________________________________________________
TH1* AlidNdPtHelper::GetContCorrHisto(TH1 *const hist) 
{
//
// get contamination histogram
//
 if(!hist) return 0;

 Int_t nbins = hist->GetNbinsX();
 TH1 *hCont = (TH1D *)hist->Clone();

 for(Int_t i=0; i<=nbins+1; i++) {
   Double_t binContent = hist->GetBinContent(i);
   Double_t binError = hist->GetBinError(i);

   hCont->SetBinContent(i,1.-binContent);
   hCont->SetBinError(i,binError);
 }

return hCont;
}


//_____________________________________________________________________________
TH1* AlidNdPtHelper::ScaleByBinWidth(TH1 *const hist) 
{
//
// scale by bin width
//
 if(!hist) return 0;

 TH1 *hScale = (TH1D *)hist->Clone();
 hScale->Scale(1.,"width");

return hScale;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::CalcRelativeDifference(const TH1 *const hist1, const TH1 *const hist2) 
{
//
// calculate rel. difference 
//

 if(!hist1) return 0;
 if(!hist2) return 0;

 TH1 *h1Clone = (TH1D *)hist1->Clone();
 h1Clone->Sumw2();

 // (rec-mc)/mc
 h1Clone->Add(hist2,-1);
 h1Clone->Divide(hist2);

return h1Clone;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::CalcRelativeDifferenceFun(const TH1 *const hist1, TF1 *const fun) 
{
//
// calculate rel. difference
// between histogram and function
//
 if(!hist1) return 0;
 if(!fun) return 0;

 TH1 *h1Clone = (TH1D *)hist1->Clone();
 h1Clone->Sumw2();

 // 
 h1Clone->Add(fun,-1);
 h1Clone->Divide(hist1);

return h1Clone;
}

//_____________________________________________________________________________
TH1* AlidNdPtHelper::NormalizeToEvent(const TH2 *const hist1, const TH1 *const hist2) 
{
// normalise to event for a given multiplicity bin
// return pt histogram 

 if(!hist1) return 0;
 if(!hist2) return 0;
 char name[256];

 Int_t nbinsX = hist1->GetNbinsX();
 //Int_t nbinsY = hist1->GetNbinsY();

 TH1D *histNorm = 0;
 for(Int_t i=0; i<=nbinsX+1; i++) {
   snprintf(name,256,"mom_%d",i);
   TH1D *hist = (TH1D*)hist1->ProjectionY(name,i+1,i+1);

   snprintf(name,256,"mom_norm");
   if(i==0) { 
     histNorm = (TH1D *)hist->Clone(name);
     histNorm->Reset();
   }

   Double_t nbEvents = hist2->GetBinContent(i);
   if(!nbEvents) { nbEvents = 1.; };

   hist->Scale(1./nbEvents);
   histNorm->Add(hist);
 }

return histNorm;
}

//_____________________________________________________________________________
//THnSparse* AlidNdPtHelper::GenerateCorrMatrix(THnSparse *const hist1, const THnSparse *const hist2, char *const name) {
THnSparse* AlidNdPtHelper::GenerateCorrMatrix(THnSparse *const hist1, const THnSparse *const hist2, const char *name) {
// generate correction matrix
if(!hist1 || !hist2) return 0; 

THnSparse *h =(THnSparse*)hist1->Clone(name);;
h->Divide(hist1,hist2,1,1,"B");

return h;
}

//_____________________________________________________________________________
//TH2* AlidNdPtHelper::GenerateCorrMatrix(TH2 *const hist1, TH2 *const hist2, char *const name) {
TH2* AlidNdPtHelper::GenerateCorrMatrix(TH2 *const hist1, TH2 *const hist2, const char *name) {
// generate correction matrix
if(!hist1 || !hist2) return 0; 

TH2D *h =(TH2D*)hist1->Clone(name);;
h->Divide(hist1,hist2,1,1,"B");

return h;
}

//_____________________________________________________________________________
//TH1* AlidNdPtHelper::GenerateCorrMatrix(TH1 *const hist1, TH1 *const hist2, char *const name) {
TH1* AlidNdPtHelper::GenerateCorrMatrix(TH1 *const hist1, TH1 *const hist2, const char* name) {
// generate correction matrix
if(!hist1 || !hist2) return 0; 

TH1D *h =(TH1D*)hist1->Clone(name);;
h->Divide(hist1,hist2,1,1,"B");

return h;
}

//_____________________________________________________________________________
//THnSparse* AlidNdPtHelper::GenerateContCorrMatrix(THnSparse *const hist1, THnSparse *const hist2, char *const name) {
THnSparse* AlidNdPtHelper::GenerateContCorrMatrix(THnSparse *const hist1, const THnSparse *const hist2, const char* name) {
// generate contamination correction matrix
if(!hist1 || !hist2) return 0; 

THnSparse *hist =  GenerateCorrMatrix(hist1, hist2, name);
if(!hist) return 0;

// only for non ZERO bins!!!!

Int_t* coord = new Int_t[hist->GetNdimensions()];
memset(coord, 0, sizeof(Int_t) * hist->GetNdimensions());

  for (Long64_t i = 0; i < hist->GetNbins(); ++i) {
    Double_t v = hist->GetBinContent(i, coord);
    hist->SetBinContent(coord, 1.0-v);
    //printf("v %f, hist->GetBinContent(i, coord) %f \n",v,hist->GetBinContent(i, coord));
    Double_t err = hist->GetBinError(coord);
    hist->SetBinError(coord, err);
  }

delete [] coord;

return hist;
}

//_____________________________________________________________________________
//TH2* AlidNdPtHelper::GenerateContCorrMatrix(TH2 *const hist1, TH2 *const hist2, char *const name) {
TH2* AlidNdPtHelper::GenerateContCorrMatrix(TH2 *const hist1, TH2 *const hist2, const char* name) {
// generate contamination correction matrix
if(!hist1 || !hist2) return 0; 

TH2 *hist = GenerateCorrMatrix(hist1, hist2, name);
if(!hist) return 0;

Int_t nBinsX = hist->GetNbinsX();
Int_t nBinsY = hist->GetNbinsY();

  for (Int_t i = 0; i < nBinsX+1; i++) {
  for (Int_t j = 0; j < nBinsY+1; j++) {
     Double_t cont = hist->GetBinContent(i,j);
     hist->SetBinContent(i,j,1.-cont);
     Double_t err = hist->GetBinError(i,j);
     hist->SetBinError(i,j,err);
  }
  }

return hist;
}

//_____________________________________________________________________________
//TH1* AlidNdPtHelper::GenerateContCorrMatrix(TH1 *const hist1, TH1 *const hist2, char *const name) {
TH1* AlidNdPtHelper::GenerateContCorrMatrix(TH1 *const hist1, TH1 *const hist2, const char* name) {
// generate contamination correction matrix
if(!hist1 || !hist2) return 0; 

TH1 *hist = GenerateCorrMatrix(hist1, hist2, name);
if(!hist) return 0;

Int_t nBinsX = hist->GetNbinsX();

  for (Int_t i = 0; i < nBinsX+1; i++) {
     Double_t cont = hist->GetBinContent(i);
     hist->SetBinContent(i,1.-cont);
     Double_t err = hist->GetBinError(i);
     hist->SetBinError(i,err);
  }

return hist;
}

//_____________________________________________________________________________
const AliESDVertex* AlidNdPtHelper::GetTPCVertexZ(const AliESDEvent* const esdEvent, const AlidNdPtEventCuts *const evtCuts, const AlidNdPtAcceptanceCuts *const accCuts, const AliESDtrackCuts *const trackCuts, Float_t fraction, Int_t ntracksMin){
  //
  // TPC Z vertexer
  //
  if(!esdEvent)
  { 
    ::Error("AlidNdPtHelper::GetTPCVertexZ()","cuts not available");
    return NULL;  
  }

  if(!evtCuts || !accCuts || !trackCuts) 
  { 
    ::Error("AlidNdPtHelper::GetTPCVertexZ()","cuts not available");
    return NULL;  
  }

  Double_t vtxpos[3]={evtCuts->GetMeanXv(),evtCuts->GetMeanYv(),evtCuts->GetMeanZv()};
  Double_t vtxsigma[3]={evtCuts->GetSigmaMeanXv(),evtCuts->GetSigmaMeanYv(),evtCuts->GetSigmaMeanZv()};
  AliESDVertex vtx0(vtxpos,vtxsigma);

  Double_t maxDCAr = accCuts->GetMaxDCAr();
  Double_t maxDCAz = accCuts->GetMaxDCAz();
  Int_t minTPCClust = trackCuts->GetMinNClusterTPC();

  //
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  TVectorD ztrack(ntracks);
  Double_t dca[2],cov[3];
  Int_t counter=0;
  for (Int_t i=0;i <ntracks; i++){
    AliESDtrack *t = esdEvent->GetTrack(i);
    if (!t) continue;
    if (!t->GetTPCInnerParam()) continue;
    if (t->GetTPCNcls()<minTPCClust) continue;
    //

    Double_t x[3]; t->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    const Double_t kMaxStep = 1;   //max step over the material

    AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(t->GetTPCInnerParam()));
    if(!tpcTrack) return 0;
    if (!tpcTrack->PropagateToDCABxByBz(&vtx0,b,kMaxStep,dca,cov)) continue;

    //
    if (TMath::Abs(dca[0])>maxDCAr) continue;
    //if (TMath::Sqrt(cov[0])>sigmaXYcut) continue;    
    if (TMath::Abs(tpcTrack->GetZ())>maxDCAz) continue;

    ztrack[counter]=tpcTrack->GetZ();
    counter++;    

    if(tpcTrack) delete tpcTrack;
  }

  //
  // Find LTM z position
  //
  Double_t mean=0, sigma=0;
  if (counter<ntracksMin) return 0;
  //
  Int_t nused = TMath::Nint(counter*fraction);
  if (nused==counter) nused=counter-1;  
  if (nused>1){
    AliMathBase::EvaluateUni(counter, ztrack.GetMatrixArray(), mean,sigma, TMath::Nint(counter*fraction));
    sigma/=TMath::Sqrt(nused);
  }else{
    mean  = TMath::Mean(counter, ztrack.GetMatrixArray());
    sigma = TMath::RMS(counter, ztrack.GetMatrixArray());
    sigma/=TMath::Sqrt(counter-1);
  }
  vtxpos[2]=mean;
  vtxsigma[2]=sigma;
  const AliESDVertex* vertex = new AliESDVertex(vtxpos, vtxsigma);
  return vertex;
}

//_____________________________________________________________________________
Int_t  AlidNdPtHelper::GetSPDMBTrackMult(const AliESDEvent* const esdEvent, Float_t deltaThetaCut, Float_t deltaPhiCut) 
{
  //
  // SPD track multiplicity
  //

  // get tracklets
  const AliMultiplicity* mult = esdEvent->GetMultiplicity();
  if (!mult)
     return 0;

  // get multiplicity from SPD tracklets
  Int_t inputCount = 0; 
  for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
  {
    //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

     Float_t phi = mult->GetPhi(i);
     if (phi < 0)
       phi += TMath::Pi() * 2;
     Float_t deltaPhi = mult->GetDeltaPhi(i);
     Float_t deltaTheta = mult->GetDeltaTheta(i);

     if (TMath::Abs(deltaPhi) > 1)
       printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

     if (deltaThetaCut > 0. && TMath::Abs(deltaTheta) > deltaThetaCut)
        continue;

     if (deltaPhiCut > 0. && TMath::Abs(deltaPhi) > deltaPhiCut)
        continue;
      
     ++inputCount;
  }

return inputCount;
}

//_____________________________________________________________________________
Int_t  AlidNdPtHelper::GetSPDMBPrimTrackMult(const AliESDEvent* const esdEvent, AliStack* const stack, Float_t deltaThetaCut, Float_t deltaPhiCut) 
{
  //
  // SPD track multiplicity
  //

  // get tracklets
  const AliMultiplicity* mult = esdEvent->GetMultiplicity();
  if (!mult)
     return 0;

  // get multiplicity from SPD tracklets
  Int_t inputCount = 0; 
  for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
  {
    //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));

     Float_t phi = mult->GetPhi(i);
     if (phi < 0)
       phi += TMath::Pi() * 2;
     Float_t deltaPhi = mult->GetDeltaPhi(i);
     Float_t deltaTheta = mult->GetDeltaTheta(i);

     if (TMath::Abs(deltaPhi) > 1)
       printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), deltaPhi);

     if (deltaThetaCut > 0. && TMath::Abs(deltaTheta) > deltaThetaCut)
        continue;

     if (deltaPhiCut > 0. && TMath::Abs(deltaPhi) > deltaPhiCut)
        continue;


     if (mult->GetLabel(i, 0) < 0 || mult->GetLabel(i, 0) != mult->GetLabel(i, 1) || 
         !stack->IsPhysicalPrimary(mult->GetLabel(i, 0)))
        continue;

      
     ++inputCount;
  }

return inputCount;
}

//_____________________________________________________________________________

THnSparse* AlidNdPtHelper::RebinTHnSparse(const THnSparse* hist1, THnSparse* hist2, const Char_t* newname, Option_t* option)
{
    THnSparse* htemp = 0;
    const THnSparse* hist = 0;
    TString opt = option;
    opt.ToLower();
    Bool_t calcErrors = kFALSE;
    Bool_t useRange = kFALSE;
    Bool_t overwrite = kFALSE;
    if (opt.Contains("e")) { calcErrors = kTRUE; } // calcluate correct errors (not implemented)
    if (opt.Contains("r")) { useRange = kTRUE; }   // use the axis range given in hist1
    if (opt.Contains("o")) { overwrite = kTRUE; }  // overwrite hist2 instead of creating a new one
    Int_t ndim  = hist1->GetNdimensions();
    if (ndim != hist2->GetNdimensions()) {
        printf("AlidNdPtHelper::RebinTHnSparse: ERROR: Histograms have different dimensions \n");
        return 0;
    }    
    Int_t* dims = new Int_t[ndim];
    for (Int_t i = 0; i < ndim; i++) { dims[i] = i; }
    if (useRange) {         
        htemp = hist1->Projection(ndim,dims,"e");
	hist = htemp; 
    } else { hist = hist1; }
    //THnSparse* hnew = hist2->Projection(ndim,dims,"o");
    //hnew->SetName(newname);
    THnSparse* hnew = 0;
    if (overwrite) { 
        hnew = hist2;
    } else {
        hnew = (THnSparse*) hist2->Clone(newname);
    }
    for (Int_t i = 0; i < ndim; i++) { hnew->GetAxis(i)->SetRange(); }
    hnew->SetTitle(hist1->GetTitle());
    hnew->Reset();
    hnew->Sumw2();
    Double_t content;
    Double_t error;
    Int_t* c = new Int_t[ndim];
    Double_t* x = new Double_t[ndim];
    Long64_t n = hist->GetNbins();
    for (Long64_t j = 0; j < n; j++) {
        content = hist->GetBinContent(j,c);
        error = hist->GetBinError(j);
        for (Int_t i = 0; i < ndim; i++) {
            x[i] = hist->GetAxis(i)->GetBinCenter(c[i]);
        }
        /* function IsInRange is protected, shit!
        if (useRange) {
            if (! hist1->IsInRange(c)) continue;
        }
        */
        if (calcErrors) {
           // implementation to be done
        } else {
           hnew->Fill(x,content);
        }
    }
    delete[] c; c=0;
    delete[] x; x=0;
    delete[] dims; dims=0;
    if (htemp) { delete htemp; htemp = 0;}
    return hnew;
}

//_____________________________________________________________________________
AliPWG0Helper::MCProcessType AlidNdPtHelper::GetEventProcessTypePA(AliHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //


  // Check for simple headers first

  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader->GenEventHeader());
  if (dpmJetGenHeader) {
    return GetDPMjetEventProcessTypePA(dpmJetGenHeader,adebug);
  }
  
  // only dpmjet currently supported
  /*
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());
  if (pythiaGenHeader) {
    return AliPWG0Helper::GetPythiaEventProcessType(pythiaGenHeader,adebug);
  }  
  */

  // check for cocktail

  AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
  if (!genCocktailHeader) {
    printf("AlidNdPtHelper::GetEventProcessTypePA : Unknown header type. \n");
    return AliPWG0Helper::kInvalidProcess;
  }

  TList* headerList = genCocktailHeader->GetHeaders();
  if (!headerList) {
    return AliPWG0Helper::kInvalidProcess;
  }

  for (Int_t i=0; i<headerList->GetEntries(); i++) {
    /*
    pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
    if (pythiaGenHeader) {
      return AliPWG0Helper::GetPythiaEventProcessType(pythiaGenHeader,adebug);
    }
    */

    dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(headerList->At(i));
    if (dpmJetGenHeader) {
      return GetDPMjetEventProcessTypePA(dpmJetGenHeader,adebug);
    }
  }
  return AliPWG0Helper::kInvalidProcess;
}


//_____________________________________________________________________________
AliPWG0Helper::MCProcessType AlidNdPtHelper::GetDPMjetEventProcessTypePA(AliGenEventHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  // here kSD means (pure) single diffractive
  // and kND means non-single-diffractive

  // can only read pythia headers, either directly or from cocktalil header
  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader);

  if (!dpmJetGenHeader) {
    printf("AlidNdPtHelper::GetDPMjetProcessTypePA : Unknown header type (not DPMjet or). \n");
    return AliPWG0Helper::kInvalidProcess;
  }
   
  Int_t nsd1=0, nsd2=0, ndd=0;
  dpmJetGenHeader->GetNDiffractive(nsd1,nsd2,ndd);
  if(adebug) {
        printf("%d+%d->%d %d\n",dpmJetGenHeader->ProjectileParticipants(),dpmJetGenHeader->TargetParticipants(),nsd1,nsd2);
  }
  if((dpmJetGenHeader->ProjectileParticipants()==nsd1) && (ndd==0)) { return AliPWG0Helper::kSD; }
  else if ((dpmJetGenHeader->ProjectileParticipants()==nsd2) && (ndd==0)) { return AliPWG0Helper::kSD; }
  else { return AliPWG0Helper::kND; }
}

