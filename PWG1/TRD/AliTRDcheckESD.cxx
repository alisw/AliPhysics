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

/////////////////////////////////////////////////////
//
// Check basic detector results at ESD level
//   - Geometrical efficiency  
//   - Tracking efficiency  
//   - PID efficiency  
//   - Refit efficiency  
//
// Author
//   Alex Bercuci <A.Bercuci@gsi.de>
//
//////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TH2I.h>
#include <TH3S.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDkink.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTrackReference.h"

#include "AliTRDcheckESD.h"

ClassImp(AliTRDcheckESD)

const Float_t AliTRDcheckESD::fgkxTPC = 290.;
const Float_t AliTRDcheckESD::fgkxTOF = 365.;
FILE* AliTRDcheckESD::fgFile = NULL;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTaskSE()
  ,fStatus(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fHistos(NULL)
  ,fResults(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("checkESD", "Check TRD @ ESD level");
  SetMC(kTRUE);
}

AliTRDcheckESD::AliTRDcheckESD(char* name):
  AliAnalysisTaskSE(name)
  ,fStatus(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fHistos(NULL)
  ,fResults(NULL)
{
  //
  // Default constructor
  //
  SetMC(kTRUE);
  DefineOutput(1, TObjArray::Class());
}

//____________________________________________________________________
AliTRDcheckESD::~AliTRDcheckESD()
{
// Destructor
  if(fHistos){
    //fHistos->Delete();
    delete fHistos;
  }
  if(fResults){
    fResults->Delete();
    delete fResults;
  }
}

//____________________________________________________________________
void AliTRDcheckESD::UserCreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
  //OpenFile(0, "RECREATE");  

  Histos();
}

//____________________________________________________________________
TGraph* AliTRDcheckESD::GetGraph(Int_t id, Option_t *opt)
{
// Retrieve graph with "id"
// Possible options are :
//   "b" - build graph if none found
//   "c" - clear existing graph

  Bool_t kBUILD = strstr(opt, "b"), // build graph if none found
         kCLEAR = strstr(opt, "c"); // clear existing graph

  const Char_t *name[] = {
    "Geo", "Trk", "Pid", "Ref", "Max06", "Mean09"
  };
  const Char_t *title[] = {
    "TRD geometrical efficiency (TRDin/TPCout)"
    ,"TRD tracking efficiency (TRDout/TRDin)"
    ,"TRD PID efficiency (TRDpid/TRDin)"
    ,"TRD refit efficiency (TRDrefit/TRDin)"
    ,"TRD Eloss (Max/90% quantile)"
    ,"TRD Eloss (Mean/60% quantile)"
  };
  const Int_t ngr = sizeof(name)/sizeof(Char_t*);
  if(ngr != kNgraphs){
    AliWarning("No of graphs defined different from definition");
    return NULL;
  }

  if(!fResults){
    fResults = new TObjArray(kNgraphs);
    fResults->SetOwner();
    fResults->SetName("results");
  }

  TGraph *g = NULL;
  if((g = dynamic_cast<TGraph*>(fResults->At(id)))){
    if(kCLEAR){ 
      for(Int_t ip=g->GetN(); ip--;) g->RemovePoint(ip);
    } else {
      PutTrendValue(name[id], g->GetMean(2));
      PutTrendValue(Form("%sRMS", name[id]), g->GetRMS(2));
    }
  } else {
    if(kBUILD){
      switch(id){
      case 0:
        g = new TGraphErrors();
        break;
      case 1:
        g = new TGraphErrors();
        break;
      case 2:
        g = new TGraphErrors();
        break;
      case 3:
        g = new TGraphErrors();
        break;
      case 4:
        g = new TGraphAsymmErrors(6);
        g->SetMarkerStyle(22);g->SetMarkerColor(kRed);
        g->SetLineColor(kBlack);g->SetLineWidth(2);
        break;
      case 5:
        g = new TGraphAsymmErrors(6);
        g->SetMarkerStyle(21);
        g->SetLineColor(kRed);g->SetLineWidth(2);
        break;
      default:
        AliWarning(Form("Graph index[%d] missing/not defined.", id));
        return NULL;
      }
      g->SetNameTitle(name[id], title[id]);
      fResults->AddAt(g, id);
    }
  }
  return g;
}

//____________________________________________________________________
void AliTRDcheckESD::UserExec(Option_t *){
  //
  // Run the Analysis
  //
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fMC = MCEvent();
  
  if(!fESD){
    AliError("ESD event missing.");
    return;
  }
  
  // Get MC information if available
  AliStack * fStack = NULL;
  if(HasMC()){
    if(!fMC){ 
      AliWarning("MC event missing");
      SetMC(kFALSE);
    } else {
      if(!(fStack = fMC->Stack())){
        AliWarning("MC stack missing");
        SetMC(kFALSE);
      }
    }
  }
  TH2 *h(NULL);
  
  AliESDtrack *esdTrack(NULL);
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    esdTrack = fESD->GetTrack(itrk);

    // track status
    ULong_t status = esdTrack->GetStatus(); //PrintStatus(status);
    if(!Bool_t(status & AliESDtrack::kTPCout)) continue;
    if(esdTrack->GetKinkIndex(0) > 0) continue;

    //Int_t nTPC(esdTrack->GetNcls(1));
    Int_t nTRD(esdTrack->GetNcls(2));
    Double_t pt(esdTrack->Pt());
    //Double_t eta(esdTrack->Eta());
    //Double_t phi(esdTrack->Phi());
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    // pid quality
    //esdTrack->GetTRDntrackletsPID();
    Bool_t kBarrel = Bool_t(status & AliESDtrack::kTRDin);

    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    const AliExternalTrackParam *ip = esdTrack->GetInnerParam();

    Double_t pt0(0.), eta0(0.), phi0(0.), ptTRD(0.); 
    // read MC info if available
    Bool_t kFOUND(kFALSE), kPhysPrim(kFALSE);
    AliMCParticle *mcParticle(NULL);
    if(HasMC()){
      AliTrackReference *ref(NULL); 
      Int_t fLabel(esdTrack->GetLabel());
      Int_t fIdx(TMath::Abs(fLabel));
      if(fIdx > fStack->GetNtrack()) continue; 
      
      // read MC particle 
      if(!(mcParticle = (AliMCParticle*) fMC->GetTrack(fIdx))) {
        AliWarning(Form("MC particle missing. Label[ %d].", fLabel));
        continue;
      }
      pt0  = mcParticle->Pt();
      eta0 = mcParticle->Eta();
      phi0 = mcParticle->Phi();
      kPhysPrim = fMC->IsPhysicalPrimary(fIdx);

      // read track references
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      if(!nRefs){
        AliWarning(Form("No TR found for track @ Label[%d].", fLabel));
        continue;
      }
      Int_t iref = 0;
      while(iref<nRefs){
        ref = mcParticle->GetTrackReference(iref);
        if(ref->LocalX() > fgkxTPC) break;
        ref=NULL; iref++;
      }
      if(ref){ 
        if(ref->LocalX() > fgkxTOF){ // track skipping TRD fiducial volume
          ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
        }
      } else { // track stopped in TPC 
        ref = mcParticle->GetTrackReference(TMath::Max(iref-1, 0));
      }
      ptTRD = ref->Pt();kFOUND=kTRUE;
    } else { // use reconstructed values
      if(op){
        Double_t x(op->GetX());
        if(x<fgkxTOF && x>fgkxTPC){
          ptTRD=op->Pt();
          kFOUND=kTRUE;
        }
      }

      if(!kFOUND && ip){
        ptTRD=ip->Pt();
        kFOUND=kTRUE;
      }
    }

    if(kFOUND){
      h = (TH2I*)fHistos->At(kTRDstat);
      if(status & AliESDtrack::kTPCout) h->Fill(ptTRD, kTPCout);
      if(status & AliESDtrack::kTRDin) h->Fill(ptTRD, kTRDin);
      if(kBarrel && (status & AliESDtrack::kTRDout)){ 
        ((TH1*)fHistos->At(kNCl))->Fill(nTRD);
        h->Fill(ptTRD, kTRDout);
      }
      if(kBarrel && (status & AliESDtrack::kTRDpid)) h->Fill(ptTRD, kTRDpid);
      if(kBarrel && (status & AliESDtrack::kTRDrefit)) h->Fill(ptTRD, kTRDref);
    }
    if(HasMC() && kBarrel && (status & AliESDtrack::kTRDout)) {
      TH3 *h3 = (TH3S*)fHistos->At(kPtRes);
      Int_t sgn = mcParticle->Charge()>0?1:-1;
      h3->Fill(pt0, 1.e2*pt/pt0-1.e2, sgn*Pdg2Idx(TMath::Abs(mcParticle->PdgCode())));
    }
    if(ip){
      h = (TH2I*)fHistos->At(kTRDmom);
      Float_t pTRD(0.);
      for(Int_t ily=6; ily--;){
        if((pTRD=esdTrack->GetTRDmomentum(ily))<0.) continue;
        h->Fill(ip->GetP()-pTRD, ily);
      }
    }
  }  
  PostData(1, fHistos);
}

//____________________________________________________________________
TObjArray* AliTRDcheckESD::Histos()
{
// Retrieve histograms array if already build or build it

  if(fHistos) return fHistos;

  fHistos = new TObjArray(kNhistos);
  //fHistos->SetOwner(kTRUE);
  
  TH1 *h = NULL;

  // clusters per tracklet
  if(!(h = (TH1I*)gROOT->FindObject("hNCl"))){
    h = new TH1I("hNCl", "Clusters per TRD track;N_{cl}^{TRD};entries", 100, 0., 200.);
  } else h->Reset();
  fHistos->AddAt(h, kNCl);

  // status bits histogram
  const Int_t kNpt(10), kNbits(5);
  Float_t Pt(0.1), Bits(.5);
  Float_t binsPt[kNpt+1], binsBits[kNbits+1];
  for(Int_t i=0;i<kNpt+1; i++,Pt+=(TMath::Exp(i*i*.015)-1.)) binsPt[i]=Pt;
  for(Int_t i=0; i<kNbits+1; i++,Bits+=1.) binsBits[i]=Bits;
  if(!(h = (TH2I*)gROOT->FindObject("hTRDstat"))){
    h = new TH2I("hTRDstat", "TRD status bits;p_{t} @ TRD [GeV/c];status;entris", kNpt, binsPt, kNbits, binsBits);
    TAxis *ay(h->GetYaxis());
    ay->SetBinLabel(1, "kTPCout");
    ay->SetBinLabel(2, "kTRDin");
    ay->SetBinLabel(3, "kTRDout");
    ay->SetBinLabel(4, "kTRDpid");
    ay->SetBinLabel(5, "kTRDrefit");
  } else h->Reset();
  fHistos->AddAt(h, kTRDstat);

  // energy loss
  if(!(h = (TH2I*)gROOT->FindObject("hTRDmom"))){
    h = new TH2I("hTRDmom", "TRD energy loss;p_{inner} - p_{ly} [GeV/c];ly;entries", 100, -1., 2., 6, -0.5, 5.5);
  } else h->Reset();
  fHistos->AddAt(h, kTRDmom);

  // pt resolution
  const Int_t kNdpt(100), kNspec(2*AliPID::kSPECIES+1);
  Float_t DPt(-3.), Spec(-AliPID::kSPECIES-0.5);
  Float_t binsDPt[kNdpt+1], binsSpec[kNspec+1];
  for(Int_t i=0; i<kNdpt+1; i++,DPt+=6.e-2) binsDPt[i]=DPt;
  for(Int_t i=0; i<kNspec+1; i++,Spec+=1.) binsSpec[i]=Spec;
  if(!(h = (TH3S*)gROOT->FindObject("hPtRes"))){
    h = new TH3S("hPtRes", "P_{t} resolution @ DCA;p_{t}^{MC} [GeV/c];#Delta p_{t}/p_{t}^{MC} [%];SPECIES", kNpt, binsPt, kNdpt, binsDPt, kNspec, binsSpec);
    TAxis *az(h->GetZaxis());
    for(Int_t i(0); i<AliPID::kSPECIES; i++){
      az->SetBinLabel(5-i, AliPID::ParticleLatexName(i));
      az->SetBinLabel(7+i, AliPID::ParticleLatexName(i));
    }
  } else h->Reset();
  fHistos->AddAt(h, kPtRes);

  return fHistos;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::Load(const Char_t *filename, const Char_t *name)
{
// Load data from performance file

  if(!TFile::Open(filename)){
    AliWarning(Form("Couldn't open file %s.", filename));
    return kFALSE;
  }
  TObjArray *o = NULL;
  if(!(o = (TObjArray*)gFile->Get(name ? name : GetName()))){
    AliWarning("Missing histogram container.");
    return kFALSE;
  }
  fHistos = (TObjArray*)o->Clone(GetName());
  gFile->Close();
  SETBIT(fStatus, kLoad);
  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDcheckESD::PutTrendValue(const Char_t *name, Double_t val)
{
// Dump trending value to default file

  if(!fgFile){
    fgFile = fopen("TRD.Performance.txt", "at");
  }
  fprintf(fgFile, "%s_%s %f\n", GetName(), name, val);
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
// Steer post-processing 
  if(!IsLoad()){
    fHistos = dynamic_cast<TObjArray *>(GetOutputData(1));
    if(!fHistos){
      AliError("Histogram container not found in output");
      return;
    }
  }

  // geometrical efficiency
  TH2I *h2 = (TH2I*)fHistos->At(kTRDstat);
  TH1 *h1[2] = {NULL, NULL};
  h1[0] = h2->ProjectionX("checkESDx0", kTPCout, kTPCout);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDin, kTRDin);
  Process(h1, (TGraphErrors*)GetGraph(0));
  delete h1[0];delete h1[1];

  // tracking efficiency
  h1[0] = h2->ProjectionX("checkESDx0", kTRDin, kTRDin);
  h1[1] = h2->ProjectionX("checkESDx1", kTRDout, kTRDout);
  Process(h1, (TGraphErrors*)GetGraph(1));
  delete h1[1];

  // PID efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDpid, kTRDpid);
  Process(h1, (TGraphErrors*)GetGraph(2));
  delete h1[1];

  // Refit efficiency
  h1[1] = h2->ProjectionX("checkESDx1", kTRDref, kTRDref);
  Process(h1, (TGraphErrors*)GetGraph(3));
  delete h1[1];
  if(!(h2 = dynamic_cast<TH2I*>(fHistos->At(kTRDmom)))) return;
 
  TGraphAsymmErrors *g06 = (TGraphAsymmErrors*)GetGraph(4), *g09 = (TGraphAsymmErrors*)GetGraph(5);
  TAxis *ax=h2->GetXaxis();
  const Int_t nq(4);
  const Double_t xq[nq] = {0.05, 0.2, 0.8, 0.95};
  Double_t yq[nq];
  for(Int_t ily=6; ily--;){
    h1[0] = h2->ProjectionX("checkESDp0", ily+1, ily+1);
    h1[0]->GetQuantiles(nq,yq,xq);
    g06->SetPoint(ily, Float_t(ily), ax->GetBinCenter(h1[0]->GetMaximumBin()));
    g06->SetPointError(ily, 0., 0., TMath::Abs(yq[0]), yq[3]);
    g09->SetPoint(ily, Float_t(ily), h1[0]->GetMean());
    g09->SetPointError(ily, 0., 0., TMath::Abs(yq[1]), yq[2]);

    //printf(" max[%f] mean[%f] q[%f %f %f %f]\n", ax->GetBinCenter(h1[0]->GetMaximumBin()), h1[0]->GetMean(), yq[0], yq[1], yq[2], yq[3]);
    delete h1[0];
  }
}

//____________________________________________________________________
Int_t AliTRDcheckESD::Pdg2Idx(Int_t pdg)
{
  switch(pdg){
  case kElectron: return AliPID::kElectron+1;  
  case kMuonMinus: return AliPID::kMuon+1;  
  case kPiPlus: return AliPID::kPion+1;  
  case kKPlus: return AliPID::kKaon+1;
  case kProton: return AliPID::kProton+1;
  } 
  return 0;
}

//____________________________________________________________________
void AliTRDcheckESD::Process(TH1 **h1, TGraphErrors *g)
{
// Generic function to process one reference plot

  Int_t n1 = 0, n2 = 0, ip=0;
  Double_t eff = 0.;

  TAxis *ax = h1[0]->GetXaxis();
  for(Int_t ib=1; ib<=ax->GetNbins(); ib++){
    if(!(n1 = (Int_t)h1[0]->GetBinContent(ib))) continue;
    n2 = (Int_t)h1[1]->GetBinContent(ib);
    eff = n2/Float_t(n1);

    ip=g->GetN();
    g->SetPoint(ip, ax->GetBinCenter(ib), eff);
    g->SetPointError(ip, 0., n2 ? eff*TMath::Sqrt(1./n1+1./n2) : 0.);
  }
}  

//____________________________________________________________________
void AliTRDcheckESD::PrintStatus(ULong_t status)
{
// Dump track status to stdout

  printf("ITS[i(%d) o(%d) r(%d)] TPC[i(%d) o(%d) r(%d) p(%d)] TRD[i(%d) o(%d) r(%d) p(%d) s(%d)] HMPID[o(%d) p(%d)]\n"
    ,Bool_t(status & AliESDtrack::kITSin)
    ,Bool_t(status & AliESDtrack::kITSout)
    ,Bool_t(status & AliESDtrack::kITSrefit)
    ,Bool_t(status & AliESDtrack::kTPCin)
    ,Bool_t(status & AliESDtrack::kTPCout)
    ,Bool_t(status & AliESDtrack::kTPCrefit)
    ,Bool_t(status & AliESDtrack::kTPCpid)
    ,Bool_t(status & AliESDtrack::kTRDin)
    ,Bool_t(status & AliESDtrack::kTRDout)
    ,Bool_t(status & AliESDtrack::kTRDrefit)
    ,Bool_t(status & AliESDtrack::kTRDpid)
    ,Bool_t(status & AliESDtrack::kTRDStop)
    ,Bool_t(status & AliESDtrack::kHMPIDout)
    ,Bool_t(status & AliESDtrack::kHMPIDpid)
  );
}

