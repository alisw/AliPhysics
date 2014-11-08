// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011

#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TTree.h>
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliDhcTask.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDMuonTrack.h"
#include "AliESDtrackCuts.h"
#include "AliPool.h"
#include "AliVParticle.h"

using std::cout;
using std::endl;

ClassImp(AliDhcTask)

//________________________________________________________________________
AliDhcTask::AliDhcTask()
: AliAnalysisTaskSE(), fVerbosity(0), fEtaMax(1), fZVtxMax(10), fPtMin(0.25), fPtMax(15),
  fTrackDepth(1000), fPoolSize(200), fTracksName(), fDoWeights(kFALSE), fFillMuons(kFALSE),
  fRequireMuon(kFALSE), fReqPtLo(0.0), fReqPtHi(1000.0),
  fPtTACrit(kTRUE), fAllTAHists(kFALSE), fMixInEtaT(kFALSE),
  fEtaTLo(-1.0), fEtaTHi(1.0), fEtaALo(-1.0), fEtaAHi(1.0), fOmitFirstEv(kTRUE),
  fDoFillSame(kFALSE), fDoMassCut(kFALSE), fCheckVertex(kTRUE), fClassName(),
  fCentMethod("V0M"), fNBdeta(20), fNBdphi(36), fTriggerMatch(kTRUE),
  fBPtT(0x0), fBPtA(0x0), fBCent(0x0), fBZvtx(0x0),
  fMixBCent(0x0), fMixBZvtx(0x0), fHEffT(0x0), fHEffA(0x0),
  fESD(0x0), fAOD(0x0), fOutputList(0x0), fHEvt(0x0), fHEvtWTr(0x0), fHTrk(0x0), fHPoolReady(0x0),
  fHPtAss(0x0), fHPtTrg(0x0), fHPtTrgEvt(0x0),
  fHPtTrgNorm1S(0x0), fHPtTrgNorm1M(0x0), fHPtTrgNorm2S(0x0), fHPtTrgNorm2M(0x0),
  fHCent(0x0), fHZvtx(0x0), fNbins(0), fHSs(0x0), fHMs(0x0), fHPts(0x0), fHSMass(0x0), fHMMass(0x0),
  fHQATp(0x0), fHQAAp(0x0), fHQATpCorr(0x0), fHQAApCorr(0x0),
  fHQATm(0x0), fHQAAm(0x0), fHQATmCorr(0x0), fHQAAmCorr(0x0),
  fHPtCentT(0x0), fHPtCentA(0x0), fIndex(0x0),
  fCentrality(99), fZVertex(99), fEsdTPCOnly(0), fPoolMgr(0),
  fUtils(0x0)
{
  // Constructor
}

//________________________________________________________________________
AliDhcTask::AliDhcTask(const char *name, Bool_t def) 
: AliAnalysisTaskSE(name), fVerbosity(0), fEtaMax(1), fZVtxMax(10), fPtMin(0.25), fPtMax(15),
  fTrackDepth(1000), fPoolSize(200), fTracksName(), fDoWeights(kFALSE), fFillMuons(kFALSE),
  fRequireMuon(kFALSE), fReqPtLo(0.0), fReqPtHi(1000.0),
  fPtTACrit(kTRUE), fAllTAHists(kFALSE), fMixInEtaT(kFALSE),
  fEtaTLo(-1.0), fEtaTHi(1.0), fEtaALo(-1.0), fEtaAHi(1.0), fOmitFirstEv(kTRUE),
  fDoFillSame(kFALSE), fDoMassCut(kFALSE), fCheckVertex(kTRUE), fClassName(),
  fCentMethod("V0M"), fNBdeta(20), fNBdphi(36), fTriggerMatch(kTRUE),
  fBPtT(0x0), fBPtA(0x0), fBCent(0x0), fBZvtx(0x0),
  fMixBCent(0x0), fMixBZvtx(0x0), fHEffT(0x0), fHEffA(0x0),
  fESD(0x0), fAOD(0x0), fOutputList(0x0), fHEvt(0x0), fHEvtWTr(0x0), fHTrk(0x0), fHPoolReady(0x0),
  fHPtAss(0x0), fHPtTrg(0x0), fHPtTrgEvt(0x0),
  fHPtTrgNorm1S(0x0), fHPtTrgNorm1M(0x0), fHPtTrgNorm2S(0x0), fHPtTrgNorm2M(0x0),
  fHCent(0x0), fHZvtx(0x0), fNbins(0), fHSs(0x0), fHMs(0x0), fHPts(0x0), fHSMass(0x0), fHMMass(0x0),
  fHQATp(0x0), fHQAAp(0x0), fHQATpCorr(0x0), fHQAApCorr(0x0),
  fHQATm(0x0), fHQAAm(0x0), fHQATmCorr(0x0), fHQAAmCorr(0x0),
  fHPtCentT(0x0), fHPtCentA(0x0), fIndex(0x0),
  fCentrality(99), fZVertex(99), fEsdTPCOnly(0), fPoolMgr(0),
  fUtils(0x0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,Tracks "
               "AOD:header,tracks,vertices,";

  if (def) {
    Double_t ptt[4] = {0.25, 1.0, 2.0, 15.0};
    fBPtT  = new TAxis(3,ptt); 
    Double_t pta[4] = {0.25, 1.0, 2.0, 15.0};
    fBPtA  = new TAxis(3,pta); 
    Double_t cent[2] = {-100.0, 100.0};
    fBCent = new TAxis(1,cent);
    Double_t zvtx[2] = {-10, 10};
    fBZvtx = new TAxis(1,zvtx);
    Double_t centmix[2] = {-100.0, 100.0};
    fMixBCent = new TAxis(1,centmix);
    Double_t zvtxmix[9] = {-10,-6,-4,-2,0,2,4,6,10};
    fMixBZvtx = new TAxis(8,zvtxmix);
  }
}

//________________________________________________________________________
void AliDhcTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (per slave on PROOF!)
  PrintDhcSettings();

  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fUtils = new AliAnalysisUtils();

  BookHistos();
  InitEventMixer(); 
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliDhcTask::PrintDhcSettings()
{
  AliInfo(Form("Dhc Task %s settings",fName.Data()));
  AliInfo(Form(" centrality estimator %s", fCentMethod.Data()));
  AliInfo(Form(" using tracks named %s", fTracksName.Data()));
  AliInfo(Form(" efficiency correct triggers? %d", fHEffT!=0));
  if (fHEffT!=0) {
    AliInfo(Form(" %d dimensions (t)", fHEffT->GetNdimensions()));
  }
  AliInfo(Form(" efficiency correct associates? %d", fHEffA!=0));
  if (fHEffA!=0) {
    AliInfo(Form(" %d dimensions (a)", fHEffA->GetNdimensions()));
  }
  AliInfo(Form(" fill muons? %d", fFillMuons));
  AliInfo(Form(" require muon even if not filling them? %d", fRequireMuon));
  if (fRequireMuon) AliInfo(Form(" require muon with %f < pt < %f", fReqPtLo, fReqPtHi));
  AliInfo(Form(" use pTT > pTA criterion? %d", fPtTACrit));
  AliInfo(Form(" create all pTT, pTA hists? %d", fAllTAHists));
  AliInfo(Form(" Mix in eta_T bins instead of z_vertex? %d", fMixInEtaT));
  AliInfo(Form(" trigger eta range %f .. %f", fEtaTLo, fEtaTHi));
  AliInfo(Form(" associate eta range %f .. %f", fEtaALo, fEtaAHi));
  AliInfo(Form(" fill same event in any case? %d", fDoFillSame));
  AliInfo(Form(" do invariant mass cut? %d", fDoMassCut));
  AliInfo(Form(" omit first event? %d\n", fOmitFirstEv));
  AliInfo(Form(" check the vertex? %d\n", fCheckVertex));
}

//________________________________________________________________________
void AliDhcTask::BookHistos()
{
  // Book histograms.

  if (fVerbosity > 1) {
    AliInfo(Form("Number of pt(t) bins: %d", fBPtT->GetNbins()));
    for (Int_t i=1; i<=fBPtT->GetNbins(); i++) {
      AliInfo(Form("pt(t) bin %d, %f to %f", i, fBPtT->GetBinLowEdge(i), fBPtT->GetBinUpEdge(i)));
    }
    AliInfo(Form("Number of pt(a) bins: %d", fBPtA->GetNbins()));
    for (Int_t i=1; i<=fBPtA->GetNbins(); i++) {
      AliInfo(Form("pt(a) bin %d, %f to %f", i, fBPtA->GetBinLowEdge(i), fBPtA->GetBinUpEdge(i)));
    }
  }

  Int_t nPtAssc=fBPtA->GetNbins();
  Int_t nPtTrig=fBPtT->GetNbins();
  Int_t nCent=fBCent->GetNbins();
  Int_t nZvtx=fBZvtx->GetNbins();
  Double_t ptt[nPtTrig+1];
  ptt[0] = fBPtT->GetBinLowEdge(1);
  for (Int_t i=1; i<=nPtTrig; i++) {
    ptt[i] = fBPtT->GetBinUpEdge(i);
  }
  Double_t pta[nPtAssc+1];
  pta[0] = fBPtA->GetBinLowEdge(1);
  for (Int_t i=1; i<=nPtAssc; i++) {
    pta[i] = fBPtA->GetBinUpEdge(i);
  }
  Double_t cent[nCent+1];
  cent[0] = fBCent->GetBinLowEdge(1);
  for (Int_t i=1; i<=nCent; i++) {
    cent[i] = fBCent->GetBinUpEdge(i);
  }
  Double_t zvtx[nZvtx+1];
  zvtx[0] = fBZvtx->GetBinLowEdge(1);
  for (Int_t i=1; i<=nZvtx; i++) {
    zvtx[i] = fBZvtx->GetBinUpEdge(i);
  }
  
  fHEvt = new TH2F("fHEvt", "Event-level variables; Zvtx; Cent", nZvtx, zvtx, nCent, cent);
  fOutputList->Add(fHEvt);
  fHEvtWTr = new TH2F("fHEvtWTr", "Events with tracks; Zvtx; Cent", nZvtx, zvtx, nCent, cent);
  fOutputList->Add(fHEvtWTr);
  fHTrk = new TH2F("fHTrk", "Track-level variables",
                   100, 0, TMath::TwoPi(), 100, -fEtaMax, +fEtaMax);
  fOutputList->Add(fHTrk);
  
  // Centrality mixing axis
  Int_t nCentBins=fMixBCent->GetNbins();
  Double_t centBins[nCentBins+1];
  centBins[0] = fMixBCent->GetBinLowEdge(1);
  for (Int_t i=1; i<=nCentBins; i++) {
    centBins[i] = fMixBCent->GetBinUpEdge(i);
  }
  // Z-vertex mixing axis
  Int_t nZvtxBins=fMixBZvtx->GetNbins();
  Double_t zvtxbin[nZvtxBins+1];
  zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
  for (Int_t i=1; i<=nZvtxBins; i++) {
    zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
  }
  fHPoolReady = new TH2F("fHPoolReady","mixing started", nZvtxBins, zvtxbin, nCentBins, centBins);
  fOutputList->Add(fHPoolReady);
  
  fHPtAss = new TH1F("fHPtAss","PtAssoc;P_{T} (GeV/c) [GeV/c]",nPtAssc,pta);
  fOutputList->Add(fHPtAss);
  fHPtTrg = new TH1F("fHPtTrg","PtTrig;P_{T} (GeV/c) [GeV/c]",nPtTrig,ptt);
  fOutputList->Add(fHPtTrg);
  fHPtTrgEvt = new TH1F("fHPtTrgEvt","PtTrig;P_{T} (GeV/c) [GeV/c]",nPtTrig,ptt);
  fHPtTrgNorm1S = new TH3F("fHPtTrgNorm1S","PtTrig;P_{T} (GeV/c) [GeV/c];centrality;z_{vtx}",nPtTrig,ptt,nCent,cent,nZvtx,zvtx);
  fOutputList->Add(fHPtTrgNorm1S);
  fHPtTrgNorm1M = (TH3*) fHPtTrgNorm1S->Clone("fHPtTrgNorm1M");
  fOutputList->Add(fHPtTrgNorm1M);
  fHPtTrgNorm2S = (TH3*) fHPtTrgNorm1S->Clone("fHPtTrgNorm2S");
  fOutputList->Add(fHPtTrgNorm2S);
  fHPtTrgNorm2M = (TH3*) fHPtTrgNorm1S->Clone("fHPtTrgNorm2M");
  fOutputList->Add(fHPtTrgNorm2M);
  
  fHCent = new TH1F("fHCent","Cent;bins",nCent,cent);
  fOutputList->Add(fHCent);
  fHZvtx = new TH1F("fHZvtx","Zvertex;bins",nZvtx,zvtx);
  fOutputList->Add(fHZvtx);
  
  fHQATp = new TH3F("fHQATp","QA trigger;p_{T} (GeV/c);#eta;#phi",
                    100,0.0,10.0,
                    40,fEtaTLo,fEtaTHi,
                    36,0.0,TMath::TwoPi());
  fOutputList->Add(fHQATp);
  fHQAAp = new TH3F("fHQAAp","QA associated;p_{T} (GeV/c);#eta;#phi",
                    100,0.0,10.0,
                    40,fEtaALo,fEtaAHi,
                    36,0.0,TMath::TwoPi());
  fOutputList->Add(fHQAAp);
  fHQATpCorr = (TH3 *) fHQATp->Clone("fHQATpCorr");
  fOutputList->Add(fHQATpCorr);
  fHQAApCorr = (TH3 *) fHQAAp->Clone("fHQAApCorr");
  fOutputList->Add(fHQAApCorr);
  fHQATm = (TH3 *) fHQATp->Clone("fHQATm");
  fOutputList->Add(fHQATm);
  fHQAAm = (TH3 *) fHQAAp->Clone("fHQAAm");
  fOutputList->Add(fHQAAm);
  fHQATmCorr = (TH3 *) fHQATm->Clone("fHQATmCorr");
  fOutputList->Add(fHQATmCorr);
  fHQAAmCorr = (TH3 *) fHQAAm->Clone("fHQAAmCorr");
  fOutputList->Add(fHQAAmCorr);

  fHPtCentT = new TH2F("fHPtCentT",Form("trigger particles;p_{T} (GeV/c);centrality (%s)",fCentMethod.Data()),
                       100,0.0,10.0,
                       100,cent[0],cent[nCent]);
  fOutputList->Add(fHPtCentT);
  fHPtCentA = new TH2F("fHPtCentA",Form("associated particles;p_{T} (GeV/c);centrality (%s)",fCentMethod.Data()),
                       100,0.0,10.0,
                       100,cent[0],cent[nCent]);
  fOutputList->Add(fHPtCentA);

  fNbins  = nPtTrig*nPtAssc*nCent*nZvtx;
  fHSs    = new TH2*[fNbins];
  fHMs    = new TH2*[fNbins];
  fHPts   = new TH1*[fNbins];
  fHSMass = new TH1*[fNbins];
  fHMMass = new TH1*[fNbins];

  fIndex = new TFormula("GlobIndex","(t-1)*[0]*[1]*[2]+(z-1)*[0]*[1]+(x-1)*[1]+(y-1)");
  fIndex->SetParameters(nPtTrig,nPtAssc,nZvtx);
  fIndex->SetParNames("NTrigBins","NAssocBins","NZvertexBins");
  fOutputList->Add(fIndex);
  
  Double_t dEtaMin = fEtaTLo - fEtaAHi;
  Double_t dEtaMax = fEtaTHi - fEtaALo;
  
  Int_t count = 0;
  for (Int_t c=1; c<=nCent; ++c) {
    for (Int_t z=1; z<=nZvtx; ++z) {
      for (Int_t t=1; t<=nPtTrig; ++t) {
        for (Int_t a=1; a<=nPtAssc; ++a) {
          fHSs[count]    = 0;
          fHMs[count]    = 0;
          fHPts[count]   = 0;
          fHSMass[count] = 0;
          fHMMass[count] = 0;
          if ((a>t)&&!fAllTAHists) {
            ++count;
            continue;
          }
          if (z==1 && t==1 && a==1) {
            TString title(Form("cen=%d (%.1f to %.1f)",
                               c, fBCent->GetBinLowEdge(c), fBCent->GetBinUpEdge(c)));
            fHSMass[count] = new TH1F(Form("hSMass%d",count), Form("Mass Same Event %s;m (GeV)",title.Data()), 10000, 0,10);
            fOutputList->Add(fHSMass[count]);
            fHMMass[count] = new TH1F(Form("hMMass%d",count), Form("Mass Mixed Event %s;m (GeV)",title.Data()), 10000, 0,10);
            fOutputList->Add(fHMMass[count]);
          }
          if (t==1 && a==1) {
            TString title(Form("cen=%d (%.1f to %.1f), zVtx=%d (%.1f to %.1f)",
                               c, fBCent->GetBinLowEdge(c), fBCent->GetBinUpEdge(c),
                               z, fBZvtx->GetBinLowEdge(z), fBZvtx->GetBinUpEdge(z)));
            fHPts[count] = new TH1F(Form("hPt%d",count), Form("Ptdist %s;p_{T} (GeV/c)",title.Data()), 200,0,20);
            fOutputList->Add(fHPts[count]);
          }
          TString title(Form("cen=%d (%.1f to %.1f), zVtx=%d (%.1f to %.1f), trig=%d (%.1f to %.1f), assoc=%d (%.1f to %.1f)",
                             c, fBCent->GetBinLowEdge(c), fBCent->GetBinUpEdge(c),
                             z, fBZvtx->GetBinLowEdge(z), fBZvtx->GetBinUpEdge(z),
                             t, fBPtT->GetBinLowEdge(t),  fBPtT->GetBinUpEdge(t),
                             a, fBPtA->GetBinLowEdge(a),  fBPtA->GetBinUpEdge(a)));
          fHSs[count] = new TH2F(Form("hS%d",count), Form("Signal %s;#Delta #varphi;#Delta #eta",title.Data()),
                                 fNBdphi,-0.5*TMath::Pi(),1.5*TMath::Pi(),fNBdeta,dEtaMin,dEtaMax);
          fHMs[count] = new TH2F(Form("hM%d",count), Form("Mixed %s;#Delta #varphi;#Delta #eta",title.Data()),
                                 fNBdphi,-0.5*TMath::Pi(),1.5*TMath::Pi(),fNBdeta,dEtaMin,dEtaMax);
          fOutputList->Add(fHSs[count]);
          fOutputList->Add(fHMs[count]);
          Int_t tcount = (Int_t)fIndex->Eval(t,a,z,c);
          if (fVerbosity>5)
            cout << count << " " << tcount << ": " << title << endl;
          if (count != tcount)
            AliFatal(Form("Index mismatch: %d %d", count, tcount));
          ++count;
        }
      }
    }
  }
}

//________________________________________________________________________
void AliDhcTask::SetAnaMode(Int_t iAna)
{
  if (iAna==kHH) {
    SetFillMuons(kFALSE);
    fEtaTLo = -1.6;
    fEtaTHi = +1.6;
    fEtaALo = -1.6;
    fEtaAHi = +1.6;
  } else if (iAna==kMuH) {
    SetFillMuons(kTRUE);
    fEtaTLo = -5.0;
    fEtaTHi = -2.0;
    fEtaALo = -1.6;
    fEtaAHi = +1.6;
  } else if (iAna==kHMu) {
    SetFillMuons(kTRUE);
    fEtaTLo = -1.6;
    fEtaTHi = +1.6;
    fEtaALo = -5.0;
    fEtaAHi = -2.0;
  } else if (iAna==kMuMu) {
    SetFillMuons(kTRUE);
    fEtaTLo = -5.0;
    fEtaTHi = -2.0;
    fEtaALo = -5.0;
    fEtaAHi = -2.0;
  } else if (iAna==kPSide) {
    SetFillMuons(kFALSE);
    fEtaTLo = -1.6;
    fEtaTHi = -0.465;
    fEtaALo = -1.6;
    fEtaAHi = +1.6;
  } else if (iAna==kASide) {
    SetFillMuons(kFALSE);
    fEtaTLo = +0.465;
    fEtaTHi = +1.6;
    fEtaALo = -1.6;
    fEtaAHi = +1.6;
  } else {
    AliInfo(Form("Unrecognized analysis option: %d", iAna));
  }
}

//________________________________________________________________________
void AliDhcTask::InitEventMixer()
{
  // The effective pool size in events is set by trackDepth, so more
  // low-mult events are required to maintain the threshold than
  // high-mult events. Centrality pools are indep. of data histogram
  // binning, no need to match.

  // Centrality pools
  Int_t nCentBins=fMixBCent->GetNbins();
  Double_t centBins[nCentBins+1];
  centBins[0] = fMixBCent->GetBinLowEdge(1);
  for (Int_t i=1; i<=nCentBins; i++) {
    centBins[i] = fMixBCent->GetBinUpEdge(i);
  }
 
  // Z-vertex pools
  Int_t nZvtxBins=fMixBZvtx->GetNbins();
  Double_t zvtxbin[nZvtxBins+1];
  zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
  for (Int_t i=1; i<=nZvtxBins; i++) {
    zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
  }

  fPoolMgr = new AliEvtPoolManager();
  fPoolMgr->SetTargetTrackDepth(fTrackDepth);
  if (fVerbosity>4)
    fPoolMgr->SetDebug(1);
  fPoolMgr->InitEventPools(fPoolSize, nCentBins, centBins, nZvtxBins, zvtxbin);
}

//________________________________________________________________________
void AliDhcTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  LoadBranches();

  if (fOmitFirstEv) {
    if (fUtils->IsFirstEventInChunk(InputEvent())) 
      return;
  }

  // Get event pointers, check for signs of life
  Int_t dType = -1;       // Will be set to kESD or kAOD.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (fESD)
    dType = kESD;
  else if (fAOD)
    dType = kAOD;
  else {
    AliError("Neither fESD nor fAOD available");
    return;
  }

  // select specific trigger classes?
  if (fClassName.Length()>0) {
    TString strFiredClass;
    if (fESD)
      strFiredClass = fESD->GetFiredTriggerClasses();
    else
      strFiredClass = fAOD->GetFiredTriggerClasses();
    
    if (fVerbosity > 10) {
      AliInfo(Form("Trigger class selection: This event has classes %s", strFiredClass.Data()));
      AliInfo(Form("Trigger class selection: selecting classes %s", fClassName.Data()));
    }

    TObjArray *arrClass = fClassName.Tokenize(",");
    Int_t nClass = arrClass->GetEntries();
    
    TString strOneClass;
    Bool_t bAccEvent = kFALSE;
    for (Int_t iClass=0; iClass<nClass; iClass++) {
      strOneClass = arrClass->At(iClass)->GetName();
      if (strFiredClass.Contains(strOneClass))
        bAccEvent = kTRUE;
    }
    
    if (!bAccEvent)
      return;
    
    if (fVerbosity > 10) {
      AliInfo(Form("After trigger class selection: This event has classes %s", strFiredClass.Data()));
    }
  }

  Bool_t mcgen = 0;
  if (fTracksName.Contains("Gen"))
    mcgen = 1;

  // Centrality, vertex, other event variables...
  if (mcgen) {
    fZVertex = 0;
    TList *list = InputEvent()->GetList();
    TClonesArray *tcaTracks = 0;
    if (list) 
      tcaTracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));
    if (tcaTracks) 
      fCentrality = tcaTracks->GetEntries();
  } else {
    if (dType == kESD) {
      const AliESDVertex* vertex = fESD->GetPrimaryVertex();
      fZVertex = vertex->GetZ();
      if (fCheckVertex && !VertexOk()) {
        if (fVerbosity > 1)
          AliInfo(Form("Event REJECTED (ESD vertex not OK). z = %.1f", fZVertex));
        return;
      }
      if(fESD->GetCentrality()) {
        fCentrality = 
          fESD->GetCentrality()->GetCentralityPercentile(fCentMethod);
      }
    } else if (dType == kAOD) {
      const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
      fZVertex = vertex->GetZ();
      if (fCheckVertex && !VertexOk()) {
        if (fVerbosity > 1)
          Info("Exec()", "Event REJECTED (AOD vertex not OK). z = %.1f", fZVertex);
        return;
      }
      AliAODHeader * header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
      if(!header) AliFatal("Not a standard AOD");

      const AliCentrality *aodCent = header->GetCentralityP();
      if (aodCent) {
        fCentrality = aodCent->GetCentralityPercentile(fCentMethod);
      }
    }
  }

  // Fill event histogram
  fHEvt->Fill(fZVertex, fCentrality);
  if (fCentrality > fBCent->GetXmax() || fCentrality < fBCent->GetXmin()) {
    if (fVerbosity > 1)
      AliInfo(Form("Event REJECTED (centrality out of range). fCentrality = %.1f", fCentrality));
    return;
  }

  // Get pool containing tracks from other events like this one
  AliEvtPool* pool = fPoolMgr->GetEventPool(fCentrality, fZVertex);
  if (!pool) {
    AliWarning(Form("No pool found. Centrality %f, ZVertex %f", 
		    fCentrality, fZVertex));
    return;
  }

  // Get array of selected tracks
  MiniEvent* sTracks = new MiniEvent(0); // deleted by pool mgr.
  if (dType == kESD) {
    GetESDTracks(sTracks);
  } else if (dType == kAOD) {
    GetAODTracks(sTracks);
  }

  if (sTracks->size()==0) {
    if (fVerbosity > 1)
      AliWarning(Form("Track array empty"));
    delete sTracks;
    return;
  }
  
  fHEvtWTr->Fill(fZVertex, fCentrality);

  if (!pool->IsReady()) {
    pool->UpdatePool(sTracks);
    if (fDoFillSame) { // fill same event right away if requested
      Correlate(*sTracks, *sTracks, kSameEvt);  
    }
    PostData(1, fOutputList);
    return;
  }

  if (pool->IsFirstReady()) {
    fHPoolReady->Fill(fZVertex,fCentrality);
    // recover events missed before the pool is ready
    Int_t nEvs = pool->GetCurrentNEvents();
    if (nEvs>1) {
      for (Int_t i=0; i<nEvs; ++i) {
        MiniEvent* evI = pool->GetEvent(i);
        for (Int_t j=0; j<nEvs; ++j) {
          MiniEvent* evJ = pool->GetEvent(j);
          if (i==j) {
            if (!fDoFillSame) { // do not fill the same events twice
              Correlate(*evI, *evJ, kSameEvt);
            }
          } else {
            Correlate(*evI, *evJ, kDiffEvt);
          }
        }
      }
    }
  } else { /* standard case: same event, then mix*/
    Correlate(*sTracks, *sTracks, kSameEvt);  
    Int_t nMix = pool->GetCurrentNEvents();
    for (Int_t jMix=0; jMix<nMix; ++jMix) {
      MiniEvent* bgTracks = pool->GetEvent(jMix);
      Correlate(*sTracks, *bgTracks, kDiffEvt);
    }
  }

  pool->UpdatePool(sTracks);
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliDhcTask::GetESDTracks(MiniEvent* miniEvt)
{
  // check / count muons
  Int_t nGoodMuons = 0;
  Bool_t bAtLeastOneMuon = kFALSE;
  if (fFillMuons || fRequireMuon) {
    for (Int_t iMu = 0; iMu<fESD->GetNumberOfMuonTracks(); iMu++) {
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iMu);
      if (muonTrack) {
        if (!IsGoodMUONtrack(*muonTrack))
          continue;
        Double_t ptMu   = muonTrack->Pt();
        if (ptMu>fReqPtLo && ptMu<fReqPtHi)
          bAtLeastOneMuon = kTRUE;
        Double_t etaMu  = muonTrack->Eta();
        if (IsTrigger(etaMu,ptMu)||IsAssociated(etaMu,ptMu))
          nGoodMuons++;
      }
    }
  }
  if (fRequireMuon && !bAtLeastOneMuon)
    return;
  
  // Loop twice: 1. Count sel. tracks. 2. Fill vector.
  if (fTracksName.IsNull()) {
    const AliESDVertex *vtxSPD = fESD->GetPrimaryVertexSPD();
    if (!vtxSPD)
      return;
    
    Int_t nTrax = fESD->GetNumberOfTracks();
    if (fVerbosity > 2)
      AliInfo(Form("%d tracks in event",nTrax));

    // Loop 1.
    Int_t nSelTrax = 0;
    TObjArray arr(nTrax);
    arr.SetOwner(1);

    if (!fEsdTPCOnly) {
      fEsdTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      //fEsdTPCOnly->SetMinNClustersTPC(70);
      fEsdTPCOnly->SetMinNCrossedRowsTPC(70);
      fEsdTPCOnly->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    }

    for (Int_t i = 0; i < nTrax; ++i) {
      AliESDtrack* esdtrack = fESD->GetTrack(i);
      if (!esdtrack) {
        AliError(Form("Couldn't get ESD track %d\n", i));
        continue;
      }
      Bool_t trkOK = fEsdTPCOnly->AcceptTrack(esdtrack);
      if (!trkOK)
        continue;
      Double_t pt = esdtrack->Pt();
      Bool_t ptOK = pt >= fPtMin && pt < fPtMax;
      if (!ptOK)
        continue;
      Double_t eta = esdtrack->Eta();
      if (TMath::Abs(eta) > fEtaMax)
        continue;

      // create a tpc only track
      AliESDtrack *newtrack = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
      if(!newtrack)
        continue;
      if (newtrack->Pt()<=0) {
        delete newtrack;
        continue;
      }

      AliExternalTrackParam exParam;
      Bool_t relate = newtrack->RelateToVertexTPC(vtxSPD,fESD->GetMagneticField(),kVeryBig,&exParam);
      if (!relate) {
        delete newtrack;
        continue;
      }

      // set the constraint parameters to the track
      newtrack->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());

      pt = newtrack->Pt();
      ptOK = pt >= fPtMin && pt < fPtMax;
      if (!ptOK) {
        delete newtrack;
        continue;
      }
      eta  = esdtrack->Eta();
      if (TMath::Abs(eta) > fEtaMax) {
        delete newtrack;
        continue;
      }
      arr.Add(newtrack);
      nSelTrax++;
    }

    for(Int_t itrack = 0; itrack < nSelTrax; itrack++) {
      AliVTrack *esdtrack = static_cast<AliESDtrack*>(arr.At(itrack));
      if(!esdtrack) {
        AliError(Form("ERROR: Could not retrieve esdtrack %d",itrack));
        continue;
      }
      Double_t pt = esdtrack->Pt();
      Double_t eta  = esdtrack->Eta();
      Double_t phi  = esdtrack->Phi();
      Int_t    sign = esdtrack->Charge() > 0 ? 1 : -1;
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        miniEvt->push_back(AliMiniTrack(pt, eta, phi, sign));
    }
  } else {
    TList *list = InputEvent()->GetList();
    TClonesArray *tcaTracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));

    if (!tcaTracks){
      AliError("Ptr to tcaTracks zero");
      return;
    }

    const Int_t ntracks = tcaTracks->GetEntries();
    Int_t nGoodTracks = 0;
    // count good tracks
    for (Int_t itrack = 0; itrack < ntracks; itrack++) {
      AliVTrack *vtrack = static_cast<AliVTrack*>(tcaTracks->At(itrack));
      if (!vtrack) {
        AliError(Form("ERROR: Could not retrieve vtrack %d",itrack));
        continue;
      }
      Double_t pt   = vtrack->Pt();
      Double_t eta  = vtrack->Eta();
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        nGoodTracks++;
    }
    if (miniEvt)
      miniEvt->reserve(nGoodTracks);
    else {
      AliError("Ptr to miniEvt zero");
      return;
    }
    // fill good tracks into minievent
    for (Int_t itrack = 0; itrack < ntracks; itrack++) {
      AliVTrack *vtrack = static_cast<AliVTrack*>(tcaTracks->At(itrack));
      if (!vtrack) {
        AliError(Form("ERROR: Could not retrieve vtrack %d",itrack));
        continue;
      }
      Double_t pt   = vtrack->Pt();
      Double_t eta  = vtrack->Eta();
      Double_t phi  = vtrack->Phi();
      Int_t    sign = vtrack->Charge() > 0 ? 1 : -1;
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        miniEvt->push_back(AliMiniTrack(pt, eta, phi, sign));
    }
  }
  
  // fill muons into minievent
  if (fFillMuons) {
    miniEvt->reserve(miniEvt->size()+nGoodMuons);
    // fill them into the mini event
    for (Int_t iMu = 0; iMu<fESD->GetNumberOfMuonTracks(); iMu++) {
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iMu);
      if (muonTrack) {
        if (!IsGoodMUONtrack(*muonTrack)) 
          continue;
        Double_t ptMu   = muonTrack->Pt();
        Double_t etaMu  = muonTrack->Eta();
        Double_t phiMu  = muonTrack->Phi();
        Int_t    signMu = muonTrack->Charge() > 0 ? 1 : -1;
        if (IsTrigger(etaMu,ptMu)||IsAssociated(etaMu,ptMu))
          miniEvt->push_back(AliMiniTrack(ptMu, etaMu, phiMu, signMu));
      }
    }
  }

}

//________________________________________________________________________
void AliDhcTask::GetAODTracks(MiniEvent* miniEvt)
{
  // check / count muons
  Int_t nGoodMuons = 0;
  Bool_t bAtLeastOneMuon = kFALSE;
  if (fFillMuons || fRequireMuon) {
    for (Int_t iMu = 0; iMu<fAOD->GetNumberOfTracks(); iMu++) {
      AliAODTrack* muonTrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iMu));
      if(!muonTrack) AliFatal("Not a standard AOD");
      if (muonTrack) {
        if (!IsGoodMUONtrack(*muonTrack))
          continue;
        Double_t ptMu   = muonTrack->Pt();
        if (ptMu>fReqPtLo && ptMu<fReqPtHi)
          bAtLeastOneMuon = kTRUE;
        Double_t etaMu  = muonTrack->Eta();
        if (IsTrigger(etaMu,ptMu)||IsAssociated(etaMu,ptMu))
          nGoodMuons++;
      }
    }
  }
  if (fRequireMuon && !bAtLeastOneMuon)
    return;
  
  // Loop twice: 1. Count sel. tracks. 2. Fill vector.
  if (fTracksName.IsNull()) {
    Int_t nTrax = fAOD->GetNumberOfTracks();
    Int_t nSelTrax = 0;

    if (fVerbosity > 2)
      AliInfo(Form("%d tracks in event",nTrax));

    // Loop 1.
    for (Int_t i = 0; i < nTrax; ++i) {
      AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!aodtrack) AliFatal("Not a standard AOD");
      if (!aodtrack) {
        AliError(Form("Couldn't get AOD track %d\n", i));
        continue;
      }
      // See $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C
      UInt_t tpcOnly = 1 << 7;
      Bool_t trkOK = aodtrack->TestFilterBit(tpcOnly);
      if (!trkOK)
        continue;
      Double_t pt = aodtrack->Pt();
      Double_t eta = aodtrack->Eta();
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        nSelTrax++;
    }

    if (miniEvt)
      miniEvt->reserve(nSelTrax);
    else {
      AliError("!miniEvt");
      return;
    }
  
    // Loop 2.  
    for (Int_t i = 0; i < nTrax; ++i) {
      AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!aodtrack) AliFatal("Not a standard AOD");
      if (!aodtrack) {
        AliError(Form("Couldn't get AOD track %d\n", i));
        continue;
      }
    
      // See $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C
      UInt_t tpcOnly = 1 << 7;
      Bool_t trkOK = aodtrack->TestFilterBit(tpcOnly);
      if (!trkOK)
        continue;
      Double_t pt = aodtrack->Pt();
      Double_t eta  = aodtrack->Eta();
      Double_t phi  = aodtrack->Phi();
      Int_t    sign = aodtrack->Charge() > 0 ? 1 : -1;
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        miniEvt->push_back(AliMiniTrack(pt, eta, phi, sign));
    }
  } else {
    TList *list = InputEvent()->GetList();
    TClonesArray *tcaTracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));

    if (!tcaTracks){
      AliError("Ptr to tcaTracks zero");
      return;
    }

    const Int_t ntracks = tcaTracks->GetEntries();
    Int_t nGoodTracks = 0;
    // count good tracks
    for (Int_t itrack = 0; itrack < ntracks; itrack++) {
      AliVTrack *vtrack = static_cast<AliVTrack*>(tcaTracks->At(itrack));
      if (!vtrack) {
        AliError(Form("ERROR: Could not retrieve vtrack %d",itrack));
        continue;
      }
      Double_t pt   = vtrack->Pt();
      Double_t eta  = vtrack->Eta();
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        nGoodTracks++;
    }
    if (miniEvt)
      miniEvt->reserve(nGoodTracks);
    else {
      AliError("Ptr to miniEvt zero");
      return;
    }
    // fill good tracks into minievent
    for (Int_t itrack = 0; itrack < ntracks; itrack++) {
      AliVTrack *vtrack = static_cast<AliVTrack*>(tcaTracks->At(itrack));
      if (!vtrack) {
        AliError(Form("ERROR: Could not retrieve vtrack %d",itrack));
        continue;
      }
      Double_t pt   = vtrack->Pt();
      Double_t eta  = vtrack->Eta();
      Double_t phi  = vtrack->Phi();
      Int_t    sign = vtrack->Charge() > 0 ? 1 : -1;
      if (IsTrigger(eta,pt)||IsAssociated(eta,pt))
        miniEvt->push_back(AliMiniTrack(pt, eta, phi, sign));
    }
  }

  // fill muons into minievent
  if (fFillMuons) {
    miniEvt->reserve(miniEvt->size()+nGoodMuons);
    // fill them into the mini event
    for (Int_t iMu = 0; iMu<fAOD->GetNumberOfTracks(); iMu++) {
      AliAODTrack* muonTrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iMu));
      if(!muonTrack) AliFatal("Not a standard AOD");
      if (muonTrack) {
        if (!IsGoodMUONtrack(*muonTrack)) 
          continue;
        Double_t ptMu   = muonTrack->Pt();
        Double_t etaMu  = muonTrack->Eta();
        Double_t phiMu  = muonTrack->Phi();
        Int_t    signMu = muonTrack->Charge() > 0 ? 1 : -1;
        if (IsTrigger(etaMu,ptMu)||IsAssociated(etaMu,ptMu))
          miniEvt->push_back(AliMiniTrack(ptMu, etaMu, phiMu, signMu));
      }
    }
  }

}

//________________________________________________________________________
Double_t AliDhcTask::DeltaPhi(Double_t phia, Double_t phib,
			      Double_t rangeMin, Double_t rangeMax) const
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();
  
  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;
  
  return dphi;
}

//________________________________________________________________________
Int_t AliDhcTask::Correlate(const MiniEvent &evt1, const MiniEvent &evt2, Int_t pairing)
{
  // Triggered angular correlations. If pairing is kSameEvt, particles
  // within evt1 are correlated. If kDiffEvt, correlate triggers from
  // evt1 with partners from evt2.

  Int_t cbin = fHCent->FindBin(fCentrality);
  if (fHCent->IsBinOverflow(cbin) ||
      fHCent->IsBinUnderflow(cbin))
    return 0;

  Int_t zbin = fHZvtx->FindBin(fZVertex);
  if (fHZvtx->IsBinOverflow(zbin) ||
      fHZvtx->IsBinUnderflow(zbin))
    return 0;

  Int_t iMax = evt1.size();
  Int_t jMax = evt2.size();

  TH2  **hist = fHMs;
  if (pairing == kSameEvt) {
    hist = fHSs;
    fHCent->Fill(fCentrality);
    fHZvtx->Fill(fZVertex);
  }

  Int_t nZvtx = fHZvtx->GetNbinsX();
  Int_t nPtTrig = fHPtTrg->GetNbinsX();
  Int_t nPtAssc = fHPtAss->GetNbinsX();

  Int_t globIndex = (cbin-1)*nZvtx*nPtTrig*nPtAssc+(zbin-1)*nPtTrig*nPtAssc;
  Int_t ptindex   = (Int_t)fIndex->Eval(1,1,zbin,cbin);
  Int_t mindex    = (Int_t)fIndex->Eval(1,1,1,cbin);


  fHPtTrgEvt->Reset();
  for (Int_t i=0; i<iMax; ++i) {
    const AliMiniTrack &a(evt1.at(i));
    Float_t pta  = a.Pt();
    fHPtTrgEvt->Fill(pta);
    if (pairing == kSameEvt) {
      fHPts[ptindex]->Fill(pta);
    }
  }

  Bool_t bCountTrg; // only count the trigger if an associated particle is found

  for (Int_t i=0; i<iMax; ++i) {
    // Trigger particles
    const AliMiniTrack &a(evt1.at(i));

    Float_t pta  = a.Pt();
    Float_t etaa = a.Eta();
    Float_t phia = a.Phi();
    Int_t   sgna = a.Sign();
    
    // brief intermezzo in the trigger particle loop: do associated particle QA here in order to avoid double counting
    if (pairing == kSameEvt) {
      if (IsAssociated(etaa,pta)) {
        Double_t aQAWght = 1.0;
        if (fHEffA) {
          const Int_t nEffDimA = fHEffA->GetNdimensions();
          Int_t effBinA[nEffDimA];
          effBinA[0] = fHEffA->GetAxis(0)->FindBin(etaa);
          effBinA[1] = fHEffA->GetAxis(1)->FindBin(pta);
          effBinA[2] = fHEffA->GetAxis(2)->FindBin(fCentrality);
          effBinA[3] = fHEffA->GetAxis(3)->FindBin(fZVertex);
          if (nEffDimA>4) {
            effBinA[4] = fHEffA->GetAxis(4)->FindBin(phia);
            if (nEffDimA>5) {
              effBinA[5] = fHEffA->GetAxis(5)->FindBin(sgna);
            }
          }
          aQAWght = fHEffA->GetBinContent(effBinA);
        }
        // fill every associated particle once unweighted, once weighted
        if (sgna>0.0) {
          fHQAAp->Fill(pta,etaa,phia);
          fHQAApCorr->Fill(pta,etaa,phia,aQAWght);
        } else {
          fHQAAm->Fill(pta,etaa,phia);
          fHQAAmCorr->Fill(pta,etaa,phia,aQAWght);
        }
        fHPtCentA->Fill(pta,fCentrality);
      }
    }
    
    // back to triggers
    if (!IsTrigger(etaa,pta))
      continue;

    Int_t abin = fHPtAss->FindBin(pta);

    // efficiency weighting
    Double_t effWtT = 1.0;
    if (fHEffT) {
      // trigger particle
      const Int_t nEffDimT = fHEffT->GetNdimensions();
      Int_t effBinT[nEffDimT];
      effBinT[0] = fHEffT->GetAxis(0)->FindBin(etaa);
      effBinT[1] = fHEffT->GetAxis(1)->FindBin(pta);
      effBinT[2] = fHEffT->GetAxis(2)->FindBin(fCentrality);
      effBinT[3] = fHEffT->GetAxis(3)->FindBin(fZVertex);
      if (nEffDimT>4) {
        effBinT[4] = fHEffT->GetAxis(4)->FindBin(phia);
        if (nEffDimT>5) {
          effBinT[5] = fHEffT->GetAxis(5)->FindBin(sgna);
        }
      }
      effWtT = fHEffT->GetBinContent(effBinT);
    }
    
    if (pairing == kSameEvt) {
      fHTrk->Fill(phia,etaa);
      if (sgna>0.0) {
        fHQATp->Fill(pta,etaa,phia);
        fHQATpCorr->Fill(pta,etaa,phia,effWtT);
      } else {
        fHQATm->Fill(pta,etaa,phia);
        fHQATmCorr->Fill(pta,etaa,phia,effWtT);
      }
      fHPtCentT->Fill(pta,fCentrality);
      fHPtTrg->Fill(pta);
      fHPtTrgNorm1S->Fill(pta,fCentrality,fZVertex,effWtT);
    } else {
      fHPtTrgNorm1M->Fill(pta,fCentrality,fZVertex,effWtT);
    }
    bCountTrg = kFALSE;

    for (Int_t j=0; j<jMax; ++j) {
      // Associated particles
      if (pairing == kSameEvt && i==j)
        continue;

      const AliMiniTrack &b(evt2.at(j));
      
      Float_t ptb  = b.Pt();
      Float_t etab = b.Eta();
      Float_t phib = b.Phi();
      Int_t   sgnb = b.Sign();
      if (fPtTACrit&&(pta < ptb)) {
        continue;
      }

      if (!IsAssociated(etab,ptb))
        continue;

      Int_t bbin = fHPtAss->FindBin(ptb);

      Float_t dphi = DeltaPhi(phia, phib);
      Float_t deta = etaa - etab;
      // invariant mass under mu mu assumption
      Float_t muMass2 = 0.1056583715*0.1056583715;
      Float_t ea2  = muMass2 + pta*pta*TMath::CosH(etaa)*TMath::CosH(etaa);
      Float_t eb2  = muMass2 + ptb*ptb*TMath::CosH(etab)*TMath::CosH(etab);
      Float_t papb = pta*TMath::Cos(phia)*ptb*TMath::Cos(phib) +
                     pta*TMath::Sin(phia)*ptb*TMath::Sin(phib) +
                     pta*TMath::SinH(etaa)*ptb*TMath::SinH(etab);
      Float_t mass = TMath::Sqrt(2.0*(muMass2 + TMath::Sqrt(ea2*eb2) - papb));
      Int_t q2 = sgna + sgnb;
      if ((q2==0) && fDoMassCut) {
        if (mass>3.0 && mass<3.2)
          continue;
        if (mass>9.2&&mass<9.8)
          continue;
      }

      Int_t index = globIndex+(abin-1)*nPtAssc+(bbin-1);
      Double_t weight = 1.0;
      // number of trigger weight event-by-event (CMS method)
      if (fDoWeights) {
        Double_t nTrgs = fHPtTrgEvt->GetBinContent(abin);
        if (nTrgs==0.0) {
          AliError(Form("Trying to work with number of triggers weight = %f",nTrgs));
          return 0;
        }
        weight *= 1./nTrgs;
      }

      // efficiency weighting
      if (fHEffT) {
        // trigger particle
        weight *= effWtT;
      }

      if (fHEffA) {
        // associated particle
        const Int_t nEffDimA = fHEffA->GetNdimensions();
        Int_t effBinA[nEffDimA];
        effBinA[0] = fHEffA->GetAxis(0)->FindBin(etab);
        effBinA[1] = fHEffA->GetAxis(1)->FindBin(ptb);
        effBinA[2] = fHEffA->GetAxis(2)->FindBin(fCentrality);
        effBinA[3] = fHEffA->GetAxis(3)->FindBin(fZVertex);
        if (nEffDimA>4) {
          effBinA[4] = fHEffA->GetAxis(4)->FindBin(phib);
          if (nEffDimA>5) {
            effBinA[5] = fHEffA->GetAxis(5)->FindBin(sgnb);
          }
        }
        weight *= fHEffA->GetBinContent(effBinA);
      }

      if (hist[index]) { // check if this histogram exists, relevant in the fPtTACrit==kFALSE case
        hist[index]->Fill(dphi,deta,weight);
        bCountTrg = kTRUE;
        if (pairing == kSameEvt) {
          fHPtAss->Fill(ptb); // fill every associated particle every time it is used
          if (q2==0)
            fHSMass[mindex]->Fill(mass);
        } else {
          if (q2==0)
            fHMMass[mindex]->Fill(mass);
        }
      }
    }
    if (bCountTrg) {
      if (pairing == kSameEvt) {
        fHPtTrgNorm2S->Fill(pta,fCentrality,fZVertex,effWtT);
      } else {
        fHPtTrgNorm2M->Fill(pta,fCentrality,fZVertex,effWtT);
      }
    }
  }

  return 1;
}

//________________________________________________________________________
void AliDhcTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}

//________________________________________________________________________
Bool_t AliDhcTask::VertexOk() const
{
  // Modified from AliAnalyseLeadingTrackUE::VertexSelection()
  
  Int_t nContributors  = 0;
  Double_t zVertex     = 999;
  TString name("");

  Int_t runno = InputEvent()->GetRunNumber();
  if (runno>=176326 && runno<=197692) { // year 12 and 13
    if (!fUtils->IsVertexSelected2013pA(InputEvent())) 
      return 0;
  }

  if (fESD) {
    const AliESDVertex* vtx = fESD->GetPrimaryVertex();
    if (!vtx)
      return 0;
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  } else {
    if (fAOD->GetNumberOfVertices() < 1)
      return 0;
    const AliAODVertex* vtx = fAOD->GetPrimaryVertex();
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
  
  // Reject if TPC-only vertex
  if (name.CompareTo("TPCVertex")==0)
    return kFALSE;
  
  // Check # contributors and range...
  if( nContributors < 1 || TMath::Abs(zVertex) > fZVtxMax ) {
    return kFALSE;
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliDhcTask::IsGoodMUONtrack(AliESDMuonTrack &track)
{
  // Applying track cuts for MUON tracks

  if (!track.ContainTrackerData()) 
    return kFALSE;
  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  if ((thetaTrackAbsEnd < 2.) || (thetaTrackAbsEnd > 10.)) 
    return kFALSE;
  Double_t eta = track.Eta();
  if ((eta < -4.) || (eta > -2.5))
    return kFALSE;
  if (fTriggerMatch) {
    if (!track.ContainTriggerData()) 
      return kFALSE;
    if (track.GetMatchTrigger() < 0.5) 
      return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliDhcTask::IsGoodMUONtrack(AliAODTrack &track)
{
  // Applying track cuts for MUON tracks

  if (!track.IsMuonTrack()) 
    return kFALSE;
  Double_t dThetaAbs = TMath::ATan(track.GetRAtAbsorberEnd()/505.)* TMath::RadToDeg();
  if ((dThetaAbs<2.) || (dThetaAbs>10.)) 
    return kFALSE;
  Double_t dEta = track.Eta();
  if ((dEta<-4.) || (dEta>-2.5)) 
    return kFALSE;
  if (fTriggerMatch) {
    if (track.GetMatchTrigger()<0.5) 
      return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliDhcTask::IsTrigger(Double_t eta, Double_t pt)
{
  // is in trigger eta,pt range?
  Int_t bin = fHPtTrg->FindBin(pt);
  if (fHPtTrg->IsBinOverflow(bin) || fHPtTrg->IsBinUnderflow(bin))
    return kFALSE;
  
  if (eta<fEtaTLo || eta>fEtaTHi)
    return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliDhcTask::IsAssociated(Double_t eta, Double_t pt)
{
  // is in associated eta,pt range?
  Int_t bin = fHPtAss->FindBin(pt);
  if (fHPtAss->IsBinOverflow(bin) || fHPtAss->IsBinUnderflow(bin))
    return kFALSE;
  
  if (eta<fEtaALo || eta>fEtaAHi)
    return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
AliDhcTask::~AliDhcTask()
{
  //Destructor
  if (fPoolMgr) 
    delete fPoolMgr;
}
