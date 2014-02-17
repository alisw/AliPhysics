//
// Jet mass analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalJetMass.h"

ClassImp(AliAnalysisTaskEmcalJetMass)

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::AliAnalysisTaskEmcalJetMass() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetMass", kTRUE),
  fContainerBase(0),
  fMinFractionShared(0),
  fJetMassAvg(0),
  fh2PtJet1VsLeadPtAllSel(0),
  fh2PtJet1VsLeadPtTagged(0),
  fh2PtVsMassJet1All(0),
  fh2PtVsMassJet1Tagged(0),
  fpPtVsMassJet1All(0),
  fpPtVsMassJet1Tagged(0),
  fh2MassVsAreaJet1All(0),
  fh2MassVsAreaJet1Tagged(0),
  fh2MassVsNConstJet1All(0),
  fh2MassVsNConstJet1Tagged(0),
  fh2EtMassOverEtRSq(0)
{
  // Default constructor.

  fh2PtJet1VsLeadPtAllSel      = new TH2F*[fNcentBins];
  fh2PtJet1VsLeadPtTagged      = new TH2F*[fNcentBins];
  fh2PtVsMassJet1All           = new TH2F*[fNcentBins];
  fh2PtVsMassJet1Tagged        = new TH2F*[fNcentBins];
  fpPtVsMassJet1All            = new TProfile*[fNcentBins];
  fpPtVsMassJet1Tagged         = new TProfile*[fNcentBins];
  fh2MassVsAreaJet1All         = new TH2F*[fNcentBins];
  fh2MassVsAreaJet1Tagged      = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1All       = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1Tagged    = new TH2F*[fNcentBins];
  fh2EtMassOverEtRSq           = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1VsLeadPtAllSel[i]     = 0;
    fh2PtJet1VsLeadPtTagged[i]     = 0;
    fh2PtVsMassJet1All[i]          = 0;
    fh2PtVsMassJet1Tagged[i]       = 0;
    fpPtVsMassJet1All[i]           = 0;
    fpPtVsMassJet1Tagged[i]        = 0;
    fh2MassVsAreaJet1All[i]        = 0;
    fh2MassVsAreaJet1Tagged[i]     = 0;
    fh2MassVsNConstJet1All[i]      = 0;
    fh2MassVsNConstJet1Tagged[i]   = 0;
    fh2EtMassOverEtRSq[i]          = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::AliAnalysisTaskEmcalJetMass(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinFractionShared(0),
  fJetMassAvg(0),
  fh2PtJet1VsLeadPtAllSel(0),
  fh2PtJet1VsLeadPtTagged(0),
  fh2PtVsMassJet1All(0),
  fh2PtVsMassJet1Tagged(0),
  fpPtVsMassJet1All(0),
  fpPtVsMassJet1Tagged(0),
  fh2MassVsAreaJet1All(0),
  fh2MassVsAreaJet1Tagged(0),
  fh2MassVsNConstJet1All(0),
  fh2MassVsNConstJet1Tagged(0),
  fh2EtMassOverEtRSq(0)
{
  // Standard constructor.

  fh2PtJet1VsLeadPtAllSel      = new TH2F*[fNcentBins];
  fh2PtJet1VsLeadPtTagged      = new TH2F*[fNcentBins];
  fh2PtVsMassJet1All           = new TH2F*[fNcentBins];
  fh2PtVsMassJet1Tagged        = new TH2F*[fNcentBins];
  fpPtVsMassJet1All            = new TProfile*[fNcentBins];
  fpPtVsMassJet1Tagged         = new TProfile*[fNcentBins];
  fh2MassVsAreaJet1All         = new TH2F*[fNcentBins];
  fh2MassVsAreaJet1Tagged      = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1All       = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1Tagged    = new TH2F*[fNcentBins];
  fh2EtMassOverEtRSq           = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1VsLeadPtAllSel[i]     = 0;
    fh2PtJet1VsLeadPtTagged[i]     = 0;
    fh2PtVsMassJet1All[i]          = 0;
    fh2PtVsMassJet1Tagged[i]       = 0;
    fpPtVsMassJet1All[i]           = 0;
    fpPtVsMassJet1Tagged[i]        = 0;
    fh2MassVsAreaJet1All[i]        = 0;
    fh2MassVsAreaJet1Tagged[i]     = 0;
    fh2MassVsNConstJet1All[i]      = 0;
    fh2MassVsNConstJet1Tagged[i]   = 0;
    fh2EtMassOverEtRSq[i]          = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::~AliAnalysisTaskEmcalJetMass()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMass::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt  = 250;
  const Double_t minPt = -50.;
  const Double_t maxPt = 200.;

  const Int_t nBinsM  = 150;
  const Double_t minM = -50.;
  const Double_t maxM = 100.;

  const Int_t nBinsArea = 100;
  const Double_t minArea = 0.;
  const Double_t maxArea = 1.;

  const Int_t nBinsNConst = 100;
  const Double_t minNConst = 0.;
  const Double_t maxNConst = 500.;

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = TString::Format("fh2PtJet1VsLeadPtAllSel_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
    fh2PtJet1VsLeadPtAllSel[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,20,0.,20.);
    fOutput->Add(fh2PtJet1VsLeadPtAllSel[i]);

    histName = TString::Format("fh2PtJet1VsLeadPtTagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
    fh2PtJet1VsLeadPtTagged[i] =  new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,20,0.,20.);
    fOutput->Add(fh2PtJet1VsLeadPtTagged[i]);

    histName = TString::Format("fh2PtVsMassJet1All_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
    fh2PtVsMassJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtVsMassJet1All[i]);

    histName = TString::Format("fh2PtVsMassJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
    fh2PtVsMassJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtVsMassJet1Tagged[i]);

    histName = TString::Format("fpPtVsMassJet1All_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1All[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsM,minM,maxM);
    fOutput->Add(fpPtVsMassJet1All[i]);

    histName = TString::Format("fpPtVsMassJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1Tagged[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsM,minM,maxM);
    fOutput->Add(fpPtVsMassJet1Tagged[i]);

    histName = TString::Format("fh2MassVsAreaJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1All[i]);

    histName = TString::Format("fh2MassVsAreaJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1Tagged[i]);

    histName = TString::Format("fh2MassVsNConstJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{N}_{constituents}",histName.Data());
    fh2MassVsNConstJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1All[i]);

    histName = TString::Format("fh2MassVsNConstJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{N}_{constituents}",histName.Data());
    fh2MassVsNConstJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1Tagged[i]);

    histName = TString::Format("fh2EtMassOverEtRSq_%d",i);
    histTitle = TString::Format("%s;#it{E}_{T};(#it{M}/(#it{E}_{T}#it{R}))^{2}",histName.Data());
    fh2EtMassOverEtRSq[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,0.,1.);
    fOutput->Add(fh2EtMassOverEtRSq[i]);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMass::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMass::FillHistograms()
{
  // Fill histograms.

  AliInfo(Form("%s",GetName()));

  AliEmcalJet* jet1 = NULL;

  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      Double_t fraction = jetCont->GetFractionSharedPt(jet1);
      if(fMinFractionShared>0. && fraction<fMinFractionShared) continue;

      Double_t ptJet1 = jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();//jetCont->GetJetPtCorr(jetCont->GetCurrentID());//
      fh2PtJet1VsLeadPtAllSel[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());
      fh2PtVsMassJet1All[fCentBin]->Fill(ptJet1,jet1->M());
      fpPtVsMassJet1All[fCentBin]->Fill(ptJet1,jet1->M());
      fh2MassVsAreaJet1All[fCentBin]->Fill(jet1->M(),jet1->Area());
      fh2MassVsNConstJet1All[fCentBin]->Fill(jet1->M(),jet1->GetNumberOfConstituents());
      
      if(jet1->GetTagStatus()<1 || !jet1->GetTaggedJet())
	continue;

      fh2PtJet1VsLeadPtTagged[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());
      fh2PtVsMassJet1Tagged[fCentBin]->Fill(ptJet1,jet1->M());
      fpPtVsMassJet1Tagged[fCentBin]->Fill(ptJet1,jet1->M());
      fh2MassVsAreaJet1Tagged[fCentBin]->Fill(jet1->M(),jet1->Area());
      fh2MassVsNConstJet1Tagged[fCentBin]->Fill(jet1->M(),jet1->GetNumberOfConstituents());
      
      Double_t Et2 = jet1->M()*jet1->M() + jet1->Pt()*jet1->Pt();
      Double_t Et = 0.;    Double_t massOverEtR = 0.;
      if(Et2>0.) Et = TMath::Sqrt(Et2);
      if((Et*jetCont->GetJetRadius())>0.) 
	massOverEtR = jet1->M()/(Et*jetCont->GetJetRadius());
      fh2EtMassOverEtRSq[fCentBin]->Fill(Et,massOverEtR*massOverEtR);

      //check if matched gen level jet exists and do analysis
      AliEmcalJet *jetGen = jet1->ClosestJet();

      Printf("AA jet");
      Printf("jet 4-vector: %f,%f,%f,%f",jet1->Px(),jet1->Py(),jet1->Pz(),jet1->E());
      Printf("pT: %f  pTcorr: %f  M: %f",jet1->Pt(),ptJet1,jet1->M());
      Printf("eta: %f  phi: %f",jet1->Eta(),jet1->Phi());
      if(jetGen) {
	Printf("gen jet");
	Printf("jet 4-vector: %f,%f,%f,%f",jetGen->Px(),jetGen->Py(),jetGen->Pz(),jetGen->E());
	Printf("pT: %f  M: %f",jetGen->Pt(),jetGen->M());
	Printf("eta: %f  phi: %f",jetGen->Eta(),jetGen->Phi());
      }

    }
  }

  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMass::GetJetMass(AliEmcalJet *jet) {

  Double_t deltaM = jet->M() - fJetMassAvg;
  Double_t scale = deltaM / jet->M();

  Double_t mt2 = jet->E()*jet->E() - jet->Pt()*jet->Pt();
  Double_t et2 = jet->M()*jet->M() + jet->Pt()*jet->Pt();

  Double_t pxScale = jet->Px()*scale;
  Double_t pyScale = jet->Py()*scale;
  Double_t pzScale = jet->Pz()*scale;
  
  // Printf("scaled jet 4-vector: %f-%f-%f-%f",pxScale,pyScale,pzScale);

  return deltaM;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMass::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;

}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetMass::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

