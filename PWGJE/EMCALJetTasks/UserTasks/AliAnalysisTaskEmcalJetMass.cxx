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
  fh2PtJet1VsLeadPtAllSel(0),
  fh2PtJet1VsLeadPtTagged(0),
  fh2PtVsMassJet1All(0),
  fh2PtVsMassJet1Tagged(0),
  fpPtVsMassJet1All(0),
  fpPtVsMassJet1Tagged(0),
  fh2MassVsAreaJet1All(0),
  fh2MassVsAreaJet1Tagged(0)
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

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1VsLeadPtAllSel[i]     = 0;
    fh2PtJet1VsLeadPtTagged[i]     = 0;
    fh2PtVsMassJet1All[i]          = 0;
    fh2PtVsMassJet1Tagged[i]       = 0;
    fpPtVsMassJet1All[i]           = 0;
    fpPtVsMassJet1Tagged[i]        = 0;
    fh2MassVsAreaJet1All[i]        = 0;
    fh2MassVsAreaJet1Tagged[i]     = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::AliAnalysisTaskEmcalJetMass(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fh2PtJet1VsLeadPtAllSel(0),
  fh2PtJet1VsLeadPtTagged(0),
  fh2PtVsMassJet1All(0),
  fh2PtVsMassJet1Tagged(0),
  fpPtVsMassJet1All(0),
  fpPtVsMassJet1Tagged(0),
  fh2MassVsAreaJet1All(0),
  fh2MassVsAreaJet1Tagged(0)
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

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1VsLeadPtAllSel[i]     = 0;
    fh2PtJet1VsLeadPtTagged[i]     = 0;
    fh2PtVsMassJet1All[i]          = 0;
    fh2PtVsMassJet1Tagged[i]       = 0;
    fpPtVsMassJet1All[i]           = 0;
    fpPtVsMassJet1Tagged[i]        = 0;
    fh2MassVsAreaJet1All[i]        = 0;
    fh2MassVsAreaJet1Tagged[i]     = 0;
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

  const Int_t nBinsArea = 100;
  const Double_t minArea = 0.;
  const Double_t maxArea = 1.;

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
    fh2PtVsMassJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtVsMassJet1All[i]);

    histName = TString::Format("fh2PtVsMassJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
    fh2PtVsMassJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtVsMassJet1Tagged[i]);

    histName = TString::Format("fpPtVsMassJet1All_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1All[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassJet1All[i]);

    histName = TString::Format("fpPtVsMassJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1Tagged[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassJet1Tagged[i]);

    histName = TString::Format("fh2MassVsAreaJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1All[i]);

    histName = TString::Format("fh2MassVsAreaJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1Tagged[i]);

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

  for(int i = 0; i < GetNJets(fContainerBase);++i) {
    AliEmcalJet *jet1 = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(i, fContainerBase));
    if(!jet1) continue;

    Double_t ptJet1 =  jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
    fh2PtJet1VsLeadPtAllSel[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());
    fh2PtVsMassJet1All[fCentBin]->Fill(ptJet1,jet1->M());
    fpPtVsMassJet1All[fCentBin]->Fill(ptJet1,jet1->M());
    fh2MassVsAreaJet1All[fCentBin]->Fill(jet1->M(),jet1->Area());    

    if(jet1->GetTagStatus()<1)
      continue;

    AliEmcalJet *jet2 = jet1->GetTaggedJet();
    if(!jet2) continue;
    
    fh2PtJet1VsLeadPtTagged[fCentBin]->Fill(ptJet1,jet1->MaxTrackPt());
    fh2PtVsMassJet1Tagged[fCentBin]->Fill(ptJet1,jet1->M());
    fpPtVsMassJet1Tagged[fCentBin]->Fill(ptJet1,jet1->M());
    fh2MassVsAreaJet1Tagged[fCentBin]->Fill(jet1->M(),jet1->Area());    
  
  }

  return kTRUE;

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

