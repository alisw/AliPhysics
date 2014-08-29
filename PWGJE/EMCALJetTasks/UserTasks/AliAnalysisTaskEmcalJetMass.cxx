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
  fJetMassType(kRaw),
  fh3PtJet1VsMassVsLeadPtAllSel(0),
  fh3PtJet1VsMassVsLeadPtTagged(0),
  fh3PtJet1VsMassVsLeadPtTaggedMatch(0),
  fpPtVsMassJet1All(0),
  fpPtVsMassJet1Tagged(0),
  fpPtVsMassJet1TaggedMatch(0),
  fh2MassVsAreaJet1All(0),
  fh2MassVsAreaJet1Tagged(0),
  fh2MassVsAreaJet1TaggedMatch(0),
  fh2MassVsNConstJet1All(0),
  fh2MassVsNConstJet1Tagged(0),
  fh2MassVsNConstJet1TaggedMatch(0),
  fh2EtMassOverEtRSq(0)
{
  // Default constructor.

  fh3PtJet1VsMassVsLeadPtAllSel        = new TH3F*[fNcentBins];
  fh3PtJet1VsMassVsLeadPtTagged        = new TH3F*[fNcentBins];
  fh3PtJet1VsMassVsLeadPtTaggedMatch   = new TH3F*[fNcentBins];
  fpPtVsMassJet1All                    = new TProfile*[fNcentBins];
  fpPtVsMassJet1Tagged                 = new TProfile*[fNcentBins];
  fpPtVsMassJet1TaggedMatch            = new TProfile*[fNcentBins];
  fh2MassVsAreaJet1All                 = new TH2F*[fNcentBins];
  fh2MassVsAreaJet1Tagged              = new TH2F*[fNcentBins];
  fh2MassVsAreaJet1TaggedMatch         = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1All               = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1Tagged            = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1TaggedMatch       = new TH2F*[fNcentBins];
  fh2EtMassOverEtRSq                   = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtJet1VsMassVsLeadPtAllSel[i]        = 0;
    fh3PtJet1VsMassVsLeadPtTagged[i]        = 0;
    fh3PtJet1VsMassVsLeadPtTaggedMatch[i]   = 0;
    fpPtVsMassJet1All[i]                    = 0;
    fpPtVsMassJet1Tagged[i]                 = 0;
    fpPtVsMassJet1TaggedMatch[i]            = 0;
    fh2MassVsAreaJet1All[i]                 = 0;
    fh2MassVsAreaJet1Tagged[i]              = 0;
    fh2MassVsAreaJet1TaggedMatch[i]         = 0;
    fh2MassVsNConstJet1All[i]               = 0;
    fh2MassVsNConstJet1Tagged[i]            = 0;
    fh2MassVsNConstJet1TaggedMatch[i]       = 0;
    fh2EtMassOverEtRSq[i]                   = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::AliAnalysisTaskEmcalJetMass(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinFractionShared(0),
  fJetMassType(kRaw),
  fh3PtJet1VsMassVsLeadPtAllSel(0),
  fh3PtJet1VsMassVsLeadPtTagged(0),
  fh3PtJet1VsMassVsLeadPtTaggedMatch(0),
  fpPtVsMassJet1All(0),
  fpPtVsMassJet1Tagged(0),
  fpPtVsMassJet1TaggedMatch(0),
  fh2MassVsAreaJet1All(0),
  fh2MassVsAreaJet1Tagged(0),
  fh2MassVsAreaJet1TaggedMatch(0),
  fh2MassVsNConstJet1All(0),
  fh2MassVsNConstJet1Tagged(0),
  fh2MassVsNConstJet1TaggedMatch(0),
  fh2EtMassOverEtRSq(0)
{
  // Standard constructor.

  fh3PtJet1VsMassVsLeadPtAllSel        = new TH3F*[fNcentBins];
  fh3PtJet1VsMassVsLeadPtTagged        = new TH3F*[fNcentBins];
  fh3PtJet1VsMassVsLeadPtTaggedMatch   = new TH3F*[fNcentBins];
  fpPtVsMassJet1All                    = new TProfile*[fNcentBins];
  fpPtVsMassJet1Tagged                 = new TProfile*[fNcentBins];
  fpPtVsMassJet1TaggedMatch            = new TProfile*[fNcentBins];
  fh2MassVsAreaJet1All                 = new TH2F*[fNcentBins];
  fh2MassVsAreaJet1Tagged              = new TH2F*[fNcentBins];
  fh2MassVsAreaJet1TaggedMatch         = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1All               = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1Tagged            = new TH2F*[fNcentBins];
  fh2MassVsNConstJet1TaggedMatch       = new TH2F*[fNcentBins];
  fh2EtMassOverEtRSq                   = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtJet1VsMassVsLeadPtAllSel[i]        = 0;
    fh3PtJet1VsMassVsLeadPtTagged[i]        = 0;
    fh3PtJet1VsMassVsLeadPtTaggedMatch[i]   = 0;
    fpPtVsMassJet1All[i]                    = 0;
    fpPtVsMassJet1Tagged[i]                 = 0;
    fpPtVsMassJet1TaggedMatch[i]            = 0;
    fh2MassVsAreaJet1All[i]                 = 0;
    fh2MassVsAreaJet1Tagged[i]              = 0;
    fh2MassVsAreaJet1TaggedMatch[i]         = 0;
    fh2MassVsNConstJet1All[i]               = 0;
    fh2MassVsNConstJet1Tagged[i]            = 0;
    fh2MassVsNConstJet1TaggedMatch[i]       = 0;
    fh2EtMassOverEtRSq[i]                   = 0;
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

  const Int_t nBinsPt  = 200;
  const Double_t minPt = -50.;
  const Double_t maxPt = 150.;

  const Int_t nBinsM  = 100;
  const Double_t minM = -20.;
  const Double_t maxM = 80.;

  const Int_t nBinsArea = 50;
  const Double_t minArea = 0.;
  const Double_t maxArea = 1.;

  const Int_t nBinsNConst = 100;
  const Double_t minNConst = 0.;
  const Double_t maxNConst = 500.;

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = TString::Format("fh3PtJet1VsMassVsLeadPtAllSel_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,lead trk}",histName.Data());
    fh3PtJet1VsMassVsLeadPtAllSel[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,20,0.,20.);
    fOutput->Add(fh3PtJet1VsMassVsLeadPtAllSel[i]);

    histName = TString::Format("fh3PtJet1VsMassVsLeadPtTagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,lead trk}",histName.Data());
    fh3PtJet1VsMassVsLeadPtTagged[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,20,0.,20.);
    fOutput->Add(fh3PtJet1VsMassVsLeadPtTagged[i]);

    histName = TString::Format("fh3PtJet1VsMassVsLeadPtTaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{p}_{T,lead trk}",histName.Data());
    fh3PtJet1VsMassVsLeadPtTaggedMatch[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,20,0.,20.);
    fOutput->Add(fh3PtJet1VsMassVsLeadPtTaggedMatch[i]);

    histName = TString::Format("fpPtVsMassJet1All_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1All[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassJet1All[i]);

    histName = TString::Format("fpPtVsMassJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1Tagged[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassJet1Tagged[i]);

    histName = TString::Format("fpPtVsMassJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}",histName.Data());
    fpPtVsMassJet1TaggedMatch[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassJet1TaggedMatch[i]);

    histName = TString::Format("fh2MassVsAreaJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1All[i]);

    histName = TString::Format("fh2MassVsAreaJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1Tagged[i]);

    histName = TString::Format("fh2MassVsAreaJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{A}",histName.Data());
    fh2MassVsAreaJet1TaggedMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2MassVsAreaJet1TaggedMatch[i]);

    histName = TString::Format("fh2MassVsNConstJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{N}_{constituents}",histName.Data());
    fh2MassVsNConstJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1All[i]);

    histName = TString::Format("fh2MassVsNConstJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{N}_{constituents}",histName.Data());
    fh2MassVsNConstJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1Tagged[i]);

    histName = TString::Format("fh2MassVsNConstJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{N}_{constituents}",histName.Data());
    fh2MassVsNConstJet1TaggedMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1TaggedMatch[i]);

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

  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {

      Double_t ptJet1 = jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
      Double_t mJet1 = GetJetMass(jet1);
      fh3PtJet1VsMassVsLeadPtAllSel[fCentBin]->Fill(ptJet1,mJet1,jet1->MaxTrackPt());
      fpPtVsMassJet1All[fCentBin]->Fill(ptJet1,mJet1);
      fh2MassVsAreaJet1All[fCentBin]->Fill(mJet1,jet1->Area());
      fh2MassVsNConstJet1All[fCentBin]->Fill(mJet1,jet1->GetNumberOfConstituents());
      
      if(jet1->GetTagStatus()<1 || !jet1->GetTaggedJet())
	continue;

      fh3PtJet1VsMassVsLeadPtTagged[fCentBin]->Fill(ptJet1,mJet1,jet1->MaxTrackPt());
      fpPtVsMassJet1Tagged[fCentBin]->Fill(ptJet1,mJet1);
      fh2MassVsAreaJet1Tagged[fCentBin]->Fill(mJet1,jet1->Area());
      fh2MassVsNConstJet1Tagged[fCentBin]->Fill(mJet1,jet1->GetNumberOfConstituents());

      Double_t fraction = jetCont->GetFractionSharedPt(jet1);
      if(fMinFractionShared>0. && fraction>fMinFractionShared) {
	fh3PtJet1VsMassVsLeadPtTaggedMatch[fCentBin]->Fill(ptJet1,mJet1,jet1->MaxTrackPt());
	fpPtVsMassJet1TaggedMatch[fCentBin]->Fill(ptJet1,mJet1);
	fh2MassVsAreaJet1TaggedMatch[fCentBin]->Fill(mJet1,jet1->Area());
	fh2MassVsNConstJet1TaggedMatch[fCentBin]->Fill(mJet1,jet1->GetNumberOfConstituents());
      }
      
      Double_t Et2 = mJet1*mJet1 + ptJet1*ptJet1;
      Double_t Et = 0.;    Double_t massOverEtR = 0.;
      if(Et2>0.) Et = TMath::Sqrt(Et2);
      if((Et*jetCont->GetJetRadius())>0.) 
	massOverEtR = mJet1/(Et*jetCont->GetJetRadius());
      fh2EtMassOverEtRSq[fCentBin]->Fill(Et,massOverEtR*massOverEtR);
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMass::GetJetMass(AliEmcalJet *jet) {
  //calc subtracted jet mass
  if(fJetMassType==kRaw)
    return jet->M();
  else if(fJetMassType==kDeriv)
    return jet->GetSecondOrderSubtracted();
  
  return 0;
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

