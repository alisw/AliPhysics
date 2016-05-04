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
#include <TProfile.h>

#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalJetMass.h"

ClassImp(AliAnalysisTaskEmcalJetMass)

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::AliAnalysisTaskEmcalJetMass() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetMass", kTRUE),
  fContainerBase(0),
  fContainerUnsub(1),
  fMinFractionShared(0),
  fUseUnsubJet(0),
  fJetMassType(kRaw),
  fUseSumw2(0),
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
  fh3PtJet1VsMassVsCentAllSel(0),
  fh3PtJet1VsMassVsCentTagged(0),
  fh3PtJet1VsMassVsCentTaggedMatch(0),
  fh3PtJet1VsRatVsLeadPtAllSel(0),
  fh3PtJet1VsRatVsLeadPtTagged(0),
  fh3PtJet1VsRatVsLeadPtTaggedMatch(0),
  fpPtVsRatJet1All(0),
  fpPtVsRatJet1Tagged(0),
  fpPtVsRatJet1TaggedMatch(0),
  fh2RatVsAreaJet1All(0),
  fh2RatVsAreaJet1Tagged(0),
  fh2RatVsAreaJet1TaggedMatch(0),
  fh2RatVsNConstJet1All(0),
  fh2RatVsNConstJet1Tagged(0),
  fh2RatVsNConstJet1TaggedMatch(0),
  fh3JetPtVsMassVsEPRelAllSel(0),
  fh3JetPtVsMassVsEPRelTagged(0),
  fh3JetPtVsMassVsEPRelTaggedMatch(0),
  fh3RhoVsLeadJetPtVsCent(0),
  fh3RhoMVsLeadJetPtVsCent(0)
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

  fh3PtJet1VsRatVsLeadPtAllSel         = new TH3F*[fNcentBins];
  fh3PtJet1VsRatVsLeadPtTagged         = new TH3F*[fNcentBins];
  fh3PtJet1VsRatVsLeadPtTaggedMatch    = new TH3F*[fNcentBins];
  fpPtVsRatJet1All                     = new TProfile*[fNcentBins];
  fpPtVsRatJet1Tagged                  = new TProfile*[fNcentBins];
  fpPtVsRatJet1TaggedMatch             = new TProfile*[fNcentBins];
  fh2RatVsAreaJet1All                  = new TH2F*[fNcentBins];
  fh2RatVsAreaJet1Tagged               = new TH2F*[fNcentBins];
  fh2RatVsAreaJet1TaggedMatch          = new TH2F*[fNcentBins];
  fh2RatVsNConstJet1All                = new TH2F*[fNcentBins];
  fh2RatVsNConstJet1Tagged             = new TH2F*[fNcentBins];
  fh2RatVsNConstJet1TaggedMatch        = new TH2F*[fNcentBins];

  fh3JetPtVsMassVsEPRelAllSel          = new TH3F*[fNcentBins];
  fh3JetPtVsMassVsEPRelTagged          = new TH3F*[fNcentBins];
  fh3JetPtVsMassVsEPRelTaggedMatch     = new TH3F*[fNcentBins];
  
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

    fh3PtJet1VsRatVsLeadPtAllSel[i]         = 0;
    fh3PtJet1VsRatVsLeadPtTagged[i]         = 0;
    fh3PtJet1VsRatVsLeadPtTaggedMatch[i]    = 0;
    fpPtVsRatJet1All[i]                     = 0;
    fpPtVsRatJet1Tagged[i]                  = 0;
    fpPtVsRatJet1TaggedMatch[i]             = 0;
    fh2RatVsAreaJet1All[i]                  = 0;
    fh2RatVsAreaJet1Tagged[i]               = 0;
    fh2RatVsAreaJet1TaggedMatch[i]          = 0;
    fh2RatVsNConstJet1All[i]                = 0;
    fh2RatVsNConstJet1Tagged[i]             = 0;
    fh2RatVsNConstJet1TaggedMatch[i]        = 0;

    fh3JetPtVsMassVsEPRelAllSel[i]          = 0;
    fh3JetPtVsMassVsEPRelTagged[i]          = 0;
    fh3JetPtVsMassVsEPRelTaggedMatch[i]     = 0;
    
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMass::AliAnalysisTaskEmcalJetMass(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fContainerUnsub(1),
  fMinFractionShared(0),
  fUseUnsubJet(0),
  fJetMassType(kRaw),
  fUseSumw2(0),
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
  fh3PtJet1VsMassVsCentAllSel(0),
  fh3PtJet1VsMassVsCentTagged(0),
  fh3PtJet1VsMassVsCentTaggedMatch(0),
  fh3PtJet1VsRatVsLeadPtAllSel(0),
  fh3PtJet1VsRatVsLeadPtTagged(0),
  fh3PtJet1VsRatVsLeadPtTaggedMatch(0),
  fpPtVsRatJet1All(0),
  fpPtVsRatJet1Tagged(0),
  fpPtVsRatJet1TaggedMatch(0),
  fh2RatVsAreaJet1All(0),
  fh2RatVsAreaJet1Tagged(0),
  fh2RatVsAreaJet1TaggedMatch(0),
  fh2RatVsNConstJet1All(0),
  fh2RatVsNConstJet1Tagged(0),
  fh2RatVsNConstJet1TaggedMatch(0),
  fh3JetPtVsMassVsEPRelAllSel(0),
  fh3JetPtVsMassVsEPRelTagged(0),
  fh3JetPtVsMassVsEPRelTaggedMatch(0),
  fh3RhoVsLeadJetPtVsCent(0),
  fh3RhoMVsLeadJetPtVsCent(0)
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

  fh3PtJet1VsRatVsLeadPtAllSel         = new TH3F*[fNcentBins];
  fh3PtJet1VsRatVsLeadPtTagged         = new TH3F*[fNcentBins];
  fh3PtJet1VsRatVsLeadPtTaggedMatch    = new TH3F*[fNcentBins];
  fpPtVsRatJet1All                     = new TProfile*[fNcentBins];
  fpPtVsRatJet1Tagged                  = new TProfile*[fNcentBins];
  fpPtVsRatJet1TaggedMatch             = new TProfile*[fNcentBins];
  fh2RatVsAreaJet1All                  = new TH2F*[fNcentBins];
  fh2RatVsAreaJet1Tagged               = new TH2F*[fNcentBins];
  fh2RatVsAreaJet1TaggedMatch          = new TH2F*[fNcentBins];
  fh2RatVsNConstJet1All                = new TH2F*[fNcentBins];
  fh2RatVsNConstJet1Tagged             = new TH2F*[fNcentBins];
  fh2RatVsNConstJet1TaggedMatch        = new TH2F*[fNcentBins];

  fh3JetPtVsMassVsEPRelAllSel          = new TH3F*[fNcentBins];
  fh3JetPtVsMassVsEPRelTagged          = new TH3F*[fNcentBins];
  fh3JetPtVsMassVsEPRelTaggedMatch     = new TH3F*[fNcentBins];
  
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

    fh3PtJet1VsRatVsLeadPtAllSel[i]         = 0;
    fh3PtJet1VsRatVsLeadPtTagged[i]         = 0;
    fh3PtJet1VsRatVsLeadPtTaggedMatch[i]    = 0;
    fpPtVsRatJet1All[i]                     = 0;
    fpPtVsRatJet1Tagged[i]                  = 0;
    fpPtVsRatJet1TaggedMatch[i]             = 0;
    fh2RatVsAreaJet1All[i]                  = 0;
    fh2RatVsAreaJet1Tagged[i]               = 0;
    fh2RatVsAreaJet1TaggedMatch[i]          = 0;
    fh2RatVsNConstJet1All[i]                = 0;
    fh2RatVsNConstJet1Tagged[i]             = 0;
    fh2RatVsNConstJet1TaggedMatch[i]        = 0;

    fh3JetPtVsMassVsEPRelAllSel[i]          = 0;
    fh3JetPtVsMassVsEPRelTagged[i]          = 0;
    fh3JetPtVsMassVsEPRelTaggedMatch[i]     = 0;
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

  const Int_t nBinsR  = 100;
  const Double_t minR = -0.2;
  const Double_t maxR = 0.8;

  const Int_t nBinsArea = 50;
  const Double_t minArea = 0.;
  const Double_t maxArea = 1.;

  const Int_t nBinsNConst = 200;
  const Double_t minNConst = 0.;
  const Double_t maxNConst = 200.;
  
  Int_t nBinsRho = 200;
  Double_t minRho = 0.;
  Double_t maxRho = 20.;
  Int_t nBinsRhom = 200;
  Double_t minRhom = 0.;
  Double_t maxRhom = 1.;

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
    fh2MassVsNConstJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1Tagged[i]);

    histName = TString::Format("fh2MassVsNConstJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1};#it{N}_{constituents}",histName.Data());
    fh2MassVsNConstJet1TaggedMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2MassVsNConstJet1TaggedMatch[i]);

    //
    histName = TString::Format("fh3PtJet1VsRatVsLeadPtAllSel_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
    fh3PtJet1VsRatVsLeadPtAllSel[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,20,0.,20.);
    fOutput->Add(fh3PtJet1VsRatVsLeadPtAllSel[i]);

    histName = TString::Format("fh3PtJet1VsRatVsLeadPtTagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
    fh3PtJet1VsRatVsLeadPtTagged[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,20,0.,20.);
    fOutput->Add(fh3PtJet1VsRatVsLeadPtTagged[i]);

    histName = TString::Format("fh3PtJet1VsRatVsLeadPtTaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}/#it{p}_{T,jet1};#it{p}_{T,lead trk}",histName.Data());
    fh3PtJet1VsRatVsLeadPtTaggedMatch[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsR,minR,maxR,20,0.,20.);
    fOutput->Add(fh3PtJet1VsRatVsLeadPtTaggedMatch[i]);

    histName = TString::Format("fpPtVsRatJet1All_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}/#it{p}_{T,jet1}",histName.Data());
    fpPtVsRatJet1All[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsRatJet1All[i]);

    histName = TString::Format("fpPtVsRatJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}/#it{p}_{T,jet1}",histName.Data());
    fpPtVsRatJet1Tagged[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsRatJet1Tagged[i]);

    histName = TString::Format("fpPtVsRatJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};Avg #it{M}_{jet1}/#it{p}_{T,jet1}",histName.Data());
    fpPtVsRatJet1TaggedMatch[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsRatJet1TaggedMatch[i]);

    histName = TString::Format("fh2RatVsAreaJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1}/#it{p}_{T,jet1};#it{A}",histName.Data());
    fh2RatVsAreaJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2RatVsAreaJet1All[i]);

    histName = TString::Format("fh2RatVsAreaJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1}/#it{p}_{T,jet1};#it{A}",histName.Data());
    fh2RatVsAreaJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2RatVsAreaJet1Tagged[i]);

    histName = TString::Format("fh2RatVsAreaJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1}/#it{p}_{T,jet1};#it{A}",histName.Data());
    fh2RatVsAreaJet1TaggedMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsArea,minArea,maxArea);
    fOutput->Add(fh2RatVsAreaJet1TaggedMatch[i]);

    histName = TString::Format("fh2RatVsNConstJet1All_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1}/#it{p}_{T,jet1};#it{N}_{constituents}",histName.Data());
    fh2RatVsNConstJet1All[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2RatVsNConstJet1All[i]);

    histName = TString::Format("fh2RatVsNConstJet1Tagged_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1}/#it{p}_{T,jet1};#it{N}_{constituents}",histName.Data());
    fh2RatVsNConstJet1Tagged[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2RatVsNConstJet1Tagged[i]);

    histName = TString::Format("fh2RatVsNConstJet1TaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{M}_{jet1}/#it{p}_{T,jet1};#it{N}_{constituents}",histName.Data());
    fh2RatVsNConstJet1TaggedMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsNConst,minNConst,maxNConst);
    fOutput->Add(fh2RatVsNConstJet1TaggedMatch[i]);

    histName = TString::Format("fh3JetPtVsMassVsEPRelAllSel_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#varphi_{jet}-#Psi",histName.Data());
    fh3JetPtVsMassVsEPRelAllSel[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,100,0.,TMath::Pi());
    fOutput->Add(fh3JetPtVsMassVsEPRelAllSel[i]);

    histName = TString::Format("fh3JetPtVsMassVsEPRelTagged_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#varphi_{jet}-#Psi",histName.Data());
    fh3JetPtVsMassVsEPRelTagged[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,100,0.,TMath::Pi());
    fOutput->Add(fh3JetPtVsMassVsEPRelTagged[i]);

    histName = TString::Format("fh3JetPtVsMassVsEPRelTaggedMatch_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#varphi_{jet}-#Psi",histName.Data());
    fh3JetPtVsMassVsEPRelTaggedMatch[i] =  new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,100,0.,TMath::Pi());
    fOutput->Add(fh3JetPtVsMassVsEPRelTaggedMatch[i]);
  }

  histName = TString::Format("fh3PtJet1VsMassVsCentAllSel");
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};centrality",histName.Data());
  fh3PtJet1VsMassVsCentAllSel = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,101,-1.,100.);
  fOutput->Add(fh3PtJet1VsMassVsCentAllSel);

  histName = TString::Format("fh3PtJet1VsMassVsCentTagged");
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};centrality",histName.Data());
  fh3PtJet1VsMassVsCentTagged = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,101,-1.,100.);
  fOutput->Add(fh3PtJet1VsMassVsCentTagged);

  histName = TString::Format("fh3PtJet1VsMassVsCentTaggedMatch");
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};centrality",histName.Data());
  fh3PtJet1VsMassVsCentTaggedMatch = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,101,-1.,100.);
  fOutput->Add(fh3PtJet1VsMassVsCentTaggedMatch);
  
  histName = "fh2RhoVsLeadJetPtVsCent";
  histTitle = Form("Rho Vs LeadJetPt Vs Cent; #rho;#it{p}_{T,lead-jet} (GeV/#it{c});centrality");
  fh3RhoVsLeadJetPtVsCent = new TH3F(histName.Data(), histTitle.Data(), nBinsRho, minRho, maxRho, nBinsPt, minPt, maxPt, 101, -1., 100.);
  fOutput->Add(fh3RhoVsLeadJetPtVsCent);
  
  histName = "fh2RhoMVsLeadJetPtVsCent";
  histTitle = Form("RhoM Vs LeadJetPt Vs Cent; #rho_{m};#it{p}_{T,lead-jet} (GeV/#it{c});centrality");
  fh3RhoMVsLeadJetPtVsCent = new TH3F(histName.Data(), histTitle.Data(), nBinsRhom, minRhom, maxRhom, nBinsPt, minPt, maxPt, 101, -1., 100.);
  fOutput->Add(fh3RhoMVsLeadJetPtVsCent);

  if(fUseSumw2) {
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
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMass::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  if(!GetNParticles()) {
     return kFALSE;
     
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMass::FillHistograms()
{
  // Fill histograms.
  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(jetCont) {
     AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(jetCont->GetRhoMassName()));
     Double_t rhomassval = 0;
     if (!rhomParam) {
      AliError(Form("%s: Could not retrieve rho_m %s!", GetName(), jetCont->GetRhoMassName().Data()));
      
    } else rhomassval = rhomParam->GetVal();

    AliEmcalJet *lj = jetCont->GetLeadingJet(); //leading pt without rho subtraction
    if(lj){
       fh3RhoVsLeadJetPtVsCent->Fill(jetCont->GetRhoVal(),lj->Pt(),fCent);
       fh3RhoMVsLeadJetPtVsCent->Fill(rhomassval,lj->Pt(),fCent);
    }                
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {

      Double_t ptJet1 = jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
      Double_t maxTrackPt = jet1->MaxTrackPt();
      Double_t mJet1 = GetJetMass(jet1);
      Double_t eJet1 = jet1->E();
      Double_t rat = -1.;
      if(eJet1<0. || eJet1>0.) rat = mJet1/eJet1;
      Double_t ep = jet1->Phi() - fEPV0;
      while (ep < 0) ep += TMath::Pi();
      while (ep >= TMath::Pi()) ep -= TMath::Pi();

      //Printf("useUnsub: %d jet: %d pT: %f M: %f E: %f",fUseUnsubJet,jetCont->GetCurrentID(),ptJet1,mJet1,eJet1);

      Double_t fraction = -1.;
      AliEmcalJet *jetUS = NULL;
      if(fUseUnsubJet) {
        AliJetContainer *jetContUS = GetJetContainer(fContainerUnsub);
        Int_t ifound = 0;
        Int_t ilab = -1;
        for(Int_t i = 0; i<jetContUS->GetNJets(); i++) {
          jetUS = jetContUS->GetJet(i);
          if(jetUS->GetLabel()==jet1->GetLabel()) {
            ifound++;
            if(ifound==1) ilab = i;
          }
        }
        if(ifound>1) AliDebug(2,Form("Found %d partners",ifound));
        if(ilab>-1) {
          jetUS = jetContUS->GetJet(ilab);
          if(!jetUS) continue;
          fraction = jetContUS->GetFractionSharedPt(jetUS);
          maxTrackPt = jetUS->MaxTrackPt();
        }
        //Printf("subtracted jet ijet: %d lab: %d eta: %f phi: %f",jetCont->GetCurrentID(),jet1->GetLabel(),jet1->Eta(),jet1->Phi());
        //Printf("unsubtracted jet ijet: %d lab: %d eta: %f phi: %f",ilab,jetUS->GetLabel(),jetUS->Eta(),jetUS->Phi());
        
        fraction = jetContUS->GetFractionSharedPt(jetUS);
        maxTrackPt = jetUS->MaxTrackPt();
      } else
	fraction = jetCont->GetFractionSharedPt(jet1);

      fh3PtJet1VsMassVsLeadPtAllSel[fCentBin]->Fill(ptJet1,mJet1,maxTrackPt);
      fpPtVsMassJet1All[fCentBin]->Fill(ptJet1,mJet1);
      fh2MassVsAreaJet1All[fCentBin]->Fill(mJet1,jet1->Area());
      fh2MassVsNConstJet1All[fCentBin]->Fill(mJet1,jet1->GetNumberOfConstituents());
      fh3PtJet1VsMassVsCentAllSel->Fill(ptJet1,mJet1,fCent);
      fh3JetPtVsMassVsEPRelAllSel[fCentBin]->Fill(ptJet1,mJet1,ep);

      fh3PtJet1VsRatVsLeadPtAllSel[fCentBin]->Fill(ptJet1,rat,maxTrackPt);
      fpPtVsRatJet1All[fCentBin]->Fill(ptJet1,rat);
      fh2RatVsAreaJet1All[fCentBin]->Fill(rat,jet1->Area());
      fh2RatVsNConstJet1All[fCentBin]->Fill(rat,jet1->GetNumberOfConstituents());
      
      if(jet1->GetTagStatus()<1 || !jet1->GetTaggedJet())
	continue;

      fh3PtJet1VsMassVsLeadPtTagged[fCentBin]->Fill(ptJet1,mJet1,maxTrackPt);
      fpPtVsMassJet1Tagged[fCentBin]->Fill(ptJet1,mJet1);
      fh2MassVsAreaJet1Tagged[fCentBin]->Fill(mJet1,jet1->Area());
      fh2MassVsNConstJet1Tagged[fCentBin]->Fill(mJet1,jet1->GetNumberOfConstituents());
      fh3PtJet1VsMassVsCentTagged->Fill(ptJet1,mJet1,fCent);
      fh3JetPtVsMassVsEPRelTagged[fCentBin]->Fill(ptJet1,mJet1,ep);

      fh3PtJet1VsRatVsLeadPtTagged[fCentBin]->Fill(ptJet1,rat,maxTrackPt);
      fpPtVsRatJet1Tagged[fCentBin]->Fill(ptJet1,rat);
      fh2RatVsAreaJet1Tagged[fCentBin]->Fill(rat,jet1->Area());
      fh2RatVsNConstJet1Tagged[fCentBin]->Fill(rat,jet1->GetNumberOfConstituents());

      //matching
      if(fMinFractionShared>0. && fraction>fMinFractionShared) {
	AliEmcalJet *jetM = NULL;
	if(fUseUnsubJet) {
	  if(jetUS) jetM = jetUS->ClosestJet();
	}
	else jetM = jet1->ClosestJet();
	if(jetM) maxTrackPt = jetM->MaxTrackPt();
	else     maxTrackPt = -1.;
	fh3PtJet1VsMassVsLeadPtTaggedMatch[fCentBin]->Fill(ptJet1,mJet1,maxTrackPt);
	fpPtVsMassJet1TaggedMatch[fCentBin]->Fill(ptJet1,mJet1);
	fh2MassVsAreaJet1TaggedMatch[fCentBin]->Fill(mJet1,jet1->Area());
	fh2MassVsNConstJet1TaggedMatch[fCentBin]->Fill(mJet1,jet1->GetNumberOfConstituents());
        fh3PtJet1VsMassVsCentTaggedMatch->Fill(ptJet1,mJet1,fCent);
	fh3JetPtVsMassVsEPRelTaggedMatch[fCentBin]->Fill(ptJet1,mJet1,ep);

	fh3PtJet1VsRatVsLeadPtTaggedMatch[fCentBin]->Fill(ptJet1,rat,maxTrackPt);
	fpPtVsRatJet1TaggedMatch[fCentBin]->Fill(ptJet1,rat);
	fh2RatVsAreaJet1TaggedMatch[fCentBin]->Fill(rat,jet1->Area());
	fh2RatVsNConstJet1TaggedMatch[fCentBin]->Fill(rat,jet1->GetNumberOfConstituents());
      }
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
    return jet->GetShapeProperties()->GetSecondOrderSubtracted();
  
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

