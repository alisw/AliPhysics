//
// Detector response jet mass analysis task.
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

#include "AliAnalysisTaskJetMassResponseDet.h"

ClassImp(AliAnalysisTaskJetMassResponseDet)

//________________________________________________________________________
AliAnalysisTaskJetMassResponseDet::AliAnalysisTaskJetMassResponseDet() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetMassResponseDet", kTRUE),
  fContainerPart(0),
  fContainerDet(0),
  fJetMassType(kRaw),
  fh2PtVsMassJetPartAll(0),
  fh2PtVsMassJetPartMatch(0),
  fh2PtVsMassJetPartTagged(0),
  fh2PtVsMassJetPartTaggedMatch(0),
  fh2PtVsMassJetDetAll(0),
  fh2PtVsMassJetDetTagged(0),
  fhnMassResponse(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskJetMassResponseDet::AliAnalysisTaskJetMassResponseDet(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerPart(0),
  fContainerDet(0),
  fJetMassType(kRaw),
  fh2PtVsMassJetPartAll(0),
  fh2PtVsMassJetPartMatch(0),
  fh2PtVsMassJetPartTagged(0),
  fh2PtVsMassJetPartTaggedMatch(0),
  fh2PtVsMassJetDetAll(0),
  fh2PtVsMassJetDetTagged(0),
  fhnMassResponse(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskJetMassResponseDet::~AliAnalysisTaskJetMassResponseDet()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetMassResponseDet::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt  = 200;
  const Double_t minPt = 0.;
  const Double_t maxPt = 200.;

  const Int_t nBinsM  = 100;
  const Double_t minM = 0.;
  const Double_t maxM = 50.;

  const Int_t nBinsMT  = 50;
  const Double_t minMT = 0.;
  const Double_t maxMT = 50.;

  const Int_t nBinsConst = 26;
  const Double_t minConst = -5.5;
  const Double_t maxConst = 20.5;

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 5;
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsMT};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt, minMT};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt, maxMT};

  //Create histograms
  TString histName = "";
  TString histTitle = "";
  
  histName = "fh2PtVsMassJetPartAll";
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
  fh2PtVsMassJetPartAll = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
  fOutput->Add(fh2PtVsMassJetPartAll);

  histName = "fh2PtVsMassJetPartMatch";
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
  fh2PtVsMassJetPartMatch = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
  fOutput->Add(fh2PtVsMassJetPartMatch);

  histName = "fh2PtVsMassJetPartTagged";
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
  fh2PtVsMassJetPartTagged = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
  fOutput->Add(fh2PtVsMassJetPartTagged);

  histName = "fh2PtVsMassJetPartTaggedMatch";
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
  fh2PtVsMassJetPartTaggedMatch = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
  fOutput->Add(fh2PtVsMassJetPartTaggedMatch);

  histName = "fh2PtVsMassJetDetAll";
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
  fh2PtVsMassJetDetAll = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
  fOutput->Add(fh2PtVsMassJetDetAll);

  histName = "fh2PtVsMassJetDetTagged";
  histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1}",histName.Data());
  fh2PtVsMassJetDetTagged = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
  fOutput->Add(fh2PtVsMassJetDetTagged);

  histName = "fhnMassResponse";
  histTitle = Form("%s;#it{M}_{det};#it{M}_{part};#it{p}_{T,det};#it{p}_{T,part};#it{M}_{det}^{tagged}",histName.Data());
  fhnMassResponse = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0,nBins0,xmin0,xmax0);
  fOutput->Add(fhnMassResponse);


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
Bool_t AliAnalysisTaskJetMassResponseDet::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetMassResponseDet::FillHistograms()
{
  // Fill histograms.

  AliJetContainer *cPart = GetJetContainer(fContainerPart);
  AliJetContainer *cDet = GetJetContainer(fContainerDet);
  AliEmcalJet* jPart = NULL;
  AliEmcalJet* jDet = NULL;

  //loop on particle level jets
  if(cPart) {
    cPart->ResetCurrentID();
    while((jPart = cPart->GetNextAcceptJet())) {
      fh2PtVsMassJetPartAll->Fill(jPart->Pt(),jPart->M());
      jDet = jPart->ClosestJet();
      if(jDet) fh2PtVsMassJetPartMatch->Fill(jPart->Pt(),jPart->M());
      if(jPart->GetTagStatus()<1 || !jPart->GetTaggedJet())
	continue;
      fh2PtVsMassJetPartTagged->Fill(jPart->Pt(),jPart->M());
      if(jDet) fh2PtVsMassJetPartTaggedMatch->Fill(jPart->Pt(),jPart->M());
    }
  }
  
  //loop on detector level jets
  if(cDet) {
    cDet->ResetCurrentID();
    while((jDet = cDet->GetNextAcceptJet())) {
      Double_t mjet = GetJetMass(jDet);     
      fh2PtVsMassJetDetAll->Fill(jDet->Pt(),mjet);
       if(jDet->GetTagStatus()>=1 && jDet->GetTaggedJet())
	 fh2PtVsMassJetDetTagged->Fill(jDet->Pt(),mjet);
       
       //fill detector response
       jPart = jDet->ClosestJet();
       if(jPart) {
	 AliEmcalJet *jDetT = jDet->GetTaggedJet();
	 Double_t mdetT = 0.;
	 if(jDetT) mdetT = jDetT->M();
	 Double_t var[5] = {GetJetMass(jDet),jPart->M(),jDet->Pt(),jPart->Pt(),mdetT};
	 fhnMassResponse->Fill(var);

	 if(jPart->Pt()>40. && jPart->Pt()<50.) {
	   if(jDet->Pt()>50.) Printf("feed-out high");
	   else if(jDet->Pt()<40.) Printf("feed-out low");
	   else Printf("correct");
	   Printf("pT Part: %f Det: %f",jPart->Pt(),jDet->Pt());
	   Printf("mass Part: %f Det: %f",jPart->M(),jDet->M());
	   Int_t nConstPart = jPart->GetNumberOfConstituents();
	   Int_t nConstDet = jDet->GetNumberOfConstituents();
	   Int_t diff = nConstPart-nConstDet;
	   Printf("nConst Part: %d  Det: %d  diff: %d",nConstPart,nConstDet,diff);
	 }
       }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskJetMassResponseDet::GetJetMass(AliEmcalJet *jet) {
  //calc subtracted jet mass
  if(fJetMassType==kRaw)
    return jet->M();
  else if(fJetMassType==kDeriv)
    return jet->GetSecondOrderSubtracted();
  
  return 0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetMassResponseDet::RetrieveEventObjects() {
  //
  // retrieve event objects
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskJetMassResponseDet::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

