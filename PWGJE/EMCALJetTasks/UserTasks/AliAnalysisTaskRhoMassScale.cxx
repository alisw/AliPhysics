//
// rho mass scale task
// task to estimate scale factor for rho_m
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

#include "AliAnalysisManager.h"
// #include "AliAODMCHeader.h"
// #include "AliMCEvent.h"
// #include "AliGenPythiaEventHeader.h"
// #include "AliAODEvent.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliEmcalParticle.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"

#include "AliAnalysisTaskRhoMassScale.h"

ClassImp(AliAnalysisTaskRhoMassScale)

//________________________________________________________________________
AliAnalysisTaskRhoMassScale::AliAnalysisTaskRhoMassScale() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskRhoMassScale", kTRUE),
  fContainerNeutral(0),
  fContainerCharged(1),
  fRhoMNeutralName(""),
  fRhoMChargedEmcalName(""),
  fRhoMCharged2xEmcalName(""),
  fRhoMNeutral(0),
  fRhoMChargedEmcal(0),
  fRhoMCharged2xEmcal(0),
  fHistScaleEmcalvsCent(0),
  fHistScale2EmcalvsCent(0),
  fHistDeltaScale2EmcalvsCent(0),
  fHistScaleEmcalvsMult(0),
  fHistScale2EmcalvsMult(0),
  fHistDeltaScale2EmcalvsMult(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskRhoMassScale::AliAnalysisTaskRhoMassScale(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerNeutral(0),
  fContainerCharged(1),
  fRhoMNeutralName(""),
  fRhoMChargedEmcalName(""),
  fRhoMCharged2xEmcalName(""),
  fRhoMNeutral(0),
  fRhoMChargedEmcal(0),
  fRhoMCharged2xEmcal(0),
  fHistScaleEmcalvsCent(0),
  fHistScale2EmcalvsCent(0),
  fHistDeltaScale2EmcalvsCent(0),
  fHistScaleEmcalvsMult(0),
  fHistScale2EmcalvsMult(0),
  fHistDeltaScale2EmcalvsMult(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskRhoMassScale::~AliAnalysisTaskRhoMassScale()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskRhoMassScale::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //Create histograms
  TString histName = "";
  TString histTitle = "";
 
  histName = "fHistScaleEmcalvsCent";
  histTitle = TString::Format("%s;Centrality;s_{EMC}",histName.Data());
  fHistScaleEmcalvsCent= new TH2F(histName.Data(),histTitle.Data(), 101, -1, 100, 500, 0, 5);
  fOutput->Add(fHistScaleEmcalvsCent);

  histName = "fHistScale2EmcalvsCent";
  histTitle = TString::Format("%s;Centrality;s_{2 #times EMC}",histName.Data());
  fHistScale2EmcalvsCent = new TH2F(histName.Data(),histTitle.Data(), 101, -1, 100, 500, 0, 5);
  fOutput->Add(fHistScale2EmcalvsCent);

  histName = "fHistDeltaScale2EmcalvsCent";
  histTitle = TString::Format("%s;Centrality;s_{2 #times EMC}-s_{EMC}",histName.Data());
  fHistDeltaScale2EmcalvsCent = new TH2F(histName.Data(),histTitle.Data(), 101, -1, 100, 500, -2.5, 2.5);
  fOutput->Add(fHistDeltaScale2EmcalvsCent);

  histName = "fHistScaleEmcalvsMult";
  histTitle = TString::Format("%s;#it{N}_{track};s_{EMC}",histName.Data());
  fHistScaleEmcalvsMult= new TH2F(histName.Data(),histTitle.Data(), 800, 0, 4000, 500, 0, 5);
  fOutput->Add(fHistScaleEmcalvsMult);

  histName = "fHistScale2EmcalvsMult";
  histTitle = TString::Format("%s;#it{N}_{track};s_{2 #times EMC}",histName.Data());
  fHistScale2EmcalvsMult = new TH2F(histName.Data(),histTitle.Data(), 800, 0, 4000, 500, 0, 5);
  fOutput->Add(fHistScale2EmcalvsMult);

  histName = "fHistDeltaScale2EmcalvsMult";
  histTitle = TString::Format("%s;#it{N}_{track};s_{2 #times EMC}-s_{EMC}",histName.Data());
  fHistDeltaScale2EmcalvsMult = new TH2F(histName.Data(),histTitle.Data(), 800, 0, 4000, 500, -2.5, 2.5);
  fOutput->Add(fHistDeltaScale2EmcalvsMult);
 
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
Bool_t AliAnalysisTaskRhoMassScale::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoMassScale::FillHistograms()
{
  // Fill histograms.

  Double_t rhomNe        = fRhoMNeutral->GetVal();
  Double_t rhomChEmcal   = fRhoMChargedEmcal->GetVal();
  Double_t rhomCh2xEmcal = fRhoMCharged2xEmcal->GetVal();

  Double_t scale = -1.; Double_t scale2 = -1.;
  if(rhomChEmcal>0.)   scale = (rhomNe+rhomChEmcal)/rhomChEmcal;
  if(rhomCh2xEmcal>0.) scale2 = (rhomNe+rhomChEmcal)/rhomCh2xEmcal;

  fHistScaleEmcalvsCent->Fill(fCent,scale);
  fHistScale2EmcalvsCent->Fill(fCent,scale2);
  fHistDeltaScale2EmcalvsCent->Fill(fCent,scale2-scale);

  Int_t mult = -1;
  if(GetParticleContainer(0))
    mult = GetParticleContainer(0)->GetNAcceptedParticles();

  fHistScaleEmcalvsMult->Fill(mult,scale);
  fHistScale2EmcalvsMult->Fill(mult,scale2);
  fHistDeltaScale2EmcalvsMult->Fill(mult,scale2-scale);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoMassScale::RetrieveEventObjects() {
  //
  // retrieve event objects
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  if (!fRhoMNeutralName.IsNull() && !fRhoMNeutral) { // get rho_m from the event
    fRhoMNeutral = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoMNeutralName));
    if (!fRhoMNeutral) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoMNeutralName.Data()));
      fLocalInitialized = kFALSE;
      return kFALSE;
    }
  }

  if (!fRhoMChargedEmcalName.IsNull() && !fRhoMChargedEmcal) { // get rho_m from the event
    fRhoMChargedEmcal = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoMChargedEmcalName));
    if (!fRhoMChargedEmcal) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoMChargedEmcalName.Data()));
      fLocalInitialized = kFALSE;
      return kFALSE;
    }
  }

  if (!fRhoMCharged2xEmcalName.IsNull() && !fRhoMCharged2xEmcal) { // get rho_m from the event
    fRhoMCharged2xEmcal = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoMCharged2xEmcalName));
    if (!fRhoMCharged2xEmcal) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoMCharged2xEmcalName.Data()));
      fLocalInitialized = kFALSE;
      return kFALSE;
    }
  }

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskRhoMassScale::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

