//
// Jet trigger rejection analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVVZERO.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskTriggerRejection.h"

ClassImp(JETriggerRejectionAna::AliAnalysisTaskTriggerRejection)

namespace JETriggerRejectionAna {
  //________________________________________________________________________
  AliAnalysisTaskTriggerRejection::AliAnalysisTaskTriggerRejection() :
    AliAnalysisTaskEmcalJet("AliAnalysisTaskTriggerRejection", kTRUE),
    fContainerFull(0),
    fContainerCharged(1),
    fMaxPatch(0),
    fhnTriggerInfo(0),
    fMainPatchType(kManual),
    fMainTrigCat(kTriggerLevel1Jet),
    fMainTrigSimple(kFALSE)
  {
    // Default constructor.
    SetMakeGeneralHistograms(kTRUE);
  }

  //________________________________________________________________________
  AliAnalysisTaskTriggerRejection::AliAnalysisTaskTriggerRejection(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fContainerFull(0),
    fContainerCharged(1),
    fMaxPatch(0),
    fhnTriggerInfo(0),
    fMainPatchType(kManual),
    fMainTrigCat(kTriggerLevel1Jet),
    fMainTrigSimple(kFALSE)
  {
    // Standard constructor.
    SetMakeGeneralHistograms(kTRUE);
  }

  //________________________________________________________________________
  AliAnalysisTaskTriggerRejection::~AliAnalysisTaskTriggerRejection()
  {
    // Destructor.
  }

  //________________________________________________________________________
  void AliAnalysisTaskTriggerRejection::UserCreateOutputObjects()
  {
    // Create user output.

    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    Int_t fgkNCentBins = 21;
    Float_t kMinCent   = 0.;
    Float_t kMaxCent   = 105.;
    
    Int_t fgkNPtBins = 170;
    Float_t kMinPt   = -50.;
    Float_t kMaxPt   = 120.;

    Int_t fgkNVZEROBins = 100;
    Float_t kMinVZERO   = 0.;
    Float_t kMaxVZERO   = 25000;

    const Int_t fgkNEPatch = 100;
    Float_t kMinEPatch = 0.;
    Float_t kMaxEPatch = 200.;

    const Int_t fgkNADC = 100;
    Float_t kMinADC = 0.;
    Float_t kMaxADC = 1500.;

    const Int_t fgkNEta = 10;
    const Int_t fgkNPhi = 10;

    const Int_t nDim = 8;//cent;V0mult;ptjet1;ptjet2;Epatch;ADCpatch;EtaPatch;PhiPatch
    const Int_t nBins[nDim] = {fgkNCentBins,fgkNVZEROBins,fgkNPtBins,fgkNPtBins,fgkNEPatch,fgkNADC,fgkNEta,fgkNPhi};
    const Double_t xmin0[nDim]  = {kMinCent,kMinVZERO,kMinPt,kMinPt,kMinEPatch,kMinADC,-0.7,1.4};
    const Double_t xmax0[nDim]  = {kMaxCent,kMaxVZERO,kMaxPt,kMaxPt,kMaxEPatch,kMaxADC, 0.7,3.14};
    fhnTriggerInfo = new THnSparseF("fhnTriggerInfo",
                                    "hnTriggerInfo;cent;V0mult;ptjet1;ptjet2;Epatch;ADCpatch;EtaPatch;PhiPatch",nDim,nBins,xmin0,xmax0);
    fOutput->Add(fhnTriggerInfo);

    // =========== Switch on Sumw2 for all histos ===========
    for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
      TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
      if (h1){
        h1->Sumw2();
        continue;
      }
      TH2 *h2 = dynamic_cast<TH2*>(fOutput->At(i));
      if (h2){
        h2->Sumw2();
        continue;
      }
      TH3 *h3 = dynamic_cast<TH3*>(fOutput->At(i));
      if (h3){
        h3->Sumw2();
        continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
      if(hn)hn->Sumw2();
    }

    TH1::AddDirectory(oldStatus);

    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  }

  //________________________________________________________________________
  void AliAnalysisTaskTriggerRejection::ExtractMainPatch() {

    //Find main trigger
    if(!fTriggerPatchInfo)
      return;

    //number of patches in event
    Int_t nPatch = fTriggerPatchInfo->GetEntriesFast();

    //extract main trigger patch
    Double_t emax = -1.;
    for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
      AliEMCALTriggerPatchInfo *patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
      if(patch->GetPatchE()>emax) {
        fMaxPatch = patch;
        emax = patch->GetPatchE();
      }
    }
  }
  //________________________________________________________________________
  Bool_t AliAnalysisTaskTriggerRejection::FillHistograms()
  {
    // Fill histograms.

    if(!GetJetContainer(fContainerFull) || !GetJetContainer(fContainerCharged) || !fMaxPatch)
      return kFALSE;

    //Get leading jets
    AliEmcalJet *leadJet1 = GetJetContainer(fContainerFull)->GetLeadingJet("rho");
    AliEmcalJet *leadJet2 = GetJetContainer(fContainerCharged)->GetLeadingJet("rho");
  
    Double_t ptLeadJet1 = -999;
    Double_t ptLeadJet2 = -999;

    if(leadJet1) ptLeadJet1 = leadJet1->Pt() - GetRhoVal(fContainerFull)*leadJet1->Area();
    if(leadJet2) ptLeadJet2 = leadJet2->Pt() - GetRhoVal(fContainerCharged)*leadJet2->Area();

    Double_t VZEROAmp = (Double_t)(InputEvent()->GetVZEROData()->GetTriggerChargeA() + InputEvent()->GetVZEROData()->GetTriggerChargeC());

    //cent;V0mult;ptjet1;ptjet2;Epatch;ADCpatch;EtaPatch;PhiPatch
    Double_t var[8] = {
      fCent,
      VZEROAmp,
      ptLeadJet1,
      ptLeadJet2,
      fMaxPatch->GetPatchE(),
      (Double_t)fMaxPatch->GetADCAmp(),
      fMaxPatch->GetEtaGeo(),
      fMaxPatch->GetPhiGeo()
    };
    fhnTriggerInfo->Fill(var);

    return kTRUE;
  }

  //________________________________________________________________________
  Bool_t AliAnalysisTaskTriggerRejection::Run()
  {
    // Run analysis code here, if needed. It will be executed before FillHistograms().

    if(fTriggerPatchInfo) { 
      if(fMainPatchType==kManual) ExtractMainPatch();
      else if(fMainPatchType==kEmcalJet) 
        fMaxPatch = GetMainTriggerPatch(fMainTrigCat,fMainTrigSimple);
    }

    return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
  }

  //_______________________________________________________________________
  void AliAnalysisTaskTriggerRejection::Terminate(Option_t *) 
  {
    // Called once at the end of the analysis.
  }
}
