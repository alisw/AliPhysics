#include <stdio.h>
#include <stdlib.h>

#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TList.h"
#include "TChain.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliAnaVZEROEPFlatenning.h"
#include "AliCentrality.h"
#include "AliEventplane.h"

// VZERO includes
#include "AliVVZERO.h"

ClassImp(AliAnaVZEROEPFlatenning)

AliAnaVZEROEPFlatenning::AliAnaVZEROEPFlatenning() 
  : AliAnalysisTaskSE("AliAnaVZEROEPFlatenning"), fEvent(0), fOutputList(0),
  fMBTrigName("CPBI"),
  fUsePhysSel(kFALSE),
  fPsiARawCentr(NULL),
  fPsiAFlatCentr(NULL),
  fPsiCRawCentr(NULL),
  fPsiCFlatCentr(NULL),
  fPsiACRawCentr(NULL),
  fPsiACFlatCentr(NULL)
{
  // Default constructor
  // Init pointers
  for(Int_t i = 0; i < 11; ++i) fX2[i] = fY2[i] = fX2Y2[i] = fCos8Psi[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fX2In[i] = fY2In[i] = fX2Y2In[i] = fCos8PsiIn[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fPsiRingRawCentr[i] = fPsiRingFlatCentr[i] = fPsiRingFlatFinalCentr[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fC2[i] = fS2[i] = fC4[i] = fS4[i] = NULL;
  for(Int_t i = 0; i < 11; ++i) fX2Corr[i] = fY2Corr[i] = fX2Y2Corr[i] = NULL;
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnaVZEROEPFlatenning::AliAnaVZEROEPFlatenning(const char *name) 
  : AliAnalysisTaskSE(name), fEvent(0), fOutputList(0),
  fMBTrigName("CPBI"),
  fUsePhysSel(kFALSE),
  fPsiARawCentr(NULL),
  fPsiAFlatCentr(NULL),
  fPsiCRawCentr(NULL),
  fPsiCFlatCentr(NULL),
  fPsiACRawCentr(NULL),
  fPsiACFlatCentr(NULL)
{
  // Constructor
  // Init pointers
  for(Int_t i = 0; i < 11; ++i) fX2[i] = fY2[i] = fX2Y2[i] = fCos8Psi[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fX2In[i] = fY2In[i] = fX2Y2In[i] = fCos8PsiIn[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fPsiRingRawCentr[i] = fPsiRingFlatCentr[i] = fPsiRingFlatFinalCentr[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fC2[i] = fS2[i] = fC4[i] = fS4[i] = NULL;
  for(Int_t i = 0; i < 11; ++i) fX2Corr[i] = fY2Corr[i] = fX2Y2Corr[i] = NULL;
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnaVZEROEPFlatenning::SetInput(const char *filename)
{
  // Read recentering
  // parameters from an input file

  AliInfo(Form("Reading input histograms from %s",filename));
  TFile *f = TFile::Open(filename);
  TList *list = (TList*)f->Get("coutput");

  for(Int_t i = 0; i < 8; ++i) {
    fX2In[i] = (TProfile*)list->FindObject(Form("fX2_%d",i))->Clone(Form("fX2In_%d",i));
    fX2In[i]->SetDirectory(0);
    fY2In[i] = (TProfile*)list->FindObject(Form("fY2_%d",i))->Clone(Form("fY2In_%d",i));
    fY2In[i]->SetDirectory(0);
    fX2Y2In[i] = (TProfile*)list->FindObject(Form("fX2Y2_%d",i))->Clone(Form("fX2Y2In_%d",i));
    fX2Y2In[i]->SetDirectory(0);
    fCos8PsiIn[i] = (TProfile*)list->FindObject(Form("fCos8Psi_%d",i))->Clone(Form("fCos8PsiIn_%d",i));
    fCos8PsiIn[i]->SetDirectory(0);
  }
  f->Close();
}

//________________________________________________________________________
void AliAnaVZEROEPFlatenning::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  for(Int_t i = 0; i < 11; ++i) {
    fX2[i] = new TProfile(Form("fX2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fX2[i]);
    fY2[i] = new TProfile(Form("fY2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fY2[i]);
    fX2Y2[i] = new TProfile(Form("fX2Y2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fX2Y2[i]);
    fCos8Psi[i] = new TProfile(Form("fCos8Psi_%d",i),"",21,0,105,"s");
    fOutputList->Add(fCos8Psi[i]);
  }
  for(Int_t i = 0; i < 11; ++i) {
    fX2Corr[i] = new TProfile(Form("fX2Corr_%d",i),"",21,0,105,"s");
    fOutputList->Add(fX2Corr[i]);
    fY2Corr[i] = new TProfile(Form("fY2Corr_%d",i),"",21,0,105,"s");
    fOutputList->Add(fY2Corr[i]);
    fX2Y2Corr[i] = new TProfile(Form("fX2Y2Corr_%d",i),"",21,0,105,"s");
    fOutputList->Add(fX2Y2Corr[i]);
  }

  for(Int_t i = 0; i < 8; ++i) {
    fPsiRingRawCentr[i] = new TH2F(Form("fPsiRingRawCentr_%d",i),"",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
    fOutputList->Add(fPsiRingRawCentr[i]);
    fPsiRingFlatCentr[i] = new TH2F(Form("fPsiRingFlatCentr_%d",i),"",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
    fOutputList->Add(fPsiRingFlatCentr[i]);
    fPsiRingFlatFinalCentr[i] = new TH2F(Form("fPsiRingFlatFinalCentr_%d",i),"",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
    fOutputList->Add(fPsiRingFlatFinalCentr[i]);
  }
  fPsiARawCentr = new TH2F("fPsiARawCentr","",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
  fOutputList->Add(fPsiARawCentr);
  fPsiAFlatCentr = new TH2F("fPsiAFlatCentr","",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
  fOutputList->Add(fPsiAFlatCentr);
  fPsiCRawCentr = new TH2F("fPsiCRawCentr","",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
  fOutputList->Add(fPsiCRawCentr);
  fPsiCFlatCentr = new TH2F("fPsiCFlatCentr","",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
  fOutputList->Add(fPsiCFlatCentr);
  fPsiACRawCentr = new TH2F("fPsiACRawCentr","",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
  fOutputList->Add(fPsiACRawCentr);
  fPsiACFlatCentr = new TH2F("fPsiACFlatCentr","",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
  fOutputList->Add(fPsiACFlatCentr);

  for(Int_t i = 0; i < 8; ++i) {
    fC2[i] = new TProfile(Form("fC2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fC2[i]);
    fS2[i] = new TProfile(Form("fS2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fS2[i]);
    fC4[i] = new TProfile(Form("fC4_%d",i),"",21,0,105,"s");
    fOutputList->Add(fC4[i]);
    fS4[i] = new TProfile(Form("fS4_%d",i),"",21,0,105,"s");
    fOutputList->Add(fS4[i]);
  }

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnaVZEROEPFlatenning::Init()
{
  // Nothing here
  // ....
}

//________________________________________________________________________
void AliAnaVZEROEPFlatenning::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event


  fEvent = InputEvent();
  if (!fEvent) {
    printf("ERROR: fEvent not available\n");
    return;
  }

  AliVVZERO* esdV0 = fEvent->GetVZEROData();
  if (!esdV0) {
    Printf("ERROR: esd/aod V0  not available");
    return;
  }

  // Phys sel
  Bool_t goodEvent = kTRUE;
  Bool_t isSelected;
  if (fUsePhysSel)
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMB | AliVEvent::kSemiCentral));
  else
    isSelected = ((esdV0->GetV0ADecision()==1) && (esdV0->GetV0CDecision()==1));

  if (!isSelected) goodEvent = kFALSE;

  AliCentrality *centrality = fEvent->GetCentrality();
  //  Float_t percentile = centrality->GetCentralityPercentile("V0M");
  Float_t spdPercentile = centrality->GetCentralityPercentile("CL1");

  TString inputDataType = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (inputDataType == "ESD") {
    AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
    // Trigger
    TString trigStr(esdEvent->GetFiredTriggerClasses());
    if (!trigStr.Contains("-B-")) return;
    if (!trigStr.Contains(fMBTrigName.Data())) return;

    const AliESDVertex *primaryVtx = esdEvent->GetPrimaryVertexSPD();
    if (!primaryVtx) goodEvent = kFALSE;
    if (!primaryVtx->GetStatus()) goodEvent = kFALSE;
    Double_t tPrimaryVtxPosition[3];
    primaryVtx->GetXYZ(tPrimaryVtxPosition);
    if (TMath::Abs(tPrimaryVtxPosition[2]) > 10.0) goodEvent = kFALSE;
  }
  else if (inputDataType == "AOD") {
    AliAODEvent* aodEvent = static_cast<AliAODEvent*>(fEvent);
    // Trigger
    TString trigStr(aodEvent->GetHeader()->GetFiredTriggerClasses());
    if (!trigStr.Contains("-B-")) return;
    if (!trigStr.Contains(fMBTrigName.Data())) return;

    const AliAODVertex *primaryVtx = aodEvent->GetPrimaryVertexSPD();
    if (!primaryVtx) goodEvent = kFALSE;
    if (primaryVtx->GetNContributors()==0) goodEvent = kFALSE;
    Double_t tPrimaryVtxPosition[3];
    primaryVtx->GetXYZ(tPrimaryVtxPosition);
    if (TMath::Abs(tPrimaryVtxPosition[2]) > 10.0) goodEvent = kFALSE;
  }
  else {
    AliFatal("Trying to get the vertex from neither ESD nor AOD event!");
    return;
  }

  if (goodEvent) {
    Double_t totMult = 0, totMultA = 0, totMultC = 0;
    Double_t c2A = 0, c2C = 0, s2A = 0, s2C = 0;
    Double_t qxA = 0, qyA = 0;
    Double_t qxC = 0, qyC = 0;
    for(Int_t iring = 0; iring < 8; ++iring) {
      Double_t ringMult = 0;
      Double_t c2 = 0, s2 = 0;
      for(Int_t iCh = iring*8; iCh < (iring+1)*8; ++iCh) {
	Double_t phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
	c2 += fEvent->GetVZEROEqMultiplicity(iCh)*TMath::Cos(2.*phi);
	s2 += fEvent->GetVZEROEqMultiplicity(iCh)*TMath::Sin(2.*phi);
	if (fEvent->GetVZEROEqMultiplicity(iCh) > 1e-6) {
	  fC2[iring]->Fill(spdPercentile,TMath::Cos(2.*phi),fEvent->GetVZEROEqMultiplicity(iCh));
	  fS2[iring]->Fill(spdPercentile,TMath::Sin(2.*phi),fEvent->GetVZEROEqMultiplicity(iCh));
	  fC4[iring]->Fill(spdPercentile,TMath::Cos(4.*phi),fEvent->GetVZEROEqMultiplicity(iCh));
	  fS4[iring]->Fill(spdPercentile,TMath::Sin(4.*phi),fEvent->GetVZEROEqMultiplicity(iCh));
	}
	ringMult += fEvent->GetVZEROEqMultiplicity(iCh);
      }
      if (ringMult < 1e-6) continue;
      totMult += ringMult;
      if (iring >= 4) {
	totMultA += ringMult;
	c2A += c2;
	s2A += s2;
      }
      else {
	totMultC += ringMult;
	c2C += c2;
	s2C += s2;
      }
      fX2[iring]->Fill(spdPercentile,c2);
      fY2[iring]->Fill(spdPercentile,s2);
      fX2Y2[iring]->Fill(spdPercentile,c2*s2);

      Double_t qxOut = 0, qyOut = 0;
      Double_t psiRingRaw = fEvent->GetEventplane()->CalculateVZEROEventPlane(fEvent,iring,iring,2,qxOut,qyOut);
      fPsiRingRawCentr[iring]->Fill(spdPercentile,psiRingRaw);
      fCos8Psi[iring]->Fill(spdPercentile, 2./4.*TMath::Cos(2.*4.*psiRingRaw));
      Double_t qxTierce = 0, qyTierce = 0;
      Double_t psiRingFlat = CalculateVZEROEventPlane(fEvent,iring,spdPercentile,qxTierce,qyTierce);
      fPsiRingFlatCentr[iring]->Fill(spdPercentile,psiRingFlat);
      Int_t ibin = fCos8PsiIn[iring]->FindBin(spdPercentile);
      Double_t psiRingFlatFinal = psiRingFlat + (fCos8PsiIn[iring]->GetBinContent(ibin)*TMath::Sin(2.*4.*psiRingFlat))/2.;
      if (psiRingFlatFinal > TMath::Pi()/2) psiRingFlatFinal -= TMath::Pi();
      if (psiRingFlatFinal <-TMath::Pi()/2) psiRingFlatFinal += TMath::Pi();
      fPsiRingFlatFinalCentr[iring]->Fill(spdPercentile,psiRingFlatFinal);
      fX2Corr[iring]->Fill(spdPercentile,qxTierce);
      fY2Corr[iring]->Fill(spdPercentile,qyTierce);
      fX2Y2Corr[iring]->Fill(spdPercentile,qxTierce*qyTierce);
      if (iring >= 4) {
	qxA += qxTierce;
	qyA += qyTierce;
      }
      else {
 	qxC += qxTierce;
	qyC += qyTierce;
     }
    }

    if (totMultA > 1e-6) {
      fX2[8]->Fill(spdPercentile,c2A);
      fY2[8]->Fill(spdPercentile,s2A);
      fX2Y2[8]->Fill(spdPercentile,c2A*s2A);
      fX2Corr[8]->Fill(spdPercentile,qxA);
      fY2Corr[8]->Fill(spdPercentile,qyA);
      fX2Y2Corr[8]->Fill(spdPercentile,qxA*qyA);
      Double_t psiA = TMath::ATan2(qyA,qxA)/2.;
      fPsiAFlatCentr->Fill(spdPercentile,psiA);
      Double_t psiAOrg = fEvent->GetEventplane()->GetEventplane("V0A",fEvent,2);
      fPsiARawCentr->Fill(spdPercentile,psiAOrg);
      fCos8Psi[8]->Fill(spdPercentile, 2./4.*TMath::Cos(2.*4.*psiAOrg));
    }
    if (totMultC > 1e-6) {
      fX2[9]->Fill(spdPercentile,c2C);
      fY2[9]->Fill(spdPercentile,s2C);
      fX2Y2[9]->Fill(spdPercentile,c2C*s2C);
      fX2Corr[9]->Fill(spdPercentile,qxC);
      fY2Corr[9]->Fill(spdPercentile,qyC);
      fX2Y2Corr[9]->Fill(spdPercentile,qxC*qyC);
      Double_t psiC = TMath::ATan2(qyC,qxC)/2.;
      fPsiCFlatCentr->Fill(spdPercentile,psiC);
      Double_t psiCOrg = fEvent->GetEventplane()->GetEventplane("V0C",fEvent,2);
      fPsiCRawCentr->Fill(spdPercentile,psiCOrg);
      fCos8Psi[9]->Fill(spdPercentile, 2./4.*TMath::Cos(2.*4.*psiCOrg));
   }
    if (totMult > 1e-6) {
      fX2[10]->Fill(spdPercentile,c2A+c2C);
      fY2[10]->Fill(spdPercentile,s2A+s2C);
      fX2Y2[10]->Fill(spdPercentile,(c2A+c2C)*(s2A+s2C));
      fX2Corr[10]->Fill(spdPercentile,qxA+qxC);
      fY2Corr[10]->Fill(spdPercentile,qyA+qyC);
      fX2Y2Corr[10]->Fill(spdPercentile,(qxA+qxC)*(qyA+qyC));
      Double_t psiAC = TMath::ATan2(qyA+qyC,qxA+qxC)/2.;
      fPsiACFlatCentr->Fill(spdPercentile,psiAC);
      Double_t psiACOrg = fEvent->GetEventplane()->GetEventplane("V0",fEvent,2);
      fPsiACRawCentr->Fill(spdPercentile,psiACOrg);
      fCos8Psi[10]->Fill(spdPercentile, 2./4.*TMath::Cos(2.*4.*psiACOrg));
    }
  }

  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnaVZEROEPFlatenning::Terminate(Option_t *) 
{
  // Check the output list and store output config file
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

Double_t AliAnaVZEROEPFlatenning::CalculateVZEROEventPlane(const AliVEvent *  event, Int_t ring, Float_t centrality, Double_t &qxTierce, Double_t &qyTierce) const
{
  // Calculate the VZERO event plane
  // taking into account the recentering/twist and rescaling
  // of the cumulants
  qxTierce = qyTierce = 0.;
  if(!event) {
    AliError("No Event received");
    return -1000.;
  }
  AliVVZERO *vzeroData = event->GetVZEROData();
  if(!vzeroData) {
    AliError("Enable to get VZERO Data");
    return -1000.;
  }

  Int_t ibin = fX2In[ring]->FindBin(centrality);
  Double_t meanX2 = fX2In[ring]->GetBinContent(ibin);
  Double_t meanY2 = fY2In[ring]->GetBinContent(ibin);
  Double_t sigmaX2 = fX2In[ring]->GetBinError(ibin);
  Double_t sigmaY2 = fY2In[ring]->GetBinError(ibin);
  Double_t rho = (fX2Y2In[ring]->GetBinContent(ibin)-meanX2*meanY2)/sigmaX2/sigmaY2;

  Double_t b = rho*sigmaX2*sigmaY2*
    TMath::Sqrt(2.*(sigmaX2*sigmaX2+sigmaY2*sigmaY2-2.*sigmaX2*sigmaY2*TMath::Sqrt(1.-rho*rho))/
		((sigmaX2*sigmaX2-sigmaY2*sigmaY2)*(sigmaX2*sigmaX2-sigmaY2*sigmaY2)+
		 4.*sigmaX2*sigmaX2*sigmaY2*sigmaY2*rho*rho));
  Double_t aPlus = TMath::Sqrt(2.*sigmaX2*sigmaX2-b*b);
  Double_t aMinus= TMath::Sqrt(2.*sigmaY2*sigmaY2-b*b);

  Double_t lambdaPlus = b/aPlus;
  Double_t lambdaMinus = b/aMinus;

  Double_t qx=0., qy=0.;
  for(Int_t iCh = ring*8; iCh < (ring+1)*8; ++iCh) {
    Double_t phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    Double_t mult = event->GetVZEROEqMultiplicity(iCh);
    qx += mult*TMath::Cos(2.*phi);
    qy += mult*TMath::Sin(2.*phi);
  }
  // Recentering
  Double_t qxPrime = qx - meanX2;
  Double_t qyPrime = qy - meanY2;
  // Twist
  Double_t trans[2][2];
  trans[0][0] = 1./(1.-lambdaPlus*lambdaMinus);
  trans[0][1] = -lambdaMinus/(1.-lambdaPlus*lambdaMinus);
  trans[1][0] = -lambdaPlus/(1.-lambdaPlus*lambdaMinus);
  trans[1][1] = 1./(1.-lambdaPlus*lambdaMinus);
  Double_t qxSeconde = trans[0][0]*qxPrime + trans[0][1]*qyPrime;
  Double_t qySeconde = trans[1][0]*qxPrime + trans[1][1]*qyPrime;
  // Rescaling
  qxTierce = qxSeconde/aPlus;
  qyTierce = qySeconde/aMinus;

  return (TMath::ATan2(qyTierce,qxTierce)/2.);
}
