#include <stdio.h>
#include <stdlib.h>

#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TList.h"
#include "TChain.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnaVZEROEPFlatenning.h"
#include "AliCentrality.h"
#include "AliEventplane.h"

// VZERO includes
#include "AliESDVZERO.h"

ClassImp(AliAnaVZEROEPFlatenning)

AliAnaVZEROEPFlatenning::AliAnaVZEROEPFlatenning() 
  : AliAnalysisTaskSE("AliAnaVZEROEPFlatenning"), fESD(0), fOutputList(0),
  fMBTrigName("CPBI"),
  fUsePhysSel(kFALSE)
{
  // Default constructor
  // Init pointers
  for(Int_t i = 0; i < 8; ++i) fX2[i] = fY2[i] = fX2Y2[i] = fCos8Psi[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fX2In[i] = fY2In[i] = fX2Y2In[i] = fCos8PsiIn[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fPsiRingRawCentr[i] = fPsiRingFlatCentr[i] = fPsiRingFlatFinalCentr[i] = NULL;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnaVZEROEPFlatenning::AliAnaVZEROEPFlatenning(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0),
  fMBTrigName("CPBI"),
  fUsePhysSel(kFALSE)
{
  // Constructor
  // Init pointers
  for(Int_t i = 0; i < 8; ++i) fX2[i] = fY2[i] = fX2Y2[i] = fCos8Psi[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fX2In[i] = fY2In[i] = fX2Y2In[i] = fCos8PsiIn[i] = NULL;
  for(Int_t i = 0; i < 8; ++i) fPsiRingRawCentr[i] = fPsiRingFlatCentr[i] = fPsiRingFlatFinalCentr[i] = NULL;

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

  for(Int_t i = 0; i < 8; ++i) {
    fX2[i] = new TProfile(Form("fX2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fX2[i]);
    fY2[i] = new TProfile(Form("fY2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fY2[i]);
    fX2Y2[i] = new TProfile(Form("fX2Y2_%d",i),"",21,0,105,"s");
    fOutputList->Add(fX2Y2[i]);
    fCos8Psi[i] = new TProfile(Form("fCos8Psi_%d",i),"",21,0,105,"s");
    fOutputList->Add(fCos8Psi[i]);
  }

  for(Int_t i = 0; i < 8; ++i) {
    fPsiRingRawCentr[i] = new TH2F(Form("fPsiRingRawCentr_%d",i),"",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
    fOutputList->Add(fPsiRingRawCentr[i]);
    fPsiRingFlatCentr[i] = new TH2F(Form("fPsiRingFlatCentr_%d",i),"",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
    fOutputList->Add(fPsiRingFlatCentr[i]);
    fPsiRingFlatFinalCentr[i] = new TH2F(Form("fPsiRingFlatFinalCentr_%d",i),"",105,0,105,100,-TMath::Pi()/2,TMath::Pi()/2);
    fOutputList->Add(fPsiRingFlatFinalCentr[i]);
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


  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  if (!esdV0) {
    Printf("ERROR: esd V0  not available");
    return;
  }

  // Trigger
  TString trigStr(fESD->GetFiredTriggerClasses());
  if (!trigStr.Contains("-B-")) return;
  if (!trigStr.Contains(fMBTrigName.Data())) return;

  // Phys sel
  Bool_t goodEvent = kTRUE;
  Bool_t isSelected;
  if (fUsePhysSel)
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  else
    isSelected = ((esdV0->GetV0ADecision()==1) && (esdV0->GetV0CDecision()==1));

  if (!isSelected) goodEvent = kFALSE;

  AliCentrality *centrality = fESD->GetCentrality();
  //  Float_t percentile = centrality->GetCentralityPercentile("V0M");
  Float_t spdPercentile = centrality->GetCentralityPercentile("CL1");

  const AliESDVertex *primaryVtx = fESD->GetPrimaryVertexSPD();
  if (!primaryVtx) goodEvent = kFALSE;
  if (!primaryVtx->GetStatus()) goodEvent = kFALSE;
  Double_t tPrimaryVtxPosition[3];
  primaryVtx->GetXYZ(tPrimaryVtxPosition);
  if (TMath::Abs(tPrimaryVtxPosition[2]) > 10.0) goodEvent = kFALSE;

  if (goodEvent) {
    for(Int_t iring = 0; iring < 8; ++iring) {
      Double_t c2 = 0;
      Double_t s2 = 0;
      Double_t c4 = 0;
      Double_t s4 = 0;
      Double_t totMult = 0;
      for(Int_t iCh = iring*8; iCh < (iring+1)*8; ++iCh) {
	Double_t phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
	c2 += fESD->GetVZEROEqMultiplicity(iCh)*TMath::Cos(2.*phi);
	s2 += fESD->GetVZEROEqMultiplicity(iCh)*TMath::Sin(2.*phi);
	c4 += fESD->GetVZEROEqMultiplicity(iCh)*TMath::Cos(4.*phi);
	s4 += fESD->GetVZEROEqMultiplicity(iCh);
	totMult += fESD->GetVZEROEqMultiplicity(iCh);
      }
      if (totMult < 1e-6) continue;

      fX2[iring]->Fill(spdPercentile,c2);
      fY2[iring]->Fill(spdPercentile,s2);
      fX2Y2[iring]->Fill(spdPercentile,c2*s2);

      Double_t psiRingRaw = fESD->GetEventplane()->CalculateVZEROEventPlane(fESD,iring,iring,2);
      fPsiRingRawCentr[iring]->Fill(spdPercentile,psiRingRaw);
      Double_t psiRingFlat2 = CalculateVZEROEventPlane(fESD,iring,spdPercentile);
      fPsiRingFlatCentr[iring]->Fill(spdPercentile,psiRingFlat2);
      fCos8Psi[iring]->Fill(spdPercentile, 2./4.*TMath::Cos(2.*4.*psiRingFlat2));
      Int_t ibin = fCos8PsiIn[iring]->FindBin(spdPercentile);
      Double_t psiRingFlatFinal = psiRingFlat2 + (fCos8PsiIn[iring]->GetBinContent(ibin)*TMath::Sin(2.*4.*psiRingFlat2))/2.;
      if (psiRingFlatFinal > TMath::Pi()/2) psiRingFlatFinal -= TMath::Pi();
      if (psiRingFlatFinal <-TMath::Pi()/2) psiRingFlatFinal += TMath::Pi();
      fPsiRingFlatFinalCentr[iring]->Fill(spdPercentile,psiRingFlatFinal);
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

Double_t AliAnaVZEROEPFlatenning::CalculateVZEROEventPlane(const AliVEvent *  event, Int_t ring, Float_t centrality) const
{
  // Calculate the VZERO event plane
  // taking into account the recentering/twist and rescaling
  // of the cumulants
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

  Double_t trans[2][2];
  trans[0][0] = 1./(1.-lambdaPlus*lambdaMinus);
  trans[0][1] = -lambdaMinus/(1.-lambdaPlus*lambdaMinus);
  trans[1][0] = -lambdaPlus/(1.-lambdaPlus*lambdaMinus);
  trans[1][1] = 1./(1.-lambdaPlus*lambdaMinus);

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
  Double_t qxSeconde = trans[0][0]*qxPrime + trans[0][1]*qyPrime;
  Double_t qySeconde = trans[1][0]*qxPrime + trans[1][1]*qyPrime;
  // Rescaling
  Double_t qxTierce = qxSeconde/aPlus;
  Double_t qyTierce = qySeconde/aMinus;

  return (TMath::ATan2(qyTierce,qxTierce)/2.);
}
