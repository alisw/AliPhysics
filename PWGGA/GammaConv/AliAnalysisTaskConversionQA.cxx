/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																								*
 * Authors: Svein Lindal, Daniel Lohner											*
 * Version 1.0																				*
 *																								*
 * Permission to use, copy, modify and distribute this software and its	 *
 * documentation strictly for non-commercial purposes is hereby granted	 *
 * without fee, provided that the above copyright notice appears in all	 *
 * copies and that both the copyright notice and this permission notice	 *
 * appear in the supporting documentation. The authors make no claims	 *
 * about the suitability of this software for any purpose. It is			*
 * provided "as is" without express or implied warranty.						*
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskConversionQA.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskConversionQA)

//________________________________________________________________________
AliAnalysisTaskConversionQA::AliAnalysisTaskConversionQA(const char *name) : AliAnalysisTaskSE(name),
   fConversionGammas(NULL),
   fConversionCuts(NULL),
   fStreamQA(NULL),
   fIsHeavyIon(kFALSE),
   fOutputList(NULL)
{
   // Default constructor

   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskConversionQA::~AliAnalysisTaskConversionQA()
{
   // default deconstructor
   if(fStreamQA){
      delete fStreamQA;
      fStreamQA = 0x0;
   }
}
//________________________________________________________________________
void AliAnalysisTaskConversionQA::UserCreateOutputObjects()
{
   // Create User Output Objects

   if(fOutputList != NULL){
      delete fOutputList;
      fOutputList = NULL;
   }
   if(fOutputList == NULL){
      fOutputList = new TList();
      fOutputList->SetOwner(kTRUE);
   }

   // V0 Reader Cuts
   TString cutnumber = fConversionCuts->GetCutNumber();

   fStreamQA = new TTreeSRedirector(Form("GammaConvV1_QA_%s.root",cutnumber.Data()));
   PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskConversionQA::UserExec(Option_t *){

   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");

   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
   if(eventQuality != 0){// Event Not Accepted
      return;
   }
   AliESDEvent* event = (AliESDEvent*) InputEvent();
   if(fIsHeavyIon && !fConversionCuts->IsCentralitySelected(event)) return;

   fConversionGammas=fV0Reader->GetReconstructedGammas();

   ProcessQA();
   PostData(1, fOutputList);
}


///________________________________________________________________________
void AliAnalysisTaskConversionQA::ProcessQA(){

   // Fill Histograms for QA and MC
   AliESDEvent* event = (AliESDEvent*) InputEvent();
   AliPIDResponse* pidResonse = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();
   for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
      AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
      if(!fConversionCuts->PhotonIsSelected(gamma,event)) continue;
      Float_t gammaPt = gamma->GetPhotonPt();
      Float_t gammaPhi = gamma->GetPhotonPhi();
      Float_t gammaTheta = gamma->Theta();
      Float_t gammaChi2NDF = gamma->GetChi2perNDF();
      Float_t gammaQt = gamma->GetArmenterosQt();
      Float_t gammaAlpha = gamma->GetArmenterosAlpha();
      Float_t gammaPsiPair = gamma->GetPsiPair();
      TVectorF daughterProp(14);
      AliVTrack * negTrack = fConversionCuts->GetTrack(event, gamma->GetTrackLabelNegative());
      AliVTrack * posTrack = fConversionCuts->GetTrack(event, gamma->GetTrackLabelPositive());
      if(!negTrack||!posTrack)return;

      daughterProp(0) = posTrack->Pt();
      daughterProp(7) = negTrack->Pt();
      daughterProp(1) = posTrack->Theta();
      daughterProp(8) = negTrack->Theta();
      daughterProp(2) = posTrack->GetTPCsignal();
      daughterProp(9) = negTrack->GetTPCsignal();
      daughterProp(3) = pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
      daughterProp(10) = pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
      if((posTrack->GetStatus() & AliESDtrack::kTOFpid) && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
         daughterProp(4) = posTrack->GetTOFsignal();
         daughterProp(5) = pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron);
      } else {
         daughterProp(4) = 20000;
         daughterProp(5) = -20;
      }
      if((negTrack->GetStatus() & AliESDtrack::kTOFpid) && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
         daughterProp(11) = negTrack->GetTOFsignal();
         daughterProp(12) = pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron);
      } else {
         daughterProp(11) = 20000;
         daughterProp(12) = -20;
      }
      daughterProp(6) = (Float_t)posTrack->GetNcls(1)/(Float_t)posTrack->GetTPCNclsF();
      daughterProp(13) = (Float_t)negTrack->GetNcls(1)/(Float_t)negTrack->GetTPCNclsF();

      if (fStreamQA){
         (*fStreamQA)<<"PhotonQA"
                     << "pt=" << gammaPt
                     << "phi=" << gammaPhi
                     << "theta=" << gammaTheta
                     << "chi2ndf=" << gammaChi2NDF
                     << "qt="<< gammaQt
                     << "alpha=" << gammaAlpha
                     << "psipair=" << gammaPsiPair
                     << "daugtherProp.=" << &daughterProp
                     << "\n";
      }
   }
}


//________________________________________________________________________
void AliAnalysisTaskConversionQA::Terminate(Option_t *)
{
   if (fStreamQA){
      fStreamQA->GetFile()->Write();
   }
}
