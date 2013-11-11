/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *				       					  *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliConversionMesonCuts.h"

#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "AliStack.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"
#include "AliAODMCParticle.h"

class iostream;

using namespace std;

ClassImp(AliConversionMesonCuts)


const char* AliConversionMesonCuts::fgkCutNames[AliConversionMesonCuts::kNCuts] = {
   "MesonKind", //0
   "BackgroundScheme", //1
   "NumberOfBGEvents", //2
   "DegreesForRotationMethod", //3
   "RapidityMesonCut", //4
   "RCut", //5
   "AlphaMesonCut", //6
   "Chi2MesonCut", //7
   "SharedElectronCuts", //8
   "RejectToCloseV0s", //9
   "UseMCPSmearing", //10
   "DcaGammaGamma", //11
   "DcaRPrimVtx", //12
   "DcaZPrimVtx" //13
};


//________________________________________________________________________
AliConversionMesonCuts::AliConversionMesonCuts(const char *name,const char *title) :
   AliAnalysisCuts(name,title),
   fHistograms(NULL),
   fMesonKind(0),
   fMaxR(200),
   fChi2CutMeson(1000),
   fAlphaMinCutMeson(0),
   fAlphaCutMeson(1),
   fRapidityCutMeson(1),
   fUseRotationMethodInBG(kFALSE),
   fDoBG(kTRUE),
   fdoBGProbability(kFALSE),
   fUseTrackMultiplicityForBG(kFALSE),
   fnDegreeRotationPMForBG(0),
   fNumberOfBGEvents(0),
   fOpeningAngle(0.005),
   fDoToCloseV0sCut(kFALSE),
   fminV0Dist(200.),
   fDoSharedElecCut(kFALSE),
   fUseMCPSmearing(kFALSE),
   fPBremSmearing(0),
   fPSigSmearing(0),
   fPSigSmearingCte(0),
   fBrem(NULL),
   fRandom(0),
   fFAlphaCut(0),
   fAlphaPtDepCut(kFALSE),
   fElectronLabelArraySize(500),
   fElectronLabelArray(NULL),
   fDCAGammaGammaCut(1000),
   fDCAZMesonPrimVtxCut(1000),
   fDCARMesonPrimVtxCut(1000),
   fBackgroundHandler(0),
   fCutString(NULL),
   hMesonCuts(NULL),
   hMesonBGCuts(NULL),
   hDCAGammaGammaMesonBefore(NULL),
   hDCAZMesonPrimVtxBefore(NULL),
   hDCARMesonPrimVtxBefore(NULL),
   hDCAGammaGammaMesonAfter(NULL),
   hDCAZMesonPrimVtxAfter(NULL),
   hDCARMesonPrimVtxAfter(NULL)

{
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
   fCutString=new TObjString((GetCutNumber()).Data());
   fElectronLabelArray = new Int_t[fElectronLabelArraySize];
   if (fBrem == NULL){
      fBrem = new TF1("fBrem","pow(-log(x),[0]/log(2.0)-1.0)/TMath::Gamma([0]/log(2.0))",0.00001,0.999999999);
      // tests done with 1.0e-14
      fBrem->SetParameter(0,fPBremSmearing);
      fBrem->SetNpx(100000);
   }

}

//________________________________________________________________________
AliConversionMesonCuts::AliConversionMesonCuts(const AliConversionMesonCuts &ref) :
   AliAnalysisCuts(ref),
   fHistograms(NULL),
   fMesonKind(ref.fMesonKind),
   fMaxR(ref.fMaxR),
   fChi2CutMeson(ref.fChi2CutMeson),
   fAlphaMinCutMeson(ref.fAlphaMinCutMeson),
   fAlphaCutMeson(ref.fAlphaCutMeson),
   fRapidityCutMeson(ref.fRapidityCutMeson),
   fUseRotationMethodInBG(ref.fUseRotationMethodInBG),
   fDoBG(ref.fDoBG),
   fdoBGProbability(ref.fdoBGProbability),
   fUseTrackMultiplicityForBG(ref.fUseTrackMultiplicityForBG),
   fnDegreeRotationPMForBG(ref.fnDegreeRotationPMForBG),
   fNumberOfBGEvents(ref. fNumberOfBGEvents),
   fOpeningAngle(ref.fOpeningAngle),
   fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
   fminV0Dist(ref.fminV0Dist),
   fDoSharedElecCut(ref.fDoSharedElecCut),
   fUseMCPSmearing(ref.fUseMCPSmearing),
   fPBremSmearing(ref.fPBremSmearing),
   fPSigSmearing(ref.fPSigSmearing),
   fPSigSmearingCte(ref.fPSigSmearingCte),
   fBrem(NULL),
   fRandom(ref.fRandom),
   fFAlphaCut(NULL),
   fAlphaPtDepCut(ref.fAlphaPtDepCut),
   fElectronLabelArraySize(ref.fElectronLabelArraySize),
   fElectronLabelArray(NULL),
   fDCAGammaGammaCut(ref.fDCAGammaGammaCut),
   fDCAZMesonPrimVtxCut(ref.fDCAZMesonPrimVtxCut),
   fDCARMesonPrimVtxCut(ref.fDCARMesonPrimVtxCut),
   fBackgroundHandler(ref.fBackgroundHandler),
   fCutString(NULL),
   hMesonCuts(NULL),
   hMesonBGCuts(NULL),
   hDCAGammaGammaMesonBefore(NULL),
   hDCAZMesonPrimVtxBefore(NULL),
   hDCARMesonPrimVtxBefore(NULL),
   hDCAGammaGammaMesonAfter(NULL),
   hDCAZMesonPrimVtxAfter(NULL),
   hDCARMesonPrimVtxAfter(NULL)
{
   // Copy Constructor
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
   fCutString=new TObjString((GetCutNumber()).Data());
   fElectronLabelArray = new Int_t[fElectronLabelArraySize];
   if (fBrem == NULL)fBrem = (TF1*)ref.fBrem->Clone("fBrem");
   // Histograms are not copied, if you need them, call InitCutHistograms
}


//________________________________________________________________________
AliConversionMesonCuts::~AliConversionMesonCuts() {
   // Destructor
   //Deleting fHistograms leads to seg fault it it's added to output collection of a task
   // if(fHistograms)
   // 	delete fHistograms;
   // fHistograms = NULL;
   if(fCutString != NULL){
      delete fCutString;
      fCutString = NULL;
   }
   if(fElectronLabelArray){
      delete fElectronLabelArray;
      fElectronLabelArray = NULL;
   }

}

//________________________________________________________________________
void AliConversionMesonCuts::InitCutHistograms(TString name, Bool_t additionalHists){

   // Initialize Cut Histograms for QA (only initialized and filled if function is called)
   TH1::AddDirectory(kFALSE);

   if(fHistograms != NULL){
      delete fHistograms;
      fHistograms=NULL;
   }
   if(fHistograms==NULL){
      fHistograms=new TList();
      fHistograms->SetOwner(kTRUE);
      if(name=="")fHistograms->SetName(Form("ConvMesonCuts_%s",GetCutNumber().Data()));
      else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
   }

   // Meson Cuts
   hMesonCuts=new TH1F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts",13,-0.5,12.5);
   hMesonCuts->GetXaxis()->SetBinLabel(1,"in");
   hMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
   hMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
   hMesonCuts->GetXaxis()->SetBinLabel(4,"opening angle");
   hMesonCuts->GetXaxis()->SetBinLabel(5,"alpha max");
   hMesonCuts->GetXaxis()->SetBinLabel(6,"alpha min");
   hMesonCuts->GetXaxis()->SetBinLabel(7,"dca gamma gamma");
   hMesonCuts->GetXaxis()->SetBinLabel(8,"dca R prim Vtx");
   hMesonCuts->GetXaxis()->SetBinLabel(9,"dca Z prim Vtx");
   hMesonCuts->GetXaxis()->SetBinLabel(10,"out");
   fHistograms->Add(hMesonCuts);

   hMesonBGCuts=new TH1F(Form("MesonBGCuts %s",GetCutNumber().Data()),"MesonBGCuts",13,-0.5,12.5);
   hMesonBGCuts->GetXaxis()->SetBinLabel(1,"in");
   hMesonBGCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
   hMesonBGCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
   hMesonBGCuts->GetXaxis()->SetBinLabel(4,"opening angle");
   hMesonBGCuts->GetXaxis()->SetBinLabel(5,"alpha max");
   hMesonBGCuts->GetXaxis()->SetBinLabel(6,"alpha min");
   hMesonBGCuts->GetXaxis()->SetBinLabel(7,"dca gamma gamma");
   hMesonBGCuts->GetXaxis()->SetBinLabel(8,"dca R prim Vtx");
   hMesonBGCuts->GetXaxis()->SetBinLabel(9,"dca Z prim Vtx");
   hMesonBGCuts->GetXaxis()->SetBinLabel(10,"out");

   fHistograms->Add(hMesonBGCuts);

   if (additionalHists){
      hDCAGammaGammaMesonBefore=new TH1F(Form("DCAGammaGammaMeson Before %s",GetCutNumber().Data()),"DCAGammaGammaMeson Before",200,0,10);
      fHistograms->Add(hDCAGammaGammaMesonBefore);

      hDCARMesonPrimVtxBefore=new TH1F(Form("DCARMesonPrimVtx Before %s",GetCutNumber().Data()),"DCARMesonPrimVtx Before",200,0,10);
      fHistograms->Add(hDCARMesonPrimVtxBefore);

      hDCAZMesonPrimVtxBefore=new TH1F(Form("DCAZMesonPrimVtx Before %s",GetCutNumber().Data()),"DCAZMesonPrimVtx Before",401,-10,10);
      fHistograms->Add(hDCAZMesonPrimVtxBefore);

   }

   hDCAGammaGammaMesonAfter=new TH1F(Form("DCAGammaGammaMeson After %s",GetCutNumber().Data()),"DCAGammaGammaMeson After",200,0,10);
   fHistograms->Add(hDCAGammaGammaMesonAfter);

   hDCAZMesonPrimVtxAfter=new TH2F(Form("InvMassDCAZMesonPrimVtx After %s",GetCutNumber().Data()),"InvMassDCAZMesonPrimVtx After",800,0,0.8,401,-10,10);
   fHistograms->Add(hDCAZMesonPrimVtxAfter);

   hDCARMesonPrimVtxAfter=new TH1F(Form("DCARMesonPrimVtx After %s",GetCutNumber().Data()),"DCARMesonPrimVtx After",200,0,10);
   fHistograms->Add(hDCARMesonPrimVtxAfter);

   TH1::AddDirectory(kTRUE);
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack, Double_t fRapidityShift){
   // Returns true for all pions within acceptance cuts for decay into 2 photons
   // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

   if(!fMCStack)return kFALSE;

   if(fMCMother->GetPdgCode()==111 || fMCMother->GetPdgCode()==221){
      if(fMCMother->R()>fMaxR)	return kFALSE; // cuts on distance from collision point

      Double_t rapidity = 10.;
      if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
         rapidity=8.-fRapidityShift;
      } else{
         rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
      }

      // Rapidity Cut
      if(abs(rapidity)>fRapidityCutMeson)return kFALSE;

      // Select only -> 2y decay channel
      if(fMCMother->GetNDaughters()!=2)return kFALSE;

      for(Int_t i=0;i<2;i++){
         TParticle *MDaughter=fMCStack->Particle(fMCMother->GetDaughter(i));
         // Is Daughter a Photon?
         if(MDaughter->GetPdgCode()!=22)return kFALSE;
         // Is Photon in Acceptance?
         //   if(bMCDaughtersInAcceptance){
         //	if(!PhotonIsSelectedMC(MDaughter,fMCStack)){return kFALSE;}
         //   }
      }
      return kTRUE;
   }
   return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMC(AliAODMCParticle *MCMother,TClonesArray *AODMCArray, Double_t fRapidityShift){
   // Returns true for all pions within acceptance cuts for decay into 2 photons
   // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

   if(!AODMCArray)return kFALSE;

   if(MCMother->GetPdgCode()==111 || MCMother->GetPdgCode()==221){
      Double_t rMeson = sqrt( (MCMother->Xv()*MCMother->Xv()) + (MCMother->Yv()*MCMother->Yv()) ) ;
      if(rMeson>fMaxR)	return kFALSE; // cuts on distance from collision point

      Double_t rapidity = 10.;
      if(MCMother->E() - MCMother->Pz() == 0 || MCMother->E() + MCMother->Pz() == 0){
         rapidity=8.-fRapidityShift;
      } else{
         rapidity = 0.5*(TMath::Log((MCMother->E()+MCMother->Pz()) / (MCMother->E()-MCMother->Pz())))-fRapidityShift;
      }

      // Rapidity Cut
      if(abs(rapidity)>fRapidityCutMeson)return kFALSE;

      // Select only -> 2y decay channel
      if(MCMother->GetNDaughters()!=2)return kFALSE;

      for(Int_t i=0;i<2;i++){
         AliAODMCParticle *MDaughter=static_cast<AliAODMCParticle*>(AODMCArray->At(MCMother->GetDaughter(i)));
         // Is Daughter a Photon?
         if(MDaughter->GetPdgCode()!=22)return kFALSE;
         // Is Photon in Acceptance?
         //   if(bMCDaughtersInAcceptance){
         //	if(!PhotonIsSelectedMC(MDaughter,fMCStack)){return kFALSE;}
         //   }
      }
      return kTRUE;
   }
   return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCDalitz(TParticle *fMCMother,AliStack *fMCStack, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma, Double_t fRapidityShift){

   // Returns true for all pions within acceptance cuts for decay into 2 photons
   // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

   if( !fMCStack )return kFALSE;

   if(	fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 ) return kFALSE;

   if(  fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

   Double_t rapidity = 10.;

   if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
      rapidity=8.-fRapidityShift;
   }
   else{
      rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
   }

   // Rapidity Cut
   if( abs(rapidity) > fRapidityCutMeson )return kFALSE;

   // Select only -> Dalitz decay channel
   if( fMCMother->GetNDaughters() != 3 )return kFALSE;

   TParticle *positron = 0x0;
   TParticle *electron = 0x0;
   TParticle    *gamma = 0x0;

   for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){

      TParticle* temp = (TParticle*)fMCStack->Particle( index );

      switch( temp->GetPdgCode() ) {
      case ::kPositron:
         positron      =  temp;
         labelpositron = index;
         break;
      case ::kElectron:
         electron      =  temp;
         labelelectron = index;
         break;
      case ::kGamma:
         gamma         =  temp;
         labelgamma    = index;
         break;
      }
   }

   if( positron && electron && gamma) return kTRUE;
   return kFALSE;


}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCChiC(TParticle *fMCMother,AliStack *fMCStack,Int_t & labelelectronChiC, Int_t & labelpositronChiC, Int_t & labelgammaChiC, Double_t fRapidityShift){
   // Returns true for all ChiC within acceptance cuts for decay into JPsi + gamma -> e+ + e- + gamma
   // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

   if(!fMCStack)return kFALSE;
   // if(fMCMother->GetPdgCode()==20443 ){
   // 	 return kFALSE;
   // }
   if(fMCMother->GetPdgCode()==10441 || fMCMother->GetPdgCode()==10443 || fMCMother->GetPdgCode()==445 ){
      if(fMCMother->R()>fMaxR)	return kFALSE; // cuts on distance from collision point

      Double_t rapidity = 10.;
      if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
         rapidity=8.-fRapidityShift;
      }
      else{
         rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
      }

      // Rapidity Cut
      if(abs(rapidity)>fRapidityCutMeson)return kFALSE;

      // Select only -> ChiC radiative (JPsi+gamma) decay channel
      if(fMCMother->GetNDaughters()!=2)return kFALSE;

      TParticle *jpsi 	= 0x0;
      TParticle *gamma 	= 0x0;
      TParticle *positron = 0x0;
      TParticle *electron = 0x0;

      Int_t labeljpsiChiC = -1;

      for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){

         TParticle* temp = (TParticle*)fMCStack->Particle( index );

         switch( temp->GetPdgCode() ) {
         case 443:
            jpsi =  temp;
            labeljpsiChiC = index;
            break;
         case 22:
            gamma    =  temp;
            labelgammaChiC = index;
            break;
         }
      }

      if ( !jpsi || ! gamma) return kFALSE;
      if(jpsi->GetNDaughters()!=2)return kFALSE;


      for(Int_t index= jpsi->GetFirstDaughter();index<= jpsi->GetLastDaughter();index++){
         TParticle* temp = (TParticle*)fMCStack->Particle( index );
         switch( temp->GetPdgCode() ) {
         case -11:
            electron =  temp;
            labelelectronChiC = index;
            break;
         case 11:
            positron =  temp;
            labelpositronChiC = index;
            break;
         }
      }
      if( !electron || !positron) return kFALSE;
      if( positron && electron && gamma) return kTRUE;
   }
   return kFALSE;
}



///________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal, Double_t fRapidityShift)
{

   // Selection of reconstructed Meson candidates
   // Use flag IsSignal in order to fill Fill different
   // histograms for Signal and Background
   TH1 *hist=0x0;

   if(IsSignal){hist=hMesonCuts;}
   else{hist=hMesonBGCuts;}

   Int_t cutIndex=0;

   if(hist)hist->Fill(cutIndex);
   cutIndex++;

   // Undefined Rapidity -> Floating Point exception
   if((pi0->E()+pi0->Pz())/(pi0->E()-pi0->Pz())<=0){
      if(hist)hist->Fill(cutIndex);
      cutIndex++;
      return kFALSE;
   }
   else{
      // PseudoRapidity Cut --> But we cut on Rapidity !!!
      cutIndex++;
      if(abs(pi0->Rapidity()-fRapidityShift)>fRapidityCutMeson){
         if(hist)hist->Fill(cutIndex);
         return kFALSE;
      }
   }
   cutIndex++;

   // Opening Angle Cut
   //fOpeningAngle=2*TMath::ATan(0.134/pi0->P());// physical minimum opening angle
   if( pi0->GetOpeningAngle() < fOpeningAngle){
      if(hist)hist->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if ( fAlphaPtDepCut == kTRUE ) {
 
	fAlphaCutMeson = fFAlphaCut->Eval( pi0->Pt() );
   }
   
   
   // Alpha Max Cut
   if(pi0->GetAlpha()>fAlphaCutMeson){
      if(hist)hist->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   // Alpha Min Cut
   if(pi0->GetAlpha()<fAlphaMinCutMeson){
      if(hist)hist->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if (hDCAGammaGammaMesonBefore)hDCAGammaGammaMesonBefore->Fill(pi0->GetDCABetweenPhotons());
   if (hDCARMesonPrimVtxBefore)hDCARMesonPrimVtxBefore->Fill(pi0->GetDCARMotherPrimVtx());

   if (pi0->GetDCABetweenPhotons() > fDCAGammaGammaCut){
      if(hist)hist->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if (pi0->GetDCARMotherPrimVtx() > fDCARMesonPrimVtxCut){
      if(hist)hist->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;


   if (hDCAZMesonPrimVtxBefore)hDCAZMesonPrimVtxBefore->Fill(pi0->GetDCAZMotherPrimVtx());

   if (abs(pi0->GetDCAZMotherPrimVtx()) > fDCAZMesonPrimVtxCut){
      if(hist)hist->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;


   if (hDCAGammaGammaMesonAfter)hDCAGammaGammaMesonAfter->Fill(pi0->GetDCABetweenPhotons());
   if (hDCARMesonPrimVtxAfter)hDCARMesonPrimVtxAfter->Fill(pi0->GetDCARMotherPrimVtx());
   if (hDCAZMesonPrimVtxAfter)hDCAZMesonPrimVtxAfter->Fill(pi0->M(),pi0->GetDCAZMotherPrimVtx());

   if(hist)hist->Fill(cutIndex);
   return kTRUE;
}



///________________________________________________________________________
///________________________________________________________________________
Bool_t AliConversionMesonCuts::UpdateCutString() {
   ///Update the cut string (if it has been created yet)

   if(fCutString && fCutString->GetString().Length() == kNCuts) {
      fCutString->SetString(GetCutNumber());
   } else {
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
   // Initialize Cuts from a given Cut string
   AliInfo(Form("Set Meson Cutnumber: %s",analysisCutSelection.Data()));
   if(analysisCutSelection.Length()!=kNCuts) {
      AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
      return kFALSE;
   }
   if(!analysisCutSelection.IsDigit()){
      AliError("Cut selection contains characters");
      return kFALSE;
   }

   const char *cutSelection = analysisCutSelection.Data();
#define ASSIGNARRAY(i)	fCuts[i] = cutSelection[i] - '0'
   for(Int_t ii=0;ii<kNCuts;ii++){
      ASSIGNARRAY(ii);
   }

   // Set Individual Cuts
   for(Int_t ii=0;ii<kNCuts;ii++){
      if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
   }

   //PrintCuts();
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetCut(cutIds cutID, const Int_t value) {
   ///Set individual cut ID

   //cout << "Updating cut  " << fgkCutNames[cutID] << " (" << cutID << ") to " << value << endl;
   switch (cutID) {
   case kMesonKind:
      if( SetMesonKind(value)) {
         fCuts[kMesonKind] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kchi2MesonCut:
      if( SetChi2MesonCut(value)) {
         fCuts[kchi2MesonCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kalphaMesonCut:
      if( SetAlphaMesonCut(value)) {
         fCuts[kalphaMesonCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kRCut:
      if( SetRCut(value)) {
         fCuts[kRCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kRapidityMesonCut:
      if( SetRapidityMesonCut(value)) {
         fCuts[kRapidityMesonCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kBackgroundScheme:
      if( SetBackgroundScheme(value)) {
         fCuts[kBackgroundScheme] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kDegreesForRotationMethod:
      if( SetNDegreesForRotationMethod(value)) {
         fCuts[kDegreesForRotationMethod] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kNumberOfBGEvents:
      if( SetNumberOfBGEvents(value)) {
         fCuts[kNumberOfBGEvents] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kuseMCPSmearing:
      if( SetMCPSmearing(value)) {
         fCuts[kuseMCPSmearing] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kElecShare:
      if( SetSharedElectronCut(value)) {
         fCuts[kElecShare] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kToCloseV0s:
      if( SetToCloseV0sCut(value)) {
         fCuts[kToCloseV0s] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kDcaGammaGamma:
      if( SetDCAGammaGammaCut(value)) {
         fCuts[kDcaGammaGamma] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kDcaZPrimVtx:
      if( SetDCAZMesonPrimVtxCut(value)) {
         fCuts[kDcaZPrimVtx] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;
   case kDcaRPrimVtx:
      if( SetDCARMesonPrimVtxCut(value)) {
         fCuts[kDcaRPrimVtx] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kNCuts:
      cout << "Error:: Cut id out of range"<< endl;
      return kFALSE;
   }

   cout << "Error:: Cut id " << cutID << " not recognized "<< endl;
   return kFALSE;

}


///________________________________________________________________________
void AliConversionMesonCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMesonKind(Int_t mesonKind){
   // Set Cut
   switch(mesonKind){
   case 0:
      fMesonKind=0;
      break;
   case 1:
      fMesonKind=1;;
      break;
   default:
      cout<<"Warning: Meson kind not defined"<<mesonKind<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetRCut(Int_t rCut){
   // Set Cut
   switch(rCut){
   case 0:
      fMaxR = 180.;
      break;
   case 1:
      fMaxR = 180.;
      break;
   case 2:
      fMaxR = 180.;
      break;
   case 3:
      fMaxR = 70.;
      break;
   case 4:
      fMaxR = 70.;
      break;
   case 5:
      fMaxR = 180.;
      break;
      // High purity cuts for PbPb
   case 9:
      fMaxR = 180.;
      break;
   default:
      cout<<"Warning: rCut not defined"<<rCut<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetChi2MesonCut(Int_t chi2MesonCut){
   // Set Cut
   switch(chi2MesonCut){
   case 0:  // 100.
      fChi2CutMeson = 100.;
      break;
   case 1:  // 50.
      fChi2CutMeson = 50.;
      break;
   case 2:  // 30.
      fChi2CutMeson = 30.;
      break;
   case 3:
      fChi2CutMeson = 200.;
      break;
   case 4:
      fChi2CutMeson = 500.;
      break;
   case 5:
      fChi2CutMeson = 1000.;
      break;
   default:
      cout<<"Warning: Chi2MesonCut not defined "<<chi2MesonCut<<endl;
      return kFALSE;
   }
   return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetAlphaMesonCut(Int_t alphaMesonCut)
{   // Set Cut
   switch(alphaMesonCut){
   case 0:	// 0- 0.7
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 0.7;
      fAlphaPtDepCut = kFALSE;
      break;
   case 1:	// Updated 31 October 2013 before 0.0 - 0.5
      if( fFAlphaCut ) delete fFAlphaCut;
      fFAlphaCut= new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
      fFAlphaCut->SetParameter(0,0.7);
      fFAlphaCut->SetParameter(1,1.2);
      fAlphaMinCutMeson	 =  0.0;
      fAlphaCutMeson	 = -1.0;
      fAlphaPtDepCut = kTRUE;
      break;
   case 2:	// Updated 31 October 2013 before 0.5-1  
      if( fFAlphaCut ) delete fFAlphaCut;
      fFAlphaCut= new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
      fFAlphaCut->SetParameter(0,0.8);
      fFAlphaCut->SetParameter(1,1.2);
      fAlphaMinCutMeson	 =  0.0;
      fAlphaCutMeson	 = -1.0;
      fAlphaPtDepCut = kTRUE;
      break;
   case 3:	// 0.0-1
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 1.;
      fAlphaPtDepCut = kFALSE;
      break;
   case 4:	// 0-0.65
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 0.65;
      fAlphaPtDepCut = kFALSE;
      break;
   case 5:	// 0-0.75
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 0.75;
      fAlphaPtDepCut = kFALSE;
      break;
   case 6:	// 0-0.8
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 0.8;
      fAlphaPtDepCut = kFALSE;
      break;
   case 7:	// 0.0-0.85
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 0.85;
      fAlphaPtDepCut = kFALSE;
      break;
   case 8:	// 0.0-0.6
      fAlphaMinCutMeson	 = 0.0;
      fAlphaCutMeson	 = 0.6;
      fAlphaPtDepCut = kFALSE;
      break;
   case 9: // Updated 11 November 2013 before 0.0 - 0.3
      if( fFAlphaCut ) delete fFAlphaCut;
      fFAlphaCut= new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
      fFAlphaCut->SetParameter(0,0.65);
      fFAlphaCut->SetParameter(1,1.2);
      fAlphaMinCutMeson  =  0.0;
      fAlphaCutMeson     = -1.0;
      fAlphaPtDepCut = kTRUE;
      break;
   default:
      cout<<"Warning: AlphaMesonCut not defined "<<alphaMesonCut<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetRapidityMesonCut(Int_t RapidityMesonCut){
   // Set Cut
   switch(RapidityMesonCut){
   case 0:  //
      fRapidityCutMeson   = 0.9;
      break;
   case 1:  //
      fRapidityCutMeson   = 0.8;
      break;
   case 2:  //
      fRapidityCutMeson   = 0.7;
      break;
   case 3:  //
      fRapidityCutMeson   = 0.6;
      break;
   case 4:  //
      fRapidityCutMeson   = 0.5;
      break;
   case 5:  //
      fRapidityCutMeson   = 0.85;
      break;
   case 6:  //
      fRapidityCutMeson   = 0.75;
      break;
   case 7:  //
      fRapidityCutMeson   = 0.3;
      break;
   case 8:  //
      fRapidityCutMeson   = 0.35;
      break;
   case 9:  //
      fRapidityCutMeson   = 0.4;
      break;
   default:
      cout<<"Warning: RapidityMesonCut not defined "<<RapidityMesonCut<<endl;
      return kFALSE;
   }
   return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetBackgroundScheme(Int_t BackgroundScheme){
   // Set Cut
   switch(BackgroundScheme){
   case 0: //Rotation
      fUseRotationMethodInBG=kTRUE;
      fdoBGProbability=kFALSE;
      break;
   case 1: // mixed event with V0 multiplicity
      fUseRotationMethodInBG=kFALSE;
      fUseTrackMultiplicityForBG=kFALSE;
      fdoBGProbability=kFALSE;
      break;
   case 2: // mixed event with track multiplicity
      fUseRotationMethodInBG=kFALSE;
      fUseTrackMultiplicityForBG=kTRUE;
      fdoBGProbability=kFALSE;
      break;
   case 3: //Rotation
      fUseRotationMethodInBG=kTRUE;
      fdoBGProbability=kTRUE;
      break;
   case 4: //No BG calculation
      cout << "no BG calculation should be done" << endl;
      fUseRotationMethodInBG=kFALSE;
      fdoBGProbability=kFALSE;
      fDoBG=kFALSE;
      break;
   case 5: //Rotation
      fUseRotationMethodInBG=kTRUE;
      fdoBGProbability=kFALSE;
      fBackgroundHandler = 1;
      break;
   case 6: // mixed event with V0 multiplicity
      fUseRotationMethodInBG=kFALSE;
      fUseTrackMultiplicityForBG=kFALSE;
      fdoBGProbability=kFALSE;
      fBackgroundHandler = 1;
      break;
   case 7: // mixed event with track multiplicity
      fUseRotationMethodInBG=kFALSE;
      fUseTrackMultiplicityForBG=kTRUE;
      fdoBGProbability=kFALSE;
      fBackgroundHandler = 1;
      break;
   case 8: //Rotation
      fUseRotationMethodInBG=kTRUE;
      fdoBGProbability=kTRUE;
      fBackgroundHandler = 1;
      break;
   default:
      cout<<"Warning: BackgroundScheme not defined "<<BackgroundScheme<<endl;
      return kFALSE;
   }
   return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetNDegreesForRotationMethod(Int_t DegreesForRotationMethod){
   // Set Cut
   switch(DegreesForRotationMethod){
   case 0:
      fnDegreeRotationPMForBG = 5;
      break;
   case 1:
      fnDegreeRotationPMForBG = 10;
      break;
   case 2:
      fnDegreeRotationPMForBG = 15;
      break;
   case 3:
      fnDegreeRotationPMForBG = 20;
      break;
   default:
      cout<<"Warning: DegreesForRotationMethod not defined "<<DegreesForRotationMethod<<endl;
      return kFALSE;
   }
   fCuts[kDegreesForRotationMethod]=DegreesForRotationMethod;
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetNumberOfBGEvents(Int_t NumberOfBGEvents){
   // Set Cut
   switch(NumberOfBGEvents){
   case 0:
      fNumberOfBGEvents = 5;
      break;
   case 1:
      fNumberOfBGEvents = 10;
      break;
   case 2:
      fNumberOfBGEvents = 15;
      break;
   case 3:
      fNumberOfBGEvents = 20;
      break;
   case 4:
      fNumberOfBGEvents = 2;
      break;
   case 5:
      fNumberOfBGEvents = 50;
      break;
   case 6:
      fNumberOfBGEvents = 80;
      break;
   case 7:
      fNumberOfBGEvents = 100;
      break;
   default:
      cout<<"Warning: NumberOfBGEvents not defined "<<NumberOfBGEvents<<endl;
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetSharedElectronCut(Int_t sharedElec) {

   switch(sharedElec){
   case 0:
      fDoSharedElecCut = kFALSE;
      break;
   case 1:
      fDoSharedElecCut = kTRUE;
      break;
   default:
      cout<<"Warning: Shared Electron Cut not defined "<<sharedElec<<endl;
      return kFALSE;
   }

   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetToCloseV0sCut(Int_t toClose) {

   switch(toClose){
   case 0:
      fDoToCloseV0sCut = kFALSE;
      fminV0Dist = 250;
      break;
   case 1:
      fDoToCloseV0sCut = kTRUE;
      fminV0Dist = 1;
      break;
   case 2:
      fDoToCloseV0sCut = kTRUE;
      fminV0Dist = 2;
      break;
   case 3:
      fDoToCloseV0sCut = kTRUE;
      fminV0Dist = 3;
      break;
   default:
      cout<<"Warning: Shared Electron Cut not defined "<<toClose<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMCPSmearing(Int_t useMCPSmearing)
{// Set Cut
   switch(useMCPSmearing){
   case 0:
      fUseMCPSmearing=0;
      fPBremSmearing=1.;
      fPSigSmearing=0.;
      fPSigSmearingCte=0.;
      break;
   case 1:
      fUseMCPSmearing=1;
      fPBremSmearing=1.0e-14;
      fPSigSmearing=0.;
      fPSigSmearingCte=0.;
      break;
   case 2:
      fUseMCPSmearing=1;
      fPBremSmearing=1.0e-15;
      fPSigSmearing=0.0;
      fPSigSmearingCte=0.;
      break;
   case 3:
      fUseMCPSmearing=1;
      fPBremSmearing=1.;
      fPSigSmearing=0.003;
      fPSigSmearingCte=0.002;
      break;
   case 4:
      fUseMCPSmearing=1;
      fPBremSmearing=1.;
      fPSigSmearing=0.003;
      fPSigSmearingCte=0.007;
      break;
   case 5:
      fUseMCPSmearing=1;
      fPBremSmearing=1.;
      fPSigSmearing=0.003;
      fPSigSmearingCte=0.016;
      break;
   case 6:
      fUseMCPSmearing=1;
      fPBremSmearing=1.;
      fPSigSmearing=0.007;
      fPSigSmearingCte=0.016;
      break;
   case 7:
      fUseMCPSmearing=1;
      fPBremSmearing=1.0e-16;
      fPSigSmearing=0.0;
      fPSigSmearingCte=0.;
      break;
   case 8:
      fUseMCPSmearing=1;
      fPBremSmearing=1.;
      fPSigSmearing=0.007;
      fPSigSmearingCte=0.014;
      break;
   case 9:
      fUseMCPSmearing=1;
      fPBremSmearing=1.;
      fPSigSmearing=0.007;
      fPSigSmearingCte=0.011;
      break;

   default:
      AliError("Warning: UseMCPSmearing not defined");
      return kFALSE;
   }
   return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetDCAGammaGammaCut(Int_t DCAGammaGamma){
   // Set Cut
   switch(DCAGammaGamma){
   case 0:  //
      fDCAGammaGammaCut   = 1000;
      break;
   case 1:  //
      fDCAGammaGammaCut   = 10;
      break;
   case 2:  //
      fDCAGammaGammaCut   = 5;
      break;
   case 3:  //
      fDCAGammaGammaCut   = 4;
      break;
   case 4:  //
      fDCAGammaGammaCut   = 3;
      break;
   case 5:  //
      fDCAGammaGammaCut   = 2.5;
      break;
   case 6:  //
      fDCAGammaGammaCut   = 2;
      break;
   case 7:  //
      fDCAGammaGammaCut   = 1.5;
      break;
   case 8:  //
      fDCAGammaGammaCut   = 1;
      break;
   case 9:  //
      fDCAGammaGammaCut   = 0.5;
      break;
   default:
      cout<<"Warning: DCAGammaGamma not defined "<<DCAGammaGamma<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetDCAZMesonPrimVtxCut(Int_t DCAZMesonPrimVtx){
   // Set Cut
   switch(DCAZMesonPrimVtx){
   case 0:  //
      fDCAZMesonPrimVtxCut   = 1000;
      break;
   case 1:  //
      fDCAZMesonPrimVtxCut   = 10;
      break;
   case 2:  //
      fDCAZMesonPrimVtxCut   = 5;
      break;
   case 3:  //
      fDCAZMesonPrimVtxCut   = 4;
      break;
   case 4:  //
      fDCAZMesonPrimVtxCut   = 3;
      break;
   case 5:  //
      fDCAZMesonPrimVtxCut   = 2.5;
      break;
   case 6:  //
      fDCAZMesonPrimVtxCut   = 2;
      break;
   case 7:  //
      fDCAZMesonPrimVtxCut   = 1.5;
      break;
   case 8:  //
      fDCAZMesonPrimVtxCut   = 1;
      break;
   case 9:  //
      fDCAZMesonPrimVtxCut   = 0.5;
      break;
   default:
      cout<<"Warning: DCAZMesonPrimVtx not defined "<<DCAZMesonPrimVtx<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::SetDCARMesonPrimVtxCut(Int_t DCARMesonPrimVtx){
   // Set Cut
   switch(DCARMesonPrimVtx){
   case 0:  //
      fDCARMesonPrimVtxCut   = 1000;
      break;
   case 1:  //
      fDCARMesonPrimVtxCut   = 10;
      break;
   case 2:  //
      fDCARMesonPrimVtxCut   = 5;
      break;
   case 3:  //
      fDCARMesonPrimVtxCut   = 4;
      break;
   case 4:  //
      fDCARMesonPrimVtxCut   = 3;
      break;
   case 5:  //
      fDCARMesonPrimVtxCut   = 2.5;
      break;
   case 6:  //
      fDCARMesonPrimVtxCut   = 2;
      break;
   case 7:  //
      fDCARMesonPrimVtxCut   = 1.5;
      break;
   case 8:  //
      fDCARMesonPrimVtxCut   = 1;
      break;
   case 9:  //
      fDCARMesonPrimVtxCut   = 0.5;
      break;
   default:
      cout<<"Warning: DCARMesonPrimVtx not defined "<<DCARMesonPrimVtx<<endl;
      return kFALSE;
   }
   return kTRUE;
}


///________________________________________________________________________
TString AliConversionMesonCuts::GetCutNumber(){
   // returns TString with current cut number
   TString a(kNCuts);
   for(Int_t ii=0;ii<kNCuts;ii++){
      a.Append(Form("%d",fCuts[ii]));
   }
   return a;
}

///________________________________________________________________________
void AliConversionMesonCuts::FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0){

   Int_t posLabel = photon->GetTrackLabelPositive();
   Int_t negLabel = photon->GetTrackLabelNegative();

   fElectronLabelArray[nV0*2] = posLabel;
   fElectronLabelArray[(nV0*2)+1] = negLabel;
}

///________________________________________________________________________
Bool_t AliConversionMesonCuts::RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s){

   Int_t posLabel = photon->GetTrackLabelPositive();
   Int_t negLabel = photon->GetTrackLabelNegative();

   for(Int_t i = 0; i<nV0s*2;i++){
      if(i==nV0*2)     continue;
      if(i==(nV0*2)+1) continue;
      if(fElectronLabelArray[i] == posLabel){
         return kFALSE;}
      if(fElectronLabelArray[i] == negLabel){
         return kFALSE;}
   }

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionMesonCuts::RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0){
   Double_t posX = photon->GetConversionX();
   Double_t posY = photon->GetConversionY();
   Double_t posZ = photon->GetConversionZ();

   for(Int_t i = 0;i<photons->GetEntries();i++){
      if(nV0 == i) continue;
      AliAODConversionPhoton *photonComp = (AliAODConversionPhoton*) photons->At(i);
      Double_t posCompX = photonComp->GetConversionX();
      Double_t posCompY = photonComp->GetConversionY();
      Double_t posCompZ = photonComp->GetConversionZ();

      Double_t dist = pow((posX - posCompX),2)+pow((posY - posCompY),2)+pow((posZ - posCompZ),2);

      if(dist < fminV0Dist*fminV0Dist){
         if(photon->GetChi2perNDF() < photonComp->GetChi2perNDF()) return kTRUE;
         else {
            return kFALSE;}
      }

   }
   return kTRUE;
}

///________________________________________________________________________
void AliConversionMesonCuts::SmearParticle(AliAODConversionPhoton* photon)
{

   if (photon==NULL) return;
   Double_t facPBrem = 1.;
   Double_t facPSig = 0.;

   Double_t phi=0.;
   Double_t theta=0.;
   Double_t P=0.;


   P=photon->P();
   phi=photon->Phi();
   if( photon->P()!=0){
      theta=acos( photon->Pz()/ photon->P());
   }

   if( fPSigSmearing != 0. || fPSigSmearingCte!=0. ){
      facPSig = TMath::Sqrt(fPSigSmearingCte*fPSigSmearingCte+fPSigSmearing*fPSigSmearing*P*P)*fRandom.Gaus(0.,1.);
   }

   if( fPBremSmearing != 1.){
      if(fBrem!=NULL){
         facPBrem = fBrem->GetRandom();
      }
   }

   photon->SetPx(facPBrem* (1+facPSig)* P*sin(theta)*cos(phi)) ;
   photon->SetPy(facPBrem* (1+facPSig)* P*sin(theta)*sin(phi)) ;
   photon->SetPz(facPBrem* (1+facPSig)* P*cos(theta)) ;
   photon->SetE(photon->P());
}
