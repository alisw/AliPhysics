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

#include "AliConversionCuts.h"

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
#include "AliLog.h"

class iostream;

using namespace std;

ClassImp(AliConversionCuts)


const char* AliConversionCuts::fgkCutNames[AliConversionCuts::kNCuts] = {
"GoodId",
"V0FinderType",
"eProbCut",
"ededxSigmaCut",
"pidedxSigmaCut",
"piMomdedxSigmaCut",
"Chi2GammaCut",
"SinglePtCut",
"ClsTPCCut",
"EtaCut",
"Chi2MesonCut",
"LowPRejectionSigmaCut",
"QtMaxCut",
"piMaxMomdedxSigmaCut",
"AlphaMesonCut",
"MinRCut",
"RapidityMesonCut",
"BackgroundScheme",
"DegreesForRotationMethod",
"NumberOfRotations",
"RemovePileUp",
"SelectV0AND",
"MultiplicityMethod",
"HeavyIon",
"CentralityMin",
"CentralityMax",
"TOFelectronPID",
"UseMCPSmearing",
"DoPhotonAsymmetryCut",
"PsiPair",
"CosinePointingAngle",
"SharedElectronCuts",
"RejectToCloseV0s",
};


//________________________________________________________________________
AliConversionCuts::AliConversionCuts(const char *name,const char *title) : AliAnalysisCuts(name,title),
    fHistograms(NULL),
    fPIDResponse(NULL),
    fEventQuality(-1),
    fMaxR(200),
    fMinR(0),
    fEtaCut(0.9),
    fEtaCutMin(-0.1),
    fPtCut(0.02),
    fSinglePtCut(0),
    fMaxZ(1000),
    fMinClsTPC(0.),
    fMinClsTPCToF(0.),
    fLineCutZRSlope(0.),
    fLineCutZValue(0),
    fLineCutZRSlopeMin(0.),
    fLineCutZValueMin(0),
    fChi2CutConversion(1000),
    fChi2CutMeson(1000),
    fPIDProbabilityCutNegativeParticle(0),
    fPIDProbabilityCutPositiveParticle(0),
    fDodEdxSigmaCut(kTRUE),
    fDoTOFsigmaCut(kFALSE), // RRnewTOF
    fPIDTRDEfficiency(1),
    fDoTRDPID(kFALSE),
    fPIDnSigmaAboveElectronLine(100),
    fPIDnSigmaBelowElectronLine(-100),
    fTofPIDnSigmaAboveElectronLine(100), // RRnewTOF
    fTofPIDnSigmaBelowElectronLine(-100), // RRnewTOF
    fPIDnSigmaAbovePionLine(0),
    fPIDnSigmaAbovePionLineHighPt(-100),
    fPIDMinPnSigmaAbovePionLine(0),
    fPIDMaxPnSigmaAbovePionLine(0),
    fDoKaonRejectionLowP(kFALSE),
    fDoProtonRejectionLowP(kFALSE),
    fDoPionRejectionLowP(kFALSE),
    fPIDnSigmaAtLowPAroundKaonLine(0),
    fPIDnSigmaAtLowPAroundProtonLine(0),
    fPIDnSigmaAtLowPAroundPionLine(0),
    fPIDMinPKaonRejectionLowP(0),
    fPIDMinPProtonRejectionLowP(0),
    fPIDMinPPionRejectionLowP(0),
    fDoQtGammaSelection(kTRUE),
    fDoHighPtQtGammaSelection(kFALSE), // RRnew
    fQtMax(100),
    fHighPtQtMax(0.), // RRnew
    fPtBorderForQt(0), // RRnew
    fXVertexCut(0.),
    fYVertexCut(0.),
    fZVertexCut(0.),
    fNSigmaMass(0.),
    fUseEtaMinCut(kFALSE),
    fUseOnFlyV0Finder(kTRUE),
    fDoPhotonAsymmetryCut(kTRUE),
    fMinPPhotonAsymmetryCut(100.),
    fMinPhotonAsymmetry(0.),
    fIsHeavyIon(kFALSE),
    fDetectorCentrality(0),
    fModCentralityClass(0),
    fMaxVertexZ(10),
    fCentralityMin(0),
    fCentralityMax(0),
    fUseCorrectedTPCClsInfo(kFALSE),
    fUseTOFpid(kFALSE),
    fAlphaMinCutMeson(0),
    fAlphaCutMeson(1),
    fRapidityCutMeson(1),
    fUseRotationMethodInBG(kFALSE),
    fdoBGProbability(kFALSE),
    fUseTrackMultiplicityForBG(kFALSE),
    fnDegreeRotationPMForBG(0),
    fnumberOfRotationEventsForBG(0),
    fUseMCPSmearing(kFALSE),
    fPBremSmearing(0),
    fPSigSmearing(0),
    fPSigSmearingCte(0),
    fBrem(NULL),
    fMultiplicityMethod(0),
    fSelectV0AND(kFALSE),
    fRemovePileUp(kFALSE),
    fOpeningAngle(0.005),
    fPsiPairCut(10000),
    fCosPAngleCut(10000),
    fDoToCloseV0sCut(kFALSE),
    fminV0Dist(200.),
    fDoSharedElecCut(kFALSE),
    fOfflineTriggerMask(0),
    fRandom(0),
    fElectronLabelArray(NULL),
    fCutString(NULL),
    hdEdxCuts(NULL),
    hTPCdEdxbefore(NULL),
    hTPCdEdxafter(NULL),
    hTOFbefore(NULL),
    hTOFafter(NULL),
    hTrackCuts(NULL),
    hPhotonCuts(NULL),
    hInvMassbefore(NULL),
    hArmenterosbefore(NULL),
    hInvMassafter(NULL),
    hArmenterosafter(NULL),
    hAcceptanceCuts(NULL),
    hCutIndex(NULL),
    hV0EventCuts(NULL),
    hCentrality(NULL),
    hVertexZ(NULL),
    hTriggerClass(NULL),
    hMesonCuts(NULL),
    hMesonBGCuts(NULL)
{
    InitPIDResponse();
    for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
    fCutString=new TObjString((GetCutNumber()).Data());

    if (fBrem == NULL){
       fBrem = new TF1("fBrem","pow(-log(x),[0]/log(2.0)-1.0)/TMath::Gamma([0]/log(2.0))",0.00001,0.999999999);
       // tests done with 1.0e-14
       fBrem->SetParameter(0,fPBremSmearing);
       fBrem->SetNpx(100000);
    }

    fElectronLabelArray = new Int_t[500];
}

//________________________________________________________________________
AliConversionCuts::~AliConversionCuts() {
    // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  // 	delete fHistograms;
  // fHistograms = NULL;
   if(fCutString != NULL){
      delete fCutString;
      fCutString = NULL;
   }
   if(fBrem){
      delete fBrem;
      fBrem = NULL;
   }
   if(fElectronLabelArray){
      delete fElectronLabelArray;
      fElectronLabelArray = NULL;
   }

}

//________________________________________________________________________
void AliConversionCuts::InitCutHistograms(TString name, Bool_t preCut){

    // Initialize Cut Histograms for QA (only initialized and filled if function is called)

    if(fHistograms != NULL){
	delete fHistograms;
	fHistograms=NULL;
    }
    if(fHistograms==NULL){
	fHistograms=new TList();
	if(name=="")fHistograms->SetName(Form("ConvCuts_%s",GetCutNumber().Data()));
	else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
    }

    // IsPhotonSelected
    hCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",10,-0.5,9.5);
    hCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
    hCutIndex->GetXaxis()->SetBinLabel(kOnFly+1,"onfly");
    hCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
    hCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"dEdx");
    hCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
    hCutIndex->GetXaxis()->SetBinLabel(kConvPointFail+1,"ConvPoint fail");
    hCutIndex->GetXaxis()->SetBinLabel(kPhotonCuts+1,"PhotonCuts");
    hCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
    fHistograms->Add(hCutIndex);

    // Track Cuts
    hTrackCuts=new TH1F(Form("TrackCuts %s",GetCutNumber().Data()),"TrackCuts",10,-0.5,9.5);
    hTrackCuts->GetXaxis()->SetBinLabel(1,"in");
    hTrackCuts->GetXaxis()->SetBinLabel(2,"likesign");
    hTrackCuts->GetXaxis()->SetBinLabel(3,"ntpccl");
    hTrackCuts->GetXaxis()->SetBinLabel(4,"acceptance");
    hTrackCuts->GetXaxis()->SetBinLabel(5,"singlept");
    hTrackCuts->GetXaxis()->SetBinLabel(6,"TPCrefit");
    hTrackCuts->GetXaxis()->SetBinLabel(7,"kink");
    hTrackCuts->GetXaxis()->SetBinLabel(8,"out");
    fHistograms->Add(hTrackCuts);

    // Photon Cuts
    hPhotonCuts=new TH1F(Form("PhotonCuts %s",GetCutNumber().Data()),"PhotonCuts",12,-0.5,11.5);
    hPhotonCuts->GetXaxis()->SetBinLabel(1,"in");
    hPhotonCuts->GetXaxis()->SetBinLabel(2,"qtcut");
    hPhotonCuts->GetXaxis()->SetBinLabel(3,"chi2");
    hPhotonCuts->GetXaxis()->SetBinLabel(4,"acceptance");
    hPhotonCuts->GetXaxis()->SetBinLabel(5,"asymmetry");
    hPhotonCuts->GetXaxis()->SetBinLabel(6,"pidprob");
    hPhotonCuts->GetXaxis()->SetBinLabel(7,"cortpcclinfo");
    hPhotonCuts->GetXaxis()->SetBinLabel(8,"PsiPair");
    hPhotonCuts->GetXaxis()->SetBinLabel(9,"CosPAngle");
    hPhotonCuts->GetXaxis()->SetBinLabel(10,"out");
    fHistograms->Add(hPhotonCuts);

    if(preCut){
       hInvMassbefore=new TH1F(Form("InvMass_before %s",GetCutNumber().Data()),"InvMass_before",100,0,0.3);
       fHistograms->Add(hInvMassbefore);
       hArmenterosbefore=new TH2F(Form("Armenteros_before %s",GetCutNumber().Data()),"Armenteros_before",200,-1,1,250,0,0.25);
       fHistograms->Add(hArmenterosbefore);
    }
    hInvMassafter=new TH1F(Form("InvMass_after %s",GetCutNumber().Data()),"InvMass_after",100,0,0.3);
    fHistograms->Add(hInvMassafter);
    hArmenterosafter=new TH2F(Form("Armenteros_after %s",GetCutNumber().Data()),"Armenteros_after",200,-1,1,250,0,0.25);
    fHistograms->Add(hArmenterosafter);

    hAcceptanceCuts=new TH1F(Form("PhotonAcceptanceCuts %s",GetCutNumber().Data()),"PhotonAcceptanceCuts",10,-0.5,9.5);
    hAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(2,"maxR");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(3,"minR");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(4,"line");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(5,"maxZ");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(6,"eta");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(7,"minpt");
    hAcceptanceCuts->GetXaxis()->SetBinLabel(8,"out");
    fHistograms->Add(hAcceptanceCuts);

    // dEdx Cuts
    hdEdxCuts=new TH1F(Form("dEdxCuts %s",GetCutNumber().Data()),"dEdxCuts",10,-0.5,9.5);
    hdEdxCuts->GetXaxis()->SetBinLabel(1,"in");
    hdEdxCuts->GetXaxis()->SetBinLabel(2,"TPCelectron");
    hdEdxCuts->GetXaxis()->SetBinLabel(3,"TPCpion");
    hdEdxCuts->GetXaxis()->SetBinLabel(4,"TPCpionhighp");
    hdEdxCuts->GetXaxis()->SetBinLabel(5,"TPCkaonlowprej");
    hdEdxCuts->GetXaxis()->SetBinLabel(6,"TPCprotonlowprej");
    hdEdxCuts->GetXaxis()->SetBinLabel(7,"TPCpionlowprej");
    hdEdxCuts->GetXaxis()->SetBinLabel(8,"TOFelectron");
    hdEdxCuts->GetXaxis()->SetBinLabel(9,"TRDelectron");
    hdEdxCuts->GetXaxis()->SetBinLabel(10,"out");
    fHistograms->Add(hdEdxCuts);
    
    TAxis *AxisBeforedEdx = NULL;
    TAxis *AxisBeforeTOF = NULL;
    if(preCut){
       hTPCdEdxbefore=new TH2F(Form("Gamma_dEdx_before %s",GetCutNumber().Data()),"dEdx Gamma before" ,150,0.05,20,400,-10,10);
       fHistograms->Add(hTPCdEdxbefore);
       AxisBeforedEdx = hTPCdEdxbefore->GetXaxis(); 

       hTOFbefore=new TH2F(Form("Gamma_TOF_before %s",GetCutNumber().Data()),"TOF Gamma before" ,150,0.05,20,400,-6,10);
       fHistograms->Add(hTOFbefore);
       AxisBeforeTOF = hTOFbefore->GetXaxis(); 
    }
    hTPCdEdxafter=new TH2F(Form("Gamma_dEdx_after %s",GetCutNumber().Data()),"dEdx Gamma after" ,150,0.05,20,400, -10,10);
    fHistograms->Add(hTPCdEdxafter);

    hTOFafter=new TH2F(Form("Gamma_TOF_after %s",GetCutNumber().Data()),"TOF Gamma after" ,150,0.05,20,400,-6,10);
    fHistograms->Add(hTOFafter);

    TAxis *AxisAfter = hTPCdEdxafter->GetXaxis(); 
    Int_t bins = AxisAfter->GetNbins();
    Double_t from = AxisAfter->GetXmin();
    Double_t to = AxisAfter->GetXmax();
    Double_t *newBins = new Double_t[bins+1];
    newBins[0] = from;
    Double_t factor = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
    AxisAfter->Set(bins, newBins);
    AxisAfter = hTOFafter->GetXaxis(); 
    AxisAfter->Set(bins, newBins);
    if(preCut){
       AxisBeforedEdx->Set(bins, newBins);
       AxisBeforeTOF->Set(bins, newBins);
    }
    delete [] newBins;
        
    // Event Cuts and Info
    if(preCut){
       hV0EventCuts=new TH1F(Form("ESD_EventCuts %s",GetCutNumber().Data()),"Event Cuts",10,-0.5,9.5);
       hV0EventCuts->GetXaxis()->SetBinLabel(1,"in");
       hV0EventCuts->GetXaxis()->SetBinLabel(2,"OfflineTrigger");
       hV0EventCuts->GetXaxis()->SetBinLabel(3,"VertexZ");
       hV0EventCuts->GetXaxis()->SetBinLabel(4,"nvtxcontr");
       hV0EventCuts->GetXaxis()->SetBinLabel(5,"pileup");
       hV0EventCuts->GetXaxis()->SetBinLabel(6,"centrsel");
       hV0EventCuts->GetXaxis()->SetBinLabel(7,"out");
       fHistograms->Add(hV0EventCuts);
       
       hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",100,0,100);
       fHistograms->Add(hCentrality);
       hVertexZ=new TH1F(Form("VertexZ %s",GetCutNumber().Data()),"VertexZ",1000,-50,50);
       fHistograms->Add(hVertexZ);
       
       hTriggerClass= new TH1F(Form("OfflineTrigger %s",GetCutNumber().Data()),"OfflineTrigger",4,-0.5,3.5);
       hTriggerClass->GetXaxis()->SetBinLabel(1,"kAny");
       hTriggerClass->GetXaxis()->SetBinLabel(2,"kMB");
       hTriggerClass->GetXaxis()->SetBinLabel(3,"kCentral");
       hTriggerClass->GetXaxis()->SetBinLabel(4,"kSemiCentral");
       fHistograms->Add(hTriggerClass);
    }
    
    // Meson Cuts
    hMesonCuts=new TH1F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts",10,-0.5,9.5);
    hMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    hMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    hMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    hMesonCuts->GetXaxis()->SetBinLabel(4,"opening angle");
    hMesonCuts->GetXaxis()->SetBinLabel(5,"alpha max");
    hMesonCuts->GetXaxis()->SetBinLabel(6,"alpha min");
    hMesonCuts->GetXaxis()->SetBinLabel(7,"out");
    fHistograms->Add(hMesonCuts);

    hMesonBGCuts=new TH1F(Form("MesonBGCuts %s",GetCutNumber().Data()),"MesonBGCuts",10,-0.5,9.5);
    hMesonBGCuts->GetXaxis()->SetBinLabel(1,"in");
    hMesonBGCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    hMesonBGCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    hMesonBGCuts->GetXaxis()->SetBinLabel(4,"opening angle");
    hMesonBGCuts->GetXaxis()->SetBinLabel(5,"alpha max");
    hMesonBGCuts->GetXaxis()->SetBinLabel(6,"alpha min");
    hMesonBGCuts->GetXaxis()->SetBinLabel(7,"out");
    fHistograms->Add(hMesonBGCuts);
}

//________________________________________________________________________
Bool_t AliConversionCuts::InitPIDResponse(){
    // Set Pointer to AliPIDResponse

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man) { 
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(fPIDResponse)return kTRUE;
    
  }
  
  return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::EventIsSelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent){
    // Process Event Selection

    Int_t cutindex=0;
    if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
    cutindex++;
    
    // Check for MC event
    if(fMCEvent){
       // Check if MC event is correctly loaded
       AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
       if (!mcHandler){
          fEventQuality = 2;
          return kFALSE;
       }
       if (!mcHandler->InitOk() ){
          fEventQuality = 2;
          return kFALSE;
       }
       if (!mcHandler->TreeK() ){
          fEventQuality = 2;
          return kFALSE;
       }
       if (!mcHandler->TreeTR() ) {
          fEventQuality = 2;
          return kFALSE;
       }
    }
    
    // Event Trigger

    if(!IsTriggerSelected()){
	if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
        fEventQuality = 3;
	return kFALSE;
    }
    cutindex++;

    // Z Vertex Position Cut
    if(!VertexZCut(fInputEvent)){
	if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
        fEventQuality = 4;
        return kFALSE;
    }
    cutindex++;

    // Number of Contributors Cut
    if(GetNumberOfContributorsVtx(fInputEvent)<=0) {
	if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
        fEventQuality = 5;
	return kFALSE;
    }
    cutindex++;

    // Pile Up Rejection

    if(fRemovePileUp){
	if(fInputEvent->IsPileupFromSPD(3,0.8,3.,2.,5.)){
	    if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
            fEventQuality = 6;
            return kFALSE;
	}
    }
    cutindex++;

    // Centrality Selection
    if(!IsCentralitySelected(fInputEvent)){
	if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
        fEventQuality = 1;
	return kFALSE;
    }
    cutindex++;

    // Fill Event Histograms
    if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
    if(hVertexZ)hVertexZ->Fill(fInputEvent->GetPrimaryVertex()->GetZ());
    if(hCentrality)hCentrality->Fill(GetCentrality(fInputEvent));

    fEventQuality = 0;
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PhotonIsSelectedMC(TParticle *particle,AliStack *fMCStack,Bool_t checkForConvertedGamma){
    // MonteCarlo Photon Selection

    if(!fMCStack)return kFALSE;

    if (particle->GetPdgCode() == 22){
              
       	if(particle->R() > fMaxR)	return kFALSE;
	if(TMath::Abs(particle->Eta())> fEtaCut || TMath::Abs(particle->Eta())< fEtaCutMin)	return kFALSE;

	if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
	    return kFALSE; // no photon as mothers!
	}

	if(particle->GetMother(0) >= fMCStack->GetNprimary()){
	    return kFALSE; // the gamma has a mother, and it is not a primary particle
	}

        if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

	// looking for conversion gammas (electron + positron from pairbuilding (= 5) )
	TParticle* ePos = NULL;
	TParticle* eNeg = NULL;

	if(particle->GetNDaughters() >= 2){
	    for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
		TParticle *tmpDaughter = fMCStack->Particle(daughterIndex);
		if(tmpDaughter->GetUniqueID() == 5){
		    if(tmpDaughter->GetPdgCode() == 11){
			eNeg = tmpDaughter;
		    } else if(tmpDaughter->GetPdgCode() == -11){
			ePos = tmpDaughter;
		    }
		}
	    }
	}

	if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
	    return kFALSE;
	}
        
        if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
           return kFALSE; // no reconstruction below the Pt cut
        }
	
        if( TMath::Abs(ePos->Eta())> fEtaCut || TMath::Abs(ePos->Eta())< fEtaCutMin || 
            TMath::Abs(eNeg->Eta())> fEtaCut || TMath::Abs(eNeg->Eta())< fEtaCutMin ) {
           return kFALSE;
        }
	
        if(ePos->R()>fMaxR){
           return kFALSE; // cuts on distance from collision point
        }

        if(TMath::Abs(ePos->Vz()) > fMaxZ){
           return kFALSE;	 // outside material
        }
        if(TMath::Abs(eNeg->Vz()) > fMaxZ){
           return kFALSE;	 // outside material
        }

        if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
           return kFALSE;  // line cut to exclude regions where we do not reconstruct
        } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
           return kFALSE;
        }
	
        if( eNeg->R() <= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
           return kFALSE; // line cut to exclude regions where we do not reconstruct
        } else if ( fEtaCutMin != -0.1 &&   eNeg->R() >= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
           return kFALSE;
        }
	
        return kTRUE;
	//if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
    }
    return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionCuts::MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack,Bool_t bMCDaughtersInAcceptance){
    // Returns true for all pions within acceptance cuts for decay into 2 photons
    // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

    if(!fMCStack)return kFALSE;
    
    if(fMCMother->GetPdgCode()==111 || fMCMother->GetPdgCode()==221){
       
       if(fMCMother->R()>fMaxR)	return kFALSE; // cuts on distance from collision point

       Double_t rapidity = 10.;
       if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
          rapidity=8.;
       }
       else{
          rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())));
       }	
       
       // Rapidity Cut
       if(TMath::Abs(rapidity)>fRapidityCutMeson)return kFALSE;

	// Select only -> 2y decay channel
	if(fMCMother->GetNDaughters()!=2)return kFALSE;

	for(Int_t i=0;i<2;i++){
	    TParticle *MDaughter=fMCStack->Particle(fMCMother->GetDaughter(i));

	    // Is Daughter a Photon?
	    if(MDaughter->GetPdgCode()!=22)return kFALSE;
            // Is Photon in Acceptance?
	    if(bMCDaughtersInAcceptance){
		if(!PhotonIsSelectedMC(MDaughter,fMCStack)){return kFALSE;}
	    }
	}
	return kTRUE;
    }
    return kFALSE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event)
{   // Specific Photon Cuts

    Int_t cutIndex = 0;
    if(hPhotonCuts)hPhotonCuts->Fill(cutIndex);
    cutIndex++;

    // Fill Histos before Cuts
    if(hInvMassbefore)hInvMassbefore->Fill(photon->GetMass());
    if(hArmenterosbefore)hArmenterosbefore->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());

    // Gamma selection based on QT from Armenteros
    if(fDoQtGammaSelection == kTRUE){
	if(!ArmenterosQtCut(photon)){
	    if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //1
	    return kFALSE;
	}
    }
    cutIndex++; //2

    // Chi Cut
    if(photon->GetChi2perNDF() > fChi2CutConversion || photon->GetChi2perNDF() <=0){
	{
	    if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //2
	    return kFALSE;
	}
    }
    cutIndex++;//3

    // Reconstruction Acceptance Cuts
    if(!AcceptanceCuts(photon)){
	if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //3
	return kFALSE;
    }

    cutIndex++; //4
    // Asymmetry Cut
    if(fDoPhotonAsymmetryCut == kTRUE){
	if(!AsymmetryCut(photon,event)){
	    if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //4
	    return kFALSE;
	}
    }

    //Check the pid probability
    cutIndex++; //5
    if(!PIDProbabilityCut(photon, event)) {
	if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //5
	return kFALSE;
    }

    cutIndex++; //6
    if(!CorrectedTPCClusterCut(photon, event)) {
	if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //6
	return kFALSE;
    }


    cutIndex++; //7
    if(!PsiPairCut(photon)) {
	  if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //7
	  return kFALSE;
    }

    cutIndex++; //8
    if(!CosinePAngleCut(photon, event)) {
	  if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //8
	  return kFALSE;
    }

    cutIndex++; //9
    if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //9

    // Histos after Cuts
    if(hInvMassafter)hInvMassafter->Fill(photon->GetMass());
    if(hArmenterosafter)hArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());


    return kTRUE;

}

///________________________________________________________________________
Bool_t AliConversionCuts::CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event)
{   //Cut on corrected TPC Cluster Info

    AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
    AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

    if(!negTrack||!posTrack)return kFALSE;

    Double_t negclsToF=0;

    if (!fUseCorrectedTPCClsInfo ){
	if(negTrack->GetTPCNclsF()!=0){
	    negclsToF = (Double_t)negTrack->GetNcls(1)/(Double_t)negTrack->GetTPCNclsF();}// Ncluster/Nfindablecluster
    }
    else {
	negclsToF = negTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
    }

    Double_t posclsToF = 0.;
    if (!fUseCorrectedTPCClsInfo ){
	if(posTrack->GetTPCNclsF()!=0	){
	    posclsToF = (Double_t)posTrack->GetNcls(1)/(Double_t)posTrack->GetTPCNclsF();
	}
    }else{
	posclsToF = posTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
    }

    if( negclsToF < fMinClsTPCToF ||	posclsToF < fMinClsTPCToF ){
    return kFALSE;
    }

    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PhotonIsSelected(AliConversionPhotonBase *photon, AliVEvent * event)
{
    //Selection of Reconstructed Photons

    FillPhotonCutIndex(kPhotonIn);

    // Get Tracks
    AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
    AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

    if(!negTrack || !posTrack) {
	FillPhotonCutIndex(kNoTracks);
	return kFALSE;
    }

    // dEdx Cuts
    if(!dEdxCuts(negTrack) || !dEdxCuts(posTrack)) {
	FillPhotonCutIndex(kdEdxCuts);
	return kFALSE;
    }

    // Track Cuts
    if(!TracksAreSelected(negTrack, posTrack)){
	FillPhotonCutIndex(kTrackCuts);
	return kFALSE;
    }

    // Photon Cuts
    if(!PhotonCuts(photon,event)){
	FillPhotonCutIndex(kPhotonCuts);
	return kFALSE;
    }

    // Photon passed cuts
    FillPhotonCutIndex(kPhotonOut);
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal)
{
    // Selection of reconstructed Meson candidates
    // Use flag IsSignal in order to fill Fill different
    // histograms for Signal and Background
    TH1 *hist=0x0;

    if(IsSignal){hist=hMesonCuts;
    }
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
        if(TMath::Abs(pi0->Rapidity())>fRapidityCutMeson){
        //if(TMath::Abs(pi0->PseudoRapidity())>fRapidityCutMeson){
	    if(hist)hist->Fill(cutIndex);
	    return kFALSE;
	}
    }
    cutIndex++;

    // Opening Angle Cut
    //fOpeningAngle=2*TMath::ATan(0.134/pi0->P());// physical minimum opening angle
    if(pi0->GetOpeningAngle()<fOpeningAngle){
	if(hist)hist->Fill(cutIndex);
	return kFALSE;
    }
    cutIndex++;

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

    if(hist)hist->Fill(cutIndex);
    return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionCuts::ArmenterosQtCut(AliConversionPhotonBase *photon)
{   // Armenteros Qt Cut

  if(fDoHighPtQtGammaSelection){
	if(photon->GetPhotonPt() < fPtBorderForQt){
	  if(photon->GetArmenterosQt()>fQtMax){
		return kFALSE;
	  }
	} else {
	  if(photon->GetArmenterosQt()>fHighPtQtMax){
		return kFALSE;
	  }
	}
  } else {

	if(photon->GetArmenterosQt()>fQtMax){
	  return kFALSE;
	}
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionCuts::AcceptanceCuts(AliConversionPhotonBase *photon) {
    // Exclude certain areas for photon reconstruction

    Int_t cutIndex=0;
    if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
    cutIndex++;

    if(photon->GetConversionRadius()>fMaxR){ // cuts on distance from collision point
	if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
	return kFALSE;
    }
    cutIndex++;

    if(photon->GetConversionRadius()<fMinR){ // cuts on distance from collision point
	if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
	return kFALSE;
    }
    cutIndex++;

    if(photon->GetConversionRadius() <= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlope)-fLineCutZValue)){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }
    else if (fUseEtaMinCut &&  photon->GetConversionRadius() >= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
	if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
	return kFALSE;
    }
    cutIndex++;

  if(TMath::Abs(photon->GetConversionZ()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
  }
    cutIndex++;


  if(TMath::Abs(photon->GetPhotonEta())> fEtaCut || TMath::Abs(photon->GetPhotonEta())< fEtaCutMin){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
  }
    cutIndex++;


  if(photon->GetPhotonPt()<fPtCut){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
  }
    cutIndex++;

  if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
 
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionCuts::SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex) {
    // Track Cuts which require AOD/ESD specific implementation

  if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  cutIndex++;

  AliAODVertex * NegVtxType=negTrack->GetProdVertex();
  AliAODVertex * PosVtxType=posTrack->GetProdVertex();
  if((NegVtxType->GetType())==AliAODVertex::kKink  || (PosVtxType->GetType())==AliAODVertex::kKink) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  return kTRUE;

}


///________________________________________________________________________
Bool_t AliConversionCuts::SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex) {
    // Track Cuts which require AOD/ESD specific implementation

   if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  cutIndex++;

  if(negTrack->GetKinkIndex(0) > 0  || posTrack->GetKinkIndex(0) > 0 ) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  return kTRUE;
}



///________________________________________________________________________
Bool_t AliConversionCuts::TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack) {
    // Track Selection for Photon Reconstruction

    Int_t cutIndex=0;
    if(hTrackCuts)hTrackCuts->Fill(cutIndex);
    cutIndex++;

  // avoid like sign
  if(negTrack->Charge() == posTrack->Charge()) {
       if(hTrackCuts)hTrackCuts->Fill(cutIndex);
        return kFALSE;
  }
  cutIndex++;

  // Number of TPC Clusters
  if( negTrack->GetNcls(1) < fMinClsTPC || posTrack->GetNcls(1) < fMinClsTPC ) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
        return kFALSE;
  }
  cutIndex++;

  // Acceptance

  if(TMath::Abs(negTrack->Eta()) > fEtaCut || TMath::Abs(negTrack->Eta()) < fEtaCutMin ||
     TMath::Abs(posTrack->Eta())> fEtaCut || TMath::Abs(posTrack->Eta())< fEtaCutMin) {
       if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  cutIndex++;

  // Single Pt Cut
  if( negTrack->Pt()< fSinglePtCut ||	posTrack->Pt()< fSinglePtCut){
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  cutIndex++;

  // AOD ESD specific cuts
  Bool_t passCuts = kTRUE;

  if(negTrack->IsA()==AliAODTrack::Class()) {
	passCuts = passCuts * SpecificTrackCuts(static_cast<AliAODTrack*>(negTrack), static_cast<AliAODTrack*>(posTrack),cutIndex);
  } else { 
	passCuts = passCuts * SpecificTrackCuts(static_cast<AliESDtrack*>(negTrack), static_cast<AliESDtrack*>(posTrack),cutIndex);
  }	

  if(!passCuts){
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
  }
  cutIndex++;

  if(hTrackCuts)hTrackCuts->Fill(cutIndex);

  return kTRUE;
		    
}

///________________________________________________________________________
Bool_t AliConversionCuts::dEdxCuts(AliVTrack *fCurrentTrack){
    // Electron Identification Cuts for Photon reconstruction

    if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
    if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error

    Int_t cutIndex=0;
    if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
    if(hTPCdEdxbefore)hTPCdEdxbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
    cutIndex++;
    

  if(fDodEdxSigmaCut == kTRUE){
      // TPC Electron Line
      if( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
		fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine){

	  if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
	  return kFALSE;
      }
      cutIndex++;

      // TPC Pion Line
	if( fCurrentTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentTrack->P()<fPIDMaxPnSigmaAbovePionLine ){
	  if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
		 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
		 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){

	      if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
	      return kFALSE;
	  }
	}
	cutIndex++;
   
	// High Pt Pion rej
	if( fCurrentTrack->P()>fPIDMaxPnSigmaAbovePionLine ){
	  if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
		 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
		 fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){

                if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
		return kFALSE;
	  }
	}
	cutIndex++;
  }
  else{cutIndex+=3;}

  if(fDoKaonRejectionLowP == kTRUE){
	if(fCurrentTrack->P()<fPIDMinPKaonRejectionLowP ){
	  if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){

	      if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
	      return kFALSE;
	  }
	}
  }
  cutIndex++;
   
  if(fDoProtonRejectionLowP == kTRUE){
	if( fCurrentTrack->P()<fPIDMinPProtonRejectionLowP ){
	  if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){

	      if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
	      return kFALSE;
	  }
	}
  }
   cutIndex++;
   
  if(fDoPionRejectionLowP == kTRUE){
	if( fCurrentTrack->P()<fPIDMinPPionRejectionLowP ){
	  if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){

	      if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
		return kFALSE;
	  }
	}
  }
  cutIndex++;
   
  if((fCurrentTrack->GetStatus() & AliESDtrack::kTOFpid) && !(fCurrentTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
     if(hTOFbefore) hTOFbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
     if(fUseTOFpid){
        if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)>fTofPIDnSigmaAboveElectronLine ||
           fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)<fTofPIDnSigmaBelowElectronLine ){
           if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
           return kFALSE;
        }
     }
     if(hTOFafter)hTOFafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
  }
     cutIndex++;
   
    // Apply TRD PID
  if(fDoTRDPID){
	if(!fPIDResponse->IdentifiedAsElectronTRD(fCurrentTrack,fPIDTRDEfficiency)){

	    if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
	    return kFALSE;
	}
  }
  cutIndex++;

  if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
  if(hTPCdEdxafter)hTPCdEdxafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::AsymmetryCut(AliConversionPhotonBase * photon,AliVEvent *event) {
    // Cut on Energy Assymetry

    for(Int_t ii=0;ii<2;ii++){

        AliVTrack *track=GetTrack(event,photon->GetTrackLabel(ii));

	if( track->P() > fMinPPhotonAsymmetryCut ){
	    Double_t trackNegAsy=0;
	    if (photon->GetPhotonP()!=0.){
		trackNegAsy= track->P()/photon->GetPhotonP();
	    }

	    if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
		return kFALSE;
	    }
	}

    }
    return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliConversionCuts::GetTrack(AliVEvent * event, Int_t label){
    //Returns pointer to the track with given ESD label
    //(Important for AOD implementation, since Track array in AOD data is different
    //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);
  if(esdEvent) {
  	if(label > event->GetNumberOfTracks() ) return NULL;
  	AliESDtrack * track = esdEvent->GetTrack(label);
  	return track;
	
  } else { 
	for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
	  AliVTrack * track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
	  
	  if(track) { 
		if(track->GetID() == label) {
		  return track;
		}
	  }
	}
  }
  
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}


///________________________________________________________________________
Bool_t AliConversionCuts::PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event){
    // Cut on Electron Probability for Photon Reconstruction

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);

  if(esdEvent){
	
	Bool_t iResult=kFALSE;
	
	Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
	Double_t *negProbArray = new Double_t[AliPID::kSPECIES];
	
	AliESDtrack* negTrack	= esdEvent->GetTrack(photon->GetTrackLabelNegative());
	AliESDtrack* posTrack	= esdEvent->GetTrack(photon->GetTrackLabelPositive());
	
	if(negProbArray && posProbArray){
	  
	  negTrack->GetTPCpid(negProbArray);
	  posTrack->GetTPCpid(posProbArray);

	  if(negProbArray[AliPID::kElectron]>=fPIDProbabilityCutNegativeParticle && posProbArray[AliPID::kElectron]>=fPIDProbabilityCutPositiveParticle){
		iResult=kTRUE;
	  }
	}
	
	delete [] posProbArray;
	delete [] negProbArray;
	return iResult;

  } else {
      ///Not possible for AODs
      return kTRUE;
  }



  
}


///________________________________________________________________________
Bool_t AliConversionCuts::AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg){
    // MC Acceptance Cuts
    //(Certain areas were excluded for photon reconstruction)

  if(particle->R()>fMaxR){
	return kFALSE;}

  if(ePos->R()>fMaxR){
	return kFALSE;
  }

  if(ePos->R()<fMinR){
	return kFALSE;
  }

  if( ePos->R() <= ((TMath::Abs(ePos->Vz())*fLineCutZRSlope)-fLineCutZValue)){
	return kFALSE;
  }
  else if (fUseEtaMinCut &&  ePos->R() >= ((TMath::Abs(ePos->Vz())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
      return kFALSE;
  }

  if(TMath::Abs(eNeg->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
	return kFALSE;
  }

  if(eNeg->Vz()!=ePos->Vz()||eNeg->R()!=ePos->R()){
	return kFALSE;
  }

  if(TMath::Abs(ePos->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
	return kFALSE;
  }

  if(TMath::Abs(particle->Eta())> fEtaCut || TMath::Abs(particle->Eta())< fEtaCutMin){
	return kFALSE;
  }

  if(TMath::Abs(ePos->Eta())> fEtaCut || TMath::Abs(ePos->Eta())< fEtaCutMin){
	return kFALSE;
  }

  if(TMath::Abs(eNeg->Eta())> fEtaCut || TMath::Abs(eNeg->Eta())< fEtaCutMin){
	return kFALSE;
  }

  if( ePos->Pt()< fSinglePtCut ||  eNeg->Pt()< fSinglePtCut){
	return kFALSE;
  }

  if(particle->Pt()<fPtCut){
	return kFALSE;
  }

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::UpdateCutString(cutIds cutID, Int_t value) {
///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
	fCutString->SetString(GetCutNumber());
  } else {
	return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
   // Initialize Cuts from a given Cut string

    AliInfo(Form("Set Cut Number: %s",analysisCutSelection.Data()));
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

  // TestFlag
  if(fCuts[0] !=9){
    AliError("Analysis Cut Selection does not start with 9");
	PrintCuts();
    return kFALSE;
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
      if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  //PrintCuts();

  // Set StandardTriggers

  if(fIsHeavyIon)SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else SelectCollisionCandidates(AliVEvent::kMB);

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID


  switch (cutID) {
  case kgoodId:
	fCuts[kgoodId] = value;
	if(value != 9) {
	  AliError("First value of cut string is wrong, aborting!!");
	  return kFALSE;
	} else {
	  return kTRUE;
	}

  case kv0FinderType:
	if( SetV0Finder(value)) {
	  fCuts[kv0FinderType] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case keProbCut:
	if( SetElectronProbCut(value)) {
	  fCuts[keProbCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kededxSigmaCut:
	if( SetTPCdEdxCutElectronLine(value)) {
	  fCuts[kededxSigmaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kpidedxSigmaCut:
	if( SetTPCdEdxCutPionLine(value)) {
	  fCuts[kpidedxSigmaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kpiMomdedxSigmaCut:
	if( SetMinMomPiondEdxCut(value)) {
	  fCuts[kpiMomdedxSigmaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kchi2GammaCut:
	if( SetChi2GammaCut(value)) {
	  fCuts[kchi2GammaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case ksinglePtCut:
	if( SetSinglePtCut(value)) {
	  fCuts[ksinglePtCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kclsTPCCut:
	if( SetTPCClusterCut(value)) {
	  fCuts[kclsTPCCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case ketaCut:
	if( SetEtaCut(value)) {
	  fCuts[ketaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kchi2MesonCut:
	if( SetChi2MesonCut(value)) {
	  fCuts[kchi2MesonCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kLowPRejectionSigmaCut:
	if( SetLowPRejectionCuts(value)) {
	  fCuts[kLowPRejectionSigmaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kQtMaxCut:
	if( SetQtMaxCut(value)) {
	  fCuts[kQtMaxCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kpiMaxMomdedxSigmaCut:
	if( SetMaxMomPiondEdxCut(value)) {
	  fCuts[kpiMaxMomdedxSigmaCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kalphaMesonCut:
	if( SetAlphaMesonCut(value)) {
	  fCuts[kalphaMesonCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kminRCut:
	if( SetRCut(value)) {
	  fCuts[kminRCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kRapidityMesonCut:
	if( SetRapidityMesonCut(value)) {
	  fCuts[kRapidityMesonCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kBackgroundScheme:
	if( SetBackgroundScheme(value)) {
	  fCuts[kBackgroundScheme] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kDegreesForRotationMethod:
	if( SetNDegreesForRotationMethod(value)) {
	  fCuts[kDegreesForRotationMethod] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kNumberOfRotations:
	if( SetNumberOfRotations(value)) {
	  fCuts[kNumberOfRotations] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kremovePileUp:
	if( SetRemovePileUp(value)) {
	  fCuts[kremovePileUp] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kselectV0AND:
	if( SetSelectV0AND(value)) {
	  fCuts[kselectV0AND] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kmultiplicityMethod:
	if( SetMultiplicityMethod(value)) {
	  fCuts[kmultiplicityMethod] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kisHeavyIon:
	if( SetIsHeavyIon(value)) {
	  fCuts[kisHeavyIon] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kCentralityMin:
	if( SetCentralityMin(value)) {
	  fCuts[kCentralityMin] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kCentralityMax:
	if( SetCentralityMax(value)) {
	  fCuts[kCentralityMax] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kTOFelectronPID:
	if( SetTOFElectronPIDCut(value)) {
	  fCuts[kTOFelectronPID] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kuseMCPSmearing:
	if( SetMCPSmearing(value)) {
	  fCuts[kuseMCPSmearing] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kdoPhotonAsymmetryCut:
	if( SetPhotonAsymmetryCut(value)) {
	  fCuts[kdoPhotonAsymmetryCut] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kPsiPair:
	if( SetPsiPairCut(value)) {
	  fCuts[kPsiPair] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kCosPAngle:
	if( SetCosPAngleCut(value)) {
	  fCuts[kCosPAngle] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;


  case kElecShare:
	if( SetSharedElectronCut(value)) {
	  fCuts[kElecShare] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;


  case kToCloseV0s:
	if( SetToCloseV0sCut(value)) {
	  fCuts[kToCloseV0s] = value;
	  UpdateCutString(cutID, value);
	  return kTRUE;
	} else return kFALSE;

  case kNCuts:
      AliError("Cut id out of range");
	return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;

}
///________________________________________________________________________
Bool_t AliConversionCuts::SetRemovePileUp(Int_t removePileUp)
{// Set Cut
    switch(removePileUp){
    case 0:
	fRemovePileUp=kFALSE;
	break;
    case 1:
	fRemovePileUp=kTRUE;
	break;
    default:
	AliError("RemovePileUpCut not defined");
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetSelectV0AND(Int_t selectV0AND)
{// Set Cut
  switch(selectV0AND){
  case 0:
      fSelectV0AND=kFALSE;
      break;
  case 1:
      fSelectV0AND=kTRUE;
      break;
  default:
      AliError("Warning: V0ANDCut not defined");
      return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetMultiplicityMethod(Int_t multiplicityMethod)
{
    // Set Cut
    fMultiplicityMethod=multiplicityMethod;

    // 0 Photon Multiplicity
    // 1 TPC Track multiplicity
    // 2 V0 Mult
    // 3 SPD Mult

    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetMCPSmearing(Int_t useMCPSmearing)
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
void AliConversionCuts::PrintCuts() {
    // Print out current Cut Selection
  for(Int_t ic = 0; ic < kNCuts; ic++) {
	printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
  }

}

///________________________________________________________________________
Bool_t AliConversionCuts::SetRCut(Int_t RCut){
    // Set Cut
    switch(RCut){
    case 0:
	fMinR=0;
	fMaxR = 180.;
	break;
    case 1:
	fMinR=2.8;
	fMaxR = 180.;
	break;
    case 2:
	fMinR=5.;
	fMaxR = 180.;
	break;
    case 3:
	fMaxR = 70.;
	fMinR = 10.;
	break;
    case 4:
	fMaxR = 70.;
	fMinR = 5.;
	break;
    case 5:
	fMaxR = 180.;
	fMinR = 10.;
	break;
    // High purity cuts for PbPb
    case 9:
	fMaxR = 180.;
	fMinR = 60.;
	break;

    default:
	AliError("RCut not defined");
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut)
{   // Set Cut
    switch(ededxSigmaCut){
    case 0: // -10,10
	fPIDnSigmaBelowElectronLine=-10;
	fPIDnSigmaAboveElectronLine=10;
	break;
    case 1: // -5,5
	fPIDnSigmaBelowElectronLine=-5;
	fPIDnSigmaAboveElectronLine=5;
	break;
    case 2: // -3,5
	fPIDnSigmaBelowElectronLine=-3;
	fPIDnSigmaAboveElectronLine=5;
	break;
    case 3: // -4,5
	fPIDnSigmaBelowElectronLine=-4;
	fPIDnSigmaAboveElectronLine=5;
	break;
    case 4: // -6,7
       fPIDnSigmaBelowElectronLine=-6;
       fPIDnSigmaAboveElectronLine=7;
       break;
    case 5: // -4,4
       fPIDnSigmaBelowElectronLine=-4;
       fPIDnSigmaAboveElectronLine=4;
       break;
    case 6: // -2.5,4
       fPIDnSigmaBelowElectronLine=-2.5;
       fPIDnSigmaAboveElectronLine=4;
       break;
    case 7: // -2,3.5
       fPIDnSigmaBelowElectronLine=-2;
       fPIDnSigmaAboveElectronLine=3.5;
       break;
    default:
	AliError("TPCdEdxCutElectronLine not defined");
	return kFALSE;
        
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut)
{   // Set Cut

    switch(pidedxSigmaCut){
    case 0:  // -10
	fPIDnSigmaAbovePionLine=-10;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 1:   // 0
	fPIDnSigmaAbovePionLine=0;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 2:  // 1
	fPIDnSigmaAbovePionLine=1;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 3:  // 1
	fPIDnSigmaAbovePionLine=-1;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 4:  // 1
	fPIDnSigmaAbovePionLine=2.5;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 5:  // 1
	fPIDnSigmaAbovePionLine=2.;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 6:  // 1
	fPIDnSigmaAbovePionLine=2.;
	fPIDnSigmaAbovePionLineHighPt=0.5;
	break;
    case 7:  // 1
	fPIDnSigmaAbovePionLine=3.5;
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    case 8:  // 1
	fPIDnSigmaAbovePionLine=2.;
	fPIDnSigmaAbovePionLineHighPt=1.;
	break;
    case 9:
	fPIDnSigmaAbovePionLine=3.0; // We need a bit less tight cut on dE/dx
	fPIDnSigmaAbovePionLineHighPt=-10;
	break;
    default:
	AliError(Form("Warning: pidedxSigmaCut not defined %d",pidedxSigmaCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetMinMomPiondEdxCut(Int_t piMomdedxSigmaCut)
{   // Set Cut
    switch(piMomdedxSigmaCut){
    case 0:  // 0.5 GeV
	fPIDMinPnSigmaAbovePionLine=0.5;
	break;
    case 1:  // 1. GeV
	fPIDMinPnSigmaAbovePionLine=1.;
	break;
    case 2:  // 1.5 GeV
	fPIDMinPnSigmaAbovePionLine=1.5;
	break;
    case 3:  // 20.0 GeV
	fPIDMinPnSigmaAbovePionLine=20.;
	break;
    case 4:  // 50.0 GeV
	fPIDMinPnSigmaAbovePionLine=50.;
	break;
    case 5:  // 0.3 GeV
	fPIDMinPnSigmaAbovePionLine=0.3;
	break;
    case 6:  // 0.25 GeV
	fPIDMinPnSigmaAbovePionLine=0.25;
	break;
    case 7:  // 0.4 GeV
	fPIDMinPnSigmaAbovePionLine=0.4;
	break;
    default:
	AliError(Form("piMomdedxSigmaCut not defined %d",piMomdedxSigmaCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetChi2GammaCut(Int_t chi2GammaCut)
{   // Set Cut

    switch(chi2GammaCut){
    case 0: // 100
	fChi2CutConversion = 100.;
	break;
    case 1:  // 50
	fChi2CutConversion = 50.;
	break;
    case 2:  // 30
	fChi2CutConversion = 30.;
	break;
    case 3:
	fChi2CutConversion = 200.;
	break;
    case 4:
	fChi2CutConversion = 500.;
	break;
    case 5:
	fChi2CutConversion = 100000.;
	break;
    case 6:
	fChi2CutConversion = 5.;
	break;
    case 7:
	fChi2CutConversion = 10.;
	break;
    case 8:
	fChi2CutConversion = 20.;
	break;
    case 9:
	fChi2CutConversion = 15.;
	break;
    default:
	AliError(Form("Warning: Chi2GammaCut not defined %d",chi2GammaCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetV0Finder(Int_t v0FinderType)
{   // Set Cut
    switch (v0FinderType){
    case 0:  // on fly V0 finder
	fUseOnFlyV0Finder=kTRUE;
	break;
    case 1:  // offline V0 finder
	fUseOnFlyV0Finder=kFALSE;
	break;
    default:
        AliError(Form(" v0FinderType not defined %d",v0FinderType));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetElectronProbCut(Int_t eProbCut)
{   // Set Cut

    switch(eProbCut){
    case 0:
	fPIDProbabilityCutNegativeParticle=0;
	fPIDProbabilityCutPositiveParticle=0;
	break;
    case 1:
	fPIDProbabilityCutNegativeParticle=0.1;
	fPIDProbabilityCutPositiveParticle=0.1;
	break;
    case 2:
	fPIDProbabilityCutNegativeParticle=0.5;
	fPIDProbabilityCutPositiveParticle=0.5;
	break;
    case 3:
	fPIDProbabilityCutNegativeParticle=0.7;
	fPIDProbabilityCutPositiveParticle=0.7;
	break;
    default:
	AliError(Form("Warning: eProbCut not defined %d",eProbCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetSinglePtCut(Int_t singlePtCut)
{   // Set Cut
    switch(singlePtCut){
    case 0: // 0.050 GeV
	fSinglePtCut = 0.050;
	break;
    case 1:  // 0.100 GeV
	fSinglePtCut = 0.100;
	break;
    case 2:  // 0.150 GeV
	fSinglePtCut = 0.150;
	break;
    case 3:  // 0.200 GeV
	fSinglePtCut = 0.200;
	break;
    case 4:  // 0.075 GeV
	fSinglePtCut = 0.075;
	break;
    case 5:  // 0.125 GeV
	fSinglePtCut = 0.125;
	break;
    default:
	AliError(Form("singlePtCut not defined %d",singlePtCut));
	return kFALSE;
    }
    return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetTPCClusterCut(Int_t clsTPCCut)
{   // Set Cut
    switch(clsTPCCut){
    case 0: // 0
	fMinClsTPC= 0.;
	break;
    case 1:  // 70
	fMinClsTPC= 70.;
	break;
    case 2:  // 80
	fMinClsTPC= 80.;
	break;
    case 3:  // 100
	fMinClsTPC= 100.;
	break;
    case 4:  // 60% of findable clusters
	fMinClsTPCToF= 0.6;
	fUseCorrectedTPCClsInfo=0;
	break;
    case 5:  // 0% of findable clusters
	fMinClsTPCToF= 0.0;
	fUseCorrectedTPCClsInfo=1;
	break;
    case 6:  // 0% of findable clusters
	fMinClsTPCToF= 0.7;
	fUseCorrectedTPCClsInfo=1;
	break;
    case 7:  // 0% of findable clusters
	fMinClsTPCToF= 0.35;
	fUseCorrectedTPCClsInfo=0;
	break;
    case 8:
	fMinClsTPCToF= 0.35;
	fUseCorrectedTPCClsInfo=1;
	break;
    case 9:
	fMinClsTPCToF= 0.6;
	fUseCorrectedTPCClsInfo=1;
	break;
    default:
	AliError(Form("Warning: clsTPCCut not defined %d",clsTPCCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetEtaCut(Int_t etaCut)
{   // Set Cut

   //Set Standard LineCutZValues
   fLineCutZValueMin = -2;
   fLineCutZValue = 7.;
   
    switch(etaCut){
    case 0: // 0.9
	fEtaCut		= 0.9;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= -0.1;
	fLineCutZRSlopeMin = 0.;
	break;
    case 1:	// 1.2
	fEtaCut		= 1.2;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= -0.1;
	fLineCutZRSlopeMin = 0.;
	break;
    case 2:	// 1.4
	fEtaCut		= 1.4;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= -0.1;
	fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCut)));
	break;
    case 3: // 0.8
	fEtaCut		= 0.8;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= -0.1;
	fLineCutZRSlopeMin = 0.;
	break;
    case 4: // 0.75
	fEtaCut		= 0.75;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= -0.1;
	fLineCutZRSlopeMin = 0.;
	break;
    case 5: // 0.9 - 1.4
	fEtaCut		= 1.4;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= 0.9;
	fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
	break;
    case 6: // 0.9 - 1.2
	fEtaCut		= 1.2;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= 0.9;
	fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
	break;
    case 7: // 0.1 - 0.8
	fEtaCut		= 0.8;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= 0.1;
	fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
	break;
    case 8: // 0.1 - 0.8
	fEtaCut		= 0.9;
	fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
	fEtaCutMin		= 0.1;
	fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
        break;
    default:
	AliError(Form(" EtaCut not defined %d",etaCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetChi2MesonCut(Int_t chi2MesonCut)
{   // Set Cut
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
	AliError(Form("Chi2MesonCut not defined %d",chi2MesonCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut)
{   // Set Cut
    switch(piMaxMomdedxSigmaCut){
    case 0:  // 100. GeV
	fPIDMaxPnSigmaAbovePionLine=100.;
	break;
    case 1:  // 5. GeV
	fPIDMaxPnSigmaAbovePionLine=5.;
	break;
    case 2:  // 4. GeV
	fPIDMaxPnSigmaAbovePionLine=4.;
	break;
    case 3:  // 3.5 GeV
	fPIDMaxPnSigmaAbovePionLine=3.5;
	break;
    case 4:  // 3. GeV
	fPIDMaxPnSigmaAbovePionLine=3.;
	break;
    default:
	AliError(Form("piMaxMomdedxSigmaCut not defined %d",piMaxMomdedxSigmaCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetIsHeavyIon(Int_t isHeavyIon)
{   // Set Cut
   switch(isHeavyIon){
		case 0:
			fIsHeavyIon=0;
			break;
		case 1:
			fIsHeavyIon=1;
			fDetectorCentrality=0;
			break;
		case 2:
			fIsHeavyIon=1;
			fDetectorCentrality=1;
			break;
		case 3: //allows to select centrality 0-45% in steps of 5% for V0 Multiplicity
			fIsHeavyIon=1;
			fDetectorCentrality=0;	  
			fModCentralityClass=1;
			break;
		case 4: //allows to select centrality 45-90% in steps of 5% for V0 Multiplicity
			fIsHeavyIon=1;
			fDetectorCentrality=0;
			fModCentralityClass=2;
			break;
		default:
			AliError(Form("SetHeavyIon not defined %d",isHeavyIon));
			return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetAlphaMesonCut(Int_t alphaMesonCut)
{   // Set Cut
    switch(alphaMesonCut){
    case 0:	// 0- 0.7
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.7;
	break;
    case 1:	// 0-0.5
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.5;
	break;
    case 2:	// 0.5-1
	fAlphaMinCutMeson	 = 0.5;
	fAlphaCutMeson	 = 1.;
	break;
    case 3:	// 0.0-1
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 1.;
	break;
    case 4:	// 0-0.65
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.65;
	break;
    case 5:	// 0-0.75
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.75;
	break;
    case 6:	// 0-0.8
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.8;
	break;
    case 7:	// 0.0-0.85
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.85;
	break;
    case 8:	// 0.0-0.6
	fAlphaMinCutMeson	 = 0.0;
	fAlphaCutMeson	 = 0.6;
	break;
    default:
        AliError(Form("AlphaMesonCut not defined %d",alphaMesonCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetRapidityMesonCut(Int_t RapidityMesonCut)
{   // Set Cut
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

    default:
        AliError(Form("RapidityMesonCut not defined %d",RapidityMesonCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut)
{   // Set Cut
    switch(LowPRejectionSigmaCut){
	case 0:  //
	fPIDnSigmaAtLowPAroundKaonLine=0;
	fPIDnSigmaAtLowPAroundProtonLine=0;
	fPIDnSigmaAtLowPAroundPionLine=0;
	break;
    case 1:  //
	fPIDnSigmaAtLowPAroundKaonLine=0.5;
	fPIDnSigmaAtLowPAroundProtonLine=0.5;
	fPIDnSigmaAtLowPAroundPionLine=0.5;
	break;
    case 2:  //
	fPIDnSigmaAtLowPAroundKaonLine=1;
	fPIDnSigmaAtLowPAroundProtonLine=1;
	fPIDnSigmaAtLowPAroundPionLine=1;
	break;
    case 3:  //
	fPIDnSigmaAtLowPAroundKaonLine=2.;
	fPIDnSigmaAtLowPAroundProtonLine=2.;
	fPIDnSigmaAtLowPAroundPionLine=2.;
	break;
    case 4:  //
	fPIDnSigmaAtLowPAroundKaonLine=0.;
	fPIDnSigmaAtLowPAroundProtonLine=0.;
	fPIDnSigmaAtLowPAroundPionLine=1;
	break;
    case 5:  //
	fPIDnSigmaAtLowPAroundKaonLine=0.;
	fPIDnSigmaAtLowPAroundProtonLine=0.;
	fPIDnSigmaAtLowPAroundPionLine=1.5;
	break;
    case 6:  //
	fPIDnSigmaAtLowPAroundKaonLine=0.;
	fPIDnSigmaAtLowPAroundProtonLine=0.;
	fPIDnSigmaAtLowPAroundPionLine=2.;
	break;
    default:
        AliError(Form("LowPRejectionSigmaCut not defined %d",LowPRejectionSigmaCut));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetTOFElectronPIDCut(Int_t TOFelectronPID){
    // Set Cut
    switch(TOFelectronPID){ // RRnewTOF start //////////////////////////////////////////////////////////////////////////
    case 0: // no cut
	fUseTOFpid = kFALSE;
	fTofPIDnSigmaBelowElectronLine=-100;
	fTofPIDnSigmaAboveElectronLine=100;
	break;
    case 1: // -7,7
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-7;
	fTofPIDnSigmaAboveElectronLine=7;
	break;
    case 2: // -5,5
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-5;
	fTofPIDnSigmaAboveElectronLine=5;
	break;
    case 3: // -3,5
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-3;
	fTofPIDnSigmaAboveElectronLine=5;
	break;
    case 4: // -2,3
	fUseTOFpid = kTRUE;
	fTofPIDnSigmaBelowElectronLine=-2;
	fTofPIDnSigmaAboveElectronLine=3;
	break;
    default:
        AliError(Form("TOFElectronCut not defined %d",TOFelectronPID));
	return kFALSE;
    } //////////////////////// RRnewTOF end //////////////////////////////////////////////////////////////////////////
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetTRDElectronCut(Int_t TRDElectronCut)
{   // Set Cut
    switch(TRDElectronCut){
    case 0:
	fDoTRDPID=kFALSE;
	break;
    case 1:
	fDoTRDPID=kTRUE;
	fPIDTRDEfficiency=0.1;
	break;
    case 8:
	fDoTRDPID=kTRUE;
	fPIDTRDEfficiency=0.8;
	break;
    case 9:
	fDoTRDPID=kTRUE;
	fPIDTRDEfficiency=0.9;
	break;
    default:
	AliError(Form("TRDElectronCut not defined %d",TRDElectronCut));
	return kFALSE;
    }

    return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetQtMaxCut(Int_t QtMaxCut)
{   // Set Cut
    switch(QtMaxCut){
    case 0: //
	fQtMax=1.;
	fDoQtGammaSelection=kFALSE;      //No Qt selection (true by default)
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break;
    case 1:
	fQtMax=0.1;
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break;
    case 2:
	fQtMax=0.07;
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break;
    case 3:
	fQtMax=0.05;
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break;
    case 4:
	fQtMax=0.03;
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break;
    case 5: // RR try to improve (get rid of) low InvMass peak in PbPb
	fQtMax=0.02;
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break; // end RR ///////////////////////////////////////////////
    case 6:  // RRnew start: pT dependent qT cut
	fQtMax=0.02;
	fDoHighPtQtGammaSelection=kTRUE;
	fHighPtQtMax=0.06;
	fPtBorderForQt=2.5;
	break; // RRnew end ////////////////////////////////////////////
    case 7:
	fQtMax=0.15;
	fDoHighPtQtGammaSelection=kFALSE; // RRnew
	fHighPtQtMax=100.;	        // RRnew
	fPtBorderForQt=100.;	        // RRnew
	break;
    default:
	AliError(Form("Warning: QtMaxCut not defined %d",QtMaxCut));
	return kFALSE;
    }
    return kTRUE;
}

//-------------------------------------------------------------
Double_t AliConversionCuts::GetCentrality(AliVEvent *event)
{   // Get Event Centrality

    AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
    if(esdEvent){
	AliCentrality *fESDCentrality=(AliCentrality*)esdEvent->GetCentrality();

	if(fDetectorCentrality==0){
	    return fESDCentrality->GetCentralityPercentile("V0M"); // default
	}
	if(fDetectorCentrality==1){
	    return fESDCentrality->GetCentralityPercentile("CL1");
	}
    }

    AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);
    if(aodEvent){
	if(aodEvent->GetHeader()){return aodEvent->GetHeader()->GetCentrality();}
    }

    return -1;
}

//-------------------------------------------------------------
Bool_t AliConversionCuts::IsCentralitySelected(AliVEvent *event)
{   // Centrality Selection
   if(!fIsHeavyIon)return kTRUE;
    
	if(fCentralityMin == 0 && fCentralityMax == 0) return kTRUE;//0-100%
	if(fCentralityMin >= fCentralityMax) return kTRUE;//0-100%
   
	Double_t centrality=GetCentrality(event);
	if(centrality<0)return kFALSE;
	
	Int_t centralityC=0;
	if (fModCentralityClass == 0){
		centralityC= Int_t(centrality/10);
		if(centralityC >= fCentralityMin && centralityC < fCentralityMax)
			return kTRUE;
		else 
			return kFALSE;
	}	else if (fModCentralityClass ==1){
		centralityC= Int_t(centrality);
		if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
			return kTRUE;
		} else { 
			return kFALSE;
		}
	} else if (fModCentralityClass ==2){
		centralityC= Int_t(centrality+1);
		if(centralityC >= (fCentralityMin*5+45) && centralityC < (fCentralityMax*5+45))
			return kTRUE;
		else 
			return kFALSE;
	}
	return kFALSE;
}

//-------------------------------------------------------------
Bool_t AliConversionCuts::SetCentralityMin(Int_t minCentrality)
{
    // Set Cut
    if(minCentrality<0||minCentrality>9){
	AliError(Form("minCentrality not defined %d",minCentrality));
	return kFALSE;
    }

    fCentralityMin=minCentrality;
    return kTRUE;
}
//-------------------------------------------------------------
Bool_t AliConversionCuts::SetCentralityMax(Int_t maxCentrality)
{
    // Set Cut
    if(maxCentrality<0||maxCentrality>9){
	AliError(Form("maxCentrality not defined %d",maxCentrality));
	return kFALSE;
    }

    fCentralityMax=maxCentrality;
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut){
    // Set Cut
    switch(doPhotonAsymmetryCut){
    case 0:
	fDoPhotonAsymmetryCut=0;
	fMinPPhotonAsymmetryCut=100.;
	fMinPhotonAsymmetry=0.;
	break;
    case 1:
	fDoPhotonAsymmetryCut=1;
	fMinPPhotonAsymmetryCut=3.5;
	fMinPhotonAsymmetry=0.04;
	break;
    case 2:
	fDoPhotonAsymmetryCut=1;
	fMinPPhotonAsymmetryCut=3.5;
	fMinPhotonAsymmetry=0.06;
	break;
    default:
	AliError(Form("PhotonAsymmetryCut not defined %d",doPhotonAsymmetryCut));
	return kFALSE;
    }
    fCuts[kdoPhotonAsymmetryCut]=doPhotonAsymmetryCut;
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetBackgroundScheme(Int_t BackgroundScheme){
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
    default:
	AliError(Form("BackgroundScheme not defined %d",BackgroundScheme));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetNDegreesForRotationMethod(Int_t DegreesForRotationMethod){
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
	AliError(Form("DegreesForRotationMethod not defined %d",DegreesForRotationMethod));
	return kFALSE;
    }
    fCuts[kDegreesForRotationMethod]=DegreesForRotationMethod;
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetNumberOfRotations(Int_t NumberOfRotations)
{   // Set Cut
    switch(NumberOfRotations){
    case 0:
	fnumberOfRotationEventsForBG = 5;
	break;
    case 1:
	fnumberOfRotationEventsForBG = 10;
	break;
    case 2:
	fnumberOfRotationEventsForBG = 15;
	break;
    case 3:
	fnumberOfRotationEventsForBG = 20;
	break;
    case 4:
	fnumberOfRotationEventsForBG = 2;
	break;
    case 5:
	fnumberOfRotationEventsForBG = 50;
	break;
    case 6:
	fnumberOfRotationEventsForBG = 80;
	break;
    case 7:
	fnumberOfRotationEventsForBG = 100;
	break;
    default:
	AliError(Form("NumberOfRotations not defined %d",NumberOfRotations));
	return kFALSE;
    }
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetPsiPairCut(Int_t psiCut) {
  

  switch(psiCut) {
  case 0:
	fPsiPairCut = 10000; // 
	break;
  case 1:
	fPsiPairCut = 0.1; // 
	break;
  case 2:
	fPsiPairCut = 0.05; // Standard
	break;
  case 3:
	fPsiPairCut = 0.035; // 
	break;
  case 4:
	fPsiPairCut = 0.15; // 
	break;
  case 5:
	fPsiPairCut = 0.2; // 
	break;
  case 6:
	fPsiPairCut = 0.03; // 
	break;
  case 7:
	fPsiPairCut = 0.025; // 
	break;
  case 8:
	fPsiPairCut = 0.01; // 
	break;
  default:
      AliError(Form("PsiPairCut not defined %d",psiCut));
      return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetCosPAngleCut(Int_t cosCut) {

    switch(cosCut){
    case 0:
	fCosPAngleCut = TMath::Pi(); //
	break;
    case 1:
	fCosPAngleCut = 0.1; //
	break;
    case 2:
	fCosPAngleCut = 0.05; //
	break;
    case 3:
	fCosPAngleCut = 0.025; // Standard
	break;
    case 4:
	fCosPAngleCut = 0.01; //
	break;
    default:
	AliError(Form("Cosine Pointing Angle cut not defined %d",cosCut));
	return kFALSE;
    }
	
    return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionCuts::VertexZCut(AliVEvent *event){
    // Cut on z position of primary vertex
     Double_t fVertexZ=event->GetPrimaryVertex()->GetZ();

     if(TMath::Abs(fVertexZ)>fMaxVertexZ)return kFALSE;
     return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetSharedElectronCut(Int_t sharedElec) {

   switch(sharedElec){
    case 0:
       fDoSharedElecCut = kFALSE;
	break;
    case 1:
       fDoSharedElecCut = kTRUE;
	break;
    default:
	AliError(Form("Shared Electron Cut not defined %d",sharedElec));
	return kFALSE;
    }
	
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetToCloseV0sCut(Int_t toClose) {

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
       AliError(Form("Shared Electron Cut not defined %d",toClose));
	return kFALSE;
   }
   return kTRUE;
}



///________________________________________________________________________

Int_t AliConversionCuts::GetNumberOfContributorsVtx(AliVEvent *event){
    // returns number of contributors to the vertex
    
    AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(event);
    if(fESDEvent){
	if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
	    return fESDEvent->GetPrimaryVertexTracks()->GetNContributors();
	}
     
	if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()<1) {
	    //		return 0;
	    //-AM test pi0s without SPD only vertex
	    if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
		return fESDEvent->GetPrimaryVertexSPD()->GetNContributors();

	    }
	    if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
		return 0;
	    }
	}
    }

    AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(event);
    if(fAODEvent){
	if(fAODEvent->GetPrimaryVertex()->GetNContributors()>0) {
	    return fAODEvent->GetPrimaryVertex()->GetNContributors();
	}
	if(fAODEvent->GetPrimaryVertex()->GetNContributors()<1) {
	    if(fAODEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
		return fAODEvent->GetPrimaryVertexSPD()->GetNContributors();
	    }
	    if(fAODEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
		AliWarning(Form("Number of contributors from bad vertex type:: %s",fAODEvent->GetPrimaryVertex()->GetName()));
		return 0;
	    }
	}
    }


  return 0;
}

///________________________________________________________________________

Bool_t AliConversionCuts::IsTriggerSelected()
{
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    UInt_t isSelected = AliVEvent::kAny;
    if( fInputHandler && fInputHandler->GetEventSelection()) {
	// Get the actual offline trigger mask for the event and AND it with the
	// requested mask. If no mask requested select by default the event.
	if (fOfflineTriggerMask)
	    isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
    }
   
    if(!isSelected)return kFALSE;

    // Fill Histogram
    if(hTriggerClass){
	if (fInputHandler->IsEventSelected() & AliVEvent::kAny)hTriggerClass->Fill(0);
	if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClass->Fill(1);
	if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)hTriggerClass->Fill(2);
	if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)hTriggerClass->Fill(3);
    }

    
    return kTRUE;

}

///________________________________________________________________________
Int_t AliConversionCuts::GetFirstTPCRow(Double_t radius){
    // Get first TPC row
	Int_t firstTPCRow=0;
	Double_t radiusI	=	84.8;
	Double_t radiusO	= 134.6;
	Double_t radiusOB = 198.;
	Double_t rSizeI	 = 0.75;
	Double_t rSizeO	 = 1.;
	Double_t rSizeOB	= 1.5;
	Int_t nClsI=63;
	Int_t nClsIO=127;

	if(radius <= radiusI){
		return firstTPCRow;
	}
	if(radius>radiusI && radius<=radiusO){
		firstTPCRow = (Int_t)((radius-radiusI)/rSizeI);
	}
	if(radius>radiusO && radius<=radiusOB){
		firstTPCRow = (Int_t)(nClsI+(radius-radiusO)/rSizeO);
	}

	if(radius>radiusOB){
		firstTPCRow =(Int_t)(nClsIO+(radius-radiusOB)/rSizeOB);
	}

	return firstTPCRow;
}

Bool_t AliConversionCuts::CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const {
  ///Check if passes cosine of pointing angle cut
   if(GetCosineOfPointingAngle(photon, event) < (TMath::Cos(fCosPAngleCut))){
      return kFALSE;
   }
   return kTRUE;
}

Double_t AliConversionCuts::GetCosineOfPointingAngle( const AliConversionPhotonBase * photon, AliVEvent * event) const{
   // calculates the pointing angle of the recalculated V0 

   Double_t momV0[3] = {0,0,0};
   if(event->IsA()==AliESDEvent::Class()){
      AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
      AliESDv0 *v0 = esdEvent->GetV0(photon->GetV0Index());
      v0->GetPxPyPz(momV0[0],momV0[1],momV0[2]);
   }
   if(event->IsA()==AliAODEvent::Class()){
      momV0[0] = photon->GetPx();
      momV0[1] = photon->GetPy();
      momV0[2] = photon->GetPz();
   }
   
   //Double_t momV0[3] = { photon->GetPx(), photon->GetPy(), photon->GetPz() }; //momentum of the V0
   Double_t PosV0[3] = { photon->GetConversionX() - event->GetPrimaryVertex()->GetX(), 
                         photon->GetConversionY() - event->GetPrimaryVertex()->GetY(), 
                         photon->GetConversionZ() - event->GetPrimaryVertex()->GetZ() }; //Recalculated V0 Position vector
   
   Double_t momV02 = momV0[0]*momV0[0] + momV0[1]*momV0[1] + momV0[2]*momV0[2];
   Double_t PosV02 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1] + PosV0[2]*PosV0[2];

   Double_t cosinePointingAngle = (PosV0[0]*momV0[0] +  PosV0[1]*momV0[1] + PosV0[2]*momV0[2] ) / TMath::Sqrt(momV02 * PosV02);
  
   return cosinePointingAngle;
}


Bool_t AliConversionCuts::PsiPairCut(const AliConversionPhotonBase * photon) const {

    if(photon->GetPsiPair() > fPsiPairCut){
	return kFALSE;}
    else{return kTRUE;}

}

///________________________________________________________________________
TString AliConversionCuts::GetCutNumber(){
    // returns TString with current cut number
  TString a(kNCuts);
  for(Int_t ii=0;ii<kNCuts;ii++){
      a.Append(Form("%d",fCuts[ii]));
  }
  return a;
}

///________________________________________________________________________
void AliConversionCuts::SmearParticle(AliAODConversionPhoton* photon)
{
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
///________________________________________________________________________
void AliConversionCuts::FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0){
   
   Int_t posLabel = photon->GetTrackLabelPositive();
   Int_t negLabel = photon->GetTrackLabelNegative();
   
   fElectronLabelArray[nV0*2] = posLabel;
   fElectronLabelArray[(nV0*2)+1] = negLabel;
}
///________________________________________________________________________
Bool_t AliConversionCuts::RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s){

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
Bool_t AliConversionCuts::RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0){


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
