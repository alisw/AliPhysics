
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *				       					 									*
 * Authors: Friederike Bock												  	*
 * Version 1.0								  								*
 *																			*
 * Permission to use, copy, modify and distribute this software and its	  	*
 * documentation strictly for non-commercial purposes is hereby granted	  	*
 * without fee, provided that the above copyright notice appears in all	  	*
 * copies and that both the copyright notice and this permission notice	  	*
 * appear in the supporting documentation. The authors make no claims	  	*
 * about the suitability of this software for any purpose. It is		  	*
 * provided "as is" without express or implied warranty.		  			*
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////


#include "AliPrimaryPionCuts.h"
#include "AliAODConversionPhoton.h"
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
#include "AliMCEvent.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "TList.h"
class iostream;

using namespace std;

ClassImp(AliPrimaryPionCuts)


const char* AliPrimaryPionCuts::fgkCutNames[AliPrimaryPionCuts::kNCuts] = {
	"kEtaCut",				// 0
	"kClsITSCut",			// 1
	"kClsTPCCut",			// 2
	"kDCAcut",				// 3
	"kPtCut",				// 4
	"kPiDedxSigmaITSCut",	// 5
	"kPiDedxSigmaTPCCut",	// 6
	"kPiTOFSigmaCut",		// 7 
	"kMassCut"				// 8
};

//________________________________________________________________________
AliPrimaryPionCuts::AliPrimaryPionCuts(const char *name,const char *title) : AliAnalysisCuts(name,title),
	fHistograms(NULL),
    fDoLightOutput(kFALSE),
	fPIDResponse(NULL),
	fEsdTrackCuts(NULL),
	fEsdTrackCutsGC(NULL),
	fEtaCut(0.9),
	fEtaShift(0.0),
	fDoEtaCut(kFALSE),
	fPtCut(0.0),
	fMinClsTPC(0), // minimum clusters in the TPC
  fChi2PerClsTPC(9999), // maximum Chi2 per cluster in the TPC
  fRequireTPCRefit(kFALSE), // require a refit in the TPC
	fMinClsTPCToF(0), // minimum clusters to findable clusters
    fMaxSharedClsTPCFrac(99), // maximum fraction of shared clusters to TPCnClus
	fMinClsITS(0), // minimum clusters to findable clusters
	fDodEdxSigmaITSCut(kFALSE),
	fDodEdxSigmaTPCCut(kTRUE),
	fDoTOFsigmaCut(kFALSE), // RRnewTOF
	fPIDnSigmaAbovePionLineITS(100),
	fPIDnSigmaBelowPionLineITS(-100),
	fPIDnSigmaAbovePionLineTPC(100),
	fPIDnSigmaBelowPionLineTPC(-100),
	fPIDnSigmaAbovePionLineTOF(100), // RRnewTOF
	fPIDnSigmaBelowPionLineTOF(-100), // RRnewTOF
	fUseCorrectedTPCClsInfo(kFALSE),
	fUseTOFpid(kFALSE),
	fRequireTOF(kFALSE),
	fDoMassCut(kFALSE),
    fDoMassCut_WithNDM(kFALSE),
	fMassCut(10),
    fMassCut_WithNDM(10),
	fUse4VecForMass(kFALSE),
	fRequireVertexConstrain(kFALSE),
	fDoWeights(kFALSE),
  fMaxDCAToVertexZ(8000),
  fMaxDCAToVertexXY(8000),
  fUsePtDepXYDCA(kFALSE),
  fUseDCAToVertex2D(kFALSE),
  fMaxDCAToVertexXYPtDep(""),
  fRunFlag(1500),
	fCutString(NULL),
  fCutStringRead(""),
	fHistCutIndex(NULL),
	fHistdEdxCuts(NULL),
	fHistITSdEdxbefore(NULL),
	fHistITSdEdxafter(NULL),
	fHistTPCdEdxbefore(NULL),
	fHistTPCdEdxafter(NULL),
	fHistTPCdEdxSignalbefore(NULL),
	fHistTPCdEdxSignalafter(NULL),
	fHistTOFbefore(NULL),
	fHistTOFafter(NULL),
	fHistTrackDCAxyPtbefore(NULL),
	fHistTrackDCAxyPtafter(NULL),
	fHistTrackDCAzPtbefore(NULL),
	fHistTrackDCAzPtafter(NULL),
	fHistTrackNFindClsPtTPCbefore(NULL),
	fHistTrackNFindClsPtTPCafter(NULL),
	fHistTrackSelectedEta(NULL),
	fHistTrackSelectedPhi(NULL),
	fHistTrackSelectedPt(NULL),
	fHistTrackSelectedPtWithoutITS(NULL),
	fStringITSClusterCut(""),
	fPeriodName("")
{
	InitPIDResponse();
	for(Int_t jj=0;jj<kNCuts;jj++){ fCuts[jj]=0; }
	fCutString=new TObjString((GetCutNumber()).Data());

	// Using standard function for setting Cuts
	if (fEsdTrackCuts==NULL) fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
}

//________________________________________________________________________
AliPrimaryPionCuts::AliPrimaryPionCuts(const AliPrimaryPionCuts &ref) : AliAnalysisCuts(ref),
	fHistograms(NULL),
    fDoLightOutput(ref.fDoLightOutput),
	fPIDResponse(NULL),
	fEsdTrackCuts(ref.fEsdTrackCuts),
	fEsdTrackCutsGC(ref.fEsdTrackCutsGC),
	fEtaCut(ref.fEtaCut),
	fEtaShift(ref.fEtaShift),
	fDoEtaCut(ref.fDoEtaCut),
	fPtCut(ref.fPtCut),
	fMinClsTPC(ref.fMinClsTPC), // minimum clusters in the TPC
    fChi2PerClsTPC(ref.fChi2PerClsTPC), // maximum Chi2 per cluster in the TPC
    fRequireTPCRefit(ref.fRequireTPCRefit), // require a refit in the TPC
	fMinClsTPCToF(ref.fMinClsTPCToF), // minimum clusters to findable clusters
    fMaxSharedClsTPCFrac(ref.fMaxSharedClsTPCFrac), // maximum fraction of shared clusters to TPCnClus
	fMinClsITS(ref.fMinClsITS), // minimum clusters to findable clusters
	fDodEdxSigmaITSCut(ref.fDodEdxSigmaITSCut),
	fDodEdxSigmaTPCCut(ref.fDodEdxSigmaTPCCut),
	fDoTOFsigmaCut(ref.fDoTOFsigmaCut), // RRnewTOF
	fPIDnSigmaAbovePionLineITS(ref.fPIDnSigmaAbovePionLineITS),
	fPIDnSigmaBelowPionLineITS(ref.fPIDnSigmaBelowPionLineITS),
	fPIDnSigmaAbovePionLineTPC(ref.fPIDnSigmaAbovePionLineTPC),
	fPIDnSigmaBelowPionLineTPC(ref.fPIDnSigmaBelowPionLineTPC),
	fPIDnSigmaAbovePionLineTOF(ref.fPIDnSigmaAbovePionLineTOF),
	fPIDnSigmaBelowPionLineTOF(ref.fPIDnSigmaBelowPionLineTOF),
	fUseCorrectedTPCClsInfo(ref.fUseCorrectedTPCClsInfo),
	fUseTOFpid(ref.fUseTOFpid),
	fRequireTOF(ref.fRequireTOF),
	fDoMassCut(ref.fDoMassCut),
    fDoMassCut_WithNDM(ref.fDoMassCut_WithNDM),
	fMassCut(ref.fMassCut),
    fMassCut_WithNDM(ref.fMassCut_WithNDM),
	fUse4VecForMass(ref.fUse4VecForMass),
	fRequireVertexConstrain(ref.fRequireVertexConstrain),
	fDoWeights(ref.fDoWeights),
    fMaxDCAToVertexZ(ref.fMaxDCAToVertexZ),
	fMaxDCAToVertexXY(ref.fMaxDCAToVertexXY),
	fUsePtDepXYDCA(ref.fUsePtDepXYDCA),
	fUseDCAToVertex2D(ref.fUseDCAToVertex2D),
	fMaxDCAToVertexXYPtDep(ref.fMaxDCAToVertexXYPtDep),
	fRunFlag(ref.fRunFlag),
	fCutString(NULL),
	fCutStringRead(""),
	fHistCutIndex(NULL),
	fHistdEdxCuts(NULL),
	fHistITSdEdxbefore(NULL),
	fHistITSdEdxafter(NULL),
	fHistTPCdEdxbefore(NULL),
	fHistTPCdEdxafter(NULL),
	fHistTPCdEdxSignalbefore(NULL),
	fHistTPCdEdxSignalafter(NULL),
	fHistTOFbefore(NULL),
	fHistTOFafter(NULL),
	fHistTrackDCAxyPtbefore(NULL),
	fHistTrackDCAxyPtafter(NULL),
	fHistTrackDCAzPtbefore(NULL),
	fHistTrackDCAzPtafter(NULL),
	fHistTrackNFindClsPtTPCbefore(NULL),
	fHistTrackNFindClsPtTPCafter(NULL),
	fHistTrackSelectedEta(NULL),
	fHistTrackSelectedPhi(NULL),
	fHistTrackSelectedPt(NULL),
	fHistTrackSelectedPtWithoutITS(NULL),
	fStringITSClusterCut(""),
	fPeriodName(ref.fPeriodName)
{
	for(Int_t jj=0;jj<kNCuts;jj++){ fCuts[jj]=0; }
	fCutString=new TObjString((GetCutNumber()).Data());

	// Using standard function for setting Cuts
	if (fEsdTrackCuts==NULL) fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
}

//________________________________________________________________________
AliPrimaryPionCuts::~AliPrimaryPionCuts() {
		// Destructor
	//Deleting fHistograms leads to seg fault it it's added to output collection of a task
	// if(fHistograms)
	// 	delete fHistograms;
	// fHistograms = NULL;

	if(fCutString != NULL){
		delete fCutString;
		fCutString = NULL;
	}
}

//________________________________________________________________________
void AliPrimaryPionCuts::InitCutHistograms(TString name, Bool_t preCut,TString cutNumber){

	// Initialize Cut Histograms for QA (only initialized and filled if function is called)

	TString cutName = "";
	
	if( cutNumber==""){
		cutName = GetCutNumber().Data();
	}
	else {
		cutName = cutNumber.Data();
	} 

	if(fHistograms != NULL){
		delete fHistograms;
		fHistograms=NULL;
	}
	if(fHistograms==NULL){
		fHistograms=new TList();
		if(name=="")fHistograms->SetName(Form("PionCuts_%s",cutName.Data()));
		else fHistograms->SetName(Form("%s_%s",name.Data(),cutName.Data()));
	}


	fHistCutIndex=new TH1F(Form("IsPionSelected %s",cutName.Data()),"IsPionSelected",10,-0.5,9.5);
	fHistCutIndex->GetXaxis()->SetBinLabel(kPionIn+1,"in");
	fHistCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
	fHistCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
	fHistCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"dEdx");
	fHistCutIndex->GetXaxis()->SetBinLabel(kPionOut+1,"out");
	fHistograms->Add(fHistCutIndex);

	// dEdx Cuts
	fHistdEdxCuts=new TH1F(Form("dEdxCuts %s",cutName.Data()),"dEdxCuts",5,-0.5,4.5);
	fHistdEdxCuts->GetXaxis()->SetBinLabel(1,"in");
	fHistdEdxCuts->GetXaxis()->SetBinLabel(2,"ITSpion");
	fHistdEdxCuts->GetXaxis()->SetBinLabel(3,"TPCpion");
	fHistdEdxCuts->GetXaxis()->SetBinLabel(4,"TOFpion");
	fHistdEdxCuts->GetXaxis()->SetBinLabel(5,"out");
	fHistograms->Add(fHistdEdxCuts);
	
	TAxis *axisBeforeITS  = NULL;
	TAxis *axisBeforedEdx = NULL;
	TAxis *axisBeforeTOF  = NULL;
	TAxis *axisBeforedEdxSignal = NULL;
    if(!fDoLightOutput){
      fHistTPCdEdxbefore=new TH2F(Form("Pion_dEdx_before %s",cutName.Data()),"dEdx pion before" ,170,0.05,50,400,-10,10);
      fHistograms->Add(fHistTPCdEdxbefore);
      axisBeforedEdx = fHistTPCdEdxbefore->GetXaxis();

      fHistTPCdEdxSignalbefore=new TH2F(Form("Pion_dEdxSignal_before %s",cutName.Data()),"dEdx pion signal before" ,170,0.05,50.0,800,0.0,200);
      fHistograms->Add(fHistTPCdEdxSignalbefore);
      axisBeforedEdxSignal = fHistTPCdEdxSignalbefore->GetXaxis();
      if(preCut){
        fHistITSdEdxbefore=new TH2F(Form("Pion_ITS_before %s",cutName.Data()),"ITS dEdx pion before" ,170,0.05,50,400,-10,10);
        fHistograms->Add(fHistITSdEdxbefore);
        axisBeforeITS = fHistITSdEdxbefore->GetXaxis();



        fHistTOFbefore=new TH2F(Form("Pion_TOF_before %s",cutName.Data()),"TOF pion before" ,170,0.05,50,400,-6,10);
        fHistograms->Add(fHistTOFbefore);
        axisBeforeTOF = fHistTOFbefore->GetXaxis();

        fHistTrackDCAxyPtbefore = new TH2F(Form("hTrack_DCAxy_Pt_before %s",cutName.Data()),"DCAxy Vs Pt of tracks before",800,-4.0,4.0,400,0.,10.);
        fHistograms->Add(fHistTrackDCAxyPtbefore);

        fHistTrackDCAzPtbefore  = new TH2F(Form("hTrack_DCAz_Pt_before %s",cutName.Data()), "DCAz  Vs Pt of tracks before",800,-4.0,4.0,400,0.,10.);
        fHistograms->Add(fHistTrackDCAzPtbefore);

        fHistTrackNFindClsPtTPCbefore = new TH2F(Form("hTrack_NFindCls_Pt_TPC_before %s",cutName.Data()),"Track: N Findable Cls TPC Vs Pt before",100,0,1,400,0.,10.);
        fHistograms->Add(fHistTrackNFindClsPtTPCbefore);
      }

      fHistITSdEdxafter=new TH2F(Form("Pion_ITS_after %s",cutName.Data()),"ITS dEdx pion after" ,170,0.05,50,400, -10,10);
      fHistograms->Add(fHistITSdEdxafter);

      fHistTPCdEdxafter=new TH2F(Form("Pion_dEdx_after %s",cutName.Data()),"dEdx pion after" ,170,0.05,50,400, -10,10);
      fHistograms->Add(fHistTPCdEdxafter);

      fHistTPCdEdxSignalafter=new TH2F(Form("Pion_dEdxSignal_after %s",cutName.Data()),"dEdx pion signal after" ,170,0.05,50.0,800,0.0,200);
      fHistograms->Add(fHistTPCdEdxSignalafter);

      fHistTOFafter=new TH2F(Form("Pion_TOF_after %s",cutName.Data()),"TOF pion after" ,170,0.05,50,400,-6,10);
      fHistograms->Add(fHistTOFafter);

      fHistTrackDCAxyPtafter  = new TH2F(Form("hTrack_DCAxy_Pt_after %s",cutName.Data()),"DCAxy Vs Pt of tracks after",800,-4.0,4.0,400,0.,10.);
      fHistograms->Add(fHistTrackDCAxyPtafter);

      fHistTrackDCAzPtafter  = new TH2F(Form("hTrack_DCAz_Pt_after %s",cutName.Data()), "DCAz Vs Pt of tracks  after",800,-4.0,4.0,400,0.,10.);
      fHistograms->Add(fHistTrackDCAzPtafter);

      fHistTrackNFindClsPtTPCafter = new TH2F(Form("hTrack_NFindCls_Pt_TPC_after %s",cutName.Data()),"Track: N Findable Cls TPC Vs Pt after",100,0,1,400,0.,10.);
      fHistograms->Add(fHistTrackNFindClsPtTPCafter);

			fHistTrackSelectedEta = new TH1F(Form("fHistTrackSelectedEta %s",cutName.Data()),"Selected Track Eta",200,-1,1);
      fHistograms->Add(fHistTrackSelectedEta);

			fHistTrackSelectedPhi = new TH1F(Form("fHistTrackSelectedPhi %s",cutName.Data()),"Selected Track Phi",200,0.,2*TMath::Pi());
      fHistograms->Add(fHistTrackSelectedPhi);

			fHistTrackSelectedPt = new TH1F(Form("fHistTrackSelectedPt %s",cutName.Data()),"Selected Track Pt",500,0.,50.);
      fHistograms->Add(fHistTrackSelectedPt);

			fHistTrackSelectedPtWithoutITS = new TH1F(Form("fHistTrackSelectedPtWithoutITS %s",cutName.Data()),"Selected Track Pt w/o ITS refit and cluster requirement",500,0.,50.);
      fHistograms->Add(fHistTrackSelectedPtWithoutITS);
    }
    if(!fDoLightOutput){
      TAxis *AxisAfter = fHistTPCdEdxafter->GetXaxis();
      Int_t bins = AxisAfter->GetNbins();
      Double_t from = AxisAfter->GetXmin();
      Double_t to = AxisAfter->GetXmax();
      Double_t *newBins = new Double_t[bins+1];
      newBins[0] = from;
      Double_t factor = TMath::Power(to/from, 1./bins);
      for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
      AxisAfter->Set(bins, newBins);
      AxisAfter = fHistTOFafter->GetXaxis();
      AxisAfter->Set(bins, newBins);
      AxisAfter = fHistITSdEdxafter->GetXaxis();
      AxisAfter->Set(bins,newBins);
      AxisAfter = fHistTPCdEdxSignalafter->GetXaxis();
      AxisAfter->Set(bins,newBins);

      if (axisBeforedEdx) axisBeforedEdx->Set(bins, newBins);
      if (axisBeforedEdxSignal) axisBeforedEdxSignal->Set(bins,newBins);
      if(preCut){
        axisBeforeITS->Set(bins, newBins);
        axisBeforeTOF->Set(bins, newBins);

      }

      delete [] newBins;
    }
	// Event Cuts and Info
}


//________________________________________________________________________
Bool_t AliPrimaryPionCuts::InitPIDResponse(){

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
Bool_t AliPrimaryPionCuts::PionIsSelectedMC(Int_t labelParticle,AliMCEvent *mcEvent){
	
    if( labelParticle < 0 || labelParticle >= mcEvent->GetNumberOfTracks() ) return kFALSE;
// 	if( mcEvent->IsPhysicalPrimary(labelParticle) == kFALSE ) return kFALSE;  // moved to actual tasks

    TParticle* particle = mcEvent->Particle(labelParticle);

	if( TMath::Abs( particle->GetPdgCode() ) != 211 )  return kFALSE;
	
	if( fDoEtaCut ){
	if( particle->Eta() > (fEtaCut + fEtaShift) || particle->Eta() < (-fEtaCut + fEtaShift) )
		return kFALSE;
	}
	
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::PionIsSelectedAODMC(Int_t labelParticle,TClonesArray *AODMCTrackArray){

    if( labelParticle < 0 || labelParticle >= AODMCTrackArray->GetSize()) return kFALSE;
// 	if( mcEvent->IsPhysicalPrimary(labelParticle) == kFALSE ) return kFALSE;  // moved to actual tasks

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelParticle));

    if( TMath::Abs( particle->GetPdgCode() ) != 211 )  return kFALSE;
    if( fDoEtaCut ){
    if( particle->Eta() > (fEtaCut + fEtaShift) || particle->Eta() < (-fEtaCut + fEtaShift) )
        return kFALSE;
    }

    return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::PionIsSelected(AliESDtrack* lTrack){
    //Selection of Reconstructed electrons
	
	Float_t b[2];
	Float_t bCov[3];
	if (lTrack == NULL){
		if (fHistCutIndex) fHistCutIndex->Fill(kNoTracks);
		return kFALSE;  
	}   
	lTrack->GetImpactParameters(b,bCov);

	if (bCov[0]<=0 || bCov[2]<=0) {
		AliDebug(1, "Estimated b resolution lower or equal zero!");
		bCov[0]=0; bCov[2]=0;
	}

	Float_t dcaToVertexXY = b[0];
	Float_t dcaToVertexZ  = b[1];
	Double_t clsToF = GetNFindableClustersTPC(lTrack);


	if (fHistCutIndex) fHistCutIndex->Fill(kPionIn);

	if (fHistTrackDCAxyPtbefore) fHistTrackDCAxyPtbefore->Fill(dcaToVertexXY,lTrack->Pt());
	if (fHistTrackDCAzPtbefore) fHistTrackDCAzPtbefore->Fill( dcaToVertexZ, lTrack->Pt());
	if (fHistTrackNFindClsPtTPCbefore) fHistTrackNFindClsPtTPCbefore->Fill( clsToF, lTrack->Pt());

	
	//if ( ! lTrack->GetConstrainedParam() ){
  //    return kFALSE;
	//}
	AliVTrack * track = dynamic_cast<AliVTrack*>(lTrack);

	// Track Cuts
	if( !TrackIsSelected(lTrack) ){
		if (fHistCutIndex) fHistCutIndex->Fill(kTrackCuts);
		return kFALSE;
	}

	// dEdx Cuts
	if( ! dEdxCuts( track ) ) {
		if(fHistCutIndex)fHistCutIndex->Fill(kdEdxCuts);
		return kFALSE;
	}

	//Pion passed the cuts
	if (fHistCutIndex) fHistCutIndex->Fill(kPionOut);
	if (fHistTrackDCAxyPtafter) fHistTrackDCAxyPtafter->Fill(dcaToVertexXY,lTrack->Pt());
	if (fHistTrackDCAzPtafter) fHistTrackDCAzPtafter->Fill(dcaToVertexZ,lTrack->Pt());
	if (fHistTrackNFindClsPtTPCafter) fHistTrackNFindClsPtTPCafter->Fill( clsToF, lTrack->Pt());

	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::PionIsSelectedAOD(AliAODTrack* lTrack){
    //Selection of Reconstructed electrons

    Float_t b[2];
    Float_t bCov[3];
    if (lTrack == NULL){
        if (fHistCutIndex) fHistCutIndex->Fill(kNoTracks);
        return kFALSE;
    }


    if (fHistCutIndex) fHistCutIndex->Fill(kPionIn);


    AliVTrack * track = dynamic_cast<AliVTrack*>(lTrack);

    // Track Cuts
    if( !TrackIsSelectedAOD(lTrack) ){
        if (fHistCutIndex) fHistCutIndex->Fill(kTrackCuts);
        return kFALSE;
    }

    // this would throw a warning for AOD tracks that are not defined at the X of vertex (mainly tracks that do
    // not pass the track selection). The plotting of DCAxy and DCAz BEFORE has therefore been moved to after the track selection
    // for AODs to avoid undefined x values.
    lTrack->GetImpactParameters(b,bCov);

    if (bCov[0]<=0 || bCov[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal zero!");
        bCov[0]=0; bCov[2]=0;
    }

    Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ  = b[1];
    Double_t clsToF = GetNFindableClustersTPC(lTrack);


    if (fHistTrackDCAxyPtbefore) fHistTrackDCAxyPtbefore->Fill(dcaToVertexXY,lTrack->Pt());
    if (fHistTrackDCAzPtbefore) fHistTrackDCAzPtbefore->Fill( dcaToVertexZ, lTrack->Pt());
    if (fHistTrackNFindClsPtTPCbefore) fHistTrackNFindClsPtTPCbefore->Fill( clsToF, lTrack->Pt());

    // dEdx Cuts
    if( ! dEdxCuts( track ) ) { // TODO: check if these still work for AODs
        if(fHistCutIndex)fHistCutIndex->Fill(kdEdxCuts);
        return kFALSE;

    }


    //Pion passed the cuts
    if (fHistCutIndex) fHistCutIndex->Fill(kPionOut);
    if (fHistTrackDCAxyPtafter) fHistTrackDCAxyPtafter->Fill(dcaToVertexXY,lTrack->Pt());
    if (fHistTrackDCAzPtafter) fHistTrackDCAzPtafter->Fill(dcaToVertexZ,lTrack->Pt());
    if (fHistTrackNFindClsPtTPCafter) fHistTrackNFindClsPtTPCafter->Fill( clsToF, lTrack->Pt());

    return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::TrackIsSelected(AliESDtrack* lTrack) {
  // Track Selection for Photon Reconstruction
  Double_t clsToF = GetNFindableClustersTPC(lTrack);
  Double_t frac_SharedClus = Double_t(lTrack->GetTPCnclsS()) / Double_t(lTrack->GetTPCncls());
  if( ! fEsdTrackCuts->AcceptTrack(lTrack) && ! fEsdTrackCutsGC->AcceptTrack(lTrack)){
    return kFALSE;
  }

	// Absolute TPC cluster cut
	// (should be already applied in fEsdTrackCuts, however it might
	// not be properly applied together with pTDependent cut )
	if(lTrack->GetTPCNcls()<fMinClsTPC) return kFALSE;

  if (frac_SharedClus > fMaxSharedClsTPCFrac){ return kFALSE; }
 
  if( fDoEtaCut ) {
    if(  lTrack->Eta() > (fEtaCut + fEtaShift) || lTrack->Eta() < (-fEtaCut + fEtaShift) ) {
      return kFALSE;
    }
  }

  if( lTrack->Pt() < fPtCut ) {
    return kFALSE;
  }

  if( clsToF < fMinClsTPCToF){
    return kFALSE;
  }

  if(!fDoLightOutput){
  	fHistTrackSelectedEta->Fill(lTrack->Eta());
  	fHistTrackSelectedPhi->Fill(lTrack->Phi());
  	fHistTrackSelectedPt->Fill(lTrack->Pt());

		if(lTrack->GetNumberOfITSClusters()==0) fHistTrackSelectedPtWithoutITS->Fill(lTrack->Pt());
	}

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliPrimaryPionCuts::TrackIsSelectedAOD(AliAODTrack* lTrack) {
  // Track Selection for Photon Reconstruction
  Double_t clsToF = GetNFindableClustersTPC(lTrack);
  Double_t frac_SharedClus = Double_t(lTrack->GetTPCnclsS()) / Double_t(lTrack->GetTPCncls());
  // apply filter bits 
  if( ! lTrack->IsHybridGlobalConstrainedGlobal()){
    return kFALSE;
  }

	if(fRequireVertexConstrain && (! lTrack->IsGlobalConstrained())){
    return kFALSE;
	}

	// since fEsdTrackCuts->AcceptTrack() is not available for AODTracks
	// the following cuts will be applied manually
	// Note that they are only effective if they are stronger 
	// than the cuts already applied on AOD refiltering level

	// Absolute TPC Cluster cut
    if (frac_SharedClus > fMaxSharedClsTPCFrac){ return kFALSE; }
	if(lTrack->GetTPCNcls()<fMinClsTPC) return kFALSE;
	if(lTrack->GetTPCchi2perCluster()>fChi2PerClsTPC) return kFALSE;
  // DCA cut 
  if(!IsDCACutAccepted(lTrack)) return kFALSE;

  // ITS Cluster Cut
	// SetClusterRequirementITS and SetRequireITSRefit can
	// not be set for AODs after filtering
	if(lTrack->GetITSNcls()<fMinClsITS) return kFALSE;

  if( fDoEtaCut ) {
    if(  lTrack->Eta() > (fEtaCut + fEtaShift) || lTrack->Eta() < (-fEtaCut + fEtaShift) ) {
      return kFALSE;
    }
  }

  if( lTrack->Pt() < fPtCut ) {
    return kFALSE;
  }

  if( clsToF < fMinClsTPCToF){
    return kFALSE;
  }

  if(!fDoLightOutput){
		fHistTrackSelectedEta->Fill(lTrack->Eta());
  	fHistTrackSelectedPhi->Fill(lTrack->Phi());
  	fHistTrackSelectedPt->Fill(lTrack->Pt());
		if(lTrack->GetITSNcls()==0) fHistTrackSelectedPtWithoutITS->Fill(lTrack->Pt());
	}

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::dEdxCuts(AliVTrack *fCurrentTrack){

	// Pion Identification Cuts for Photon reconstruction

	if(!fPIDResponse){  InitPIDResponse();  }// Try to reinitialize PID Response
	if(!fPIDResponse){  AliError("No PID Response"); return kFALSE;}// if still missing fatal error

	Int_t cutIndex=0;

	if(fHistdEdxCuts)fHistdEdxCuts->Fill(cutIndex);
	if(fHistITSdEdxbefore)fHistITSdEdxbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kPion));
	if(fHistTPCdEdxbefore)fHistTPCdEdxbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kPion));
	if(fHistTPCdEdxSignalbefore)fHistTPCdEdxSignalbefore->Fill(fCurrentTrack->P(),TMath::Abs(fCurrentTrack->GetTPCsignal()));

	cutIndex++;

	if( fDodEdxSigmaITSCut == kTRUE ){
		if( fPIDResponse->NumberOfSigmasITS(fCurrentTrack,AliPID::kPion)<fPIDnSigmaBelowPionLineITS ||
				fPIDResponse->NumberOfSigmasITS(fCurrentTrack,AliPID::kPion)> fPIDnSigmaAbovePionLineITS ){	
			if(fHistdEdxCuts)fHistdEdxCuts->Fill(cutIndex);
			return kFALSE;
		}
	}
		
	if(fHistITSdEdxafter)fHistITSdEdxafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kPion));
	
	
	cutIndex++;
		
		
	if(fDodEdxSigmaTPCCut == kTRUE){
		// TPC Pion Line
		if( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaBelowPionLineTPC ||
			fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)>fPIDnSigmaAbovePionLineTPC){
			if(fHistdEdxCuts)fHistdEdxCuts->Fill(cutIndex);
			return kFALSE;
		}
		cutIndex++;
	} else { cutIndex+=1; }
	
	if( ( fCurrentTrack->GetStatus() & AliESDtrack::kTOFpid ) && ( !( fCurrentTrack->GetStatus() & AliESDtrack::kTOFmismatch) ) ){
		if(fHistTOFbefore) fHistTOFbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kPion));
		if(fUseTOFpid){
			if( fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kPion)>fPIDnSigmaAbovePionLineTOF ||
				fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kPion)<fPIDnSigmaBelowPionLineTOF ){
				if(fHistdEdxCuts)fHistdEdxCuts->Fill(cutIndex);
				return kFALSE;
			}
		}
		if(fHistTOFafter)fHistTOFafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kPion));
	} else if ( fRequireTOF == kTRUE ) {
		if(fHistdEdxCuts)fHistdEdxCuts->Fill(cutIndex);
		return kFALSE;
	}
	cutIndex++;
		
	if(fHistdEdxCuts)fHistdEdxCuts->Fill(cutIndex);
	if(fHistTPCdEdxafter)fHistTPCdEdxafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kPion));
	if(fHistTPCdEdxSignalafter)fHistTPCdEdxSignalafter->Fill(fCurrentTrack->P(),TMath::Abs(fCurrentTrack->GetTPCsignal()));
	
	return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliPrimaryPionCuts::GetTrack(AliVEvent * event, Int_t label){
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
	cout << "track not found " << label << " " << event->GetNumberOfTracks() << endl;
	return NULL;
}

///________________________________________________________________________
Double_t AliPrimaryPionCuts::GetNFindableClustersTPC(AliVTrack* lTrack){

	Double_t clsToF=0;
	if ( !fUseCorrectedTPCClsInfo ){
		if(lTrack->GetTPCNclsF()!=0){
			clsToF = (Double_t)lTrack->GetNcls(1)/(Double_t)lTrack->GetTPCNclsF();
		}
	} else {
		clsToF = lTrack->GetTPCClusterInfo(2,0); //NOTE ask friederike	
	}
	return clsToF;

}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::UpdateCutString() {
///Update the cut string (if it has been created yet)

	if(fCutString && fCutString->GetString().Length() == kNCuts) {
		fCutString->SetString(GetCutNumber());
	} else {
		return kFALSE;
	}
	return kTRUE;

}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());
  
  // Set basic cuts for AOD compability
	SetHybridTrackCutsAODFiltering(fRunFlag);
	AliInfo(Form("Presetting ESD cuts with prefiltering for runflag %d",fRunFlag));

	// Initialize Cuts from a given Cut string

	AliInfo(Form("Set PionCuts Number: %s",analysisCutSelection.Data()));
	
	if(analysisCutSelection.Length()!=kNCuts) {
		AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
		return kFALSE;
	}
	if(!analysisCutSelection.IsAlnum()){
		AliError("Cut selection is not alphanumeric");
		return kFALSE;
	}
	
  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
	const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
	for(Int_t ii=0;ii<kNCuts;ii++){
		ASSIGNARRAY(ii);
	}

	// Set Individual Cuts
	for(Int_t ii=0;ii<kNCuts;ii++){
		if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
	}

	fEsdTrackCutsGC = (AliESDtrackCuts*) fEsdTrackCuts->Clone();
	
	if(fRunFlag==1500){
		fEsdTrackCutsGC->SetRequireITSRefit(kTRUE);
	} else{
		fEsdTrackCutsGC->SetRequireITSRefit(kFALSE);
	}
	fEsdTrackCutsGC->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);


	PrintCutsWithValues();
	return kTRUE;
}
///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetCut(cutIds cutID, const Int_t value) {
	///Set individual cut ID

	//cout << "Updating cut  " << fgkCutNames[cutID] << " (" << cutID << ") to " << value << endl;

	switch (cutID) {
		case kEtaCut:
			if( SetEtaCut(value)) {
				fCuts[kEtaCut] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kClsITSCut:
			if( SetITSClusterCut(value) ) {
				fCuts[kClsITSCut] = value;
				UpdateCutString();
				return kTRUE;		 	
			} else return kFALSE;
		case kClsTPCCut:
			if( SetTPCClusterCut(value)) {
				fCuts[kClsTPCCut] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kDCACut:
			if( SetDCACut(value)) {
				fCuts[kDCACut] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kPtCut: 	
			if( SetPtCut(value)) {
				fCuts[kPtCut] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kPidedxSigmaITSCut:
			if( SetITSdEdxCutPionLine(value)) { 
				fCuts[kPidedxSigmaITSCut] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kPidedxSigmaTPCCut:
			if( SetTPCdEdxCutPionLine(value)) { 
				fCuts[kPidedxSigmaTPCCut] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kPiTOFSigmaPID:
			if( SetTOFPionPIDCut(value)) {
				fCuts[kPiTOFSigmaPID] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		case kMassCut:
			if( SetMassCut(value)) {
				fCuts[kMassCut] = value;
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
void AliPrimaryPionCuts::PrintCuts() {
    // Print out current Cut Selection
	for(Int_t ic = 0; ic < kNCuts; ic++) {
		printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
	}
}

///________________________________________________________________________
void AliPrimaryPionCuts::PrintCutsWithValues() {
   // Print out current Cut Selection with value
	printf("\nCharged Pion cutnumber \n");
	for(Int_t ic = 0; ic < kNCuts; ic++) {
		printf("%d",fCuts[ic]);
	}
	printf("\n\n");

	printf("Acceptance cuts \n");
	if (fDoEtaCut) printf("\t |eta_{pi+-}| < %3.2f  \n", fEtaCut);
	else 	printf("\t none \n");
	printf("Track cuts \n");
	printf("\t %s \n", fStringITSClusterCut.Data());
	printf("\t min N cluster TPC > %3.2f \n", fMinClsTPC);
	printf("\t min N cluster TPC/ findable > %3.2f \n", fMinClsTPCToF);

    printf("\t max Chi2 per cluster TPC < %3.2f \n", fChi2PerClsTPC);
    printf("\t require TPC refit ? %d \n", fRequireTPCRefit);
// 	printf("\t dca > %3.2f \n", fMinClsTPCToF);
// 	"kDCAcut",				// 3
	printf("\t min pT > %3.2f \n", fPtCut);
	printf("PID cuts \n");
	if (fDodEdxSigmaITSCut)printf("\t %3.2f < ITS n_sigma pi < %3.2f \n", fPIDnSigmaBelowPionLineITS, fPIDnSigmaAbovePionLineITS );
	if (fDodEdxSigmaTPCCut)printf("\t %3.2f < TPC n_sigma pi < %3.2f \n", fPIDnSigmaBelowPionLineTPC, fPIDnSigmaAbovePionLineTPC );
	if (fDoTOFsigmaCut)printf("\t %3.2f < TOF n_sigma pi < %3.2f \n", fPIDnSigmaBelowPionLineTOF, fPIDnSigmaAbovePionLineTOF );
	if (fDoMassCut) printf("two-pion mass cut < %3.2f \n", fMassCut);
	printf("\n\n");
}



///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetITSdEdxCutPionLine(Int_t ededxSigmaCut){ 
	switch(ededxSigmaCut){
		case 0: 
			fDodEdxSigmaITSCut = kFALSE;
			fPIDnSigmaBelowPionLineITS=-100;
			fPIDnSigmaAbovePionLineITS= 100;
			break;
		case 1: // -10,10
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-10;
			fPIDnSigmaAbovePionLineITS=10;
			break;
		case 2: // -6,7
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-6;
			fPIDnSigmaAbovePionLineITS=7;
			break;
		case 3: // -5,5
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-5;
			fPIDnSigmaAbovePionLineITS=5;
			break;
		case 4: // -4,5
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-4;
			fPIDnSigmaAbovePionLineITS=5;
			break;
		case 5: // -3,5
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-3;
			fPIDnSigmaAbovePionLineITS=5;
			break;
		case 6: // -4,4
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-4;
			fPIDnSigmaAbovePionLineITS=4;
			break;
		case 7: // -2.5,4
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-2.5;
			fPIDnSigmaAbovePionLineITS=4;
			break;
		case 8: // -2,3.5
			fDodEdxSigmaITSCut = kTRUE;
			fPIDnSigmaBelowPionLineITS=-2;
			fPIDnSigmaAbovePionLineITS=3.5;
			break;
		default:
			cout<<"Warning: ITSdEdxCutPionLine not defined"<<ededxSigmaCut<<endl;
			return kFALSE;
		
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetTPCdEdxCutPionLine(Int_t ededxSigmaCut){
	switch(ededxSigmaCut){
		case 0: 
			fDodEdxSigmaTPCCut = kFALSE;
			fPIDnSigmaBelowPionLineTPC=-10;
			fPIDnSigmaAbovePionLineTPC=10;
			break;
		case 1: // -10,10
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-10;
			fPIDnSigmaAbovePionLineTPC=10;
			break;
		case 2: // -6,7
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-6;
			fPIDnSigmaAbovePionLineTPC=7;
			break;
		case 3: // -5,5
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-5;
			fPIDnSigmaAbovePionLineTPC=5;
			break;
		case 4: // -4,5
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-4;
			fPIDnSigmaAbovePionLineTPC=5;
			break;	
		case 5: // -4,4
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-4;
			fPIDnSigmaAbovePionLineTPC=4;
			break;
		case 6: // -3,4
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-3;
			fPIDnSigmaAbovePionLineTPC=4;
			break;
		case 7: // -3,3
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-3;
			fPIDnSigmaAbovePionLineTPC=3;
			break;
		case 8: // -2,3.
			fDodEdxSigmaTPCCut = kTRUE;
			fPIDnSigmaBelowPionLineTPC=-2;
			fPIDnSigmaAbovePionLineTPC=3.;
			break;
        case 9: // -2.5, 2.5
            fDodEdxSigmaTPCCut = kTRUE;
            fPIDnSigmaBelowPionLineTPC=-2;
            fPIDnSigmaAbovePionLineTPC=3.;
            break;
        case 10: //a -3.5, 3.5
            fDodEdxSigmaTPCCut = kTRUE;
            fPIDnSigmaBelowPionLineTPC=-2;
            fPIDnSigmaAbovePionLineTPC=3.;
            break;
        case 11: //b -2., 2.
            fDodEdxSigmaTPCCut = kTRUE;
            fPIDnSigmaBelowPionLineTPC=-2;
            fPIDnSigmaAbovePionLineTPC=3.;
            break;
		default:
			cout<<"Warning: TPCdEdxCutPionLine not defined"<<ededxSigmaCut<<endl;
			return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetITSClusterCut(Int_t clsITSCut){
	if( !fEsdTrackCuts ) {
		cout<<"Warning: AliESDtrackCut is not initialized "<<endl;
		return kFALSE;
	}

	switch(clsITSCut){
		case 0: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
			fStringITSClusterCut= "no SPD cluster requirement";
			break;
		case 1: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
			fStringITSClusterCut= "first SPD cluster required";
			break;  //1 hit first layer of SPD
		case 2: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetRequireITSRefit(kTRUE);
			fStringITSClusterCut= "first or second SPD cluster required";
			break; //1 hit in any layer of SPD
		case 3: 
		  fMinClsITS = 4;
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
			fEsdTrackCuts->SetMinNClustersITS(fMinClsITS);
			fStringITSClusterCut= "first SPD cluster required, min number of ITS clusters = 4";
			// 4 hits in total in the ITS. At least 1 hit in the first layer of SPD  
			break;
		case 4: 
		  fMinClsITS = 3;
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(fMinClsITS);
			fStringITSClusterCut= "first or second SPD cluster required, min number of ITS clusters = 3";
			// 3 hits in total in the ITS. At least 1 hit in any layer of SPD
			break;
		case 5: 
		 	fMinClsITS = 4;
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(fMinClsITS);
			fStringITSClusterCut= "first or second SPD cluster required, min number of ITS clusters = 4";
			// 4 hits in total in the ITS. At least 1 hit in any layer of SPD
			break;
		case 6: 
			fMinClsITS = 5;
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(fMinClsITS);
			fStringITSClusterCut= "first or second SPD cluster required, min number of ITS clusters = 5";
			// 5 hits in total in the ITS. At least 1 hit in any layer of SPD
			break;
		case 7: 
			fMinClsITS = 1;
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(fMinClsITS);
			fStringITSClusterCut= "first or second SPD cluster required, min number of ITS clusters = 5";
			// 5 hits in total in the ITS. At least 1 hit in any layer of SPD
			break;
		default:
			cout<<"Warning: clsITSCut not defined "<<clsITSCut<<endl;
			return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetTPCClusterCut(Int_t clsTPCCut){  
	switch(clsTPCCut){
		case 0: // 0
			fMinClsTPC= 0.;
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			break;
		case 1:  // 70
			fMinClsTPC= 70.;
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			break;
		case 2:  // 80
			fMinClsTPC= 80.;
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			break;
		case 3:  // 100
			fMinClsTPC= 100.;
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			break;
		case 4:  // 0% of findable clusters
				fMinClsTPC= 70.;  
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			fMinClsTPCToF= 0.0;
			fUseCorrectedTPCClsInfo=0;
			break;
		case 5:  // 35% of findable clusters
			fMinClsTPC = 70.;  
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			fMinClsTPCToF= 0.35;
			fUseCorrectedTPCClsInfo=0;
			break;
		case 6:  // 60% of findable clusters
			fMinClsTPC= 70.;  
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			fMinClsTPCToF= 0.6;
			fUseCorrectedTPCClsInfo=0;
			break;
		case 7:  // 70% of findable clusters
			fMinClsTPC= 70.;  
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			fMinClsTPCToF= 0.7;
			fUseCorrectedTPCClsInfo=0;
			break;
		case 8: fMinClsTPC = 0.;  
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			fMinClsTPCToF= 0.35;
			fUseCorrectedTPCClsInfo=0;
			break;
		case 9:  // 35% of findable clusters
			fMinClsTPC = 70.;  
			fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
			fMinClsTPCToF= 0.35;
			fUseCorrectedTPCClsInfo=1;
			break;
        case 10: //a
            fMinClsTPC     = 80.;
            fChi2PerClsTPC = 4;
            fRequireTPCRefit    = kTRUE;
            fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
            // Other Cuts concerning TPC
            fEsdTrackCuts->SetMaxChi2PerClusterTPC(fChi2PerClsTPC);
            fEsdTrackCuts->SetRequireTPCRefit(fRequireTPCRefit);
            break;
        case 11: //b settings as in PHOS public omega
            fMinClsTPC     = 70.;
            fChi2PerClsTPC = 4;
            fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
            fEsdTrackCuts->SetMaxChi2PerClusterTPC(fChi2PerClsTPC);
            break;
        case 12:  //c 80 + refit
            fMinClsTPC= 80.;
            fRequireTPCRefit    = kTRUE;
            fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
            break;
        case 13:  //d 80 + refit + vertex constrain (only for AOD)
				    fRequireVertexConstrain = kTRUE;
            fMinClsTPC= 80.;
            fRequireTPCRefit    = kTRUE;
            fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
            break;
        case 14:  //e 80 + refit, Shared cluster Fraction =0
            fMinClsTPC= 80.;
            fRequireTPCRefit    = kTRUE;
            fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
            fMaxSharedClsTPCFrac=0.;
            break;
        case 15:  //f 80 + refit, Shared cluster Fraction <=0.4
            fMinClsTPC= 80.;
            fRequireTPCRefit    = kTRUE;
            fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
            fMaxSharedClsTPCFrac=0.4;
            break;
				default:
						cout<<"Warning: clsTPCCut not defined "<<clsTPCCut<<endl;
						return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetEtaCut(Int_t etaCut){ 
	// Set eta Cut
	switch(etaCut){
		case 0: 
			fEtaCut = 100.;
			fDoEtaCut = kFALSE;
			break;
		case 1:	// 1.4
			fEtaCut	= 1.4;
			fDoEtaCut = kTRUE;
			break;
		case 2:	// 1.2
			fEtaCut	= 1.2;
			fDoEtaCut = kTRUE;
			break;
		case 3: // 0.9
			fEtaCut	= 0.9;
			fDoEtaCut = kTRUE;
			break;
		case 4: // 0.8
			fEtaCut	= 0.8;
			fDoEtaCut = kTRUE;
			break;
		case 5: // 0.75
			fEtaCut	= 0.75;
			fDoEtaCut = kTRUE;
			break;
		case 6: //0.6
			fEtaCut = 0.6; //changed from 0.4 to 0.6 2013.06.10
			fDoEtaCut = kTRUE;
			break;
		case 7: //0.5
			fEtaCut = 0.5; //changed from 0.3 to 0.5 2013.06.10
			fDoEtaCut = kTRUE;
			break;
		case 8: 
			fEtaCut = 0.4;
			fDoEtaCut = kTRUE;
			break;
		default:
			cout<<"Warning: EtaCut not defined "<<etaCut<<endl;
			return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetPtCut(Int_t ptCut){ 
	switch(ptCut){	  
		case 0: 
			fPtCut = 0.075;		
			break;
		case 1:	 // 0.1
			fPtCut	= 0.1; 	
			break;
		case 2:	 // 0.125 GeV
			fPtCut	= 0.125;		
			break;
		case 3: // 0.15 GeV
			fPtCut	= 0.15;
			break;
        case 4: // 0.40 GeV
            fPtCut  = 0.40;
            break;
		default:
			cout<<"Warning: PtCut not defined "<<ptCut<<endl;
			return kFALSE;
	}
	return kTRUE;
}


///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetDCACut(Int_t dcaCut)
{ 
	// Set DCA Cut
	if( !fEsdTrackCuts ) {
		cout<<"Warning: AliESDtrackCut is not initialized "<<endl;
		return kFALSE;
	}
  
	switch(dcaCut){	  
		case 0: 
			//Open cuts//
			fMaxDCAToVertexZ=1000;
			fMaxDCAToVertexXY=1000;
			fEsdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
			fEsdTrackCuts->SetMaxDCAToVertexXY(fMaxDCAToVertexXY);
			fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
			break;
    case 1:
			fUsePtDepXYDCA=kTRUE;
			fMaxDCAToVertexXYPtDep = "0.0182+0.0350/pt^1.01";
			fEsdTrackCuts->SetMaxDCAToVertexXYPtDep(fMaxDCAToVertexXYPtDep.Data());
			fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
			break;
		case 2: 
		  fMaxDCAToVertexZ=2;
			fMaxDCAToVertexXY=1;
			fEsdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
			fEsdTrackCuts->SetMaxDCAToVertexXY(fMaxDCAToVertexXY);
			fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
			break; 
        case 3:
            fMaxDCAToVertexZ = 3.0;
						fUsePtDepXYDCA=kTRUE;
						fMaxDCAToVertexXYPtDep = "0.0182+0.0350/pt^1.01";
            fEsdTrackCuts->SetMaxDCAToVertexXYPtDep(fMaxDCAToVertexXYPtDep.Data());
            fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
            fEsdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
            break;
        case 4:
				    fMaxDCAToVertexZ=3;
			      fMaxDCAToVertexXY=0.5;
            fEsdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
            fEsdTrackCuts->SetMaxDCAToVertexXY(fMaxDCAToVertexXY);
            fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
            break;
				case 5:
					 fMaxDCAToVertexZ=3.2;
			     fMaxDCAToVertexXY=2.4;
					 fUseDCAToVertex2D=kTRUE;
				   fEsdTrackCuts->SetMaxDCAToVertexXY(fMaxDCAToVertexXY);
           fEsdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
           fEsdTrackCuts->SetDCAToVertex2D(fUseDCAToVertex2D);
					 break;
				case 6: // temp
					 fMaxDCAToVertexZ=0.5;
			     fMaxDCAToVertexXY=0.5;
					 fUseDCAToVertex2D=kTRUE;
				   fEsdTrackCuts->SetMaxDCAToVertexXY(fMaxDCAToVertexXY);
           fEsdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAToVertexZ);
					 break;
		default:
			cout<<"Warning: dcaCut not defined "<<dcaCut<<endl;
			return kFALSE;
	}
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetTOFPionPIDCut(Int_t TOFelectronPID){
    // Set Cut
	switch(TOFelectronPID){ 
		case 0: // no cut
			fRequireTOF = kFALSE;
			fUseTOFpid = kFALSE;
			fPIDnSigmaBelowPionLineTOF=-100;
			fPIDnSigmaAbovePionLineTOF=100;
			break;
		case 1: // -7,7
			fRequireTOF = kFALSE;
			fUseTOFpid = kTRUE;
			fPIDnSigmaBelowPionLineTOF=-7;
			fPIDnSigmaAbovePionLineTOF=7;
			break;
		case 2: // -5,5
			fRequireTOF = kFALSE;
			fUseTOFpid = kTRUE;
			fPIDnSigmaBelowPionLineTOF=-5;
			fPIDnSigmaAbovePionLineTOF=5;
			break;
		case 3: // -3,5
			fRequireTOF = kFALSE;
			fUseTOFpid = kTRUE;
			fPIDnSigmaBelowPionLineTOF=-3;
			fPIDnSigmaAbovePionLineTOF=5;
			break;
		case 4: // -2,3
			fRequireTOF = kFALSE;
			fUseTOFpid = kTRUE;
			fPIDnSigmaBelowPionLineTOF=-2;
			fPIDnSigmaAbovePionLineTOF=3;
			break;
		case 5: // -3, 3 TOF mandatory
			fRequireTOF = kTRUE;
			fUseTOFpid  = kTRUE;
			fPIDnSigmaBelowPionLineTOF= -3;
			fPIDnSigmaAbovePionLineTOF=  3;
			break;
		default:
			cout<<"Warning: TOFPionCut not defined "<<TOFelectronPID<<endl;
		return kFALSE;
    } 
    return kTRUE;
}

///________________________________________________________________________
Bool_t AliPrimaryPionCuts::SetMassCut(Int_t massCut){
    // Set Cut
	switch(massCut){ 
		case 0: // no cut
			fDoMassCut = kFALSE;
			fMassCut = 10;
			break;
		case 1: // cut at 1 GeV/c^2
			fDoMassCut = kTRUE;
			fMassCut = 1;
			break;
		case 2: // cut at 0.7 GeV/c^2
			fDoMassCut = kTRUE;
			fMassCut = 0.75;
			break;
		case 3: // cut at 0.6 GeV/c^2
			fDoMassCut = kTRUE;
			fMassCut = 0.6;
			break;
		case 4: // cut at eta mass
			fDoMassCut = kTRUE;
			fMassCut = 0.547853;
			break;
		case 5: // cut at 0.5 GeV/c^2
			fDoMassCut = kTRUE;
			fMassCut = 0.5;
			break;
        case 6: // cut at 0.65 GeV/c^2
            fDoMassCut = kTRUE;
            fMassCut = 0.65;
            break;
        case 7: // cut at 0.7 GeV/c^2
            fDoMassCut = kTRUE;
            fMassCut = 0.7;
            break;
        case 8: // cut at 0.85 GeV/c^2
            fDoMassCut = kTRUE;
            fMassCut = 0.85;
            break;
        case 9: // cut at 1.5 GeV/c^2
            fDoMassCut = kTRUE;
            fMassCut = 1.5;
            break;
        case 10: //a overload mass cut for chi2 of vParticle
            fUse4VecForMass = kTRUE;
            fDoMassCut = kTRUE;
            fMassCut = 0.85;
            break;
        case 11: //b overload mass cut for chi2 of vParticle
            fUse4VecForMass = kTRUE;
            fDoMassCut = kTRUE;
            fDoMassCut_WithNDM= kTRUE;
            fMassCut = 0.600;
            fMassCut_WithNDM = 0.600;
            break;
        case 12: //c overload mass cut for chi2 of vParticle
            fUse4VecForMass = kTRUE;
            fDoMassCut = kTRUE;
            fDoMassCut_WithNDM= kTRUE;
            fMassCut = 0.85;
            fMassCut_WithNDM = 1.0;
            break;
        case 13: //d overload mass cut for chi2 of vParticle
            fUse4VecForMass = kTRUE;
            fDoMassCut = kTRUE;
            fDoMassCut_WithNDM= kTRUE;
            fMassCut = 0.85;
            fMassCut_WithNDM = 0.85;
            break;
        case 14: //e overload mass cut for chi2 of vParticle
            fUse4VecForMass = kTRUE;
            fDoMassCut = kTRUE;
            fDoMassCut_WithNDM= kTRUE;
            fMassCut = 0.600;
            fMassCut_WithNDM = 0.600;
            break;

		default:
			cout<<"Warning: MassCut not defined "<<massCut<<endl;
		return kFALSE;
    } 
    return kTRUE;
}


///________________________________________________________________________
TString AliPrimaryPionCuts::GetCutNumber(){
	// returns TString with current cut number
	return fCutStringRead;
}


///________________________________________________________________________
AliPrimaryPionCuts* AliPrimaryPionCuts::GetStandardCuts2010PbPb(){
    //Create and return standard 2010 PbPb cuts
    AliPrimaryPionCuts *cuts=new AliPrimaryPionCuts("StandardCuts2010PbPb","StandardCuts2010PbPb");
    if(!cuts->InitializeCutsFromCutString("000000400")){
		cout<<"Warning: Initialization of Standardcuts2010PbPb failed"<<endl;	
	}
    return cuts;
}

///________________________________________________________________________
AliPrimaryPionCuts* AliPrimaryPionCuts::GetStandardCuts2010pp(){
    //Create and return standard 2010 PbPb cuts
    AliPrimaryPionCuts *cuts=new AliPrimaryPionCuts("StandardCuts2010pp","StandardCuts2010pp");
                                          
    if(!cuts->InitializeCutsFromCutString("000000400")){
		cout<<"Warning: Initialization of Standardcuts2010pp failed"<<endl;
	}
    return cuts;
}

///________________________________________________________________________
void AliPrimaryPionCuts::SetHybridTrackCutsAODFiltering(Int_t runflag= 1000){
   // As preselection apply all cuts that are applied in AOD filtering
	 // so that ESD results are comparable to
	 // SetHybridFilterMaskGlobalConstrainedGlobal
    if(runflag == 0){
        fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    } else if(runflag==1){
		fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
	} else if(runflag==1500){
		fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
		fEsdTrackCuts->SetMaxDCAToVertexXY(2.4);
		fEsdTrackCuts->SetMaxDCAToVertexZ(3.2);
		fEsdTrackCuts->SetDCAToVertex2D(kTRUE);
		fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
		fEsdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
	} else if(runflag==1000){
		fEsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
		TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
		fEsdTrackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
		fEsdTrackCuts->SetMinNClustersTPC(70);
		fEsdTrackCuts->SetMaxChi2PerClusterTPC(4);
		fEsdTrackCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
		fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
		fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
		fEsdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
		// ITS
		fEsdTrackCuts->SetRequireITSRefit(kTRUE);
		//accept secondaries
		fEsdTrackCuts->SetMaxDCAToVertexXY(2.4);
		fEsdTrackCuts->SetMaxDCAToVertexZ(3.2);
		fEsdTrackCuts->SetDCAToVertex2D(kTRUE);
		//reject fakes
		fEsdTrackCuts->SetMaxChi2PerClusterITS(36);
		fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);	

		fEsdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	
		fEsdTrackCuts->SetEtaRange(-0.9,0.9);
		fEsdTrackCuts->SetPtRange(0.15, 1E+15);

	 	fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	} else{
		AliFatal(Form("Runflag %d is an invalid option for track preselection! ",runflag));
	}
}

//--------------------------------------------------------------------------
void AliPrimaryPionCuts::SetPtDepDCACuts(Double_t pt) {
  /// set the pt-dependent DCA cuts
  TString tmp = fMaxDCAToVertexXYPtDep;
  tmp.ReplaceAll("pt","x");
  TFormula CutMaxDCAToVertexXYPtDep("CutMaxDCAToVertexXYPtDep",tmp.Data());
   
  fMaxDCAToVertexXY=CutMaxDCAToVertexXYPtDep.Eval(pt);

  return;
}

//--------------------------------------------------------------------------
Bool_t AliPrimaryPionCuts::IsDCACutAccepted(AliAODTrack* lTrack) {
if(fUsePtDepXYDCA) SetPtDepDCACuts(lTrack->Pt());
  
	Float_t b[2];
  Float_t bCov[3];
  lTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }

  Float_t dcaToVertexXY = b[0];
  Float_t dcaToVertexZ = b[1];

  Float_t dcaToVertex = -1;
 
  if (fUseDCAToVertex2D){
    dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY/fMaxDCAToVertexXY/fMaxDCAToVertexXY + dcaToVertexZ*dcaToVertexZ/fMaxDCAToVertexZ/fMaxDCAToVertexZ);
	}
  else{
    dcaToVertex = TMath::Sqrt(dcaToVertexXY*dcaToVertexXY + dcaToVertexZ*dcaToVertexZ);
	}

	if (fUseDCAToVertex2D && dcaToVertex > 1)
    return kFALSE;
  if (!fUseDCAToVertex2D && TMath::Abs(dcaToVertexXY) > fMaxDCAToVertexXY)
    return kFALSE;
  if (!fUseDCAToVertex2D && TMath::Abs(dcaToVertexZ) > fMaxDCAToVertexZ)
    return kFALSE;

	return kTRUE;
}

