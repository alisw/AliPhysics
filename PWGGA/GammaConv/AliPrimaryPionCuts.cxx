
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
	fPIDResponse(NULL),
	fEsdTrackCuts(NULL),
	fEtaCut(0.9),
	fEtaShift(0.0),
	fDoEtaCut(kFALSE),
	fPtCut(0.0),
	fMinClsTPC(0), // minimum clusters in the TPC
    fChi2PerClsTPC(0), // maximum Chi2 per cluster in the TPC
    fRequireTPCRefit(kFALSE), // require a refit in the TPC
	fMinClsTPCToF(0), // minimum clusters to findable clusters
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
	fMassCut(10),
	fDoWeights(kFALSE),
    fMaxDCAToVertexZ(8000),
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
	fStringITSClusterCut("")
{
	InitPIDResponse();
	for(Int_t jj=0;jj<kNCuts;jj++){ fCuts[jj]=0; }
	fCutString=new TObjString((GetCutNumber()).Data());

	// Using standard function for setting Cuts
	Bool_t selectPrimaries=kFALSE;
	if (fEsdTrackCuts==NULL)fEsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
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

	if(preCut){
		fHistITSdEdxbefore=new TH2F(Form("Pion_ITS_before %s",cutName.Data()),"ITS dEdx pion before" ,150,0.05,20,400,-10,10);
		fHistograms->Add(fHistITSdEdxbefore);
		axisBeforeITS = fHistITSdEdxbefore->GetXaxis();

		fHistTPCdEdxbefore=new TH2F(Form("Pion_dEdx_before %s",cutName.Data()),"dEdx pion before" ,150,0.05,20,400,-10,10);
		fHistograms->Add(fHistTPCdEdxbefore);
		axisBeforedEdx = fHistTPCdEdxbefore->GetXaxis();

		fHistTPCdEdxSignalbefore=new TH2F(Form("Pion_dEdxSignal_before %s",cutName.Data()),"dEdx pion signal before" ,150,0.05,20.0,800,0.0,200);
		fHistograms->Add(fHistTPCdEdxSignalbefore);
		axisBeforedEdxSignal = fHistTPCdEdxSignalbefore->GetXaxis();

		fHistTOFbefore=new TH2F(Form("Pion_TOF_before %s",cutName.Data()),"TOF pion before" ,150,0.05,20,400,-6,10);
		fHistograms->Add(fHistTOFbefore);
		axisBeforeTOF = fHistTOFbefore->GetXaxis();
		
		fHistTrackDCAxyPtbefore = new TH2F(Form("hTrack_DCAxy_Pt_before %s",cutName.Data()),"DCAxy Vs Pt of tracks before",800,-4.0,4.0,400,0.,10.);
		fHistograms->Add(fHistTrackDCAxyPtbefore); 	
		
		fHistTrackDCAzPtbefore  = new TH2F(Form("hTrack_DCAz_Pt_before %s",cutName.Data()), "DCAz  Vs Pt of tracks before",800,-4.0,4.0,400,0.,10.);
		fHistograms->Add(fHistTrackDCAzPtbefore); 
		
		fHistTrackNFindClsPtTPCbefore = new TH2F(Form("hTrack_NFindCls_Pt_TPC_before %s",cutName.Data()),"Track: N Findable Cls TPC Vs Pt before",100,0,1,400,0.,10.);
		fHistograms->Add(fHistTrackNFindClsPtTPCbefore); 
	}

	fHistITSdEdxafter=new TH2F(Form("Pion_ITS_after %s",cutName.Data()),"ITS dEdx pion after" ,150,0.05,20,400, -10,10);
	fHistograms->Add(fHistITSdEdxafter);

	fHistTPCdEdxafter=new TH2F(Form("Pion_dEdx_after %s",cutName.Data()),"dEdx pion after" ,150,0.05,20,400, -10,10);
	fHistograms->Add(fHistTPCdEdxafter);

	fHistTPCdEdxSignalafter=new TH2F(Form("Pion_dEdxSignal_after %s",cutName.Data()),"dEdx pion signal after" ,150,0.05,20.0,800,0.0,200);
	fHistograms->Add(fHistTPCdEdxSignalafter);

	fHistTOFafter=new TH2F(Form("Pion_TOF_after %s",cutName.Data()),"TOF pion after" ,150,0.05,20,400,-6,10);
	fHistograms->Add(fHistTOFafter);
	
	fHistTrackDCAxyPtafter  = new TH2F(Form("hTrack_DCAxy_Pt_after %s",cutName.Data()),"DCAxy Vs Pt of tracks after",800,-4.0,4.0,400,0.,10.);
	fHistograms->Add(fHistTrackDCAxyPtafter); 
	
	fHistTrackDCAzPtafter  = new TH2F(Form("hTrack_DCAz_Pt_after %s",cutName.Data()), "DCAz Vs Pt of tracks  after",800,-4.0,4.0,400,0.,10.);
	fHistograms->Add(fHistTrackDCAzPtafter); 
	
	fHistTrackNFindClsPtTPCafter = new TH2F(Form("hTrack_NFindCls_Pt_TPC_after %s",cutName.Data()),"Track: N Findable Cls TPC Vs Pt after",100,0,1,400,0.,10.);
	fHistograms->Add(fHistTrackNFindClsPtTPCafter); 
	
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
	
	if(preCut){
		axisBeforeITS->Set(bins, newBins);
		axisBeforedEdx->Set(bins, newBins);
		axisBeforedEdxSignal->Set(bins,newBins);
		axisBeforeTOF->Set(bins, newBins);
	
	}
	delete [] newBins;
	
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

	
	if ( ! lTrack->GetConstrainedParam() ){
		return kFALSE;
	}
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
Bool_t AliPrimaryPionCuts::TrackIsSelected(AliESDtrack* lTrack) {
	// Track Selection for Photon Reconstruction
	
	Double_t clsToF = GetNFindableClustersTPC(lTrack);
	
	if( ! fEsdTrackCuts->AcceptTrack(lTrack) ){
		return kFALSE;
	}
		
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
Double_t AliPrimaryPionCuts::GetNFindableClustersTPC(AliESDtrack* lTrack){

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
			fStringITSClusterCut= "first or second SPD cluster required";
			break; //1 hit in any layer of SPD
		case 3: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
			fEsdTrackCuts->SetMinNClustersITS(4);
			fStringITSClusterCut= "first SPD cluster required, min number of ITS clusters = 4";
			// 4 hits in total in the ITS. At least 1 hit in the first layer of SPD  
			break;
		case 4: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(3);
			fStringITSClusterCut= "first or second SPD cluster required, min number of ITS clusters = 3";
			// 3 hits in total in the ITS. At least 1 hit in any layer of SPD
			break;
		case 5: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(4);
			fStringITSClusterCut= "first or second SPD cluster required, min number of ITS clusters = 4";
			// 4 hits in total in the ITS. At least 1 hit in any layer of SPD
			break;
		case 6: 
			fEsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
			fEsdTrackCuts->SetMinNClustersITS(5);
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
        case 10:
             fMinClsTPC     = 80.;
             fChi2PerClsTPC = 4;
             fRequireTPCRefit    = kTRUE;
             fEsdTrackCuts->SetMinNClustersTPC(fMinClsTPC);
             // Other Cuts concerning TPC
             fEsdTrackCuts->SetMaxChi2PerClusterTPC(fChi2PerClsTPC);
             fEsdTrackCuts->SetRequireTPCRefit(fRequireTPCRefit);
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
			fEsdTrackCuts->SetMaxDCAToVertexZ(1000);
			fEsdTrackCuts->SetMaxDCAToVertexXY(1000);
			fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
			break;
        case 1:
			fEsdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
			fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
			break;
		case 2: 
			fEsdTrackCuts->SetMaxDCAToVertexZ(2);
			fEsdTrackCuts->SetMaxDCAToVertexXY(1);
			fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
			break; 
        case 3:
            fMaxDCAToVertexZ = 3.0;
            fEsdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
            fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
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

