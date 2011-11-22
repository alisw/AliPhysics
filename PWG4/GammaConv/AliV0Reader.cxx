/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																																				*
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt												*
 * Version 1.0																														*
 *																																				*
 * Permission to use, copy, modify and distribute this software and its	 *
 * documentation strictly for non-commercial purposes is hereby granted	 *
 * without fee, provided that the above copyright notice appears in all	 *
 * copies and that both the copyright notice and this permission notice	 *
 * appear in the supporting documentation. The authors make no claims		 *
 * about the suitability of this software for any purpose. It is					*
 * provided "as is" without express or implied warranty.									*
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include <TMath.h>

//---- ANALYSIS system ----
#include "AliV0Reader.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliKFVertex.h"

#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliGammaConversionBGHandler.h"
#include "AliESDtrackCuts.h"
#include "TRandom3.h"

class iostream;
class AliESDv0;
class TFormula;
class TRandom3;

using namespace std;

ClassImp(AliV0Reader)


AliESDpid* AliV0Reader::fgESDpid = 0x0;

AliV0Reader::AliV0Reader() :
	TObject(),
	fMCStack(NULL),
	// fMCTruth(NULL),
	fMCEvent(NULL),		// for CF
	fChain(NULL),
	// fESDHandler(NULL),
	fESDEvent(NULL),
	fCFManager(NULL),
	//fESDpid(NULL),
	fHistograms(NULL),
	fCurrentV0IndexNumber(0),
	fCurrentV0(NULL),
	fCurrentNegativeKFParticle(NULL),
	fCurrentPositiveKFParticle(NULL),
	fCurrentMotherKFCandidate(NULL),
	fCurrentNegativeESDTrack(NULL),
	fCurrentPositiveESDTrack(NULL),
	fNegativeTrackLorentzVector(NULL),
	fPositiveTrackLorentzVector(NULL),
	fMotherCandidateLorentzVector(NULL),
	fCurrentXValue(0),
	fCurrentYValue(0),
	fCurrentZValue(0),
	fPositiveTrackPID(0),
	fNegativeTrackPID(0),
	fNegativeMCParticle(NULL),
	fPositiveMCParticle(NULL),
	fMotherMCParticle(NULL),
	fMotherCandidateKFMass(0),
	fMotherCandidateKFWidth(0),
	fUseKFParticle(kTRUE),
	fUseESDTrack(kFALSE),
	fDoMC(kFALSE),
	fMaxVertexZ(100.),// 100 cm(from the 0)
	fMaxR(10000),// 100 meter(outside of ALICE)
	fMinR(0),// 100 meter(outside of ALICE)
	fEtaCut(0.),
	fEtaCutMin(-0.1),
	fRapidityMesonCut(0.),
	fPtCut(0.),
	fSinglePtCut(0.),
	fMaxZ(0.),
	fMinClsTPC(0.),
	fMinClsTPCToF(0.),
	fLineCutZRSlope(0.),
	fLineCutZValue(0.),
	fLineCutZRSlopeMin(0.),
	fLineCutZValueMin(0.),
	fChi2CutConversion(0.),
	fChi2CutMeson(0.),
	fAlphaCutMeson(1.),
	fAlphaMinCutMeson(0.),
	fPIDProbabilityCutNegativeParticle(0),
	fPIDProbabilityCutPositiveParticle(0),
	fDodEdxSigmaCut(kFALSE),
	fDoTOFsigmaCut(kFALSE), // RRnewTOF
	fPIDnSigmaAboveElectronLine(100),
	fPIDnSigmaBelowElectronLine(-100),
	fTofPIDnSigmaAboveElectronLine(100), // RRnewTOF
	fTofPIDnSigmaBelowElectronLine(-100), // RRnewTOF
	fPIDnSigmaAbovePionLine(-100),
	fPIDnSigmaAbovePionLineHighPt(-100),
	fPIDMinPnSigmaAbovePionLine(100), 
	fPIDMaxPnSigmaAbovePionLine(100), 
	fDoKaonRejectionLowP(kFALSE),
	fDoProtonRejectionLowP(kFALSE),
	fDoPionRejectionLowP(kFALSE),
	fPIDnSigmaAtLowPAroundKaonLine(0),
	fPIDnSigmaAtLowPAroundProtonLine(0),
	fPIDnSigmaAtLowPAroundPionLine(0),
	fPIDMinPKaonRejectionLowP(0),
	fPIDMinPProtonRejectionLowP(0),
	fPIDMinPPionRejectionLowP(0),
	fDoQtGammaSelection(kFALSE),
	fDoHighPtQtGammaSelection(kFALSE), // RRnew
	fQtMax(100.),
	fHighPtQtMax(100.), // RRnew
	fPtBorderForQt(100.), // RRnew
	fXVertexCut(0.),
	fYVertexCut(0.),
	fZVertexCut(0.),
	fNSigmaMass(0.),
	fUseImprovedVertex(kFALSE),
	fUseOwnXYZCalculation(kFALSE),
	fUseConstructGamma(kFALSE),
	fDoCF(kFALSE),
	fUseEtaMinCut(kFALSE),
	fUseOnFlyV0Finder(kTRUE),
	fUpdateV0AlreadyCalled(kFALSE),
	fCurrentEventGoodV0s(NULL),
	fV0Pindex(),
	fV0Nindex(),
//	fPreviousEventGoodV0s(),
	fCalculateBackground(kFALSE),
	fBGEventHandler(NULL),
	fBGEventInitialized(kFALSE),
	fEsdTrackCuts(NULL),
	fNumberOfESDTracks(0),
	fNEventsForBGCalculation(20),
	fUseChargedTrackMultiplicityForBG(kTRUE),
	fNumberOfGoodV0s(0),
	fIsHeavyIon(0),
	fUseCorrectedTPCClsInfo(kFALSE),
	fUseMCPSmearing(kTRUE),
	fPBremSmearing(1.),
	fPSigSmearing(0.),
	fPSigSmearingCte(0.),
	fRandom(0),
	fBrem(NULL),
	fDoPhotonAsymmetryCut(0),
	fMinPPhotonAsymmetryCut(100.),
	fMinPhotonAsymmetry(0.)
{
	//fESDpid = new AliESDpid;	
}


AliV0Reader::AliV0Reader(const AliV0Reader & original) :
	TObject(original),
	fMCStack(original.fMCStack),
	// fMCTruth(original.fMCTruth),
	fMCEvent(original.fMCEvent),	// for CF
	fChain(original.fChain),
	//	fESDHandler(original.fESDHandler),
	fESDEvent(original.fESDEvent),
	fCFManager(original.fCFManager),
	// fESDpid(original.fESDpid),
	fHistograms(original.fHistograms),
	fCurrentV0IndexNumber(original.fCurrentV0IndexNumber),
	fCurrentV0(original.fCurrentV0),
	fCurrentNegativeKFParticle(original.fCurrentNegativeKFParticle),
	fCurrentPositiveKFParticle(original.fCurrentPositiveKFParticle),
	fCurrentMotherKFCandidate(original.fCurrentMotherKFCandidate),
	fCurrentNegativeESDTrack(original.fCurrentNegativeESDTrack),
	fCurrentPositiveESDTrack(original.fCurrentPositiveESDTrack),
	fNegativeTrackLorentzVector(original.fNegativeTrackLorentzVector),
	fPositiveTrackLorentzVector(original.fPositiveTrackLorentzVector),
	fMotherCandidateLorentzVector(original.fMotherCandidateLorentzVector),
	fCurrentXValue(original.fCurrentXValue),
	fCurrentYValue(original.fCurrentYValue),
	fCurrentZValue(original.fCurrentZValue),
	fPositiveTrackPID(original.fPositiveTrackPID),
	fNegativeTrackPID(original.fNegativeTrackPID),
	fNegativeMCParticle(original.fNegativeMCParticle),
	fPositiveMCParticle(original.fPositiveMCParticle),
	fMotherMCParticle(original.fMotherMCParticle),
	fMotherCandidateKFMass(original.fMotherCandidateKFMass),
	fMotherCandidateKFWidth(original.fMotherCandidateKFWidth),
	fUseKFParticle(kTRUE),
	fUseESDTrack(kFALSE),
	fDoMC(kFALSE),
	fMaxVertexZ(original.fMaxVertexZ),
	fMaxR(original.fMaxR),
	fMinR(original.fMinR),
	fEtaCut(original.fEtaCut),
	fEtaCutMin(original.fEtaCutMin),
	fRapidityMesonCut(original.fRapidityMesonCut),
	fPtCut(original.fPtCut),
	fSinglePtCut(original.fSinglePtCut),
	fMaxZ(original.fMaxZ),
	fMinClsTPC(original.fMinClsTPC),
	fMinClsTPCToF(original.fMinClsTPCToF),
	fLineCutZRSlope(original.fLineCutZRSlope),
	fLineCutZValue(original.fLineCutZValue),
	fLineCutZRSlopeMin(original.fLineCutZRSlopeMin),
	fLineCutZValueMin(original.fLineCutZValueMin),
	fChi2CutConversion(original.fChi2CutConversion),
	fChi2CutMeson(original.fChi2CutMeson),
	fAlphaCutMeson(original.fAlphaCutMeson),
	fAlphaMinCutMeson(original.fAlphaMinCutMeson),
	fPIDProbabilityCutNegativeParticle(original.fPIDProbabilityCutNegativeParticle),
	fPIDProbabilityCutPositiveParticle(original.fPIDProbabilityCutPositiveParticle),
	fDodEdxSigmaCut(original.fDodEdxSigmaCut),
	fDoTOFsigmaCut(original.fDoTOFsigmaCut), // RRnewTOF
	fPIDnSigmaAboveElectronLine(original.fPIDnSigmaAboveElectronLine),
	fPIDnSigmaBelowElectronLine(original.fPIDnSigmaBelowElectronLine),
	fTofPIDnSigmaAboveElectronLine(original.fTofPIDnSigmaAboveElectronLine), // RRnewTOF
	fTofPIDnSigmaBelowElectronLine(original.fTofPIDnSigmaBelowElectronLine), // RRnewTOF
	fPIDnSigmaAbovePionLine(original.fPIDnSigmaAbovePionLine), 
	fPIDnSigmaAbovePionLineHighPt(original.fPIDnSigmaAbovePionLineHighPt), 
	fPIDMinPnSigmaAbovePionLine(original.fPIDMinPnSigmaAbovePionLine), 
	fPIDMaxPnSigmaAbovePionLine(original.fPIDMaxPnSigmaAbovePionLine), 
	fDoKaonRejectionLowP(original.fDoKaonRejectionLowP),
	fDoProtonRejectionLowP(original.fDoProtonRejectionLowP),
	fDoPionRejectionLowP(original.fDoPionRejectionLowP),
	fPIDnSigmaAtLowPAroundKaonLine(original.fPIDnSigmaAtLowPAroundKaonLine),
	fPIDnSigmaAtLowPAroundProtonLine(original.fPIDnSigmaAtLowPAroundProtonLine),
	fPIDnSigmaAtLowPAroundPionLine(original.fPIDnSigmaAtLowPAroundPionLine),
	fPIDMinPKaonRejectionLowP(original.fPIDMinPKaonRejectionLowP),
	fPIDMinPProtonRejectionLowP(original.fPIDMinPProtonRejectionLowP),
	fPIDMinPPionRejectionLowP(original.fPIDMinPPionRejectionLowP),
	fDoQtGammaSelection(original.fDoQtGammaSelection),
	fDoHighPtQtGammaSelection(original.fDoHighPtQtGammaSelection), // RRnew
	fQtMax(original.fQtMax),
	fHighPtQtMax(original.fHighPtQtMax), // RRnew
	fPtBorderForQt(original.fPtBorderForQt), // RRnew
	fXVertexCut(original.fXVertexCut),
	fYVertexCut(original.fYVertexCut),
	fZVertexCut(original.fZVertexCut),
	fNSigmaMass(original.fNSigmaMass),
	fUseImprovedVertex(original.fUseImprovedVertex),
	fUseOwnXYZCalculation(original.fUseOwnXYZCalculation),
	fUseConstructGamma(original.fUseConstructGamma),
	fDoCF(original.fDoCF),
	fUseEtaMinCut(original.fUseEtaMinCut),
	fUseOnFlyV0Finder(original.fUseOnFlyV0Finder),
	fUpdateV0AlreadyCalled(original.fUpdateV0AlreadyCalled),
	fCurrentEventGoodV0s(original.fCurrentEventGoodV0s),
	fV0Pindex(original.fV0Pindex),
	fV0Nindex(original.fV0Nindex),
	//	fPreviousEventGoodV0s(original.fPreviousEventGoodV0s),
	fCalculateBackground(original.fCalculateBackground),
	fBGEventHandler(original.fBGEventHandler),
	fBGEventInitialized(original.fBGEventInitialized),
	fEsdTrackCuts(original.fEsdTrackCuts),
	fNumberOfESDTracks(original.fNumberOfESDTracks),
	fNEventsForBGCalculation(original.fNEventsForBGCalculation),
	fUseChargedTrackMultiplicityForBG(original.fUseChargedTrackMultiplicityForBG),
	fNumberOfGoodV0s(original.fNumberOfGoodV0s),
	fIsHeavyIon(original.fIsHeavyIon),
	fUseCorrectedTPCClsInfo(original.fUseCorrectedTPCClsInfo),
	fUseMCPSmearing(original.fUseMCPSmearing),
	fPBremSmearing(original.fPBremSmearing),
	fPSigSmearing(original.fPSigSmearing),
	fPSigSmearingCte(original.fPSigSmearingCte),
	fRandom(original.fRandom),
	fBrem(original.fBrem),
	fDoPhotonAsymmetryCut(original.fDoPhotonAsymmetryCut),
	fMinPPhotonAsymmetryCut(original.fMinPPhotonAsymmetryCut),
	fMinPhotonAsymmetry(original.fMinPhotonAsymmetry)
{
	
}


AliV0Reader & AliV0Reader::operator = (const AliV0Reader & /*source*/)
{
	// assignment operator
	return *this;
}
AliV0Reader::~AliV0Reader()
{
	//	if(fESDpid){
	// delete fESDpid;
	//}
}

//____________________________________________________________________________
void AliV0Reader::SetInputAndMCEvent(AliVEvent* esd, AliMCEvent* mc) {
	// Connect the data pointers

	SetInputEvent(esd);
	SetMC(mc);

}


void AliV0Reader::Initialize(){
	//see header file for documentation

	fUpdateV0AlreadyCalled = kFALSE;	

	/*
	// Get the input handler from the manager
	fESDHandler = (AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if(fESDHandler == NULL){
		//print warning here
	}

	// Get pointer to esd event from input handler
	fESDEvent = fESDHandler->GetEvent();
	if(fESDEvent == NULL){
		//print warning here
	}

	//Get pointer to MCTruth
	fMCTruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
	*/



	//	fMCTruth = mcH->MCEvent();
	//	fMC = mcH->MCEvent();
	// stack = fMC->Stack();


	//if(fMCTruth == NULL){
		//print warning here
	// fDoMC = kFALSE;
	//}

	if(fMCEvent == NULL){
		fDoMC = kFALSE;
	}

	//Get pointer to the mc stack
	//	if(fMCTruth){
	if(fMCEvent){
		fMCStack = fMCEvent->Stack();
		if(fMCStack == NULL){
			//print warning here
		}
		// Better parameters for MonteCarlo from A. Kalweit 2010/01/8
//		 fESDpid->GetTPCResponse().SetBetheBlochParameters( 2.15898e+00/50.,
// 							1.75295e+01,
// 							3.40030e-09,
// 							1.96178e+00,
// 							3.91720e+00);
	}
	else{
		// Better parameters for data from A. Kalweit 2010/01/8
 //		fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,
// 						 2.63394e+01,
// 						 5.04114e-11,
// 						 2.12543e+00,
// 						 4.88663e+00);
	}
	
	// for CF
	//Get pointer to the mc event
	if(fDoCF && fDoMC){
		//fMCEvent = fMCTruth->MCEvent();
		if(fMCEvent == NULL){
			//print warning here
			fDoCF = kFALSE;
		}	
	}
	
	fUseEtaMinCut = kFALSE;
	if ( fEtaCutMin != -0.1) {
		fUseEtaMinCut = kTRUE;
	}


	AliKFParticle::SetField(fESDEvent->GetMagneticField());

	//	fCurrentEventGoodV0s = new TClonesArray("TClonesArray", 0);
	if(fCurrentEventGoodV0s == NULL){
		fCurrentEventGoodV0s = new TClonesArray("AliKFParticle", 0);
	}

	fV0Pindex.clear();
	fV0Nindex.clear();

	if(gRandom != NULL){
		delete gRandom;
		gRandom= new TRandom3(0);
	}


	if (fBrem == NULL){
		fBrem = new TF1("fBrem","pow(-log(x),[0]/log(2.0)-1.0)/TMath::Gamma([0]/log(2.0))",0.00001,0.999999999);
		// tests done with 1.0e-14
		fBrem->SetParameter(0,fPBremSmearing);
		fBrem->SetNpx(100000);
	}

	if(fCalculateBackground == kTRUE){
		if(fBGEventInitialized == kFALSE){

			
			Double_t *zBinLimitsArray = new Double_t[9];
			zBinLimitsArray[0] = -50.00;
			zBinLimitsArray[1] = -3.375;
			zBinLimitsArray[2] = -1.605;
			zBinLimitsArray[3] = -0.225;
			zBinLimitsArray[4] = 1.065;
			zBinLimitsArray[5] = 2.445;
			zBinLimitsArray[6] = 4.245;
			zBinLimitsArray[7] = 50.00;
			zBinLimitsArray[8] = 1000.00;
			
			Double_t *multiplicityBinLimitsArray= new Double_t[6];
			if(fUseChargedTrackMultiplicityForBG == kTRUE){
				multiplicityBinLimitsArray[0] = 0;
				multiplicityBinLimitsArray[1] = 8.5;
				multiplicityBinLimitsArray[2] = 16.5;
				multiplicityBinLimitsArray[3] = 27.5;
				multiplicityBinLimitsArray[4] = 41.5;
				multiplicityBinLimitsArray[5] = 100.;
				if(fIsHeavyIon){
					multiplicityBinLimitsArray[0] = 0;
					multiplicityBinLimitsArray[1] = 200.;
					multiplicityBinLimitsArray[2] = 500.;
					multiplicityBinLimitsArray[3] = 1000.;
					multiplicityBinLimitsArray[4] = 1500.;
					multiplicityBinLimitsArray[5] = 3000.;
				}
				fBGEventHandler = new AliGammaConversionBGHandler(9,6,fNEventsForBGCalculation);
			} else {
				multiplicityBinLimitsArray[0] = 2;
				multiplicityBinLimitsArray[1] = 3;
				multiplicityBinLimitsArray[2] = 4;
				multiplicityBinLimitsArray[3] = 5;
				multiplicityBinLimitsArray[4] = 9999;
				if(fIsHeavyIon){
					multiplicityBinLimitsArray[0] = 2;
					multiplicityBinLimitsArray[1] = 10;
					multiplicityBinLimitsArray[2] = 30;
					multiplicityBinLimitsArray[3] = 50;
					multiplicityBinLimitsArray[4] = 9999;
				}

				fBGEventHandler = new AliGammaConversionBGHandler(9,5,fNEventsForBGCalculation);
			}


			
			/*
			// ---------------------------------
			Double_t *zBinLimitsArray = new Double_t[1];
			zBinLimitsArray[0] = 999999.00;

			Double_t *multiplicityBinLimitsArray= new Double_t[1];
			multiplicityBinLimitsArray[0] = 99999999.00;
			fBGEventHandler = new AliGammaConversionBGHandler(1,1,10);
			// ---------------------------------
			*/
			fBGEventHandler->Initialize(zBinLimitsArray, multiplicityBinLimitsArray);
			fBGEventInitialized = kTRUE;
		}
	}
}

AliESDv0* AliV0Reader::GetV0(Int_t index){
	//see header file for documentation
	fCurrentV0 = fESDEvent->GetV0(index);
	UpdateV0Information();
	return fCurrentV0;
}

Int_t AliV0Reader::GetNumberOfContributorsVtx(){
	// returns number of contributors to the vertex
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
			cout<<"number of contributors from bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
			return 0;
		}
	}
	return 0;
}

Bool_t AliV0Reader::CheckForPrimaryVertex(){
	//see headerfile for documentation
	if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
		return 1;
	}
	if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()<1) {
	// SPD vertex
		if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
			// return 0;
			//-AM test pi0s without SPD only vertex
			//cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
			return 1;
		}
		if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
			//			cout<<"bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
			return 0;
		}
	}
	return 0;
	//	return fESDEvent->GetPrimaryVertex()->GetNContributors()>0;
}

Bool_t AliV0Reader::CheckForPrimaryVertexZ(){
	//see headerfile for documentation
	if(TMath::Abs(fESDEvent->GetPrimaryVertex()->GetZ())<GetMaxVertexZ()){
		return kTRUE;
	}else{
		return kFALSE;
	}
	return kTRUE;
}

Bool_t AliV0Reader::CheckV0FinderStatus(Int_t index){
	// see headerfile for documentation
	if(fUseOnFlyV0Finder){
		if(!GetV0(index)->GetOnFlyStatus()){
			return kFALSE;
		}
	}
	if(!fUseOnFlyV0Finder){
		if(GetV0(index)->GetOnFlyStatus()){
			return kFALSE;
		}
	}
	return kTRUE;
}



Bool_t AliV0Reader::NextV0(){
	//see header file for documentation
	Bool_t iResult=kFALSE;
	while(fCurrentV0IndexNumber<fESDEvent->GetNumberOfV0s()){
		fCurrentV0 = fESDEvent->GetV0(fCurrentV0IndexNumber);

		fUpdateV0AlreadyCalled=kFALSE;

		if(fHistograms != NULL){
			fHistograms->FillHistogram("ESD_AllV0s_InvMass",GetMotherCandidateMass());
		}
		
		// moved it up here so that the correction framework can access pt and eta information
		if(UpdateV0Information() == kFALSE){
			fCurrentV0IndexNumber++;
			continue;
		}
 
		Double_t containerInput[3];
		if(fDoCF){
			containerInput[0] = GetMotherCandidatePt();
			containerInput[1] = GetMotherCandidateEta();
			containerInput[2] = GetMotherCandidateMass();
		}
		/*
		if(fDoCF){
			containerInput[0] = GetMotherCandidatePt();
			containerInput[1] = GetMotherCandidateEta();
			containerInput[2] = GetMotherCandidateMass();
			
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepLikeSign);		// for CF	
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCRefit);		// for CF	
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepKinks);		// for CF	
		}
		*/

		//checks if on the fly mode is set
		if ( !CheckV0FinderStatus(fCurrentV0IndexNumber) ){
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutGetOnFly_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepGetOnFly);		// for CF	
		}

		if(fHistograms != NULL){
			fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_InvMass",GetMotherCandidateMass());
		}
 
		Double_t armenterosQtAlfa[2];
		Double_t armenterosQtAlphaESDMC[2];
		Double_t armenterosQtAlphaMC[2];
		GetArmenterosQtAlfa(GetNegativeKFParticle(), 
			GetPositiveKFParticle(), 
			GetMotherCandidateKFCombination(),
			armenterosQtAlfa);
	 
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_alfa_qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
 
	 
		if(fCurrentNegativeESDTrack->Charge() == fCurrentPositiveESDTrack->Charge()){						 // avoid like sign
			//	iResult=kFALSE;
			if(fHistograms != NULL ){
				fHistograms->FillHistogram("ESD_CutLikeSign_InvMass",GetMotherCandidateMass());
				// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
				// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
				//	fUpdateV0AlreadyCalled = kTRUE;
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepLikeSign);		// for CF	
		}
 	
		if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kTPCrefit) || 
			!(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kTPCrefit) ){
			//	if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kITSrefit) || 
			//			!(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kITSrefit) ){
			//	iResult=kFALSE;
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutRefit_InvMass",GetMotherCandidateMass());
				// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
				// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
				//fUpdateV0AlreadyCalled = kTRUE;
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCRefit);		// for CF	
		}
 	


		if( fCurrentNegativeESDTrack->GetKinkIndex(0) > 0 || 
			fCurrentPositiveESDTrack->GetKinkIndex(0) > 0) {			
			//iResult=kFALSE;
			if(fHistograms != NULL ){
				fHistograms->FillHistogram("ESD_CutKink_InvMass",GetMotherCandidateMass());
				// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
				// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
				//fUpdateV0AlreadyCalled = kTRUE;
			}
			fCurrentV0IndexNumber++;
			continue;
		}

		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepKinks);		// for CF	
		}
 	
		if(GetXYRadius()>fMaxR){ // cuts on distance from collision point
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutR_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}	
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepR);			// for CF
		}
		if(GetXYRadius()<fMinR){ // cuts on distance from collision point
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutMinR_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}
				
		//if((TMath::Abs(fCurrentZValue)*fLineCutZRSlope)-fLineCutZValue > GetXYRadius() ) { // cuts out regions where we do not reconstruct
		if( GetXYRadius() <= ((TMath::Abs(fCurrentZValue)*fLineCutZRSlope)-fLineCutZValue)){
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutLine_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		} else if (fUseEtaMinCut &&  GetXYRadius() >= ((TMath::Abs(fCurrentZValue)*fLineCutZRSlopeMin)-fLineCutZValueMin )){ 
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutLine_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}

		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepLine);			// for CF
		}
		
		if(TMath::Abs(fCurrentZValue) > fMaxZ ){ // cuts out regions where we do not reconstruct
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutZ_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepZ);		// for CF	
		}

		if(fUseKFParticle){
			if(TMath::Abs(fMotherCandidateLorentzVector->Eta())> fEtaCut || TMath::Abs(fMotherCandidateLorentzVector->Eta())< fEtaCutMin){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutEta_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}

			if(TMath::Abs(fCurrentNegativeKFParticle->GetEta())> fEtaCut || TMath::Abs(fCurrentNegativeKFParticle->GetEta())< fEtaCutMin){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutEta_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}

			if(TMath::Abs(fCurrentPositiveKFParticle->GetEta())> fEtaCut || TMath::Abs(fCurrentPositiveKFParticle->GetEta())< fEtaCutMin){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutEta_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}
		}

fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_Pt_Qt",GetMotherCandidatePt(),armenterosQtAlfa[0]);

		if(fDoMC){
			if ( HasSameMCMother() == kTRUE){ 
				GetArmenterosQtAlfa(fNegativeMCParticle, 
					fPositiveMCParticle, 
					fMotherMCParticle,
					armenterosQtAlphaMC);
			}
			GetArmenterosQtAlfa(fNegativeMCParticle, 
				fPositiveMCParticle, 
				GetMotherCandidateKFCombination(),
				armenterosQtAlphaESDMC );
		}
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_goodtracks_alfa_qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
		if( fCurrentNegativeKFParticle->GetPt()> 0.150 &&	fCurrentPositiveKFParticle->GetPt()> 0.150){
			fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_minPt_GT_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
		}
		if(fDoMC){
			fHistograms->FillHistogram("ESD_TrueConvAllV0s_ESDMother_Alpha_Qt",armenterosQtAlphaESDMC[1],armenterosQtAlphaESDMC[0]);
			if ( HasSameMCMother() == kTRUE){ 
				fHistograms->FillHistogram("ESD_TrueConvSameMother_ESDMother_Alpha_Qt",armenterosQtAlphaESDMC[1],armenterosQtAlphaESDMC[0]);
				fHistograms->FillHistogram("ESD_TrueConvSameMother_MCMother_Alpha_Qt",armenterosQtAlphaMC[1],armenterosQtAlphaMC[0]);
				if (fMotherMCParticle->GetPdgCode() == 22 ){
					fHistograms->FillHistogram("ESD_TrueConvGamma_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
					fHistograms->FillHistogram("ESD_TrueConvGamma_Pt_Qt",GetMotherCandidatePt(),armenterosQtAlfa[0]);	
				} else if ( fMotherMCParticle->GetPdgCode() == 310 ){
					fHistograms->FillHistogram("ESD_TrueConvK0s_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				} else if ( fMotherMCParticle->GetPdgCode() == 113 ){
					fHistograms->FillHistogram("ESD_TrueConvRho0_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				} else if ( fMotherMCParticle->GetPdgCode() == 333 ){
					fHistograms->FillHistogram("ESD_TrueConvPhi_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				} else if ( (fMotherMCParticle->GetPdgCode() == 3122 || fMotherMCParticle->GetPdgCode() == -3122) ){
					fHistograms->FillHistogram("ESD_TrueConvLambda_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				} else if ( (fMotherMCParticle->GetPdgCode() == 2114 || fMotherMCParticle->GetPdgCode() == -2114) ){
					fHistograms->FillHistogram("ESD_TrueConvDelta_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				} else if ( (fMotherMCParticle->GetPdgCode() == 313 || 
								fMotherMCParticle->GetPdgCode() == 323 || 
								fMotherMCParticle->GetPdgCode() == -323 ) ){
					fHistograms->FillHistogram("ESD_TrueConvKStar_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
				} else {
					fHistograms->FillHistogram("ESD_TrueConvUnknown_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
					fHistograms->FillHistogram("ESD_TrueConvUnknown_Qt_PDG",fMotherMCParticle->GetPdgCode());
				//	cout << "unidentfied mother: pdg-C mother " << fMotherMCParticle->GetPdgCode() << " daughters " << fNegativeMCParticle->GetPdgCode() << "\t" << fPositiveMCParticle->GetPdgCode() << endl;
				}
			}	else {
				fHistograms->FillHistogram("ESD_TrueConvComb_Alpha_Qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
			}
		}

		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_E_dEdxP",fCurrentNegativeESDTrack->P(),fCurrentNegativeESDTrack->GetTPCsignal());
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_P_dEdxP",fCurrentPositiveESDTrack->P(),fCurrentPositiveESDTrack->GetTPCsignal());
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_E_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_P_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_PiPl_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_PiMi_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_KPl_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kKaon));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_KMi_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kKaon));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_PPl_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kProton));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_PMi_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kProton));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_MuPl_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kMuon));
		fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_MuMi_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kMuon));
 
		if(fDodEdxSigmaCut == kTRUE){
			if( fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
			fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine ||
			fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
			fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine ){
				//iResult=kFALSE;
				if(fHistograms != NULL ){
					fHistograms->FillHistogram("ESD_CutdEdxSigmaElectronLine_InvMass",GetMotherCandidateMass());
					// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
					// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
					//fUpdateV0AlreadyCalled = kTRUE;
				}
				fCurrentV0IndexNumber++;
				continue;
			}
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepdEdxElectronselection);							 // for CF
			}

			fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_PiPl_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion));
			fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_PiMi_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion));
			fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_MuPl_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kMuon));
			fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_MuMi_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kMuon));

			
			if( fCurrentPositiveESDTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentPositiveESDTrack->P()<fPIDMaxPnSigmaAbovePionLine ){
				if(fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
					fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
					fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
					//		iResult=kFALSE;
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutdEdxSigmaPionLine_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
			
			if( fCurrentNegativeESDTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentNegativeESDTrack->P()<fPIDMaxPnSigmaAbovePionLine){
				if(fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
					fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
					fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
					//		iResult=kFALSE;
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutdEdxSigmaPionLine_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}

			// High Pt
			if( fCurrentPositiveESDTrack->P()>fPIDMaxPnSigmaAbovePionLine ){
				if(fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
					fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
					fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){
					//		iResult=kFALSE;
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutdEdxSigmaPionLine_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
			
			if( fCurrentNegativeESDTrack->P()>fPIDMaxPnSigmaAbovePionLine){
				if(fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
					fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
					fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){
			//		iResult=kFALSE;
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutdEdxSigmaPionLine_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}


			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepdEdxPionrejection);							 // for CF
			}

		}

		fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_KPl_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kKaon));
		fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_KMi_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kKaon));
		
		if(fDoKaonRejectionLowP == kTRUE){
			if( fCurrentNegativeESDTrack->P()<fPIDMinPKaonRejectionLowP ){
				if( TMath::Abs(fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutKaonRejectionLowP_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
			if( fCurrentPositiveESDTrack->P()<fPIDMinPKaonRejectionLowP ){
				if( TMath::Abs(fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutKaonRejectionLowP_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
		}

		fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_PPl_SigdEdxP",fCurrentPositiveESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kProton));
		fHistograms->FillHistogram("ESD_ConvGammaBeforeCorresCut_PMi_SigdEdxP",fCurrentNegativeESDTrack->P(),fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kKaon));
		
		if(fDoProtonRejectionLowP == kTRUE){
			if( fCurrentNegativeESDTrack->P()<fPIDMinPProtonRejectionLowP){
				if( TMath::Abs(fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutProtonRejectionLowP_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
			if( fCurrentPositiveESDTrack->P()<fPIDMinPProtonRejectionLowP ){
				if( TMath::Abs(fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutProtonRejectionLowP_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
		}

		if(fDoPionRejectionLowP == kTRUE){
			if( fCurrentNegativeESDTrack->P()<fPIDMinPPionRejectionLowP ){
				if( TMath::Abs(fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutPionRejectionLowP_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
			if( fCurrentPositiveESDTrack->P()<fPIDMinPPionRejectionLowP ){
				if( TMath::Abs(fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutPionRejectionLowP_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
		}


		if( fDoTOFsigmaCut == kTRUE ){ // RRnewTOF start ///////////////////////////////////////////////////////////////////////////// 
			Bool_t PosTrackNotTOFelec = kFALSE;
			Bool_t NegTrackNotTOFelec = kFALSE;
			if( (fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kTOFpid) && !(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kTOFmismatch) ){
				Double_t t0pos = fgESDpid->GetTOFResponse().GetStartTime(fCurrentPositiveESDTrack->P());
				Double_t fnSigmaPos = fgESDpid->NumberOfSigmasTOF(fCurrentPositiveESDTrack, AliPID::kElectron, t0pos);
				if( (fnSigmaPos>fTofPIDnSigmaAboveElectronLine) || (fnSigmaPos<fTofPIDnSigmaBelowElectronLine) ) PosTrackNotTOFelec = kTRUE;
			}
			if( (fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kTOFpid) && !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kTOFmismatch) ){
				Double_t t0neg = fgESDpid->GetTOFResponse().GetStartTime(fCurrentNegativeESDTrack->P());
				Double_t fnSigmaNeg = fgESDpid->NumberOfSigmasTOF(fCurrentNegativeESDTrack, AliPID::kElectron, t0neg);
				if( (fnSigmaNeg>fTofPIDnSigmaAboveElectronLine) || (fnSigmaNeg<fTofPIDnSigmaBelowElectronLine) ) NegTrackNotTOFelec = kTRUE;	
			}
			if( (PosTrackNotTOFelec==kTRUE) || (NegTrackNotTOFelec==kTRUE) ){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutTOFsigmaElec_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}
		} /////////////////////////////// RRnewTOF end ///////////////////////////////////////////////////////////////////////////////


		// Gamma selection based on QT from Armenteros
		if(fDoQtGammaSelection == kTRUE){ // RRnew start : apply different qT-cut above/below
			if(fDoHighPtQtGammaSelection){
				if(GetMotherCandidatePt() < fPtBorderForQt){
					if(armenterosQtAlfa[0]>fQtMax){
						if(fHistograms != NULL){
							fHistograms->FillHistogram("ESD_CutQt_InvMass",GetMotherCandidateMass());
						}
						fCurrentV0IndexNumber++;
						continue;
					}
				} else {
					if(armenterosQtAlfa[0]>fHighPtQtMax)	{
						if(fHistograms != NULL){
							fHistograms->FillHistogram("ESD_CutQt_InvMass",GetMotherCandidateMass());
						}
						fCurrentV0IndexNumber++;
						continue;
					}
				}
			} else {
				if(armenterosQtAlfa[0]>fQtMax){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutQt_InvMass",GetMotherCandidateMass());
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
		}	 // RRnew end

		if(fDoPhotonAsymmetryCut == kTRUE){
			if( fNegativeTrackLorentzVector->P()>fMinPPhotonAsymmetryCut ){
				Double_t trackNegAsy=0;
				if (fCurrentMotherKFCandidate->GetP()!=0.){
					trackNegAsy= fNegativeTrackLorentzVector->P()/fMotherCandidateLorentzVector->P();
				}
				if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutPhotonAsymmetry_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}

			if( fPositiveTrackLorentzVector->P()>fMinPPhotonAsymmetryCut ){
				Double_t trackPosAsy=0;
				if (fCurrentMotherKFCandidate->GetP()!=0.){
					trackPosAsy= fPositiveTrackLorentzVector->P()/fMotherCandidateLorentzVector->P();
				}
				if( trackPosAsy<fMinPhotonAsymmetry ||trackPosAsy>(1.- fMinPhotonAsymmetry)){
					if(fHistograms != NULL){
						fHistograms->FillHistogram("ESD_CutPhotonAsymmetry_InvMass",GetMotherCandidateMass());
						// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
						// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
						//fUpdateV0AlreadyCalled = kTRUE;
					}
					fCurrentV0IndexNumber++;
					continue;
				}
			}
		}
		//checks if we have a prim vertex
		//if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) { 
		if(GetNumberOfContributorsVtx()<=0) { 
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutNContributors_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepNContributors);		// for CF	
		}
		
		//Check the pid probability
		if(CheckPIDProbability(fPIDProbabilityCutNegativeParticle,fPIDProbabilityCutPositiveParticle)==kFALSE){
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutPIDProb_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCPID);			// for CF
		}
		
		
		/* Moved further up so corr framework can work
			 if(UpdateV0Information() == kFALSE){
			 fCurrentV0IndexNumber++;
			 continue;
			 }
		*/
		if(fCurrentNegativeESDTrack->GetNcls(1) < fMinClsTPC ||	fCurrentPositiveESDTrack->GetNcls(1) < fMinClsTPC ){
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutMinNClsTPC_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}
		if(fDoCF){
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepMinClsTPC);		// for CF	
		}
		Double_t negclsToF = 0.;
		if (!fUseCorrectedTPCClsInfo ){
			if(fCurrentNegativeESDTrack->GetTPCNclsF()!=0	){
				negclsToF = (Double_t)fCurrentNegativeESDTrack->GetNcls(1)/(Double_t)fCurrentNegativeESDTrack->GetTPCNclsF();
			}
		} else {
			negclsToF = fCurrentNegativeESDTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(GetXYRadius()));
		}

		Double_t posclsToF = 0.;
		if (!fUseCorrectedTPCClsInfo ){
			if(fCurrentPositiveESDTrack->GetTPCNclsF()!=0	){
				posclsToF = (Double_t)fCurrentPositiveESDTrack->GetNcls(1)/(Double_t)fCurrentPositiveESDTrack->GetTPCNclsF();
			}
		}else{
			posclsToF = fCurrentPositiveESDTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(GetXYRadius()));
		}

		if( negclsToF < fMinClsTPCToF ||	posclsToF < fMinClsTPCToF ){
			if(fHistograms != NULL){
				fHistograms->FillHistogram("ESD_CutMinNClsTPCToF_InvMass",GetMotherCandidateMass());
			}
			fCurrentV0IndexNumber++;
			continue;
		}



		
		if(fUseKFParticle){


			if( fCurrentNegativeKFParticle->GetPt()< fSinglePtCut ||	fCurrentPositiveKFParticle->GetPt()< fSinglePtCut){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutSinglePt_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepSinglePt);		// for CF	
			}


			if(fCurrentMotherKFCandidate->GetNDF()<=0){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutNDF_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepNDF);		// for CF	
			}
			
			Double_t chi2V0 = fCurrentMotherKFCandidate->GetChi2()/fCurrentMotherKFCandidate->GetNDF();
			if(chi2V0 > fChi2CutConversion || chi2V0 <=0){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutChi2_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepChi2);			// for CF
			}
		
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepEta);			// for CF
			}
			
			if(fMotherCandidateLorentzVector->Pt()<fPtCut){
				if(fHistograms != NULL){
					fHistograms->FillHistogram("ESD_CutPt_InvMass",GetMotherCandidateMass());
				}
				fCurrentV0IndexNumber++;
				continue;
			}
			if(fDoCF){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepPt);			// for CF
			}
			
		}
		else if(fUseESDTrack){
			//TODO
		}

		if(fHistograms != NULL){
			fHistograms->FillHistogram("ESD_GoodV0s_InvMass",GetMotherCandidateMass());
		}

		//		fCurrentEventGoodV0s.push_back(*fCurrentMotherKFCandidate);

		if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
			fCurrentMotherKFCandidate->E()=fCurrentMotherKFCandidate->GetP();
		}

		if(fDoMC&& fUseMCPSmearing>0){
			SmearKFParticle(fCurrentMotherKFCandidate);

		}		

		new((*fCurrentEventGoodV0s)[fCurrentEventGoodV0s->GetEntriesFast()])	AliKFParticle(*fCurrentMotherKFCandidate);
		fV0Pindex.push_back(fCurrentV0->GetPindex());
		fV0Nindex.push_back(fCurrentV0->GetNindex());

		iResult=kTRUE;//means we have a v0 who survived all the cuts applied
		
		fNumberOfGoodV0s++;

		fCurrentV0IndexNumber++;
		
		break;
	}
	return iResult; 
}

Bool_t AliV0Reader::UpdateV0Information(){
	//see header file for documentation
	
	Bool_t iResult=kTRUE;		 				// for taking out not refitted, kinks and like sign tracks 
	
	Bool_t switchTracks = kFALSE;
	
	fCurrentNegativeESDTrack = fESDEvent->GetTrack(fCurrentV0->GetNindex());
	fCurrentPositiveESDTrack = fESDEvent->GetTrack(fCurrentV0->GetPindex());


	if(fCurrentPositiveESDTrack->GetSign() == -1 && fCurrentNegativeESDTrack->GetSign() == 1){	// switch wrong signed tracks
		fCurrentNegativeESDTrack = fESDEvent->GetTrack(fCurrentV0->GetPindex());
		fCurrentPositiveESDTrack = fESDEvent->GetTrack(fCurrentV0->GetNindex());
		switchTracks = kTRUE;
	}

	
	if(fCurrentNegativeKFParticle != NULL){
		delete fCurrentNegativeKFParticle;
	}
	if(switchTracks == kFALSE){
		fCurrentNegativeKFParticle = new AliKFParticle(*(fCurrentV0->GetParamN()),fNegativeTrackPID);
	}
	else{
		fCurrentNegativeKFParticle = new AliKFParticle(*(fCurrentV0->GetParamP()),fNegativeTrackPID);
	}
	
	if(fCurrentPositiveKFParticle != NULL){
		delete fCurrentPositiveKFParticle;
	}
	if(switchTracks == kFALSE){
		fCurrentPositiveKFParticle = new AliKFParticle(*(fCurrentV0->GetParamP()),fPositiveTrackPID);
	}
	else{
		fCurrentPositiveKFParticle = new AliKFParticle(*(fCurrentV0->GetParamN()),fPositiveTrackPID);
	}
		
	if(fCurrentMotherKFCandidate != NULL){
		delete fCurrentMotherKFCandidate;
	}

	if(fUseConstructGamma==kTRUE){
		fCurrentMotherKFCandidate = new AliKFParticle;//(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
		fCurrentMotherKFCandidate->ConstructGamma(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
	}else{
		fCurrentMotherKFCandidate = new AliKFParticle(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
		if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
			fCurrentMotherKFCandidate->SetMassConstraint(0,fNSigmaMass);
		}
	}
	if(fUseImprovedVertex == kTRUE){
		AliKFVertex primaryVertexImproved(*GetPrimaryVertex());
		primaryVertexImproved+=*fCurrentMotherKFCandidate;
		fCurrentMotherKFCandidate->SetProductionVertex(primaryVertexImproved);
	}
	
	fCurrentMotherKFCandidate->GetMass(fMotherCandidateKFMass,fMotherCandidateKFWidth);
		
	if(fNegativeTrackLorentzVector != NULL){
		delete fNegativeTrackLorentzVector;
	}
	if(fUseKFParticle){
		fNegativeTrackLorentzVector = new TLorentzVector(fCurrentNegativeKFParticle->Px(),fCurrentNegativeKFParticle->Py(),fCurrentNegativeKFParticle->Pz());
	}
	else	{ //if(fUseESDTrack){
		fNegativeTrackLorentzVector = new TLorentzVector(fCurrentNegativeESDTrack->Px(),fCurrentNegativeESDTrack->Py(),fCurrentNegativeESDTrack->Pz());
	}
	
	if(fPositiveTrackLorentzVector != NULL){
		delete fPositiveTrackLorentzVector;
	}
	if(fUseKFParticle){
		fPositiveTrackLorentzVector = new TLorentzVector(fCurrentPositiveKFParticle->Px(),fCurrentPositiveKFParticle->Py(),fCurrentPositiveKFParticle->Pz());
	}
	else	{ // if(fUseESDTrack){	fPositiveTrackLorentzVector must be reinitialized, so assuming use ESD if not kfparticle Svein. 
		fPositiveTrackLorentzVector = new TLorentzVector(fCurrentPositiveESDTrack->Px(),fCurrentPositiveESDTrack->Py(),fCurrentPositiveESDTrack->Pz());
	}
	
	if(fMotherCandidateLorentzVector != NULL){
		delete fMotherCandidateLorentzVector;
	}

	fMotherCandidateLorentzVector = new TLorentzVector(*fNegativeTrackLorentzVector + *fPositiveTrackLorentzVector);
	
	if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
		fMotherCandidateLorentzVector->SetXYZM(fMotherCandidateLorentzVector->Px() ,fMotherCandidateLorentzVector->Py(),fMotherCandidateLorentzVector->Pz(),0.); 
	}
		
	
	if(fDoMC == kTRUE){
		fMotherMCParticle= NULL;
		if(switchTracks == kFALSE){
			fNegativeMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetNindex())->GetLabel()));
			fPositiveMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetPindex())->GetLabel()));
		}else{
			fNegativeMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetPindex())->GetLabel()));
			fPositiveMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetNindex())->GetLabel()));
		}

		if(fPositiveMCParticle->GetMother(0)>-1){
			fMotherMCParticle = fMCStack->Particle(fPositiveMCParticle->GetMother(0));
		}
	}



	

	// for CF
//	 Double_t containerInput[3];
//	 if(fDoCF){
//		 containerInput[0] = GetMotherCandidatePt();
//		 containerInput[1] = GetMotherCandidateEta();
//		 containerInput[2] = GetMotherCandidateMass();
		
//		 fCFManager->GetParticleContainer()->Fill(containerInput,kStepLikeSign);		// for CF	
//		 fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCRefit);		// for CF	
//		 fCFManager->GetParticleContainer()->Fill(containerInput,kStepKinks);		// for CF	
//	 }
	

	if(fUseOwnXYZCalculation == kFALSE){
		if(fUseConstructGamma == kFALSE){
			fCurrentV0->GetXYZ(fCurrentXValue,fCurrentYValue,fCurrentZValue);
		}else{
			fCurrentXValue=GetMotherCandidateKFCombination()->GetX();
			fCurrentYValue=GetMotherCandidateKFCombination()->GetY();
			fCurrentZValue=GetMotherCandidateKFCombination()->GetZ();
		}
	}
	else{
		Double_t convpos[2]; //Double_t convpos[3];
		convpos[0]=0;
		convpos[1]=0;
//     convpos[2]=0;

		GetConvPosXY(GetPositiveESDTrack(),GetNegativeESDTrack(),GetMagneticField(),convpos);
		
		/*if(switchTracks == kFALSE){
			GetConversionPoint((fCurrentV0->GetParamP()),(fCurrentV0->GetParamN()),convpos);
		}else{
			GetConversionPoint((fCurrentV0->GetParamN()),(fCurrentV0->GetParamP()),convpos);
		}*/
		
		fCurrentXValue = convpos[0];
		fCurrentYValue = convpos[1];
// 		fCurrentZValue = convpos[2];
		fCurrentZValue = GetConvPosZ(GetPositiveESDTrack(),GetNegativeESDTrack(),GetMagneticField());
	}
	/*
	if(fCurrentNegativeESDTrack->GetSign() == fCurrentPositiveESDTrack->GetSign()){						 // avoid like sign
		iResult=kFALSE;
		if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE && doFillHistos == kTRUE){
			fHistograms->FillHistogram("ESD_CutLikeSign_InvMass",GetMotherCandidateMass());
			// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
			// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
			fUpdateV0AlreadyCalled = kTRUE;
		}
	}
	
	if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kTPCrefit) || 
			!(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kTPCrefit) ){
		//	if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kITSrefit) || 
		//			!(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kITSrefit) ){
		iResult=kFALSE;
		if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE	&& doFillHistos == kTRUE){
			fHistograms->FillHistogram("ESD_CutRefit_InvMass",GetMotherCandidateMass());
			// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
			// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
			fUpdateV0AlreadyCalled = kTRUE;
		}
	}
	
	if( fCurrentNegativeESDTrack->GetKinkIndex(0) > 0 || 
			fCurrentPositiveESDTrack->GetKinkIndex(0) > 0) {			
		
		iResult=kFALSE;
		if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE && doFillHistos == kTRUE ){
			fHistograms->FillHistogram("ESD_CutKink_InvMass",GetMotherCandidateMass());
			// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
			// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
			fUpdateV0AlreadyCalled = kTRUE;
		}
	}

	if(fDodEdxSigmaCut == kTRUE){

		if( fESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
	fESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine ||
	fESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
	fESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine ){
			iResult=kFALSE;
			if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE	&& doFillHistos == kTRUE){
	fHistograms->FillHistogram("ESD_CutdEdxSigmaElectronLine_InvMass",GetMotherCandidateMass());
	// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
	// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
	fUpdateV0AlreadyCalled = kTRUE;
			}
		}
		if( fCurrentPositiveESDTrack->P()>fPIDMinPnSigmaAbovePionLine){
			if(fESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
	 fESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
	 fESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
	iResult=kFALSE;
	if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE	&& doFillHistos == kTRUE){
		fHistograms->FillHistogram("ESD_CutdEdxSigmaPionLine_InvMass",GetMotherCandidateMass());
		// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
		// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
		fUpdateV0AlreadyCalled = kTRUE;
	}
			}
		}

		if( fCurrentNegativeESDTrack->P()>fPIDMinPnSigmaAbovePionLine){
			if(fESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
	 fESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
	 fESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
	iResult=kFALSE;
	if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE && doFillHistos == kTRUE ){
		fHistograms->FillHistogram("ESD_CutdEdxSigmaPionLine_InvMass",GetMotherCandidateMass());
		// to avoid filling the other cut histograms. So in this case fUpdateV0AlreadyCalled also serves as a flag for the histogram filling
		// it will anyway be set to true at the end of the UpdateV0Information function, and there are no return until the end
		fUpdateV0AlreadyCalled = kTRUE;
	}
			}
		}
	}
	*/
	fUpdateV0AlreadyCalled = kTRUE;

	return iResult;
}



Bool_t AliV0Reader::HasSameMCMother(){
	//see header file for documentation
	
	Bool_t iResult = kFALSE;
	if(fDoMC == kTRUE){
		if(fNegativeMCParticle != NULL && fPositiveMCParticle != NULL){
			if(fNegativeMCParticle->GetMother(0) == fPositiveMCParticle->GetMother(0))
	if(fMotherMCParticle){
		iResult = kTRUE;
	}
		}
	}
	return iResult;
}

Bool_t AliV0Reader::CheckPIDProbability(Double_t negProbCut, Double_t posProbCut){
	//see header file for documentation
	
	Bool_t iResult=kFALSE;
	
	//	Double_t *posProbArray = new Double_t[10];
	//	Double_t *negProbArray = new Double_t[10];
	//-AM The TPCpid method expects an array of length kSPECIES that is 5 not 10 

	Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
	Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

	AliESDtrack* negTrack	= GetNegativeESDTrack();
	AliESDtrack* posTrack	= GetPositiveESDTrack();
	//fESDEvent->GetTrack(fCurrentV0->GetNindex());
		//fESDEvent->GetTrack(fCurrentV0->GetPindex());
	//-AM for switchtracks==true the above is a bug

	if(negProbArray && posProbArray){

		negTrack->GetTPCpid(negProbArray);
		posTrack->GetTPCpid(posProbArray);
		
		//	if(negProbArray != NULL && posProbArray != NULL){ // this is not allowed anymore for some reason(RC19)
		if(negProbArray[GetSpeciesIndex(-1)]>=negProbCut && posProbArray[GetSpeciesIndex(1)]>=posProbCut){
			iResult=kTRUE;
		}
	}
	delete [] posProbArray;
	delete [] negProbArray;
	return iResult;
}

void AliV0Reader::GetPIDProbability(Double_t &negPIDProb,Double_t & posPIDProb){
	// see header file for documentation

	//Double_t *posProbArray = new Double_t[10];
	// Double_t *negProbArray = new Double_t[10];
	//-AM The TPCpid method expects an array of length kSPECIES that is 5 not 10 
	Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
	Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

//	 AliESDtrack* negTrack	= fESDEvent->GetTrack(fCurrentV0->GetNindex());
//	 AliESDtrack* posTrack	= fESDEvent->GetTrack(fCurrentV0->GetPindex());
	//-AM for switchtracks the above is a bug
	AliESDtrack* negTrack	= GetNegativeESDTrack();
	AliESDtrack* posTrack	= GetPositiveESDTrack();

	if(negProbArray && posProbArray){
		negTrack->GetTPCpid(negProbArray);
		posTrack->GetTPCpid(posProbArray);
		
		//	if(negProbArray!=NULL && posProbArray!=NULL){ // this is not allowed anymore for some reason(RC19)
		negPIDProb = negProbArray[GetSpeciesIndex(-1)];
		posPIDProb = posProbArray[GetSpeciesIndex(1)];
	}
	delete [] posProbArray;
	delete [] negProbArray;
}

void AliV0Reader::GetPIDProbabilityMuonPion(Double_t &negPIDProb,Double_t & posPIDProb){
	// see header file for documentation


	Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
	Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

	// AliESDtrack* negTrack	= fESDEvent->GetTrack(fCurrentV0->GetNindex());
	// AliESDtrack* posTrack	= fESDEvent->GetTrack(fCurrentV0->GetPindex());
	//-AM for switchtracks the above is a bug

	AliESDtrack* negTrack	= GetNegativeESDTrack();
	AliESDtrack* posTrack	= GetPositiveESDTrack();

	if(negProbArray && posProbArray){
		negTrack->GetTPCpid(negProbArray);
		posTrack->GetTPCpid(posProbArray);
		
		//	if(negProbArray!=NULL && posProbArray!=NULL){ // this is not allowed anymore for some reason(RC19)

		negPIDProb = negProbArray[1]+negProbArray[2];
		posPIDProb = posProbArray[1]+posProbArray[2];
	}
	delete [] posProbArray;
	delete [] negProbArray;
}

void AliV0Reader::UpdateEventByEventData(){
	//see header file for documentation
	if(fCurrentEventGoodV0s->GetEntriesFast() >0 ){
		if(fCalculateBackground){
			if(fUseChargedTrackMultiplicityForBG == kTRUE){
	fBGEventHandler->AddEvent(fCurrentEventGoodV0s,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),CountESDTracks());
	//filling z and multiplicity histograms
	fHistograms->FillHistogram("ESD_Z_distribution",fESDEvent->GetPrimaryVertex()->GetZ());
	fHistograms->FillHistogram("ESD_multiplicity_distribution",CountESDTracks());
	fHistograms->FillHistogram("ESD_ZvsMultiplicity",fESDEvent->GetPrimaryVertex()->GetZ(),CountESDTracks());
			}
			else{ // means we use #V0s for multiplicity
	fBGEventHandler->AddEvent(fCurrentEventGoodV0s,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfGoodV0s);
	//filling z and multiplicity histograms
	fHistograms->FillHistogram("ESD_Z_distribution",fESDEvent->GetPrimaryVertex()->GetZ());
	fHistograms->FillHistogram("ESD_multiplicity_distribution",fNumberOfGoodV0s);
	fHistograms->FillHistogram("ESD_ZvsMultiplicity",fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfGoodV0s);
			}
		}
	}
	fCurrentEventGoodV0s->Delete();
	fCurrentV0IndexNumber=0;
	fNumberOfESDTracks=0;

	fV0Pindex.clear();
	fV0Nindex.clear();
	


	//	fBGEventHandler->PrintBGArray(); // for debugging
}


Double_t AliV0Reader::GetNegativeTrackPhi() const{
	//see header file for documentation
	
	Double_t offset=0;
	if(fNegativeTrackLorentzVector->Phi()> TMath::Pi()){
		offset = -2*TMath::Pi();
	}
	return fNegativeTrackLorentzVector->Phi()+offset;
}

Double_t AliV0Reader::GetPositiveTrackPhi() const{
	//see header file for documentation
	
	Double_t offset=0;
	if(fPositiveTrackLorentzVector->Phi()> TMath::Pi()){
		offset = -2*TMath::Pi();
	}
	return fPositiveTrackLorentzVector->Phi()+offset;
}

Double_t AliV0Reader::GetMotherCandidatePhi() const{
	//see header file for documentation
	
	Double_t offset=0;
	if(fMotherCandidateLorentzVector->Phi()> TMath::Pi()){
		offset = -2*TMath::Pi();
	}
	return fMotherCandidateLorentzVector->Phi()+offset;
}


Double_t AliV0Reader::GetMotherCandidateRapidity() const{
	//see header file for documentation
	
	Double_t rapidity=0;
	if(fMotherCandidateLorentzVector->Energy() - fMotherCandidateLorentzVector->Pz() == 0 || fMotherCandidateLorentzVector->Energy() + fMotherCandidateLorentzVector->Pz() == 0) rapidity=0;
	else rapidity = 0.5*(TMath::Log((fMotherCandidateLorentzVector->Energy() + fMotherCandidateLorentzVector->Pz()) / (fMotherCandidateLorentzVector->Energy()-fMotherCandidateLorentzVector->Pz())));
	return rapidity;
	
}





Int_t AliV0Reader::GetSpeciesIndex(Int_t chargeOfTrack){
	//see header file for documentation
	
	Int_t iResult = 10; // Unknown particle
	
	if(chargeOfTrack==-1){ //negative track
		switch(abs(fNegativeTrackPID)){
		case 11:			 //electron
			iResult = 0;
			break;
		case 13:			 //muon
			iResult = 1;
			break;
		case 211:			//pion
			iResult = 2;
			break;
		case 321:			//kaon
			iResult = 3;
			break;
		case 2212:		 //proton
			iResult = 4;
			break;
		case 22:			 //photon
			iResult = 5;
			break;
		case 111:			//pi0
			iResult = 6;
			break;
		case 2112:		 //neutron
			iResult = 7;
			break;
		case 311:			//K0
			iResult = 8;
			break;
				
			//Put in here for kSPECIES::kEleCon	????
		}
	}
	else if(chargeOfTrack==1){ //positive track
		switch(abs(fPositiveTrackPID)){
		case 11:			 //electron
			iResult = 0;
			break;
		case 13:			 //muon
			iResult = 1;
			break;
		case 211:			//pion
			iResult = 2;
			break;
		case 321:			//kaon
			iResult = 3;
			break;
		case 2212:		 //proton
			iResult = 4;
			break;
		case 22:			 //photon
			iResult = 5;
			break;
		case 111:			//pi0
			iResult = 6;
			break;
		case 2112:		 //neutron
			iResult = 7;
			break;
		case 311:			//K0
			iResult = 8;
			break;
				
			//Put in here for kSPECIES::kEleCon	????
		}
	}
	else{
		//Wrong parameter.. Print warning
	}
	return iResult;
}

Bool_t AliV0Reader::GetHelixCenter(AliESDtrack* track, Double_t b,Int_t charge, Double_t center[2]){ 
// Bool_t AliV0Reader::GetHelixCenter(const AliExternalTrackParam *track, Double_t b,Int_t charge, Double_t center[2]){
	// see header file for documentation
		
	Double_t pi = 3.14159265358979323846;
	
	Double_t	helix[6];
	track->GetHelixParameters(helix,b);
	
	Double_t xpos =	helix[5];
	Double_t ypos =	helix[0];
	Double_t radius = TMath::Abs(1./helix[4]);
	Double_t phi = helix[2];

	if(phi < 0){
		phi = phi + 2*pi;
	}

	phi -= pi/2.;
	Double_t xpoint =	radius * TMath::Cos(phi);
	Double_t ypoint =	radius * TMath::Sin(phi);

	if(b<0){
	  if(charge > 0){
		xpoint = - xpoint;
		ypoint = - ypoint;
	  }
	
	} else if(b>0){
	  if(charge < 0){
		xpoint = - xpoint;
		ypoint = - ypoint;
	  }
	}

	center[0] =	xpos + xpoint;
	center[1] =	ypos + ypoint;

	return 1;
}

/*void AliV0Reader::GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3]){

	Double_t helixcenterpos[2];
	GetHelixCenter(pparam,GetMagneticField(),pparam->Charge(),helixcenterpos);
	
	Double_t helixcenterneg[2];
	GetHelixCenter(nparam,GetMagneticField(),nparam->Charge(),helixcenterneg);
	
	Double_t helixpos[6];
	pparam->GetHelixParameters(helixpos,GetMagneticField());
	Double_t posradius = TMath::Abs(1./helixpos[4]);
	
	Double_t helixneg[6];
	nparam->GetHelixParameters(helixneg,GetMagneticField());
	Double_t negradius = TMath::Abs(1./helixneg[4]);
	
	// Calculate xy-position
	
	Double_t xpos = helixcenterpos[0];
	Double_t ypos = helixcenterpos[1];
	Double_t xneg = helixcenterneg[0];
	Double_t yneg = helixcenterneg[1];
	
	convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
	convpos[1] = (ypos*negradius + yneg*posradius)/(negradius+posradius);
	
	
	// Calculate z-position
	
	Double_t deltaXPos = convpos[0] -       xpos;
	Double_t deltaYPos = convpos[1] -      ypos;
	
	Double_t deltaXNeg = convpos[0] -      xneg;
	Double_t deltaYNeg = convpos[1] -      yneg;
	
	Double_t alphaPos =    TMath::Pi() + TMath::ATan2(-deltaYPos,-deltaXPos);
	Double_t alphaNeg =    TMath::Pi() + TMath::ATan2(-deltaYNeg,-deltaXNeg);
	Double_t vertexXNeg =  xneg +  TMath::Abs(negradius)*TMath::Cos(alphaNeg);
	Double_t vertexYNeg =  yneg +  TMath::Abs(negradius)*TMath::Sin(alphaNeg);
	
	Double_t vertexXPos =  xpos +  TMath::Abs(posradius)*TMath::Cos(alphaPos);
	Double_t vertexYPos =  ypos +  TMath::Abs(posradius)*TMath::Sin(alphaPos);
	
	Double_t x0neg =        helixneg[5];
	Double_t y0neg =        helixneg[0];
	
	Double_t x0pos =        helixpos[5];
	Double_t y0pos =        helixpos[0];
	
	Double_t dNeg = TMath::Sqrt((vertexXNeg - x0neg)*(vertexXNeg - x0neg)+(vertexYNeg - y0neg)*(vertexYNeg - y0neg));
	
	Double_t dPos = TMath::Sqrt((vertexXPos - x0pos)*(vertexXPos - x0pos)+(vertexYPos - y0pos)*(vertexYPos - y0pos));
	
	Double_t rNeg = TMath::Sqrt(negradius*negradius - dNeg*dNeg/4.);

	Double_t rPos = TMath::Sqrt(posradius*posradius - dPos*dPos/4.);
	
	Double_t deltabetaNeg = 2*(TMath::Pi() + TMath::ATan2(-dNeg/2.,-rNeg));
	Double_t deltabetaPos = 2*(TMath::Pi() + TMath::ATan2(-dPos/2.,-rPos));
	
	Double_t deltaUNeg = negradius*deltabetaNeg;
	Double_t deltaUPos = posradius*deltabetaPos;

	Double_t zphaseNeg = nparam->GetZ() +  deltaUNeg * nparam->GetTgl();
	Double_t zphasePos = pparam->GetZ() +  deltaUPos * pparam->GetTgl();

	convpos[2] = (zphasePos*negradius+zphaseNeg*posradius)/(negradius+posradius);

  

}*/

Bool_t AliV0Reader::GetConvPosXY(AliESDtrack* ptrack, AliESDtrack* ntrack, Double_t b, Double_t convpos[2]){
	//see header file for documentation

	Double_t helixcenterpos[2];
	GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

	Double_t helixcenterneg[2];
	GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

	Double_t	poshelix[6];
	ptrack->GetHelixParameters(poshelix,b);
	Double_t posradius = TMath::Abs(1./poshelix[4]);

	Double_t	neghelix[6];
	ntrack->GetHelixParameters(neghelix,b);
	Double_t negradius = TMath::Abs(1./neghelix[4]);

	Double_t xpos = helixcenterpos[0];
	Double_t ypos = helixcenterpos[1];
	Double_t xneg = helixcenterneg[0];
	Double_t yneg = helixcenterneg[1];

	convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
	convpos[1] = (ypos*negradius+	yneg*posradius)/(negradius+posradius);

	return 1;
}



Double_t AliV0Reader::GetConvPosZ(AliESDtrack* ptrack,AliESDtrack* ntrack, Double_t b){
	//see header file for documentation

	Double_t	helixpos[6];
	ptrack->GetHelixParameters(helixpos,b);

	Double_t	helixneg[6];
	ntrack->GetHelixParameters(helixneg,b);

	Double_t negtrackradius =	TMath::Abs(1./helixneg[4]);
	Double_t postrackradius =	TMath::Abs(1./helixpos[4]);

	Double_t pi = 3.14159265358979323846;

	Double_t convpos[2];
	GetConvPosXY(ptrack,ntrack,b,convpos);

	 Double_t convposx = convpos[0];
	 Double_t convposy = convpos[1];

	 Double_t helixcenterpos[2];
	 GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

	 Double_t helixcenterneg[2];
	 GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

	 Double_t xpos = helixcenterpos[0];
	 Double_t ypos = helixcenterpos[1];
	 Double_t xneg = helixcenterneg[0];
	 Double_t yneg = helixcenterneg[1];

	 Double_t deltaXPos = convposx -	xpos;
	 Double_t deltaYPos = convposy -	ypos;

	 Double_t deltaXNeg = convposx -	xneg;
	 Double_t deltaYNeg = convposy -	yneg;

	 Double_t alphaPos =	pi + TMath::ATan2(-deltaYPos,-deltaXPos);
	 Double_t alphaNeg =	pi + TMath::ATan2(-deltaYNeg,-deltaXNeg);

	 Double_t vertexXNeg =	xneg +	TMath::Abs(negtrackradius)*
	 TMath::Cos(alphaNeg);
	 Double_t vertexYNeg =	yneg +	TMath::Abs(negtrackradius)*
	 TMath::Sin(alphaNeg);

	 Double_t vertexXPos =	xpos +	TMath::Abs(postrackradius)*
	 TMath::Cos(alphaPos);
	 Double_t vertexYPos =	ypos +	TMath::Abs(postrackradius)*
	 TMath::Sin(alphaPos);

	 Double_t x0neg =	 helixneg[5];
	 Double_t y0neg =	 helixneg[0];

	 Double_t x0pos =	 helixpos[5];
	 Double_t y0pos =	 helixpos[0];

	 Double_t dNeg = TMath::Sqrt((vertexXNeg -	x0neg)*(vertexXNeg - x0neg)
															 +(vertexYNeg -	y0neg)*(vertexYNeg - y0neg));

	 Double_t dPos = TMath::Sqrt((vertexXPos -	x0pos)*(vertexXPos - x0pos)
															 +(vertexYPos -	y0pos)*(vertexYPos - y0pos));

	 Double_t rNeg =	TMath::Sqrt(negtrackradius*negtrackradius -
	 dNeg*dNeg/4.);

	 Double_t rPos = TMath::Sqrt(postrackradius*postrackradius -
	 dPos*dPos/4.);

	 Double_t deltabetaNeg =	2*(pi +	 TMath::ATan2(-dNeg/2.,-rNeg));
	 Double_t deltabetaPos = 2*(pi + TMath::ATan2(-dPos/2.,-rPos));

	 Double_t deltaUNeg = negtrackradius*deltabetaNeg;
	 Double_t deltaUPos = postrackradius*deltabetaPos;

	 Double_t zphaseNeg = ntrack->GetZ() +	deltaUNeg * ntrack->GetTgl();
	 Double_t zphasePos = ptrack->GetZ() +	deltaUPos * ptrack->GetTgl();

	 Double_t convposz = (zphasePos*negtrackradius+zphaseNeg*postrackradius)/(negtrackradius+postrackradius);

	 return convposz;
}

AliGammaConversionKFVector* AliV0Reader::GetBGGoodV0s(Int_t /*event*/) const{
	/*
	if(fUseChargedTrackMultiplicityForBG == kTRUE){
		return fBGEventHandler->GetBGGoodV0s(event,fESDEvent->GetPrimaryVertex()->GetZ(),CountESDTracks());
	}
	else{ // means we use #v0s as multiplicity
		return fBGEventHandler->GetBGGoodV0s(event,fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfGoodV0s);
	}
	*/
	return NULL;
}

Int_t AliV0Reader::CountESDTracks(){
	// see header file for documentation
	if(fNumberOfESDTracks == 0){ // count the good esd tracks
		for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);			
			if(!curTrack){
	continue;
			}
			if(fEsdTrackCuts->AcceptTrack(curTrack) ){
	fNumberOfESDTracks++;
			}
		}
	}

	return fNumberOfESDTracks;
}

Bool_t AliV0Reader::CheckIfPi0IsMother(Int_t label){
	// see headerfile for documentation
	Bool_t iResult=kFALSE;
	//	cout<<"Checking particle label, particle is: "<<fMCStack->Particle(TMath::Abs(label))->GetName()<<endl;
	if(fMCStack->Particle(TMath::Abs(label))->GetPdgCode() == 111){
		iResult=kTRUE;
	}
	return iResult;
}

Bool_t AliV0Reader::CheckIfEtaIsMother(Int_t label){
	// see headerfile for documentation
	Bool_t iResult=kFALSE;
	//	cout<<"Checking particle label, particle is: "<<fMCStack->Particle(TMath::Abs(label))->GetName()<<endl;
	if(fMCStack->Particle(TMath::Abs(label))->GetPdgCode() == 221){
		iResult=kTRUE;
	}
	return iResult;
}


Bool_t AliV0Reader::GetArmenterosQtAlfa(const AliKFParticle* negativeKFParticle, const AliKFParticle * positiveKFParticle, const AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] ){
	//see header file for documentation

	TVector3 momentumVectorPositiveKF(positiveKFParticle->GetPx(),positiveKFParticle->GetPy(),positiveKFParticle->GetPz());
	TVector3 momentumVectorNegativeKF(negativeKFParticle->GetPx(),negativeKFParticle->GetPy(),negativeKFParticle->GetPz());
	TVector3 vecV0(gammaKFCandidate->GetPx(),gammaKFCandidate->GetPy(),gammaKFCandidate->GetPz());

	Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
	Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));
	
	Float_t alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
	

	Float_t qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);
			
	armenterosQtAlfa[0]=qt;
	armenterosQtAlfa[1]=alfa;

	return 1;

}

Bool_t AliV0Reader::GetArmenterosQtAlfa(const TParticle* negativeParticle, const TParticle * positiveParticle, const AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] ){
	//see header file for documentation

	TVector3 momentumVectorPositiveKF(positiveParticle->Px(),positiveParticle->Py(),positiveParticle->Pz());
	TVector3 momentumVectorNegativeKF(negativeParticle->Px(),negativeParticle->Py(),negativeParticle->Pz());
	TVector3 vecV0(gammaKFCandidate->GetPx(),gammaKFCandidate->GetPy(),gammaKFCandidate->GetPz());

	Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
	Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));
	
	Float_t alfa;
	Float_t qt;
	if ( positiveParticle->GetPdgCode() == 11 || positiveParticle->GetPdgCode() == 13 || positiveParticle->GetPdgCode() == 15){
		alfa =((momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)-(momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorNegativeKF.Mag()*TMath::Sin(thetaV0neg);
	} else if ( negativeParticle->GetPdgCode() == -11 || negativeParticle->GetPdgCode() == -13 || negativeParticle->GetPdgCode() == -15){
		alfa =((momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)-(momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorNegativeKF.Mag()*TMath::Sin(thetaV0neg);
	} else if (positiveParticle->GetPdgCode() < 0 && positiveParticle->GetPdgCode() != -11 && positiveParticle->GetPdgCode() != -13 && positiveParticle->GetPdgCode() != -15 ){
		alfa =((momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)-(momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorNegativeKF.Mag()*TMath::Sin(thetaV0neg);
	} else {
		alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);
	}
		
	armenterosQtAlfa[0]=qt;
	armenterosQtAlfa[1]=alfa;
	return 1;

}

Bool_t AliV0Reader::GetArmenterosQtAlfa(const TParticle* negativeParticle, const TParticle * positiveParticle, const TParticle * gammaCandidate, Double_t armenterosQtAlfa[2] ){
	//see header file for documentation

	TVector3 momentumVectorPositiveKF(positiveParticle->Px(),positiveParticle->Py(),positiveParticle->Pz());
	TVector3 momentumVectorNegativeKF(negativeParticle->Px(),negativeParticle->Py(),negativeParticle->Pz());
	TVector3 vecV0(gammaCandidate->Px(),gammaCandidate->Py(),gammaCandidate->Pz());

	Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
	Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));
	
	Float_t alfa;
	Float_t qt;
	if ( positiveParticle->GetPdgCode() == 11 || positiveParticle->GetPdgCode() == 13 || positiveParticle->GetPdgCode() == 15){
		alfa =((momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)-(momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorNegativeKF.Mag()*TMath::Sin(thetaV0neg);
	} else if ( negativeParticle->GetPdgCode() == -11 || negativeParticle->GetPdgCode() == -13 || negativeParticle->GetPdgCode() == -15){
		alfa =((momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)-(momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorNegativeKF.Mag()*TMath::Sin(thetaV0neg);
	} else if (positiveParticle->GetPdgCode() < 0 && positiveParticle->GetPdgCode() != -11 && positiveParticle->GetPdgCode() != -13 && positiveParticle->GetPdgCode() != -15 ){
		alfa =((momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)-(momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorNegativeKF.Mag()*TMath::Sin(thetaV0neg);
	} else {
		alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
		((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
		qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);
	}
			
	armenterosQtAlfa[0]=qt;
	armenterosQtAlfa[1]=alfa;
	return 1;

}


Int_t AliV0Reader::GetFirstTPCRow(Double_t radius){


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
void AliV0Reader::SmearKFParticle(AliKFParticle * kfParticle)
{
	Double_t facPBrem = 1.;
	Double_t facPSig = 0.;

	Double_t phi=0.;
	Double_t theta=0.;
	Double_t P=0.;

	P=kfParticle->GetP();
	phi=kfParticle->GetPhi();
	if( kfParticle->GetP()!=0){
		theta=acos( kfParticle->Pz()/ kfParticle->GetP());
	}

	if( fPSigSmearing != 0. || fPSigSmearingCte!=0. ){ 
		facPSig = TMath::Sqrt(fPSigSmearingCte*fPSigSmearingCte+fPSigSmearing*fPSigSmearing*P*P)*fRandom.Gaus(0.,1.);
	}
	
	if( fPBremSmearing != 1.){
		if(fBrem!=NULL){
			facPBrem = fBrem->GetRandom();
		}
	}

	kfParticle->Px() = facPBrem* (1+facPSig)* P*sin(theta)*cos(phi) ;
	kfParticle->Py() = facPBrem* (1+facPSig)* P*sin(theta)*sin(phi) ;
	kfParticle->Pz() = facPBrem* (1+facPSig)* P*cos(theta) ;
	kfParticle->E() = kfParticle->GetP();
}
