/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
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
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliKFVertex.h"

#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliGammaConversionBGHandler.h"
#include "AliESDtrackCuts.h"


class iostream;
class AliESDv0;
class TFormula;

using namespace std;

ClassImp(AliV0Reader)


AliESDpid* AliV0Reader::fgESDpid = 0x0;

AliV0Reader::AliV0Reader() :
  TObject(),
  fMCStack(NULL),
  // fMCTruth(NULL),
  fMCEvent(NULL),    // for CF
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
  fMaxR(10000),// 100 meter(outside of ALICE)
  fEtaCut(0.),
  fPtCut(0.),
  fSinglePtCut(0.),
  fMaxZ(0.),
  fMinClsTPC(0.),
  fLineCutZRSlope(0.),
  fLineCutZValue(0.),
  fChi2CutConversion(0.),
  fChi2CutMeson(0.),
  fAlphaCutMeson(1.),
  fPIDProbabilityCutNegativeParticle(0),
  fPIDProbabilityCutPositiveParticle(0),
  fDodEdxSigmaCut(kFALSE),
  fPIDnSigmaAboveElectronLine(100),
  fPIDnSigmaBelowElectronLine(-100),
  fPIDnSigmaAbovePionLine(-100), 
  fPIDMinPnSigmaAbovePionLine(100), 
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
  fQtMax(100.),
  fXVertexCut(0.),
  fYVertexCut(0.),
  fZVertexCut(0.),
  fNSigmaMass(0.),
  fUseImprovedVertex(kFALSE),
  fUseOwnXYZCalculation(kFALSE),
  fDoCF(kFALSE),
  fUseOnFlyV0Finder(kTRUE),
  fUpdateV0AlreadyCalled(kFALSE),
  fCurrentEventGoodV0s(NULL),
//  fPreviousEventGoodV0s(),
  fCalculateBackground(kFALSE),
  fBGEventHandler(NULL),
  fBGEventInitialized(kFALSE),
  fEsdTrackCuts(NULL),
  fNumberOfESDTracks(0),
  nEventsForBGCalculation(10)
{
  //fESDpid = new AliESDpid;	
}


AliV0Reader::AliV0Reader(const AliV0Reader & original) :
  TObject(original),
  fMCStack(original.fMCStack),
  // fMCTruth(original.fMCTruth),
  fMCEvent(original.fMCEvent),  // for CF
  fChain(original.fChain),
  //  fESDHandler(original.fESDHandler),
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
  fMaxR(original.fMaxR),
  fEtaCut(original.fEtaCut),
  fPtCut(original.fPtCut),
  fSinglePtCut(original.fSinglePtCut),
  fMaxZ(original.fMaxZ),
  fMinClsTPC(original.fMinClsTPC),
  fLineCutZRSlope(original.fLineCutZRSlope),
  fLineCutZValue(original.fLineCutZValue),
  fChi2CutConversion(original.fChi2CutConversion),
  fChi2CutMeson(original.fChi2CutMeson),
  fAlphaCutMeson(original.fAlphaCutMeson),
  fPIDProbabilityCutNegativeParticle(original.fPIDProbabilityCutNegativeParticle),
  fPIDProbabilityCutPositiveParticle(original.fPIDProbabilityCutPositiveParticle),
  fDodEdxSigmaCut(original.fDodEdxSigmaCut),
  fPIDnSigmaAboveElectronLine(original.fPIDnSigmaAboveElectronLine),
  fPIDnSigmaBelowElectronLine(original.fPIDnSigmaBelowElectronLine),
  fPIDnSigmaAbovePionLine(original.fPIDnSigmaAbovePionLine), 
  fPIDMinPnSigmaAbovePionLine(original.fPIDMinPnSigmaAbovePionLine), 
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
  fQtMax(original.fQtMax),
  fXVertexCut(original.fXVertexCut),
  fYVertexCut(original.fYVertexCut),
  fZVertexCut(original.fZVertexCut),
  fNSigmaMass(original.fNSigmaMass),
  fUseImprovedVertex(original.fUseImprovedVertex),
  fUseOwnXYZCalculation(original.fUseOwnXYZCalculation),
  fDoCF(original.fDoCF),
  fUseOnFlyV0Finder(original.fUseOnFlyV0Finder),
  fUpdateV0AlreadyCalled(original.fUpdateV0AlreadyCalled),
  fCurrentEventGoodV0s(original.fCurrentEventGoodV0s),
  //  fPreviousEventGoodV0s(original.fPreviousEventGoodV0s),
  fCalculateBackground(original.fCalculateBackground),
  fBGEventHandler(original.fBGEventHandler),
  fBGEventInitialized(original.fBGEventInitialized),
  fEsdTrackCuts(original.fEsdTrackCuts),
  fNumberOfESDTracks(original.fNumberOfESDTracks),
  nEventsForBGCalculation(original.nEventsForBGCalculation)
{
	
}


AliV0Reader & AliV0Reader::operator = (const AliV0Reader & /*source*/)
{
  // assignment operator
  return *this;
}
AliV0Reader::~AliV0Reader()
{
  //  if(fESDpid){
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



  //  fMCTruth = mcH->MCEvent();
  //  fMC = mcH->MCEvent();
  // stack = fMC->Stack();


  //if(fMCTruth == NULL){
    //print warning here
  // fDoMC = kFALSE;
  //}

  if(fMCEvent == NULL){
   fDoMC = kFALSE;
  }

  //Get pointer to the mc stack
  //  if(fMCTruth){
  if(fMCEvent){
    fMCStack = fMCEvent->Stack();
    if(fMCStack == NULL){
      //print warning here
    }
    // Better parameters for MonteCarlo from A. Kalweit 2010/01/8
//     fESDpid->GetTPCResponse().SetBetheBlochParameters( 2.15898e+00/50.,
// 				      1.75295e+01,
// 				      3.40030e-09,
// 				      1.96178e+00,
// 				      3.91720e+00);
  }
  else{
    // Better parameters for data from A. Kalweit 2010/01/8
 //    fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,
// 				     2.63394e+01,
// 				     5.04114e-11,
// 				     2.12543e+00,
// 				     4.88663e+00);
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
	
  AliKFParticle::SetField(fESDEvent->GetMagneticField());

  //  fCurrentEventGoodV0s = new TClonesArray("TClonesArray", 0);
  if(fCurrentEventGoodV0s == NULL){
    fCurrentEventGoodV0s = new TClonesArray("AliKFParticle", 0);
  }

  if(fCalculateBackground == kTRUE){
    if(fBGEventInitialized == kFALSE){

      
      Double_t *zBinLimitsArray = new Double_t[8];
      zBinLimitsArray[0] = -50.00;
      zBinLimitsArray[1] = -4.07;
      zBinLimitsArray[2] = -2.17;
      zBinLimitsArray[3] = -0.69;
      zBinLimitsArray[4] = 0.69;
      zBinLimitsArray[5] = 2.17;
      zBinLimitsArray[6] = 4.11;
      zBinLimitsArray[7] = 50.00;
      
      Double_t *multiplicityBinLimitsArray= new Double_t[5];
      multiplicityBinLimitsArray[0] = 0;
      multiplicityBinLimitsArray[1] = 8.5;
      multiplicityBinLimitsArray[2] = 16.5;
      multiplicityBinLimitsArray[3] = 27.5;
      multiplicityBinLimitsArray[4] = 41.5;
          
      fBGEventHandler = new AliGammaConversionBGHandler(8,5,nEventsForBGCalculation);
      
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

Bool_t AliV0Reader::CheckForPrimaryVertex(){
  //see headerfile for documentation
  return fESDEvent->GetPrimaryVertex()->GetNContributors()>0;
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
    GetArmenterosQtAlfa(GetNegativeKFParticle(), 
			GetPositiveKFParticle(), 
			GetMotherCandidateKFCombination(),
			armenterosQtAlfa);
   
    fHistograms->FillHistogram("ESD_AllV0sCurrentFinder_alfa_qt",armenterosQtAlfa[1],armenterosQtAlfa[0]);
 
   
    if(fCurrentNegativeESDTrack->GetSign() == fCurrentPositiveESDTrack->GetSign()){             // avoid like sign
      //  iResult=kFALSE;
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
      //  if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kITSrefit) || 
      //      !(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kITSrefit) ){
      //  iResult=kFALSE;
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
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepdEdx_electronselection);               // for CF
      }

      if( fCurrentPositiveESDTrack->P()>fPIDMinPnSigmaAbovePionLine){
	if(fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
	   fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
	   fgESDpid->NumberOfSigmasTPC(fCurrentPositiveESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
	  //	  iResult=kFALSE;
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
      
      if( fCurrentNegativeESDTrack->P()>fPIDMinPnSigmaAbovePionLine){
	if(fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
	   fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
	   fgESDpid->NumberOfSigmasTPC(fCurrentNegativeESDTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
	  //	  iResult=kFALSE;
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
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepdEdx_pionrejection);               // for CF
      }

    }
    //    Float_t fPIDMinPKaonRejectionLowP=1.5;
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
    //    Float_t fPIDMinPProtonRejection=2;
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
    //    Float_t fPIDMinPPionRejection=0.3;
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

    // Gamma selection based on QT from Armenteros
    if(fDoQtGammaSelection == kTRUE){
      if(armenterosQtAlfa[0]>fQtMax){
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("ESD_CutQt_InvMass",GetMotherCandidateMass());
	}
	fCurrentV0IndexNumber++;
	continue;
      }
    }

    //checks if we have a prim vertex
    if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) { 
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
		
		
    if((TMath::Abs(fCurrentZValue)*fLineCutZRSlope)-fLineCutZValue > GetXYRadius() ){ // cuts out regions where we do not reconstruct
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
		
    /* Moved further up so corr framework can work
       if(UpdateV0Information() == kFALSE){
       fCurrentV0IndexNumber++;
       continue;
       }
    */
    if(fCurrentNegativeESDTrack->GetNcls(1) < fMinClsTPC ||  fCurrentPositiveESDTrack->GetNcls(1) < fMinClsTPC ){
      if(fHistograms != NULL){
	fHistograms->FillHistogram("ESD_CutMinNClsTPC_InvMass",GetMotherCandidateMass());
      }
      fCurrentV0IndexNumber++;
      continue;
    }
    if(fDoCF){
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepMinClsTPC);		// for CF	
    }

		
    if(fUseKFParticle){


      if( fCurrentNegativeKFParticle->GetPt()< fSinglePtCut ||  fCurrentPositiveKFParticle->GetPt()< fSinglePtCut){
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
			
      if(TMath::Abs(fMotherCandidateLorentzVector->Eta())> fEtaCut){
	if(fHistograms != NULL){
	  fHistograms->FillHistogram("ESD_CutEta_InvMass",GetMotherCandidateMass());
	}
	fCurrentV0IndexNumber++;
	continue;
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

    //    fCurrentEventGoodV0s.push_back(*fCurrentMotherKFCandidate);

    new((*fCurrentEventGoodV0s)[fCurrentEventGoodV0s->GetEntriesFast()])  AliKFParticle(*fCurrentMotherKFCandidate);

    iResult=kTRUE;//means we have a v0 who survived all the cuts applied
		
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

  if(fCurrentPositiveESDTrack->GetSign() == -1 && fCurrentNegativeESDTrack->GetSign() == 1){  // switch wrong signed tracks
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
  fCurrentMotherKFCandidate = new AliKFParticle(*fCurrentNegativeKFParticle,*fCurrentPositiveKFParticle);
	
	
  if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
    fCurrentMotherKFCandidate->SetMassConstraint(0,fNSigmaMass);
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
  else if(fUseESDTrack){
    fNegativeTrackLorentzVector = new TLorentzVector(fCurrentNegativeESDTrack->Px(),fCurrentNegativeESDTrack->Py(),fCurrentNegativeESDTrack->Pz());
  }
	
  if(fPositiveTrackLorentzVector != NULL){
    delete fPositiveTrackLorentzVector;
  }
  if(fUseKFParticle){
    fPositiveTrackLorentzVector = new TLorentzVector(fCurrentPositiveKFParticle->Px(),fCurrentPositiveKFParticle->Py(),fCurrentPositiveKFParticle->Pz());
  }
  else if(fUseESDTrack){
    fPositiveTrackLorentzVector = new TLorentzVector(fCurrentPositiveESDTrack->Px(),fCurrentPositiveESDTrack->Py(),fCurrentPositiveESDTrack->Pz());
  }
	
  if(fMotherCandidateLorentzVector != NULL){
    delete fMotherCandidateLorentzVector;
  }
  if(fUseKFParticle){
    fMotherCandidateLorentzVector = new TLorentzVector(*fNegativeTrackLorentzVector + *fPositiveTrackLorentzVector);
  }
  else if(fUseESDTrack){
    fMotherCandidateLorentzVector = new TLorentzVector(*fNegativeTrackLorentzVector + *fPositiveTrackLorentzVector);
  }
	
  if(fPositiveTrackPID==-11 && fNegativeTrackPID==11){
    fMotherCandidateLorentzVector->SetXYZM(fMotherCandidateLorentzVector->Px() ,fMotherCandidateLorentzVector->Py(),fMotherCandidateLorentzVector->Pz(),0.); 
  }
    
	
  if(fDoMC == kTRUE){
    fMotherMCParticle= NULL;
    fNegativeMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetNindex())->GetLabel()));
    fPositiveMCParticle = fMCStack->Particle(TMath::Abs(fESDEvent->GetTrack(fCurrentV0->GetPindex())->GetLabel()));
    if(fPositiveMCParticle->GetMother(0)>-1){
      fMotherMCParticle = fMCStack->Particle(fPositiveMCParticle->GetMother(0));
    }
  }
	

  // for CF
//   Double_t containerInput[3];
//   if(fDoCF){
//     containerInput[0] = GetMotherCandidatePt();
//     containerInput[1] = GetMotherCandidateEta();
//     containerInput[2] = GetMotherCandidateMass();
    
//     fCFManager->GetParticleContainer()->Fill(containerInput,kStepLikeSign);		// for CF	
//     fCFManager->GetParticleContainer()->Fill(containerInput,kStepTPCRefit);		// for CF	
//     fCFManager->GetParticleContainer()->Fill(containerInput,kStepKinks);		// for CF	
//   }
  

  if(fUseOwnXYZCalculation == kFALSE){
    fCurrentV0->GetXYZ(fCurrentXValue,fCurrentYValue,fCurrentZValue);
  }
  else{
    Double_t convpos[2];
    convpos[0]=0;
    convpos[1]=0;

    GetConvPosXY(GetPositiveESDTrack(),GetNegativeESDTrack(),GetMagneticField(),convpos);
    fCurrentXValue = convpos[0];
    fCurrentYValue = convpos[1];
    fCurrentZValue = GetConvPosZ(GetPositiveESDTrack(),GetNegativeESDTrack(),GetMagneticField());
  }
  /*
  if(fCurrentNegativeESDTrack->GetSign() == fCurrentPositiveESDTrack->GetSign()){             // avoid like sign
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
    //  if( !(fCurrentNegativeESDTrack->GetStatus() & AliESDtrack::kITSrefit) || 
    //      !(fCurrentPositiveESDTrack->GetStatus() & AliESDtrack::kITSrefit) ){
    iResult=kFALSE;
    if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE  && doFillHistos == kTRUE){
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
      if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE  && doFillHistos == kTRUE){
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
	if(fHistograms != NULL && fUpdateV0AlreadyCalled == kFALSE  && doFillHistos == kTRUE){
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
	
  //  Double_t *posProbArray = new Double_t[10];
  //  Double_t *negProbArray = new Double_t[10];
  //-AM The TPCpid method expects an array of length kSPECIES that is 5 not 10 

  Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
  Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

  AliESDtrack* negTrack  = GetNegativeESDTrack();
  AliESDtrack* posTrack  = GetPositiveESDTrack();
  //fESDEvent->GetTrack(fCurrentV0->GetNindex());
    //fESDEvent->GetTrack(fCurrentV0->GetPindex());
  //-AM for switchtracks==true the above is a bug

  negTrack->GetTPCpid(negProbArray);
  posTrack->GetTPCpid(posProbArray);
	
  //  if(negProbArray != NULL && posProbArray != NULL){ // this is not allowed anymore for some reason(RC19)
  if(negProbArray && posProbArray){
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

//   AliESDtrack* negTrack  = fESDEvent->GetTrack(fCurrentV0->GetNindex());
//   AliESDtrack* posTrack  = fESDEvent->GetTrack(fCurrentV0->GetPindex());
  //-AM for switchtracks the above is a bug
  AliESDtrack* negTrack  = GetNegativeESDTrack();
  AliESDtrack* posTrack  = GetPositiveESDTrack();


  negTrack->GetTPCpid(negProbArray);
  posTrack->GetTPCpid(posProbArray);
	
  //  if(negProbArray!=NULL && posProbArray!=NULL){ // this is not allowed anymore for some reason(RC19)
  if(negProbArray && posProbArray){
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

  // AliESDtrack* negTrack  = fESDEvent->GetTrack(fCurrentV0->GetNindex());
  // AliESDtrack* posTrack  = fESDEvent->GetTrack(fCurrentV0->GetPindex());
  //-AM for switchtracks the above is a bug

  AliESDtrack* negTrack  = GetNegativeESDTrack();
  AliESDtrack* posTrack  = GetPositiveESDTrack();

  negTrack->GetTPCpid(negProbArray);
  posTrack->GetTPCpid(posProbArray);
	
  //  if(negProbArray!=NULL && posProbArray!=NULL){ // this is not allowed anymore for some reason(RC19)
  if(negProbArray && posProbArray){
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
      fBGEventHandler->AddEvent(fCurrentEventGoodV0s,fESDEvent->GetPrimaryVertex()->GetZ(),CountESDTracks());
      //filling z and multiplicity histograms
      fHistograms->FillHistogram("ESD_Z_distribution",fESDEvent->GetPrimaryVertex()->GetZ());
      fHistograms->FillHistogram("ESD_multiplicity_distribution",CountESDTracks());
      fHistograms->FillHistogram("ESD_ZvsMultiplicity",fESDEvent->GetPrimaryVertex()->GetZ(),CountESDTracks());
    }
  }
  fCurrentEventGoodV0s->Delete();
  fCurrentV0IndexNumber=0;
  fNumberOfESDTracks=0;
  //  fBGEventHandler->PrintBGArray(); // for debugging
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
    case 11:       //electron
      iResult = 0;
      break;
    case 13:       //muon
      iResult = 1;
      break;
    case 211:      //pion
      iResult = 2;
      break;
    case 321:      //kaon
      iResult = 3;
      break;
    case 2212:     //proton
      iResult = 4;
      break;
    case 22:       //photon
      iResult = 5;
      break;
    case 111:      //pi0
      iResult = 6;
      break;
    case 2112:     //neutron
      iResult = 7;
      break;
    case 311:      //K0
      iResult = 8;
      break;
				
      //Put in here for kSPECIES::kEleCon  ????
    }
  }
  else if(chargeOfTrack==1){ //positive track
    switch(abs(fPositiveTrackPID)){
    case 11:       //electron
      iResult = 0;
      break;
    case 13:       //muon
      iResult = 1;
      break;
    case 211:      //pion
      iResult = 2;
      break;
    case 321:      //kaon
      iResult = 3;
      break;
    case 2212:     //proton
      iResult = 4;
      break;
    case 22:       //photon
      iResult = 5;
      break;
    case 111:      //pi0
      iResult = 6;
      break;
    case 2112:     //neutron
      iResult = 7;
      break;
    case 311:      //K0
      iResult = 8;
      break;
				
      //Put in here for kSPECIES::kEleCon  ????
    }
  }
  else{
    //Wrong parameter.. Print warning
  }
  return iResult;
}

Bool_t AliV0Reader::GetHelixCenter(AliESDtrack* track, Double_t b,Int_t charge, Double_t center[2]){
  // see header file for documentation
  
  Double_t pi = 3.14159265358979323846;
  
  Double_t  helix[6];
  track->GetHelixParameters(helix,b);
  
  Double_t xpos =  helix[5];
  Double_t ypos =  helix[0];
  Double_t radius = TMath::Abs(1./helix[4]);
  Double_t phi = helix[2];

  if(phi < 0){
    phi = phi + 2*pi;
  }

  phi -= pi/2.;
  Double_t xpoint =  radius * TMath::Cos(phi);
  Double_t ypoint =  radius * TMath::Sin(phi);

  if(b<0){
    if(charge > 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }

    if(charge < 0){
      xpoint =  xpoint;
      ypoint =  ypoint;
    }
  }
  if(b>0){
    if(charge > 0){
      xpoint =  xpoint;
      ypoint =  ypoint;
    }

    if(charge < 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
  }
  center[0] =  xpos + xpoint;
  center[1] =  ypos + ypoint;

  return 1;
}

Bool_t AliV0Reader::GetConvPosXY(AliESDtrack* ptrack, AliESDtrack* ntrack, Double_t b, Double_t convpos[2]){
  //see header file for documentation

  Double_t helixcenterpos[2];
  GetHelixCenter(ptrack,b,ptrack->Charge(),helixcenterpos);

  Double_t helixcenterneg[2];
  GetHelixCenter(ntrack,b,ntrack->Charge(),helixcenterneg);

  Double_t  poshelix[6];
  ptrack->GetHelixParameters(poshelix,b);
  Double_t posradius = TMath::Abs(1./poshelix[4]);

  Double_t  neghelix[6];
  ntrack->GetHelixParameters(neghelix,b);
  Double_t negradius = TMath::Abs(1./neghelix[4]);

  Double_t xpos = helixcenterpos[0];
  Double_t ypos = helixcenterpos[1];
  Double_t xneg = helixcenterneg[0];
  Double_t yneg = helixcenterneg[1];

  convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
  convpos[1] = (ypos*negradius+  yneg*posradius)/(negradius+posradius);

  return 1;
}



Double_t AliV0Reader::GetConvPosZ(AliESDtrack* ptrack,AliESDtrack* ntrack, Double_t b){
  //see header file for documentation

  Double_t  helixpos[6];
  ptrack->GetHelixParameters(helixpos,b);

  Double_t  helixneg[6];
  ntrack->GetHelixParameters(helixneg,b);

  Double_t negtrackradius =  TMath::Abs(1./helixneg[4]);
  Double_t postrackradius =  TMath::Abs(1./helixpos[4]);

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

   Double_t deltaXPos = convposx -  xpos;
   Double_t deltaYPos = convposy -  ypos;

   Double_t deltaXNeg = convposx -  xneg;
   Double_t deltaYNeg = convposy -  yneg;

   Double_t alphaPos =  pi + TMath::ATan2(-deltaYPos,-deltaXPos);
   Double_t alphaNeg =  pi + TMath::ATan2(-deltaYNeg,-deltaXNeg);

   Double_t vertexXNeg =  xneg +  TMath::Abs(negtrackradius)*
   TMath::Cos(alphaNeg);
   Double_t vertexYNeg =  yneg +  TMath::Abs(negtrackradius)*
   TMath::Sin(alphaNeg);

   Double_t vertexXPos =  xpos +  TMath::Abs(postrackradius)*
   TMath::Cos(alphaPos);
   Double_t vertexYPos =  ypos +  TMath::Abs(postrackradius)*
   TMath::Sin(alphaPos);

   Double_t x0neg =   helixneg[5];
   Double_t y0neg =   helixneg[0];

   Double_t x0pos =   helixpos[5];
   Double_t y0pos =   helixpos[0];

   Double_t dNeg = TMath::Sqrt((vertexXNeg -  x0neg)*(vertexXNeg - x0neg)
                               +(vertexYNeg -  y0neg)*(vertexYNeg - y0neg));

   Double_t dPos = TMath::Sqrt((vertexXPos -  x0pos)*(vertexXPos - x0pos)
                               +(vertexYPos -  y0pos)*(vertexYPos - y0pos));

   Double_t rNeg =  TMath::Sqrt(negtrackradius*negtrackradius -
   dNeg*dNeg/4.);

   Double_t rPos = TMath::Sqrt(postrackradius*postrackradius -
   dPos*dPos/4.);

   Double_t deltabetaNeg =  2*(pi +   TMath::ATan2(-dNeg/2.,-rNeg));
   Double_t deltabetaPos = 2*(pi + TMath::ATan2(-dPos/2.,-rPos));

   Double_t deltaUNeg = negtrackradius*deltabetaNeg;
   Double_t deltaUPos = postrackradius*deltabetaPos;

   Double_t zphaseNeg = ntrack->GetZ() +  deltaUNeg * ntrack->GetTgl();
   Double_t zphasePos = ptrack->GetZ() +  deltaUPos * ptrack->GetTgl();

   Double_t convposz = (zphasePos*negtrackradius+zphaseNeg*postrackradius)/(negtrackradius+postrackradius);

   return convposz;
}

AliGammaConversionKFVector* AliV0Reader::GetBGGoodV0s(Int_t event){

  return fBGEventHandler->GetBGGoodV0s(event,fESDEvent->GetPrimaryVertex()->GetZ(),CountESDTracks());
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
  //  cout<<"Checking particle label, particle is: "<<fMCStack->Particle(TMath::Abs(label))->GetName()<<endl;
  if(fMCStack->Particle(TMath::Abs(label))->GetPdgCode() == 111){
    iResult=kTRUE;
  }
  return iResult;
}


Bool_t AliV0Reader::GetArmenterosQtAlfa(AliKFParticle* positiveKFParticle, AliKFParticle * negativeKFParticle, AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] ){
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


