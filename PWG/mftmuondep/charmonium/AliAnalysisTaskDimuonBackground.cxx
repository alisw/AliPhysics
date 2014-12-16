#include "TH1D.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskDimuonBackground.h"
#include "TDatabasePDG.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliAODDimuon.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"
#include "TClonesArray.h"
#include "AliMFTConstants.h"
#include "AliMFTAnalysisTools.h"
#include "TRandom.h"
#include "TList.h"

ClassImp(AliAnalysisTaskDimuonBackground)

//====================================================================================================================================================

AliAnalysisTaskDimuonBackground::AliAnalysisTaskDimuonBackground() : 
  AliAnalysisTaskSE(),
  fVertexMode(0),
  fMinTriggerMatch(0),
  fSingleMuonMinEta(-9999), 
  fSingleMuonMaxEta(9999), 
  fSingleMuonMinPt(0),
  fSingleMuonMaxChi2(9999),
  fHistogramList(0),
  fMainInputHandler(0),
  fMixingInputHandler(0)
{
  
  // Default constructor

  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut = 0;
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut = 0;

  for (Int_t i=0; i<2; i++) fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] = 0;
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances = 0;
  
  for (Int_t i=0; i<3; i++) { fPrimaryVertex[i]=0; fPrimaryVertexMixEv[i]=0; fPrimaryVertexTrue[i]=0; fPrimaryVertexMixEvTrue[i]=0; }
  fVtxResolutionITS[0] = 5.e-4;
  fVtxResolutionITS[1] = 5.e-4;
  fVtxResolutionITS[2] = 4.e-4;

  fMassJpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();

}

//====================================================================================================================================================

AliAnalysisTaskDimuonBackground::AliAnalysisTaskDimuonBackground(const Char_t *name) : 
  AliAnalysisTaskSE(name),
  fVertexMode(0),
  fMinTriggerMatch(0),
  fSingleMuonMinEta(-9999), 
  fSingleMuonMaxEta(9999), 
  fSingleMuonMinPt(0),
  fSingleMuonMaxChi2(9999),
  fHistogramList(0),
  fMainInputHandler(0),
  fMixingInputHandler(0) 
{

  // Constructor

  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut = 0;
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut = 0;

  for (Int_t i=0; i<2; i++) fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] = 0;
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances = 0;

  for (Int_t i=0; i<3; i++) { fPrimaryVertex[i]=0; fPrimaryVertexMixEv[i]=0; fPrimaryVertexTrue[i]=0; fPrimaryVertexMixEvTrue[i]=0; }
  fVtxResolutionITS[0] = 5.e-4;
  fVtxResolutionITS[1] = 5.e-4;
  fVtxResolutionITS[2] = 4.e-4;

  fMassJpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();

  // Define input and output slots here
  DefineOutput(1, TList::Class());

}

//====================================================================================================================================================

void AliAnalysisTaskDimuonBackground::UserCreateOutputObjects() {

  // Called once

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);

  // Max 100 dimensions
  Int_t  nBins[100] = {0};
  Double_t min[100] = {0};
  Double_t max[100] = {0};

  TString tag[2] = {"SingleEvents", "MixedEvents"};

  //-------------------------------------------------------------------------------

  nBins[0] = 100;  
  nBins[1] = 1000;
  nBins[2] = 1000;
  nBins[3] = 100;
  nBins[4] = 100;

  min[0] = 0.;
  min[1] = -10.;
  min[2] = 0.;
  min[3] = 0.;
  min[4] = -4.5;

  max[0] = 1.;
  max[1] = 10.;
  max[2] = 10.;
  max[3] = 10.;
  max[4] = -2.0;

  for (Int_t i=0; i<2; i++) {
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] = new THnSparseD(Form("fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_%s",tag[i].Data()), "", 5, nBins, min, max);
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] -> Sumw2(); 
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] -> GetAxis(0)->SetTitle("PCA Quality");
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] -> GetAxis(1)->SetTitle("t_{z}  [ps]");
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] -> GetAxis(2)->SetTitle("Mass  [GeV/c^{2}]");
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] -> GetAxis(3)->SetTitle("p_{T}  [GeV/c]");
    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i] -> GetAxis(4)->SetTitle("Rapidity");
    fHistogramList->Add(fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[i]);
  }

  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances = new THnSparseD("fHistPCAQualityVsPPDecayTimeVsMassVsPt_Resonances", "", 5, nBins, min, max);
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances -> Sumw2(); 
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances -> GetAxis(0)->SetTitle("PCA Quality");
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances -> GetAxis(1)->SetTitle("t_{z}  [ps]");
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances -> GetAxis(2)->SetTitle("Mass  [GeV/c^{2}]");
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances -> GetAxis(3)->SetTitle("p_{T}  [GeV/c]");
  fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances -> GetAxis(4)->SetTitle("Rapidity");
  fHistogramList->Add(fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances);

  nBins[0] = 200;  
  nBins[1] = 100;
  nBins[2] = 100;

  min[0] = 0.;
  min[1] = 0.;
  min[2] = -4.5;

  max[0] = 10.;
  max[1] = 10.;
  max[2] = -2.0;
  
  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut = new THnSparseD("fHistSingleMuonsChi2VsPtVsRapidityBeforeCut", "", 3, nBins, min, max);
  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut -> Sumw2(); 
  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut -> GetAxis(0)->SetTitle("p_{T}  [GeV/c]");
  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut -> GetAxis(1)->SetTitle("#chi^{2}");
  fHistSingleMuonsChi2VsPtVsRapidityBeforeCut -> GetAxis(2)->SetTitle("Rapidity");
  fHistogramList->Add(fHistSingleMuonsChi2VsPtVsRapidityBeforeCut);
  
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut = new THnSparseD("fHistSingleMuonsChi2VsPtVsRapidityAfterCut", "", 3, nBins, min, max);
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut -> Sumw2(); 
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut -> GetAxis(0)->SetTitle("p_{T}  [GeV/c]");
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut -> GetAxis(1)->SetTitle("#chi^{2}");
  fHistSingleMuonsChi2VsPtVsRapidityAfterCut -> GetAxis(2)->SetTitle("Rapidity");
  fHistogramList->Add(fHistSingleMuonsChi2VsPtVsRapidityAfterCut);

  //-------------------------------------------------------------------------------

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  SetMainInputHandler(mgr);
  if (fMainInputHandler) SetMixingInputHandler(fMainInputHandler);

  PostData(1, fHistogramList);

}

//====================================================================================================================================================

void AliAnalysisTaskDimuonBackground::UserExec(Option_t *) {

  // Main loop
  // Called for each event

  AliAODEvent *aodEv = dynamic_cast<AliAODEvent*>(GetMainEvent());

  if (!aodEv) return;

  if (!(AliMUONTrackExtrap::IsFieldON())) {
    ((AliAODHeader*) aodEv->GetHeader())->InitMagneticField();
    AliMUONTrackExtrap::SetField();
  }

  // Getting primary vertex, either from the generation or from the reconstruction -------------------

  AliAODMCHeader *mcHeader = (AliAODMCHeader*) (aodEv->GetList()->FindObject(AliAODMCHeader::StdBranchName()));

  if (fVertexMode == kGenerated) {
    mcHeader->GetVertex(fPrimaryVertex);
    for (Int_t i=0; i<3; i++) fPrimaryVertex[i] = gRandom->Gaus(fPrimaryVertex[i], fVtxResolutionITS[i]);
  }
  else if (fVertexMode == kReconstructed) {
    aodEv->GetPrimaryVertex()->GetXYZ(fPrimaryVertex);
  }

  TClonesArray *stackMC = (TClonesArray*) (aodEv->GetList()->FindObject(AliAODMCParticle::StdBranchName()));

  AliAODTrack *recMuon[2] = {0};
  AliAODMCParticle *mcMuon[2]={0}, *mcMother[2]={0};
  Double_t var[100] = {0};

  //--------------------------------------------------------------------------------

  for (Int_t iTrack=0; iTrack<aodEv->GetNumberOfTracks(); iTrack++) { 

    recMuon[0] = (AliAODTrack*) aodEv->GetTrack(iTrack);
    if (!(recMuon[0]->IsMuonGlobalTrack()))                                  continue;
    if (AliMFTAnalysisTools::IsTrackInjected(recMuon[0], mcHeader, stackMC)) continue;    // Only HIJING tracks are considered to build the background

    var[0] = recMuon[0]->Pt();
    var[1] = recMuon[0]->Chi2perNDF();
    var[2] = recMuon[0]->Y();

    fHistSingleMuonsChi2VsPtVsRapidityBeforeCut -> Fill(var);

    if (!IsSingleMuonCutPassed(recMuon[0]))                                  continue;

    fHistSingleMuonsChi2VsPtVsRapidityAfterCut -> Fill(var);

    mcMuon[0]   = NULL;
    mcMother[0] = NULL;

    if (recMuon[0]->GetLabel()>=0)              mcMuon[0]   = (AliAODMCParticle*) stackMC->At(recMuon[0]->GetLabel());
    if (mcMuon[0] && mcMuon[0]->GetMother()>=0) mcMother[0] = (AliAODMCParticle*) stackMC->At(mcMuon[0]->GetMother());

    for (Int_t jTrack=0; jTrack<iTrack; jTrack++) { 

      recMuon[1] = (AliAODTrack*) aodEv->GetTrack(jTrack);
      if (!(recMuon[1]->IsMuonGlobalTrack()))                                  continue;
      if (AliMFTAnalysisTools::IsTrackInjected(recMuon[1], mcHeader, stackMC)) continue;    // Only HIJING tracks are considered to build the background
      if (!IsSingleMuonCutPassed(recMuon[1]))                                  continue;

      mcMuon[1]   = NULL;
      mcMother[1] = NULL;
      
      if (recMuon[1]->GetLabel()>=0)              mcMuon[1]   = (AliAODMCParticle*) stackMC->At(recMuon[1]->GetLabel());
      if (mcMuon[1] && mcMuon[1]->GetMother()>=0) mcMother[1] = (AliAODMCParticle*) stackMC->At(mcMuon[1]->GetMother());
      
      TObjArray *dimuon = new TObjArray();
      dimuon -> Add(recMuon[0]);
      dimuon -> Add(recMuon[1]);	
      
      if (recMuon[0]->Charge() + recMuon[1]->Charge()) {
	delete dimuon;
	continue;
      }

      Bool_t isDimuonResonance = kFALSE;
      if (mcMother[0] && mcMother[1] && mcMuon[0]->GetMother()==mcMuon[1]->GetMother()) {
	if (AliMFTAnalysisTools::IsPDGResonance(mcMother[0]->GetPdgCode())) isDimuonResonance = kTRUE;
      }

      Double_t pca[3]={0};
      Double_t pcaQuality=0;
      TLorentzVector kinem(0,0,0,0);
      if (!AliMFTAnalysisTools::CalculatePCA(dimuon, pca, pcaQuality, kinem)) {
	delete dimuon;
	continue;
      }

      Double_t ppDecayTime = AliMFTAnalysisTools::GetPseudoProperDecayTimeZ(fPrimaryVertex[2], pca[2], fMassJpsi, kinem.Pz());

      var[0] = pcaQuality;
      var[1] = ppDecayTime;
      var[2] = kinem.M();
      var[3] = kinem.Pt();
      var[4] = kinem.Rapidity();

      if (!isDimuonResonance) fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[kSingleEvents] -> Fill(var);
      else                    fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity_Resonances     -> Fill(var);

      delete dimuon;

    }   // end of loop on 2nd muon

  }   // end of loop on 1st muon

  //--------------------------------------------------------------------------------
  
  PostData(1, fHistogramList);

}

//====================================================================================================================================================

void AliAnalysisTaskDimuonBackground::UserExecMix(Option_t *) {

  // Performing Mixing

  if (!fMixingInputHandler) return;

  Int_t bufferSize  = fMixingInputHandler->BufferSize();

  AliAODEvent *aodEv = dynamic_cast<AliAODEvent*>(GetMainEvent());

  if (!aodEv) return;

  if (!(AliMUONTrackExtrap::IsFieldON())) {
    ((AliAODHeader*) aodEv->GetHeader())->InitMagneticField();
    AliMUONTrackExtrap::SetField();
  }

  // Getting primary vertex, either from the generation or from the reconstruction -------------------

  AliAODMCHeader *mcHeader = (AliAODMCHeader*) (aodEv->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  mcHeader->GetVertex(fPrimaryVertexTrue);

  if (fVertexMode == kGenerated) {
    mcHeader->GetVertex(fPrimaryVertex);
    for (Int_t i=0; i<3; i++) fPrimaryVertex[i] = gRandom->Gaus(fPrimaryVertex[i], fVtxResolutionITS[i]);
  }
  else if (fVertexMode == kReconstructed) {
    aodEv->GetPrimaryVertex()->GetXYZ(fPrimaryVertex);
  }

  TClonesArray *stackMC = (TClonesArray*) (aodEv->GetList()->FindObject(AliAODMCParticle::StdBranchName()));

  // -------------------------------------------------------------------------------------------------

  AliAODTrack *recMuon[2] = {0};
  Double_t var[100] = {0};

  for (Int_t iBuffer=0; iBuffer<bufferSize; iBuffer++) {

    AliAODEvent *aodEvMix = dynamic_cast<AliAODEvent *>(GetMixedEvent(iBuffer));
    if (!aodEvMix) continue;

    // Getting primary vertex, either from the generation or from the reconstruction -------------------
    
    AliAODMCHeader *mcHeaderMixEv = (AliAODMCHeader*) (aodEvMix->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    mcHeaderMixEv->GetVertex(fPrimaryVertexMixEvTrue);
 
    if (fVertexMode == kGenerated) {
      mcHeaderMixEv->GetVertex(fPrimaryVertexMixEv);
      for (Int_t i=0; i<3; i++) fPrimaryVertexMixEv[i] = gRandom->Gaus(fPrimaryVertexMixEv[i], fVtxResolutionITS[i]);
    }
    else if (fVertexMode == kReconstructed) {
      aodEvMix->GetPrimaryVertex()->GetXYZ(fPrimaryVertexMixEv);
    }

    TClonesArray *stackMCMixEv = (TClonesArray*) (aodEvMix->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  
    // Loop over MUON+MFT muons of 1st and 2nd events
    
    //--------------------------------------------------------------------------------------------------
    
    for (Int_t iTrack=0; iTrack<aodEv->GetNumberOfTracks(); iTrack++) { 
      
      recMuon[0] = (AliAODTrack*) aodEv->GetTrack(iTrack);
      if (!(recMuon[0]->IsMuonGlobalTrack()))                                  continue;
      if (AliMFTAnalysisTools::IsTrackInjected(recMuon[0], mcHeader, stackMC)) continue;    // Only HIJING tracks are considered to build the background
      if (!IsSingleMuonCutPassed(recMuon[0]))                                  continue;

      for (Int_t jTrack=0; jTrack<aodEvMix->GetNumberOfTracks(); jTrack++) { 

	recMuon[1] = (AliAODTrack*) aodEvMix->GetTrack(jTrack);
	if (!(recMuon[1]->IsMuonGlobalTrack()))                                            continue;
	if (AliMFTAnalysisTools::IsTrackInjected(recMuon[1], mcHeaderMixEv, stackMCMixEv)) continue;    // Only HIJING tracks are considered to build the background
	if (!IsSingleMuonCutPassed(recMuon[1]))                                            continue;

	// ----------- Translating muons to a common primary vertex (the origin). Preserve original tracks in case one wants to use them

	AliAODTrack *recMuonTranslated[2] = {0};
	for (Int_t i=0; i<2; i++) recMuonTranslated[i] = new AliAODTrack(*(recMuon[i]));

	if (!(AliMFTAnalysisTools::TranslateMuonToOrigin(recMuonTranslated[0], fPrimaryVertexTrue)) ||
	    !(AliMFTAnalysisTools::TranslateMuonToOrigin(recMuonTranslated[1], fPrimaryVertexMixEvTrue))) {
	  printf("Error: tracks cannot be translated!!!\n");
	  delete recMuonTranslated[0];
	  delete recMuonTranslated[1];
	  continue;
	}

	TObjArray *dimuon = new TObjArray();
	dimuon -> Add(recMuonTranslated[0]);
	dimuon -> Add(recMuonTranslated[1]);	

	if (recMuonTranslated[0]->Charge() + recMuonTranslated[1]->Charge()) {
	  delete dimuon;
	  delete recMuonTranslated[0];
	  delete recMuonTranslated[1];
	  continue;
	}

	Double_t pca[3]={0};
	Double_t pcaQuality=0;
	TLorentzVector kinem(0,0,0,0);
	if (!AliMFTAnalysisTools::CalculatePCA(dimuon, pca, pcaQuality, kinem)) {
	  delete dimuon;
	  delete recMuonTranslated[0];
	  delete recMuonTranslated[1];
	  continue;
	}
	
	Double_t ppDecayTime = AliMFTAnalysisTools::GetPseudoProperDecayTimeZ(fPrimaryVertex[2]-fPrimaryVertexTrue[2], pca[2], fMassJpsi, kinem.Pz());
	
	var[0] = pcaQuality;
	var[1] = ppDecayTime;
	var[2] = kinem.M();
	var[3] = kinem.Pt();
	var[4] = kinem.Rapidity();

	fHistPCAQualityVsPPDecayTimeVsMassVsPtVsRapidity[kMixedEvents] -> Fill(var);
	
	delete dimuon;
	delete recMuonTranslated[0];
	delete recMuonTranslated[1];
	
      }   // end of loop on 2nd muon
      
    }   // end of loop on 1st muon
    
    //--------------------------------------------------------------------------------------------------
   
  }   // end of loop on events in buffer

  PostData(1, fHistogramList);

}

//====================================================================================================================================================

void AliAnalysisTaskDimuonBackground::Terminate(Option_t *) {

  // Draw result to the screen
  // Called once at the end of the query
  
}

//====================================================================================================================================================

AliVEvent *AliAnalysisTaskDimuonBackground::GetMainEvent() {

  // Access to MainEvent
  
  AliMultiInputEventHandler *inEvHMainMulti = fMainInputHandler;
  if (inEvHMainMulti) {
    AliInputEventHandler *inEvMain = dynamic_cast<AliInputEventHandler *>(inEvHMainMulti->GetFirstInputEventHandler());
    if (inEvMain) return inEvMain->GetEvent();
  }
  
  return 0;

}

//====================================================================================================================================================

AliVEvent *AliAnalysisTaskDimuonBackground::GetMixedEvent(Int_t buffId) {

  // Access to Mixed event with buffer id
  
  AliMultiInputEventHandler *inEvHMain = fMainInputHandler;
  
  if (inEvHMain) {
    
    AliMixInputEventHandler *mixIH = fMixingInputHandler;
    if (!mixIH) return 0;

    if (mixIH->CurrentBinIndex() < 0) {
      AliDebug(AliLog::kDebug + 1, "Current event mixEH->CurrentEntry() == -1");
      return 0;
    }
    
    AliMultiInputEventHandler *inEvHMixedCurrent = dynamic_cast<AliMultiInputEventHandler *>(mixIH->InputEventHandler(buffId));
    if (!inEvHMixedCurrent) return 0;

    AliInputEventHandler *ihMixedCurrent = inEvHMixedCurrent->GetFirstInputEventHandler();
    if (ihMixedCurrent) return ihMixedCurrent->GetEvent();

  }
  
  return 0;

}

//====================================================================================================================================================

AliMultiInputEventHandler *AliAnalysisTaskDimuonBackground::SetMainInputHandler(AliAnalysisManager *mgr) {
  
  // Sets main input handler
  
  if (!fMainInputHandler) fMainInputHandler = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
  
  return fMainInputHandler;

}

//====================================================================================================================================================

AliMixInputEventHandler *AliAnalysisTaskDimuonBackground::SetMixingInputHandler(AliMultiInputEventHandler *mainIH) {

  // Sets mixing input handler
  
  if (!fMixingInputHandler) fMixingInputHandler = dynamic_cast<AliMixInputEventHandler *>(mainIH->GetFirstMultiInputHandler());
  
  return fMixingInputHandler;

}

//====================================================================================================================================================

Bool_t AliAnalysisTaskDimuonBackground::IsSingleMuonCutPassed(AliAODTrack *mu) {

  if (mu->GetMatchTrigger() < fMinTriggerMatch)                   return kFALSE;
  if (mu->Eta()<fSingleMuonMinEta || mu->Eta()>fSingleMuonMaxEta) return kFALSE;
  if (mu->Pt() < fSingleMuonMinPt)                                return kFALSE;
  if (mu->Chi2perNDF() > fSingleMuonMaxChi2)                      return kFALSE;

  return kTRUE;

}

//====================================================================================================================================================

