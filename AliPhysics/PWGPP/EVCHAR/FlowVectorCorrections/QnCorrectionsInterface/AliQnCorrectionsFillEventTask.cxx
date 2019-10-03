/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
/***********************************************************
 Variable fill for Flow Qn vector corrections framework
 Based on work of Ionut-Cristian Arsene
 ***********************************************************/

#include <Riostream.h>

#include "AliQnCorrectionsFillEventTask.h"
#include "AliQnCorrectionsVarManagerTask.h"

#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsDetector.h"
#include "AliQnCorrectionsManager.h"

#include "AliQnCorrectionsHistos.h"

#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>

#include <AliInputEventHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliMultSelection.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODForwardMult.h>
#include <AliForwardUtil.h>

#include <AliVEvent.h>
#include <AliVZDC.h>
#include <AliVVZERO.h>
#include <AliVParticle.h>

#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include <AliESDFMD.h>

#include <AliAODInputHandler.h>
#include <AliAODEvent.h>
#include <AliAODHeader.h>
#include <AliAODTrack.h>

#include <AliLog.h>


ClassImp(AliQnCorrectionsFillEventTask)

const Float_t AliQnCorrectionsFillEventTask::fVZEROSignalThreshold = 0.01;
const Float_t AliQnCorrectionsFillEventTask::fTZEROSignalThreshold = 0.01;
const Float_t AliQnCorrectionsFillEventTask::fZDCSignalThreshold = 100.0;
const Float_t AliQnCorrectionsFillEventTask::fFMDSignalThreshold = 0.01;

AliQnCorrectionsFillEventTask::AliQnCorrectionsFillEventTask() :
AliQnCorrectionsVarManagerTask(),
fEvent(NULL),
fAliQnCorrectionsManager(NULL),
fEventHistos(NULL),
fDataBank(NULL),
fUseOnlyCentCalibEvents(kTRUE),
fUseTPCStandaloneTracks(kFALSE),
fFillVZERO(kFALSE),
fFillTPC(kFALSE),
fFillZDC(kFALSE),
fFillTZERO(kFALSE),
fFillFMD(kFALSE),
fFillRawFMD(kFALSE),
fFillSPD(kFALSE),
fIsAOD(kFALSE),
fIsESD(kFALSE)
{
  //
  // Default constructor
  //

}

AliQnCorrectionsFillEventTask::AliQnCorrectionsFillEventTask(const char *name) :
AliQnCorrectionsVarManagerTask(name),
fEvent(NULL),
fAliQnCorrectionsManager(NULL),
fEventHistos(NULL),
fDataBank(NULL),
fUseOnlyCentCalibEvents(kTRUE),
fUseTPCStandaloneTracks(kFALSE),
fFillVZERO(kFALSE),
fFillTPC(kFALSE),
fFillZDC(kFALSE),
fFillTZERO(kFALSE),
fFillFMD(kFALSE),
fFillRawFMD(kFALSE),
fFillSPD(kFALSE),
fIsAOD(kFALSE),
fIsESD(kFALSE)
{
  //
  // Default constructor
  //

}


//_____________________________________________________________________________
AliQnCorrectionsFillEventTask::~AliQnCorrectionsFillEventTask()
{
  //
  // Destructor
  //
}



//__________________________________________________________________
void AliQnCorrectionsFillEventTask::SetDetectors() {
  //
  // determine which detectors are used (to call only the necessary detector fill functions)
  if(fAliQnCorrectionsManager->FindDetector(kTPC  ) != NULL)    fFillTPC = kTRUE;
  if(fAliQnCorrectionsManager->FindDetector(kVZERO) != NULL)  fFillVZERO = kTRUE;
  if(fAliQnCorrectionsManager->FindDetector(kTZERO) != NULL)  fFillTZERO = kTRUE;
  if(fAliQnCorrectionsManager->FindDetector(kZDC  ) != NULL)    fFillZDC = kTRUE;
  if(fAliQnCorrectionsManager->FindDetector(kFMD  ) != NULL)    fFillFMD = kTRUE;
  if(fAliQnCorrectionsManager->FindDetector(kFMDraw) != NULL)fFillRawFMD = kTRUE;
  if(fAliQnCorrectionsManager->FindDetector(kSPD  ) != NULL)    fFillSPD = kTRUE;
}


//__________________________________________________________________
void AliQnCorrectionsFillEventTask::FillEventData() {

  TString aod = "AOD";
  TString esd = "ESD";

  fIsAOD = ( aod.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );
  fIsESD = ( esd.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );

  FillEventInfo();
  FillDetectors();
}

//__________________________________________________________________
void AliQnCorrectionsFillEventTask::FillEventInfo() {
  //
  // fill event info
  //


  fDataBank[kRunNo]       = fEvent->GetRunNumber();
  fDataBank[kVtxX]        = -999.;
  fDataBank[kVtxY]        = -999.;
  fDataBank[kVtxZ]        = -999.;
  const AliVVertex *primVtx = fEvent->GetPrimaryVertex();
  if (primVtx){
    fDataBank[kVtxX]        = primVtx->GetX();
    fDataBank[kVtxY]        = primVtx->GetY();
    fDataBank[kVtxZ]        = primVtx->GetZ();
    fDataBank[kNVtxContributors]    = primVtx->GetNContributors();
  }

  AliMultSelection *MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
  if(MultSelection) fDataBank[kVZEROMultPercentile] = MultSelection->GetMultiplicityPercentile("V0M", fUseOnlyCentCalibEvents);

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  AliCentrality* cent = esdEvent->GetCentrality();
  if(cent){
    fDataBank[kCentVZERO]   = cent->GetCentralityPercentile("V0M");
    fDataBank[kCentSPD]     = cent->GetCentralityPercentile("CL1");
    fDataBank[kCentTPC]     = cent->GetCentralityPercentile("TRK");
    fDataBank[kCentQuality] = cent->GetQuality();
  }


  AliVVZERO* vzero = fEvent->GetVZEROData();
  fDataBank[kVZEROATotalMult]     = vzero->GetMTotV0A();
  fDataBank[kVZEROCTotalMult]     = vzero->GetMTotV0C();
  fDataBank[kVZEROTotalMult]      = fDataBank[kVZEROATotalMult]+fDataBank[kVZEROCTotalMult];

  AliMultiplicity* spdmult = (AliMultiplicity*) fEvent->GetMultiplicity();
  fDataBank[kSPDntracklets]      = spdmult->GetNumberOfTracklets();
  fDataBank[kSPDnSingleClusters] = spdmult->GetNumberOfSingleClusters();
}



//__________________________________________________________________
void AliQnCorrectionsFillEventTask::FillTrackInfo(AliVParticle* particle) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;

  fDataBank[kPx]        = particle->Px();
  fDataBank[kPy]        = particle->Py();
  fDataBank[kPz]        = particle->Pz();
  fDataBank[kPt]        = particle->Pt();
  fDataBank[kP]         = particle->P();
  fDataBank[kPhi]       = particle->Phi();
  fDataBank[kTheta]     = particle->Theta();
  fDataBank[kEta]       = particle->Eta();
  fDataBank[kCharge]    = particle->Charge();
  fDataBank[kDcaXY]     = dcaxy;
  fDataBank[kDcaZ]      = dcaz;

  AliAODTrack* aodTrack=static_cast<AliAODTrack*>(particle);

  //fDataBank[VAR::kITSncls]       = particle->GetNcls(0);
  fDataBank[kTPCncls]       = aodTrack->GetTPCNcls();
  fDataBank[kTPCchi2]       = aodTrack->Chi2perNDF();
  fDataBank[kTPCsignal]     = aodTrack->GetTPCsignal();
  for(Int_t ibit=0; ibit<9; ibit++) fDataBank[kFilterBit+ibit]     = aodTrack->TestFilterBit(BIT(ibit));
  fDataBank[kFilterBitMask768]  = (aodTrack->TestFilterBit(BIT(8))||aodTrack->TestFilterBit(BIT(9)));

}


//__________________________________________________________________
void AliQnCorrectionsFillEventTask::FillTrackInfo(AliESDtrack* particle) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;
  particle->GetImpactParameters(dcaxy,dcaz);

  fDataBank[kPx]        = particle->Px();
  fDataBank[kPy]        = particle->Py();
  fDataBank[kPz]        = particle->Pz();
  fDataBank[kPt]        = particle->Pt();
  fDataBank[kP]         = particle->P();
  fDataBank[kPhi]       = particle->Phi();
  fDataBank[kTheta]     = particle->Theta();
  fDataBank[kEta]       = particle->Eta();
  fDataBank[kCharge]    = particle->Charge();
  fDataBank[kDcaXY]     = dcaxy;
  fDataBank[kDcaZ]      = dcaz;

  fDataBank[kTPCncls]       = particle->GetTPCNcls();
  fDataBank[kTPCnclsIter1]  = particle->GetTPCNclsIter1();
  fDataBank[kTPCchi2]       = fDataBank[kTPCncls]>0 ? particle->GetTPCchi2()/fDataBank[kTPCncls] : 0.0;
  fDataBank[kTPCchi2Iter1]  = fDataBank[kTPCnclsIter1]>0 ? particle->GetTPCchi2Iter1()/fDataBank[kTPCnclsIter1] : 0.0;
  fDataBank[kTPCsignal]     = particle->GetTPCsignal();



}

//_________________________________
void AliQnCorrectionsFillEventTask::FillDetectors(){

  if(fFillTPC)   FillTPC();
  if(fFillVZERO) FillVZERO();
  if(fFillZDC)   FillZDC();
  if(fFillTZERO) FillTZERO();
  if(fFillFMD)   FillFMD();
  if(fFillRawFMD)FillRawFMD();
  if(fFillSPD) FillSPDTracklets();
}


//_________________________________
void AliQnCorrectionsFillEventTask::FillTPC(){
  //
  // fill TPC info
  //

  if(fIsAOD) FillAodTPC();
  if(fIsESD) FillEsdTPC();
}


//_________________________________
void AliQnCorrectionsFillEventTask::FillAodTPC(){
  //
  // fill AOD TPC info
  //

  AliVParticle* vTrack;

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    vTrack = fEvent->GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if (!vTrack) continue;

    FillTrackInfo(vTrack);
    fEventHistos->FillHistClass("TrackQA_NoCuts", fDataBank);

    Int_t nNoOfAcceptedConf = fAliQnCorrectionsManager->AddDataVector(kTPC, vTrack->Phi());

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
        fEventHistos->FillHistClass(Form("TrackQA_%s",
            fAliQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(kTPC, conf)),
            fDataBank);
    }
  }
}




//_________________________________
void AliQnCorrectionsFillEventTask::FillEsdTPC(){
  //
  // fill ESD TPC info
  //

  AliESDtrack* esdTrack;

  const AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  const AliESDEvent& esd = *esdEvent;

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliESDtrack* track = NULL;
    esdTrack = esd.GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if(fUseTPCStandaloneTracks) track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(&esd),esdTrack->GetID());
    else track = esdTrack;
    if (!track) continue;

    FillTrackInfo(track);
    fEventHistos->FillHistClass("TrackQA_NoCuts", fDataBank);

    Int_t nNoOfAcceptedConf = fAliQnCorrectionsManager->AddDataVector(kTPC, track->Phi());

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
        fEventHistos->FillHistClass(Form("TrackQA_%s",
            fAliQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(kTPC, conf)),
            fDataBank);
    }

    if(fUseTPCStandaloneTracks) delete track;
  }
}



//_________________________________________________________________________________
void AliQnCorrectionsFillEventTask::FillSPDTracklets() {
  //
  // fill SPD info
  //

  Int_t nTracklets = 0;

  AliMultiplicity* mult = (AliMultiplicity*) fEvent->GetMultiplicity();
  nTracklets = mult->GetNumberOfTracklets();
  for(Int_t iTracklet=0; iTracklet<nTracklets; ++iTracklet) {
    fDataBank[kSPDtrackletEta]    = mult->GetEta(iTracklet);
    fDataBank[kSPDtrackletPhi]    = mult->GetPhi(iTracklet);

    Int_t nNoOfAcceptedConf = fAliQnCorrectionsManager->AddDataVector(kSPD, fDataBank[kSPDtrackletPhi]);

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
      fEventHistos->FillHistClass(Form("TrackletQA_%s",
          fAliQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(kSPD, conf)),
          fDataBank);
    }
  }
}

void AliQnCorrectionsFillEventTask::FillVZERO(){
  //
  // fill VZERO info
  //

  Double_t weight=0.;
  static const Double_t phi[8] = {1*TMath::Pi()/8.0, 3*TMath::Pi()/8.0, 5*TMath::Pi()/8.0, 7*TMath::Pi()/8.0,
      9*TMath::Pi()/8.0, 11*TMath::Pi()/8.0, 13*TMath::Pi()/8.0, 15*TMath::Pi()/8.0};

  AliVVZERO* vzero = fEvent->GetVZEROData();

  for(Int_t ich=0; ich<64; ich++){
    weight=vzero->GetMultiplicity(ich);
    if(weight > fVZEROSignalThreshold) {
      fAliQnCorrectionsManager->AddDataVector(kVZERO, phi[ich%8], weight, ich);   // 1st ich is position in array, 2nd ich is channel id
    }
  }
}



void AliQnCorrectionsFillEventTask::FillTZERO(){
  //
  // fill ESD TZERO info
  //

  Double_t weight=0.0;
  const Double_t kX[24] = {/* Cside */ 0.905348,0.571718,0.0848977,-0.424671,-0.82045,-0.99639,-0.905348,-0.571718,-0.0848977,0.424671,0.82045,0.99639,
                           /* Aside */ 0.99995,0.870982,0.508635,0.00999978,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635,-0.0100001,0.491315,0.860982};
  const Double_t kY[24] = {/* Cside */ 0.424671,0.82045,0.99639,0.905348,0.571718,0.0848976,-0.424671,-0.82045,-0.99639,-0.905348,-0.571719,-0.0848975,
                           /* Aside */ -0.00999983,0.491315,0.860982,0.99995,0.870982,0.508635,0.00999974,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635};
  const Double_t phi[24] = {TMath::ATan2(kY[0],kX[0]), TMath::ATan2(kY[1],kX[1]), TMath::ATan2(kY[2],kX[2]),
      TMath::ATan2(kY[3],kX[3]), TMath::ATan2(kY[4],kX[4]), TMath::ATan2(kY[5],kX[5]),
      TMath::ATan2(kY[6],kX[6]), TMath::ATan2(kY[7],kX[7]),
      TMath::ATan2(kY[8],kX[8]), TMath::ATan2(kY[9],kX[9]), TMath::ATan2(kY[10],kX[10]),
      TMath::ATan2(kY[11],kX[11]), TMath::ATan2(kY[12],kX[12]), TMath::ATan2(kY[13],kX[13]),
      TMath::ATan2(kY[14],kX[14]), TMath::ATan2(kY[15],kX[15]),
      TMath::ATan2(kY[16],kX[16]), TMath::ATan2(kY[17],kX[17]), TMath::ATan2(kY[18],kX[18]),
      TMath::ATan2(kY[19],kX[19]), TMath::ATan2(kY[20],kX[20]), TMath::ATan2(kY[21],kX[21]),
      TMath::ATan2(kY[22],kX[22]), TMath::ATan2(kY[23],kX[23]) };

  if (fIsESD) {
    const AliESDTZERO* esdT0 = dynamic_cast<AliESDEvent*>(fEvent)->GetESDTZERO();
    if (esdT0 != NULL) {
      for(Int_t ich=0; ich<24; ich++){
        weight=esdT0->GetT0amplitude()[ich];
        if(weight > fTZEROSignalThreshold) {
          fAliQnCorrectionsManager->AddDataVector(kTZERO, phi[ich], weight, ich);   // 1st ich is position in array, 2nd ich is channel id
        }
      }
    }
    else {
      AliError("AliESDTZERO not available");
    }
  }
  else {
    const AliAODTZERO* aodT0 = dynamic_cast<AliAODEvent*>(fEvent)->GetTZEROData();
    if (aodT0 != NULL) {
      for(Int_t ich=0; ich<24; ich++){
        weight=aodT0->GetAmp(ich);
        if(weight > fTZEROSignalThreshold) {
          fAliQnCorrectionsManager->AddDataVector(kTZERO, phi[ich], weight, ich);   // 1st ich is position in array, 2nd ich is channel id
        }
      }
    }
    else {
      AliError("AliAODTZERO not available");
    }
  }
}




//_________________________________
void AliQnCorrectionsFillEventTask::FillZDC(){
  //
  // fill ZDC info
  //


  Double_t weight=0.0;
  const Double_t kX[10] = { /* Cside */ 0.0,  -1.75,  1.75, -1.75, 1.75,
                            /* Aside */  0.0,  1.75, -1.75, 1.75, -1.75  };
  const Double_t kY[10] = { /* Cside */ 0.0,  -1.75, -1.75,  1.75, 1.75,
                            /* Aside */  0.0, -1.75, -1.75, 1.75,  1.75  };
  const Double_t phi[10] = {0.0, TMath::ATan2(kY[1],kX[1]), TMath::ATan2(kY[2],kX[2]),
      TMath::ATan2(kY[3],kX[3]), TMath::ATan2(kY[4],kX[4]), 0.0,
      TMath::ATan2(kY[6],kX[6]), TMath::ATan2(kY[7],kX[7]),
      TMath::ATan2(kY[8],kX[8]), TMath::ATan2(kY[9],kX[9]) };


  AliVZDC* zdc = (AliVZDC*) fEvent->GetZDCData();

  Double_t ZDCenergy[10];
  for(Int_t i=0; i<5; ++i)    ZDCenergy[i]  = zdc->GetZNCTowerEnergy()[i];
  for(Int_t i=5; i<10; ++i)   ZDCenergy[i]  = zdc->GetZNATowerEnergy()[i-5];

  for(Int_t ich=1; ich<10; ich++){
    if(ich==5) continue;
    weight=ZDCenergy[ich];
    if(weight > fZDCSignalThreshold) {
      fAliQnCorrectionsManager->AddDataVector(kZDC, phi[ich], weight, ich);   // 1st ich is position in array, 2nd ich is channel id
    }
  }
}


void AliQnCorrectionsFillEventTask::FillFMD()
{
  //
  // fill FMD info
  //

  Float_t m,phi;

  AliAODEvent* aodEvent=0x0;
  if(fIsAOD) aodEvent = (AliAODEvent*) fEvent;
  if(fIsESD) aodEvent = AliForwardUtil::GetAODEvent(this);


  if (!aodEvent) {
    AliFatal("Didn't get AOD event. Aborting! Check the AOD event handler presence.\n");
    return;
  }


  TObject* obj = aodEvent->FindListObject("Forward");
  if (!obj) {
    AliError("Didn't get the AOD Forward multiplicity object instance\n");
    return;
  }

  AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(obj);

  const TH2D& d2Ndetadphi = aodForward->GetHistogram();

  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();


  // Loop over eta 
  Int_t nFMD=-1;
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
    if (!valid) continue; // No data expected for this eta 

    // Loop over phi 
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
      m     =  d2Ndetadphi.GetBinContent(iEta, iPhi);
      if(m > fFMDSignalThreshold) {
        nFMD++;
        fAliQnCorrectionsManager->AddDataVector(kFMD, phi, m, iEta*nPhi+iPhi);   // 1st ich is position in array, 2nd ich is channel id
      }
    }
  }
}



//_________________________________
void AliQnCorrectionsFillEventTask::FillRawFMD()
{
  //
  // fill Raw FMD info
  //
  Bool_t isESD = (fEvent->IsA()==AliESDEvent::Class());
  if(!isESD) return;

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);

  AliESDFMD* esdFmd = esdEvent->GetFMDData();
/*
  Int_t id=-1;
  Int_t maxDet=3;
  Int_t maxRing=2;
  Int_t maxSector;
  Int_t maxStrip;
  Float_t m=0.0;
  Double_t phi,eta;
  Char_t ring;
*/
  Int_t nNoOfDetectors       = 3;           ///< the number of FMD detectors
  Int_t detectorNumber[]     = {1,2,3};     ///< the number of the FMD detector
  Int_t nNoOfRings[]         = {1,2,2};     ///< the number of rings for each detector
  Char_t ringId[]            = {'I','O'};   ///< ring identity
  Int_t ringNoOfSectors[]    = {20,40};   ///< ring number of sectors
  Int_t ringNoOfStrips[]     = {512,256};     ///< ring number of strips per sector
  Int_t nSectorId = 0;


  /**
   * Get the azimuthal angle of
   * @f$ \text{FMD}\langle detector\rangle\lange ring\rangle_{\langle
   * sector\rangle\langle strip\rangle}@f$
   *
   * @param detector Detector number (1-3)
   * @param ring     Ring identifier ('I' or 'O')
   * @param sector   Sector number (0-19, or 0-39)
   * @param strip    Strip number (0-511, or 0-255)
   *
   * @return Azimuthal angle
   */

  for(Int_t detector = 0; detector < nNoOfDetectors; detector++) {
    for(Int_t ring = 0; ring < nNoOfRings[detector]; ring++) {
      for(Int_t sector = 0; sector < ringNoOfSectors[ring]; sector++) {
        Double_t phi  =  esdFmd->Phi(detectorNumber[detector], ringId[ring], sector, 0) / 180. * TMath::Pi();
        for(Int_t strip = 0; strip < ringNoOfStrips[ring]; strip++) {
          Double_t eta  =  esdFmd->Eta(detectorNumber[detector], ringId[ring], sector, strip);
          Float_t m = esdFmd->Multiplicity(detectorNumber[detector], ringId[ring], sector, strip);
          if(m !=  AliESDFMD::kInvalidMult) {
            fDataBank[kFMDEta] = eta;
            fAliQnCorrectionsManager->AddDataVector(kFMDraw, phi, m, nSectorId);   // 1st ich is position in array, 2nd ich is channel id
          }
        }  // end loop over strips
        nSectorId++;
      }  // end loop over sectors      
    }  // end loop over rings
  } // end loop over detectors
}


