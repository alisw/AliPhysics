/**************************************************************************
 *  Authors : Domenico Colella                                            *
 *                                                                        *
 *                                                                        *
 * Derived from the:                                                      *
 *  - Original AliAnalysisTaskCheckCascade (A. Maire, G. Hippolyte)       *
 *  - Adapted to PbPb analysis (M. Nicassio)                              *
 *  - Adapted to work on all collisidng systems (pp, PbPb and pPb),       *
 *    ESD/AOD and experimental/MC data (D. Colella)                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;


#include <Riostream.h>
#include "THnSparse.h"
#include "TVector3.h"
#include "TMath.h"

#include "AliLog.h"
#include "AliCentrality.h"
#include "AliHeader.h"   //for MC
#include "AliMCEvent.h"  //for MC
#include "AliStack.h"    //for MC
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h" 
#include "AliAODInputHandler.h"
#include "AliMultiplicity.h"

#include "AliGenEventHeader.h"         //for MC
#include "AliGenCocktailEventHeader.h" //for MC
#include "AliGenHijingEventHeader.h"   //for MC
#include "AliAODMCParticle.h"          //for MC

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDpid.h"

#include "AliMultiplicity.h"
#include "AliMultSelection.h"

#include "AliAnalysisTaskQAMultistrange.h"

ClassImp(AliAnalysisTaskQAMultistrange)



//------------------------------------------------------------
AliAnalysisTaskQAMultistrange::AliAnalysisTaskQAMultistrange() 
  : AliAnalysisTaskSE           (), 
    fisMC                       (kFALSE),
    fAnalysisType               ("ESD"), 
    fPIDResponse                (0),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCutnTPCcls         (kTRUE),
    fMinnTPCcls                 (70),  
    fMinPtCutOnDaughterTracks   (0.),
    fEtaCutOnDaughterTracks     (0.8),


    fListHistMultistrangeQA(0),
      fHistEventSel(0),
      fHistCascadeMultiplicityXiMinus(0), fHistCascadeMultiplicityXiPlus(0), fHistCascadeMultiplicityOmegaMinus(0), fHistCascadeMultiplicityOmegaPlus(0),
      fHistVarDcaCascDaughtXiMinus(0), fHistVarDcaCascDaughtXiPlus(0), fHistVarDcaCascDaughtOmegaMinus(0), fHistVarDcaCascDaughtOmegaPlus(0),
      fHistVarDcaBachToPrimVertexXiMinus(0), fHistVarDcaBachToPrimVertexXiPlus(0), fHistVarDcaBachToPrimVertexOmegaMinus(0), fHistVarDcaBachToPrimVertexOmegaPlus(0),
      fHistVarCascCosineOfPointingAngleXiMinus(0), fHistVarCascCosineOfPointingAngleXiPlus(0), fHistVarCascCosineOfPointingAngleOmegaMinus(0), fHistVarCascCosineOfPointingAngleOmegaPlus(0),
      fHistVarCascRadiusXiMinus(0), fHistVarCascRadiusXiPlus(0), fHistVarCascRadiusOmegaMinus(0), fHistVarCascRadiusOmegaPlus(0),
      fHistVarInvMassLambdaAsCascDghterXiMinus(0), fHistVarInvMassLambdaAsCascDghterXiPlus(0), fHistVarInvMassLambdaAsCascDghterOmegaMinus(0), fHistVarInvMassLambdaAsCascDghterOmegaPlus(0),
      fHistVarDcaV0DaughtersXiMinus(0), fHistVarDcaV0DaughtersXiPlus(0), fHistVarDcaV0DaughtersOmegaMinus(0), fHistVarDcaV0DaughtersOmegaPlus(0),
      fHistVarV0CosineOfPAToCascVertexXiMinus(0), fHistVarV0CosineOfPAToCascVertexXiPlus(0), fHistVarV0CosineOfPAToCascVertexOmegaMinus(0), fHistVarV0CosineOfPAToCascVertexOmegaPlus(0),
      fHistVarV0RadiusXiMinus(0), fHistVarV0RadiusXiPlus(0), fHistVarV0RadiusOmegaMinus(0), fHistVarV0RadiusOmegaPlus(0),
      fHistVarDcaV0ToPrimVertexXiMinus(0), fHistVarDcaV0ToPrimVertexXiPlus(0), fHistVarDcaV0ToPrimVertexOmegaMinus(0), fHistVarDcaV0ToPrimVertexOmegaPlus(0),
      fHistVarDcaPosToPrimVertexXiMinus(0), fHistVarDcaPosToPrimVertexXiPlus(0), fHistVarDcaPosToPrimVertexOmegaMinus(0), fHistVarDcaPosToPrimVertexOmegaPlus(0),
      fHistVarDcaNegToPrimVertexXiMinus(0), fHistVarDcaNegToPrimVertexXiPlus(0), fHistVarDcaNegToPrimVertexOmegaMinus(0), fHistVarDcaNegToPrimVertexOmegaPlus(0),
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      fHistVarTransvMomentumXiMinus(0), fHistVarTransvMomentumXiPlus(0), fHistVarTransvMomentumOmegaMinus(0), fHistVarTransvMomentumOmegaPlus(0),
      fHistVarRapidityXiMinus(0), fHistVarRapidityXiPlus(0), fHistVarRapidityOmegaMinus(0), fHistVarRapidityOmegaPlus(0),
      fHistVarCascProperLengthXiMinus(0), fHistVarCascProperLengthXiPlus(0), fHistVarCascProperLengthOmegaMinus(0), fHistVarCascProperLengthOmegaPlus(0),
      fHistVarV0ProperLengthXiMinus(0), fHistVarV0ProperLengthXiPlus(0), fHistVarV0ProperLengthOmegaMinus(0), fHistVarV0ProperLengthOmegaPlus(0),
      fHistGenVarTotMomXiMinus(0),fHistGenVarTotMomXiPlus(0), fHistGenVarTotMomOmegaMinus(0), fHistGenVarTotMomOmegaPlus(0),
      fHistGenVarTransvMomXiMinus(0), fHistGenVarTransvMomXiPlus(0), fHistGenVarTransvMomOmegaMinus (0), fHistGenVarTransvMomOmegaPlus(0),
      fHistGenVarYXiMinus(0), fHistGenVarYXiPlus(0), fHistGenVarYOmegaMinus(0), fHistGenVarYOmegaPlus(0),
      fHistGenVarEtaXiMinus(0), fHistGenVarEtaXiPlus(0), fHistGenVarEtaOmegaMinus(0), fHistGenVarEtaOmegaPlus(0),
      fHistGenVarThetaXiMinus(0), fHistGenVarThetaXiPlus(0), fHistGenVarThetaOmegaMinus(0), fHistGenVarThetaOmegaPlus(0),
      fHistGenVarPhiXiMinus(0), fHistGenVarPhiXiPlus(0), fHistGenVarPhiOmegaMinus(0), fHistGenVarPhiOmegaPlus(0)



{
  // Dummy Constructor
}


//----------------------------------------------------------------------------
AliAnalysisTaskQAMultistrange::AliAnalysisTaskQAMultistrange(const char *name) 
  : AliAnalysisTaskSE           (name),
    fisMC                       (kFALSE),
    fAnalysisType               ("ESD"), 
    fPIDResponse                (0),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCutnTPCcls         (kTRUE),
    fMinnTPCcls                 (70),
    fMinPtCutOnDaughterTracks   (0.),
    fEtaCutOnDaughterTracks     (0.8),


    fListHistMultistrangeQA(0),
      fHistEventSel(0),
      fHistCascadeMultiplicityXiMinus(0), fHistCascadeMultiplicityXiPlus(0), fHistCascadeMultiplicityOmegaMinus(0), fHistCascadeMultiplicityOmegaPlus(0),
      fHistVarDcaCascDaughtXiMinus(0), fHistVarDcaCascDaughtXiPlus(0), fHistVarDcaCascDaughtOmegaMinus(0), fHistVarDcaCascDaughtOmegaPlus(0),
      fHistVarDcaBachToPrimVertexXiMinus(0), fHistVarDcaBachToPrimVertexXiPlus(0), fHistVarDcaBachToPrimVertexOmegaMinus(0), fHistVarDcaBachToPrimVertexOmegaPlus(0),
      fHistVarCascCosineOfPointingAngleXiMinus(0), fHistVarCascCosineOfPointingAngleXiPlus(0), fHistVarCascCosineOfPointingAngleOmegaMinus(0), fHistVarCascCosineOfPointingAngleOmegaPlus(0),
      fHistVarCascRadiusXiMinus(0), fHistVarCascRadiusXiPlus(0), fHistVarCascRadiusOmegaMinus(0), fHistVarCascRadiusOmegaPlus(0),
      fHistVarInvMassLambdaAsCascDghterXiMinus(0), fHistVarInvMassLambdaAsCascDghterXiPlus(0), fHistVarInvMassLambdaAsCascDghterOmegaMinus(0), fHistVarInvMassLambdaAsCascDghterOmegaPlus(0),
      fHistVarDcaV0DaughtersXiMinus(0), fHistVarDcaV0DaughtersXiPlus(0), fHistVarDcaV0DaughtersOmegaMinus(0), fHistVarDcaV0DaughtersOmegaPlus(0),
      fHistVarV0CosineOfPAToCascVertexXiMinus(0), fHistVarV0CosineOfPAToCascVertexXiPlus(0), fHistVarV0CosineOfPAToCascVertexOmegaMinus(0), fHistVarV0CosineOfPAToCascVertexOmegaPlus(0),
      fHistVarV0RadiusXiMinus(0), fHistVarV0RadiusXiPlus(0), fHistVarV0RadiusOmegaMinus(0), fHistVarV0RadiusOmegaPlus(0),
      fHistVarDcaV0ToPrimVertexXiMinus(0), fHistVarDcaV0ToPrimVertexXiPlus(0), fHistVarDcaV0ToPrimVertexOmegaMinus(0), fHistVarDcaV0ToPrimVertexOmegaPlus(0),
      fHistVarDcaPosToPrimVertexXiMinus(0), fHistVarDcaPosToPrimVertexXiPlus(0), fHistVarDcaPosToPrimVertexOmegaMinus(0), fHistVarDcaPosToPrimVertexOmegaPlus(0),
      fHistVarDcaNegToPrimVertexXiMinus(0), fHistVarDcaNegToPrimVertexXiPlus(0), fHistVarDcaNegToPrimVertexOmegaMinus(0), fHistVarDcaNegToPrimVertexOmegaPlus(0),
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      fHistVarTransvMomentumXiMinus(0), fHistVarTransvMomentumXiPlus(0), fHistVarTransvMomentumOmegaMinus(0), fHistVarTransvMomentumOmegaPlus(0),
      fHistVarRapidityXiMinus(0), fHistVarRapidityXiPlus(0), fHistVarRapidityOmegaMinus(0), fHistVarRapidityOmegaPlus(0),
      fHistVarCascProperLengthXiMinus(0), fHistVarCascProperLengthXiPlus(0), fHistVarCascProperLengthOmegaMinus(0), fHistVarCascProperLengthOmegaPlus(0),
      fHistVarV0ProperLengthXiMinus(0), fHistVarV0ProperLengthXiPlus(0), fHistVarV0ProperLengthOmegaMinus(0), fHistVarV0ProperLengthOmegaPlus(0),
      fHistGenVarTotMomXiMinus(0),fHistGenVarTotMomXiPlus(0), fHistGenVarTotMomOmegaMinus(0), fHistGenVarTotMomOmegaPlus(0),
      fHistGenVarTransvMomXiMinus(0), fHistGenVarTransvMomXiPlus(0), fHistGenVarTransvMomOmegaMinus (0), fHistGenVarTransvMomOmegaPlus(0),
      fHistGenVarYXiMinus(0), fHistGenVarYXiPlus(0), fHistGenVarYOmegaMinus(0), fHistGenVarYOmegaPlus(0),
      fHistGenVarEtaXiMinus(0), fHistGenVarEtaXiPlus(0), fHistGenVarEtaOmegaMinus(0), fHistGenVarEtaOmegaPlus(0),
      fHistGenVarThetaXiMinus(0), fHistGenVarThetaXiPlus(0), fHistGenVarThetaOmegaMinus(0), fHistGenVarThetaOmegaPlus(0),
      fHistGenVarPhiXiMinus(0), fHistGenVarPhiXiPlus(0), fHistGenVarPhiOmegaMinus(0), fHistGenVarPhiOmegaPlus(0)



{
  // Constructor
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());

  AliLog::SetClassDebugLevel("AliAnalysisTaskQAMultistrange",1);
}


//-------------------------------------------------------------
AliAnalysisTaskQAMultistrange::~AliAnalysisTaskQAMultistrange() 

{
  //
  // Destructor
  //
  // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
  // They will be deleted when fListCascade is deleted by the TSelector dtor
  // Because of TList::SetOwner() ...
  if (fListHistMultistrangeQA && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fListHistMultistrangeQA; fListHistMultistrangeQA = 0x0; }
}



//-----------------------------------------------------------
void AliAnalysisTaskQAMultistrange::UserCreateOutputObjects() 

{
  // Create histograms
  // Called once

  //____________________________________________________________________________
  // Option for AliLog: to suppress the extensive info prompted by a run with MC
  if (fisMC) AliLog::SetGlobalLogLevel(AliLog::kError);

  //_______________________________________________
  // Particle Identification Setup (new PID object)
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  //_________________
  // Define the TList
  fListHistMultistrangeQA = new TList();
  fListHistMultistrangeQA->SetOwner();

  //__________________
  // Define the Histos
  // -- Event selection distribution
  if (! fHistEventSel) {
        fHistEventSel = new TH1F("fHistEventSel", "Event selection;Evt. Sel. Step;Count",2, 0, 2);
        fHistEventSel->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventSel->GetXaxis()->SetBinLabel(2, "Selected");
        fListHistMultistrangeQA->Add(fHistEventSel);
  }
  // -- Cascade multiplicity distributions (for MC generated particles)
  if (! fHistCascadeMultiplicityXiMinus) {
     fHistCascadeMultiplicityXiMinus = new TH1F("fHistCascadeMultiplicityMCXiMinus","Xi Minus per event;Nbr of Xi Minus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityXiMinus);
  } 
  if (! fHistCascadeMultiplicityXiPlus) {
     fHistCascadeMultiplicityXiPlus = new TH1F("fHistCascadeMultiplicityMCXiPlus","Xi Plus per event;Nbr of Xi Plus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityXiPlus);
  }
  if (! fHistCascadeMultiplicityOmegaMinus) {
     fHistCascadeMultiplicityOmegaMinus = new TH1F("fHistCascadeMultiplicityMCOmegaMinus","Omega Minus per event;Nbr of Omega Minus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityOmegaMinus);
  }
  if (! fHistCascadeMultiplicityOmegaPlus) {
     fHistCascadeMultiplicityOmegaPlus = new TH1F("fHistCascadeMultiplicityMCOmegaPlus","Omega Plus per event;Nbr of Omega Plus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityOmegaPlus);
  }
  // -- Cascade variable distributions for reconstructed particles (in case of MC the association is required)
  // --- DcaCascDaught  (25 bins, [0.0,2.5])
  if (!fHistVarDcaCascDaughtXiMinus) {
     fHistVarDcaCascDaughtXiMinus = new TH2F("fHistVarDcaCascDaughtXiMinus",";#it{p_{T}} (GeV/#it{c});DCA cascade daughters (cm)", 250, 0.0, 25.0, 25, 0.0, 2.5);
     fListHistMultistrangeQA->Add(fHistVarDcaCascDaughtXiMinus);
  }
  if (!fHistVarDcaCascDaughtXiPlus) {
     fHistVarDcaCascDaughtXiPlus = new TH2F("fHistVarDcaCascDaughtXiPlus",";#it{p_{T}} (GeV/#it{c});DCA cascade daughters (cm)", 250, 0.0, 25.0, 25, 0.0, 2.5);
     fListHistMultistrangeQA->Add(fHistVarDcaCascDaughtXiPlus);
  }
  if (!fHistVarDcaCascDaughtOmegaMinus) {
     fHistVarDcaCascDaughtOmegaMinus = new TH2F("fHistVarDcaCascDaughtOmegaMinus",";#it{p_{T}} (GeV/#it{c});DCA cascade daughters (cm)", 250, 0.0, 25.0, 25, 0.0, 2.5);
     fListHistMultistrangeQA->Add(fHistVarDcaCascDaughtOmegaMinus);
  }
  if (!fHistVarDcaCascDaughtOmegaPlus) {
     fHistVarDcaCascDaughtOmegaPlus = new TH2F("fHistVarDcaCascDaughtOmegaPlus",";#it{p_{T}} (GeV/#it{c});DCA cascade daughters (cm)", 250, 0.0, 25.0, 25, 0.0, 2.5);
     fListHistMultistrangeQA->Add(fHistVarDcaCascDaughtOmegaPlus);
  }
  // --- DcaBachToPrimVertex (25 bins, [0.0,0.25])
  if (!fHistVarDcaBachToPrimVertexXiMinus) {
     fHistVarDcaBachToPrimVertexXiMinus = new TH2F("fHistVarDcaBachToPrimVertexXiMinus",";#it{p_{T}} (GeV/#it{c});DCA bachelor to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaBachToPrimVertexXiMinus);
  }
  if (!fHistVarDcaBachToPrimVertexXiPlus) {
     fHistVarDcaBachToPrimVertexXiPlus = new TH2F("fHistVarDcaBachToPrimVertexXiPlus",";#it{p_{T}} (GeV/#it{c});DCA bachelor to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaBachToPrimVertexXiPlus);
  }
  if (!fHistVarDcaBachToPrimVertexOmegaMinus) {
     fHistVarDcaBachToPrimVertexOmegaMinus = new TH2F("fHistVarDcaBachToPrimVertexOmegaMinus",";#it{p_{T}} (GeV/#it{c});DCA bachelor to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaBachToPrimVertexOmegaMinus);
  }
  if (!fHistVarDcaBachToPrimVertexOmegaPlus) {
     fHistVarDcaBachToPrimVertexOmegaPlus = new TH2F("fHistVarDcaBachToPrimVertexOmegaPlus",";#it{p_{T}} (GeV/#it{c});DCA bachelor to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaBachToPrimVertexOmegaPlus);
  }
  // --- CascCosineOfPointingAngle (30 bins, [0.97,1.0])
  if (!fHistVarCascCosineOfPointingAngleXiMinus) {
     fHistVarCascCosineOfPointingAngleXiMinus = new TH2F("fHistVarCascCosineOfPointingAngleXiMinus",";#it{p_{T}} (GeV/#it{c});Cascade cosine of PA", 250, 0.0, 25.0, 120, 0.97, 1.0);
     fListHistMultistrangeQA->Add(fHistVarCascCosineOfPointingAngleXiMinus);
  }
  if (!fHistVarCascCosineOfPointingAngleXiPlus) {
     fHistVarCascCosineOfPointingAngleXiPlus = new TH2F("fHistVarCascCosineOfPointingAngleXiPlus",";#it{p_{T}} (GeV/#it{c});Cascade cosine of PA", 250, 0.0, 25.0, 120, 0.97, 1.0);
     fListHistMultistrangeQA->Add(fHistVarCascCosineOfPointingAngleXiPlus);
  }
  if (!fHistVarCascCosineOfPointingAngleOmegaMinus) {
     fHistVarCascCosineOfPointingAngleOmegaMinus = new TH2F("fHistVarCascCosineOfPointingAngleOmegaMinus",";#it{p_{T}} (GeV/#it{c});Cascade cosine of PA", 250, 0.0, 25.0, 120, 0.97, 1.0);
     fListHistMultistrangeQA->Add(fHistVarCascCosineOfPointingAngleOmegaMinus);
  }
  if (!fHistVarCascCosineOfPointingAngleOmegaPlus) {
     fHistVarCascCosineOfPointingAngleOmegaPlus = new TH2F("fHistVarCascCosineOfPointingAngleOmegaPlus",";#it{p_{T}} (GeV/#it{c});Cascade cosine of PA", 250, 0.0, 25.0, 120, 0.97, 1.0);
     fListHistMultistrangeQA->Add(fHistVarCascCosineOfPointingAngleOmegaPlus);
  }
  // --- CascRadius (40 bins, [0.0,4.0)
  if (!fHistVarCascRadiusXiMinus) {
     fHistVarCascRadiusXiMinus = new TH2F("fHistVarCascRadiusXiMinus",";#it{p_{T}} (GeV/#it{c});Cascade fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarCascRadiusXiMinus);
  }
  if (!fHistVarCascRadiusXiPlus) {
     fHistVarCascRadiusXiPlus = new TH2F("fHistVarCascRadiusXiPlus",";#it{p_{T}} (GeV/#it{c});Cascade fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarCascRadiusXiPlus);
  }
  if (!fHistVarCascRadiusOmegaMinus) {
     fHistVarCascRadiusOmegaMinus = new TH2F("fHistVarCascRadiusOmegaMinus",";#it{p_{T}} (GeV/#it{c});Cascade fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarCascRadiusOmegaMinus);
  }
  if (!fHistVarCascRadiusOmegaPlus) {
     fHistVarCascRadiusOmegaPlus = new TH2F("fHistVarCascRadiusOmegaPlus",";#it{p_{T}} (GeV/#it{c});Cascade fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarCascRadiusOmegaPlus);
  }
  // --- InvMassLambdaAsCascDghter (30 bins, [1.1,1.3])
  if (!fHistVarInvMassLambdaAsCascDghterXiMinus) {
     fHistVarInvMassLambdaAsCascDghterXiMinus = new TH2F("fHistVarInvMassLambdaAsCascDghterXiMinus",";#it{p_{T}} (GeV/#it{c});V^{0} invariant mass (GeV/c^{2})", 250, 0.0, 25.0, 30, 1.10, 1.13);
     fListHistMultistrangeQA->Add(fHistVarInvMassLambdaAsCascDghterXiMinus);
  }
  if (!fHistVarInvMassLambdaAsCascDghterXiPlus) {
     fHistVarInvMassLambdaAsCascDghterXiPlus = new TH2F("fHistVarInvMassLambdaAsCascDghterXiPlus",";#it{p_{T}} (GeV/#it{c});V^{0} invariant mass (GeV/c^{2})", 250, 0.0, 25.0, 30, 1.10, 1.13);
     fListHistMultistrangeQA->Add(fHistVarInvMassLambdaAsCascDghterXiPlus);
  }
  if (!fHistVarInvMassLambdaAsCascDghterOmegaMinus) {
     fHistVarInvMassLambdaAsCascDghterOmegaMinus = new TH2F("fHistVarInvMassLambdaAsCascDghterOmegaMinus",";#it{p_{T}} (GeV/#it{c});V^{0} invariant mass (GeV/c^{2})", 250, 0.0, 25.0, 30, 1.10, 1.13);
     fListHistMultistrangeQA->Add(fHistVarInvMassLambdaAsCascDghterOmegaMinus);
  }
  if (!fHistVarInvMassLambdaAsCascDghterOmegaPlus) {
     fHistVarInvMassLambdaAsCascDghterOmegaPlus = new TH2F("fHistVarInvMassLambdaAsCascDghterOmegaPlus",";#it{p_{T}} (GeV/#it{c});V^{0} invariant mass (GeV/c^{2})", 250, 0.0, 25.0, 30, 1.10, 1.13);
     fListHistMultistrangeQA->Add(fHistVarInvMassLambdaAsCascDghterOmegaPlus);
  }
  // --- DcaV0Daughters (20 bins, [0.0,2.0])
  if (!fHistVarDcaV0DaughtersXiMinus) {
     fHistVarDcaV0DaughtersXiMinus = new TH2F("fHistVarDcaV0DaughtersXiMinus",";#it{p_{T}} (GeV/#it{c});DCA V^{0} daughters (cm)", 250, 0.0, 25.0, 20, 0.0, 2.0);
     fListHistMultistrangeQA->Add(fHistVarDcaV0DaughtersXiMinus);
  }
  if (!fHistVarDcaV0DaughtersXiPlus) {
     fHistVarDcaV0DaughtersXiPlus = new TH2F("fHistVarDcaV0DaughtersXiPlus",";#it{p_{T}} (GeV/#it{c});DCA V^{0} daughters (cm)", 250, 0.0, 25.0, 20, 0.0, 2.0);
     fListHistMultistrangeQA->Add(fHistVarDcaV0DaughtersXiPlus);
  }
  if (!fHistVarDcaV0DaughtersOmegaMinus) {
     fHistVarDcaV0DaughtersOmegaMinus = new TH2F("fHistVarDcaV0DaughtersOmegaMinus",";#it{p_{T}} (GeV/#it{c});DCA V^{0} daughters (cm)", 250, 0.0, 25.0, 20, 0.0, 2.0);
     fListHistMultistrangeQA->Add(fHistVarDcaV0DaughtersOmegaMinus);
  }
  if (!fHistVarDcaV0DaughtersOmegaPlus) {
     fHistVarDcaV0DaughtersOmegaPlus = new TH2F("fHistVarDcaV0DaughtersOmegaPlus",";#it{p_{T}} (GeV/#it{c});DCA V^{0} daughters (cm)", 250, 0.0, 25.0, 20, 0.0, 2.0);
     fListHistMultistrangeQA->Add(fHistVarDcaV0DaughtersOmegaPlus);
  }
  // --- V0CosineOfPointingAngleToCascVertex (100 bins, [0.9,1.0])
  if (!fHistVarV0CosineOfPAToCascVertexXiMinus) {
     fHistVarV0CosineOfPAToCascVertexXiMinus = new TH2F("fHistVarV0CosineOfPAToCascVertexXiMinus",";#it{p_{T}} (GeV/#it{c});V^{0} cosine of PA (to cascade vertex)", 250, 0.0, 25.0, 100, 0.9, 1.0);
     fListHistMultistrangeQA->Add(fHistVarV0CosineOfPAToCascVertexXiMinus);
  }
  if (!fHistVarV0CosineOfPAToCascVertexXiPlus) {
     fHistVarV0CosineOfPAToCascVertexXiPlus = new TH2F("fHistVarV0CosineOfPAToCascVertexXiPlus",";#it{p_{T}} (GeV/#it{c});V^{0} cosine of PA (to cascade vertex)", 250, 0.0, 25.0, 100, 0.9, 1.0);
     fListHistMultistrangeQA->Add(fHistVarV0CosineOfPAToCascVertexXiPlus);
  }
  if (!fHistVarV0CosineOfPAToCascVertexOmegaMinus) {
     fHistVarV0CosineOfPAToCascVertexOmegaMinus = new TH2F("fHistVarV0CosineOfPAToCascVertexOmegaMinus",";#it{p_{T}} (GeV/#it{c});V^{0} cosine of PA (to cascade vertex)", 250, 0.0, 25.0, 100, 0.9, 1.0);
     fListHistMultistrangeQA->Add(fHistVarV0CosineOfPAToCascVertexOmegaMinus);
  }
  if (!fHistVarV0CosineOfPAToCascVertexOmegaPlus) {
     fHistVarV0CosineOfPAToCascVertexOmegaPlus = new TH2F("fHistVarV0CosineOfPAToCascVertexOmegaPlus",";#it{p_{T}} (GeV/#it{c});V^{0} cosine of PA (to cascade vertex)", 250, 0.0, 25.0, 100, 0.9, 1.0);
     fListHistMultistrangeQA->Add(fHistVarV0CosineOfPAToCascVertexOmegaPlus);
  }
  // --- V0Radius (40 bins, [0.0,4.0])
  if (!fHistVarV0RadiusXiMinus) {
     fHistVarV0RadiusXiMinus = new TH2F("fHistVarV0RadiusXiMinus",";#it{p_{T}} (GeV/#it{c});V^{0} fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarV0RadiusXiMinus);
  }
  if (!fHistVarV0RadiusXiPlus) {
     fHistVarV0RadiusXiPlus = new TH2F("fHistVarV0RadiusXiPlus",";#it{p_{T}} (GeV/#it{c});V^{0} fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarV0RadiusXiPlus);
  }
  if (!fHistVarV0RadiusOmegaMinus) {
     fHistVarV0RadiusOmegaMinus = new TH2F("fHistVarV0RadiusOmegaMinus",";#it{p_{T}} (GeV/#it{c});V^{0} fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarV0RadiusOmegaMinus);
  }
  if (!fHistVarV0RadiusOmegaPlus) {
     fHistVarV0RadiusOmegaPlus = new TH2F("fHistVarV0RadiusOmegaPlus",";#it{p_{T}} (GeV/#it{c});V^{0} fiducial volume radius (cm)", 250, 0.0, 25.0, 40, 0.0, 4.0);
     fListHistMultistrangeQA->Add(fHistVarV0RadiusOmegaPlus);
  }
  // --- DcaV0ToPrimVertex (40 bins, [0.0,0.4])
  if (!fHistVarDcaV0ToPrimVertexXiMinus) {
     fHistVarDcaV0ToPrimVertexXiMinus = new TH2F("fHistVarDcaV0ToPrimVertexXiMinus",";#it{p_{T}} (GeV/#it{c});V^{0} DCA to PV (cm)", 250, 0.0, 25.0, 40, 0.0, 0.4);
     fListHistMultistrangeQA->Add(fHistVarDcaV0ToPrimVertexXiMinus);
  }
  if (!fHistVarDcaV0ToPrimVertexXiPlus) {
     fHistVarDcaV0ToPrimVertexXiPlus = new TH2F("fHistVarDcaV0ToPrimVertexXiPlus",";#it{p_{T}} (GeV/#it{c});V^{0} DCA to PV (cm)", 250, 0.0, 25.0, 40, 0.0, 0.4);
     fListHistMultistrangeQA->Add(fHistVarDcaV0ToPrimVertexXiPlus);
  }
  if (!fHistVarDcaV0ToPrimVertexOmegaMinus) {
     fHistVarDcaV0ToPrimVertexOmegaMinus = new TH2F("fHistVarDcaV0ToPrimVertexOmegaMinus",";#it{p_{T}} (GeV/#it{c});V^{0} DCA to PV (cm)", 250, 0.0, 25.0, 40, 0.0, 0.4);
     fListHistMultistrangeQA->Add(fHistVarDcaV0ToPrimVertexOmegaMinus);
  }
  if (!fHistVarDcaV0ToPrimVertexOmegaPlus) {
     fHistVarDcaV0ToPrimVertexOmegaPlus = new TH2F("fHistVarDcaV0ToPrimVertexOmegaPlus",";#it{p_{T}} (GeV/#it{c});V^{0} DCA to PV (cm)", 250, 0.0, 25.0, 40, 0.0, 0.4);
     fListHistMultistrangeQA->Add(fHistVarDcaV0ToPrimVertexOmegaPlus);
  }
  // --- DcaPosToPrimVertex (25 bins, [0.0,0.25])
  if (!fHistVarDcaPosToPrimVertexXiMinus) {
     fHistVarDcaPosToPrimVertexXiMinus = new TH2F("fHistVarDcaPosToPrimVertexXiMinus",";#it{p_{T}} (GeV/#it{c});Positive V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaPosToPrimVertexXiMinus);
  }
  if (!fHistVarDcaPosToPrimVertexXiPlus) {
     fHistVarDcaPosToPrimVertexXiPlus = new TH2F("fHistVarDcaPosToPrimVertexXiPlus",";#it{p_{T}} (GeV/#it{c});Positive V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaPosToPrimVertexXiPlus);
  }
  if (!fHistVarDcaPosToPrimVertexOmegaMinus) {
     fHistVarDcaPosToPrimVertexOmegaMinus = new TH2F("fHistVarDcaPosToPrimVertexOmegaMinus",";#it{p_{T}} (GeV/#it{c});Positive V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaPosToPrimVertexOmegaMinus);
  }
  if (!fHistVarDcaPosToPrimVertexOmegaPlus) {
     fHistVarDcaPosToPrimVertexOmegaPlus = new TH2F("fHistVarDcaPosToPrimVertexOmegaPlus",";#it{p_{T}} (GeV/#it{c});Positive V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaPosToPrimVertexOmegaPlus);
  }  
  // --- DcaNegToPrimVertex (25 bins, [0.0,0.25])
  if (!fHistVarDcaNegToPrimVertexXiMinus) {
     fHistVarDcaNegToPrimVertexXiMinus = new TH2F("fHistVarDcaNegToPrimVertexXiMinus",";#it{p_{T}} (GeV/#it{c});Negative V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaNegToPrimVertexXiMinus);
  }
  if (!fHistVarDcaNegToPrimVertexXiPlus) {
     fHistVarDcaNegToPrimVertexXiPlus = new TH2F("fHistVarDcaNegToPrimVertexXiPlus",";#it{p_{T}} (GeV/#it{c});Negative V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaNegToPrimVertexXiPlus);
  }
  if (!fHistVarDcaNegToPrimVertexOmegaMinus) {
     fHistVarDcaNegToPrimVertexOmegaMinus = new TH2F("fHistVarDcaNegToPrimVertexOmegaMinus",";#it{p_{T}} (GeV/#it{c});Negative V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaNegToPrimVertexOmegaMinus);
  }
  if (!fHistVarDcaNegToPrimVertexOmegaPlus) {
     fHistVarDcaNegToPrimVertexOmegaPlus = new TH2F("fHistVarDcaNegToPrimVertexOmegaPlus",";#it{p_{T}} (GeV/#it{c});Negative V^{0} daughter DCA to PV (cm)", 250, 0.0, 25.0, 25, 0.0, 0.25);
     fListHistMultistrangeQA->Add(fHistVarDcaNegToPrimVertexOmegaPlus);
  }  
  // -- Invariant mass distributions (150 bins, [1.25,1.40] for Xi and [1.62,1.74] for Omega)
  if (! fHistMassXiMinus) {
     fHistMassXiMinus = new TH1F("fHistMassXiMinus", "#Xi^{-} candidates;M(#Lambda,#pi^{-}) (GeV/c^{2}); Counts", 150, 1.25, 1.40);
     fListHistMultistrangeQA->Add(fHistMassXiMinus);
  }
  if (! fHistMassXiPlus) {
     fHistMassXiPlus = new TH1F("fHistMassXiPlus", "#Xi^{+} candidates; M(#bar{#Lambda}^{0},#pi^{+}) (GeV/c^{2}); Counts", 150, 1.25, 1.40);
     fListHistMultistrangeQA->Add(fHistMassXiPlus);
  }
  if (! fHistMassOmegaMinus) {
     fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus", "#Omega^{-} candidates; M(#Lambda,K^{-}) (GeV/c^{2}); Counts", 120, 1.62, 1.74);
     fListHistMultistrangeQA->Add(fHistMassOmegaMinus);
  }
  if (! fHistMassOmegaPlus) {
     fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus", "#Omega^{+} candidates;M(#bar{#Lambda}^{0},K^{+}) (GeV/c^{2}); Counts", 120, 1.62, 1.74);
     fListHistMultistrangeQA->Add(fHistMassOmegaPlus);
  }
  // --- CascTransvMom (250, [0.0,25.0])
  if (!fHistVarTransvMomentumXiMinus) {
     fHistVarTransvMomentumXiMinus = new TH1F("fHistVarTransvMomentumXiMinus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistVarTransvMomentumXiMinus);
  }
  if (!fHistVarTransvMomentumXiPlus) {
     fHistVarTransvMomentumXiPlus = new TH1F("fHistVarTransvMomentumXiPlus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistVarTransvMomentumXiPlus);
  }
  if (!fHistVarTransvMomentumOmegaMinus) {
     fHistVarTransvMomentumOmegaMinus = new TH1F("fHistVarTransvMomentumOmegaMinus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistVarTransvMomentumOmegaMinus);
  }
  if (!fHistVarTransvMomentumOmegaPlus) {
     fHistVarTransvMomentumOmegaPlus = new TH1F("fHistVarTransvMomentumOmegaPlus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistVarTransvMomentumOmegaPlus);
  }
  // --- Y (110 bins, [-1.1, 1.1])
  if (!fHistVarRapidityXiMinus) {
     fHistVarRapidityXiMinus = new TH1F("fHistVarRapidityXiMinus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistVarRapidityXiMinus);
  }
  if (!fHistVarRapidityXiPlus) {
     fHistVarRapidityXiPlus = new TH1F("fHistVarRapidityXiPlus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistVarRapidityXiPlus);
  }
  if (!fHistVarRapidityOmegaMinus) {
     fHistVarRapidityOmegaMinus = new TH1F("fHistVarRapidityOmegaMinus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistVarRapidityOmegaMinus);
  }
  if (!fHistVarRapidityOmegaPlus) {
     fHistVarRapidityOmegaPlus = new TH1F("fHistVarRapidityOmegaPlus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistVarRapidityOmegaPlus);
  }
  // --- CascadeProperLength (100 bins, [0.0, 100.])
  if (!fHistVarCascProperLengthXiMinus) {
     fHistVarCascProperLengthXiMinus =  new TH2F("fHistVarCascProperLengthXiMinus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarCascProperLengthXiMinus);
  }
  if (!fHistVarCascProperLengthXiPlus) {
     fHistVarCascProperLengthXiPlus =  new TH2F("fHistVarCascProperLengthXiPlus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarCascProperLengthXiPlus);
  }
  if (!fHistVarCascProperLengthOmegaMinus) {
     fHistVarCascProperLengthOmegaMinus =  new TH2F("fHistVarCascProperLengthOmegaMinus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarCascProperLengthOmegaMinus);
  }
  if (!fHistVarCascProperLengthOmegaPlus) {
     fHistVarCascProperLengthOmegaPlus =  new TH2F("fHistVarCascProperLengthOmegaPlus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarCascProperLengthOmegaPlus);
  }
  // --- V0ProperLength (100 bins, [0.0, 100.])
  if (!fHistVarV0ProperLengthXiMinus) {
     fHistVarV0ProperLengthXiMinus =  new TH2F("fHistVarV0ProperLengthXiMinus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarV0ProperLengthXiMinus);
  }
  if (!fHistVarV0ProperLengthXiPlus) {
     fHistVarV0ProperLengthXiPlus =  new TH2F("fHistVarV0ProperLengthXiPlus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarV0ProperLengthXiPlus);
  }
  if (!fHistVarV0ProperLengthOmegaMinus) {
     fHistVarV0ProperLengthOmegaMinus =  new TH2F("fHistVarV0ProperLengthOmegaMinus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarV0ProperLengthOmegaMinus);
  }
  if (!fHistVarV0ProperLengthOmegaPlus) {
     fHistVarV0ProperLengthOmegaPlus =  new TH2F("fHistVarV0ProperLengthOmegaPlus",";#it{p_{T}} (GeV/#it{c});mL/p (cm)", 250, 0.0, 25.0, 100, 0.0, 100.);
     fListHistMultistrangeQA->Add(fHistVarV0ProperLengthOmegaPlus);
  }  
  // -- Cascade variable distributions for generated particles  
  // --- Total Momentum [250 bins, (0.0,25.0)]
  if (!fHistGenVarTotMomXiMinus) {
     fHistGenVarTotMomXiMinus = new TH1F("fHistGenVarTotMomXiMinus",";#it{p} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTotMomXiMinus);
  }
  if (!fHistGenVarTotMomXiPlus) {
     fHistGenVarTotMomXiPlus = new TH1F("fHistGenVarTotMomXiPlus",";#it{p} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTotMomXiPlus);
  }
  if (!fHistGenVarTotMomOmegaMinus) {
     fHistGenVarTotMomOmegaMinus = new TH1F("fHistGenVarTotMomOmegaMinus",";#it{p} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTotMomOmegaMinus);
  }
  if (!fHistGenVarTotMomOmegaPlus) {
     fHistGenVarTotMomOmegaPlus = new TH1F("fHistGenVarTotMomOmegaPlus",";#it{p} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTotMomOmegaPlus);
  }
  // --- Transverse Momentum [250 bins, (0.0,25.0)]
  if (!fHistGenVarTransvMomXiMinus) {
     fHistGenVarTransvMomXiMinus = new TH1F("fHistGenVarTransvMomXiMinus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTransvMomXiMinus);
  }
  if (!fHistGenVarTransvMomXiPlus) {
     fHistGenVarTransvMomXiPlus = new TH1F("fHistGenVarTransvMomXiPlus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTransvMomXiPlus);
  }
  if (!fHistGenVarTransvMomOmegaMinus) {
     fHistGenVarTransvMomOmegaMinus = new TH1F("fHistGenVarTransvMomOmegaMinus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTransvMomOmegaMinus);
  }
  if (!fHistGenVarTransvMomOmegaPlus) {
     fHistGenVarTransvMomOmegaPlus = new TH1F("fHistGenVarTransvMomOmegaPlus",";#it{p_{T}} (GeV/#it{c});Counts", 250, 0.0, 25.0);
     fListHistMultistrangeQA->Add(fHistGenVarTransvMomOmegaPlus);
  }
  // --- Y [110 bins, (-1.1,1.1)]
  if (!fHistGenVarYXiMinus) {
     fHistGenVarYXiMinus = new TH1F("fHistGenVarYXiMinus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistGenVarYXiMinus);
  }
  if (!fHistGenVarYXiPlus) {
     fHistGenVarYXiPlus = new TH1F("fHistGenVarYXiPlus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistGenVarYXiPlus);
  }
  if (!fHistGenVarYOmegaMinus) {
     fHistGenVarYOmegaMinus = new TH1F("fHistGenVarYOmegaMinus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistGenVarYOmegaMinus);
  }
  if (!fHistGenVarYOmegaPlus) {
     fHistGenVarYOmegaPlus = new TH1F("fHistGenVarYOmegaPlus",";Y;Counts", 110, -1.1, 1.1);
     fListHistMultistrangeQA->Add(fHistGenVarYOmegaPlus);
  }
  // --- Eta [200 bins, (-10.0,10.0)]
  if (!fHistGenVarEtaXiMinus) {
     fHistGenVarEtaXiMinus = new TH1F("fHistGenVarEtaXiMinus",";#eta;Counts", 200, -10.0, 10.0);
     fListHistMultistrangeQA->Add(fHistGenVarEtaXiMinus);
  }
  if (!fHistGenVarEtaXiPlus) {
     fHistGenVarEtaXiPlus = new TH1F("fHistGenVarEtaXiPlus",";#eta;Counts", 200, -10.0, 10.0);
     fListHistMultistrangeQA->Add(fHistGenVarEtaXiPlus);
  }
  if (!fHistGenVarEtaOmegaMinus) {
     fHistGenVarEtaOmegaMinus = new TH1F("fHistGenVarEtaOmegaMinus",";#eta;Counts", 200, -10.0, 10.0);
     fListHistMultistrangeQA->Add(fHistGenVarEtaOmegaMinus);
  }
  if (!fHistGenVarEtaOmegaPlus) {
     fHistGenVarEtaOmegaPlus = new TH1F("fHistGenVarEtaOmegaPlus",";#eta;Counts", 200, -10.0, 10.0);
     fListHistMultistrangeQA->Add(fHistGenVarEtaOmegaPlus);
  }
  // --- Theta [200 bins, (-10.0,190.0)]
  if (!fHistGenVarThetaXiMinus) {
     fHistGenVarThetaXiMinus = new TH1F("fHistGenVarThetaXiMinus",";#theta;Counts", 200, -10.0, 190.0);
     fListHistMultistrangeQA->Add(fHistGenVarThetaXiMinus);
  }
  if (!fHistGenVarThetaXiPlus) {
     fHistGenVarThetaXiPlus = new TH1F("fHistGenVarThetaXiPlus",";#theta;Counts", 200, -10.0, 190.0);
     fListHistMultistrangeQA->Add(fHistGenVarThetaXiPlus);
  }
  if (!fHistGenVarThetaOmegaMinus) {
     fHistGenVarThetaOmegaMinus = new TH1F("fHistGenVarThetaOmegaMinus",";#theta;Counts", 200, -10.0, 190.0);
     fListHistMultistrangeQA->Add(fHistGenVarThetaOmegaMinus);
  }
  if (!fHistGenVarThetaOmegaPlus) {
     fHistGenVarThetaOmegaPlus = new TH1F("fHistGenVarThetaOmegaPlus",";#theta;Counts", 200, -10.0, 190.0);
     fListHistMultistrangeQA->Add(fHistGenVarThetaOmegaPlus);
  }
  // --- Phi [180 bins, (0.0,360.0)]
  if (!fHistGenVarPhiXiMinus) {
     fHistGenVarPhiXiMinus = new TH1F("fHistGenVarPhiXiMinus",";#phi;Counts", 180, 0.0, 360.0);
     fListHistMultistrangeQA->Add(fHistGenVarPhiXiMinus);
  }
  if (!fHistGenVarPhiXiPlus) {
     fHistGenVarPhiXiPlus = new TH1F("fHistGenVarPhiXiPlus",";#phi;Counts", 180, 0.0, 360.0);
     fListHistMultistrangeQA->Add(fHistGenVarPhiXiPlus);
  }
  if (!fHistGenVarPhiOmegaMinus) {
     fHistGenVarPhiOmegaMinus = new TH1F("fHistGenVarPhiOmegaMinus",";#phi;Counts", 180, 0.0, 360.0);
     fListHistMultistrangeQA->Add(fHistGenVarPhiOmegaMinus);
  }
  if (!fHistGenVarPhiOmegaPlus) {
     fHistGenVarPhiOmegaPlus = new TH1F("fHistGenVarPhiOmegaPlus",";#phi;Counts", 180, 0.0, 360.0);
     fListHistMultistrangeQA->Add(fHistGenVarPhiOmegaPlus);
  }


  PostData(1, fListHistMultistrangeQA);


}// end UserCreateOutputObjects


//------------------------------------------------------
void AliAnalysisTaskQAMultistrange::UserExec(Option_t *) 

{


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Main loop (called for each event)
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   //________________________
   // Define ESD/AOD handlers
   AliESDEvent  *lESDevent = 0x0;
   AliAODEvent  *lAODevent = 0x0;
   AliMCEvent   *lMCevent  = 0x0; //for MC
   AliStack     *lMCstack  = 0x0; //for MC
   TClonesArray *arrayMC   = 0;   //for MC  

   //______________________
   // Check the PIDresponse
   if(!fPIDResponse) { AliError("Cannot get pid response");  return; }
   
   //_____________________
   // Check the InputEvent 
   if (fAnalysisType == "ESD") {
       lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
       if (!lESDevent) { AliWarning("ERROR: lESDevent not available \n");  return; }
       if (fisMC) {
           lMCevent = MCEvent();
           if (!lMCevent) { AliWarning("ERROR: Could not retrieve MC event \n");  return; }
           lMCstack = lMCevent->Stack();
           if (!lMCstack) { AliWarning("ERROR: Could not retrieve MC stack \n");  return; }
       }
   } else if (fAnalysisType == "AOD") {
       lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
       if (!lAODevent) { AliWarning("ERROR: lAODevent not available \n");  return; }
       if (fisMC) {
           arrayMC = (TClonesArray*) lAODevent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
           if (!arrayMC) { AliWarning("ERROR: MC particles branch not found!\n"); return; }
       }
   } else {
     AliWarning("Analysis type (ESD or AOD) not specified \n");
     return;
   }

   //_________________________________________________
   // - Fill the event plot before any event selection 
   fHistEventSel->Fill(0.5);

   //______________________________________________________________________________________________
   // - Perform the event selection (via AliPPVsMultUtils) and acquire the multiplicity information
   Int_t lEvSelCode = 100;
   AliMultSelection *MultSelection = 0x0;
   if      (fAnalysisType == "ESD")  MultSelection = (AliMultSelection*) lESDevent->FindListObject("MultSelection");
   else if (fAnalysisType == "AOD")  MultSelection = (AliMultSelection*) lAODevent->FindListObject("MultSelection");
   if (!MultSelection) {
          AliWarning("AliMultSelection object not found!");  //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
          PostData(1, fListHistMultistrangeQA);
          return;
   } else {
          AliWarning("AliMultSelection object found!");
          lEvSelCode = MultSelection->GetEvSelCode();  //Event Selection Code
   }
   // - Remove events
   if (lEvSelCode != 0) {
          AliWarning(Form("lEvSelCode value = %i. Run Not good! REMOVE",lEvSelCode));
          PostData(1, fListHistMultistrangeQA);
          return;
   }

   //_________________________________________________
   // - Fill the event plot after all event selections 
   fHistEventSel->Fill(1.5); 

   //_________________
   // - Magnetic field
   Double_t lMagneticField = -10.;
   if      (fAnalysisType == "ESD") lMagneticField = lESDevent->GetMagneticField();
   else if (fAnalysisType == "AOD") lMagneticField = lAODevent->GetMagneticField();

   //__________________________________________________________
   // - Get Vertex and fill the plot before any event selection
   Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
   if (fAnalysisType == "ESD") { 
       const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();
       if (!lPrimaryBestESDVtx) { 
             AliWarning("No prim. vertex in ESD... return!");
             PostData(1, fListHistMultistrangeQA);
             return;
       }
       lPrimaryBestESDVtx->GetXYZ(lBestPrimaryVtxPos);
   } else if (fAnalysisType == "AOD") {
       const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
       if (!lPrimaryBestAODVtx) {
             AliWarning("No prim. vertex in AOD... return!");
             PostData(1, fListHistMultistrangeQA);
             return;
       }
       lPrimaryBestAODVtx->GetXYZ(lBestPrimaryVtxPos);
   }


  ////////////////////////////               
  // MC GENERATED CASCADE PART
  ////////////////////////////

  //%%%%%%%%%%%%%%%%%
  // Gen cascade loop
  if (fisMC) {

      Int_t lNbMCPrimary = 0;
      if      (fAnalysisType == "ESD") lNbMCPrimary = lMCstack->GetNtrack();   
      else if (fAnalysisType == "AOD") lNbMCPrimary = arrayMC->GetEntries();
      Int_t ngenximinus    = 0;
      Int_t ngenxiplus     = 0;
      Int_t ngenomegaminus = 0;
      Int_t ngenomegaplus  = 0;


      for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) {

           Float_t partP      = 0.;
           Float_t partPt     = 0.;
           Float_t partEta    = 0.;
           Float_t partTheta  = 0.;
           Float_t partPhi    = 0.;
           Float_t partRap    = 0.;
           Float_t partEnergy = 0.; //for Rapidity
           Float_t partPz     = 0.; //for Rapidity
           Int_t    PDGcode    = 0;

           if ( fAnalysisType == "ESD" ) {
               TParticle* lCurrentParticlePrimary = 0x0;
               lCurrentParticlePrimary = lMCstack->Particle( iCurrentLabelStack );
               if (!lCurrentParticlePrimary) {
                   AliWarning("Cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...");
                   continue;
               }
               if (!lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) continue;
               TParticle* cascMC = 0x0;
               cascMC = lCurrentParticlePrimary;
               if (!cascMC) {
                   AliWarning("MC TParticle pointer to Cascade = 0x0 ! Skip ...");
                   continue;
               }
               partP      = cascMC->P();
               partPt     = cascMC->Pt();
               partEta    = cascMC->Eta();
               partTheta  = cascMC->Theta()*180.0/TMath::Pi();
               partPhi    = cascMC->Phi()*180.0/TMath::Pi();
               partEnergy = cascMC->Energy(); 
               partPz     = cascMC->Pz();
               PDGcode    = lCurrentParticlePrimary->GetPdgCode();
           } else if ( fAnalysisType == "AOD" ) {
               AliAODMCParticle *lCurrentParticleaod = (AliAODMCParticle*) arrayMC->At(iCurrentLabelStack);
               if (!lCurrentParticleaod) {
                   AliWarning("Cascade loop - MC TParticle pointer to current stack particle = 0x0 ! Skip ...");
                   continue;
               }
               if (!lCurrentParticleaod->IsPhysicalPrimary()) continue;
               partP      = lCurrentParticleaod->P();
               partPt     = lCurrentParticleaod->Pt();
               partEta    = lCurrentParticleaod->Eta(); 
               partTheta  = lCurrentParticleaod->Theta()*180.0/TMath::Pi();
               partPhi    = lCurrentParticleaod->Phi()*180.0/TMath::Pi();
               partEnergy = lCurrentParticleaod->E();
               partPz     = lCurrentParticleaod->Pz();
               PDGcode    = lCurrentParticleaod->GetPdgCode();
           }
           partRap = 0.5*TMath::Log((partEnergy + partPz) / (partEnergy - partPz + 1.e-13));

           if (PDGcode == 3312) {
               fHistGenVarTotMomXiMinus->Fill(partP);
               fHistGenVarTransvMomXiMinus->Fill(partPt);
               fHistGenVarYXiMinus->Fill(partRap);
               fHistGenVarEtaXiMinus->Fill(partEta);
               fHistGenVarThetaXiMinus->Fill(partTheta);
               fHistGenVarPhiXiMinus->Fill(partPhi);
               ngenximinus++;
           } else if (PDGcode == -3312) {
               fHistGenVarTotMomXiPlus->Fill(partP);
               fHistGenVarTransvMomXiPlus->Fill(partPt);
               fHistGenVarYXiPlus->Fill(partRap);
               fHistGenVarEtaXiPlus->Fill(partEta);
               fHistGenVarThetaXiPlus->Fill(partTheta);
               fHistGenVarPhiXiPlus->Fill(partPhi);
               ngenxiplus++;
           } else if (PDGcode == 3334)  {
               fHistGenVarTotMomOmegaMinus->Fill(partP);
               fHistGenVarTransvMomOmegaMinus->Fill(partPt);
               fHistGenVarYOmegaMinus->Fill(partRap);
               fHistGenVarEtaOmegaMinus->Fill(partEta);
               fHistGenVarThetaOmegaMinus->Fill(partTheta);
               fHistGenVarPhiOmegaMinus->Fill(partPhi);
               ngenomegaminus++;
           } else if (PDGcode == -3334) {
               fHistGenVarTotMomOmegaPlus->Fill(partP);
               fHistGenVarTransvMomOmegaPlus->Fill(partPt);
               fHistGenVarYOmegaPlus->Fill(partRap);
               fHistGenVarEtaOmegaPlus->Fill(partEta);
               fHistGenVarThetaOmegaPlus->Fill(partTheta);
               fHistGenVarPhiOmegaPlus->Fill(partPhi);
               ngenomegaplus++;
           }

      }
      fHistCascadeMultiplicityXiMinus->Fill(ngenximinus);
      fHistCascadeMultiplicityXiPlus->Fill(ngenxiplus);
      fHistCascadeMultiplicityOmegaMinus->Fill(ngenomegaminus);
      fHistCascadeMultiplicityOmegaPlus->Fill(ngenomegaplus);
  }        

  //////////////////////////////
  // CASCADE RECONSTRUCTION PART
  //////////////////////////////

  //%%%%%%%%%%%%%
  // Cascade loop
  Int_t ncascades = -22;
  if      (fAnalysisType == "ESD") ncascades = lESDevent->GetNumberOfCascades();
  else if (fAnalysisType == "AOD") ncascades = lAODevent->GetNumberOfCascades();

  for (Int_t iXi = 0; iXi < ncascades; iXi++) {// This is the begining of the Cascade loop (ESD or AOD)

    // -------------------------------------
    // - Initialisation of the local variables that will be needed for ESD/AOD
    // -- Container variables (1st round)
    Float_t lDcaXiDaughters              = -1. ;                   //[Container]
    Float_t lXiCosineOfPointingAngle     = -1. ;                   //[Container]
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };             //Useful to define other variables: radius fid. vol., ctau, etc. for cascade
    Float_t lXiRadius                    = -1000. ;                //[Container]
    UShort_t lPosTPCClusters              = -1;                     //To check the quality of the tracks. For ESD only ...
    UShort_t lNegTPCClusters              = -1;                     //To check the quality of the tracks. For ESD only ...
    UShort_t lBachTPCClusters             = -1;                     //To check the quality of the tracks. For ESD only ...
    Float_t lInvMassLambdaAsCascDghter   = 0.;                     //[Container]
    Float_t lDcaV0DaughtersXi            = -1.;                    //[Container]
    Float_t lDcaBachToPrimVertexXi       = -1.;                    //[Container]
    Float_t lDcaV0ToPrimVertexXi         = -1.;                    //[Container]
    Float_t lDcaPosToPrimVertexXi        = -1.;                    //[Container]
    Float_t lDcaNegToPrimVertexXi        = -1.;                    //[Container]
    Float_t lV0CosineOfPointingAngle     = -1.;                    //[Container]
    Float_t lV0toXiCosineOfPointingAngle = -1.;                    //[Container] 
    Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. };             //Useful to define other variables: radius fid. vol., ctau, etc. for VO 
    Float_t lV0RadiusXi                  = -1000.0;                //[Container]
    Double_t lV0quality                   = 0.;                     //  ??
    Float_t lInvMassXiMinus              = 0.;                     //[Container]
    Float_t lInvMassXiPlus               = 0.;                     //[Container]
    Float_t lInvMassOmegaMinus           = 0.;                     //[Container]
    Float_t lInvMassOmegaPlus            = 0.;                     //[Container]
    // -- PID treatment
    Bool_t   lIsBachelorKaonForTPC = kFALSE; 
    Bool_t   lIsBachelorPionForTPC = kFALSE; 
    Bool_t   lIsNegPionForTPC      = kFALSE; 
    Bool_t   lIsPosPionForTPC      = kFALSE; 
    Bool_t   lIsNegProtonForTPC    = kFALSE; 
    Bool_t   lIsPosProtonForTPC    = kFALSE; 
    // -- MC Association
    Int_t    lblPosV0Dghter         = 0; 
    Int_t    lblNegV0Dghter         = 0;
    Int_t    lblMotherPosV0Dghter   = 0;
    Int_t    lblMotherNegV0Dghter   = 0; 
    Int_t    lblBach                = 0;
    Int_t    lblGdMotherPosV0Dghter = 0;
    Int_t    lblGdMotherNegV0Dghter = 0; 
    Int_t    lblMotherBach          = 0;
    Bool_t   lAssoXiMinus    = kFALSE;
    Bool_t   lAssoXiPlus     = kFALSE;
    Bool_t   lAssoOmegaMinus = kFALSE;
    Bool_t   lAssoOmegaPlus  = kFALSE;
    TParticle *mcPosV0Dghter         = 0x0;
    TParticle *mcNegV0Dghter         = 0x0; 
    TParticle *mcMotherPosV0Dghter   = 0x0; 
    TParticle *mcMotherNegV0Dghter   = 0x0; 
    TParticle *mcBach                = 0x0;
    TParticle *mcGdMotherPosV0Dghter = 0x0; 
    TParticle *mcGdMotherNegV0Dghter = 0x0; 
    TParticle *mcMotherBach          = 0x0;
    AliAODMCParticle *mcPosV0Dghteraod         = 0x0;
    AliAODMCParticle *mcNegV0Dghteraod         = 0x0;
    AliAODMCParticle *mcMotherPosV0Dghteraod   = 0x0;
    AliAODMCParticle *mcMotherNegV0Dghteraod   = 0x0;
    AliAODMCParticle *mcBachaod                = 0x0;
    AliAODMCParticle *mcGdMotherPosV0Dghteraod = 0x0;
    AliAODMCParticle *mcGdMotherNegV0Dghteraod = 0x0;
    AliAODMCParticle *mcMotherBachaod          = 0x0;
    // -- More container variables and quality checks
    Double_t lXiMomX          = 0.;                               //Useful to define other variables: lXiTransvMom, lXiTotMom
    Double_t lXiMomY          = 0.;                               //Useful to define other variables: lXiTransvMom, lXiTotMom
    Double_t lXiMomZ          = 0.;                               //Useful to define other variables: lXiTransvMom, lXiTotMom
    Float_t lXiTransvMom      = 0.;                               //
    Float_t lXiTotMom         = 0.;                               //Useful to define other variables: cTau
    Double_t lV0PMomX         = 0.;                               //Useful to define other variables: lV0TotMom, lpTrackTransvMom
    Double_t lV0PMomY         = 0.;                               //Useful to define other variables: lV0TotMom, lpTrackTransvMom
    Double_t lV0PMomZ         = 0.;                               //Useful to define other variables: lV0TotMom, lpTrackTransvMom
    Double_t lV0NMomX         = 0.;                               //Useful to define other variables: lV0TotMom, lnTrackTransvMom
    Double_t lV0NMomY         = 0.;                               //Useful to define other variables: lV0TotMom, lnTrackTransvMom
    Double_t lV0NMomZ         = 0.;                               //Useful to define other variables: lV0TotMom, lnTrackTransvMom
    Float_t lV0TotMom         = 0.;                               //Useful to define other variables: lctauV0
    Double_t lBachMomX        = 0.;                               //Useful to define other variables: lBachTransvMom
    Double_t lBachMomY        = 0.;                               //Useful to define other variables: lBachTransvMom
    Double_t lBachMomZ        = 0.;                               //Useful to define other variables: lBachTransvMom
    Float_t lBachTransvMom    = 0.;                               //Selection on the min bachelor pT
    Float_t lpTrackTransvMom  = 0.;                               //Selection on the min bachelor pT
    Float_t lnTrackTransvMom  = 0.;                               //Selection on the min bachelor pT
    Short_t  lChargeXi        = -2;                               //Useful to select the particles based on the charge
    Float_t lRapXi            = -20.0;                            //
    Float_t lRapOmega         = -20.0;                            //
    Float_t  etaBach          = 0.;                               //Selection on the eta range
    Float_t  etaPos           = 0.;                               //Selection on the eta range
    Float_t  etaNeg           = 0.;                               //Selection on the eta range
    Float_t cascadeMass       = 0.;
    // --  variables for the cascade cut optimisation: ESD and AOD 
    if (fAnalysisType == "ESD") { 
  
          // -------------------------------------
          // - Load the cascades from the handler 
          AliESDcascade *xi = lESDevent->GetCascade(iXi);
          if (!xi) { AliWarning("ERROR: Cascade not found!"); continue; }
        
          // ---------------------------------------------------------------------------
          // - Assigning the necessary variables for specific AliESDcascade data members 	
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay
          lDcaXiDaughters	   = xi->GetDcaXiDaughters();
          lXiCosineOfPointingAngle = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
		                       // Take care : the best available vertex should be used (like in AliCascadeVertexer)
          xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
          lXiRadius        	= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
		
          // -------------------------------------------------------------------------------------------------------------------------------
          // - Around the tracks : Bach + V0 (ESD). Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
          UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
          UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
          UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                                    // Care track label can be negative in MC production (linked with the track quality)
                                    // However = normally, not the case for track index ...
          // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
          if (lBachIdx == lIdxNegXi) { AliWarning("ERROR: this track has been already used!"); continue; }   
          if (lBachIdx == lIdxPosXi) { AliWarning("ERROR: this track has been already used!"); continue; }
          // - Get the track for the daughters
          AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
          AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
          AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
          if (!pTrackXi || !nTrackXi || !bachTrackXi )  { AliWarning("ERROR: one of the daughter track do not exist!"); continue; }
          // - Get the TPCnumber of cluster for the daughters
          lPosTPCClusters   = pTrackXi->GetTPCNcls();
          lNegTPCClusters   = nTrackXi->GetTPCNcls();
          lBachTPCClusters  = bachTrackXi->GetTPCNcls();
      
          // ------------------------------------
          // - Rejection of a poor quality tracks
          if (fkQualityCutTPCrefit) {
                // 1 - Poor quality related to TPCrefit
                ULong_t pStatus    = pTrackXi->GetStatus();
                ULong_t nStatus    = nTrackXi->GetStatus();
                ULong_t bachStatus = bachTrackXi->GetStatus();
                if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning(" V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning(" V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning(" Bach.   track has no TPCrefit ... continue!"); continue; }
          }
          if (fkQualityCutnTPCcls) {
                // 2 - Poor quality related to TPC clusters
                if (lPosTPCClusters  < fMinnTPCcls) { AliWarning(" V0 Pos. track has less than minn TPC clusters ... continue!"); continue; }
                if (lNegTPCClusters  < fMinnTPCcls) { AliWarning(" V0 Neg. track has less than minn TPC clusters ... continue!"); continue; }
                if (lBachTPCClusters < fMinnTPCcls) { AliWarning(" Bach.   track has less than minn TPC clusters ... continue!"); continue; }
          }

          // ------------------------------
          etaPos  = pTrackXi->Eta();             
          etaNeg  = nTrackXi->Eta();
          etaBach = bachTrackXi->Eta();
          lInvMassLambdaAsCascDghter = xi->GetEffMass();
          lDcaV0DaughtersXi          = xi->GetDcaV0Daughters(); 
          lV0CosineOfPointingAngle   = xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
          lDcaV0ToPrimVertexXi       = xi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );	
          lDcaBachToPrimVertexXi     = TMath::Abs( bachTrackXi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1],	lMagneticField) ); 
          xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
          lV0RadiusXi		     = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );
          lDcaPosToPrimVertexXi      = TMath::Abs( pTrackXi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); 
          lDcaNegToPrimVertexXi      = TMath::Abs( nTrackXi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); 
	
           //----------------------------------------------------------------------------------------------------       
           // - Around effective masses. Change mass hypotheses to cover all the possibilities:  Xi-/+, Omega -/+
           if (bachTrackXi->Charge() < 0) {
                // Calculate the effective mass of the Xi- candidate. pdg code 3312 = Xi-
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3312); 	
                lInvMassXiMinus = xi->GetEffMassXi();
                // Calculate the effective mass of the Xi- candidate. pdg code 3334 = Omega-
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3334); 	
                lInvMassOmegaMinus = xi->GetEffMassXi();
		// Back to default hyp.			
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3312);
           }// end if negative bachelor
           if ( bachTrackXi->Charge() >  0 ) {
                // Calculate the effective mass of the Xi+ candidate. pdg code -3312 = Xi+
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3312); 	
                lInvMassXiPlus = xi->GetEffMassXi();
                // Calculate the effective mass of the Xi+ candidate. pdg code -3334  = Omega+
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3334); 	
                lInvMassOmegaPlus = xi->GetEffMassXi();
		// Back to "default" hyp.
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3312); 
           }// end if positive bachelor

           // ----------------------------------------------	
	   // - TPC PID : 3-sigma bands on Bethe-Bloch curve
           // Bachelor
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
           // Negative V0 daughter
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
           // Positive V0 daughter
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;

           // -----------------------------------------
           // - MC Association in case of MC production 
           if (fisMC) {
             lblPosV0Dghter = 0;  lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );
             lblNegV0Dghter = 0;  lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
             mcPosV0Dghter = 0x0; mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
             mcNegV0Dghter = 0x0; mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );               
             lblMotherPosV0Dghter = 0.;  lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother();
             lblMotherNegV0Dghter = 0.;  lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
             if (lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // must have same mother
             if (lblMotherPosV0Dghter < 0) continue;                     // this particle is primary, no mother   
             if (lblMotherNegV0Dghter < 0) continue;                     // this particle is primary, no mother
             mcMotherPosV0Dghter = 0x0; mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
             mcMotherNegV0Dghter = 0x0; mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );  
             lblBach = 0; lblBach = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
             mcBach = 0x0; mcBach = lMCstack->Particle( lblBach );
             lblGdMotherPosV0Dghter = 0; lblGdMotherPosV0Dghter = mcMotherPosV0Dghter->GetFirstMother() ;
             lblGdMotherNegV0Dghter = 0; lblGdMotherNegV0Dghter = mcMotherNegV0Dghter->GetFirstMother() ;
             lblMotherBach = 0; lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother() );
             if(lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter) continue; // must have same grand-mother
             if(lblGdMotherPosV0Dghter < 0) continue;                       // primary lambda ...   
             if(lblGdMotherNegV0Dghter < 0) continue;                       // primary lambda ...                            
             if(lblMotherBach != lblGdMotherPosV0Dghter) continue;          // must have same mother bach and V0 daughters
             mcGdMotherPosV0Dghter = 0x0; mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
             mcGdMotherNegV0Dghter = 0x0; mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );
             mcMotherBach          = 0x0; mcMotherBach          = lMCstack->Particle( lblMotherBach );
             if (!(lMCstack->IsPhysicalPrimary(lblMotherBach))) continue;
             if      (mcMotherBach->GetPdgCode() == 3312  && mcGdMotherPosV0Dghter->GetPdgCode() == 3312  && mcGdMotherNegV0Dghter->GetPdgCode() == 3312)  {lAssoXiMinus    = kTRUE; cascadeMass = 1.321;}
             else if (mcMotherBach->GetPdgCode() == -3312 && mcGdMotherPosV0Dghter->GetPdgCode() == -3312 && mcGdMotherNegV0Dghter->GetPdgCode() == -3312) {lAssoXiPlus     = kTRUE; cascadeMass = 1.321;}
             else if (mcMotherBach->GetPdgCode() == 3334  && mcGdMotherPosV0Dghter->GetPdgCode() == 3334  && mcGdMotherNegV0Dghter->GetPdgCode() == 3334)  {lAssoOmegaMinus = kTRUE; cascadeMass = 1.672;}
             else if (mcMotherBach->GetPdgCode() == -3334 && mcGdMotherPosV0Dghter->GetPdgCode() == -3334 && mcGdMotherNegV0Dghter->GetPdgCode() == -3334) {lAssoOmegaPlus  = kTRUE; cascadeMass = 1.672;}
           }
        
           // ------------------------------
           // - Miscellaneous pieces of info that may help regarding data quality assessment.
           xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
           lXiTransvMom = TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY );
           lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ );
           xi->GetNPxPyPz(lV0NMomX,lV0NMomY,lV0NMomZ);
           xi->GetPPxPyPz(lV0PMomX,lV0PMomY,lV0PMomZ);
           lV0TotMom = TMath::Sqrt(TMath::Power(lV0NMomX+lV0PMomX,2)+TMath::Power(lV0NMomY+lV0PMomY,2)+TMath::Power(lV0NMomZ+lV0PMomZ,2));	
           xi->GetBPxPyPz( lBachMomX, lBachMomY, lBachMomZ );
           lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY );
           lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX + lV0NMomY*lV0NMomY );
           lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX + lV0PMomY*lV0PMomY );
           lChargeXi = xi->Charge();
           lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
           lRapXi    = xi->RapXi();
           lRapOmega = xi->RapOmega();
	
    } else if (fAnalysisType == "AOD") {

           // -------------------------------------
           // - Load the cascades from the handler	
           const AliAODcascade *xi = lAODevent->GetCascade(iXi);
           if (!xi) continue;
		
	   //----------------------------------------------------------------------------        
           // - Assigning the necessary variables for specific AliESDcascade data members  
           lDcaXiDaughters		= xi->DcaXiDaughters();
           lXiCosineOfPointingAngle 	= xi->CosPointingAngleXi( lBestPrimaryVtxPos[0], 
							  lBestPrimaryVtxPos[1], 
							  lBestPrimaryVtxPos[2] );
           lPosXi[0] = xi->DecayVertexXiX();
           lPosXi[1] = xi->DecayVertexXiY();
           lPosXi[2] = xi->DecayVertexXiZ();
           lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );	

           //-------------------------------------------------------------------------------------------------------------------------------
           // - Around the tracks: Bach + V0 (ESD). Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
           AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
           AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
           AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );
           if (!pTrackXi || !nTrackXi || !bachTrackXi ) continue;
           UInt_t lIdxPosXi  = (UInt_t) TMath::Abs( pTrackXi->GetID() );  
           UInt_t lIdxNegXi  = (UInt_t) TMath::Abs( nTrackXi->GetID() );
           UInt_t lBachIdx   = (UInt_t) TMath::Abs( bachTrackXi->GetID() );
                                // Care track label can be negative in MC production (linked with the track quality)
                                // However = normally, not the case for track index ...

           // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
           if (lBachIdx == lIdxNegXi) continue; 
           if (lBachIdx == lIdxPosXi) continue;
           // - Get the TPCnumber of cluster for the daughters
           lPosTPCClusters   = pTrackXi->GetTPCNcls(); 
           lNegTPCClusters   = nTrackXi->GetTPCNcls();
           lBachTPCClusters  = bachTrackXi->GetTPCNcls();

           // ------------------------------------
           // - Rejection of a poor quality tracks
           if (fkQualityCutTPCrefit) {
                // - Poor quality related to TPCrefit
                if (!(pTrackXi->IsOn(AliAODTrack::kTPCrefit)))    { AliWarning(" V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if (!(nTrackXi->IsOn(AliAODTrack::kTPCrefit)))    { AliWarning(" V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if (!(bachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning(" Bach.   track has no TPCrefit ... continue!"); continue; }
           }
           if (fkQualityCutnTPCcls) {
                // - Poor quality related to TPC clusters
                if (lPosTPCClusters  < fMinnTPCcls) continue; 
                if (lNegTPCClusters  < fMinnTPCcls) continue; 
                if (lBachTPCClusters < fMinnTPCcls) continue; 
           }
  
           // ------------------------------------------------------------------------------------------------------------------------------
           // - Around the tracks: Bach + V0 (AOD). Necessary variables for AODcascade data members coming from the AODv0 part (inheritance)
           etaPos  = pTrackXi->Eta();
           etaNeg  = nTrackXi->Eta();
           etaBach = bachTrackXi->Eta();
           lChargeXi 			= xi->ChargeXi();
           if ( lChargeXi < 0) 	lInvMassLambdaAsCascDghter	= xi->MassLambda();
           else                 lInvMassLambdaAsCascDghter	= xi->MassAntiLambda();
           lDcaV0DaughtersXi 		= xi->DcaV0Daughters(); 
           lDcaV0ToPrimVertexXi 		= xi->DcaV0ToPrimVertex();
           lDcaBachToPrimVertexXi 		= xi->DcaBachToPrimVertex(); 
           lPosV0Xi[0] = xi->DecayVertexV0X();
           lPosV0Xi[1] = xi->DecayVertexV0Y();
           lPosV0Xi[2] = xi->DecayVertexV0Z(); 
           lV0RadiusXi	= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
           lV0CosineOfPointingAngle     = xi->CosPointingAngle( lBestPrimaryVtxPos ); 
           lDcaPosToPrimVertexXi	= xi->DcaPosToPrimVertex(); 
           lDcaNegToPrimVertexXi	= xi->DcaNegToPrimVertex(); 

           // ---------------------------------------------------------------------------------------------------       
           // - Around effective masses. Change mass hypotheses to cover all the possibilities:  Xi-/+, Omega -/+
           if ( lChargeXi < 0 )	lInvMassXiMinus	   = xi->MassXi();
           if ( lChargeXi > 0 )	lInvMassXiPlus 	   = xi->MassXi();
           if ( lChargeXi < 0 )	lInvMassOmegaMinus = xi->MassOmega();
           if ( lChargeXi > 0 )	lInvMassOmegaPlus  = xi->MassOmega();

           // ----------------------------------------------
           // - TPC PID : 3-sigma bands on Bethe-Bloch curve
           // Bachelor
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
           // Negative V0 daughter
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
           // Positive V0 daughter
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;

           // -----------------------------------------
           // - MC Association in case of MC production 
           if (fisMC) {
             lblPosV0Dghter = 0; lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );
             lblNegV0Dghter = 0; lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
             mcPosV0Dghteraod = 0x0; mcPosV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblPosV0Dghter );
             mcNegV0Dghteraod = 0x0; mcNegV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblNegV0Dghter );
             lblMotherPosV0Dghter = 0;  lblMotherPosV0Dghter = mcPosV0Dghteraod->GetMother();
             lblMotherNegV0Dghter = 0;  lblMotherNegV0Dghter = mcNegV0Dghteraod->GetMother();
             if (lblMotherPosV0Dghter != lblMotherNegV0Dghter) continue; // must have same mother
             if (lblMotherPosV0Dghter < 0 ) continue;                    // this particle is primary, no mother
             if (lblMotherNegV0Dghter < 0 ) continue;                    // this particle is primary, no mother
             mcMotherPosV0Dghteraod = 0x0; mcMotherPosV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblMotherPosV0Dghter );
             mcMotherNegV0Dghteraod = 0x0; mcMotherNegV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblMotherNegV0Dghter );
             lblBach = 0;  lblBach = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
             mcBachaod = 0x0; mcBachaod = (AliAODMCParticle*) arrayMC->At( lblBach );
             lblGdMotherPosV0Dghter = 0; lblGdMotherPosV0Dghter = mcMotherPosV0Dghteraod->GetMother() ;
             lblGdMotherNegV0Dghter = 0; lblGdMotherNegV0Dghter = mcMotherNegV0Dghteraod->GetMother() ;
             lblMotherBach = 0; lblMotherBach = (Int_t) TMath::Abs( mcBachaod->GetMother() );
             if (lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter ) continue;
             if (lblGdMotherPosV0Dghter < 0 ) continue;                    // primary lambda ...
             if (lblGdMotherNegV0Dghter < 0 ) continue;                    // primary lambda ...
             if (lblMotherBach != lblGdMotherPosV0Dghter ) continue;       //same mother for bach and V0 daughters
             mcGdMotherPosV0Dghteraod = 0x0; mcGdMotherPosV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblGdMotherPosV0Dghter );
             mcGdMotherNegV0Dghteraod = 0x0; mcGdMotherNegV0Dghteraod = (AliAODMCParticle*) arrayMC->At( lblGdMotherNegV0Dghter );
             mcMotherBachaod          = 0x0; mcMotherBachaod          = (AliAODMCParticle*) arrayMC->At( lblMotherBach );
             // - Check if cascade is primary
             if (!(mcMotherBachaod->IsPhysicalPrimary())) continue;
             // - Manage boolean for association
             if      (mcMotherBachaod->GetPdgCode() == 3312  && mcGdMotherPosV0Dghteraod->GetPdgCode() == 3312  && mcGdMotherNegV0Dghteraod->GetPdgCode() == 3312 ) {lAssoXiMinus = kTRUE;    cascadeMass = 1.321;}
             else if (mcMotherBachaod->GetPdgCode() == -3312 && mcGdMotherPosV0Dghteraod->GetPdgCode() == -3312 && mcGdMotherNegV0Dghteraod->GetPdgCode() == -3312) {lAssoXiPlus = kTRUE;     cascadeMass = 1.321;}
             else if (mcMotherBachaod->GetPdgCode() == 3334  && mcGdMotherPosV0Dghteraod->GetPdgCode() == 3334  && mcGdMotherNegV0Dghteraod->GetPdgCode() == 3334 ) {lAssoOmegaMinus = kTRUE; cascadeMass = 1.672;}
             else if (mcMotherBachaod->GetPdgCode() == -3334 && mcGdMotherPosV0Dghteraod->GetPdgCode() == -3334 && mcGdMotherNegV0Dghteraod->GetPdgCode() == -3334) {lAssoOmegaPlus = kTRUE;  cascadeMass = 1.672;}
           }

           // ---------------------------------
           // - Miscellaneous pieces of info that may help regarding data quality assessment.
           lXiMomX = xi->MomXiX();
           lXiMomY = xi->MomXiY();
           lXiMomZ = xi->MomXiZ();
           lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
           lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
           Double_t lV0MomX = xi->MomV0X();
           Double_t lV0MomY = xi->MomV0Y();
           Double_t lV0MomZ = xi->MomV0Z();
           lV0TotMom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2)+TMath::Power(lV0MomZ,2));
           lBachMomX = xi->MomBachX();
           lBachMomY = xi->MomBachY();
           lBachMomZ = xi->MomBachZ();		
           lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
           lV0NMomX = xi->MomNegX();
           lV0NMomY = xi->MomNegY();
           lV0PMomX = xi->MomPosX();
           lV0PMomY = xi->MomPosY(); 
           lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
           lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
           lV0toXiCosineOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
           lRapXi    = xi->RapXi();
           lRapOmega = xi->RapOmega();

    }// end of AOD treatment

    // ---------------------------------------
    // Cut on pt of the three daughter tracks
    if (lBachTransvMom<fMinPtCutOnDaughterTracks)   { AliWarning("ERROR: bachelor pT < lowlimit");          continue; }
    if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) { AliWarning("ERROR: positive daughter pT < lowlimit"); continue; }
    if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) { AliWarning("ERROR: negative daughter pT < lowlimit"); continue; }

    // ---------------------------------------------------
    // Cut on pseudorapidity of the three daughter tracks
    if (TMath::Abs(etaBach) > 0.8) { AliWarning("ERROR: bachelor eta > maxlimit");          continue; }
    if (TMath::Abs(etaPos)  > 0.8) { AliWarning("ERROR: positive daughter eta > maxlimit"); continue; }
    if (TMath::Abs(etaNeg)  > 0.8) { AliWarning("ERROR: negative daughter eta > maxlimit"); continue; }

    // ----------------------------------
    // Calculate proper time for cascade
    if (!fisMC) {
      if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
           ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) cascadeMass = 1.321;
      if ( ( (lChargeXi<0) && lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
           ( (lChargeXi>0) && lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) cascadeMass = 1.672; 
    }
    Double_t lctau = TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power(( lPosXi[2]-lBestPrimaryVtxPos[2]),2));
    if (lXiTotMom != 0) lctau = lctau*cascadeMass/lXiTotMom;
    else                lctau = -1.;
    // Calculate proper time for Lambda (reconstructed)
    Float_t lambdaMass = 1.115683; // PDG mass
    Float_t distV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2));
    Float_t lctauV0 = -1.;
    if (lV0TotMom != 0) lctauV0 = distV0Xi*lambdaMass/lV0TotMom;


    // ------------------- 
    // Fill the Histograms
    if (lChargeXi < 0) { 
        if ((!fisMC && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) || (fisMC && lAssoXiMinus)) {
            fHistVarDcaCascDaughtXiMinus->Fill(lXiTransvMom,lDcaXiDaughters);
            fHistVarDcaBachToPrimVertexXiMinus->Fill(lXiTransvMom,lDcaBachToPrimVertexXi);
            fHistVarCascCosineOfPointingAngleXiMinus->Fill(lXiTransvMom,lXiCosineOfPointingAngle);
            fHistVarCascRadiusXiMinus->Fill(lXiTransvMom,lXiRadius);
            fHistVarInvMassLambdaAsCascDghterXiMinus->Fill(lXiTransvMom,lInvMassLambdaAsCascDghter);
            fHistVarDcaV0DaughtersXiMinus->Fill(lXiTransvMom,lDcaV0DaughtersXi);
            fHistVarV0CosineOfPAToCascVertexXiMinus->Fill(lXiTransvMom,lV0toXiCosineOfPointingAngle);
            fHistVarV0RadiusXiMinus->Fill(lXiTransvMom,lV0RadiusXi);
            fHistVarDcaV0ToPrimVertexXiMinus->Fill(lXiTransvMom,lDcaV0ToPrimVertexXi);
            fHistVarDcaPosToPrimVertexXiMinus->Fill(lXiTransvMom,lDcaPosToPrimVertexXi);
            fHistVarDcaNegToPrimVertexXiMinus->Fill(lXiTransvMom,lDcaNegToPrimVertexXi);
            fHistMassXiMinus->Fill( lInvMassXiMinus );
            fHistVarTransvMomentumXiMinus->Fill(lXiTransvMom);
            fHistVarRapidityXiMinus->Fill(lRapXi);
            fHistVarCascProperLengthXiMinus->Fill(lXiTransvMom,lctau);
            fHistVarV0ProperLengthXiMinus->Fill(lXiTransvMom,lctauV0);
        }
        if ((!fisMC && lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) || (fisMC && lAssoOmegaMinus)) {
            fHistVarDcaCascDaughtOmegaMinus->Fill(lXiTransvMom,lDcaXiDaughters);
            fHistVarDcaBachToPrimVertexOmegaMinus->Fill(lXiTransvMom,lDcaBachToPrimVertexXi);
            fHistVarCascCosineOfPointingAngleOmegaMinus->Fill(lXiTransvMom,lXiCosineOfPointingAngle);
            fHistVarCascRadiusOmegaMinus->Fill(lXiTransvMom,lXiRadius);
            fHistVarInvMassLambdaAsCascDghterOmegaMinus->Fill(lXiTransvMom,lInvMassLambdaAsCascDghter);
            fHistVarDcaV0DaughtersOmegaMinus->Fill(lXiTransvMom,lDcaV0DaughtersXi);
            fHistVarV0CosineOfPAToCascVertexOmegaMinus->Fill(lXiTransvMom,lV0toXiCosineOfPointingAngle);
            fHistVarV0RadiusOmegaMinus->Fill(lXiTransvMom,lV0RadiusXi);
            fHistVarDcaV0ToPrimVertexOmegaMinus->Fill(lXiTransvMom,lDcaV0ToPrimVertexXi);
            fHistVarDcaPosToPrimVertexOmegaMinus->Fill(lXiTransvMom,lDcaPosToPrimVertexXi);
            fHistVarDcaNegToPrimVertexOmegaMinus->Fill(lXiTransvMom,lDcaNegToPrimVertexXi);
            fHistMassOmegaMinus->Fill( lInvMassOmegaMinus );
            fHistVarTransvMomentumOmegaMinus->Fill(lXiTransvMom);
            fHistVarRapidityOmegaMinus->Fill(lRapXi);
            fHistVarCascProperLengthOmegaMinus->Fill(lXiTransvMom,lctau);
            fHistVarV0ProperLengthOmegaMinus->Fill(lXiTransvMom,lctauV0); 
        }
    } else {
        if ((!fisMC && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) || (fisMC && lAssoXiPlus)) {
            fHistVarDcaCascDaughtXiPlus->Fill(lXiTransvMom,lDcaXiDaughters);
            fHistVarDcaBachToPrimVertexXiPlus->Fill(lXiTransvMom,lDcaBachToPrimVertexXi);
            fHistVarCascCosineOfPointingAngleXiPlus->Fill(lXiTransvMom,lXiCosineOfPointingAngle);
            fHistVarCascRadiusXiPlus->Fill(lXiTransvMom,lXiRadius);
            fHistVarInvMassLambdaAsCascDghterXiPlus->Fill(lXiTransvMom,lInvMassLambdaAsCascDghter);
            fHistVarDcaV0DaughtersXiPlus->Fill(lXiTransvMom,lDcaV0DaughtersXi);
            fHistVarV0CosineOfPAToCascVertexXiPlus->Fill(lXiTransvMom,lV0toXiCosineOfPointingAngle);
            fHistVarV0RadiusXiPlus->Fill(lXiTransvMom,lV0RadiusXi);
            fHistVarDcaV0ToPrimVertexXiPlus->Fill(lXiTransvMom,lDcaV0ToPrimVertexXi);
            fHistVarDcaPosToPrimVertexXiPlus->Fill(lXiTransvMom,lDcaPosToPrimVertexXi);
            fHistVarDcaNegToPrimVertexXiPlus->Fill(lXiTransvMom,lDcaNegToPrimVertexXi);
            fHistMassXiPlus->Fill( lInvMassXiPlus );
            fHistVarTransvMomentumXiPlus->Fill(lXiTransvMom);
            fHistVarRapidityXiPlus->Fill(lRapXi);
            fHistVarCascProperLengthXiPlus->Fill(lXiTransvMom,lctau);
            fHistVarV0ProperLengthXiPlus->Fill(lXiTransvMom,lctauV0);
        }
        if ((!fisMC && lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) || (fisMC && lAssoOmegaPlus)) {
            fHistVarDcaCascDaughtOmegaPlus->Fill(lXiTransvMom,lDcaXiDaughters);
            fHistVarDcaBachToPrimVertexOmegaPlus->Fill(lXiTransvMom,lDcaBachToPrimVertexXi);
            fHistVarCascCosineOfPointingAngleOmegaPlus->Fill(lXiTransvMom,lXiCosineOfPointingAngle);
            fHistVarCascRadiusOmegaPlus->Fill(lXiTransvMom,lXiRadius);
            fHistVarInvMassLambdaAsCascDghterOmegaPlus->Fill(lXiTransvMom,lInvMassLambdaAsCascDghter);
            fHistVarDcaV0DaughtersOmegaPlus->Fill(lXiTransvMom,lDcaV0DaughtersXi);
            fHistVarV0CosineOfPAToCascVertexOmegaPlus->Fill(lXiTransvMom,lV0toXiCosineOfPointingAngle);
            fHistVarV0RadiusOmegaPlus->Fill(lXiTransvMom,lV0RadiusXi);
            fHistVarDcaV0ToPrimVertexOmegaPlus->Fill(lXiTransvMom,lDcaV0ToPrimVertexXi);
            fHistVarDcaPosToPrimVertexOmegaPlus->Fill(lXiTransvMom,lDcaPosToPrimVertexXi);
            fHistVarDcaNegToPrimVertexOmegaPlus->Fill(lXiTransvMom,lDcaNegToPrimVertexXi);
            fHistMassOmegaPlus->Fill( lInvMassOmegaPlus );
            fHistVarTransvMomentumOmegaPlus->Fill(lXiTransvMom);
            fHistVarRapidityOmegaPlus->Fill(lRapXi);
            fHistVarCascProperLengthOmegaPlus->Fill(lXiTransvMom,lctau);
            fHistVarV0ProperLengthOmegaPlus->Fill(lXiTransvMom,lctauV0);
        }
    }
    
  }// end of the Cascade loop (ESD or AOD)
    
  
  // Post output data.
  PostData(1, fListHistMultistrangeQA);

}// End UserExec


//-------------------------------------------------------
void AliAnalysisTaskQAMultistrange::Terminate(Option_t *) 

{

}
