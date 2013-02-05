/**************************************************************************
 *  Authors : Antonin Maire, Boris Hippolyte                              *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//            AliAnalysisTaskCheckCascadepp276 class
//
//            Origin AliAnalysisTaskCheckCascade which has four roles :
//              1. QAing the Cascades from ESD and AOD
//                 Origin:  AliAnalysisTaskESDCheckV0 by Boris Hippolyte Nov2007, hippolyt@in2p3.fr
//              2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//              3. Supply an AliCFContainer meant to define the optimised topological selections
//              4. Rough azimuthal correlation study (Eta, Phi)
//              Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//              Modified :           A.Maire Mar2010 
//
//              Adapted to PbPb analysis: M. Nicassio, maria.nicassio@ba.infn.it
//               Feb-August2011
//                - Physics selection moved to the run.C macro
//                - Centrality selection added (+ setters) and histos
//                - flag and setters added (CF container usage, vertex range)
//                - histo added and histo/container binning changed 
//                - protection in the destructor for CAF usage          
//                - AliWarning disabled
//                - number of tracklets from AOD also          
//                - automatic settings for PID
//               September2011
//                - proper time histos/container added (V0 and Cascades)
//                - cosine PA V0 wrt Xi vertex in the container  
//               November2011
//                - re-run V0's and cascade's vertexers (SetCuts instead SetDefaultCuts!!)
//                - problems of libraries on Grid --> code copied in the task (from AliRoot v5-10-AN
//                  where new pt dependent V0's cosPA cut implemented by Iouri) 
//                - AOD analysis part completed 
//
//
//              Adapted to pp 2.76 analysis: D. Colella, domenico.colella@ba.infn.it
//               Gen-now 2012
//                - Physics selection re-moved here (mainly for normalization in the efficiency calcuation)
//                - Centrality selection deleted
//                - 
//
//-----------------------------------------------------------------

class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"

#include "AliESDVZERO.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h" 
#include "AliAODInputHandler.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAODTrack.h"

#include "AliAnalysisTaskCheckCascadepp276.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCheckCascadepp276)



//________________________________________________________________________
AliAnalysisTaskCheckCascadepp276::AliAnalysisTaskCheckCascadepp276() 
  : AliAnalysisTaskSE(), 
    fAnalysisType               ("ESD"),
    fESDtrackCuts               (0),
    fPIDResponse                (0),
    fkRerunV0CascVertexers      (0),
    fkSDDSelectionOn            (kTRUE),
    fkQualityCutZprimVtxPos     (kTRUE),
    fkQualityCutNoTPConlyPrimVtx(kTRUE),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCutnTPCcls         (kTRUE),
    fkQualityCutPileup          (kTRUE),
    fwithSDD                    (kTRUE),
    fMinnTPCcls                 (0),
    fkExtraSelections           (0),
    fVtxRange                   (0),
    fVtxRangeMin                (0),
    fMinPtCutOnDaughterTracks   (0),
    fEtaCutOnDaughterTracks     (0),

    // - Plots initialisation
    fListHistCascade(0),
      // Cascades multiplicity plots
      fHistCascadeMultiplicityBeforeAnySel(0),
      fHistCascadeMultiplicityAfterSDDSel(0),
      fHistCascadeMultiplicityAfterPhysicsSel(0),
      fHistCascadeMultiplicityForSelEvtNoTPCOnly(0),
      fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
      fHistCascadeMultiplicityAfterVertexCutSel(0),
      // Tracks multiplicity plots
      fHistTrackMultiplicityBeforeAnySel(0),
      fHistTrackMultiplicityAfterSDDSel(0),
      fHistTrackMultiplicityAfterPhysicsSel(0),
      fHistTrackMultiplicityForSelEvtNoTPCOnly(0),
      fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
      fHistTrackMultiplicityAfterVertexCutSel(0),
      // Vertex position plots (BestVertex)
      fHistPVx(0), fHistPVy(0), fHistPVz(0),
      fHistPVxAnalysis(0), fHistPVyAnalysis(0), fHistPVzAnalysis(0),    
      // TPC cluster distributions for daughters
      fHistPosV0TPCClusters(0), 
      fHistNegV0TPCClusters(0), 
      fHistBachTPCClusters(0),
      // Cut's variables distributions
      fHistEffMassXi(0), 
      fHistDcaXiDaughters(0), 
      fHistDcaBachToPrimVertex(0), 
      fHistXiCosineOfPointingAngle(0), 
      fHistXiRadius(0),
      fHistMassLambdaAsCascDghter(0),
      fHistDcaV0DaughtersXi(0),
      fHistDcaV0ToPrimVertexXi(0), 
      fHistV0CosineOfPointingAngleXi(0),
      fHistV0RadiusXi(0),
      fHistDcaPosToPrimVertexXi(0), 
      fHistDcaNegToPrimVertexXi(0), 
      // Invariant mass distributions
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      // Transverse and total momentum distributions
      fHistXiTransvMom(0), fHistXiTotMom(0), fHistBachTransvMomXi(0), fHistBachTotMomXi(0),
      // Others QA plots
      fHistChargeXi(0),
      fHistV0toXiCosineOfPointingAngle(0),
      fHistRapXi(0), fHistRapOmega(0), 
      fHistEtaXi(0), fHistEtaBachXi(0), fHistEtaPosXi(0), fHistEtaNegXi(0),
      fHistThetaXi(0), 
      fHistPhiXi(0),
      f2dHistArmenteros(0),			
      f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
      f2dHistEffMassLambdaVsEffMassXiPlus(0),  f2dHistEffMassXiVsEffMassOmegaPlus(0),
      f2dHistXiRadiusVsEffMassXiMinus(0),      f2dHistXiRadiusVsEffMassXiPlus(0),
      f2dHistXiRadiusVsEffMassOmegaMinus(0),   f2dHistXiRadiusVsEffMassOmegaPlus(0),
      f2dHistTPCdEdxOfCascDghters(0),
      f2dHistDcaXiDaughtersvsInvMass(0), 
      f2dHistDcaBachToPrimVertexvsInvMass(0), 
      f2dHistXiCosineOfPointingAnglevsInvMass(0),
      f2dHistMassLambdaAsCascDghtervsInvMass(0),
      f2dHistDcaV0DaughtersXivsInvMass(0),
      f2dHistDcaV0ToPrimVertexXivsInvMass(0),
      // Containers for cuts study 
      fCFContCascadePIDXiMinus(0),
      fCFContCascadePIDXiPlus(0),
      fCFContCascadePIDOmegaMinus(0),
      fCFContCascadePIDOmegaPlus(0),
      fCFContCascadeCuts(0)
    
    {
     // Dummy Constructor
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascSels        [iCascSelIdx ] = -1.; }
    }


//________________________________________________________________________
AliAnalysisTaskCheckCascadepp276::AliAnalysisTaskCheckCascadepp276(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAnalysisType               ("ESD"), 
    fESDtrackCuts               (0), 
    fPIDResponse                (0),
    fkRerunV0CascVertexers      (0),
    fkSDDSelectionOn            (kTRUE),
    fkQualityCutZprimVtxPos     (kTRUE),
    fkQualityCutNoTPConlyPrimVtx(kTRUE),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCutnTPCcls         (kTRUE),
    fkQualityCutPileup          (kTRUE),
    fwithSDD                    (kTRUE),
    fMinnTPCcls                 (0),
    fkExtraSelections           (0),
    fVtxRange                   (0),
    fVtxRangeMin                (0),
    fMinPtCutOnDaughterTracks   (0),
    fEtaCutOnDaughterTracks     (0),
     
    // - Plots initialisation
    fListHistCascade(0),

      // Cascades multiplicity plots
      fHistCascadeMultiplicityBeforeAnySel(0),
      fHistCascadeMultiplicityAfterSDDSel(0),
      fHistCascadeMultiplicityAfterPhysicsSel(0),
      fHistCascadeMultiplicityForSelEvtNoTPCOnly(0),
      fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
      fHistCascadeMultiplicityAfterVertexCutSel(0),
      // Tracks multiplicity plots
      fHistTrackMultiplicityBeforeAnySel(0),
      fHistTrackMultiplicityAfterSDDSel(0),
      fHistTrackMultiplicityAfterPhysicsSel(0),
      fHistTrackMultiplicityForSelEvtNoTPCOnly(0),
      fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup(0),
      fHistTrackMultiplicityAfterVertexCutSel(0),
      // Vertex position plots (BestVertex)
      fHistPVx(0), fHistPVy(0), fHistPVz(0),
      fHistPVxAnalysis(0), fHistPVyAnalysis(0), fHistPVzAnalysis(0),
      // TPC cluster distributions for daughters
      fHistPosV0TPCClusters(0), fHistNegV0TPCClusters(0), fHistBachTPCClusters(0),
      // Cut's variables distributions
      fHistEffMassXi(0),
      fHistDcaXiDaughters(0),
      fHistDcaBachToPrimVertex(0),
      fHistXiCosineOfPointingAngle(0),
      fHistXiRadius(0),
      fHistMassLambdaAsCascDghter(0),
      fHistDcaV0DaughtersXi(0),
      fHistDcaV0ToPrimVertexXi(0),
      fHistV0CosineOfPointingAngleXi(0),
      fHistV0RadiusXi(0),
      fHistDcaPosToPrimVertexXi(0),
      fHistDcaNegToPrimVertexXi(0),
      // Invariant mass distributions
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      // Transverse and total momentum distributions
      fHistXiTransvMom(0), fHistXiTotMom(0), fHistBachTransvMomXi(0), fHistBachTotMomXi(0),
      // Others QA plots
      fHistChargeXi(0),
      fHistV0toXiCosineOfPointingAngle(0),
      fHistRapXi(0), fHistRapOmega(0),
      fHistEtaXi(0), fHistEtaBachXi(0), fHistEtaPosXi(0), fHistEtaNegXi(0),
      fHistThetaXi(0),
      fHistPhiXi(0),
      f2dHistArmenteros(0),
      f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
      f2dHistEffMassLambdaVsEffMassXiPlus(0),  f2dHistEffMassXiVsEffMassOmegaPlus(0),
      f2dHistXiRadiusVsEffMassXiMinus(0),      f2dHistXiRadiusVsEffMassXiPlus(0),
      f2dHistXiRadiusVsEffMassOmegaMinus(0),   f2dHistXiRadiusVsEffMassOmegaPlus(0),
      f2dHistTPCdEdxOfCascDghters(0),
      f2dHistDcaXiDaughtersvsInvMass(0),
      f2dHistDcaBachToPrimVertexvsInvMass(0),
      f2dHistXiCosineOfPointingAnglevsInvMass(0),
      f2dHistMassLambdaAsCascDghtervsInvMass(0),
      f2dHistDcaV0DaughtersXivsInvMass(0),
      f2dHistDcaV0ToPrimVertexXivsInvMass(0),
      // Containers for cuts study 
      fCFContCascadePIDXiMinus(0),
      fCFContCascadePIDXiPlus(0),
      fCFContCascadePIDOmegaMinus(0),
      fCFContCascadePIDOmegaPlus(0),
      fCFContCascadeCuts(0)
    
    //_____Costructor____
    {
     // Define input and output slots here
     // Input slot #0 works with a TChain
     // DefineInput(0, TChain::Class());
     // Output slot #1 writes into a TList container (cascade)
        // default p-p values
        fV0Sels[0] =  33.  ;     // max allowed chi2
        fV0Sels[1] =   0.073;    // min allowed impact parameter for the 1st daughter 
        fV0Sels[2] =   0.073;    // min allowed impact parameter for the 2nd daughter 
        fV0Sels[3] =   1.18;     // max allowed DCA between the daughter tracks       
        fV0Sels[4] =    .983;    // min allowed cosine of V0's pointing angle         
        fV0Sels[5] =   2.67;     // min radius of the fiducial volume                 
        fV0Sels[6] = 100.;       // max radius of the fiducial volume                 

        fCascSels[0] =  33.;     // max allowed chi2 (same as PDC07)
        fCascSels[1] =   0.03;   // min allowed V0 impact parameter                    
        fCascSels[2] =   0.008;  // "window" around the Lambda mass                    
        fCascSels[3] =   0.0204; // min allowed bachelor's impact parameter          
        fCascSels[4] =   1.68;   // max allowed DCA between the V0 and the bachelor    
        fCascSels[5] =   0.9826; // min allowed cosine of the cascade pointing angle   
        fCascSels[6] =   0.38;   // min radius of the fiducial volume                  
        fCascSels[7] = 100.;     // max radius of the fiducial volume                  

     // Output slot #0 writes into a TList container (Cascade)
     DefineOutput(1, TList::Class());
     DefineOutput(2, AliCFContainer::Class());
     DefineOutput(3, AliCFContainer::Class());
     DefineOutput(4, AliCFContainer::Class());
     DefineOutput(5, AliCFContainer::Class());
     DefineOutput(6, AliCFContainer::Class());
     AliLog::SetClassDebugLevel("AliAnalysisTaskCheckCascadepp276",1);
    } 


    //_____Destructor_____
    AliAnalysisTaskCheckCascadepp276::~AliAnalysisTaskCheckCascadepp276() {
      // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
      // They will be deleted when fListCascade is deleted by the TSelector dtor
      // Because of TList::SetOwner() ...   
       if (fListHistCascade && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())           { delete fListHistCascade; fListHistCascade = 0x0; }
       if (fCFContCascadePIDXiMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())   { delete fCFContCascadePIDXiMinus; fCFContCascadePIDXiMinus = 0x0; }
       if (fCFContCascadePIDXiPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())    { delete fCFContCascadePIDXiPlus; fCFContCascadePIDXiPlus = 0x0; }
       if (fCFContCascadePIDOmegaMinus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){ delete fCFContCascadePIDOmegaMinus; fCFContCascadePIDOmegaMinus = 0x0; }
       if (fCFContCascadePIDOmegaPlus && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadePIDOmegaPlus; fCFContCascadePIDOmegaPlus = 0x0; }  
       if (fCFContCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())         { delete fCFContCascadeCuts; fCFContCascadeCuts = 0x0; }
       if (fESDtrackCuts)                                                                          { delete fESDtrackCuts; fESDtrackCuts = 0x0; }
    }


//________________________________________________________________________
void AliAnalysisTaskCheckCascadepp276::UserCreateOutputObjects() {
  // Create histograms
  // Called once


 fListHistCascade = new TList();
 fListHistCascade->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

 //-----------------------------------------------
 // Particle Identification Setup (new PID object)
 //-----------------------------------------------
 AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
 AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
 fPIDResponse = inputHandler->GetPIDResponse();

 // Only used to get the number of primary reconstructed tracks
 if (fAnalysisType == "ESD" && (! fESDtrackCuts )){
   fESDtrackCuts = new AliESDtrackCuts();
 }

 //---------------------------------------------------
 // Initialize cuts to re-run V0 and cascade vertexers
 //---------------------------------------------------
 // Not validated; to be checked
 fV0Sels[0] =  33.  ;     // max allowed chi2
 fV0Sels[1] =   0.073;    // min allowed impact parameter for the 1st daughter 
 fV0Sels[2] =   0.073;    // min allowed impact parameter for the 2nd daughter 
 fV0Sels[3] =   1.18;     // max allowed DCA between the daughter tracks       
 fV0Sels[4] =    .983;    // min allowed cosine of V0's pointing angle         
 fV0Sels[5] =   2.67;     // min radius of the fiducial volume                 
 fV0Sels[6] = 100.;       // max radius of the fiducial volume                 

 fCascSels[0] =  33.;     // max allowed chi2 (same as PDC07)
 fCascSels[1] =   0.03;   // min allowed V0 impact parameter                    
 fCascSels[2] =   0.008;  // "window" around the Lambda mass                    
 fCascSels[3] =   0.0204; // min allowed bachelor's impact parameter           //check cuts 
 fCascSels[4] =   1.68;   // max allowed DCA between the V0 and the bachelor    
 fCascSels[5] =   0.9826; // min allowed cosine of the cascade pointing angle   
 fCascSels[6] =   0.38;   // min radius of the fiducial volume                  
 fCascSels[7] = 100.;     // max radius of the fiducial volume 

 //----------------------
 // Initialize the histos
 //----------------------
 
 // - Cascades multiplicity plots 
 if(! fHistCascadeMultiplicityBeforeAnySel) {
        fHistCascadeMultiplicityBeforeAnySel = new TH1F("fHistCascadeMultiplicityBeforeAnySel",
                        "Cascades per event (before any selections);Nbr of Cascades/Evt;Events",                
                        50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityBeforeAnySel);
 }
 if(! fHistCascadeMultiplicityAfterSDDSel) {
        fHistCascadeMultiplicityAfterSDDSel = new TH1F("fHistCascadeMultiplicityAfterSDDSel", 
			"Cascades per event (after the SDD selection);Nbr of Cascades/Evt;Events", 
			50, 0, 50); 		
	fListHistCascade->Add(fHistCascadeMultiplicityAfterSDDSel);
 }
 if(! fHistCascadeMultiplicityAfterPhysicsSel) {
        fHistCascadeMultiplicityAfterPhysicsSel = new TH1F("fHistCascadeMultiplicityAfterPhysicsSel",
                        "Cascades per event (after physics selection);Nbr of Cascades/Evt;Events",
                        50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterPhysicsSel);
 }
 if(! fHistCascadeMultiplicityForSelEvtNoTPCOnly) {
        fHistCascadeMultiplicityForSelEvtNoTPCOnly = new TH1F("fHistCascadeMultiplicityForSelEvtNoTPCOnly",
                        "Cascades per event (for selected events with well-established PV);Nbr of Cascades/Evt;Events",
                        50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityForSelEvtNoTPCOnly);
 }
 if(! fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup) {
        fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup",
                        "Cascades per event (for selected events with well-establisched PV and no pile-up);Nbr of Cascades/Evt;Events",
                        50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup);
 }
 if(! fHistCascadeMultiplicityAfterVertexCutSel) {
        fHistCascadeMultiplicityAfterVertexCutSel = new TH1F("fHistCascadeMultiplicityAfterVertexCutSel",
                                                             "Cascades per event (after vertex cut selection);Nbr of Cascades/Evt;Events",
                                                             50, 0, 50);
        fListHistCascade->Add(fHistCascadeMultiplicityAfterVertexCutSel);
 }
 // - Tracks multiplicity plots 
 if(! fHistTrackMultiplicityBeforeAnySel) {
	fHistTrackMultiplicityBeforeAnySel = new TH1F("fHistTrackMultiplicityBeforeAnySel", 
			"Tracks per event (before any selections);Nbr of Cascades/Evt;Events", 
			200, 0, 200); 		
	fListHistCascade->Add(fHistTrackMultiplicityBeforeAnySel);
 } 
 if(! fHistTrackMultiplicityAfterSDDSel) {
        fHistTrackMultiplicityAfterSDDSel = new TH1F("fHistTrackMultiplicityAfterSDDSel",                  
                        "Tracks per event (after the SDD selection);Nbr of Cascades/Evt;Events",
                        200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityAfterSDDSel);
 }
 if(! fHistTrackMultiplicityAfterPhysicsSel) {
        fHistTrackMultiplicityAfterPhysicsSel = new TH1F("fHistTrackMultiplicityAfterPhysicsSel",
                        "Tracks per event (after physics selection);Nbr of Cascades/Evt;Events",
                        200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityAfterPhysicsSel);
 }
 if(! fHistTrackMultiplicityForSelEvtNoTPCOnly) {
        fHistTrackMultiplicityForSelEvtNoTPCOnly = new TH1F("fHistTrackMultiplicityForSelEvtNoTPCOnly",
                        "Tracks per event (for selected events with well-established PV);Nbr of Cascades/Evt;Events",
                        200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityForSelEvtNoTPCOnly);
 }
 if(! fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup) {
        fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup",
                        "Tracks per event (for selected events with well-establisched PV and no pile-up);Nbr of Cascades/Evt;Events",
                        200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup);
 }
 if(! fHistTrackMultiplicityAfterVertexCutSel) {
        fHistTrackMultiplicityAfterVertexCutSel = new TH1F("fHistTrackMultiplicityAfterVertexCutSel",
                                                           "Tracks per event (after vertex cut selection);Nbr of Cascades/Evt;Events",
                                                           200, 0, 200);
        fListHistCascade->Add(fHistTrackMultiplicityAfterVertexCutSel);
 }
 // - Vertex position plots
 if(! fHistPVx ){
	fHistPVx = new TH1F("fHistPVx", "Best PV position in x; x (cm); Events", 2000, -0.5, 0.5);
	fListHistCascade->Add(fHistPVx);
 }
 if(! fHistPVy ){
        fHistPVy = new TH1F("fHistPVy", "Best PV position in y; y (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVy);
 }
 if(! fHistPVz ){
        fHistPVz = new TH1F("fHistPVz", "Best PV position in z; z (cm); Events", 400, -20, 20);
        fListHistCascade->Add(fHistPVz);
 } 
 if(! fHistPVxAnalysis ){
        fHistPVxAnalysis = new TH1F("fHistPVxAnalysis", "Best PV position in x (after events selections); x (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVxAnalysis);
 }
 if(! fHistPVyAnalysis ){
        fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", "Best PV position in y (after events selections); y (cm); Events", 2000, -0.5, 0.5);
        fListHistCascade->Add(fHistPVyAnalysis);
 }
 if(! fHistPVzAnalysis ){
        fHistPVzAnalysis = new TH1F("fHistPVzAnalysis", "Best PV position in z (after events selections); z (cm); Events", 400, -20, 20);
        fListHistCascade->Add(fHistPVzAnalysis);
 }
 // - TPC clusetr sdistributions for daughters (histos for events containing at least ONE CASCADE)
 if(! fHistPosV0TPCClusters ){
        fHistPosV0TPCClusters = new TH1F("fHistPosV0TPCClusters", "TPC clusters for Pos. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Pos.); Track counts", 165, 0.0, 165.0);
        fListHistCascade->Add(fHistPosV0TPCClusters);
 }
 if(! fHistNegV0TPCClusters ){
        fHistNegV0TPCClusters = new TH1F("fHistNegV0TPCClusters", "TPC clusters for Neg. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Neg.); Track counts", 165, 0.0, 165.0);
        fListHistCascade->Add(fHistNegV0TPCClusters);
 }
 if(! fHistBachTPCClusters ){
        fHistBachTPCClusters = new TH1F("fHistBachTPCClusters", "TPC clusters for Bachelor track; Nbr of TPC clusters (Bach); Track counts", 165, 0.0, 165.0);
        fListHistCascade->Add(fHistBachTPCClusters);
 }
 // - Cut's variables distributions (typical histos for cascades): as example only for the Xi (both particle and anti-particle)
 if(! fHistEffMassXi) {
     fHistEffMassXi = new TH1F("fHistEffMassXi", "Xi candidates; Invariant Mass (GeV/c^{2}); Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistEffMassXi);
 }  
 if(! fHistDcaXiDaughters ){
     fHistDcaXiDaughters = new TH1F("fHistDcaXiDaughters", "DCA between Xi daughters; DCA (cm); Counts", 210, 0., 2.1);
     fListHistCascade->Add(fHistDcaXiDaughters);
 }
 if(! fHistDcaBachToPrimVertex) {
     fHistDcaBachToPrimVertex = new TH1F("fHistDcaBachToPrimVertex", "Impact parameter of Bach. to Prim. Vertex; DCA (cm); Counts", 250, 0., 0.25);
     fListHistCascade->Add(fHistDcaBachToPrimVertex);
 }
 if(! fHistXiCosineOfPointingAngle) {
     fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl); Counts", 301, 0.97, 1.0001);
     fListHistCascade->Add(fHistXiCosineOfPointingAngle);
 }
 if(! fHistXiRadius ){
     fHistXiRadius = new TH1F("fHistXiRadius", "Cascade decay transv. radius; r (cm); Counts" , 1050, 0., 105.0);
     fListHistCascade->Add(fHistXiRadius);
 }
 if(! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter", "#Lambda associated to cascade candidates; Eff. Mass (GeV/c^{2}); Counts", 300, 1.0, 1.3);
     fListHistCascade->Add(fHistMassLambdaAsCascDghter);
 }
 if(! fHistDcaV0DaughtersXi) {
     fHistDcaV0DaughtersXi = new TH1F("fHistDcaV0DaughtersXi", "DCA between V0 daughters, in cascade; DCA (cm); Counts", 320, 0., 1.6);
     fListHistCascade->Add(fHistDcaV0DaughtersXi);
 }
 if(! fHistDcaV0ToPrimVertexXi) {
     fHistDcaV0ToPrimVertexXi = new TH1F("fHistDcaV0ToPrimVertexXi", "Impact parameter of V0  to Prim. Vertex, in cascade; DCA (cm); Counts", 200, 0., 1.);
     fListHistCascade->Add(fHistDcaV0ToPrimVertexXi);
 }
 if(! fHistV0CosineOfPointingAngleXi) {
     fHistV0CosineOfPointingAngleXi = new TH1F("fHistV0CosineOfPointingAngleXi", "Cosine of V0 Pointing Angle, in cascade; Cos(V0 Point. Angl); Counts", 201, 0.8, 1.001);
     fListHistCascade->Add(fHistV0CosineOfPointingAngleXi);
 }
 if(! fHistV0RadiusXi) {
     fHistV0RadiusXi = new TH1F("fHistV0RadiusXi", "V0 decay radius, in cascade; radius (cm); Counts", 1050, 0., 105.0);
     fListHistCascade->Add(fHistV0RadiusXi);
 }
 if(! fHistDcaPosToPrimVertexXi) {
     fHistDcaPosToPrimVertexXi = new TH1F("fHistDcaPosToPrimVertexXi", "Impact parameter of V0 pos daughter to Prim. Vertex; DCA (cm); Counts", 300, 0, 3);
     fListHistCascade->Add(fHistDcaPosToPrimVertexXi);
 }
 if(! fHistDcaNegToPrimVertexXi) {
     fHistDcaNegToPrimVertexXi = new TH1F("fHistDcaNegToPrimVertexXi", "Impact parameter of V0 neg daughter to Prim. Vertex; DCA (cm); Counts", 300, 0, 3);
     fListHistCascade->Add(fHistDcaNegToPrimVertexXi);
 }
 // - Effective mass histos for cascades.
    //By cascade hyp  
 if(! fHistMassXiMinus) {
     fHistMassXiMinus = new TH1F("fHistMassXiMinus", "#Xi^{-} candidates; M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistMassXiMinus);
 } 
 if(! fHistMassXiPlus) {
     fHistMassXiPlus = new TH1F("fHistMassXiPlus", "#Xi^{+} candidates; M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistMassXiPlus);
 }
 if(! fHistMassOmegaMinus) {
     fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus", "#Omega^{-} candidates; M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500, 1.5, 2.5);
     fListHistCascade->Add(fHistMassOmegaMinus);
 }
 if(! fHistMassOmegaPlus) {
     fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus", "#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2}); Counts", 500, 1.5, 2.5);
     fListHistCascade->Add(fHistMassOmegaPlus);
 }
 // - Transverse and total momentum distributions
 if(! fHistXiTransvMom ){
     fHistXiTransvMom = new TH1F("fHistXiTransvMom", "#Xi transverse momentum (cand. around the mass peak); p_{t}(#Xi) (GeV/c); Counts", 100, 0.0, 10.0);
     fListHistCascade->Add(fHistXiTransvMom);
 }
 if(! fHistXiTotMom ){
     fHistXiTotMom = new TH1F("fHistXiTotMom", "#Xi momentum norm (cand. around the mass peak); p_{tot}(#Xi) (GeV/c); Counts", 150, 0.0, 15.0);
     fListHistCascade->Add(fHistXiTotMom);
 }
 if(! fHistBachTransvMomXi ){
     fHistBachTransvMomXi = new TH1F("fHistBachTransvMomXi", "#Xi Bach. transverse momentum (cand. around the mass peak); p_{t}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
     fListHistCascade->Add(fHistBachTransvMomXi);
 }
 if(! fHistBachTotMomXi ){
     fHistBachTotMomXi = new TH1F("fHistBachTotMomXi", "#Xi Bach. momentum norm (cand. around the mass peak); p_{tot}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
     fListHistCascade->Add(fHistBachTotMomXi);
 }
 // - Others QA plots
    //TH1
 if(! fHistChargeXi ){
     fHistChargeXi = new TH1F("fHistChargeXi", "Charge of Xi candidates; Sign; Counts", 5, -2.0, 3.0);
     fListHistCascade->Add(fHistChargeXi);
 }
 if(! fHistV0toXiCosineOfPointingAngle) {
     fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ; Cos(V0 Point. Angl / Xi vtx); Counts", 1101, 0.89, 1.0001);
     fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
 }
 if(! fHistRapXi ){
     fHistRapXi = new TH1F("fHistRapXi", "Rapidity of #Xi candidates (around the mass peak); y; Counts", 20, -1.0, 1.0);
     fListHistCascade->Add(fHistRapXi);
 }
 if(! fHistRapOmega ){
     fHistRapOmega = new TH1F("fHistRapOmega", "Rapidity of #Omega candidates (around the mass peak); y; Counts", 20, -1.0, 1.0);
     fListHistCascade->Add(fHistRapOmega);
 }
 if(! fHistEtaXi ){
     fHistEtaXi = new TH1F("fHistEtaXi", "Pseudo-rap. of #Xi candidates (around the mass peak); #eta; Counts", 20, -1.0, 1.0);
     fListHistCascade->Add(fHistEtaXi);
 }
 if(! fHistEtaBachXi){
     fHistEtaBachXi = new TH1F("fHistEtaBachXi", "Pseudo-rap. of #Xi bachelor; #eta; Counts", 40, -2.0, 2.0);
     fListHistCascade->Add(fHistEtaBachXi);
 }
 if(! fHistEtaPosXi){
     fHistEtaPosXi = new TH1F("fHistEtaPosXi", "Pseudo-rap. of #Xi positive meson daughter; #eta; Counts", 40, -2.0, 2.0);
     fListHistCascade->Add(fHistEtaPosXi);
 }
 if(! fHistEtaNegXi){
     fHistEtaNegXi = new TH1F("fHistEtaNegXi", "Pseudo-rap. of #Xi negative meson daughter; #eta; Counts", 40, -2.0, 2.0);
     fListHistCascade->Add(fHistEtaNegXi);
 }
 if(! fHistThetaXi ){
     fHistThetaXi = new TH1F("fHistThetaXi", "#theta of #Xi candidates (around the mass peak); #theta (deg); Counts", 180, 0., 180.0);
     fListHistCascade->Add(fHistThetaXi);
 }
 if(! fHistPhiXi ){
     fHistPhiXi = new TH1F("fHistPhiXi", "#phi of #Xi candidates (around the mass peak); #phi (deg); Counts", 360, 0., 360.);
     fListHistCascade->Add(fHistPhiXi);
 }
 if(! f2dHistArmenteros) {
     f2dHistArmenteros = new TH2F("f2dHistArmenteros", "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm}; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
     fListHistCascade->Add(f2dHistArmenteros);
 }
    //TH2
 if(! f2dHistEffMassLambdaVsEffMassXiMinus) {
     f2dHistEffMassLambdaVsEffMassXiMinus = new TH2F("f2dHistEffMassLambdaVsEffMassXiMinus", "M_{#Lambda} Vs M_{#Xi^{-} candidates}; Inv. M_{#Lambda^{0}} (GeV/c^{2}); M( #Lambda , #pi^{-} ) (GeV/c^{2})", 300, 1.1, 1.13, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiMinus);
 }
 if(! f2dHistEffMassXiVsEffMassOmegaMinus) {
     f2dHistEffMassXiVsEffMassOmegaMinus = new TH2F("f2dHistEffMassXiVsEffMassOmegaMinus", "M_{#Xi^{-} candidates} Vs M_{#Omega^{-} candidates}; M( #Lambda , #pi^{-} ) (GeV/c^{2}); M( #Lambda , K^{-} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaMinus);
 }
 if(! f2dHistEffMassLambdaVsEffMassXiPlus) {
     f2dHistEffMassLambdaVsEffMassXiPlus = new TH2F("f2dHistEffMassLambdaVsEffMassXiPlus", "M_{#Lambda} Vs M_{#Xi^{+} candidates}; Inv. M_{#Lambda^{0}} (GeV/c^{2}); M( #Lambda , #pi^{+} ) (GeV/c^{2})", 300, 1.1, 1.13, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiPlus);
 }
 if(! f2dHistEffMassXiVsEffMassOmegaPlus) {
     f2dHistEffMassXiVsEffMassOmegaPlus = new TH2F("f2dHistEffMassXiVsEffMassOmegaPlus", "M_{#Xi^{+} candidates} Vs M_{#Omega^{+} candidates}; M( #Lambda , #pi^{+} ) (GeV/c^{2}); M( #Lambda , K^{+} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaPlus);
 }
 if(! f2dHistXiRadiusVsEffMassXiMinus) {
     f2dHistXiRadiusVsEffMassXiMinus = new TH2F("f2dHistXiRadiusVsEffMassXiMinus", "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2})", 450, 0., 45.0, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiMinus);
 }
 if(! f2dHistXiRadiusVsEffMassXiPlus) {
     f2dHistXiRadiusVsEffMassXiPlus = new TH2F("f2dHistXiRadiusVsEffMassXiPlus", "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2})", 450, 0., 45.0, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiPlus);
 }
 if(! f2dHistXiRadiusVsEffMassOmegaMinus) {
     f2dHistXiRadiusVsEffMassOmegaMinus = new TH2F("f2dHistXiRadiusVsEffMassOmegaMinus", "Transv. R_{Xi Decay} Vs M_{#Omega^{-} candidates}; r_{cascade} (cm); M( #Lambda , K^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaMinus);
 }
 if(! f2dHistXiRadiusVsEffMassOmegaPlus) {
     f2dHistXiRadiusVsEffMassOmegaPlus = new TH2F("f2dHistXiRadiusVsEffMassOmegaPlus", "Transv. R_{Xi Decay} Vs M_{#Omega^{+} candidates}; r_{cascade} (cm); M( #Lambda , K^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
     fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaPlus);
 }
 if(! f2dHistTPCdEdxOfCascDghters){
     f2dHistTPCdEdxOfCascDghters = new TH2F("f2dHistTPCdEdxOfCascDghters", "TPC dE/dx of the cascade daughters; charge x || #vec{p}_{TPC inner wall}(Casc. daughter) || (GeV/c); TPC signal (ADC)", 2000, -10.0, 10.0, 450, 0., 900.);
     fListHistCascade->Add(f2dHistTPCdEdxOfCascDghters);
 }
 if(! f2dHistDcaXiDaughtersvsInvMass){
     f2dHistDcaXiDaughtersvsInvMass = new TH2F("f2dHistDcaXiDaughtersvsInvMass", "DCA between Xi Daughters; DCA (cm); Number of Cascades", 100, 0., 0.5, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaXiDaughtersvsInvMass);
 }
 if(! f2dHistDcaBachToPrimVertexvsInvMass) {
     f2dHistDcaBachToPrimVertexvsInvMass = new TH2F("f2dHistDcaBachToPrimVertexvsInvMass", "DCA of Bach. to Prim. Vertex; DCA (cm); Number of Cascades", 250, 0., 0.25, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaBachToPrimVertexvsInvMass);
 }
 if(! f2dHistXiCosineOfPointingAnglevsInvMass){
     f2dHistXiCosineOfPointingAnglevsInvMass = new TH2F("f2dHistXiCosineOfPointingAnglevsInvMass", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl); Number of Xis", 200, 0.99, 1.0, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistXiCosineOfPointingAnglevsInvMass);
 }
 if(! f2dHistMassLambdaAsCascDghtervsInvMass){ 
     f2dHistMassLambdaAsCascDghtervsInvMass = new TH2F("f2dHistMassLambdaAsCascDghtervsInvMass","#Lambda associated to Casc. candidates; Eff. Mass (GeV/c^{2}); Counts", 300, 1.00, 1.3, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistMassLambdaAsCascDghtervsInvMass);
 }
 if(! f2dHistDcaV0DaughtersXivsInvMass){
     f2dHistDcaV0DaughtersXivsInvMass = new TH2F("f2dHistDcaV0DaughtersXivsInvMass", "DCA between V0 daughters, in cascade; DCA (cm); Number of V0s", 120, 0., 0.6, 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaV0DaughtersXivsInvMass);
 }
 if(! f2dHistDcaV0ToPrimVertexXivsInvMass){
     f2dHistDcaV0ToPrimVertexXivsInvMass = new TH2F("f2dHistDcaV0ToPrimVertexXivsInvMass", "DCA of V0 to Prim. Vertex, in cascade; DCA (cm); Number of Cascades", 200, 0., 1., 400, 1.2, 2.0);
     fListHistCascade->Add(f2dHistDcaV0ToPrimVertexXivsInvMass);
 }
 // - CFContainer PID study Xi minus
 if(!fCFContCascadePIDXiMinus)  {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 800;
   lNbBinsPerVar[2] = 22;
   if (fkSDDSelectionOn) {
        if (fwithSDD) fCFContCascadePIDXiMinus = new AliCFContainer(Form("fCFContCascadePIDXiMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDon",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
        else if (!fwithSDD) fCFContCascadePIDXiMinus = new AliCFContainer(Form("fCFContCascadePIDXiMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDoff",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar ); 
   } else if (!fkSDDSelectionOn) fCFContCascadePIDXiMinus = new AliCFContainer(Form("fCFContCascadePIDXiMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_woSDD",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
   fCFContCascadePIDXiMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
     //Setting the step title : one per PID case
   fCFContCascadePIDXiMinus->SetStepTitle(0, "No PID");
   fCFContCascadePIDXiMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDXiMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDXiMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");  
     //Setting the variable title, per axis
   fCFContCascadePIDXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
   fCFContCascadePIDXiMinus->SetVarTitle(2, "Y_{#Xi}");
   fListHistCascade->Add(fCFContCascadePIDXiMinus);
 }
 // - CFContainer PID study Xi plus
 if (!fCFContCascadePIDXiPlus) {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 800;
   lNbBinsPerVar[2] = 22;
   if (fkSDDSelectionOn) {
        if (fwithSDD) fCFContCascadePIDXiPlus = new AliCFContainer(Form("fCFContCascadePIDXiPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDon",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
        else if (!fwithSDD) fCFContCascadePIDXiPlus = new AliCFContainer(Form("fCFContCascadePIDXiPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDoff",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
   } else if (!fkSDDSelectionOn) fCFContCascadePIDXiPlus = new AliCFContainer(Form("fCFContCascadePIDXiPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_woSDD",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
   fCFContCascadePIDXiPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
     //Setting the step title : one per PID case
   fCFContCascadePIDXiPlus->SetStepTitle(0, "No PID");
   fCFContCascadePIDXiPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDXiPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDXiPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
     //Setting the variable title, per axis
   fCFContCascadePIDXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
   fCFContCascadePIDXiPlus->SetVarTitle(2, "Y_{#Xi}");
   fListHistCascade->Add(fCFContCascadePIDXiPlus);
 }
 // - CFContainer PID study Omega minus
 if(!fCFContCascadePIDOmegaMinus)  {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 1000;
   lNbBinsPerVar[2] = 22;
   if (fkSDDSelectionOn) {
        if (fwithSDD) fCFContCascadePIDOmegaMinus = new AliCFContainer(Form("fCFContCascadePIDOmegaMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDon",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
        else if (!fwithSDD) fCFContCascadePIDOmegaMinus = new AliCFContainer(Form("fCFContCascadePIDOmegaMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDoff",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
   } else if (!fkSDDSelectionOn) fCFContCascadePIDOmegaMinus = new AliCFContainer(Form("fCFContCascadePIDOmegaMinus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_woSDD",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
   fCFContCascadePIDOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
     //Setting the step title : one per PID case
   fCFContCascadePIDOmegaMinus->SetStepTitle(0, "No PID");
   fCFContCascadePIDOmegaMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDOmegaMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDOmegaMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
     //Setting the variable title, per axis
   fCFContCascadePIDOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
   fCFContCascadePIDOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
   fListHistCascade->Add(fCFContCascadePIDOmegaMinus);
 }
 // - CFContainer PID study Omega plus
 if(!fCFContCascadePIDOmegaPlus)  {
   const Int_t  lNbSteps      =  7 ;
   const Int_t  lNbVariables  =  3 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[3] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 1000;
   lNbBinsPerVar[2] = 22; 
   if (fkSDDSelectionOn) {
        if (fwithSDD) fCFContCascadePIDOmegaPlus = new AliCFContainer(Form("fCFContCascadePIDOmegaPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDon",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
        else if (!fwithSDD) fCFContCascadePIDOmegaPlus = new AliCFContainer(Form("fCFContCascadePIDOmegaPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDoff",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
   } else if (!fkSDDSelectionOn) fCFContCascadePIDOmegaPlus = new AliCFContainer(Form("fCFContCascadePIDOmegaPlus_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_woSDD",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
     //Setting the bin limits 
   fCFContCascadePIDOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
   fCFContCascadePIDOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
   fCFContCascadePIDOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity 
     //Setting the step title : one per PID case
   fCFContCascadePIDOmegaPlus->SetStepTitle(0, "No PID");
   fCFContCascadePIDOmegaPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
   fCFContCascadePIDOmegaPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
   fCFContCascadePIDOmegaPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
   fCFContCascadePIDOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
   fCFContCascadePIDOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
   fCFContCascadePIDOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
     //Setting the variable title, per axis
   fCFContCascadePIDOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
   fCFContCascadePIDOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
   fCFContCascadePIDOmegaPlus->SetVarTitle(2, "Y_{#Omega}");  
   fListHistCascade->Add(fCFContCascadePIDOmegaPlus);
 }
 // - CFContainer: towards the optimisation of topological selections
 if(! fCFContCascadeCuts) {
	// Container meant to store all the relevant distributions corresponding to the cut variables.
        // NB: overflow/underflow of variables on which we want to cut later should be 0!!! 
   const Int_t  lNbSteps      =  4 ;
   const Int_t  lNbVariables  =  19 ;
     //Array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[lNbVariables] = {0};
   lNbBinsPerVar[0]  = 25;     //DcaCascDaughters             :  [0.0,2.4,3.0]       -> Rec.Cut = 2.0;
   lNbBinsPerVar[1]  = 25;     //DcaBachToPrimVertex          :  [0.0,0.24,100.0]    -> Rec.Cut = 0.01; 
   lNbBinsPerVar[2]  = 30;     //CascCosineOfPointingAngle    :  [0.97,1.0]          -> Rec.Cut = 0.98;
   lNbBinsPerVar[3]  = 40;     //CascRadius                   :  [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
   lNbBinsPerVar[4]  = 30;     //InvMassLambdaAsCascDghter    :  [1.1,1.3]           -> Rec.Cut = 0.008;
   lNbBinsPerVar[5]  = 20;     //DcaV0Daughters               :  [0.0,2.0]           -> Rec.Cut = 1.5;
   lNbBinsPerVar[6]  = 201;    //V0CosineOfPointingAngle      :  [0.89,1.0]          -> Rec.Cut = 0.9;
   lNbBinsPerVar[7]  = 40;     //V0Radius                     :  [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
   lNbBinsPerVar[8]  = 40;     //DcaV0ToPrimVertex            :  [0.0,0.39,110.0]    -> Rec.Cut = 0.01;  
   lNbBinsPerVar[9]  = 25;     //DcaPosToPrimVertex           :  [0.0,0.24,100.0]    -> Rec.Cut = 0.05;
   lNbBinsPerVar[10] = 25;     //DcaNegToPrimVertex           :  [0.0,0.24,100.0]    -> Rec.Cut = 0.05
   lNbBinsPerVar[11] = 150;    //InvMassXi                    :   2-MeV/c2 bins
   lNbBinsPerVar[12] = 120;    //InvMassOmega                 :   2-MeV/c2 bins
   lNbBinsPerVar[13] = 100;    //XiTransvMom                  :  [0.0,10.0]
   lNbBinsPerVar[14] = 110;    //Y(Xi)                        :   0.02 in rapidity units
   lNbBinsPerVar[15] = 110;    //Y(Omega)                     :   0.02 in rapidity units
   lNbBinsPerVar[16] = 112;    //Proper lenght of cascade       
   lNbBinsPerVar[17] = 112;    //Proper lenght of V0
   lNbBinsPerVar[18] = 112;    //Distance V0-Xi in transverse plane
   if (fkSDDSelectionOn) {
        if (fwithSDD) fCFContCascadeCuts = new AliCFContainer(Form("fCFContCascadeCuts_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDon",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar);
        else if (!fwithSDD) fCFContCascadeCuts = new AliCFContainer(Form("fCFContCascadeCuts_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_wSDDoff",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar);
   } else if (!fkSDDSelectionOn) fCFContCascadeCuts = new AliCFContainer(Form("fCFContCascadeCuts_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f_woSDD",fMinnTPCcls,fVtxRange,fVtxRangeMin,fMinPtCutOnDaughterTracks,fEtaCutOnDaughterTracks),"Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar);
     //Setting the bin limits 
     //0 -  DcaXiDaughters
   Double_t *lBinLim0  = new Double_t[ lNbBinsPerVar[0] + 1 ];
        for(Int_t i=0; i< lNbBinsPerVar[0]; i++) lBinLim0[i] = (Double_t)0.0 + (2.4 - 0.0)/(lNbBinsPerVar[0] - 1) * (Double_t)i;
        lBinLim0[ lNbBinsPerVar[0] ] = 3.0;
   fCFContCascadeCuts -> SetBinLimits(0, lBinLim0);  
   delete [] lBinLim0;
     //1 - DcaToPrimVertexXi
   Double_t *lBinLim1  = new Double_t[ lNbBinsPerVar[1] + 1 ];
        for(Int_t i=0; i<lNbBinsPerVar[1]; i++) lBinLim1[i] = (Double_t)0.0 + (0.24  - 0.0)/(lNbBinsPerVar[1] - 1) * (Double_t)i;
        lBinLim1[ lNbBinsPerVar[1] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(1, lBinLim1);  
   delete [] lBinLim1;    
     //2 - CascCosineOfPointingAngle 
   fCFContCascadeCuts->SetBinLimits(2, 0.97, 1.);
     //3 - CascRadius
   Double_t *lBinLim3  = new Double_t[ lNbBinsPerVar[3]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[3]; i++)   lBinLim3[i]  = (Double_t)0.0   + (3.9  - 0.0 )/(lNbBinsPerVar[3] - 1)  * (Double_t)i ;
        lBinLim3[ lNbBinsPerVar[3] ] = 1000.0;
   fCFContCascadeCuts -> SetBinLimits(3,  lBinLim3 );        
   delete [] lBinLim3;
     //4 - InvMassLambdaAsCascDghter
   fCFContCascadeCuts->SetBinLimits(4, 1.1, 1.13);
     //5 - DcaV0Daughters
   fCFContCascadeCuts -> SetBinLimits(5, 0., 2.);
     //6 - V0CosineOfPointingAngle
   fCFContCascadeCuts -> SetBinLimits(6, 0.8, 1.001);
     //7 - V0Radius
   Double_t *lBinLim7 = new Double_t[ lNbBinsPerVar[7] + 1];
        for(Int_t i=0; i< lNbBinsPerVar[7];i++) lBinLim7[i] = (Double_t)0.0 + (3.9 - 0.0)/(lNbBinsPerVar[7] - 1) * (Double_t)i;
        lBinLim7[ lNbBinsPerVar[7] ] = 1000.0;
   fCFContCascadeCuts -> SetBinLimits(7, lBinLim7); 
   delete [] lBinLim7;
     //8 - DcaV0ToPrimVertex
   Double_t *lBinLim8  = new Double_t[ lNbBinsPerVar[8]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[8];i++)   lBinLim8[i]  = (Double_t)0.0   + (0.39  - 0.0 )/(lNbBinsPerVar[8]-1)  * (Double_t)i ;
        lBinLim8[ lNbBinsPerVar[8]  ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(8,  lBinLim8 );    
   delete [] lBinLim8;
     //9 - DcaPosToPrimVertex
   Double_t *lBinLim9  = new Double_t[ lNbBinsPerVar[9]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[9];i++)   lBinLim9[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[9]-1)  * (Double_t)i ;
        lBinLim9[ lNbBinsPerVar[9]  ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(9,  lBinLim9 );        
   delete [] lBinLim9;
     //10 - DcaNegToPrimVertex
   Double_t *lBinLim10  = new Double_t[ lNbBinsPerVar[10]+1 ];
        for(Int_t i=0; i< lNbBinsPerVar[10];i++)   lBinLim10[i]  = (Double_t)0.0   + (0.24  - 0.0 )/(lNbBinsPerVar[10]-1)  * (Double_t)i ;
        lBinLim10[ lNbBinsPerVar[10]  ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(10,  lBinLim10 );            // DcaPosToPrimVertexXi : 0.0 to 0.25 
   delete [] lBinLim10;
     //11 - InvMassXi
   fCFContCascadeCuts->SetBinLimits(11, 1.25, 1.40);
     //12 - InvMassOmega
   fCFContCascadeCuts->SetBinLimits(12, 1.62, 1.74);
     //13 - XiTransvMom
   fCFContCascadeCuts->SetBinLimits(13, 0.0, 10.0); 
     //14 - Y(Xi)
   fCFContCascadeCuts->SetBinLimits(14, -1.1, 1.1);
     //15 - Y(Omega)
   fCFContCascadeCuts->SetBinLimits(15, -1.1, 1.1);
     //16 - Proper time of cascade
   Double_t *lBinLim16  = new Double_t[ lNbBinsPerVar[16]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[16];i++) lBinLim16[i] = (Double_t) -1. + (110. + 1.0 ) / (lNbBinsPerVar[16] - 1) * (Double_t) i;
   lBinLim16[ lNbBinsPerVar[16] ] = 2000.0;
   fCFContCascadeCuts->SetBinLimits(16, lBinLim16);
     //17 - Proper time of V0
   fCFContCascadeCuts->SetBinLimits(17, lBinLim16);
     //18 - Distance V0-Xi in transverse plane
   fCFContCascadeCuts->SetBinLimits(18, lBinLim16);
     // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
   fCFContCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates");
   fCFContCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
   fCFContCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates");
   fCFContCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
     // Setting the variable title, per axis
   fCFContCascadeCuts->SetVarTitle(0,  "Dca(cascade daughters) (cm)");
   fCFContCascadeCuts->SetVarTitle(1,  "ImpactParamToPV(bachelor) (cm)");
   fCFContCascadeCuts->SetVarTitle(2,  "cos(cascade PA)");
   fCFContCascadeCuts->SetVarTitle(3,  "R_{2d}(cascade decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(as casc dghter) (GeV/c^{2})");
   fCFContCascadeCuts->SetVarTitle(5,  "Dca(V0 daughters) in Xi (cm)");
   fCFContCascadeCuts->SetVarTitle(6,  "cos(V0 PA) in cascade");
   fCFContCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(8,  "ImpactParamToPV(V0) (cm)");
   fCFContCascadeCuts->SetVarTitle(9,  "ImpactParamToPV(Pos) (cm)");
   fCFContCascadeCuts->SetVarTitle(10, "ImpactParamToPV(Neg) (cm)");
   fCFContCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
   fCFContCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
   fCFContCascadeCuts->SetVarTitle(13, "pt(cascade) (GeV/c)");
   fCFContCascadeCuts->SetVarTitle(14, "Y(Xi)");
   fCFContCascadeCuts->SetVarTitle(15, "Y(Omega)");
   fCFContCascadeCuts->SetVarTitle(16, "mL/p (cascade) (cm)");
   fCFContCascadeCuts->SetVarTitle(17, "mL/p (V0) (cm)");
   fCFContCascadeCuts->SetVarTitle(18, "Distance V0-Cascade in transverse plane (cm)");
   fListHistCascade->Add(fCFContCascadeCuts);
 }

 PostData(1, fListHistCascade);
 PostData(2, fCFContCascadePIDXiMinus);
 PostData(3, fCFContCascadePIDXiPlus);
 PostData(4, fCFContCascadePIDOmegaMinus);
 PostData(5, fCFContCascadePIDOmegaPlus);
 PostData(6, fCFContCascadeCuts);
} // end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskCheckCascadepp276::UserExec(Option_t *) {

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Main loop (called for each event)
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  //----------------
  //Define variables 
  AliESDEvent *lESDevent = 0x0;
  AliAODEvent *lAODevent = 0x0;

  //---------------------
  //Check the PIDresponse
  if(!fPIDResponse) {
       AliError("Cannot get pid response");
       return;
  }

  ///////////////////
  // EVENT SELECTIONS
  ///////////////////
  // In order:
  // 1) SDD selection
  // 2) Physics selection
  // 3) Select only looking at events with well-established PV
  // 4) Pileup selection
  // 5) |Z| < 10 cm

  //----------------------
  // Before any selections
  //----------------------
  //- Define the variables
  Int_t ncascadesBeforeAnySel = 0;
  Int_t nTrackMultiplicityBeforeAnySel = 0;
  if (fAnalysisType == "ESD") {
      // - Load the InputEvent and check
      lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
      if (!lESDevent) {
         AliWarning("ERROR: lESDevent not available \n");
         return;
      }
      // - Take the number of cascades and tracks before any events selection
      ncascadesBeforeAnySel = lESDevent->GetNumberOfCascades();
      nTrackMultiplicityBeforeAnySel = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
  } else if (fAnalysisType == "AOD") {
      // - Load the InputEvent and check
      lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
      if (!lAODevent) {
          AliWarning("ERROR: lAODevent not available \n");
          return;
      }
      // - Take the number of cascades and tracks before any events selection
      ncascadesBeforeAnySel  = lAODevent->GetNumberOfCascades();
      nTrackMultiplicityBeforeAnySel = -100;  //FIXME: I can't find the equivalent method for the AOD  
  } else {
      Printf("Analysis type (ESD or AOD) not specified \n");
      return;
  }
  // - Fill the plots
  fHistCascadeMultiplicityBeforeAnySel->Fill(ncascadesBeforeAnySel);
  fHistTrackMultiplicityBeforeAnySel->Fill(nTrackMultiplicityBeforeAnySel);
 
  //--------------
  // SDD selection
  //--------------
  // - Define the variables
  Int_t ncascadesAfterSDDSel = 0;
  Int_t nTrackMultiplicityAfterSDDSel = 0;
  // - Selection for ESD and AOD
  if (fAnalysisType == "ESD") {
     if (fkSDDSelectionOn) {
        TString trcl = lESDevent->GetFiredTriggerClasses();
        //cout<<"Fired Trigger Classes: "<<trcl<<endl;
        if (fwithSDD){
          if(!(trcl.Contains("ALLNOTRD"))) {
               cout<<"We are selecting events with SDD turn ON. This event has the SDD turn OFF. =>  RETURN!! (Exclude it)..."<<endl;
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
          }
        } else if (!fwithSDD){
          if((trcl.Contains("ALLNOTRD"))) {
               cout<<"We are selecting events with SDD turn OFF. This event has the SDD turn ON. =>  RETURN!! (Exclude it)..."<<endl;
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
          }
        }
     }
     // - Take the number of cascades and tracks after the SDD selection
     ncascadesAfterSDDSel = lESDevent->GetNumberOfCascades();
     nTrackMultiplicityAfterSDDSel = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
  } else if (fAnalysisType == "AOD") {
     if (fkSDDSelectionOn) {
        TString trcl = lAODevent->GetFiredTriggerClasses();
        if (fwithSDD){
           if(!(trcl.Contains("ALLNOTRD"))) {
                PostData(1, fListHistCascade);
                PostData(2, fCFContCascadePIDXiMinus);
                PostData(3, fCFContCascadePIDXiPlus);
                PostData(4, fCFContCascadePIDOmegaMinus);
                PostData(5, fCFContCascadePIDOmegaPlus);
                PostData(6, fCFContCascadeCuts);
                cout<<"We are selecting events with SDD turn ON. This event has the SDD turn OFF. =>  RETURN!! (Exclude it)..."<<endl;
                return;
           }
        } else if (!fwithSDD) {
           if((trcl.Contains("ALLNOTRD"))) {
                PostData(1, fListHistCascade);
                PostData(2, fCFContCascadePIDXiMinus);
                PostData(3, fCFContCascadePIDXiPlus);
                PostData(4, fCFContCascadePIDOmegaMinus);
                PostData(5, fCFContCascadePIDOmegaPlus);
                PostData(6, fCFContCascadeCuts);
                cout<<"We are selecting events with SDD turn OFF. This event has the SDD turn ON. =>  RETURN!! (Exclude it)..."<<endl;
                return;
           }
        }
     }
     // - Take the number of cascades and tracks after the SDD selection
     ncascadesAfterSDDSel = lAODevent->GetNumberOfCascades();
     nTrackMultiplicityAfterSDDSel = -100; //FIXME: I can't find the equivalent method for the AOD
  }
  // - Fill the plots
  fHistCascadeMultiplicityAfterSDDSel->Fill(ncascadesAfterSDDSel);
  fHistTrackMultiplicityAfterSDDSel->Fill(nTrackMultiplicityAfterSDDSel);

  //----------------------------------------------
  // Physics selection (+ re-vertexer for the ESD)
  //----------------------------------------------
  // - Define the variables
  Int_t ncascadesAfterPhysicsSel = 0;
  Int_t nTrackMultiplicityAfterPhysicsSel = 0;
  // - Selection for ESD and AOD
  if (fAnalysisType == "ESD") {
      UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
      Bool_t isSelected = 0;
      isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
      if(! isSelected){
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          cout<<"We are selecting the events that past tha Physics Selection. This event does not pass the Physics Selection. =>  RETURN!! (Exclude it)..."<<endl;
          return;
      }
      // - Take the number of cascades and tracks after physics selection
      ncascadesAfterPhysicsSel = lESDevent->GetNumberOfCascades();    
      nTrackMultiplicityAfterPhysicsSel = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);  
      // - Cascade vertexer (ESD)
      // Relaunch V0 and Cascade vertexers
      if (fkRerunV0CascVertexers) { 
            lESDevent->ResetCascades();
            lESDevent->ResetV0s();
            //AliV0vertexer *lV0vtxer = new AliV0vertexer();
            //AliCascadeVertexer *lCascVtxer = new AliCascadeVertexer();
            //lV0vtxer->GetCuts(fV0Sels);
            //lCascVtxer->GetCuts(fCascSels);
            //lV0vtxer->SetCuts(fV0Sels);      // NB don't use SetDefaultCuts!! because it acts on static variables 
            //lCascVtxer->SetCuts(fCascSels);
            //lV0vtxer->Tracks2V0vertices(lESDevent);
            //lCascVtxer->V0sTracks2CascadeVertices(lESDevent);
            //delete lV0vtxer;
            //delete lCascVtxer;
      }         
  } else if (fAnalysisType == "AOD") {
      UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
      Bool_t isSelected = 0;
      isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
      if(! isSelected){
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          cout<<"We are selecting the events that past tha Physics Selection. This event does not pass the Physics Selection. =>  RETURN!! (Exclude it)..."<<endl;
          return;
      }    
      // - Take the number of cascades and tracks after the physics selection
      ncascadesAfterPhysicsSel = lAODevent->GetNumberOfCascades();
      nTrackMultiplicityAfterPhysicsSel = -100;  //FIXME: I can't find the equivalent method for the AOD    
  } 
  // - Fill the plots
  fHistCascadeMultiplicityAfterPhysicsSel->Fill(ncascadesAfterPhysicsSel);
  fHistTrackMultiplicityAfterPhysicsSel->Fill(nTrackMultiplicityAfterPhysicsSel);

  //------------------------------
  // Well-established PV selection
  //------------------------------
  // - Define variables
  Int_t ncascadesForSelEvtNoTPCOnly = 0;
  Int_t nTrackMultiplicityForSelEvtNoTPCOnly = 0;
  // - Selection for ESD and AOD
  if (fAnalysisType == "ESD") {
      // - Vertex coordinates: get the PVs stored in the ESD found with tracks and SPD
      const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
      const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
      // - Select only looking at events with well-established PV
      if (fkQualityCutNoTPConlyPrimVtx) {
          if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus() ){
              AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
          }
      }
      // - Take the number of cascades and tracks after TPConly selection
      ncascadesForSelEvtNoTPCOnly = lESDevent->GetNumberOfCascades();
      nTrackMultiplicityForSelEvtNoTPCOnly = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
  } else if (fAnalysisType == "AOD") {
      // - Vertex coordinates: get the PVs stored in the AOD found with tracks and SPD
      const AliAODVertex *lPrimarySPDVtx = lAODevent->GetPrimaryVertexSPD();
      const AliAODVertex *lPrimaryTrackingAODVtx = lAODevent->GetPrimaryVertex();
      // - Select only looking at events with well-established PV
      if (fkQualityCutNoTPConlyPrimVtx) {
          if (!lPrimarySPDVtx && !lPrimaryTrackingAODVtx) {
              AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
          }
      }
      // - Take the number of cascades and tracks after TPConly selection
      ncascadesForSelEvtNoTPCOnly = lAODevent->GetNumberOfCascades();
      nTrackMultiplicityForSelEvtNoTPCOnly = -100;  //FIXME: I can't find the equivalent method for the AOD
  }
  // - Fill the plots
  fHistCascadeMultiplicityForSelEvtNoTPCOnly->Fill(ncascadesForSelEvtNoTPCOnly);
  fHistTrackMultiplicityForSelEvtNoTPCOnly->Fill(nTrackMultiplicityForSelEvtNoTPCOnly);
    
  //----------------
  // Pilup selection
  //----------------
  // - Define variables
  Int_t ncascadesForSelEvtNoTPCOnlyNoPileup = 0;
  Int_t nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = 0;
  // - Selection for ESD and AOD
  if (fAnalysisType == "ESD") {
      if (fkQualityCutPileup) {
          if(lESDevent->IsPileupFromSPD()){
              AliWarning("Pb / Pile-up event ... return!");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
          }
      }
      // - Take the number of cascades and tracks after Pileup selection
      ncascadesForSelEvtNoTPCOnlyNoPileup = lESDevent->GetNumberOfCascades();
      nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5);
  } else if (fAnalysisType == "AOD") {
      if (fkQualityCutPileup) {
          if(lAODevent->IsPileupFromSPD()){
              AliWarning("Pb / Pile-up event ... return!");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
          }
      }
      // - Take the number of cascades and tracks after Pileup selection
      ncascadesForSelEvtNoTPCOnlyNoPileup = lAODevent->GetNumberOfCascades();
      nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup = -100;   //FIXME: I can't find the equivalent method for the AOD
  }
  fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup->Fill(ncascadesForSelEvtNoTPCOnlyNoPileup);
  fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup->Fill(nTrackMultiplicityForSelEvtNoTPCOnlyNoPileup);
    
  //----------------------------------------------------
  // Vertex Z position selection (+ magnetic field info)
  //----------------------------------------------------
  // - Define variables
  Double_t lBestPrimaryVtxPos[3]  = {-100.0, -100.0, -100.0};
  Double_t lMagneticField         = -10.;
  Double_t tPrimaryVtxPosition[3] = {-100.0, -100.0, -100.0};
  Int_t ncascadesAfterVertexSel          = 0;
  Int_t nTrackMultiplicityAfterVertexSel = 0; 
  // - Selection for ESD and AOD
  if (fAnalysisType == "ESD") {
      // - Vertex coordinates: get the best primary vertex available for the event 
      const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex(); 
      if (!lPrimaryBestESDVtx){
          AliWarning("No prim. vertex in ESD... return!");
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          return;
      }
      lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
      // - Fill the vertex plots before any event selection on vertex position
      const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
      tPrimaryVtxPosition[0] = primaryVtx->GetX();
      tPrimaryVtxPosition[1] = primaryVtx->GetY();
      tPrimaryVtxPosition[2] = primaryVtx->GetZ();
	  fHistPVx->Fill( tPrimaryVtxPosition[0] );
	  fHistPVy->Fill( tPrimaryVtxPosition[1] );
	  fHistPVz->Fill( tPrimaryVtxPosition[2] );	  
      // - Get magnetic filed info
      lMagneticField = lESDevent->GetMagneticField();
      //if(TMath::Abs(lMagneticField ) < 10e-6) continue;
      // - Selection on the primary vertex Z position 
      if (fkQualityCutZprimVtxPos) {
          if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange || TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin) {
               AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
               PostData(1, fListHistCascade);
               PostData(2, fCFContCascadePIDXiMinus);
               PostData(3, fCFContCascadePIDXiPlus);
               PostData(4, fCFContCascadePIDOmegaMinus);
               PostData(5, fCFContCascadePIDOmegaPlus);
               PostData(6, fCFContCascadeCuts);
               return;
          }
      }
      // - Take the number of cascades and tracks after vertex Z position selection
      ncascadesAfterVertexSel = lESDevent->GetNumberOfCascades();
      nTrackMultiplicityAfterVertexSel = fESDtrackCuts->GetReferenceMultiplicity(lESDevent,AliESDtrackCuts::kTrackletsITSTPC,0.5); 
  } else if (fAnalysisType == "AOD") {
      // - Vertex coordinates: get the best primary vertex available for the event
      const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
      if (!lPrimaryBestAODVtx){
          AliWarning("No prim. vertex in AOD... return!");
          PostData(1, fListHistCascade);
          PostData(2, fCFContCascadePIDXiMinus);
          PostData(3, fCFContCascadePIDXiPlus);
          PostData(4, fCFContCascadePIDOmegaMinus);
          PostData(5, fCFContCascadePIDOmegaPlus);
          PostData(6, fCFContCascadeCuts);
          return;
      }
      lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
      // - Fill the vertex plots before any event selection on vertex position
      const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
      tPrimaryVtxPosition[0] = primaryVtx->GetX();
      tPrimaryVtxPosition[1] = primaryVtx->GetY();
      tPrimaryVtxPosition[2] = primaryVtx->GetZ();
	  fHistPVx->Fill( tPrimaryVtxPosition[0] );
	  fHistPVy->Fill( tPrimaryVtxPosition[1] );
	  fHistPVz->Fill( tPrimaryVtxPosition[2] );
      // - Get magnetic filed info
      lMagneticField = lAODevent->GetMagneticField();
      //if(TMath::Abs(lMagneticField ) < 10e-6) continue;
      // - Selection on the primary vertex Z position 
      if (fkQualityCutZprimVtxPos) {
          if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange || TMath::Abs(lBestPrimaryVtxPos[2]) < fVtxRangeMin) {
              AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
              PostData(1, fListHistCascade);
              PostData(2, fCFContCascadePIDXiMinus);
              PostData(3, fCFContCascadePIDXiPlus);
              PostData(4, fCFContCascadePIDOmegaMinus);
              PostData(5, fCFContCascadePIDOmegaPlus);
              PostData(6, fCFContCascadeCuts);
              return;
          }
      }
      // - Take the number of cascades and tracks after vertex Z position selection
      ncascadesAfterVertexSel = lAODevent->GetNumberOfCascades();
      nTrackMultiplicityAfterVertexSel = -100; //FIXME: I can't find the equivalent method for the AOD
  } 
  // - Fill the plots
  fHistCascadeMultiplicityAfterVertexCutSel->Fill(ncascadesAfterVertexSel);
  fHistTrackMultiplicityAfterVertexCutSel->Fill(nTrackMultiplicityAfterVertexSel);

  // - Vertex position plots: after any event selections
  tPrimaryVtxPosition[0] = -100.0;
  tPrimaryVtxPosition[1] = -100.0;
  tPrimaryVtxPosition[2] = -100.0;
  if (fAnalysisType == "ESD") {
       const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
       tPrimaryVtxPosition[0] = primaryVtx->GetX();
       tPrimaryVtxPosition[1] = primaryVtx->GetY();
       tPrimaryVtxPosition[2] = primaryVtx->GetZ();
  } else if (fAnalysisType == "AOD") {
       const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
       tPrimaryVtxPosition[0] = primaryVtx->GetX();
       tPrimaryVtxPosition[1] = primaryVtx->GetY();
       tPrimaryVtxPosition[2] = primaryVtx->GetZ();
  }
  fHistPVxAnalysis->Fill( tPrimaryVtxPosition[0] );
  fHistPVyAnalysis->Fill( tPrimaryVtxPosition[1] );
  fHistPVzAnalysis->Fill( tPrimaryVtxPosition[2] );
  


  //////////////////////////////
  // CASCADE RECONSTRUCTION PART
  //////////////////////////////
  
  //%%%%%%%%%%%%%
  // Cascade loop
  Int_t ncascades = 0;
  if      (fAnalysisType == "ESD") ncascades = lESDevent->GetNumberOfCascades();
  else if (fAnalysisType == "AOD") ncascades = lAODevent->GetNumberOfCascades();

  for (Int_t iXi = 0; iXi < ncascades; iXi++) {// This is the begining of the Cascade loop (ESD or AOD)
	   
      // -----------------------------------------------------------------------
      // - Initialisation of the local variables that will be needed for ESD/AOD

      // - 0th part of initialisation : around primary vertex ...
      //Double_t lBestPrimaryVtxRadius3D     = -500.0;
      // - 1st part of initialisation : variables needed to store AliESDCascade data members
      Double_t lEffMassXi                  = 0.;
      Double_t lDcaXiDaughters             = -1.;
      Double_t lXiCosineOfPointingAngle    = -1.;
      Double_t lPosXi[3]                   = { -1000.0, -1000.0, -1000.0 };
      Double_t lXiRadius                   = -1000. ;
      // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
      UShort_t lPosTPCClusters             = -1; // For ESD only ...
      UShort_t lNegTPCClusters             = -1; // For ESD only ...
      UShort_t lBachTPCClusters            = -1; // For ESD only ...
      Double_t lInnerWallMomCascDghters[3] = {-100., -100., -100.};
      Double_t lTPCSignalCascDghters   [3] = {-100., -100., -100.};
      // - 3rd part of initialisation : about V0 part in cascades
      Double_t lInvMassLambdaAsCascDghter  = 0.;
      Double_t lDcaV0DaughtersXi           = -1.;
      Double_t lDcaBachToPrimVertexXi      = -1.;
      Double_t lDcaV0ToPrimVertexXi        = -1.;
      Double_t lDcaPosToPrimVertexXi       = -1.;
      Double_t lDcaNegToPrimVertexXi       = -1.;
      Double_t lV0CosineOfPointingAngleXi  = -1. ;
      Double_t lPosV0Xi[3]                 = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
      Double_t lV0RadiusXi                 = -1000.0;
      Double_t lV0quality                  = 0.;
      // - 4th part of initialisation : Effective masses
      Double_t lInvMassXiMinus             = 0.;
      Double_t lInvMassXiPlus              = 0.;
      Double_t lInvMassOmegaMinus          = 0.;
      Double_t lInvMassOmegaPlus           = 0.;
      // - 5th part of initialisation : PID treatment
      Bool_t   lIsPosInXiProton            = kFALSE;
      Bool_t   lIsPosInXiPion              = kFALSE;
      Bool_t   lIsPosInOmegaProton         = kFALSE;
      Bool_t   lIsPosInOmegaPion           = kFALSE;
      Bool_t   lIsNegInXiProton            = kFALSE;
      Bool_t   lIsNegInXiPion              = kFALSE;
      Bool_t   lIsNegInOmegaProton         = kFALSE;
      Bool_t   lIsNegInOmegaPion           = kFALSE;
      Bool_t   lIsBachelorKaon             = kFALSE;
      Bool_t   lIsBachelorPion             = kFALSE; 
      Bool_t   lIsBachelorKaonForTPC       = kFALSE; 
      Bool_t   lIsBachelorPionForTPC       = kFALSE; 
      Bool_t   lIsNegPionForTPC            = kFALSE; 
      Bool_t   lIsPosPionForTPC            = kFALSE; 
      Bool_t   lIsNegProtonForTPC          = kFALSE; 
      Bool_t   lIsPosProtonForTPC          = kFALSE; 
      // - 6th part of initialisation : extra info for QA
      Double_t lXiMomX         = 0.;
      Double_t lXiMomY         = 0.;
      Double_t lXiMomZ         = 0.;
      Double_t lXiTransvMom    = 0.;
      Double_t lXiTotMom       = 0.;
      Double_t lV0PMomX        = 0.;
      Double_t lV0PMomY        = 0.;
      Double_t lV0PMomZ        = 0.;
      Double_t lV0NMomX        = 0.;
      Double_t lV0NMomY        = 0.;
      Double_t lV0NMomZ        = 0.;
      Double_t lV0TotMom       = 0.;
      Double_t lBachMomX       = 0.;
      Double_t lBachMomY       = 0.;
      Double_t lBachMomZ       = 0.;
      Double_t lBachTransvMom  = 0.;
      Double_t lBachTotMom     = 0.;
      Double_t lpTrackTransvMom  = 0.;
      Double_t lnTrackTransvMom  = 0.;
      Short_t  lChargeXi       = -2;
      Double_t lV0toXiCosineOfPointingAngle = 0.;
      Double_t lRapXi   = -20.0, lRapOmega = -20.0, lEta = -20.0, lTheta = 360., lPhi = 720.;
      Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
      Double_t etaPos = 0., etaNeg = 0., etaBach = 0.;
  
      if (fAnalysisType == "ESD") { 
  
          // -------------------------------------
          // - Load the cascades from the handler
          AliESDcascade *xi = lESDevent->GetCascade(iXi);
          if (!xi) continue;

          //----------------------------------------------------------------------------        
          // - Assigning the necessary variables for specific AliESDcascade data members	
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis: cascade = Xi-decay
          lEffMassXi  	           = xi->GetEffMassXi();
          lDcaXiDaughters	   = xi->GetDcaXiDaughters();
          lXiCosineOfPointingAngle = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
	                                                               //Take care : the best available vertex should be used (like in AliCascadeVertexer)
          xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
          lXiRadius        	   = TMath::Sqrt( lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1] );
		
          //-------------------------------------------------------------------------------------------------------------------------------
	  // - Around the tracks: Bach + V0 (ESD). Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
          UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
          UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
          UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                                      //Care track label can be negative in MC production (linked with the track quality)
                                      //However = normally, not the case for track index ...
          // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
          if (lBachIdx == lIdxNegXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue; }
          if (lBachIdx == lIdxPosXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue; }
          // - Get the track for the daughters
          AliESDtrack *pTrackXi	   = lESDevent->GetTrack( lIdxPosXi );
          AliESDtrack *nTrackXi	   = lESDevent->GetTrack( lIdxNegXi );
          AliESDtrack *bachTrackXi = lESDevent->GetTrack( lBachIdx );
          if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
               AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
               continue;
          }
          // - Get the TPCnumber of cluster for the daughters
          lPosTPCClusters   = pTrackXi->GetTPCNcls();
          lNegTPCClusters   = nTrackXi->GetTPCNcls();
          lBachTPCClusters  = bachTrackXi->GetTPCNcls();

          //-------------------------------------
          // - Rejection of a poor quality tracks
          if (fkQualityCutTPCrefit) {
                // - Poor quality related to TPCrefit
                ULong_t pStatus    = pTrackXi->GetStatus();
                ULong_t nStatus    = nTrackXi->GetStatus();
                ULong_t bachStatus = bachTrackXi->GetStatus();
                if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
          }
          if (fkQualityCutnTPCcls) {
                // - Poor quality related to TPC clusters
                if (lPosTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
                if (lNegTPCClusters  < fMinnTPCcls) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
                if (lBachTPCClusters < fMinnTPCcls) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
          }

          //-----------------------------------
          const AliExternalTrackParam *pExtTrack    = pTrackXi->GetInnerParam();
          const AliExternalTrackParam *nExtTrack    = nTrackXi->GetInnerParam();
          const AliExternalTrackParam *bachExtTrack = bachTrackXi->GetInnerParam();
          if (pExtTrack) {
                lInnerWallMomCascDghters[0] = pExtTrack->GetP() * pExtTrack->Charge();
                lTPCSignalCascDghters   [0] = pTrackXi->GetTPCsignal();
          }
          if (nExtTrack) {
                lInnerWallMomCascDghters[1] = nExtTrack->GetP() * nExtTrack->Charge();
                lTPCSignalCascDghters   [1] = nTrackXi->GetTPCsignal();
           }
           if (bachExtTrack) {
                lInnerWallMomCascDghters[2] = bachExtTrack->GetP() * bachExtTrack->Charge();
                lTPCSignalCascDghters   [2] = bachTrackXi->GetTPCsignal();
           }
           etaPos  = pTrackXi->Eta();
           etaNeg  = nTrackXi->Eta();
           etaBach = bachTrackXi->Eta();
           lInvMassLambdaAsCascDghter = xi->GetEffMass(); //This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
           lDcaV0DaughtersXi 	      = xi->GetDcaV0Daughters(); 
           lV0CosineOfPointingAngleXi = xi->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],
                                                                       lBestPrimaryVtxPos[1],
                                                                       lBestPrimaryVtxPos[2] );
           lDcaV0ToPrimVertexXi = xi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
           lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); //Note: AliExternalTrackParam::GetD returns an algebraic value ...
           xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
           lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );	
           lDcaPosToPrimVertexXi = TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); 
           lDcaNegToPrimVertexXi = TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lMagneticField ) ); 

           //-----------------------------------------
	   // - Extra-selection for cascade candidates
           if (fkExtraSelections) { //in AliCascadeVertexer
               if (lDcaXiDaughters > 0.3) continue;  
               if (lXiCosineOfPointingAngle < 0.999 ) continue; 
               if (lDcaV0ToPrimVertexXi < 0.05) continue; 
               if (lDcaBachToPrimVertexXi < 0.03) continue; 
               //if (TMath::Abs(lInvMassLambdaAsCascDghter-1.11568) > 0.006 ) continue;   
               if (lDcaV0DaughtersXi > 1.) continue;  
               if (lV0CosineOfPointingAngleXi < 0.998) continue; 
               if (lDcaPosToPrimVertexXi < 0.1) continue; 
               if (lDcaNegToPrimVertexXi < 0.1) continue; 
	       if (lXiRadius < .9) continue; 
               //if (lXiRadius > 100) continue; 
	       if (lV0RadiusXi < 0.9) continue;   
               //if (lV0RadiusXi > 100) continue; 
           }

           //----------------------------------------------------------------------------------------------------	
           // - Around effective masses. Change mass hypotheses to cover all the possibilities:  Xi-/+, Omega -/+
           if ( bachTrackXi->Charge() < 0 )	{
                //Calculate the effective mass of the Xi- candidate: Xi- hyp. (pdg code 3312)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3312); 	
                lInvMassXiMinus = xi->GetEffMassXi();
                //Calculate the effective mass of the Xi- candidate: Omega- hyp. (pdg code 3334)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3334); 	
                lInvMassOmegaMinus = xi->GetEffMassXi();
                //Back to "default" hyp. (Xi-)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , 3312); 
           } // end if negative bachelor
           if ( bachTrackXi->Charge() >  0 ) {
                //Calculate the effective mass of the Xi- candidate: Xi+ hyp. (pdg code -3312)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3312); 	
                lInvMassXiPlus = xi->GetEffMassXi();
                //Calculate the effective mass of the Xi- candidate: Omega+ hyp. (pdg code -3334)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3334); 	
                lInvMassOmegaPlus = xi->GetEffMassXi();
                //Back to "default" hyp. (Xi-)
                lV0quality = 0.;
                xi->ChangeMassHypothesis(lV0quality , -3312); 	
           } // end if positive bachelor

           //--------------------------------
           // - PID on the daughter tracks
	   // - Combined PID ->  removed, add when will be used

	   // - TPC PID : 3-sigma bands on Bethe-Bloch curve
           //Bachelor
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
           //Negative V0 daughter
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
           //Positive V0 daughter
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
           if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
           /*
           const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXi    ->GetInnerParam();
           const AliExternalTrackParam *nInnerWallTrackXi    = nTrackXi    ->GetInnerParam();
           const AliExternalTrackParam *bachInnerWallTrackXi = bachTrackXi ->GetInnerParam();
           if (pInnerWallTrackXi && nInnerWallTrackXi && bachInnerWallTrackXi ) { 
                Double_t pMomInnerWall    = pInnerWallTrackXi   ->GetP();
                Double_t nMomInnerWall    = nInnerWallTrackXi   ->GetP();
                Double_t bachMomInnerWall = bachInnerWallTrackXi->GetP();
                //Bachelor
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 3)                              lIsBachelorPionForTPC = kTRUE;
                if (bachMomInnerWall < 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 5) lIsBachelorKaonForTPC = kTRUE;
                if (bachMomInnerWall > 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;                
                //Negative V0 daughter
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 3  )                           lIsNegPionForTPC   = kTRUE;
                if (nMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 5 )   lIsNegProtonForTPC = kTRUE;
                if (nMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 3 )   lIsNegProtonForTPC = kTRUE;       
                //Positive V0 daughter
                if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 3 )                            lIsPosPionForTPC   = kTRUE;
                if (pMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 5)     lIsPosProtonForTPC = kTRUE;
                if (pMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 3)     lIsPosProtonForTPC = kTRUE;
           }*/
		
           //---------------------------------
           // - Extra info for QA (ESD)
           // Miscellaneous pieces of info that may help regarding data quality assessment.
           // Cascade transverse and total momentum
           xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
           lXiTransvMom	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY );
           lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ );
           // V0 total momentum
           xi->GetNPxPyPz(lV0NMomX,lV0NMomY,lV0NMomZ);
           xi->GetPPxPyPz(lV0PMomX,lV0PMomY,lV0PMomZ);
           lV0TotMom = TMath::Sqrt(TMath::Power(lV0NMomX+lV0PMomX,2) + TMath::Power(lV0NMomY+lV0PMomY,2) + TMath::Power(lV0NMomZ+lV0PMomZ,2));
           // Bachelor total momentum
           xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
           lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY );
           lBachTotMom     = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY + lBachMomZ*lBachMomZ );
           lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
           lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
           lChargeXi       = xi->Charge();
           lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
           lRapXi    = xi->RapXi();
           lRapOmega = xi->RapOmega();
           lEta      = xi->Eta();
           lTheta    = xi->Theta()*180.0/TMath::Pi();
           lPhi      = xi->Phi()*180.0/TMath::Pi();
           lAlphaXi  = xi->AlphaXi();
           lPtArmXi  = xi->PtArmXi();
	   // Extra-cut = Anti-splitting cut for lambda daughters
           Bool_t kAntiSplittingLambda = kFALSE;	
           if (kAntiSplittingLambda) { // not used
               Double_t lNMomX = 0., lNMomY = 0., lNMomZ = 0.;
               Double_t lPMomX = 0., lPMomY = 0., lPMomZ = 0.;
               xi->GetPPxPyPz(lPMomX, lPMomY, lPMomZ); 
               xi->GetNPxPyPz(lNMomX, lNMomY, lNMomZ); 
               if ( xi->Charge() < 0) {// Xi- or Omega-
                   if (TMath::Abs(lBachTransvMom - TMath::Sqrt( lNMomX*lNMomX + lNMomY*lNMomY )  ) < 0.075) continue;
	       } else {                //Xi+ or Omega+
	           if(TMath::Abs(lBachTransvMom - TMath::Sqrt( lPMomX*lPMomX + lPMomY*lPMomY ) ) < 0.075) continue;
	       }
           }

    } // end of ESD treatment
 
    else if (fAnalysisType == "AOD") {

           // -------------------------------------
           // - Load the cascades from the handler
           const AliAODcascade *xi = lAODevent->GetCascade(iXi);
           if (!xi) continue;
		
           //----------------------------------------------------------------------------        
           // - Assigning the necessary variables for specific AliESDcascade data members		
           lEffMassXi        	   = xi->MassXi(); // default working hypothesis : cascade = Xi- decay
           lDcaXiDaughters	   = xi->DcaXiDaughters();
           lXiCosineOfPointingAngle = xi->CosPointingAngleXi( lBestPrimaryVtxPos[0], 
		  					      lBestPrimaryVtxPos[1], 
		 					      lBestPrimaryVtxPos[2] );
           lPosXi[0] = xi->DecayVertexXiX();
           lPosXi[1] = xi->DecayVertexXiY();
           lPosXi[2] = xi->DecayVertexXiZ();
           lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1] );

           //-------------------------------------------------------------------------------------------------------------------------------
           // - Around the tracks: Bach + V0 (ESD). Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
           AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
           AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
           AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );
           if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
                AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
                continue;
           }
           UInt_t lIdxPosXi  = (UInt_t) TMath::Abs( pTrackXi->GetID() );  
           UInt_t lIdxNegXi  = (UInt_t) TMath::Abs( nTrackXi->GetID() );
           UInt_t lBachIdx   = (UInt_t) TMath::Abs( bachTrackXi->GetID() );
                                                // Care track label can be negative in MC production (linked with the track quality)
                                                // However = normally, not the case for track index ...
           // - Rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
           if (lBachIdx == lIdxNegXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue; }
           if (lBachIdx == lIdxPosXi) { AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue; }
           // - Get the TPCnumber of cluster for the daughters
           lPosTPCClusters   = pTrackXi->GetTPCNcls(); // FIXME: Is this ok? or something like in LambdaK0PbPb task AOD?
           lNegTPCClusters   = nTrackXi->GetTPCNcls();
           lBachTPCClusters  = bachTrackXi->GetTPCNcls();

           //-------------------------------------
           // - Rejection of a poor quality tracks
           if (fkQualityCutTPCrefit) {
                // - Poor quality related to TPCrefit
                if (!(pTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
                if (!(nTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
                if (!(bachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
           }
           if (fkQualityCutnTPCcls) {
                // - Poor quality related to TPC clusters
                if (lPosTPCClusters  < fMinnTPCcls) { //AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!");
                    continue; }
                if (lNegTPCClusters  < fMinnTPCcls) { //AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!");
                    continue; }
                if (lBachTPCClusters < fMinnTPCcls) { //AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!");
                    continue; }
           }

           //---------------------------------------
           // - Around the tracks: Bach + V0 (AOD). Necessary variables for AODcascade data members coming from the AODv0 part (inheritance)
           etaPos  = pTrackXi->Eta();
           etaNeg  = nTrackXi->Eta();
           etaBach = bachTrackXi->Eta();
           lChargeXi = xi->ChargeXi();
           if ( lChargeXi < 0) lInvMassLambdaAsCascDghter = xi->MassLambda();
           else 	       lInvMassLambdaAsCascDghter = xi->MassAntiLambda();
           lDcaV0DaughtersXi  	  = xi->DcaV0Daughters(); 
           lDcaV0ToPrimVertexXi   = xi->DcaV0ToPrimVertex();
           lDcaBachToPrimVertexXi = xi->DcaBachToPrimVertex(); 
           lPosV0Xi[0] = xi->DecayVertexV0X();
           lPosV0Xi[1] = xi->DecayVertexV0Y();
           lPosV0Xi[2] = xi->DecayVertexV0Z(); 
           lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );
           lV0CosineOfPointingAngleXi = xi->CosPointingAngle( lBestPrimaryVtxPos ); 
           lDcaPosToPrimVertexXi      = xi->DcaPosToPrimVertex(); 
           lDcaNegToPrimVertexXi      = xi->DcaNegToPrimVertex(); 

           //----------------------------------------------------------------------------------------------------       
           // - Around effective masses. Change mass hypotheses to cover all the possibilities:  Xi-/+, Omega -/+
           if ( lChargeXi < 0 )	lInvMassXiMinus	   = xi->MassXi();
           if ( lChargeXi > 0 )	lInvMassXiPlus 	   = xi->MassXi();
           if ( lChargeXi < 0 )	lInvMassOmegaMinus = xi->MassOmega();
           if ( lChargeXi > 0 )	lInvMassOmegaPlus  = xi->MassOmega();

           //--------------------------------
           // - PID on the daughter tracks
           // - Combined PID ->  removed, add when will be used

           // - TPC PID : 3-sigma bands on Bethe-Bloch curve
           //Bachelor
           if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
           if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
           //Negative V0 daughter
           if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
           if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
           //Positive V0 daughter
           if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
           if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;

           //---------------------------------
           // - Extra info for QA (AOD)
           // Miscellaneous pieces of info that may help regarding data quality assessment.
           // Cascade transverse and total momentum	
           lXiMomX = xi->MomXiX();
           lXiMomY = xi->MomXiY();
           lXiMomZ = xi->MomXiZ();
           lXiTransvMom = TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY );
           lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ );
           Double_t lV0MomX = xi->MomV0X();
           Double_t lV0MomY = xi->MomV0Y();
           Double_t lV0MomZ = xi->MomV0Z();
           lV0TotMom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2)+TMath::Power(lV0MomZ,2));
           lBachMomX = xi->MomBachX();
           lBachMomY = xi->MomBachY();
           lBachMomZ = xi->MomBachZ();		
           lBachTransvMom = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY );
           lBachTotMom    = TMath::Sqrt( lBachMomX*lBachMomX + lBachMomY*lBachMomY + lBachMomZ*lBachMomZ );
           lV0NMomX = xi->MomNegX();
           lV0NMomY = xi->MomNegY();
           lV0PMomX = xi->MomPosX();
           lV0PMomY = xi->MomPosY();
           lnTrackTransvMom = TMath::Sqrt( lV0NMomX*lV0NMomX   + lV0NMomY*lV0NMomY );
           lpTrackTransvMom = TMath::Sqrt( lV0PMomX*lV0PMomX   + lV0PMomY*lV0PMomY );
           lV0toXiCosineOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );
           lRapXi    = xi->RapXi();
           lRapOmega = xi->RapOmega();
           lEta      = xi->Eta();	              	// Will not work ! need a method Pz(), Py() Px() 
           lTheta    = xi->Theta() *180.0/TMath::Pi();  // in AODcascade.
           lPhi      = xi->Phi()   *180.0/TMath::Pi();  // Here, we will get eta, theta, phi for the V0 ...
           lAlphaXi  = xi->AlphaXi();
           lPtArmXi  = xi->PtArmXi();

    } // end of AOD treatment

    // Cut on pt of the three daughter tracks
    if (lBachTransvMom<fMinPtCutOnDaughterTracks) continue;
    if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
    if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
      
    // Cut on pseudorapidity of the three daughter tracks
    if (TMath::Abs(etaBach)>fEtaCutOnDaughterTracks) continue;
    if (TMath::Abs(etaPos)>fEtaCutOnDaughterTracks) continue;
    if (TMath::Abs(etaNeg)>fEtaCutOnDaughterTracks) continue;
      
      
    //----------------------------------
    // Calculate proper lenght for cascade
    Double_t cascadeMass = 0.;
    if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) cascadeMass = 1.321;
    if ( ( (lChargeXi<0) && lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )  ) cascadeMass = 1.672;
    Double_t lctau =  TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power(( lPosXi[2]-lBestPrimaryVtxPos[2]),2));
    if (lXiTotMom!=0) lctau = lctau*cascadeMass/lXiTotMom;
    else lctau = -1.;

    //-------------------------------------------------
    // Calculate proper lenght for Lambda (reconstructed)
    Float_t lambdaMass = 1.115683; // PDG mass
    Float_t distV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2));
    Float_t lctauV0 = -1.;
    if (lV0TotMom!=0) lctauV0 = distV0Xi*lambdaMass/lV0TotMom;
    Float_t distTV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));

    //--------------
    /*// For AliEVE      
         if(lChargeXi < 0&& lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) {
             if (lXiTransvMom>2.&&lXiTransvMom<4.&&(lInvMassXiMinus<1.322&&lInvMassXiMinus>1.320)&&(lXiRadius<8.&&lXiRadius>3.)) {
                         // FIXME : Just to know which file is currently open : locate the file containing Xi
                  cout << "Name of the file containing Xi candidate(s) :" 
                       << CurrentFileName() 
                       << " / entry: "     << Entry()
                       << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
                       << " AliESDcascade number " << iXi 
                       << " : mass(Xi-) = " << lInvMassXiMinus
                       << " / charge = "   << lChargeXi
                       << " / pt(Casc) = " << lXiTransvMom
                       << " / Decay 2d R(Xi) = " << lXiRadius 
                       << endl;
             }
         }
         if(lChargeXi < 0&& lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) {
             if (lXiTransvMom>2&&lXiTransvMom<4&&(lInvMassOmegaMinus<1.674&&lInvMassOmegaMinus>1.670)&&(lXiRadius<8.&&lXiRadius>3.)) {
                  cout << "Name of the file containing Omega candidate(s) :"
                       << CurrentFileName()
                       << " / entry: "     << Entry()
                       << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
                       << " AliESDcascade number " << iXi 
                       << " : mass(Omega-) = " << lInvMassOmegaMinus
                       << " / charge = "   << lChargeXi
                       << " / pt(Casc) = " << lXiTransvMom
                       << " / Decay 2d R(Xi) = " << lXiRadius
                       << endl;

             }
         }*/
          

    // - 
    fHistPosV0TPCClusters->Fill( lPosTPCClusters );
    fHistNegV0TPCClusters->Fill( lNegTPCClusters );
    fHistBachTPCClusters->Fill( lBachTPCClusters );
    f2dHistTPCdEdxOfCascDghters->Fill( lInnerWallMomCascDghters[0] , lTPCSignalCascDghters[0] );
    f2dHistTPCdEdxOfCascDghters->Fill( lInnerWallMomCascDghters[1] , lTPCSignalCascDghters[1] );
    f2dHistTPCdEdxOfCascDghters->Fill( lInnerWallMomCascDghters[2] , lTPCSignalCascDghters[2] );

    //----------------
    //Plot with PID on  
    if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) {
                                       // NOTE : 
                                       // with this condition, it could happen that a cascade candidate satisfies the wrong requirement,
                                       // e.g. one looks at a Xi- candidate for which lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC = kFALSE
                                       // Expectation: it should be excluded, but lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC = kTRUE
                                       // then this bad Xi-candidate will contribute anyway (OR condition).
                                       // Hence: the extra condition on the sign of the Cascade
           //if (TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.010 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.010)
           fHistEffMassXi->Fill( lEffMassXi );
           fHistDcaXiDaughters->Fill( lDcaXiDaughters );         	    // Flag CascadeVtxer: Cut Variable e 
           fHistDcaBachToPrimVertex->Fill( lDcaBachToPrimVertexXi   );	    // Flag CascadeVtxer: Cut Variable d
           fHistXiCosineOfPointingAngle->Fill( lXiCosineOfPointingAngle );  // Flag CascadeVtxer: Cut Variable f
           fHistXiRadius->Fill( lXiRadius );	                            // Flag CascadeVtxer: Cut Variable g+h
           fHistMassLambdaAsCascDghter->Fill( lInvMassLambdaAsCascDghter ); // Flag CascadeVtxer: Cut Variable c
           fHistDcaV0DaughtersXi->Fill( lDcaV0DaughtersXi );
           fHistV0CosineOfPointingAngleXi->Fill( lV0CosineOfPointingAngleXi ); 
           fHistV0RadiusXi->Fill( lV0RadiusXi );
           fHistDcaV0ToPrimVertexXi->Fill( lDcaV0ToPrimVertexXi );	    // Flag CascadeVtxer: Cut Variable b
           fHistDcaPosToPrimVertexXi->Fill( lDcaPosToPrimVertexXi );
           fHistDcaNegToPrimVertexXi->Fill( lDcaNegToPrimVertexXi );
           fHistChargeXi->Fill( lChargeXi );
           fHistV0toXiCosineOfPointingAngle->Fill( lV0toXiCosineOfPointingAngle );
           if ( TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.012 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.012) { // One InvMass should be different from 0
                fHistXiTransvMom->Fill( lXiTransvMom );
                fHistXiTotMom->Fill( lXiTotMom );
                fHistBachTransvMomXi->Fill( lBachTransvMom );
                fHistBachTotMomXi->Fill( lBachTotMom );
                fHistRapXi->Fill( lRapXi );
                fHistEtaXi->Fill( lEta );
                if (lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) {
                      fHistEtaBachXi->Fill( etaBach );
                      fHistEtaPosXi->Fill( etaPos );
                      fHistEtaNegXi->Fill( etaNeg );
                }
                fHistThetaXi->Fill( lTheta );
                fHistPhiXi->Fill( lPhi );
           }
           if ( TMath::Abs( lInvMassOmegaMinus-1.672 ) < 0.012 || TMath::Abs( lInvMassOmegaPlus-1.672 ) < 0.012 ) { // One InvMass should be different from 0
                fHistRapOmega->Fill( lRapOmega ); 
           }
           f2dHistArmenteros->Fill( lAlphaXi, lPtArmXi );
    } // end with PID ...

    //-----------------------
    // - Invariant mass plots
    //Plots 1D
    if ( lChargeXi < 0 ) {
         fHistMassXiMinus->Fill( lInvMassXiMinus );
         fHistMassOmegaMinus->Fill( lInvMassOmegaMinus );
         f2dHistDcaXiDaughtersvsInvMass->Fill(lDcaXiDaughters,lInvMassXiMinus);
         f2dHistDcaBachToPrimVertexvsInvMass->Fill(lDcaBachToPrimVertexXi,lInvMassXiMinus); 
         f2dHistXiCosineOfPointingAnglevsInvMass->Fill(lXiCosineOfPointingAngle,lInvMassXiMinus);
         f2dHistMassLambdaAsCascDghtervsInvMass->Fill(lInvMassLambdaAsCascDghter,lInvMassXiMinus);
         f2dHistDcaV0DaughtersXivsInvMass->Fill(lDcaV0DaughtersXi,lInvMassXiMinus);
         f2dHistDcaV0ToPrimVertexXivsInvMass->Fill(lDcaV0ToPrimVertexXi,lInvMassXiMinus);
    }
    if ( lChargeXi > 0 ) {
      fHistMassXiPlus->Fill( lInvMassXiPlus );
      fHistMassOmegaPlus->Fill( lInvMassOmegaPlus );
    }
    //Plots 2D, 3D
    if ( lChargeXi < 0 ) {
      f2dHistEffMassLambdaVsEffMassXiMinus->Fill( lInvMassLambdaAsCascDghter, lInvMassXiMinus ); 
      f2dHistEffMassXiVsEffMassOmegaMinus ->Fill( lInvMassXiMinus, lInvMassOmegaMinus );
      f2dHistXiRadiusVsEffMassXiMinus     ->Fill( lXiRadius, lInvMassXiMinus );
      f2dHistXiRadiusVsEffMassOmegaMinus  ->Fill( lXiRadius, lInvMassOmegaMinus );
    } else {
      f2dHistEffMassLambdaVsEffMassXiPlus ->Fill( lInvMassLambdaAsCascDghter, lInvMassXiPlus );
      f2dHistEffMassXiVsEffMassOmegaPlus  ->Fill( lInvMassXiPlus, lInvMassOmegaPlus );
      f2dHistXiRadiusVsEffMassXiPlus      ->Fill( lXiRadius, lInvMassXiPlus);
      f2dHistXiRadiusVsEffMassOmegaPlus   ->Fill( lXiRadius, lInvMassOmegaPlus );
    }

    //---------------------------------------------	
    // - Filling the AliCFContainers related to PID
    Double_t lContainerPIDVars[3] = {0.0};
    // Xi Minus		
    if ( lChargeXi < 0 ) {
          lContainerPIDVars[0] = lXiTransvMom;
          lContainerPIDVars[1] = lInvMassXiMinus;
          lContainerPIDVars[2] = lRapXi;
          //No PID
          fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 0); // No PID
	  //TPC PID
          if ( lIsBachelorPionForTPC )                                           fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track		
          if ( lIsBachelorPionForTPC && lIsPosProtonForTPC )                     fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
          if ( lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	  //Combined PID
          if ( lIsBachelorPion )                                      fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor	
          if ( lIsBachelorPion && lIsPosInXiProton )                  fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
          if (lIsBachelorPion && lIsPosInXiProton && lIsNegInXiPion ) fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 	
    // Xi Plus		
    if ( lChargeXi > 0 ) {
          lContainerPIDVars[0] = lXiTransvMom;
          lContainerPIDVars[1] = lInvMassXiPlus;
          lContainerPIDVars[2] = lRapXi;
	  //No PID
          fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 0); // No PID
          //TPC PID
          if ( lIsBachelorPionForTPC )                                           fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
          if ( lIsBachelorPionForTPC && lIsNegProtonForTPC )                     fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
          if ( lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC ) fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	  //Combined PID
          if ( lIsBachelorPion )                                      fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
          if ( lIsBachelorPion && lIsNegInXiProton )                  fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
          if (lIsBachelorPion && lIsNegInXiProton && lIsPosInXiPion ) fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; 
    // Omega Minus		
    if ( lChargeXi < 0 ) {
          lContainerPIDVars[0] = lXiTransvMom;
          lContainerPIDVars[1] = lInvMassOmegaMinus;
          lContainerPIDVars[2] = lRapOmega;
	  //No PID
          fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
       	  //TPC PID
          if ( lIsBachelorKaonForTPC )                                           fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
          if ( lIsBachelorKaonForTPC && lIsPosProtonForTPC )                     fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
          if ( lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	  //Combined PID
          if ( lIsBachelorKaon )                                            fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
          if ( lIsBachelorKaon && lIsPosInOmegaProton )     	            fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
          if (lIsBachelorKaon && lIsPosInOmegaProton && lIsNegInOmegaPion ) fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.;
    // Omega Plus		
    if ( lChargeXi > 0 ) {
      lContainerPIDVars[0] = lXiTransvMom;
      lContainerPIDVars[1] = lInvMassOmegaPlus;
      lContainerPIDVars[2] = lRapOmega;
	// No PID
      fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 0); // No PID
	// TPC PID
      if ( lIsBachelorKaonForTPC  )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
      if( lIsBachelorKaonForTPC && 
	  lIsNegProtonForTPC     )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
      if ( lIsBachelorKaonForTPC && 
	   lIsNegProtonForTPC    && 
	   lIsPosPionForTPC       )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
	// Combined PID
      if ( lIsBachelorKaon        )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
      if ( lIsBachelorKaon       && 
           lIsNegInOmegaProton    )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
      if (lIsBachelorKaon     && 
	  lIsNegInOmegaProton && 
	  lIsPosInOmegaPion    )
	fCFContCascadePIDOmegaPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
        	
    //--------------------------------------------------------------------
    // Filling the AliCFContainer (optimisation of topological selections)
    Double_t lContainerCutVars[19] = {0.0};
                        
    lContainerCutVars[0]  = lDcaXiDaughters;
    lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
    lContainerCutVars[2]  = lXiCosineOfPointingAngle;
    lContainerCutVars[3]  = lXiRadius;
    lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
    lContainerCutVars[5]  = lDcaV0DaughtersXi;
    lContainerCutVars[6]  = lV0toXiCosineOfPointingAngle;
    lContainerCutVars[7]  = lV0RadiusXi;
    lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
    lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
    lContainerCutVars[10] = lDcaNegToPrimVertexXi;
    lContainerCutVars[13] = lXiTransvMom;
    lContainerCutVars[16] = lctau;
    lContainerCutVars[17] = lctauV0;
    lContainerCutVars[18] = distTV0Xi;
 
    if ( lChargeXi < 0 ) {
         lContainerCutVars[11] = lInvMassXiMinus;
         lContainerCutVars[12] = lInvMassOmegaMinus;
         lContainerCutVars[14] = lRapXi;
         lContainerCutVars[15] = -1.;
         if (lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
         lContainerCutVars[11] = lInvMassXiMinus;
         lContainerCutVars[12] = lInvMassOmegaMinus;
         lContainerCutVars[14] = -1.;
         lContainerCutVars[15] = lRapOmega;
         if (lIsBachelorKaonForTPC && lIsPosProtonForTPC && lIsNegPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,2); // for Omega-
    } else {
         lContainerCutVars[11] = lInvMassXiPlus;
         lContainerCutVars[12] = lInvMassOmegaPlus; 
         lContainerCutVars[14] = lRapXi;
         lContainerCutVars[15] = -1.; 
         if (lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,1); // for Xi+
         lContainerCutVars[11] = lInvMassXiPlus;
         lContainerCutVars[12] = lInvMassOmegaPlus;
         lContainerCutVars[14] = -1.;
         lContainerCutVars[15] = lRapOmega;
         if (lIsBachelorKaonForTPC && lIsNegProtonForTPC && lIsPosPionForTPC) fCFContCascadeCuts->Fill(lContainerCutVars,3); // for Omega+ 
    }                 
  } //end of the Cascade loop (ESD or AOD)
    
  // Post output data.
 PostData(1, fListHistCascade);
 PostData(2, fCFContCascadePIDXiMinus);
 PostData(3, fCFContCascadePIDXiPlus);
 PostData(4, fCFContCascadePIDOmegaMinus);
 PostData(5, fCFContCascadePIDOmegaPlus);
 PostData(6, fCFContCascadeCuts);
}

//________________________________________________________________________
Int_t AliAnalysisTaskCheckCascadepp276::DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent) {
    // Checking the number of tracks with TPCrefit for each event
    // Needed for a rough assessment of the event multiplicity
        
    Int_t nTrackWithTPCrefitMultiplicity = 0;
    for (Int_t iTrackIdx = 0; iTrackIdx < (InputEvent())->GetNumberOfTracks(); iTrackIdx++) {
      AliESDtrack *esdTrack	= 0x0;
      esdTrack	= lESDevent->GetTrack( iTrackIdx );
      if (!esdTrack) { AliWarning("Pb / Could not retrieve one track within the track loop for TPCrefit check ..."); continue; }

      ULong_t lTrackStatus    = esdTrack->GetStatus();
      if ((lTrackStatus&AliESDtrack::kTPCrefit)    == 0) continue;
      else nTrackWithTPCrefitMultiplicity++;
        // FIXME :
        // The goal here is to get a better assessment of the event multiplicity.
        // (InputEvent())->GetNumberOfTracks() takes into account ITS std alone tracks + global tracks
        // This may introduce a bias. Hence the number of TPC refit tracks.
        // Note : the event multiplicity = analysis on its own... See Jacek's or Jan Fiete's analysis on dN/d(eta)

    } // end loop over all event tracks
    return  nTrackWithTPCrefitMultiplicity;
}


//________________________________________________________________________
void AliAnalysisTaskCheckCascadepp276::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

/*  TList *cRetrievedList = 0x0;
         cRetrievedList = (TList*)GetOutputData(1);
	if(!cRetrievedList){
		AliWarning("ERROR - AliAnalysisTaskCheckCascadepp276: ouput data container list not available\n"); return;
	}
  fHistTrackMultiplicity = dynamic_cast<TH1F*> (   cRetrievedList->FindObject("fHistTrackMultiplicity") );
  if (!fHistTrackMultiplicity) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadepp276: fHistTrackMultiplicity not available\n"); return;
	}
  fHistMassXiMinus    = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiMinus") );	
	if (!fHistMassXiMinus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadepp276: fHistMassXiMinus not available\n"); return;
	}
  fHistMassXiPlus     = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiPlus") );
	if (!fHistMassXiPlus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadepp276: fHistMassXiPlus not available\n"); return;
	}	
  fHistMassOmegaMinus = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaMinus") );
	if (!fHistMassOmegaMinus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadepp276: fHistMassOmegaMinus not available\n"); return;
	}
  fHistMassOmegaPlus  = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaPlus") );	
	if (!fHistMassOmegaPlus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadepp276: fHistMassOmegaPlus not available\n"); return;
	}
  
  TCanvas *canCheckCascade = new TCanvas("AliAnalysisTaskCheckCascadep276","CheckCascade overview",10,10,1010,660);
  canCheckCascade->Divide(2,2);
  
  canCheckCascade->cd(1);
  canCheckCascade->cd(1)->SetLogy();
  fHistTrackMultiplicity->SetMarkerStyle(kFullStar);  
  fHistTrackMultiplicity->GetXaxis()->SetLabelFont(42);
  fHistTrackMultiplicity->GetYaxis()->SetLabelFont(42);
  fHistTrackMultiplicity->SetTitleFont(42, "xy");
  fHistTrackMultiplicity->GetXaxis()->SetTitleOffset(1.1);
  fHistTrackMultiplicity->DrawCopy("H");
  
  canCheckCascade->cd(2);  
  fHistMassXiMinus ->SetMarkerStyle(kFullCircle);
  fHistMassXiMinus ->SetMarkerSize(0.5);
  fHistMassXiMinus ->GetXaxis()->SetLabelFont(42);
  fHistMassXiMinus ->GetYaxis()->SetLabelFont(42);
  fHistMassXiMinus ->SetTitleFont(42, "xy");
  fHistMassXiMinus ->GetXaxis()->SetTitleOffset(1.1);
  fHistMassXiMinus ->GetYaxis()->SetTitleOffset(1.3);
  //fHistMassXiMinus->Rebin(2);
  fHistMassXiMinus ->GetXaxis()->SetRangeUser(1.24, 1.42);
  fHistMassXiMinus ->DrawCopy("E");
  
  fHistMassXiPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassXiPlus ->SetMarkerColor(kRed+2);
  fHistMassXiPlus ->SetLineColor(kRed+2);
  fHistMassXiPlus ->SetMarkerSize(0.5);
  //fHistMassXiPlus ->Rebin(2);
  fHistMassXiPlus ->DrawCopy("ESAME");
  
  
  TLegend *legendXi =new TLegend(0.67,0.34,0.97,0.54);
 		legendXi->SetTextFont(42);
 		legendXi->SetTextSize(0.05);
 		legendXi->SetFillColor(kWhite);
 		legendXi->AddEntry( fHistMassXiMinus,"#Xi^{-} candidates","lp");
 		legendXi->AddEntry( fHistMassXiPlus,"#Xi^{+} candidates","lp");
 		legendXi->Draw();
  
  
  canCheckCascade->cd(3);  
  fHistMassOmegaPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassOmegaPlus ->SetMarkerColor(kRed+2);
  fHistMassOmegaPlus ->SetLineColor(kRed+2);
  fHistMassOmegaPlus ->SetMarkerSize(0.5);
  fHistMassOmegaPlus ->GetXaxis()->SetLabelFont(42);
  fHistMassOmegaPlus ->GetYaxis()->SetLabelFont(42);
  fHistMassOmegaPlus ->SetTitleFont(42, "xy");
  fHistMassOmegaPlus ->GetXaxis()->SetTitleOffset(1.1);
  fHistMassOmegaPlus ->GetYaxis()->SetTitleOffset(1.25);
  //fHistMassOmegaPlus ->Rebin(2);
  fHistMassOmegaPlus ->GetXaxis()->SetRangeUser(1.6, 1.84);
  fHistMassOmegaPlus ->DrawCopy("E");
  
  fHistMassOmegaMinus->SetMarkerStyle(kFullCircle);
  fHistMassOmegaMinus->SetMarkerSize(0.5);
  //fHistMassOmegaMinus->Rebin(2);
  fHistMassOmegaMinus->DrawCopy("ESAME");

  
   TLegend *legendOmega = new TLegend(0.67,0.34,0.97,0.54);
 		legendOmega->SetTextFont(42);
 		legendOmega->SetTextSize(0.05);
 		legendOmega->SetFillColor(kWhite);
 		legendOmega->AddEntry( fHistMassOmegaMinus,"#Omega^{-} candidates","lp");
 		legendOmega->AddEntry( fHistMassOmegaPlus,"#Omega^{+} candidates","lp");
 		legendOmega->Draw();
     */
}
