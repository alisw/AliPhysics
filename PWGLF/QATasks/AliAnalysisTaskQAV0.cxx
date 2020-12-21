/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// QA Task designed to investigate V0 characteristics in any data
// 
// Stripped down and adapted from AliAnalysisTaskExtractV0: 
// --- smaller, more convenient output for debugging
// --- optional TH3 list of histos to enable a light-weight analysis
// --- This code is under development, so...
//
//  Please Report Any Bugs! 
//
//   --- David Dobrigkeit Chinellato
//        (david.chinellato@gmail.com)
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

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
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliLightV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliESDHeader.h"

#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskQAV0.h"
#include "AliMultSelectionTask.h"

//debugging purposes
#include "TObjectTable.h"

ClassImp(AliAnalysisTaskQAV0)

AliAnalysisTaskQAV0::AliAnalysisTaskQAV0() 
  : AliAnalysisTaskSE(), 
  //Output lists
  fOutput(0), 

  //Histos
  fHistEvent(0),
  fHistTopDCANegToPV(0),
  fHistTopDCAPosToPV(0),
  fHistTopDCAV0Daughters(0),
  fHistTopCosinePA(0),
  fHistTopV0Radius(0),
  fHistSelectedTopDCANegToPV(0),
  fHistSelectedTopDCAPosToPV(0),
  fHistSelectedTopDCAV0Daughters(0),
  fHistSelectedTopCosinePA(0),
  fHistSelectedTopV0Radius(0),

  f2dHistInvMassK0Short(0),
  f2dHistInvMassLambda(0),
  f2dHistInvMassAntiLambda(0),

  f2dHistInvMassWithdEdxK0Short(0),
  f2dHistInvMassWithdEdxLambda(0),
  f2dHistInvMassWithdEdxAntiLambda(0),

  f2dHistInvMassWithdEdxCowboysK0Short(0),
  f2dHistInvMassWithdEdxCowboysLambda(0),
  f2dHistInvMassWithdEdxCowboysAntiLambda(0),

  f2dHistInvMassWithdEdxSailorsK0Short(0),
  f2dHistInvMassWithdEdxSailorsLambda(0),
  f2dHistInvMassWithdEdxSailorsAntiLambda(0),

  f2dHistResponseNegativeAsPion(0),
  f2dHistResponseNegativeAsProton(0),
  f2dHistResponsePositiveAsPion(0),
  f2dHistResponsePositiveAsProton(0),

  f2dHistdEdxSignalPionFromLambda(0),
  f2dHistdEdxSignalProtonFromLambda(0),
  f2dHistResponsePionFromLambda(0),
  f2dHistResponseProtonFromLambda(0),

fTrigType(AliVEvent::kINT7),

  //Task Control / Utils
  fPIDResponse(0),
  fkRunV0Vertexer ( kFALSE ),
  fkUseLightV0Vertexer ( kFALSE ),
  fdEdxCut (3) 
{
  // Dummy Constructor
  for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
}

AliAnalysisTaskQAV0::AliAnalysisTaskQAV0(const char *name) 
  : AliAnalysisTaskSE(name), 
  //Output lists
  fOutput(0), 

  //Histos
  fHistEvent(0),
  fHistTopDCANegToPV(0),
  fHistTopDCAPosToPV(0),
  fHistTopDCAV0Daughters(0),
  fHistTopCosinePA(0),
  fHistTopV0Radius(0),
  fHistSelectedTopDCANegToPV(0),
  fHistSelectedTopDCAPosToPV(0),
  fHistSelectedTopDCAV0Daughters(0),
  fHistSelectedTopCosinePA(0),
  fHistSelectedTopV0Radius(0),

  f2dHistInvMassK0Short(0),
  f2dHistInvMassLambda(0),
  f2dHistInvMassAntiLambda(0),

  f2dHistInvMassWithdEdxK0Short(0),
  f2dHistInvMassWithdEdxLambda(0),
  f2dHistInvMassWithdEdxAntiLambda(0),

  f2dHistInvMassWithdEdxCowboysK0Short(0),
  f2dHistInvMassWithdEdxCowboysLambda(0),
  f2dHistInvMassWithdEdxCowboysAntiLambda(0),

  f2dHistInvMassWithdEdxSailorsK0Short(0),
  f2dHistInvMassWithdEdxSailorsLambda(0),
  f2dHistInvMassWithdEdxSailorsAntiLambda(0),

  f2dHistResponseNegativeAsPion(0),
  f2dHistResponseNegativeAsProton(0),
  f2dHistResponsePositiveAsPion(0),
  f2dHistResponsePositiveAsProton(0),

  f2dHistdEdxSignalPionFromLambda(0),
  f2dHistdEdxSignalProtonFromLambda(0),
  f2dHistResponsePionFromLambda(0),
  f2dHistResponseProtonFromLambda(0),

fTrigType(AliVEvent::kINT7),

  //Task Control / Utils
  fPIDResponse(0),
  fkRunV0Vertexer ( kFALSE ),
  fkUseLightV0Vertexer ( kFALSE ),
  fdEdxCut (3) 
{
  // Constructor
  //VERTEXER CUTS
  // REALLY LOOSE? Be careful when attempting to run over PbPb if fkRunV0Vertexer is set! 
  fV0VertexerSels[0] =  33.  ;  // max allowed chi2
  fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
  fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
  fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
  fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
  fV0VertexerSels[5] =   0.5 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
  fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)

  //SELECTION CUTS
  // REALLY LOOSE? Be careful when attempting to run over PbPb if fkRunV0Vertexer is set! 
  fV0Sels[0] =  33.  ;  // max allowed chi2
  fV0Sels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
  fV0Sels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
  fV0Sels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
  fV0Sels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
  fV0Sels[5] =   0.5 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
  fV0Sels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)

  // Output slot #0 writes into a TList container (Lambda Histos and fTree)
  DefineOutput(1, TList::Class());
}


AliAnalysisTaskQAV0::~AliAnalysisTaskQAV0()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fOutput){
      delete fOutput;
      fOutput = 0x0;
   }
}

//________________________________________________________________________
void AliAnalysisTaskQAV0::UserCreateOutputObjects()
{
  //Define Output Lists
  fOutput = new TList();
  fOutput->SetOwner(); 

  //Histogram Output: Event-by-Event
  fHistEvent = new TH1D( "fHistEvent", ";Evt. Sel. Step;Count",4,0,4); 
  fHistEvent->GetXaxis()->SetBinLabel(1, "Processed");
  fHistEvent->GetXaxis()->SetBinLabel(2, "Phys-Sel");  
  fHistEvent->GetXaxis()->SetBinLabel(3, "Has Vtx");  
  fHistEvent->GetXaxis()->SetBinLabel(4, "Vtx |z|<10cm");  
  fOutput->Add(fHistEvent); 

  //Topological Selection Histograms, 1D 
  fHistTopDCANegToPV        = new TH1D( "fHistTopDCANegToPV",";DCA Neg. Daughter to PV (cm);Counts",200,0,1); 
  fHistTopDCAPosToPV        = new TH1D( "fHistTopDCAPosToPV",";DCA Pos. Daughter to PV (cm);Counts",200,0,1); 
  fHistTopDCAV0Daughters    = new TH1D( "fHistTopDCAV0Daughters",";DCA V0 Daughters (#sigma);Counts",200,0,2); 
  fHistTopCosinePA          = new TH1D( "fHistTopCosinePA",";Cosine of PA;Counts",10010,-1.001,1.001); 
  fHistTopV0Radius          = new TH1D( "fHistTopV0Radius",";Decay Radius (cm);Counts",200,0.,10);

  fOutput->Add( fHistTopDCANegToPV     );
  fOutput->Add( fHistTopDCAPosToPV     );
  fOutput->Add( fHistTopDCAV0Daughters );
  fOutput->Add( fHistTopCosinePA       );
  fOutput->Add( fHistTopV0Radius       );

  //Zoomed In 
  fHistSelectedTopDCANegToPV        = new TH1D( "fHistSelectedTopDCANegToPV",";DCA Neg. Daughter to PV (cm);Counts",200,fV0Sels[1],1); 
  fHistSelectedTopDCAPosToPV        = new TH1D( "fHistSelectedTopDCAPosToPV",";DCA Pos. Daughter to PV (cm);Counts",200,fV0Sels[2],1); 
  fHistSelectedTopDCAV0Daughters    = new TH1D( "fHistSelectedTopDCAV0Daughters",";DCA V0 Daughters (#sigma);Counts",200,0,fV0Sels[3]); 
  fHistSelectedTopCosinePA          = new TH1D( "fHistSelectedTopCosinePA",";Cosine of PA;Counts",400,fV0Sels[4],1); 
  fHistSelectedTopV0Radius          = new TH1D( "fHistSelectedTopV0Radius",";Decay Radius (cm);Counts",200,fV0Sels[5],10);

  fOutput->Add( fHistSelectedTopDCANegToPV     );
  fOutput->Add( fHistSelectedTopDCAPosToPV     );
  fOutput->Add( fHistSelectedTopDCAV0Daughters );
  fOutput->Add( fHistSelectedTopCosinePA       );
  fOutput->Add( fHistSelectedTopV0Radius       );

  //Invariant Mass Plots
  f2dHistInvMassK0Short     = new TH2D( "f2dHistInvMassK0Short"   , ";p_{T};M(#pi^{+},#pi^{-})"   ,250,0,25,500,0.25,0.75);
  f2dHistInvMassLambda      = new TH2D( "f2dHistInvMassLambda"    , ";p_{T};M(p,#pi^{-})"         ,250,0,25,300,1.07,1.115+0.255);
  f2dHistInvMassAntiLambda  = new TH2D( "f2dHistInvMassAntiLambda", ";p_{T};M(#pi^{+},#bar{p})"   ,250,0,25,300,1.07,1.115+0.255);

  fOutput->Add( f2dHistInvMassK0Short     );
  fOutput->Add( f2dHistInvMassLambda      );
  fOutput->Add( f2dHistInvMassAntiLambda  );

  //Invariant Mass Plots, with dEdx
  f2dHistInvMassWithdEdxK0Short     = new TH2D( "f2dHistInvMassWithdEdxK0Short"   , ";p_{T};M(#pi^{+},#pi^{-})"   ,250,0,25,500,0.25,0.75);
  f2dHistInvMassWithdEdxLambda      = new TH2D( "f2dHistInvMassWithdEdxLambda"    , ";p_{T};M(p,#pi^{-})"         ,250,0,25,300,1.07,1.115+0.255);
  f2dHistInvMassWithdEdxAntiLambda  = new TH2D( "f2dHistInvMassWithdEdxAntiLambda", ";p_{T};M(#pi^{+},#bar{p})"   ,250,0,25,300,1.07,1.115+0.255);

  fOutput->Add( f2dHistInvMassWithdEdxK0Short     );
  fOutput->Add( f2dHistInvMassWithdEdxLambda      );
  fOutput->Add( f2dHistInvMassWithdEdxAntiLambda  );
    
    //Invariant Mass Plots, with dEdx + are cowboys
    f2dHistInvMassWithdEdxCowboysK0Short     = new TH2D( "f2dHistInvMassWithdEdxCowboysK0Short"   , ";p_{T};M(#pi^{+},#pi^{-})"   ,250,0,25,500,0.25,0.75);
    f2dHistInvMassWithdEdxCowboysLambda      = new TH2D( "f2dHistInvMassWithdEdxCowboysLambda"    , ";p_{T};M(p,#pi^{-})"         ,250,0,25,300,1.07,1.115+0.255);
    f2dHistInvMassWithdEdxCowboysAntiLambda  = new TH2D( "f2dHistInvMassWithdEdxCowboysAntiLambda", ";p_{T};M(#pi^{+},#bar{p})"   ,250,0,25,300,1.07,1.115+0.255);
    
    fOutput->Add( f2dHistInvMassWithdEdxCowboysK0Short     );
    fOutput->Add( f2dHistInvMassWithdEdxCowboysLambda      );
    fOutput->Add( f2dHistInvMassWithdEdxCowboysAntiLambda  );
    
    //Invariant Mass Plots, with dEdx + are sailors
    f2dHistInvMassWithdEdxSailorsK0Short     = new TH2D( "f2dHistInvMassWithdEdxSailorsK0Short"   , ";p_{T};M(#pi^{+},#pi^{-})"   ,250,0,25,500,0.25,0.75);
    f2dHistInvMassWithdEdxSailorsLambda      = new TH2D( "f2dHistInvMassWithdEdxSailorsLambda"    , ";p_{T};M(p,#pi^{-})"         ,250,0,25,300,1.07,1.115+0.255);
    f2dHistInvMassWithdEdxSailorsAntiLambda  = new TH2D( "f2dHistInvMassWithdEdxSailorsAntiLambda", ";p_{T};M(#pi^{+},#bar{p})"   ,250,0,25,300,1.07,1.115+0.255);
    
    fOutput->Add( f2dHistInvMassWithdEdxSailorsK0Short     );
    fOutput->Add( f2dHistInvMassWithdEdxSailorsLambda      );
    fOutput->Add( f2dHistInvMassWithdEdxSailorsAntiLambda  );

  //dE/dx QA for main analysis and for calibration check 
  f2dHistResponseNegativeAsPion    = new TH2D( "f2dHistResponseNegativeAsPion", ";p_{T}^{V0};N#sigma",500,0,5,400,-20,20);
  f2dHistResponseNegativeAsProton  = new TH2D( "f2dHistResponseNegativeAsProton", ";p_{T}^{V0};N#sigma",500,0,5,400,-20,20);
  f2dHistResponsePositiveAsPion    = new TH2D( "f2dHistResponsePositiveAsPion", ";p_{T}^{V0};N#sigma",500,0,5,400,-20,20);
  f2dHistResponsePositiveAsProton  = new TH2D( "f2dHistResponsePositiveAsProton", ";p_{T}^{V0};N#sigma",500,0,5,400,-20,20);

  //Clean Signal Check from Lambdas: stricter cuts, raw signal check 
  f2dHistdEdxSignalPionFromLambda    = new TH2D( "f2dHistdEdxSignalPionFromLambda", ";p_{T}^{V0};TPC Signal",500,0,5,8000,0,800);
  f2dHistdEdxSignalProtonFromLambda  = new TH2D( "f2dHistdEdxSignalProtonFromLambda", ";p_{T}^{V0};TPC Signal",500,0,5,8000,0,800);
  f2dHistResponsePionFromLambda      = new TH2D( "f2dHistResponsePionFromLambda", ";p_{T}^{V0};N#sigma",500,0,5,400,-20,20);
  f2dHistResponseProtonFromLambda    = new TH2D( "f2dHistResponseProtonFromLambda", ";p_{T}^{V0};N#sigma",500,0,5,400,-20,20);


  fOutput->Add( f2dHistResponseNegativeAsPion        );
  fOutput->Add( f2dHistResponseNegativeAsProton      );
  fOutput->Add( f2dHistResponsePositiveAsPion        );
  fOutput->Add( f2dHistResponsePositiveAsProton      );

  fOutput->Add( f2dHistdEdxSignalPionFromLambda        );
  fOutput->Add( f2dHistdEdxSignalProtonFromLambda      );
  fOutput->Add( f2dHistResponsePionFromLambda          );
  fOutput->Add( f2dHistResponseProtonFromLambda        );

//------------------------------------------------
// Particle Identification Setup
//------------------------------------------------

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
   //Regular output: Histograms
   PostData(1, fOutput);
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskQAV0::UserExec(Option_t *) 
{
   // Main loop
   // Called for each event
   //gObjectTable->Print();
   AliESDEvent *lESDevent = 0x0;

   //AliAODEvent *lAODevent = 0x0;
   Int_t    nV0s                        = -1;

   Double_t lTrkgPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lMagneticField                 = -10.;
	
   // Connect to the InputEvent	
   // After these lines, we should have an ESD/AOD event + the number of cascades in it.

   // Appropriate for ESD analysis! 
   lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
   if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
   }
    
    //------------------------------------------------
    // Rerun V0 vertexer, if asked for
    // --- WARNING: Be careful when using in PbPb
    //------------------------------------------------
    if( fkRunV0Vertexer ){
        if(!fkUseLightV0Vertexer){
            lESDevent->ResetV0s();
            AliV0vertexer lV0vtxer;
            lV0vtxer.SetCuts(fV0VertexerSels);
            lV0vtxer.Tracks2V0vertices(lESDevent);
        }
        if(fkUseLightV0Vertexer){
            lESDevent->ResetV0s();
            AliLightV0vertexer lV0vtxer;
            lV0vtxer.SetCuts(fV0VertexerSels);
            lV0vtxer.Tracks2V0vertices(lESDevent);
        }
    }
    
    fHistEvent->Fill(0.5);

//------------------------------------------------
// Physics Selection
//------------------------------------------------
  
  //Standard Min-Bias Selection
    if ( ! AliMultSelectionTask::IsSelectedTrigger(lESDevent, fTrigType) ) {
    PostData(1, fOutput);
    return;
  }
  fHistEvent->Fill(1.5); 
  
//------------------------------------------------
// After Trigger Selection
//------------------------------------------------


	nV0s = lESDevent->GetNumberOfV0s();
  
//------------------------------------------------
// Getting: Primary Vertex + MagField Info
//------------------------------------------------

	const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
	// get the vtx stored in ESD found with tracks
	lPrimaryTrackingESDVtx->GetXYZ( lTrkgPrimaryVtxPos );
        
	const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();	
  // get the best primary vertex available for the event
	// As done in AliCascadeVertexer, we keep the one which is the best one available.
	// between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
	// This one will be used for next calculations (DCA essentially)
   lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

   Double_t tPrimaryVtxPosition[3];
   const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
   tPrimaryVtxPosition[0] = primaryVtx->GetX();
   tPrimaryVtxPosition[1] = primaryVtx->GetY();
   tPrimaryVtxPosition[2] = primaryVtx->GetZ();

  //------------------------------------------------
  // Primary Vertex Requirements Section:
  //  ---> pp and PbPb: Only requires |z|<10cm
  //  ---> pPb: all requirements checked at this stage
  //------------------------------------------------
  
  //Roberto's PV selection criteria, implemented 17th April 2013
  
  // vertex selection 
  Bool_t fHasVertex = kFALSE;
  const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    vertex = lESDevent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
    else fHasVertex = kTRUE;
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vertex->IsFromVertexerZ() && (zRes>0.25)) fHasVertex = kFALSE;
  }
  else fHasVertex = kTRUE;
  
  //Is First event in chunk rejection: Still present!
  if(fHasVertex == kFALSE) {
    AliWarning("Pb / | PV does not satisfy selection criteria!");
    PostData(1, fOutput);
    return;
  }
  
  fHistEvent->Fill(2.5); 

  //17 April Fix: Always do primary vertex Z selection, after pA vertex selection from Roberto
  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0) {
    AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
    PostData(1, fOutput);
    return;
  }

  fHistEvent->Fill(3.5); 
  
  lMagneticField = lESDevent->GetMagneticField( );

//------------------------------------------------
// Only look at events with well-established PV
//------------------------------------------------
	
   //const AliESDVertex *lPrimaryTrackingESDVtxCheck = lESDevent->GetPrimaryVertexTracks();
   const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
   //if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtxCheck->GetStatus() ){
   //   AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
   //   PostData(1, fOutput);
   //   return;
   //}
  


//------------------------------------------------
// MAIN LAMBDA LOOP STARTS HERE
//------------------------------------------------

//Variable definition
  Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
  Double_t lChi2V0 = 0;
  Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
  Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
  Double_t lV0CosineOfPointingAngle = 0;
  Double_t lV0Radius = 0, lPt = 0;
  Double_t lRapK0Short = 0, lRapLambda = 0;
  Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lAlphaV0 = 0, lPtArmV0 = 0;
  Double_t lNegEta = -100, lPosEta = -100;
  Double_t fMinV0Pt = 0; 
  Double_t fMaxV0Pt = 100; 

   Int_t nv0s = 0;
   nv0s = lESDevent->GetNumberOfV0s();

   //for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
   for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
   {// This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
      if (!v0) continue;

      //Only use Offline Candidates for QA       
      lOnFlyStatus = v0->GetOnFlyStatus();
      if( lOnFlyStatus == kTRUE ) continue; 

      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 

      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] ); 
      Double_t lV0TotalMomentum = TMath::Sqrt(
      tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

      lPt = v0->Pt();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();
      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

      UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
      UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

      Double_t lMomPos[3]; v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3]; v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

       //Provisions for cowboy/sailor check
       Double_t lModp1 = TMath::Sqrt( lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1] );
       Double_t lModp2 = TMath::Sqrt( lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1] );
       
       //Calculate vec prod with momenta projected to xy plane
       Double_t lVecProd = (lMomPos[0]*lMomNeg[1] - lMomPos[1]*lMomNeg[0]) / (lModp1*lModp2);
       
       if ( lMagneticField < 0 ) lVecProd *= -1; //invert sign

       Bool_t lIsCowboy = kFALSE;
       if (lVecProd < 0) lIsCowboy = kTRUE;
       
      AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
      AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
         Printf("ERROR: Could not retreive one of the daughter track");
         continue;
      }

      //Daughter Eta for Eta selection, afterwards
      lNegEta = nTrack->Eta();
      lPosEta = pTrack->Eta();

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->GetSign() == nTrack->GetSign()){
         continue;
      } 

      //________________________________________________________________________
      // Track quality cuts 
      Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
      Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
      Int_t lLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
      if( lNegTrackCrossedRows < lLeastNbrCrossedRows )
         lLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
     
      //Get status flags
      //lPosTrackStatus = pTrack->GetStatus();
      //lNegTrackStatus = nTrack->GetStatus();
     
      if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
	
      //GetKinkIndex condition
      if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

      //Findable clusters > 0 condition
      if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;

      //Compute ratio Crossed Rows / Findable clusters
      //Note: above test avoids division by zero! 
      Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF())); 
      Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF())); 

      Float_t lLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
      if( lNegTrackCrossedRowsOverFindable < lLeastRatioCrossedRowsOverFindable )
         lLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

      //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
      if ( lLeastRatioCrossedRowsOverFindable < 0.8 ) continue;

      //End track Quality Cuts
      //________________________________________________________________________

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(tPrimaryVtxPosition[0],
							tPrimaryVtxPosition[1],
							lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(tPrimaryVtxPosition[0],
							tPrimaryVtxPosition[1],
							lMagneticField) );


      lChi2V0 = v0->GetChi2V0();
      lDcaV0Daughters = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
      lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
      //fTreeVariableV0CosineOfPointingAngle=lV0CosineOfPointingAngle;

      // Getting invariant mass infos directly from ESD
      v0->ChangeMassHypothesis(310);
      lInvMassK0s = v0->GetEffMass();
      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();
      lAlphaV0 = v0->AlphaV0();
      lPtArmV0 = v0->PtArmV0();
      
      //Official means of acquiring N-sigmas 
      Float_t lNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
      Float_t lNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
      Float_t lNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
      Float_t lNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );

//This requires an Invariant Mass Hypothesis afterwards
      Float_t lDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
					);
      lDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

//------------------------------------------------
// Fill Main Output Histograms
//------------------------------------------------

  //Topological Variable Checks, One-Dimensional      
  fHistTopDCANegToPV      -> Fill( lDcaNegToPrimVertex      ) ; 
  fHistTopDCAPosToPV      -> Fill( lDcaPosToPrimVertex      ) ; 
  fHistTopDCAV0Daughters  -> Fill( lDcaV0Daughters          ) ; 
  fHistTopCosinePA        -> Fill( lV0CosineOfPointingAngle ) ; 
  fHistTopV0Radius        -> Fill( lV0Radius                ) ; 


  if( lDcaNegToPrimVertex > fV0Sels[1] && lDcaPosToPrimVertex > fV0Sels[2]      && 
      lDcaV0Daughters     < fV0Sels[3] && lV0CosineOfPointingAngle > fV0Sels[4] &&
      lV0Radius           > fV0Sels[5] && lV0Radius < fV0Sels [6] ){ 

    //Topological Variables zoomed in at selection level (whatever that may be)
    //May be slightly redundant if no specific extra configuration was done 
    fHistSelectedTopDCANegToPV      -> Fill( lDcaNegToPrimVertex      ) ; 
    fHistSelectedTopDCAPosToPV      -> Fill( lDcaPosToPrimVertex      ) ; 
    fHistSelectedTopDCAV0Daughters  -> Fill( lDcaV0Daughters          ) ; 
    fHistSelectedTopCosinePA        -> Fill( lV0CosineOfPointingAngle ) ; 
    fHistSelectedTopV0Radius        -> Fill( lV0Radius                ) ; 

    //Specific fV0Sel selection level, but no dEdx applied 
    f2dHistInvMassK0Short     -> Fill ( lPt , lInvMassK0s        )   ; 
    f2dHistInvMassLambda      -> Fill ( lPt , lInvMassLambda     )   ;
    f2dHistInvMassAntiLambda  -> Fill ( lPt , lInvMassAntiLambda )   ;

    //General dE/dx QA 
    f2dHistResponseNegativeAsPion     -> Fill( lPt, lNSigmasNegPion      ); 
    f2dHistResponseNegativeAsProton   -> Fill( lPt, lNSigmasNegProton    ); 
    f2dHistResponsePositiveAsPion     -> Fill( lPt, lNSigmasPosPion      ); 
    f2dHistResponsePositiveAsProton   -> Fill( lPt, lNSigmasPosProton    ); 

    //Clean Sample From Lambdas
    //Very strict cuts to ensure dealing with good Lambdas 
    if ( lDcaV0Daughters < 1.0 && lV0CosineOfPointingAngle > 0.999 && TMath::Abs( lInvMassK0s - 0.497614 ) > 0.012 
          && TMath::Abs( lInvMassAntiLambda - 1.115683) > 0.08 && TMath::Abs( lInvMassLambda - 1.115683) < 0.002 ) { 
      
      f2dHistdEdxSignalPionFromLambda     -> Fill( lPt, nTrack-> GetTPCsignal() );
      f2dHistdEdxSignalProtonFromLambda   -> Fill( lPt, pTrack-> GetTPCsignal() );
      f2dHistResponsePionFromLambda     -> Fill( lPt, lNSigmasNegPion   );
      f2dHistResponseProtonFromLambda   -> Fill( lPt, lNSigmasPosProton );      
    }

    //Specific fV0Sel selection level, dE/dx applied
      if ( TMath::Abs(lNSigmasPosPion)   < fdEdxCut && TMath::Abs(lNSigmasNegPion)   < fdEdxCut ) {
          f2dHistInvMassWithdEdxK0Short       -> Fill ( lPt , lInvMassK0s        )   ;
          if (lIsCowboy) f2dHistInvMassWithdEdxCowboysK0Short       -> Fill ( lPt , lInvMassK0s        )   ;
          if (!lIsCowboy) f2dHistInvMassWithdEdxSailorsK0Short       -> Fill ( lPt , lInvMassK0s        )   ;
      }
      if ( TMath::Abs(lNSigmasPosProton) < fdEdxCut && TMath::Abs(lNSigmasNegPion)   < fdEdxCut ) {
          f2dHistInvMassWithdEdxLambda        -> Fill ( lPt , lInvMassLambda     )   ;
          if (lIsCowboy) f2dHistInvMassWithdEdxCowboysLambda        -> Fill ( lPt , lInvMassLambda     )   ;
          if (!lIsCowboy) f2dHistInvMassWithdEdxSailorsLambda        -> Fill ( lPt , lInvMassLambda     )   ;
      }
      if ( TMath::Abs(lNSigmasPosPion)   < fdEdxCut && TMath::Abs(lNSigmasNegProton) < fdEdxCut ) {
          f2dHistInvMassWithdEdxAntiLambda    -> Fill ( lPt , lInvMassAntiLambda )   ;
          if (lIsCowboy) f2dHistInvMassWithdEdxCowboysAntiLambda    -> Fill ( lPt , lInvMassAntiLambda )   ;
          if (!lIsCowboy) f2dHistInvMassWithdEdxSailorsAntiLambda    -> Fill ( lPt , lInvMassAntiLambda )   ;
      }

  }


  //Topological Variable Checks, Two-Dimensional 

//------------------------------------------------
// End Filling of main Histograms
//------------------------------------------------

  }// This is the end of the V0 loop


  // Post output data.
  PostData(1, fOutput);

}

//________________________________________________________________________
void AliAnalysisTaskQAV0::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query
   // This will draw the V0 candidate multiplicity, whose 
   // number of entries corresponds to the number of triggered events. 
   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskQAV0 : ouput data container list not available\n");
      return;
   }		
   fHistEvent = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEvent")  );
   if (!fHistEvent) {
      Printf("ERROR - AliAnalysisTaskQAV0 : fHistEvent not available");
      return;
   }
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskQAV0","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();
   fHistEvent->SetMarkerStyle(22);
   fHistEvent->DrawCopy("E");
}
