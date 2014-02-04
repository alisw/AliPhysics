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
//            AliQAProdMultistrange class
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
//                - Centrality selection added (+ setters) 
//                - setters added (vertex range)
//                - histo added and histo/container binning changed 
//                - protection in the destructor for CAF usage          
//                - AliWarning disabled
//                - automatic settings for PID
//               September2011
//                - proper time histos/container added (V0 and Cascades)
//               November2011
//                - re-run V0's and cascade's vertexers (SetCuts instead of SetDefaultCuts!!)
//               Genuary2012 
//                - AOD analysis part completed 
//               March2012
//                - min number of TPC clusters for track selection as a parameter       
//               July2012
//                - cut on min pt for daughter tracks as a parameter (+control histos)
//               August2012 
//                - cut on pseudorapidity for daughter tracks as a parameter (+control histos for Xi-)
//-----------------------------------------------------------------

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
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h" 
#include "AliAODInputHandler.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAODTrack.h"

#include "AliQAProdMultistrange.h"

ClassImp(AliQAProdMultistrange)



//________________________________________________________________________
AliQAProdMultistrange::AliQAProdMultistrange() 
  : AliAnalysisTaskSE(), 
    fAnalysisType               ("ESD"), 
    fESDtrackCuts               (0),
    fCollidingSystem            ("PbPb"),
    fPIDResponse                (0),
    fkSDDSelectionOn            (kTRUE),
    fkQualityCutZprimVtxPos     (kTRUE),
    fkQualityCutNoTPConlyPrimVtx(kTRUE),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCutnTPCcls         (kTRUE),
    fkQualityCutPileup          (kTRUE),
    fwithSDD                    (kTRUE),
    fMinnTPCcls                 (0),  
    fCentrLowLim                (0),
    fCentrUpLim                 (0),
    fCentrEstimator             (0),
    fkUseCleaning               (0),
    fVtxRange                   (0),
    fMinPtCutOnDaughterTracks   (0),
    fEtaCutOnDaughterTracks     (0),

    
    fCFContCascadeCuts(0)
    

{
  // Dummy Constructor
}


//________________________________________________________________________
AliQAProdMultistrange::AliQAProdMultistrange(const char *name) 
  : AliAnalysisTaskSE(name),
    fAnalysisType               ("ESD"), 
    fESDtrackCuts               (0),
    fCollidingSystem            ("PbPb"),
    fPIDResponse                (0),
    fkSDDSelectionOn            (kTRUE),
    fkQualityCutZprimVtxPos     (kTRUE),
    fkQualityCutNoTPConlyPrimVtx(kTRUE),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCutnTPCcls         (kTRUE),
    fkQualityCutPileup          (kTRUE),
    fwithSDD                    (kTRUE),
    fMinnTPCcls                 (0),
    fCentrLowLim                (0),
    fCentrUpLim                 (0),
    fCentrEstimator             (0),
    fkUseCleaning               (0),
    fVtxRange                   (0),
    fMinPtCutOnDaughterTracks   (0),
    fEtaCutOnDaughterTracks     (0),


    fCFContCascadeCuts(0)
    

{
  // Constructor
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, AliCFContainer::Class());

  AliLog::SetClassDebugLevel("AliQAProdMultistrange",1); // MN this should (?) enable only AliFatal
}


AliQAProdMultistrange::~AliQAProdMultistrange() {
  //
  // Destructor
  //
  // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
  // They will be deleted when fListCascade is deleted by the TSelector dtor
  // Because of TList::SetOwner() ...
  if (fCFContCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadeCuts;     fCFContCascadeCuts = 0x0;  }
  if (fESDtrackCuts)         { delete fESDtrackCuts;        fESDtrackCuts = 0x0; }
}



//________________________________________________________________________
void AliQAProdMultistrange::UserCreateOutputObjects() {
  // Create histograms
  // Called once

 //-----------------------------------------------
 // Particle Identification Setup (new PID object)
 //-----------------------------------------------
 AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
 AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
 fPIDResponse = inputHandler->GetPIDResponse();


 // Only used to get the number of primary reconstructed tracks
 if (fAnalysisType == "ESD" && (! fESDtrackCuts )){
   fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
   //Printf("CheckCascade - ESDtrackCuts set up to 2010 std ITS-TPC cuts...");
 }


 //---------------------------------------------------
 // Define the container for the topological variables
 //---------------------------------------------------
  if(! fCFContCascadeCuts) {
      // Container meant to store all the relevant distributions corresponding to the cut variables.
      // NB: overflow/underflow of variables on which we want to cut later should be 0!!! 
      const Int_t  lNbSteps      =  4 ;
      const Int_t  lNbVariables  =  21 ;
      //Array for the number of bins in each dimension :
      Int_t lNbBinsPerVar[lNbVariables] = {0};
      lNbBinsPerVar[0]  = 25;     //DcaCascDaughters             :  [0.0,2.4,3.0]       -> Rec.Cut = 2.0;
      lNbBinsPerVar[1]  = 25;     //DcaBachToPrimVertex          :  [0.0,0.24,100.0]    -> Rec.Cut = 0.01; 
      lNbBinsPerVar[2]  = 30;     //CascCosineOfPointingAngle    :  [0.97,1.0]          -> Rec.Cut = 0.98;
      lNbBinsPerVar[3]  = 40;     //CascRadius                   :  [0.0,3.9,1000.0]    -> Rec.Cut = 0.2;
      lNbBinsPerVar[4]  = 30;     //InvMassLambdaAsCascDghter    :  [1.1,1.3]           -> Rec.Cut = 0.008;
      lNbBinsPerVar[5]  = 20;     //DcaV0Daughters               :  [0.0,2.0]           -> Rec.Cut = 1.5;
      lNbBinsPerVar[6]  = 201;    //V0CosineOfPointingAngleToPV  :  [0.89,1.0]          -> Rec.Cut = 0.9;
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
      lNbBinsPerVar[18] = 201;    //V0CosineOfPointingAngleToXiV
      lNbBinsPerVar[19] = 11;     //Centrality
      lNbBinsPerVar[20] = 100;    //ESD track multiplicity
      //define the container
      fCFContCascadeCuts = new AliCFContainer("fCFContCascadeCuts","Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar );
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
      //6 - V0CosineOfPointingAngleToPV
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
      fCFContCascadeCuts -> SetBinLimits(10,  lBinLim10 );     
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
      //18 - V0CosineOfPointingAngleToXiV
      fCFContCascadeCuts -> SetBinLimits(18, 0.8, 1.001);
      //19
      Double_t *lBinLim19  = new Double_t[ lNbBinsPerVar[19]+1 ];
         for(Int_t i=3; i< lNbBinsPerVar[19]+1;i++)   lBinLim19[i]  = (Double_t)(i-1)*10.;
         lBinLim19[0] = 0.0; 
         lBinLim19[1] = 5.0; 
         lBinLim19[2] = 10.0;
      fCFContCascadeCuts->SetBinLimits(19,  lBinLim19 );     
      delete [] lBinLim19;
      //20
      fCFContCascadeCuts->SetBinLimits(20, 0.0, 6000.0);
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
      fCFContCascadeCuts->SetVarTitle(6,  "cos(V0 PA) in cascade to PV");
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
      fCFContCascadeCuts->SetVarTitle(18,  "cos(V0 PA) in cascade to XiV");
      fCFContCascadeCuts->SetVarTitle(19, "Centrality");
      fCFContCascadeCuts->SetVarTitle(20, "ESD track multiplicity");
  }

PostData(1, fCFContCascadeCuts);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliQAProdMultistrange::UserExec(Option_t *) {

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Main loop (called for each event)
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //-----------------------
  //Define ESD/AOD handlers 
  AliESDEvent *lESDevent = 0x0;
  AliAODEvent *lAODevent = 0x0;

  //---------------------
  //Check the PIDresponse
  if(!fPIDResponse) {
       AliError("Cannot get pid response");
       return;
  }

  //__________________________________________________
  // After these lines we should have an ESD/AOD event

  //---------------------------------------------------------
  //Load the InputEvent and check, before any event selection	
  //---------------------------------------------------------
  Float_t  lPrimaryTrackMultiplicity = -1.;
  AliCentrality* centrality = 0;
  if (fAnalysisType == "ESD") {
      lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
      if (!lESDevent) {
          AliWarning("ERROR: lESDevent not available \n");
          return;
      }
      if (fCollidingSystem == "PbPb") lPrimaryTrackMultiplicity = fESDtrackCuts->CountAcceptedTracks(lESDevent);
      if (fCollidingSystem == "PbPb") centrality = lESDevent->GetCentrality();
      
  } else if (fAnalysisType == "AOD") {
      lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
      if (!lAODevent) {
          AliWarning("ERROR: lAODevent not available \n");
          return;
      }
      if (fCollidingSystem == "PbPb") {
          lPrimaryTrackMultiplicity = 0;
          Int_t    nTrackMultiplicity = (InputEvent())->GetNumberOfTracks();
          for (Int_t itrack = 0; itrack < nTrackMultiplicity; itrack++) {
               AliAODTrack* track = lAODevent->GetTrack(itrack);
               if (track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) lPrimaryTrackMultiplicity++; 
          }
      }
      if (fCollidingSystem == "PbPb") centrality = lAODevent->GetCentrality();
  } else {
    Printf("Analysis type (ESD or AOD) not specified \n");
    return;
  }

  //-----------------------------------------
  // Centrality selection for PbPb collisions
  //-----------------------------------------
  Float_t lcentrality = 0.;
  if (fCollidingSystem == "PbPb") { 
       if (fkUseCleaning) lcentrality = centrality->GetCentralityPercentile(fCentrEstimator.Data());
       else {
           lcentrality = centrality->GetCentralityPercentileUnchecked(fCentrEstimator.Data());
           if (centrality->GetQuality()>1) {
               PostData(1, fCFContCascadeCuts);
               return;
           }
       }
       if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim) { 
           PostData(1, fCFContCascadeCuts);
           return;
       }
  } else if (fCollidingSystem == "pp") lcentrality = 0.;


  //----------------------------------------
  // SDD selection for pp@2.76TeV collisions
  //----------------------------------------
  if (fCollidingSystem == "pp") {
      if (fAnalysisType == "ESD") {
          if (fkSDDSelectionOn) {
              TString trcl = lESDevent->GetFiredTriggerClasses();
              if      (fwithSDD) { if(!(trcl.Contains("ALLNOTRD"))) { PostData(1, fCFContCascadeCuts); return; } }
              else if (!fwithSDD){ if((trcl.Contains("ALLNOTRD")))  { PostData(1, fCFContCascadeCuts); return; } }
          }
      } else if (fAnalysisType == "AOD") {
          if (fkSDDSelectionOn) {
              TString trcl = lAODevent->GetFiredTriggerClasses();
              if      (fwithSDD)  { if(!(trcl.Contains("ALLNOTRD"))) { PostData(1, fCFContCascadeCuts); return; } }
              else if (!fwithSDD) { if((trcl.Contains("ALLNOTRD")))  { PostData(1, fCFContCascadeCuts); return; } }
          }
      }
  }

  //--------------------------------------------
  // Physics selection for pp@2.76TeV collisions
  //--------------------------------------------
  // - moved to the runGrid for the PbPb collisions
  if (fCollidingSystem == "pp") {
      if (fAnalysisType == "ESD") {
          UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
          Bool_t isSelected = 0;
          isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
          if(! isSelected){ PostData(1, fCFContCascadeCuts); return; }
      } else if (fAnalysisType == "AOD") {
          UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
          Bool_t isSelected = 0;
          isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
          if(! isSelected){ PostData(1, fCFContCascadeCuts); return; }
      }
  }

  //------------------------------
  // Well-established PV selection
  //------------------------------
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0}; 
  Double_t lMagneticField        = -10.;
  if (fAnalysisType == "ESD") {
       const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();   
       const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();	
       if (fkQualityCutNoTPConlyPrimVtx) {
           const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
           if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus() ){
               AliWarning(" No SPD prim. vertex nor prim. Tracking vertex ... return !");
               PostData(1, fCFContCascadeCuts);
               return;
           }
       }
       lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
       lMagneticField = lESDevent->GetMagneticField( );
  } else if (fAnalysisType == "AOD") {
       const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
       if (!lPrimaryBestAODVtx){
           AliWarning("No prim. vertex in AOD... return!");
           PostData(1, fCFContCascadeCuts);
           return;
       }
       lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
       lMagneticField = lAODevent->GetMagneticField();  
  }

  //------------------------------------------
  // Pilup selection for pp@2.76TeV collisions
  //------------------------------------------
  if (fCollidingSystem == "pp") { 
      if (fAnalysisType == "ESD") {
          if (fkQualityCutPileup) { if(lESDevent->IsPileupFromSPD()){ PostData(1, fCFContCascadeCuts); return; } }
      } else if (fAnalysisType == "AOD") {
          if (fkQualityCutPileup) { if(lAODevent->IsPileupFromSPD()){ PostData(1, fCFContCascadeCuts); return; } }
      }
  }

  //----------------------------
  // Vertex Z position selection
  //----------------------------
  if (fkQualityCutZprimVtxPos) {
      if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange ) {
          PostData(1, fCFContCascadeCuts);
          return;
      }
  }



  //////////////////////////////
  // CASCADE RECONSTRUCTION PART
  //////////////////////////////

  //%%%%%%%%%%%%%
  // Cascade loop
  Int_t ncascades = 0;
  if      (fAnalysisType == "ESD") ncascades = lESDevent->GetNumberOfCascades();
  else if (fAnalysisType == "AOD") ncascades = lAODevent->GetNumberOfCascades();

  for (Int_t iXi = 0; iXi < ncascades; iXi++) {// This is the begining of the Cascade loop (ESD or AOD)
	   
    // -------------------------------------
    // - Initialisation of the local variables that will be needed for ESD/AOD
    // -- Container variables (1st round)
    Double_t lDcaXiDaughters              = -1. ;                   //[Container]
    Double_t lXiCosineOfPointingAngle     = -1. ;                   //[Container]
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };             //Useful to define other variables: radius fid. vol., ctau, etc. for cascade
    Double_t lXiRadius                    = -1000. ;                //[Container]
    UShort_t lPosTPCClusters              = -1;                     //To check the quality of the tracks. For ESD only ...
    UShort_t lNegTPCClusters              = -1;                     //To check the quality of the tracks. For ESD only ...
    UShort_t lBachTPCClusters             = -1;                     //To check the quality of the tracks. For ESD only ...
    Double_t lInvMassLambdaAsCascDghter   = 0.;                     //[Container]
    Double_t lDcaV0DaughtersXi            = -1.;                    //[Container]
    Double_t lDcaBachToPrimVertexXi       = -1.;                    //[Container]
    Double_t lDcaV0ToPrimVertexXi         = -1.;                    //[Container]
    Double_t lDcaPosToPrimVertexXi        = -1.;                    //[Container]
    Double_t lDcaNegToPrimVertexXi        = -1.;                    //[Container]
    Double_t lV0CosineOfPointingAngle     = -1.;                    //[Container]
    Double_t lV0toXiCosineOfPointingAngle = -1.;                    //[Container] 
    Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. };             //Useful to define other variables: radius fid. vol., ctau, etc. for VO 
    Double_t lV0RadiusXi                  = -1000.0;                //[Container]
    Double_t lV0quality                   = 0.;                     //  ??
    Double_t lInvMassXiMinus              = 0.;                     //[Container]
    Double_t lInvMassXiPlus               = 0.;                     //[Container]
    Double_t lInvMassOmegaMinus           = 0.;                     //[Container]
    Double_t lInvMassOmegaPlus            = 0.;                     //[Container]
    // -- PID treatment
    Bool_t   lIsBachelorKaonForTPC = kFALSE; 
    Bool_t   lIsBachelorPionForTPC = kFALSE; 
    Bool_t   lIsNegPionForTPC      = kFALSE; 
    Bool_t   lIsPosPionForTPC      = kFALSE; 
    Bool_t   lIsNegProtonForTPC    = kFALSE; 
    Bool_t   lIsPosProtonForTPC    = kFALSE; 
    // -- More container variables and quality checks
    Double_t lXiMomX           = 0.;                               //Useful to define other variables: lXiTransvMom, lXiTotMom
    Double_t lXiMomY           = 0.;                               //Useful to define other variables: lXiTransvMom, lXiTotMom
    Double_t lXiMomZ           = 0.;                               //Useful to define other variables: lXiTransvMom, lXiTotMom
    Double_t lXiTransvMom      = 0.;                               //[Container]
    Double_t lXiTotMom         = 0.;                               //Useful to define other variables: cTau
    Double_t lV0PMomX          = 0.;                               //Useful to define other variables: lV0TotMom, lpTrackTransvMom
    Double_t lV0PMomY          = 0.;                               //Useful to define other variables: lV0TotMom, lpTrackTransvMom
    Double_t lV0PMomZ          = 0.;                               //Useful to define other variables: lV0TotMom, lpTrackTransvMom
    Double_t lV0NMomX          = 0.;                               //Useful to define other variables: lV0TotMom, lnTrackTransvMom
    Double_t lV0NMomY          = 0.;                               //Useful to define other variables: lV0TotMom, lnTrackTransvMom
    Double_t lV0NMomZ          = 0.;                               //Useful to define other variables: lV0TotMom, lnTrackTransvMom
    Double_t lV0TotMom         = 0.;                               //Useful to define other variables: lctauV0
    Double_t lBachMomX         = 0.;                               //Useful to define other variables: lBachTransvMom
    Double_t lBachMomY         = 0.;                               //Useful to define other variables: lBachTransvMom
    Double_t lBachMomZ         = 0.;                               //Useful to define other variables: lBachTransvMom
    Double_t lBachTransvMom    = 0.;                               //Selection on the min bachelor pT
    Double_t lpTrackTransvMom  = 0.;                               //Selection on the min bachelor pT
    Double_t lnTrackTransvMom  = 0.;                               //Selection on the min bachelor pT
    Short_t  lChargeXi         = -2;                               //Useful to select the particles based on the charge
    Double_t lRapXi            = -20.0;                            //[Container]
    Double_t lRapOmega         = -20.0;                            //[Container]
    Float_t  etaBach           = 0.;                               //Selection on the eta range
    Float_t  etaPos            = 0.;                               //Selection on the eta range
    Float_t  etaNeg            = 0.;                               //Selection on the eta range
    // --  variables for the AliCFContainer dedicated to cascade cut optmisiation: ESD and AOD 
    if (fAnalysisType == "ESD") { 
  
          // -------------------------------------
          // - Load the cascades from the handler 
          AliESDcascade *xi = lESDevent->GetCascade(iXi);
          if (!xi) continue;
        
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
          if (lBachIdx == lIdxNegXi) continue;    
          if (lBachIdx == lIdxPosXi) continue; 
          // - Get the track for the daughters
          AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
          AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
          AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
          if (!pTrackXi || !nTrackXi || !bachTrackXi )  continue;
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
           lPosTPCClusters   = pTrackXi->GetTPCNcls(); // FIXME: Is this ok? or something like in LambdaK0PbPb task AOD?
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
           if ( lChargeXi < 0 )		lInvMassXiMinus 	= xi->MassXi();
           if ( lChargeXi > 0 )		lInvMassXiPlus 		= xi->MassXi();
           if ( lChargeXi < 0 )		lInvMassOmegaMinus 	= xi->MassOmega();
           if ( lChargeXi > 0 )		lInvMassOmegaPlus 	= xi->MassOmega();

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

           //---------------------------------
           // - Extra info for QA (AOD)
           // Miscellaneous pieces of info that may help regarding data quality assessment.
           // Cascade transverse and total momentum     
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

    //---------------------------------------
    // Cut on pt of the three daughter tracks
    if (lBachTransvMom<fMinPtCutOnDaughterTracks) continue;
    if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
    if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) continue;
 
    //---------------------------------------------------
    // Cut on pseudorapidity of the three daughter tracks
    if (TMath::Abs(etaBach)>fEtaCutOnDaughterTracks) continue;
    if (TMath::Abs(etaPos)>fEtaCutOnDaughterTracks) continue;
    if (TMath::Abs(etaNeg)>fEtaCutOnDaughterTracks) continue;

    //----------------------------------
    // Calculate proper time for cascade
    Double_t cascadeMass = 0.;
    if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  ) cascadeMass = 1.321;
    if ( ( (lChargeXi<0) && lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )  ) cascadeMass = 1.672; 
    Double_t lctau =  TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power(( lPosXi[2]-lBestPrimaryVtxPos[2]),2));
    if (lXiTotMom!=0)         lctau = lctau*cascadeMass/lXiTotMom;
    else lctau = -1.;
    // Calculate proper time for Lambda (reconstructed)
    Float_t lambdaMass = 1.115683; // PDG mass
    Float_t distV0Xi =  TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2));
    Float_t lctauV0 = -1.;
    if (lV0TotMom!=0) lctauV0 = distV0Xi*lambdaMass/lV0TotMom;

        	
    // Fill the AliCFContainer (optimisation of topological selections)
    Double_t lContainerCutVars[21] = {0.0};
    lContainerCutVars[0]  = lDcaXiDaughters;
    lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
    lContainerCutVars[2]  = lXiCosineOfPointingAngle;
    lContainerCutVars[3]  = lXiRadius;
    lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
    lContainerCutVars[5]  = lDcaV0DaughtersXi;
    lContainerCutVars[6]  = lV0CosineOfPointingAngle;
    lContainerCutVars[7]  = lV0RadiusXi;
    lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;
    lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
    lContainerCutVars[10] = lDcaNegToPrimVertexXi;
    lContainerCutVars[13] = lXiTransvMom;
    lContainerCutVars[16] = lctau;
    lContainerCutVars[17] = lctauV0;
    lContainerCutVars[18] = lV0toXiCosineOfPointingAngle;
    lContainerCutVars[19] = lcentrality;
    lContainerCutVars[20] = lPrimaryTrackMultiplicity;
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


  }// end of the Cascade loop (ESD or AOD)
    
  
  // Post output data.
  PostData(1, fCFContCascadeCuts); 
}

//________________________________________________________________________
void AliQAProdMultistrange::Terminate(Option_t *) 
{

}
