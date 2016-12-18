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
#include "AliCFContainer.h"
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
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      fHistCascadeMultiplicityXiMinus(0), fHistCascadeMultiplicityXiPlus(0), fHistCascadeMultiplicityOmegaMinus(0), fHistCascadeMultiplicityOmegaPlus(0),
      fCFContCascadeCuts(0),
      fCFContCascadeMCgen(0)

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
      fHistMassXiMinus(0), fHistMassXiPlus(0), fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
      fHistCascadeMultiplicityXiMinus(0), fHistCascadeMultiplicityXiPlus(0), fHistCascadeMultiplicityOmegaMinus(0), fHistCascadeMultiplicityOmegaPlus(0),
      fCFContCascadeCuts(0),
      fCFContCascadeMCgen(0) 

{
  // Constructor
  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliCFContainer::Class());
  DefineOutput(3, AliCFContainer::Class());

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
  if (fCFContCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())      { delete fCFContCascadeCuts;      fCFContCascadeCuts      = 0x0; }
  if (fCFContCascadeMCgen && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())     { delete fCFContCascadeMCgen;     fCFContCascadeMCgen     = 0x0; }
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
  if(! fHistEventSel) {
        fHistEventSel = new TH1F("fHistEventSel", "Event selection;Evt. Sel. Step;Count",2, 0, 2);
        fHistEventSel->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventSel->GetXaxis()->SetBinLabel(2, "Selected");
        fListHistMultistrangeQA->Add(fHistEventSel);
  }
  if(! fHistMassXiMinus) {
     fHistMassXiMinus = new TH1F("fHistMassXiMinus", "#Xi^{-} candidates;M(#Lambda,#pi^{-}) (GeV/c^{2}); Counts", 150, 1.25, 1.40);
     fListHistMultistrangeQA->Add(fHistMassXiMinus);
  } 
  if(! fHistMassXiPlus) {
     fHistMassXiPlus = new TH1F("fHistMassXiPlus", "#Xi^{+} candidates; M(#bar{#Lambda}^{0},#pi^{+}) (GeV/c^{2}); Counts", 150, 1.25, 1.40);
     fListHistMultistrangeQA->Add(fHistMassXiPlus);
  }
  if(! fHistMassOmegaMinus) {
     fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus", "#Omega^{-} candidates; M(#Lambda,K^{-}) (GeV/c^{2}); Counts", 120, 1.62, 1.74);
     fListHistMultistrangeQA->Add(fHistMassOmegaMinus);
  }
  if(! fHistMassOmegaPlus) {
     fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus", "#Omega^{+} candidates;M(#bar{#Lambda}^{0},K^{+}) (GeV/c^{2}); Counts", 120, 1.62, 1.74); 
     fListHistMultistrangeQA->Add(fHistMassOmegaPlus);                                                                                                    
  }
  if(! fHistCascadeMultiplicityXiMinus) {
     fHistCascadeMultiplicityXiMinus = new TH1F("fHistCascadeMultiplicityXiMinus","Xi Minus per event;Nbr of Xi Minus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityXiMinus);
  } 
  if(! fHistCascadeMultiplicityXiPlus) {
     fHistCascadeMultiplicityXiPlus = new TH1F("fHistCascadeMultiplicityXiPlus","Xi Plus per event;Nbr of Xi Plus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityXiPlus);
  }
  if(! fHistCascadeMultiplicityOmegaMinus) {
     fHistCascadeMultiplicityOmegaMinus = new TH1F("fHistCascadeMultiplicityOmegaMinus","Omega Minus per event;Nbr of Omega Minus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityOmegaMinus);
  }
  if(! fHistCascadeMultiplicityOmegaPlus) {
     fHistCascadeMultiplicityOmegaPlus = new TH1F("fHistCascadeMultiplicityOmegaPlus","Omega Plus per event;Nbr of Omega Plus/Evt;Events", 50, 0, 50);
     fListHistMultistrangeQA->Add(fHistCascadeMultiplicityOmegaPlus);
  }


  //___________________________________________________
  // Define the container for the topological variables
  if(! fCFContCascadeCuts) {
      // NB: overflow/underflow of variables on which we want to cut later should be 0!!! 
      const Int_t  lNbSteps      =  4 ;
      const Int_t  lNbVariables  =  20; 
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
      lNbBinsPerVar[13] = 150;    //XiTransvMom                  :  [0.0,15.0]
      lNbBinsPerVar[14] = 110;    //Y(Xi)                        :   0.02 in rapidity units
      lNbBinsPerVar[15] = 110;    //Y(Omega)                     :   0.02 in rapidity units
      lNbBinsPerVar[16] = 112;    //Proper lenght of cascade       
      lNbBinsPerVar[17] = 112;    //Proper lenght of V0
      lNbBinsPerVar[18] = 201;    //V0CosineOfPointingAngleToXiV
      lNbBinsPerVar[19] = 11;     //Centrality
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
      fCFContCascadeCuts->SetBinLimits(13, 0.0, 15.0);
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
       //19 - Centrality
      Double_t *lBinLim19  = new Double_t[ lNbBinsPerVar[19]+1 ];
         for(Int_t i=3; i< lNbBinsPerVar[19]+1;i++)   lBinLim19[i]  = (Double_t)(i-1)*10.;
         lBinLim19[0] = 0.0; 
         lBinLim19[1] = 5.0; 
         lBinLim19[2] = 10.0;
      fCFContCascadeCuts->SetBinLimits(19,  lBinLim19 );     
      delete [] lBinLim19;
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
      fCFContCascadeCuts->SetVarTitle(18, "cos(V0 PA) in cascade to XiV");
      fCFContCascadeCuts->SetVarTitle(19, "Centrality");
  }


  //_______________________________________________
  // Define the Container for the MC generated info
  if(! fCFContCascadeMCgen) {
      // NB: overflow/underflow of variables on which we want to cut later should be 0!!! 
      const Int_t  lNbStepsMC      =  4 ; 
      const Int_t  lNbVariablesMC  =  7; 
      //Array for the number of bins in each dimension :
      Int_t lNbBinsPerVarMC[lNbVariablesMC] = {0};
      lNbBinsPerVarMC[0] = 200;    //Total momentum        : [0.0,20.0]
      lNbBinsPerVarMC[1] = 200;    //Transverse momentum   : [0.0,20.0]
      lNbBinsPerVarMC[2] = 110;    //Y                     : [-1.1,1.1]  
      lNbBinsPerVarMC[3] = 200;    //eta                   : [-10, 10]
      lNbBinsPerVarMC[4] = 200;    //theta                 : [-10, 190] 
      lNbBinsPerVarMC[5] = 360;    //Phi                   : [0., 360.]
      lNbBinsPerVarMC[6] = 11;     //Centrality
      //define the container
      fCFContCascadeMCgen = new AliCFContainer("fCFContCascadeMCgen","Container for MC gen cascade ", lNbStepsMC, lNbVariablesMC, lNbBinsPerVarMC );
      //Setting the bin limits 
       //0 - Total Momentum
      fCFContCascadeMCgen->SetBinLimits(0, 0.0, 20.0);
       //1 - Transverse Momentum 
      fCFContCascadeMCgen->SetBinLimits(1, 0.0, 20.0);
       //2 - Y
      fCFContCascadeMCgen->SetBinLimits(2, -1.1, 1.1);
       //3 - Eta
      fCFContCascadeMCgen->SetBinLimits(3, -10, 10);
       //4 - Theta
      fCFContCascadeMCgen->SetBinLimits(4, -10, 190);
       //5 - Phi
      fCFContCascadeMCgen->SetBinLimits(5, 0.0, 360.0);
       //6 - Centrality
      Double_t *lBinLim6MC  = new Double_t[ lNbBinsPerVarMC[6]+1 ];
         for(Int_t i=3; i< lNbBinsPerVarMC[6]+1;i++)   lBinLim6MC[i]  = (Double_t)(i-1)*10.;
         lBinLim6MC[0] = 0.0;
         lBinLim6MC[1] = 5.0;
         lBinLim6MC[2] = 10.0;
      fCFContCascadeMCgen->SetBinLimits(6,  lBinLim6MC);
      delete [] lBinLim6MC;
      // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
      fCFContCascadeMCgen->SetStepTitle(0, "#Xi^{-} candidates");
      fCFContCascadeMCgen->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
      fCFContCascadeMCgen->SetStepTitle(2, "#Omega^{-} candidates");
      fCFContCascadeMCgen->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
      // Setting the variable title, per axis
      fCFContCascadeMCgen->SetVarTitle(0,  "MC gen p_tot (GeV/c)");
      fCFContCascadeMCgen->SetVarTitle(1,  "MC gen p_T (GeV/c)");
      fCFContCascadeMCgen->SetVarTitle(2,  "MC gen Rapidity");
      fCFContCascadeMCgen->SetVarTitle(3,  "MC gen Pseudo-rapidity");
      fCFContCascadeMCgen->SetVarTitle(4,  "MC gen Theta");
      fCFContCascadeMCgen->SetVarTitle(5,  "MC gen Phi");
      fCFContCascadeMCgen->SetVarTitle(6,  "MC gen Centrality");
  }
 
  PostData(1, fListHistMultistrangeQA);
  PostData(2, fCFContCascadeCuts);
  PostData(3, fCFContCascadeMCgen);


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
   Float_t lPercentile = 500;
   Float_t lPercentileEmbeddedSelection = 500;
   Int_t lEvSelCode = 100;
   AliMultSelection *MultSelection = 0x0;
   if      (fAnalysisType == "ESD")  MultSelection = (AliMultSelection*) lESDevent->FindListObject("MultSelection");
   else if (fAnalysisType == "AOD")  MultSelection = (AliMultSelection*) lAODevent->FindListObject("MultSelection");
   if (!MultSelection) {
          AliWarning("AliMultSelection object not found!");  //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
          PostData(1, fListHistMultistrangeQA);
          PostData(2, fCFContCascadeCuts);
          PostData(3, fCFContCascadeMCgen);
          return;
   } else {
          AliWarning("AliMultSelection object found!");
          lPercentile = MultSelection->GetMultiplicityPercentile("V0M"); //V0M Multiplicity Percentile
          lPercentileEmbeddedSelection = MultSelection->GetMultiplicityPercentile("V0M", kTRUE );
          lEvSelCode = MultSelection->GetEvSelCode();  //Event Selection Code
   }
   // - Remove events
   if (lEvSelCode != 0) {
          AliWarning(Form("lEvSelCode value = %i. Run Not good! REMOVE",lEvSelCode));
          PostData(1, fListHistMultistrangeQA);
          PostData(2, fCFContCascadeCuts);
          PostData(3, fCFContCascadeMCgen);
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
             PostData(2, fCFContCascadeCuts);
             PostData(3, fCFContCascadeMCgen);
             return;
       }
       lPrimaryBestESDVtx->GetXYZ(lBestPrimaryVtxPos);
   } else if (fAnalysisType == "AOD") {
       const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();
       if (!lPrimaryBestAODVtx) {
             AliWarning("No prim. vertex in AOD... return!");
             PostData(1, fListHistMultistrangeQA);
             PostData(2, fCFContCascadeCuts);
             PostData(3, fCFContCascadeMCgen);
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

           Double_t partP      = 0.;
           Double_t partPt     = 0.;
           Double_t partEta    = 0.;
           Double_t partTheta  = 0.;
           Double_t partPhi    = 0.;
           Double_t partRap    = 0.;
           Double_t partEnergy = 0.; //for Rapidity
           Double_t partPz     = 0.; //for Rapidity
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

           Double_t lContainerCutVarsMC[7] = {0.0};
           lContainerCutVarsMC[0]  = partP;
           lContainerCutVarsMC[1]  = partPt;
           lContainerCutVarsMC[2]  = partRap;
           lContainerCutVarsMC[3]  = partEta;
           lContainerCutVarsMC[4]  = partTheta;
           lContainerCutVarsMC[5]  = partPhi;
           lContainerCutVarsMC[6]  = lPercentileEmbeddedSelection;
           if (PDGcode == 3312)  {fCFContCascadeMCgen->Fill(lContainerCutVarsMC,0); ngenximinus++;} // for Xi-
           if (PDGcode == -3312) {fCFContCascadeMCgen->Fill(lContainerCutVarsMC,1); ngenxiplus++;} // for Xi+
           if (PDGcode == 3334)  {fCFContCascadeMCgen->Fill(lContainerCutVarsMC,2); ngenomegaminus++;} // for Omega-
           if (PDGcode == -3334) {fCFContCascadeMCgen->Fill(lContainerCutVarsMC,3); ngenomegaplus++;} // for Omega+   
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
    if (lBachTransvMom<fMinPtCutOnDaughterTracks)   { AliWarning("ERROR: bachelor pT < lowlimit");          continue; }
    if (lpTrackTransvMom<fMinPtCutOnDaughterTracks) { AliWarning("ERROR: positive daughter pT < lowlimit"); continue; }
    if (lnTrackTransvMom<fMinPtCutOnDaughterTracks) { AliWarning("ERROR: negative daughter pT < lowlimit"); continue; }


    //---------------------------------------------------
    // Cut on pseudorapidity of the three daughter tracks
    if (TMath::Abs(etaBach) > 0.8) { AliWarning("ERROR: bachelor eta > maxlimit");          continue; }
    if (TMath::Abs(etaPos)  > 0.8) { AliWarning("ERROR: positive daughter eta > maxlimit"); continue; }
    if (TMath::Abs(etaNeg)  > 0.8) { AliWarning("ERROR: negative daughter eta > maxlimit"); continue; }


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
    // Fill the TH1F without PID info
    if        ( lChargeXi < 0 ) {
         fHistMassXiMinus->Fill( lInvMassXiMinus );
         fHistMassOmegaMinus->Fill( lInvMassOmegaMinus );
    } else if ( lChargeXi > 0 ) {
      fHistMassXiPlus->Fill( lInvMassXiPlus );
      fHistMassOmegaPlus->Fill( lInvMassOmegaPlus );
    }
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
    if (fisMC) lContainerCutVars[19] = lPercentileEmbeddedSelection;
    else       lContainerCutVars[19] = lPercentile;
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
  PostData(1, fListHistMultistrangeQA);
  PostData(2, fCFContCascadeCuts); 
  PostData(3, fCFContCascadeMCgen);

}// End UserExec


//-------------------------------------------------------
void AliAnalysisTaskQAMultistrange::Terminate(Option_t *) 

{

}