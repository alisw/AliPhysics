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
//            AliAnalysisTaskCheckCascadePbPb class
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
//              Adapted to PbPb analysis: M. Nicassio Feb2011, maria.nicassio@ba.infn.it
//              - Physics selection moved to the run.C macro
//              - Centrality selection added (+ setters)
//              - flag and setters added (CF container usage, vertex range)
//              - histo added and histo/container binning changed 
//              - protection in the destructor for CAF usage          
//              - AliWarning disabled
//              - number of tracklets from AOD also           
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
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
//      #include "AliV0vertexer.h"
//      #include "AliCascadeVertexer.h"
#include "AliESDpid.h"

#include "AliESDVZERO.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"

#include "AliAnalysisTaskCheckCascadePbPb.h"

ClassImp(AliAnalysisTaskCheckCascadePbPb)



//________________________________________________________________________
AliAnalysisTaskCheckCascadePbPb::AliAnalysisTaskCheckCascadePbPb() 
  : AliAnalysisTaskSE(), fAnalysisType("ESD"), fCollidingSystems(0), fESDpid(0), /*fPaveTextBookKeeping(0),*/
    fkRerunV0CascVertexers      (0),
    fkQualityCutZprimVtxPos     (kTRUE),
    fkQualityCutNoTPConlyPrimVtx(kTRUE),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCut80TPCcls        (kTRUE),
    fkExtraSelections           (0),
    fCentrLowLim                (0),
    fCentrUpLim                 (0),
    fCentrEstimator             (0),
    fVtxRange                   (0),
    fUseCFContCascadeCuts       (0),


    	// - Cascade part initialisation
    fListHistCascade(0),
    fHistCascadeMultiplicityBeforeEvSel(0),
    fHistCascadeMultiplicityForCentrEvt(0), fHistTrackMultiplicityForCentrEvt(0), fHistTPCrefitTrackMultiplicityForCentrEvt(0),
    fHistCascadeMultiplicityForTrigEvtAndZprimVtx(0),
    fHistCascadeMultiplicityForSelEvt(0),
    fHistPosBestPrimaryVtxXForSelEvt(0), fHistPosBestPrimaryVtxYForSelEvt(0), fHistPosBestPrimaryVtxZForSelEvt(0),
    fHistTPCrefitTrackMultiplicityForCascadeEvt(0),
    fHistPosV0TPCClusters(0), fHistNegV0TPCClusters(0), fHistBachTPCClusters(0),
    fHistVtxStatus(0),

    fHistPosTrkgPrimaryVtxXForCascadeEvt(0), fHistPosTrkgPrimaryVtxYForCascadeEvt(0), fHistPosTrkgPrimaryVtxZForCascadeEvt(0), fHistTrkgPrimaryVtxRadius(0),
    fHistPosBestPrimaryVtxXForCascadeEvt(0), fHistPosBestPrimaryVtxYForCascadeEvt(0), fHistPosBestPrimaryVtxZForCascadeEvt(0), fHistBestPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsBestPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),

    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2Xi(0), 
    fHistDcaV0DaughtersXi(0),
    fHistDcaV0ToPrimVertexXi(0), 
    fHistV0CosineOfPointingAngleXi(0),
    fHistV0RadiusXi(0),
    fHistDcaPosToPrimVertexXi(0), fHistDcaNegToPrimVertexXi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMomXi(0),   fHistBachTotMomXi(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEtaXi(0),
    fHistThetaXi(0), fHistPhiXi(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0),
    f2dHistXiRadiusVsEffMassOmegaMinus(0), f2dHistXiRadiusVsEffMassOmegaPlus(0),
    
    f2dHistTPCdEdxOfCascDghters(0),
    
//    f3dHistXiPtVsEffMassVsYXiMinus(0), f3dHistXiPtVsEffMassVsYXiPlus(0),
//    f3dHistXiPtVsEffMassVsYOmegaMinus(0), f3dHistXiPtVsEffMassVsYOmegaPlus(0),
    
    fCFContCascadePIDXiMinus(0),
    fCFContCascadePIDXiPlus(0),
    fCFContCascadePIDOmegaMinus(0),
    fCFContCascadePIDOmegaPlus(0),
    fCFContCascadeCuts(0),
    
//    fHnSpAngularCorrXiMinus(0), fHnSpAngularCorrXiPlus(0), 
//    fHnSpAngularCorrOmegaMinus(0), fHnSpAngularCorrOmegaPlus(0),
    fV0Ampl(0),

    fHistDcaXiDaughtersvsInvMass(0), fHistDcaBachToPrimVertexvsInvMass(0), fHistXiCosineOfPointingAnglevsInvMass(0),
    fHistMassLambdaAsCascDghtervsInvMass(0),fHistDcaV0DaughtersXivsInvMass(0),fHistDcaV0ToPrimVertexXivsInvMass(0)


{
  // Dummy Constructor
        for(Int_t iAlephIdx   = 0; iAlephIdx   < 5; iAlephIdx++   ) { fAlephParameters [iAlephIdx]    = -1.; }
        for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
        for(Int_t iCascSelIdx = 0; iCascSelIdx < 8; iCascSelIdx++ ) { fCascSels        [iCascSelIdx ] = -1.; }
}


//________________________________________________________________________
AliAnalysisTaskCheckCascadePbPb::AliAnalysisTaskCheckCascadePbPb(const char *name) 
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fCollidingSystems(0), fESDpid(0), /*fPaveTextBookKeeping(0),*/
    fkRerunV0CascVertexers      (0),
    fkQualityCutZprimVtxPos     (kTRUE),
    fkQualityCutNoTPConlyPrimVtx(kTRUE),
    fkQualityCutTPCrefit        (kTRUE),
    fkQualityCut80TPCcls        (kTRUE),
    fkExtraSelections           (0),
    fCentrLowLim                (0),
    fCentrUpLim                 (0),
    fCentrEstimator             (0),
    fVtxRange                   (0),
    fUseCFContCascadeCuts       (0),
     
    	// - Cascade part initialisation
    fListHistCascade(0),
    fHistCascadeMultiplicityBeforeEvSel(0),
    fHistCascadeMultiplicityForCentrEvt(0), fHistTrackMultiplicityForCentrEvt(0), fHistTPCrefitTrackMultiplicityForCentrEvt(0),
    fHistCascadeMultiplicityForTrigEvtAndZprimVtx(0),
    fHistCascadeMultiplicityForSelEvt(0),
    fHistPosBestPrimaryVtxXForSelEvt(0), fHistPosBestPrimaryVtxYForSelEvt(0), fHistPosBestPrimaryVtxZForSelEvt(0),
    fHistTPCrefitTrackMultiplicityForCascadeEvt(0),
    fHistPosV0TPCClusters(0), fHistNegV0TPCClusters(0), fHistBachTPCClusters(0),
    fHistVtxStatus(0),

    fHistPosTrkgPrimaryVtxXForCascadeEvt(0), fHistPosTrkgPrimaryVtxYForCascadeEvt(0), fHistPosTrkgPrimaryVtxZForCascadeEvt(0), fHistTrkgPrimaryVtxRadius(0),
    fHistPosBestPrimaryVtxXForCascadeEvt(0), fHistPosBestPrimaryVtxYForCascadeEvt(0), fHistPosBestPrimaryVtxZForCascadeEvt(0), fHistBestPrimaryVtxRadius(0),
    f2dHistTrkgPrimVtxVsBestPrimVtx(0),

    fHistEffMassXi(0),  fHistChi2Xi(0),  
    fHistDcaXiDaughters(0), fHistDcaBachToPrimVertex(0), fHistXiCosineOfPointingAngle(0), fHistXiRadius(0),

    fHistMassLambdaAsCascDghter(0),
    fHistV0Chi2Xi(0),
    fHistDcaV0DaughtersXi(0),
    fHistDcaV0ToPrimVertexXi(0), 
    fHistV0CosineOfPointingAngleXi(0),
    fHistV0RadiusXi(0),
    fHistDcaPosToPrimVertexXi(0), fHistDcaNegToPrimVertexXi(0), 

    fHistMassXiMinus(0), fHistMassXiPlus(0),
    fHistMassOmegaMinus(0), fHistMassOmegaPlus(0),
    fHistMassWithCombPIDXiMinus(0), fHistMassWithCombPIDXiPlus(0),
    fHistMassWithCombPIDOmegaMinus(0), fHistMassWithCombPIDOmegaPlus(0),

    fHistXiTransvMom(0),    fHistXiTotMom(0),
    fHistBachTransvMomXi(0),   fHistBachTotMomXi(0),

    fHistChargeXi(0),
    fHistV0toXiCosineOfPointingAngle(0),

    fHistRapXi(0), fHistRapOmega(0), fHistEtaXi(0),
    fHistThetaXi(0), fHistPhiXi(0),

    f2dHistArmenteros(0),			
    f2dHistEffMassLambdaVsEffMassXiMinus(0), f2dHistEffMassXiVsEffMassOmegaMinus(0),
    f2dHistEffMassLambdaVsEffMassXiPlus(0), f2dHistEffMassXiVsEffMassOmegaPlus(0),
    f2dHistXiRadiusVsEffMassXiMinus(0), f2dHistXiRadiusVsEffMassXiPlus(0),
    f2dHistXiRadiusVsEffMassOmegaMinus(0), f2dHistXiRadiusVsEffMassOmegaPlus(0),
    
    f2dHistTPCdEdxOfCascDghters(0),
    
//    f3dHistXiPtVsEffMassVsYXiMinus(0), f3dHistXiPtVsEffMassVsYXiPlus(0),
//    f3dHistXiPtVsEffMassVsYOmegaMinus(0), f3dHistXiPtVsEffMassVsYOmegaPlus(0),
    
    fCFContCascadePIDXiMinus(0),
    fCFContCascadePIDXiPlus(0),
    fCFContCascadePIDOmegaMinus(0),
    fCFContCascadePIDOmegaPlus(0),
    fCFContCascadeCuts(0),
    
//    fHnSpAngularCorrXiMinus(0), fHnSpAngularCorrXiPlus(0), 
//    fHnSpAngularCorrOmegaMinus(0), fHnSpAngularCorrOmegaPlus(0),
    fV0Ampl(0),

    fHistDcaXiDaughtersvsInvMass(0), fHistDcaBachToPrimVertexvsInvMass(0), fHistXiCosineOfPointingAnglevsInvMass(0),
    fHistMassLambdaAsCascDghtervsInvMass(0),fHistDcaV0DaughtersXivsInvMass(0),fHistDcaV0ToPrimVertexXivsInvMass(0)


{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TList container (cascade)
        
        for(Int_t iAlephIdx   = 0; iAlephIdx   < 5; iAlephIdx++   ) { fAlephParameters [iAlephIdx]    = -1.; }
        
        fV0Sels[0] =  33.  ;  // max allowed chi2
        fV0Sels[1] =   0.001; // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        fV0Sels[2] =   0.001; // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        fV0Sels[3] =   0.2 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
        fV0Sels[4] =   0.0 ;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
        fV0Sels[5] =   0.1 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
        fV0Sels[6] = 100.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)

        fCascSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
        fCascSels[1] =   0.001;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
        fCascSels[2] =   0.008;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
        fCascSels[3] =   0.03;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
        fCascSels[4] =   0.3  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
        fCascSels[5] =   0.9998 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
        fCascSels[6] =   0.1  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
        fCascSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )

  // Output slot #0 writes into a TList container (Cascade)
  DefineOutput(1, TList::Class());
  /*DefineOutput(2, TPaveText::Class());*/
  AliLog::SetClassDebugLevel("AliAnalysisTaskCheckCascadePbPb",1); // MN this should (?) enable only AliFatal
}


AliAnalysisTaskCheckCascadePbPb::~AliAnalysisTaskCheckCascadePbPb()
{
  //
  // Destructor
  //

  // For all TH1, 2, 3 HnSparse and CFContainer are in the fListCascade TList.
  // They will be deleted when fListCascade is deleted by the TSelector dtor
  // Because of TList::SetOwner() ...
        
  if (fListHistCascade && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())         { delete fListHistCascade;     fListHistCascade = 0x0;    }
  if (fESDpid)                  { delete fESDpid;              fESDpid = 0x0;} // fESDpid is not stored into the TList
  //if (fPaveTextBookKeeping)     { delete fPaveTextBookKeeping; fPaveTextBookKeeping = 0x0;} // fPaveTextBookKeeping is not strored in the TList
}



//________________________________________________________________________
void AliAnalysisTaskCheckCascadePbPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once



 fListHistCascade = new TList();
 fListHistCascade->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
 
 
if(! fESDpid){

  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>( (AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler() );
  
  if( !lmcEvtHandler ) { // !0x0 = real data or !1 = there is an MC handler available (useMC = kTRUE in AnalysisTrainNew), so = data from MC
        
        // Reasonable parameters extracted for real p-p event (Dec 2009 - GSI Pass5) - A.Kalweit
        fAlephParameters[0] = 0.0283086;        // No extra-division to apply in SetBlochParam
        fAlephParameters[1] = 2.63394e+01;
        fAlephParameters[2] = 5.04114e-11;
        fAlephParameters[3] = 2.12543e+00;
        fAlephParameters[4] = 4.88663e+00; 
        Printf("CheckCascade - Check Aleph Param in case of REAL Data (fAlephParameters[1] = %f)\n",  fAlephParameters[1]);
  } else {
        // Reasonable parameters extracted for p-p simulation (LHC09a4) - A.Kalweit
         // fAlephParameters[0] = 4.23232575531564326e+00;//50*0.76176e-1; // do not forget to divide this value by 50 in SetBlochParam !
         // fAlephParameters[1] = 8.68482806165147636e+00;//10.632; 
         // fAlephParameters[2] = 1.34000000000000005e-05;//0.13279e-4;
         // fAlephParameters[3] = 2.30445734159456084e+00;//1.8631;
         // fAlephParameters[4] = 2.25624744086878559e+00;//1.9479;  
          
        // Reasonable parameters extracted for MC LHC09d10 event (Jan 2010) - A.Kalweit
        fAlephParameters[0] = 2.15898e+00/50.;
        fAlephParameters[1] = 1.75295e+01;
        fAlephParameters[2] = 3.40030e-09;
        fAlephParameters[3] = 1.96178e+00;
        fAlephParameters[4] = 3.91720e+00;
        Printf("CheckCascade - Check Aleph Param win case MC data (fAlephParameters[1] = %f)\n",  fAlephParameters[1]);
  }

  fESDpid = new AliESDpid();
  fESDpid->GetTPCResponse().SetBetheBlochParameters( fAlephParameters[0],
                                                     fAlephParameters[1],
                                                     fAlephParameters[2],
                                                     fAlephParameters[3],
                                                     fAlephParameters[4] );
}
 
 
/* 
if( !fPaveTextBookKeeping){
        // FIXME : prepare a field with the AliRoot+root distrib ...

        fPaveTextBookKeeping = new TPaveText(0.1, 0.1, 0.9, 0.9,"NDC");
        fPaveTextBookKeeping->SetName("fPaveTextBookKeeping");
        fPaveTextBookKeeping->SetBorderSize(0);
        fPaveTextBookKeeping->SetTextAlign(12);
        fPaveTextBookKeeping->SetFillColor(kWhite);
        fPaveTextBookKeeping->SetTextFont(42);        // regular Arial or Helvetica,
        fPaveTextBookKeeping->SetTextColor(kBlue+4);
        
        fPaveTextBookKeeping->AddText( "Task CHECK CASCADE analysis" );
        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
        fPaveTextBookKeeping->AddText( Form("AnalysisType : %s ", fAnalysisType.Data() ));
        if(!fCollidingSystems)  fPaveTextBookKeeping->AddText("Colliding system : p-p collisions ");
        else                    fPaveTextBookKeeping->AddText("Colliding system : A-A collisions ");

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
    
        if(fkRerunV0CascVertexers){
                fPaveTextBookKeeping->AddText("A.1. With V0 vertexer : ");
                fPaveTextBookKeeping->AddText( Form("  - V0 #chi^{2} _________________ <  %.3f ",               fV0Sels[0] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ 1^{st} daughter) ___ >  %.3f     cm ",  fV0Sels[1] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ 2^{nd} daughter) __  >  %.3f     cm",   fV0Sels[2] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA between V0 daughters ___ <  %.3f      cm",         fV0Sels[3] ));
                fPaveTextBookKeeping->AddText( Form("  - cos(V0 pointing angle) _______ >  %.3f ",              fV0Sels[4] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(V0 decay) ________ >  %.3f             cm", fV0Sels[5] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(V0 decay) ________ <  %.3f         cm",     fV0Sels[6] ));
                
                fPaveTextBookKeeping->AddText(" "); 
                
                fPaveTextBookKeeping->AddText("A.2. With Casc. vertexer : ");
                fPaveTextBookKeeping->AddText( Form("  - Casc. #chi^{2} ______________  <  %.3f ",                               fCascSels[0] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ V0) _________ >  %.3f    cm",                            fCascSels[1] ));
                fPaveTextBookKeeping->AddText( Form("  - | M_{#Lambda}(reco) - M_{#Lambda}(pdg) | _______ <  %.3f    GeV/c^{2}", fCascSels[2] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA(prim. Vtx/ Bach) _______ >  %.3f    cm",                            fCascSels[3] ));
                fPaveTextBookKeeping->AddText( Form("  - DCA between Bach/ #Lambda ______ <  %.3f    cm",                        fCascSels[4] ));
                fPaveTextBookKeeping->AddText( Form("  - cos(Casc. pointing angle) ____ >  %.3f ",                               fCascSels[5] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(Casc. decay) ______ >  %.3f         cm",                     fCascSels[6] ));
                fPaveTextBookKeeping->AddText( Form("  - R_{transv}(Casc. decay) ______ <  %.3f     cm",                         fCascSels[7] ));
        }
        else{   fPaveTextBookKeeping->AddText("A. No rerunning of the V0/Casc. vertexers ... See std cuts in (AliRoot+Rec.C) used for this prod. cycle");}

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");
        
        if(fkQualityCutZprimVtxPos)      fPaveTextBookKeeping->AddText("B. Quality Cut(prim. Vtx z-Pos)    = ON  ");
        else                             fPaveTextBookKeeping->AddText("B. Quality Cut(prim. Vtx z-Pos)    = Off ");
        if(fkQualityCutNoTPConlyPrimVtx) fPaveTextBookKeeping->AddText("C. Quality Cut(No TPC prim. vtx) = ON  ");
        else                             fPaveTextBookKeeping->AddText("C. Quality Cut(No TPC prim. vtx) = Off ");
        if(fkQualityCutTPCrefit)         fPaveTextBookKeeping->AddText("D. Quality Cut(TPCrefit)               = ON  ");
        else                             fPaveTextBookKeeping->AddText("D. Quality Cut(TPCrefit)               = Off ");
        if(fkQualityCut80TPCcls)         fPaveTextBookKeeping->AddText("E. Quality Cut(80 TPC clusters)   = ON  ");
        else                             fPaveTextBookKeeping->AddText("E. Quality Cut(80 TPC clusters)   = Off ");
        if(fkExtraSelections)            fPaveTextBookKeeping->AddText("F. Extra Analysis Selections         = ON  ");
        else                             fPaveTextBookKeeping->AddText("F. Extra Analysis Selections         = Off ");

        fPaveTextBookKeeping->AddText("- - - - - - - - - - - ");

        fPaveTextBookKeeping->AddText("G. TPC Aleph Param : ");
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [0] =  %.5g", fAlephParameters[0] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [1] =  %.5g", fAlephParameters[1] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [2] =  %.5g", fAlephParameters[2] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [3] =  %.5g", fAlephParameters[3] ));
        fPaveTextBookKeeping->AddText( Form("   - fAlephParam [4] =  %.5g", fAlephParameters[4] ));
}       
*/ 
 
	// - General histos
	//--------------
 
if(! fHistCascadeMultiplicityBeforeEvSel) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicityBeforeEvSel = new TH1F("fHistCascadeMultiplicityBeforeEvSel", 
			"Cascades per event (before vertex and centrality selections);Nbr of Cascades/Evt;Events", 
			100, 0, 200000); 		
	else // pp collisions
	  	fHistCascadeMultiplicityBeforeEvSel = new TH1F("fHistCascadeMultiplicityBeforeEvSel", 
			"Cascades per event (before vertex and centrality selections);Nbr of Cascades/Evt;Events", 
			25, 0, 25);
	fListHistCascade->Add(fHistCascadeMultiplicityBeforeEvSel);
}

        // - Histos for events passing the trigger selection
        //--------------
        
if(! fHistCascadeMultiplicityForCentrEvt) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicityForCentrEvt = new TH1F("fHistCascadeMultiplicityForCentrEvt", 
			"Cascades per event (for triggered evt);Nbr of Cascades/Evt;Events", 
			100, 0, 200000); 		
	else // pp collisions
	  	fHistCascadeMultiplicityForCentrEvt = new TH1F("fHistCascadeMultiplicityForCentrEvt", 
			"Cascades per event (for triggered evt);Nbr of Cascades/Evt;Events", 
			25, 0, 25);
	fListHistCascade->Add(fHistCascadeMultiplicityForCentrEvt);
}

 
if(! fHistTrackMultiplicityForCentrEvt) {	
	if(fCollidingSystems)// AA collisions	
		fHistTrackMultiplicityForCentrEvt = new TH1F("fHistTrackMultiplicityForCentrEvt", 
			"Track Multiplicity (for triggered evt);Nbr of tracks/Evt;Events", 
			200, 0, 12000); 		
	else // pp collisions
		fHistTrackMultiplicityForCentrEvt = new TH1F("fHistTrackMultiplicityForCentrEvt", 
			"Track Multiplicity (for triggered evt);Nbr of tracks/Evt;Events", 
			300, 0, 300);
	fListHistCascade->Add(fHistTrackMultiplicityForCentrEvt);
}

if(! fHistTPCrefitTrackMultiplicityForCentrEvt) {	
	if(fCollidingSystems)// AA collisions	
		fHistTPCrefitTrackMultiplicityForCentrEvt = new TH1F("fHistTPCrefitTrackMultiplicityForCentrEvt", 
			"TPCrefit track Multiplicity (for triggered evt);Nbr of TPCrefit tracks/Evt;Events", 
			200, 0, 12000); 		
	else // pp collisions
		fHistTPCrefitTrackMultiplicityForCentrEvt = new TH1F("fHistTPCrefitTrackMultiplicityForCentrEvt", 
			"TPCrefit track Multiplicity (for triggered evt);Nbr of TPCrefit tracks/Evt;Events", 
			300, 0, 300);
	fListHistCascade->Add(fHistTPCrefitTrackMultiplicityForCentrEvt);
}




        // - Histos for events passing the trigger selection + |z(prim. vertex)| < XX cm
        //--------------
        
if(! fHistCascadeMultiplicityForTrigEvtAndZprimVtx) {
        if(fCollidingSystems)// AA collisions
                fHistCascadeMultiplicityForTrigEvtAndZprimVtx = new TH1F("fHistCascadeMultiplicityForTrigEvtAndZprimVtx", 
                        "Cascades per event;Nbr of Cascades/Evt;Events", 
                        100, 0, 200000);           
        else // pp collisions
                fHistCascadeMultiplicityForTrigEvtAndZprimVtx = new TH1F("fHistCascadeMultiplicityForTrigEvtAndZprimVtx", 
                        "Cascades per event;Nbr of Cascades/Evt;Events", 
                        25, 0, 25);
        fListHistCascade->Add(fHistCascadeMultiplicityForTrigEvtAndZprimVtx);
}


        // - Histos for events passing the event selection at the analysis level
        //--------------
        
if(! fHistCascadeMultiplicityForSelEvt) {
	if(fCollidingSystems)// AA collisions
		fHistCascadeMultiplicityForSelEvt = new TH1F("fHistCascadeMultiplicityForSelEvt", 
			"Cascades per event;Nbr of Cascades/Evt;Events", 
			100, 0, 200000); 		
	else // pp collisions
	  	fHistCascadeMultiplicityForSelEvt = new TH1F("fHistCascadeMultiplicityForSelEvt", 
			"Cascades per event;Nbr of Cascades/Evt;Events", 
			25, 0, 25);
	fListHistCascade->Add(fHistCascadeMultiplicityForSelEvt);
}

if(! fHistPosBestPrimaryVtxXForSelEvt ){
	fHistPosBestPrimaryVtxXForSelEvt   = new TH1F( "fHistPosBestPrimaryVtxXForSelEvt" , "Best Prim. Vertex Position in x; x (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxXForSelEvt);
}

if(! fHistPosBestPrimaryVtxYForSelEvt){
	fHistPosBestPrimaryVtxYForSelEvt   = new TH1F( "fHistPosBestPrimaryVtxYForSelEvt" , "Best Prim. Vertex Position in y; y (cm); Events" , 360, -0.9, 0.9 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxYForSelEvt);
}

if(! fHistPosBestPrimaryVtxZForSelEvt ){
	fHistPosBestPrimaryVtxZForSelEvt   = new TH1F( "fHistPosBestPrimaryVtxZForSelEvt" , "Best Prim. Vertex Position in z; z (cm); Events" , 300, -30.0, 30.0 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxZForSelEvt);
}




        // - Histos for events containing at least ONE CASCADE
        //--------------
        
if(! fHistTPCrefitTrackMultiplicityForCascadeEvt) {
	if(fCollidingSystems)// AA collisions	
		fHistTPCrefitTrackMultiplicityForCascadeEvt = new TH1F("fHistTPCrefitTrackMultiplicityForCascadeEvt", 
			"TPCrefit track Multiplicity (for evt with Casc.);Nbr of TPCrefit tracks/Evt with cascade(s);Events", 
			200, 0, 20000); 		
	else // pp collisions
		fHistTPCrefitTrackMultiplicityForCascadeEvt = new TH1F("fHistTPCrefitTrackMultiplicityForCascadeEvt", 
			"TPCrefit track Multiplicity (for evt with Casc.);Nbr of TPCrefit tracks/Evt with cascade(s);Events", 
			300, 0, 300);
	fListHistCascade->Add(fHistTPCrefitTrackMultiplicityForCascadeEvt);
}

if(! fHistPosV0TPCClusters ){
        fHistPosV0TPCClusters = new TH1F("fHistPosV0TPCClusters", "TPC clusters for Pos. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Pos.); Track counts", 165, 0.0 ,165.0);
        fListHistCascade->Add(fHistPosV0TPCClusters);
}

if(! fHistNegV0TPCClusters ){
        fHistNegV0TPCClusters = new TH1F("fHistNegV0TPCClusters", "TPC clusters for Neg. V0 daughter track, in Casc; Nbr of TPC clusters (V0 Neg.); Track counts", 165, 0.0 ,165.0);
        fListHistCascade->Add(fHistNegV0TPCClusters);
}

if(! fHistBachTPCClusters ){
        fHistBachTPCClusters = new TH1F("fHistBachTPCClusters", "TPC clusters for Bachelor track; Nbr of TPC clusters (Bach); Track counts", 165, 0.0 ,165.0);
        fListHistCascade->Add(fHistBachTPCClusters);
}



if(! fHistVtxStatus ){
	fHistVtxStatus   = new TH1F( "fHistVtxStatus" , "Does a Trckg Prim.vtx exist ?; true=1 or false=0; Nb of Events" , 4, -1.0, 3.0 );
	fListHistCascade->Add(fHistVtxStatus);
}


	// - Vertex Positions
  
if(! fHistPosTrkgPrimaryVtxXForCascadeEvt ){
	fHistPosTrkgPrimaryVtxXForCascadeEvt   = new TH1F( "fHistPosTrkgPrimaryVtxXForCascadeEvt" , "Trkg Prim. Vertex Position in x; x (cm); Events" , 120, -0.6, 0.6 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxXForCascadeEvt);
}


if(! fHistPosTrkgPrimaryVtxYForCascadeEvt){
	fHistPosTrkgPrimaryVtxYForCascadeEvt   = new TH1F( "fHistPosTrkgPrimaryVtxYForCascadeEvt" , "Trkg Prim. Vertex Position in y; y (cm); Events" , 120, -0.6, 0.6 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxYForCascadeEvt);
}

if(! fHistPosTrkgPrimaryVtxZForCascadeEvt ){
	fHistPosTrkgPrimaryVtxZForCascadeEvt   = new TH1F( "fHistPosTrkgPrimaryVtxZForCascadeEvt" , "Trkg Prim. Vertex Position in z; z (cm); Events" , 200, -20.0, 20.0 );
	fListHistCascade->Add(fHistPosTrkgPrimaryVtxZForCascadeEvt);
}

if(! fHistTrkgPrimaryVtxRadius ){
	fHistTrkgPrimaryVtxRadius  = new TH1F( "fHistTrkgPrimaryVtxRadius",  "Trkg Prim. Vertex radius; r (cm); Events" , 150, 0., 15.0 );
	fListHistCascade->Add(fHistTrkgPrimaryVtxRadius);
}




if(! fHistPosBestPrimaryVtxXForCascadeEvt ){
	fHistPosBestPrimaryVtxXForCascadeEvt   = new TH1F( "fHistPosBestPrimaryVtxXForCascadeEvt" , "Best Prim. Vertex Position in x; x (cm); Events" , 120, -0.6, 0.6 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxXForCascadeEvt);
}

if(! fHistPosBestPrimaryVtxYForCascadeEvt){
	fHistPosBestPrimaryVtxYForCascadeEvt   = new TH1F( "fHistPosBestPrimaryVtxYForCascadeEvt" , "Best Prim. Vertex Position in y; y (cm); Events" , 120, -0.6, 0.6 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxYForCascadeEvt);
}

if(! fHistPosBestPrimaryVtxZForCascadeEvt ){
	fHistPosBestPrimaryVtxZForCascadeEvt   = new TH1F( "fHistPosBestPrimaryVtxZForCascadeEvt" , "Best Prim. Vertex Position in z; z (cm); Events" , 200, -20.0, 20.0 );
	fListHistCascade->Add(fHistPosBestPrimaryVtxZForCascadeEvt);
}

if(! fHistBestPrimaryVtxRadius ){
	fHistBestPrimaryVtxRadius  = new TH1F( "fHistBestPrimaryVtxRadius",  "Best Prim.  vertex radius; r (cm); Events" , 150, 0., 15.0 );
	fListHistCascade->Add(fHistBestPrimaryVtxRadius);
}

if(! f2dHistTrkgPrimVtxVsBestPrimVtx) {
	f2dHistTrkgPrimVtxVsBestPrimVtx = new TH2F( "f2dHistTrkgPrimVtxVsBestPrimVtx", "r_{Trck Prim. Vtx} Vs r_{Best Prim. Vtx}; r_{Track Vtx} (cm); r_{Best Vtx} (cm)", 150, 0., 15.0, 150, 0., 15.);
	fListHistCascade->Add(f2dHistTrkgPrimVtxVsBestPrimVtx);
}




// - Typical histos for cascades


if(! fHistEffMassXi) {
     fHistEffMassXi = new TH1F("fHistEffMassXi", "Cascade candidates ; Invariant Mass (GeV/c^{2}) ; Counts", 400, 1.2, 2.0);
     fListHistCascade->Add(fHistEffMassXi);
}
   
if(! fHistChi2Xi ){
	fHistChi2Xi = new TH1F("fHistChi2Xi", "Cascade #chi^{2}; #chi^{2}; Number of Cascades", 160, 0, 40);
	fListHistCascade->Add(fHistChi2Xi);
}
  
if(! fHistDcaXiDaughters ){
	fHistDcaXiDaughters = new TH1F( "fHistDcaXiDaughters",  "DCA between Xi Daughters; DCA (cm) ; Number of Cascades", 100, 0., 0.5);
	fListHistCascade->Add(fHistDcaXiDaughters);
}

if(! fHistDcaBachToPrimVertex) {
	fHistDcaBachToPrimVertex = new TH1F("fHistDcaBachToPrimVertex", "DCA of Bach. to Prim. Vertex;DCA (cm);Number of Cascades", 250, 0., 0.25);
	fListHistCascade->Add(fHistDcaBachToPrimVertex);
}

if(! fHistXiCosineOfPointingAngle) {
	fHistXiCosineOfPointingAngle = new TH1F("fHistXiCosineOfPointingAngle", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 200, 0.99, 1.0);
	fListHistCascade->Add(fHistXiCosineOfPointingAngle);
}

if(! fHistXiRadius ){
	fHistXiRadius  = new TH1F( "fHistXiRadius",  "Casc. decay transv. radius; r (cm); Counts" , 1050, 0., 105.0 );
	fListHistCascade->Add(fHistXiRadius);
}


// - Histos about ~ the "V0 part" of the cascade,  coming by inheritance from AliESDv0



if (! fHistMassLambdaAsCascDghter) {
     fHistMassLambdaAsCascDghter = new TH1F("fHistMassLambdaAsCascDghter","#Lambda associated to Casc. candidates;Eff. Mass (GeV/c^{2});Counts", 300,1.00,1.3);
    fListHistCascade->Add(fHistMassLambdaAsCascDghter);
}

if (! fHistV0Chi2Xi) {
	fHistV0Chi2Xi = new TH1F("fHistV0Chi2Xi", "V0 #chi^{2}, in cascade; #chi^{2};Counts", 160, 0, 40);
	fListHistCascade->Add(fHistV0Chi2Xi);
}

if (! fHistDcaV0DaughtersXi) {
	fHistDcaV0DaughtersXi = new TH1F("fHistDcaV0DaughtersXi", "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 0.6);
	fListHistCascade->Add(fHistDcaV0DaughtersXi);
}

if (! fHistDcaV0ToPrimVertexXi) {
	fHistDcaV0ToPrimVertexXi = new TH1F("fHistDcaV0ToPrimVertexXi", "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
	fListHistCascade->Add(fHistDcaV0ToPrimVertexXi);
}

if (! fHistV0CosineOfPointingAngleXi) {
	fHistV0CosineOfPointingAngleXi = new TH1F("fHistV0CosineOfPointingAngleXi", "Cosine of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts", 200, 0.98, 1.0);
	fListHistCascade->Add(fHistV0CosineOfPointingAngleXi);
}

if (! fHistV0RadiusXi) {
	fHistV0RadiusXi = new TH1F("fHistV0RadiusXi", "V0 decay radius, in cascade; radius (cm); Counts", 1050, 0., 105.0);
	fListHistCascade->Add(fHistV0RadiusXi);
}

if (! fHistDcaPosToPrimVertexXi) {
	fHistDcaPosToPrimVertexXi = new TH1F("fHistDcaPosToPrimVertexXi", "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHistCascade->Add(fHistDcaPosToPrimVertexXi);
}

if (! fHistDcaNegToPrimVertexXi) {
	fHistDcaNegToPrimVertexXi = new TH1F("fHistDcaNegToPrimVertexXi", "DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts", 300, 0, 3);
	fListHistCascade->Add(fHistDcaNegToPrimVertexXi);
}




	// - Effective mass histos for cascades.
// By cascade hyp  
if (! fHistMassXiMinus) {
    fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
    fListHistCascade->Add(fHistMassXiMinus);
}
  
if (! fHistMassXiPlus) {
    fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassXiPlus);
}

if (! fHistMassOmegaMinus) {
	fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassOmegaMinus);
}
 
if (! fHistMassOmegaPlus) {
	fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassOmegaPlus);
}

// By cascade hyp + bachelor PID
if (! fHistMassWithCombPIDXiMinus) {
    fHistMassWithCombPIDXiMinus = new TH1F("fHistMassWithCombPIDXiMinus","#Xi^{-} candidates, with Bach. comb. PID;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiMinus);
}
  
if (! fHistMassWithCombPIDXiPlus) {
    fHistMassWithCombPIDXiPlus = new TH1F("fHistMassWithCombPIDXiPlus","#Xi^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , #pi^{+} ) (GeV/c^{2});Counts",400,1.2,2.0);
    fListHistCascade->Add(fHistMassWithCombPIDXiPlus);
}

if (! fHistMassWithCombPIDOmegaMinus) {
	fHistMassWithCombPIDOmegaMinus = new TH1F("fHistMassWithCombPIDOmegaMinus","#Omega^{-} candidates, with Bach. comb. PID;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaMinus);
}
 
if (! fHistMassWithCombPIDOmegaPlus) {
	fHistMassWithCombPIDOmegaPlus = new TH1F("fHistMassWithCombPIDOmegaPlus","#Omega^{+} candidates, with Bach. comb. PID;M( #bar{#Lambda}^{0} , K^{+} ) (GeV/c^{2});Counts", 500,1.5,2.5);
    fListHistCascade->Add(fHistMassWithCombPIDOmegaPlus);
}



	// - Complements for QA

if(! fHistXiTransvMom ){
	fHistXiTransvMom  = new TH1F( "fHistXiTransvMom" , "#Xi transverse momentum (cand. around the mass peak) ; p_{t}(#Xi) (GeV/c); Counts", 100, 0.0, 10.0);
	fListHistCascade->Add(fHistXiTransvMom);
}

if(! fHistXiTotMom ){
	fHistXiTotMom  = new TH1F( "fHistXiTotMom" , "#Xi momentum norm (cand. around the mass peak); p_{tot}(#Xi) (GeV/c); Counts", 150, 0.0, 15.0);
	fListHistCascade->Add(fHistXiTotMom);
}


if(! fHistBachTransvMomXi ){
	fHistBachTransvMomXi  = new TH1F( "fHistBachTransvMomXi" , "#Xi Bach. transverse momentum (cand. around the mass peak) ; p_{t}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
	fListHistCascade->Add(fHistBachTransvMomXi);
}

if(! fHistBachTotMomXi ){
	fHistBachTotMomXi  = new TH1F( "fHistBachTotMomXi" , "#Xi Bach. momentum norm (cand. around the mass peak); p_{tot}(Bach.) (GeV/c); Counts", 100, 0.0, 5.0);
	fListHistCascade->Add(fHistBachTotMomXi);
}


if(! fHistChargeXi ){
	fHistChargeXi  = new TH1F( "fHistChargeXi" , "Charge of casc. candidates ; Sign ; Counts", 5, -2.0, 3.0);
	fListHistCascade->Add(fHistChargeXi);
}


if (! fHistV0toXiCosineOfPointingAngle) {
	fHistV0toXiCosineOfPointingAngle = new TH1F("fHistV0toXiCosineOfPointingAngle", "Cos. of V0 Ptng Angl / Xi vtx ;Cos(V0 Point. Angl / Xi vtx); Counts", 100, 0.99, 1.0);
	fListHistCascade->Add(fHistV0toXiCosineOfPointingAngle);
}


if(! fHistRapXi ){
	fHistRapXi  = new TH1F( "fHistRapXi" , "Rapidity of #Xi candidates (around the mass peak); y ; Counts", 20, -1.0, 1.0);
	fListHistCascade->Add(fHistRapXi);
}

if(! fHistRapOmega ){
	fHistRapOmega  = new TH1F( "fHistRapOmega" , "Rapidity of #Omega candidates (around the mass peak); y ; Counts", 20, -1.0, 1.0);
	fListHistCascade->Add(fHistRapOmega);
}

if(! fHistEtaXi ){
	fHistEtaXi  = new TH1F( "fHistEtaXi" , "Pseudo-rap. of #Xi candidates (around the mass peak) ; #eta ; Counts", 20, -1.0, 1.0);
	fListHistCascade->Add(fHistEtaXi);
}

if(! fHistThetaXi ){
	fHistThetaXi  = new TH1F( "fHistThetaXi" , "#theta of #Xi candidates (around the mass peak); #theta (deg) ; Counts", 180, 0., 180.0);
	fListHistCascade->Add(fHistThetaXi);
}

if(! fHistPhiXi ){
	fHistPhiXi  = new TH1F( "fHistPhiXi" , "#phi of #Xi candidates (around the mass peak); #phi (deg) ; Counts", 360, 0., 360.);
	fListHistCascade->Add(fHistPhiXi);
}


if(! f2dHistArmenteros) {
	f2dHistArmenteros = new TH2F( "f2dHistArmenteros", "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
	fListHistCascade->Add(f2dHistArmenteros);
}

//-------

if(! f2dHistEffMassLambdaVsEffMassXiMinus) {
	f2dHistEffMassLambdaVsEffMassXiMinus = new TH2F( "f2dHistEffMassLambdaVsEffMassXiMinus", "M_{#Lambda} Vs M_{#Xi^{-} candidates} ; Inv. M_{#Lambda^{0}} (GeV/c^{2}) ; M( #Lambda , #pi^{-} ) (GeV/c^{2})", 300, 1.1,1.13, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiMinus);
}

if(! f2dHistEffMassXiVsEffMassOmegaMinus) {
	f2dHistEffMassXiVsEffMassOmegaMinus = new TH2F( "f2dHistEffMassXiVsEffMassOmegaMinus", "M_{#Xi^{-} candidates} Vs M_{#Omega^{-} candidates} ; M( #Lambda , #pi^{-} ) (GeV/c^{2}) ; M( #Lambda , K^{-} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaMinus);
}

if(! f2dHistEffMassLambdaVsEffMassXiPlus) {
	f2dHistEffMassLambdaVsEffMassXiPlus = new TH2F( "f2dHistEffMassLambdaVsEffMassXiPlus", "M_{#Lambda} Vs M_{#Xi^{+} candidates} ; Inv. M_{#Lambda^{0}} (GeV/c^{2}) ; M( #Lambda , #pi^{+} ) (GeV/c^{2})", 300, 1.1,1.13, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistEffMassLambdaVsEffMassXiPlus);
}

if(! f2dHistEffMassXiVsEffMassOmegaPlus) {
	f2dHistEffMassXiVsEffMassOmegaPlus = new TH2F( "f2dHistEffMassXiVsEffMassOmegaPlus", "M_{#Xi^{+} candidates} Vs M_{#Omega^{+} candidates} ; M( #Lambda , #pi^{+} ) (GeV/c^{2}) ; M( #Lambda , K^{+} ) (GeV/c^{2})", 400, 1.2, 2.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistEffMassXiVsEffMassOmegaPlus);
}

//-------

if(! f2dHistXiRadiusVsEffMassXiMinus) {
	f2dHistXiRadiusVsEffMassXiMinus = new TH2F( "f2dHistXiRadiusVsEffMassXiMinus", "Transv. R_{Xi Decay} Vs M_{#Xi^{-} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiMinus);
}

if(! f2dHistXiRadiusVsEffMassXiPlus) {
	f2dHistXiRadiusVsEffMassXiPlus = new TH2F( "f2dHistXiRadiusVsEffMassXiPlus", "Transv. R_{Xi Decay} Vs M_{#Xi^{+} candidates}; r_{cascade} (cm); M( #Lambda , #pi^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 400, 1.2, 2.0);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassXiPlus);
}

if(! f2dHistXiRadiusVsEffMassOmegaMinus) {
	f2dHistXiRadiusVsEffMassOmegaMinus = new TH2F( "f2dHistXiRadiusVsEffMassOmegaMinus", "Transv. R_{Xi Decay} Vs M_{#Omega^{-} candidates}; r_{cascade} (cm); M( #Lambda , K^{-} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaMinus);
}

if(! f2dHistXiRadiusVsEffMassOmegaPlus) {
	f2dHistXiRadiusVsEffMassOmegaPlus = new TH2F( "f2dHistXiRadiusVsEffMassOmegaPlus", "Transv. R_{Xi Decay} Vs M_{#Omega^{+} candidates}; r_{cascade} (cm); M( #Lambda , K^{+} ) (GeV/c^{2}) ", 450, 0., 45.0, 500, 1.5, 2.5);
	fListHistCascade->Add(f2dHistXiRadiusVsEffMassOmegaPlus);
}

//------

if(! f2dHistTPCdEdxOfCascDghters){
        f2dHistTPCdEdxOfCascDghters = new TH2F( "f2dHistTPCdEdxOfCascDghters", "TPC dE/dx of the cascade daughters; charge x || #vec{p}_{TPC inner wall}(Casc. daughter) || (GeV/c); TPC signal (ADC) ", 2000, -10.0, 10.0, 450, 0., 900.);
	fListHistCascade->Add(f2dHistTPCdEdxOfCascDghters);
}



// Part 2 : Raw material for yield extraction -------
/*
if(! f3dHistXiPtVsEffMassVsYXiMinus) {
	f3dHistXiPtVsEffMassVsYXiMinus = new TH3F( "f3dHistXiPtVsEffMassVsYXiMinus", "Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{-} ) (GeV/c^{2}) ;Y_{#Xi} ", 100, 0., 10.0, 400, 1.2, 2.0, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYXiMinus);
}

if(! f3dHistXiPtVsEffMassVsYXiPlus) {
	f3dHistXiPtVsEffMassVsYXiPlus = new TH3F( "f3dHistXiPtVsEffMassVsYXiPlus", "Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}; Pt_{cascade} (GeV/c); M( #Lambda , #pi^{+} ) (GeV/c^{2}); Y_{#Xi}", 100, 0., 10.0, 400, 1.2, 2.0, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYXiPlus);
}

if(! f3dHistXiPtVsEffMassVsYOmegaMinus) {
	f3dHistXiPtVsEffMassVsYOmegaMinus = new TH3F( "f3dHistXiPtVsEffMassVsYOmegaMinus", "Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{-} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYOmegaMinus);
}

if(! f3dHistXiPtVsEffMassVsYOmegaPlus) {
	f3dHistXiPtVsEffMassVsYOmegaPlus = new TH3F( "f3dHistXiPtVsEffMassVsYOmegaPlus", "Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}; Pt_{cascade} (GeV/c); M( #Lambda , K^{+} ) (GeV/c^{2}); Y_{#Omega}", 100, 0., 10.0, 500, 1.5, 2.5, 44, -1.1,1.1);
	fListHistCascade->Add(f3dHistXiPtVsEffMassVsYOmegaPlus);
}
*/
//--
if(!fCFContCascadePIDXiMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 22;//44;
  lNbBinsPerVar[3] = 100; //250;
   
  
  fCFContCascadePIDXiMinus = new AliCFContainer("fCFContCascadePIDXiMinus","Pt_{cascade} Vs M_{#Xi^{-} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDXiMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDXiMinus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDXiMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDXiMinus->SetBinLimits(3, 0.0, 20000.0  );    // SPD tracklets Multiplicity
  else
        fCFContCascadePIDXiMinus->SetBinLimits(3, 0.0, 250.0  );     // SPD tracklets Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDXiMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDXiMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDXiMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDXiMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDXiMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDXiMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDXiMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDXiMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDXiMinus->SetVarTitle(1, "M( #Lambda , #pi^{-} ) (GeV/c^{2})");
  fCFContCascadePIDXiMinus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDXiMinus->SetVarTitle(3, "SPD tracklets Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDXiMinus);
  
}

if (!fCFContCascadePIDXiPlus) {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 400;
  lNbBinsPerVar[2] = 22; //44;
  lNbBinsPerVar[3] = 100; //250;
   
  
  fCFContCascadePIDXiPlus = new AliCFContainer("fCFContCascadePIDXiPlus","Pt_{cascade} Vs M_{#Xi^{+} candidates} Vs Y_{#Xi}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDXiPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDXiPlus->SetBinLimits(1,   1.2  ,   2.0 );	// Xi Effective mass
  fCFContCascadePIDXiPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDXiPlus->SetBinLimits(3, 0.0, 20000.0  );    // SPD tracklets Multiplicity
  else
        fCFContCascadePIDXiPlus->SetBinLimits(3, 0.0, 250.0  );     // SPD tracklets Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDXiPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDXiPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDXiPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDXiPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDXiPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDXiPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDXiPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDXiPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDXiPlus->SetVarTitle(1, "M( #Lambda , #pi^{+} ) (GeV/c^{2})");
  fCFContCascadePIDXiPlus->SetVarTitle(2, "Y_{#Xi}");
  fCFContCascadePIDXiPlus->SetVarTitle(3, "SPD tracklets Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDXiPlus);
  
}


if(!fCFContCascadePIDOmegaMinus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 22; //44;
  lNbBinsPerVar[3] = 100; //250;
   
  
  fCFContCascadePIDOmegaMinus = new AliCFContainer("fCFContCascadePIDOmegaMinus","Pt_{cascade} Vs M_{#Omega^{-} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDOmegaMinus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDOmegaMinus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDOmegaMinus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDOmegaMinus->SetBinLimits(3, 0.0, 20000.0  );    //SPD tracklets Multiplicity
  else
        fCFContCascadePIDOmegaMinus->SetBinLimits(3, 0.0, 250.0  );     // SPD tracklets Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDOmegaMinus->SetStepTitle(0, "No PID");
  fCFContCascadePIDOmegaMinus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDOmegaMinus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDOmegaMinus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDOmegaMinus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDOmegaMinus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDOmegaMinus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDOmegaMinus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDOmegaMinus->SetVarTitle(1, "M( #Lambda , K^{-} ) (GeV/c^{2})");
  fCFContCascadePIDOmegaMinus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDOmegaMinus->SetVarTitle(3, "SPD tracklets Multiplicity");
  
  fListHistCascade->Add(fCFContCascadePIDOmegaMinus);
  
}

if(!fCFContCascadePIDOmegaPlus)  {
  const	Int_t  lNbSteps      =  7 ;
  const Int_t  lNbVariables  =  4 ;

  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[4] = {0};
  lNbBinsPerVar[0] = 100;
  lNbBinsPerVar[1] = 500;
  lNbBinsPerVar[2] = 22; //44;
  lNbBinsPerVar[3] = 100; // 250
   
  
  fCFContCascadePIDOmegaPlus = new AliCFContainer("fCFContCascadePIDOmegaPlus","Pt_{cascade} Vs M_{#Omega^{+} candidates} Vs Y_{#Omega}", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits (valid  for v4-18-10-AN)
  fCFContCascadePIDOmegaPlus->SetBinLimits(0,   0.0  ,  10.0 );	// Pt(Cascade)
  fCFContCascadePIDOmegaPlus->SetBinLimits(1,   1.5  ,   2.5 );	// Omega Effective mass
  fCFContCascadePIDOmegaPlus->SetBinLimits(2,  -1.1  ,   1.1 );	// Rapidity
  if(fCollidingSystems) 
        fCFContCascadePIDOmegaPlus->SetBinLimits(3, 0.0, 20000.0  );    // SPD tracklets Multiplicity
  else
        fCFContCascadePIDOmegaPlus->SetBinLimits(3, 0.0, 250.0  );     // SPD tracklets Multiplicity
  
  // Setting the step title : one per PID case
  fCFContCascadePIDOmegaPlus->SetStepTitle(0, "No PID");
  fCFContCascadePIDOmegaPlus->SetStepTitle(1, "TPC PID / 4-#sigma cut on Bachelor track");
  fCFContCascadePIDOmegaPlus->SetStepTitle(2, "TPC PID / 4-#sigma cut on Bachelor+Baryon tracks");
  fCFContCascadePIDOmegaPlus->SetStepTitle(3, "TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks");
  fCFContCascadePIDOmegaPlus->SetStepTitle(4, "Comb. PID / Bachelor");
  fCFContCascadePIDOmegaPlus->SetStepTitle(5, "Comb. PID / Bachelor+Baryon");
  fCFContCascadePIDOmegaPlus->SetStepTitle(6, "Comb. PID / Bachelor+Baryon+Meson");
  
  // Setting the variable title, per axis
  fCFContCascadePIDOmegaPlus->SetVarTitle(0, "Pt_{cascade} (GeV/c)");
  fCFContCascadePIDOmegaPlus->SetVarTitle(1, "M( #Lambda , K^{+} ) (GeV/c^{2})");
  fCFContCascadePIDOmegaPlus->SetVarTitle(2, "Y_{#Omega}");
  fCFContCascadePIDOmegaPlus->SetVarTitle(3, "SPD tracklets Multiplicity"); 
  
  fListHistCascade->Add(fCFContCascadePIDOmegaPlus);
  
}


// Part 3 : Towards the optimisation of topological selections -------
if(! fCFContCascadeCuts&&fUseCFContCascadeCuts) {
   
	// Container meant to store all the relevant distributions corresponding to the cut variables.
	// So far, 20 variables have been identified.
	// The following will be done in quite a brut force way ... 
        //          - Define a user binning to have less bins in each dimension? 
  const	Int_t  lNbSteps      =  4 ;
  const Int_t  lNbVariables  =  20 ;
  
  //array for the number of bins in each dimension :
  Int_t lNbBinsPerVar[20] = {0};
  lNbBinsPerVar[0]  = 100;
  lNbBinsPerVar[1]  = 125;
  lNbBinsPerVar[2]  = 100;
  lNbBinsPerVar[3]  = 440;
  lNbBinsPerVar[4]  = 30;
  lNbBinsPerVar[5]  = 50;
  
  lNbBinsPerVar[6]  = 100;
  lNbBinsPerVar[7]  = 40;
  lNbBinsPerVar[8]  = 100;
  lNbBinsPerVar[9]  = 25;
  lNbBinsPerVar[10] = 25;
  
  lNbBinsPerVar[11] = 75;  // 2-MeV/c2 bins
  lNbBinsPerVar[12] = 60;  // 2-MeV/c2 bins
  
  lNbBinsPerVar[13] = 100;
  
  lNbBinsPerVar[14] = 44; // 0.05 in rapidity units
  lNbBinsPerVar[15] = 44; // 0.05 in rapidity units
  
  lNbBinsPerVar[16] = 20;
 
  lNbBinsPerVar[17] = 50;
  lNbBinsPerVar[18] = 100;
  lNbBinsPerVar[19] = 24;
    
 fCFContCascadeCuts = new AliCFContainer("fCFContCascadeCuts","Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar );
  
  
  //setting the bin limits 
   
  //0
  fCFContCascadeCuts->SetBinLimits(0, 0., 2.);                 // DcaXiDaughters : 0.0 to 2.0
  //1
  fCFContCascadeCuts->SetBinLimits(1, 0., .5);                 // DcaBachToPrimVertexXi : 0.0 to 0.5
  //2 
  fCFContCascadeCuts->SetBinLimits(2, .99, 1.);                // XiCosineOfPointingAngle : 0.99 to 1.0	
  //3
  fCFContCascadeCuts->SetBinLimits(3, 0., 5.);                 // XiRadius : 0.0 to 5.0	
  //4
  fCFContCascadeCuts->SetBinLimits(4, 1.1, 1.13 );             // InvMassLambdaAsCascDghter
  //5
  fCFContCascadeCuts->SetBinLimits(5, 0., 2.);                 // DcaV0DaughtersXi : 0.0 to 2.0	
  //6 
  fCFContCascadeCuts->SetBinLimits(6, .99, 1.);                // V0CosineOfPointingAngleXi : 0.99 to 1.0	
  //7
  fCFContCascadeCuts->SetBinLimits(7, 0., 20.);                // V0RadiusXi : 0.0 to 20.0	
  //8
  fCFContCascadeCuts->SetBinLimits(8, 0., .4);                 // DcaV0ToPrimVertexXi : 0.0 to 0.4	
  //9
  fCFContCascadeCuts->SetBinLimits(9, 0., 0.25);               // DcaPosToPrimVertexXi : 0.0 to 0.25	
  //10
  fCFContCascadeCuts->SetBinLimits(10, 0., 0.25);              // DcaPosToPrimVertexXi : 0.0 to 0.25	
  //11
  fCFContCascadeCuts->SetBinLimits(11, 1.25, 1.40);            // InvMassXi
  fCFContCascadeCuts->SetBinLimits(12, 1.62, 1.74);            // InvMassOmega
  fCFContCascadeCuts->SetBinLimits(13, 0.0, 10.0);             // XiTransvMom 
  fCFContCascadeCuts->SetBinLimits(14, -1.1, 1.1);             // Y(Xi)
  fCFContCascadeCuts->SetBinLimits(15, -1.1, 1.1);             // Y(Omega)
  fCFContCascadeCuts->SetBinLimits(16, -10.0, 10.0);           // BestPrimaryVtxPosZ
  if (fCollidingSystems) {
    fCFContCascadeCuts->SetBinLimits(17, 0.0, 10000.0);        // TPCrefitTrackMultiplicity
    fCFContCascadeCuts->SetBinLimits(18, 0.0, 10000.0);        // SPDTrackletsMultiplicity
  } else {  
    fCFContCascadeCuts->SetBinLimits(17, 0.0, 250.0);          // TPCrefitTrackMultiplicity
    fCFContCascadeCuts->SetBinLimits(18, 0.0, 200.0);          // SPDTrackletsMultiplicity
  }
  fCFContCascadeCuts->SetBinLimits(19, 25.0, 165.0);           // BachTPCClusters
  
  // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
  fCFContCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates");
  fCFContCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
  fCFContCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates");
  fCFContCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
  
  // Setting the variable title, per axis
  // fCFContCascadeCuts->SetVarTitle(40,  "Chi2Xi");
  fCFContCascadeCuts->SetVarTitle(0,  "Dca(XiDaughters) (cm)");
  fCFContCascadeCuts->SetVarTitle(1,  "Dca(Bach/PrimVertex) (cm)");
  fCFContCascadeCuts->SetVarTitle(2,  "cos(Xi pointing angle)");
  fCFContCascadeCuts->SetVarTitle(3,  "R_{2d}(Xi decay) (cm)");
  fCFContCascadeCuts->SetVarTitle(4,  "M_{#Lambda}(As Casc Dghter) (GeV/c^{2})");
  fCFContCascadeCuts->SetVarTitle(5,  "Dca(V0 Daughters) in Xi (cm)");
  
  fCFContCascadeCuts->SetVarTitle(6,  "cos(V0 pointing Angle) in Casc");
  fCFContCascadeCuts->SetVarTitle(7,  "R_{2d}(V0 decay) (cm)");
  fCFContCascadeCuts->SetVarTitle(8,  "Dca(V0/PrimVertex) (cm)");
  fCFContCascadeCuts->SetVarTitle(9,  "Dca(Pos/PrimVertex) (cm)");
  fCFContCascadeCuts->SetVarTitle(10, "Dca(Neg/PrimVertex) (cm)");
  
  fCFContCascadeCuts->SetVarTitle(11, "Inv. Mass(Xi) (GeV/c^{2})");
  fCFContCascadeCuts->SetVarTitle(12, "Inv. Mass(Omega) (GeV/c^{2})");
  
  fCFContCascadeCuts->SetVarTitle(13, "pt(Casc.) (GeV/c)");
  
  fCFContCascadeCuts->SetVarTitle(14, "Y(Xi)");
  fCFContCascadeCuts->SetVarTitle(15, "Y(Omega)");
  
  fCFContCascadeCuts->SetVarTitle(16, "Z-position(BestPrimVtx) (cm)");
  
  fCFContCascadeCuts->SetVarTitle(17, "TPCrefit track Multiplicity");
  fCFContCascadeCuts->SetVarTitle(18, "SPD tracklets Multiplicity");
  fCFContCascadeCuts->SetVarTitle(19, "Bach.TPC Clusters");
  
  fListHistCascade->Add(fCFContCascadeCuts);
}

  fV0Ampl = new TH1F("fV0Ampl","",500,0.,30000);
  fListHistCascade->Add(fV0Ampl);


  fHistDcaXiDaughtersvsInvMass = new TH2F( "fHistDcaXiDaughtersvsInvMass",  "DCA between Xi Daughters; DCA (cm) ; Number of Cascades", 100, 0., 0.5,400,1.2,2.0);
        fListHistCascade->Add(fHistDcaXiDaughtersvsInvMass); 
  fHistDcaBachToPrimVertexvsInvMass = new TH2F("fHistDcaBachToPrimVertexvsInvMass", "DCA of Bach. to Prim. Vertex;DCA (cm);Number of Cascades", 250, 0., 0.25,400,1.2,2.0);
        fListHistCascade->Add(fHistDcaBachToPrimVertexvsInvMass); 
  fHistXiCosineOfPointingAnglevsInvMass= new TH2F("fHistXiCosineOfPointingAnglevsInvMass", "Cosine of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis", 200, 0.99, 1.0,400,1.2,2.0);
        fListHistCascade->Add(fHistXiCosineOfPointingAnglevsInvMass); 
  fHistMassLambdaAsCascDghtervsInvMass= new TH2F("fHistMassLambdaAsCascDghtervsInvMass","#Lambda associated to Casc. candidates;Eff. Mass (GeV/c^{2});Counts", 300,1.00,1.3,400,1.2,2.0);
    fListHistCascade->Add(fHistMassLambdaAsCascDghtervsInvMass);
  fHistDcaV0DaughtersXivsInvMass = new TH2F("fHistDcaV0DaughtersXivsInvMass", "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 0.6,400,1.2,2.0);
        fListHistCascade->Add(fHistDcaV0DaughtersXivsInvMass);
  fHistDcaV0ToPrimVertexXivsInvMass = new TH2F("fHistDcaV0ToPrimVertexXivsInvMass", "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.,400,1.2,2.0);
        fListHistCascade->Add(fHistDcaV0ToPrimVertexXivsInvMass);



// Part 4 : Angular correlation study -------
/*
if(! fHnSpAngularCorrXiMinus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
	
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.30};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.34};
   fHnSpAngularCorrXiMinus = new THnSparseF("fHnSpAngularCorrXiMinus", "Angular Correlation for #Xi^{-}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrXiMinus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrXiMinus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrXiMinus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrXiMinus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrXiMinus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrXiMinus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrXiMinus);
}

if(! fHnSpAngularCorrXiPlus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.30};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.34};
   fHnSpAngularCorrXiPlus = new THnSparseF("fHnSpAngularCorrXiPlus", "Angular Correlation for #Xi^{+}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrXiPlus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrXiPlus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrXiPlus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrXiPlus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrXiPlus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrXiPlus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrXiPlus);
}

if(! fHnSpAngularCorrOmegaMinus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
	
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.65};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.69};
   fHnSpAngularCorrOmegaMinus = new THnSparseF("fHnSpAngularCorrOmegaMinus", "Angular Correlation for #Omega^{-}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrOmegaMinus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrOmegaMinus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrOmegaMinus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrOmegaMinus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrOmegaMinus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrOmegaMinus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrOmegaMinus);
}

if(! fHnSpAngularCorrOmegaPlus){
	// Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks
	// Delta Phi  = 360 bins de -180., 180.
	// Delta Eta  = 120 bins de -3.0, 3.0
	// Pt Cascade = 100 bins de 0., 10.0,
	// Pt track = 150 bins de 0., 15.0
   Int_t    bins[5] = { 360, 120, 100, 150, 40};
   Double_t xmin[5] = {-50., -3.,  0.,  0., 1.65};
   Double_t xmax[5] = { 310., 3., 10., 15., 1.69};
   fHnSpAngularCorrOmegaPlus = new THnSparseF("fHnSpAngularCorrOmegaPlus", "Angular Correlation for #Omega^{+}:", 5, bins, xmin, xmax);
	fHnSpAngularCorrOmegaPlus->GetAxis(0)->SetTitle(" #Delta#phi(Casc,Track) (deg)");
	fHnSpAngularCorrOmegaPlus->GetAxis(1)->SetTitle(" #Delta#eta(Casc,Track)");
	fHnSpAngularCorrOmegaPlus->GetAxis(2)->SetTitle(" Pt_{Casc} (GeV/c)");
	fHnSpAngularCorrOmegaPlus->GetAxis(3)->SetTitle(" Pt_{any track} (GeV/c)");
	fHnSpAngularCorrOmegaPlus->GetAxis(4)->SetTitle(" Eff. Inv Mass (GeV/c^{2})");
	fHnSpAngularCorrOmegaPlus->Sumw2();
   fListHistCascade->Add(fHnSpAngularCorrOmegaPlus);
}
*/

PostData(1, fListHistCascade);
/* PostData(2, fPaveTextBookKeeping);*/
}// end UserCreateOutputObjects






//________________________________________________________________________
void AliAnalysisTaskCheckCascadePbPb::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  AliESDEvent *lESDevent = 0x0;
  AliAODEvent *lAODevent = 0x0;
  Int_t    ncascades                      = -1;
  Int_t    nTrackMultiplicity             = -1;
  Int_t    nTrackWithTPCrefitMultiplicity =  0;

  Short_t  lStatusTrackingPrimVtx         = -2;
  Double_t lTrkgPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
  Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
  Double_t lMagneticField                 = -10.;

  // Connect to the InputEvent	
  // After these lines, we should have an ESD/AOD event + the number of cascades in it.
		
  if (fAnalysisType == "ESD") {
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
    }
        
    fHistCascadeMultiplicityBeforeEvSel->Fill ( lESDevent->GetNumberOfCascades() );  
        
    //-------------------------------------------------

    // 0 - Trigger managment
    // NOTE : Check the availability of the proper trigger  --> MN Physics selection moved to runProof  macro

    //-------------------------------------------------

    // 1 - Cascade vertexer (ESD)
        
    if (fkRerunV0CascVertexers) { // FIXME : relaunch V0 and Cascade vertexers
//      if (fAnalysisType == "ESD" ) {
//        lESDevent->ResetCascades();
//        lESDevent->ResetV0s();
//                 
//        AliV0vertexer lV0vtxer;
//        AliCascadeVertexer lCascVtxer;
//                 
//        lV0vtxer.SetDefaultCuts(fV0Sels);
//        lCascVtxer.SetDefaultCuts(fCascSels);
//                 
//        lV0vtxer.Tracks2V0vertices(lESDevent);
//        lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
//      }
    }// end if(RelaunchV0CascVertexers)
        
    //-------------------------------------------------
    ncascades                      = lESDevent->GetNumberOfCascades();
    nTrackWithTPCrefitMultiplicity = DoESDTrackWithTPCrefitMultiplicity(lESDevent);  

  }//if (fAnalysisType == "ESD")
 
  // MN: variable for centrality selection
  AliESDVZERO* esdV0 = lESDevent->GetVZEROData();
  Float_t multV0A=esdV0->GetMTotV0A();
  Float_t multV0C=esdV0->GetMTotV0C();
//  fV0Ampl->Fill(multV0A+multV0C);

  AliCentrality *centrality = lESDevent->GetCentrality();
  if (!(centrality->IsEventInCentralityClass(fCentrLowLim,fCentrUpLim,fCentrEstimator.Data()))) {
    PostData(1, fListHistCascade);
    return; 
  }

  fV0Ampl->Fill(multV0A+multV0C);
  
  if (fAnalysisType == "AOD") {  
    lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() ); 
    if (!lAODevent) {
      AliWarning("ERROR: lAODevent not available \n");
      return;
    }
    ncascades                      = lAODevent->GetNumberOfCascades();
    nTrackWithTPCrefitMultiplicity = -1;
        
  }

  // For AOD or ESD ...
  nTrackMultiplicity = (InputEvent())->GetNumberOfTracks();

  //-------------------------------------------------

  fHistTrackMultiplicityForCentrEvt         ->Fill( nTrackMultiplicity             );
  fHistTPCrefitTrackMultiplicityForCentrEvt ->Fill( nTrackWithTPCrefitMultiplicity );
  fHistCascadeMultiplicityForCentrEvt       ->Fill( ncascades                      );

  // ---------------------------------------------------------------

  // I - Global characteristics of the events + general histos (filled for any selected events and/or for the analysed events)

       // - I.Step 1 : Characteristics of the event : prim. Vtx + magnetic field (ESD)
       //-------------

  if (fAnalysisType == "ESD") {
    const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    // get the vtx stored in ESD found with tracks
    lPrimaryTrackingESDVtx->GetXYZ( lTrkgPrimaryVtxPos );
        
    const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();	
    // get the best primary vertex available for the event
    // As done in AliCascadeVertexer, we keep the one which is the best one available.
    // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
    // This one will be used for next calculations (DCA essentially)
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    lStatusTrackingPrimVtx  = lPrimaryTrackingESDVtx->GetStatus();

    // FIXME : quality cut on the z-position of the prim vertex.
    if (fkQualityCutZprimVtxPos) {
      if (TMath::Abs(lBestPrimaryVtxPos[2]) > fVtxRange ) { 
        AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
        PostData(1, fListHistCascade); 
        return; 
      }
    }
        
    fHistCascadeMultiplicityForTrigEvtAndZprimVtx->Fill( ncascades  );
        
    // FIXME : remove TPC-only primary vertex : retain only events with tracking + SPD vertex
    if (fkQualityCutNoTPConlyPrimVtx) {
      const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
      if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus() ){
        AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
        PostData(1, fListHistCascade); 
        return;
      }
    }
        
        // NOTE : For older evts
        
        // As previously done in AliCascadeVertexer, we keep, between both retrieved vertices (SPD or Tracking) 
        // the one which is the best one available.
        // This one will be used for next calculations (DCA essentially)
        // At that time, the TPC-only primary vertex was not considered
        
        
    lMagneticField = lESDevent->GetMagneticField( );
    // FIXME if(TMath::Abs(lMagneticField ) < 10e-6) continue;
        
  }// end if(ESD)
        
  if (fAnalysisType == "AOD") {
    // To be developed
    const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();	
    // get the best primary vertex available for the event
    // We may keep the one which is the best one available = GetVertex(0)
    // Pb with pile-up to expect
    // This one will be used for next calculations (DCA essentially)
    lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
        
    lStatusTrackingPrimVtx  = -1;
    lTrkgPrimaryVtxPos[0]   = -100.0;
    lTrkgPrimaryVtxPos[1]   = -100.0;
    lTrkgPrimaryVtxPos[2]   = -100.0;
    lMagneticField = 0.;
  }

       // - I.Step 2 : Filling histos that characterize the selected event : x,y,z prim. Vtx distrib. (ESD)
       //-------------

  fHistCascadeMultiplicityForSelEvt ->Fill( ncascades );
  fHistPosBestPrimaryVtxXForSelEvt  ->Fill( lBestPrimaryVtxPos[0] );
  fHistPosBestPrimaryVtxYForSelEvt  ->Fill( lBestPrimaryVtxPos[1] );
  fHistPosBestPrimaryVtxZForSelEvt  ->Fill( lBestPrimaryVtxPos[2] );
  
  // ---------------------------------------------------------------
  // II - Calcultaion Part dedicated to Xi vertices
  
  for (Int_t iXi = 0; iXi < ncascades; iXi++) {// This is the begining of the Cascade loop (ESD or AOD)
	   
    // -------------------------------------
    // II.Init - Initialisation of the local variables that will be needed for ESD/AOD

  
        // - 0th part of initialisation : around primary vertex ...
	
    Double_t lTrkgPrimaryVtxRadius3D = -500.0;
    Double_t lBestPrimaryVtxRadius3D = -500.0;

        // - 1st part of initialisation : variables needed to store AliESDCascade data members
    Double_t lEffMassXi      = 0. ;
    Double_t lChi2Xi         = -1. ;
    Double_t lDcaXiDaughters = -1. ;
    Double_t lXiCosineOfPointingAngle = -1. ;
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
    Double_t lXiRadius = -1000. ;
        
        // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
    Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
    Int_t    lNegTPCClusters    = -1; // For ESD only ...
    Int_t    lBachTPCClusters   = -1; // For ESD only ...
        
    Double_t lInnerWallMomCascDghters[3] = {-100., -100., -100.};
    Double_t lTPCSignalCascDghters   [3] = {-100., -100., -100.};
       
        
        // - 3rd part of initialisation : about V0 part in cascades
    Double_t lInvMassLambdaAsCascDghter = 0.;
    Double_t lV0Chi2Xi         = -1. ;
    Double_t lDcaV0DaughtersXi = -1.;
		
    Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
    Double_t lDcaPosToPrimVertexXi  = -1.;
    Double_t lDcaNegToPrimVertexXi  = -1.;
    Double_t lV0CosineOfPointingAngleXi = -1. ;
    Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
    Double_t lV0RadiusXi = -1000.0;
    Double_t lV0quality  = 0.;

	
	// - 4th part of initialisation : Effective masses
    Double_t lInvMassXiMinus    = 0.;
    Double_t lInvMassXiPlus     = 0.;
    Double_t lInvMassOmegaMinus = 0.;
    Double_t lInvMassOmegaPlus  = 0.;
  
	// - 5th part of initialisation : PID treatment
    Bool_t   lIsPosInXiProton      = kFALSE;
    Bool_t   lIsPosInXiPion        = kFALSE;
    Bool_t   lIsPosInOmegaProton   = kFALSE;
    Bool_t   lIsPosInOmegaPion     = kFALSE;
			
    Bool_t   lIsNegInXiProton      = kFALSE;
    Bool_t   lIsNegInXiPion        = kFALSE;
    Bool_t   lIsNegInOmegaProton   = kFALSE;
    Bool_t   lIsNegInOmegaPion     = kFALSE;

    Bool_t   lIsBachelorKaon       = kFALSE;
    Bool_t   lIsBachelorPion       = kFALSE; 
	
    Bool_t   lIsBachelorKaonForTPC = kFALSE; // For ESD only ...//FIXME : wait for availability in AOD
    Bool_t   lIsBachelorPionForTPC = kFALSE; // For ESD only ...
    Bool_t   lIsNegPionForTPC      = kFALSE; // For ESD only ...
    Bool_t   lIsPosPionForTPC      = kFALSE; // For ESD only ...
    Bool_t   lIsNegProtonForTPC    = kFALSE; // For ESD only ...
    Bool_t   lIsPosProtonForTPC    = kFALSE; // For ESD only ...

	// - 6th part of initialisation : extra info for QA
    Double_t lXiMomX       = 0., lXiMomY = 0., lXiMomZ = 0.;
    Double_t lXiTransvMom  = 0. ;
    Double_t lXiTotMom     = 0. ;
		
    Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
    Double_t lBachTransvMom  = 0.;
    Double_t lBachTotMom     = 0.;
	
    Short_t  lChargeXi = -2;
    Double_t lV0toXiCosineOfPointingAngle = 0. ;
	
    Double_t lRapXi   = -20.0, lRapOmega = -20.0,  lEta = -20.0, lTheta = 360., lPhi = 720. ;
    Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
	
  	// - 7th part of initialisation : variables for the AliCFContainer dedicated to cascade cut optmisiation
    Int_t    lSPDTrackletsMultiplicity = -1;
  
  	// - 8th part of initialisation : variables needed for Angular correlations
    TVector3 lTVect3MomXi(0.,0.,0.);
    Int_t    lArrTrackID[3] = {-1, -1, -1};

	  
    if (fAnalysisType == "ESD") { 
  
    // -------------------------------------
    // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
  
      AliESDcascade *xi = lESDevent->GetCascade(iXi);
      if (!xi) continue;
        
        
                // - II.Step 1 : around primary vertex
                //-------------
      lTrkgPrimaryVtxRadius3D = TMath::Sqrt(  lTrkgPrimaryVtxPos[0] * lTrkgPrimaryVtxPos[0] +
                                              lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
                                              lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );

      lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
                                              lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
                                              lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );


	
		// - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)	
		//-------------
      lV0quality = 0.;
      xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay

      lEffMassXi  	        = xi->GetEffMassXi();
      lChi2Xi 			= xi->GetChi2Xi();
      lDcaXiDaughters		= xi->GetDcaXiDaughters();
      lXiCosineOfPointingAngle 	= xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                       lBestPrimaryVtxPos[1],
                                                                       lBestPrimaryVtxPos[2] );
		// Take care : the best available vertex should be used (like in AliCascadeVertexer)
	
      xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
      lXiRadius        	= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
		
		

		// - II.Step 3 : around the tracks : Bach + V0 (ESD)
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
      UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
      UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
      UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                // Care track label can be negative in MC production (linked with the track quality)
                // However = normally, not the case for track index ...
        
                // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
      if (lBachIdx == lIdxNegXi) {
        AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
      }
      if (lBachIdx == lIdxPosXi) {
        AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
      }
        
      AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
      AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
      AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
      if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
        AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
        continue;
      }
        
        
      lPosTPCClusters   = pTrackXi->GetTPCNcls();
      lNegTPCClusters   = nTrackXi->GetTPCNcls();
      lBachTPCClusters  = bachTrackXi->GetTPCNcls();
      
                // FIXME : rejection of a poor quality tracks
      if (fkQualityCutTPCrefit) {
                // 1 - Poor quality related to TPCrefit
        ULong_t pStatus    = pTrackXi->GetStatus();
        ULong_t nStatus    = nTrackXi->GetStatus();
        ULong_t bachStatus = bachTrackXi->GetStatus();
        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
      }
      if (fkQualityCut80TPCcls) {
                // 2 - Poor quality related to TPC clusters
        if (lPosTPCClusters  < 80) { AliWarning("Pb / V0 Pos. track has less than 80 TPC clusters ... continue!"); continue; }
        if (lNegTPCClusters  < 80) { AliWarning("Pb / V0 Neg. track has less than 80 TPC clusters ... continue!"); continue; }
        if(lBachTPCClusters < 80) { AliWarning("Pb / Bach.   track has less than 80 TPC clusters ... continue!"); continue; }
      }
        
      const AliExternalTrackParam *pExtTrack    = pTrackXi    ->GetInnerParam();
      const AliExternalTrackParam *nExtTrack    = nTrackXi    ->GetInnerParam();
      const AliExternalTrackParam *bachExtTrack = bachTrackXi ->GetInnerParam();
       
      if (pExtTrack) {
        lInnerWallMomCascDghters[0] = pExtTrack ->GetP() * pExtTrack ->Charge();
        lTPCSignalCascDghters   [0] = pTrackXi  ->GetTPCsignal();
      }
      if (nExtTrack) {
        lInnerWallMomCascDghters[1] = nExtTrack ->GetP() * nExtTrack ->Charge();
        lTPCSignalCascDghters   [1] = nTrackXi  ->GetTPCsignal();
      }
      if (bachExtTrack) {
        lInnerWallMomCascDghters[2] = bachExtTrack ->GetP() * bachExtTrack ->Charge();
        lTPCSignalCascDghters   [2] = bachTrackXi  ->GetTPCsignal();
      }


      lInvMassLambdaAsCascDghter	= xi->GetEffMass();
		// This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
      lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters(); 
      lV0Chi2Xi 			= xi->GetChi2V0();
	
      lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
									  lBestPrimaryVtxPos[1],
									  lBestPrimaryVtxPos[2] );

      lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0], 
						    lBestPrimaryVtxPos[1], 
						    lBestPrimaryVtxPos[2] );
		
      lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0], 
						     		lBestPrimaryVtxPos[1], 
						     		lMagneticField  ) ); 
	        // Note : AliExternalTrackParam::GetD returns an algebraic value ...
		
      xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
      lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
	
      lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
						 		lBestPrimaryVtxPos[1], 
						 		lMagneticField  )     ); 
	
      lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
					    			lBestPrimaryVtxPos[1], 
					    			lMagneticField  )     ); 
		
		// - II.Step 3' : extra-selection for cascade candidates
		// Towards optimisation of AA selection
        // FIXME
      if (fkExtraSelections) {
       // if(lChi2Xi > 2000) continue; // in AliCascadeVertexer   
       // if(lV0Chi2Xi > 2000) continue; // in AliV0vertexer
		
        if (lDcaXiDaughters > 0.4) continue; // in AliCascadeVertexer 
        if (lXiCosineOfPointingAngle < 0.999 ) continue; // in AliCascadeVertexer
        if (lDcaV0ToPrimVertexXi < 0.05) continue; // in AliCascadeVertexer
        if (lDcaBachToPrimVertexXi < 0.03) continue; // in AliCascadeVertexer
////      if (TMath::Abs(lInvMassLambdaAsCascDghter-1.11568) > 0.006 ) continue;  // in AliCascadeVertexer 

        if (lDcaV0DaughtersXi > 1.) continue; // in AliV0vertexer 
        if (lV0CosineOfPointingAngleXi < 0.998) continue; // in AliV0vertexer
        if (lDcaPosToPrimVertexXi < 0.1) continue; // in AliV0vertexer
        if (lDcaNegToPrimVertexXi < 0.1) continue; // in AliV0vertexer


//	  if(lXiRadius < 1.0) continue; // in AliCascadeVertexer
//	  if(lXiRadius > 100) continue; // in AliCascadeVertexer
//	  if(lV0RadiusXi < 1.0) continue; // in AliV0vertexer  
//	  if(lV0RadiusXi > 100) continue; // in AliV0vertexer

      }
	
	
	
		// - II.Step 4 : around effective masses (ESD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

	
      if ( bachTrackXi->Charge() < 0 )	{
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3312); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3312 = Xi-
        lInvMassXiMinus = xi->GetEffMassXi();
		
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3334); 	
			// Calculate the effective mass of the Xi- candidate. 
			// pdg code 3334 = Omega-
        lInvMassOmegaMinus = xi->GetEffMassXi();
					
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
      }// end if negative bachelor
	
	
      if ( bachTrackXi->Charge() >  0 ) {
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , -3312); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3312 = Xi+
        lInvMassXiPlus = xi->GetEffMassXi();
		
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , -3334); 	
			// Calculate the effective mass of the Xi+ candidate. 
			// pdg code -3334  = Omega+
        lInvMassOmegaPlus = xi->GetEffMassXi();
		
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
      }// end if positive bachelor
	
	
	
		// - II.Step 5 : PID on the daughter tracks
		//-------------
	
	// A - Combined PID
	// Reasonable guess for the priors for the cascade track sample (e-, mu, pi, K, p)
      Double_t lPriorsGuessXi[5]    = {0, 0, 2, 0, 1};
      Double_t lPriorsGuessOmega[5] = {0, 0, 1, 1, 1};
	
	// Combined VO-positive-daughter PID
      AliPID pPidXi;		pPidXi.SetPriors(    lPriorsGuessXi    );
      AliPID pPidOmega;	pPidOmega.SetPriors( lPriorsGuessOmega );
		
      if ( pTrackXi->IsOn(AliESDtrack::kESDpid) ) {  // Combined PID exists
        Double_t r[10] = {0.}; pTrackXi->GetESDpid(r);
        pPidXi.SetProbabilities(r);
        pPidOmega.SetProbabilities(r);
		
		// Check if the V0 positive track is a proton (case for Xi-)
        Double_t pproton = pPidXi.GetProbability(AliPID::kProton);
        if (pproton > pPidXi.GetProbability(AliPID::kElectron) &&
            pproton > pPidXi.GetProbability(AliPID::kMuon)     &&
            pproton > pPidXi.GetProbability(AliPID::kPion)     &&
            pproton > pPidXi.GetProbability(AliPID::kKaon)     )     lIsPosInXiProton = kTRUE;
		
		// Check if the V0 positive track is a pi+ (case for Xi+)
        Double_t ppion = pPidXi.GetProbability(AliPID::kPion);
	if (ppion > pPidXi.GetProbability(AliPID::kElectron) &&
	    ppion > pPidXi.GetProbability(AliPID::kMuon)     &&
	    ppion > pPidXi.GetProbability(AliPID::kKaon)     &&
	    ppion > pPidXi.GetProbability(AliPID::kProton)   )     lIsPosInXiPion = kTRUE;
		
		
		// Check if the V0 positive track is a proton (case for Omega-)
	pproton = 0.;
	pproton = pPidOmega.GetProbability(AliPID::kProton);
	if (pproton > pPidOmega.GetProbability(AliPID::kElectron) &&
	    pproton > pPidOmega.GetProbability(AliPID::kMuon)     &&
	    pproton > pPidOmega.GetProbability(AliPID::kPion)     &&
	    pproton > pPidOmega.GetProbability(AliPID::kKaon)     )  lIsPosInOmegaProton = kTRUE;
		
		// Check if the V0 positive track is a pi+ (case for Omega+)
	ppion = 0.;
	ppion = pPidOmega.GetProbability(AliPID::kPion);
	if (ppion > pPidOmega.GetProbability(AliPID::kElectron) &&
	    ppion > pPidOmega.GetProbability(AliPID::kMuon)     &&
	    ppion > pPidOmega.GetProbability(AliPID::kKaon)     &&
	    ppion > pPidOmega.GetProbability(AliPID::kProton)   )    lIsPosInOmegaPion = kTRUE;
		
      }// end if V0 positive track with existing combined PID	
	
	
	// Combined VO-negative-daughter PID
      AliPID nPidXi;		nPidXi.SetPriors(    lPriorsGuessXi    );
      AliPID nPidOmega;	nPidOmega.SetPriors( lPriorsGuessOmega );
		
      if ( nTrackXi->IsOn(AliESDtrack::kESDpid) ) {  // Combined PID exists
        Double_t r[10] = {0.}; nTrackXi->GetESDpid(r);
        nPidXi.SetProbabilities(r);
        nPidOmega.SetProbabilities(r);
		
		// Check if the V0 negative track is a pi- (case for Xi-)
        Double_t ppion = nPidXi.GetProbability(AliPID::kPion);
        if (ppion > nPidXi.GetProbability(AliPID::kElectron) &&
            ppion > nPidXi.GetProbability(AliPID::kMuon)     &&
            ppion > nPidXi.GetProbability(AliPID::kKaon)     &&
            ppion > nPidXi.GetProbability(AliPID::kProton)   )     lIsNegInXiPion = kTRUE;

		// Check if the V0 negative track is an anti-proton (case for Xi+)
        Double_t pproton = nPidXi.GetProbability(AliPID::kProton);
	if (pproton > nPidXi.GetProbability(AliPID::kElectron) &&
	    pproton > nPidXi.GetProbability(AliPID::kMuon)     &&
	    pproton > nPidXi.GetProbability(AliPID::kPion)     &&
	    pproton > nPidXi.GetProbability(AliPID::kKaon)     )     lIsNegInXiProton = kTRUE;
		
		// Check if the V0 negative track is a pi- (case for Omega-)
	ppion = 0.;
	ppion = nPidOmega.GetProbability(AliPID::kPion);
	if (ppion > nPidOmega.GetProbability(AliPID::kElectron) &&
	    ppion > nPidOmega.GetProbability(AliPID::kMuon)     &&
	    ppion > nPidOmega.GetProbability(AliPID::kKaon)     &&
	    ppion > nPidOmega.GetProbability(AliPID::kProton)   )    lIsNegInOmegaPion = kTRUE;
		
		// Check if the V0 negative track is an anti-proton (case for Omega+)
	pproton = 0.;
	pproton = nPidOmega.GetProbability(AliPID::kProton);
	if (pproton > nPidOmega.GetProbability(AliPID::kElectron) &&
	    pproton > nPidOmega.GetProbability(AliPID::kMuon)     &&
	    pproton > nPidOmega.GetProbability(AliPID::kPion)     &&
	    pproton > nPidOmega.GetProbability(AliPID::kKaon)     )  lIsNegInOmegaProton = kTRUE;
      }// end if V0 negative track with existing combined PID	
	
		
	// Combined bachelor PID
      AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
      AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
      if ( bachTrackXi->IsOn(AliESDtrack::kESDpid) ) {  // Combined PID exists
        Double_t r[10] = {0.}; bachTrackXi->GetESDpid(r);
        bachPidXi.SetProbabilities(r);
        bachPidOmega.SetProbabilities(r);
		// Check if the bachelor track is a pion
        Double_t ppion = bachPidXi.GetProbability(AliPID::kPion);
        if (ppion > bachPidXi.GetProbability(AliPID::kElectron) &&
          ppion > bachPidXi.GetProbability(AliPID::kMuon)     &&
          ppion > bachPidXi.GetProbability(AliPID::kKaon)     &&
          ppion > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
		// Check if the bachelor track is a kaon
        Double_t pkaon = bachPidOmega.GetProbability(AliPID::kKaon);
	if (pkaon > bachPidOmega.GetProbability(AliPID::kElectron) &&
	    pkaon > bachPidOmega.GetProbability(AliPID::kMuon)     &&
	    pkaon > bachPidOmega.GetProbability(AliPID::kPion)     &&
	    pkaon > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;	
      }// end if bachelor track with existing combined PID
	
	
	// B - TPC PID : 3-sigma bands on Bethe-Bloch curve
	
        // Bachelor
      if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 4) lIsBachelorKaonForTPC = kTRUE;
      if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 4) lIsBachelorPionForTPC = kTRUE;
      
        // Negative V0 daughter
      if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 4) lIsNegPionForTPC   = kTRUE;
      if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < 4) lIsNegProtonForTPC = kTRUE;
       
        // Positive V0 daughter
      if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 4) lIsPosPionForTPC   = kTRUE;
      if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 4) lIsPosProtonForTPC = kTRUE;
        
      /*
      const AliExternalTrackParam *pInnerWallTrackXi    = pTrackXi    ->GetInnerParam();
      const AliExternalTrackParam *nInnerWallTrackXi    = nTrackXi    ->GetInnerParam();
      const AliExternalTrackParam *bachInnerWallTrackXi = bachTrackXi ->GetInnerParam();
      if (pInnerWallTrackXi && nInnerWallTrackXi && bachInnerWallTrackXi ) {
                
        Double_t pMomInnerWall    = pInnerWallTrackXi   ->GetP();
        Double_t nMomInnerWall    = nInnerWallTrackXi   ->GetP();
        Double_t bachMomInnerWall = bachInnerWallTrackXi->GetP();
              
        // Bachelor
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < 3)                              lIsBachelorPionForTPC = kTRUE;
        if (bachMomInnerWall < 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 5) lIsBachelorKaonForTPC = kTRUE;
        if (bachMomInnerWall > 0.350  && TMath::Abs(fESDpid->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < 3) lIsBachelorKaonForTPC = kTRUE;
                
        // Negative V0 daughter
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < 3  )                           lIsNegPionForTPC   = kTRUE;
        if (nMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 5 )   lIsNegProtonForTPC = kTRUE;
        if (nMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( nTrackXi,AliPID::kProton ) ) < 3 )   lIsNegProtonForTPC = kTRUE;
                
        // Positive V0 daughter
        if (TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < 3 )                            lIsPosPionForTPC   = kTRUE;
        if (pMomInnerWall < 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 5)     lIsPosProtonForTPC = kTRUE;
        if (pMomInnerWall > 0.6  && TMath::Abs(fESDpid->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < 3)     lIsPosProtonForTPC = kTRUE;
      }
      */
		
		// - II.Step 6 : extra info for QA (ESD)
		// miscellaneous pieces of info that may help regarding data quality assessment.
		//-------------

      xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
      lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
      lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
		
      xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
      lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
      lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

      lChargeXi = xi->Charge();

      lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
	
      lRapXi    = xi->RapXi();
      lRapOmega = xi->RapOmega();
      lEta      = xi->Eta();
      lTheta    = xi->Theta() *180.0/TMath::Pi();
      lPhi      = xi->Phi()   *180.0/TMath::Pi();
      lAlphaXi  = xi->AlphaXi();
      lPtArmXi  = xi->PtArmXi();
	
	
	//FIXME : Extra-cut = Anti-splitting cut for lambda daughters
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
	
      //FIXME : Just to know which file is currently open : locate the file containing Xi
      // cout << "Name of the file containing Xi candidate(s) :" 
      //        << CurrentFileName() 
      //        << " / entry: "     << Entry()
      //        << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
      //        << " : mass(Xi) = " << xi->GetEffMassXi() 
      //        << " / charge = "   << lChargeXi
      //        << " / pt(Casc) = " << lXiTransvMom
      //        << " / Decay 2d R(Xi) = " << lXiRadius 
      //        << " / Track Index(Pos)  = " << lIdxPosXi << "/ Nb(TPC clusters) = " << lPosTPCClusters 
      //        << " / Track Index(Neg)  = " << lIdxNegXi << "/ Nb(TPC clusters) = " << lNegTPCClusters 
      //        << " / Track Index(Bach) = " << lBachIdx  << "/ Nb(TPC clusters) = " << lBachTPCClusters 
      //        << endl;

	
		// II.Step 7 - Complementary info for monitoring the cascade cut variables
	
      const AliMultiplicity *lAliMult = lESDevent->GetMultiplicity();
      if(fAnalysisType == "ESD") lSPDTrackletsMultiplicity       = lAliMult->GetNumberOfTracklets();
      else if(fAnalysisType == "AOD") lSPDTrackletsMultiplicity = lAODevent->GetTracklets()->GetNumberOfTracklets();
	
        	// II.Step 8 - Azimuthal correlation study
		//-------------
	
      lTVect3MomXi.SetXYZ( lXiMomX, lXiMomY, lXiMomZ );
      lArrTrackID[0] = pTrackXi   ->GetID();
      lArrTrackID[1] = nTrackXi   ->GetID();
      lArrTrackID[2] = bachTrackXi->GetID();
	
	
    }// end of ESD treatment
  
 
    if (fAnalysisType == "AOD") {
	
	// -------------------------------------
	// II.AOD - Calcultaion Part dedicated to Xi vertices (ESD)
	
      const AliAODcascade *xi = lAODevent->GetCascade(iXi);
      if (!xi) continue;
		
	// Just to know which file is currently open : locate the file containing Xi
	// cout << "Name of the file containing Xi candidate(s) :" <<  fesdH->GetTree()->GetCurrentFile()->GetName() << endl;
	
	
		// - II.Step 1 : Characteristics of the event : prim. Vtx + magnetic field (AOD)
		//-------------
	

      lTrkgPrimaryVtxRadius3D = -500. ;
      // FIXME : We don't have the different prim. vertex at the AOD level -> nothing to do.

      lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
                                              lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
					      lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );
		
	
		// - II.Step 2 : Assigning the necessary variables for specific AliAODcascade data members (AOD)	
		//-------------
	
      lEffMassXi  			= xi->MassXi(); // default working hypothesis : cascade = Xi- decay
      lChi2Xi 			= xi->Chi2Xi();
      lDcaXiDaughters			= xi->DcaXiDaughters();
      lXiCosineOfPointingAngle 	= xi->CosPointingAngleXi( lBestPrimaryVtxPos[0], 
								  lBestPrimaryVtxPos[1], 
								  lBestPrimaryVtxPos[2] );
					// Take care : 
					// the best available vertex should be used (like in AliCascadeVertexer)

      lPosXi[0] = xi->DecayVertexXiX();
      lPosXi[1] = xi->DecayVertexXiY();
      lPosXi[2] = xi->DecayVertexXiZ();
      lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
		

		// - II.Step 3 : around the tracks : Bach + V0 (AOD)
		// ~ Necessary variables for AODcascade data members coming from the AODv0 part (inheritance)
		//-------------
		
      lChargeXi 			= xi->ChargeXi();
	
      if ( lChargeXi < 0)  	
	lInvMassLambdaAsCascDghter	= xi->MassLambda();
      else 			
        lInvMassLambdaAsCascDghter	= xi->MassAntiLambda();

      lDcaV0DaughtersXi 		= xi->DcaV0Daughters(); 
      lV0Chi2Xi 			= xi->Chi2V0();

      lV0CosineOfPointingAngleXi 	= xi->CosPointingAngle( lBestPrimaryVtxPos );
      lDcaV0ToPrimVertexXi 		= xi->DcaV0ToPrimVertex();
	
      lDcaBachToPrimVertexXi 		= xi->DcaBachToPrimVertex(); 
	
	
      lPosV0Xi[0] = xi->DecayVertexV0X();
      lPosV0Xi[1] = xi->DecayVertexV0Y();
      lPosV0Xi[2] = xi->DecayVertexV0Z(); 
      lV0RadiusXi	= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );

      lDcaPosToPrimVertexXi		= xi->DcaPosToPrimVertex(); 
      lDcaNegToPrimVertexXi		= xi->DcaNegToPrimVertex(); 

		// - II.Step 4 : around effective masses (AOD)
		// ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		//-------------

      if ( lChargeXi < 0 )		lInvMassXiMinus 	= xi->MassXi();
      if ( lChargeXi > 0 )		lInvMassXiPlus 		= xi->MassXi();
      if ( lChargeXi < 0 )		lInvMassOmegaMinus 	= xi->MassOmega();
      if ( lChargeXi > 0 )		lInvMassOmegaPlus 	= xi->MassOmega();

	
		// - II.Step 5 : PID on the daughters (To be developed ...)
		//-------------
	
	// Combined PID
	
	/* 
	// Reasonable guess for the priors for the cascade track sample
	Double_t lPriorsGuessXi[5]    = {0.0, 0.0, 2, 0, 1};
	Double_t lPriorsGuessOmega[5] = {0.0, 0.0, 1, 1, 1};
	AliPID bachPidXi;	bachPidXi.SetPriors(    lPriorsGuessXi    );
	AliPID bachPidOmega;	bachPidOmega.SetPriors( lPriorsGuessOmega );
	
	const AliAODTrack *bachTrackXi = lAODevent->GetTrack( xi->GetBachID() ); // FIXME : GetBachID not implemented ?
	
     if( bachTrackXi->IsOn(AliESDtrack::kESDpid) ){  // Combined PID exists, the AOD flags = a copy of the ESD ones
		Double_t r[10]; bachTrackXi->GetPID(r);
		bachPidXi.SetProbabilities(r);
		bachPidOmega.SetProbabilities(r);
		// Check if the bachelor track is a pion
		Double_t ppion = bachPidXi.GetProbability(AliPID::kPion);
		if (ppion > bachPidXi.GetProbability(AliPID::kElectron) &&
		    ppion > bachPidXi.GetProbability(AliPID::kMuon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kKaon)     &&
		    ppion > bachPidXi.GetProbability(AliPID::kProton)   )     lIsBachelorPion = kTRUE;
		// Check if the bachelor track is a kaon
		Double_t pkaon = bachPidOmega.GetProbability(AliPID::kKaon);
		if (pkaon > bachPidOmega.GetProbability(AliPID::kElectron) &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kMuon)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kPion)     &&
		    pkaon > bachPidOmega.GetProbability(AliPID::kProton)   )  lIsBachelorKaon = kTRUE;
		
      }// end if bachelor track with existing combined PID
      */
	
	// TPC PID
	
		// - II.Step 6 : extra info for QA (AOD)
		// miscellaneous pieces onf info that may help regarding data quality assessment.
		//-------------
      lXiMomX = xi->MomXiX();
      lXiMomY = xi->MomXiY();
      lXiMomZ = xi->MomXiZ();
      lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
      lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
	
      lBachMomX = xi->MomBachX();
      lBachMomY = xi->MomBachY();
      lBachMomZ = xi->MomBachZ();		
      lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
      lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	
      lV0toXiCosineOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );

      lRapXi    = xi->RapXi();
      lRapOmega = xi->RapOmega();
      lEta      = xi->Eta();				// Will not work ! need a method Pz(), Py() Px() 
      lTheta    = xi->Theta() *180.0/TMath::Pi();     // in AODcascade.
      lPhi      = xi->Phi()   *180.0/TMath::Pi();     // Here, we will get eta, theta, phi for the V0 ...
      lAlphaXi  = xi->AlphaXi();
      lPtArmXi  = xi->PtArmXi();

		// II.Step 7 - Complementary info for monitoring the cascade cut variables
	//FIXME : missing for AOD : Tacklet Multiplicity + TPCCluster
	
		// II.Step 8 - Azimuthal correlation study
		//-------------
	
      lTVect3MomXi.SetXYZ( lXiMomX, lXiMomY, lXiMomZ );
	
      AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
      AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
      AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );	
      if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
	AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
	continue;
      }
		
      lArrTrackID[0] = pTrackXi   ->GetID();
      lArrTrackID[1] = nTrackXi   ->GetID();
      lArrTrackID[2] = bachTrackXi->GetID();
	
    }// end of AOD treatment


    // -------------------------------------
    // II.Fill - Filling the TH1,2,3Fs, HnSparses, CFContainers, FOR events with CASCADES !
        
        // if( lIsBachelorKaonForTPC )
        //                 // FIXME : Just to know which file is currently open : locate the file containing Xi
        // cout << "Name of the file containing Omega candidate(s) :" 
        //         << CurrentFileName() 
        //         << " / entry: "     << Entry()
        //         << " / in file: "   << lESDevent->GetEventNumberInFile()   // <- Cvetan / From Mihaela: AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetReadEntry();
        //         << " : mass(Omega+) = " << lInvMassOmegaPlus 
        //         << " : mass(Omega-) = " << lInvMassOmegaMinus
        //         << " / charge = "   << lChargeXi
        //         << " / pt(Casc) = " << lXiTransvMom
        //         << " / Decay 2d R(Xi) = " << lXiRadius 
        //         << endl;
  

	// - II.Fill.Step 1	 : primary vertex
  
    fHistTPCrefitTrackMultiplicityForCascadeEvt->Fill( nTrackWithTPCrefitMultiplicity );
        
    fHistPosV0TPCClusters           ->Fill( lPosTPCClusters  );
    fHistNegV0TPCClusters           ->Fill( lNegTPCClusters  );
    fHistBachTPCClusters            ->Fill( lBachTPCClusters );
       
    f2dHistTPCdEdxOfCascDghters     ->Fill( lInnerWallMomCascDghters[0] , lTPCSignalCascDghters[0]  );
    f2dHistTPCdEdxOfCascDghters     ->Fill( lInnerWallMomCascDghters[1] , lTPCSignalCascDghters[1]  );
    f2dHistTPCdEdxOfCascDghters     ->Fill( lInnerWallMomCascDghters[2] , lTPCSignalCascDghters[2]  );
        
    fHistVtxStatus                  ->Fill( lStatusTrackingPrimVtx   );  // 1 if tracking vtx = ok

    if ( lStatusTrackingPrimVtx ) {
      fHistPosTrkgPrimaryVtxXForCascadeEvt  ->Fill( lTrkgPrimaryVtxPos[0]    );
      fHistPosTrkgPrimaryVtxYForCascadeEvt  ->Fill( lTrkgPrimaryVtxPos[1]    );
      fHistPosTrkgPrimaryVtxZForCascadeEvt  ->Fill( lTrkgPrimaryVtxPos[2]    );
      fHistTrkgPrimaryVtxRadius             ->Fill( lTrkgPrimaryVtxRadius3D );
    }

    fHistPosBestPrimaryVtxXForCascadeEvt   ->Fill( lBestPrimaryVtxPos[0]    );
    fHistPosBestPrimaryVtxYForCascadeEvt   ->Fill( lBestPrimaryVtxPos[1]    );
    fHistPosBestPrimaryVtxZForCascadeEvt   ->Fill( lBestPrimaryVtxPos[2]    );
    fHistBestPrimaryVtxRadius              ->Fill( lBestPrimaryVtxRadius3D  );
	
    f2dHistTrkgPrimVtxVsBestPrimVtx->Fill( lTrkgPrimaryVtxRadius3D, lBestPrimaryVtxRadius3D );

        // **************************** With PID on ? ... for the signal region ? ************FIXME**************************************
    if ( ( (lChargeXi<0) && lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) ||
         ( (lChargeXi>0) && lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )  )
                // NOTE : 
                // with this condition, it could happen that a cascade candidate satisfies the wrong requirement,
                // e.g. one looks at a Xi- candidate for which lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC = kFALSE
                //      Expectation: it should be excluded.
                //      but lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC = kTRUE
                //      then this bad Xi-candidate will contribute anyway (OR condition).
                // Hence : the extra condition on the sign of the Cascade
                                                                                                 {
                //         if( TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.010 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.010)
                        
                // II.Fill.Step 2
      fHistEffMassXi			->Fill( lEffMassXi               );
      fHistChi2Xi			->Fill( lChi2Xi                  );	// Flag CascadeVtxer: Cut Variable a
      fHistDcaXiDaughters		->Fill( lDcaXiDaughters          );	// Flag CascadeVtxer: Cut Variable e 
      fHistDcaBachToPrimVertex  	->Fill( lDcaBachToPrimVertexXi   );	// Flag CascadeVtxer: Cut Variable d
      fHistXiCosineOfPointingAngle	->Fill( lXiCosineOfPointingAngle );	// Flag CascadeVtxer: Cut Variable f
      fHistXiRadius			->Fill( lXiRadius                );	// Flag CascadeVtxer: Cut Variable g+h
                
                
                // II.Fill.Step 3
      fHistMassLambdaAsCascDghter	->Fill( lInvMassLambdaAsCascDghter );	// Flag CascadeVtxer: Cut Variable c
      fHistV0Chi2Xi			->Fill( lV0Chi2Xi                  );	
      fHistDcaV0DaughtersXi		->Fill( lDcaV0DaughtersXi          );
      fHistV0CosineOfPointingAngleXi	->Fill( lV0CosineOfPointingAngleXi ); 
      fHistV0RadiusXi			->Fill( lV0RadiusXi                );
                
      fHistDcaV0ToPrimVertexXi	->Fill( lDcaV0ToPrimVertexXi       );	// Flag CascadeVtxer: Cut Variable b
      fHistDcaPosToPrimVertexXi	->Fill( lDcaPosToPrimVertexXi      );
      fHistDcaNegToPrimVertexXi	->Fill( lDcaNegToPrimVertexXi      );
                
        
                // II.Fill.Step 4 : extra QA info
                
      fHistChargeXi                   ->Fill( lChargeXi      );
      fHistV0toXiCosineOfPointingAngle->Fill( lV0toXiCosineOfPointingAngle );
        
      if ( TMath::Abs( lInvMassXiMinus-1.3217 ) < 0.012 || TMath::Abs( lInvMassXiPlus-1.3217 ) < 0.012) {// One InvMass should be different from 0
        fHistXiTransvMom        ->Fill( lXiTransvMom   );
        fHistXiTotMom           ->Fill( lXiTotMom      );
        
        fHistBachTransvMomXi    ->Fill( lBachTransvMom );
        fHistBachTotMomXi       ->Fill( lBachTotMom    );
        
        fHistRapXi              ->Fill( lRapXi         );
        fHistEtaXi              ->Fill( lEta           );
        fHistThetaXi            ->Fill( lTheta         );
        fHistPhiXi              ->Fill( lPhi           );
      }

      if ( TMath::Abs( lInvMassOmegaMinus-1.672 ) < 0.012 || TMath::Abs( lInvMassOmegaPlus-1.672 ) < 0.012 ) {// One InvMass should be different from 0
        fHistRapOmega           ->Fill( lRapOmega            ); 
      }
        
      f2dHistArmenteros	        ->Fill( lAlphaXi, lPtArmXi   );
    }// end with PID ...
	
	// II.Fill.Step 5 : inv mass plots 1D
    if ( lChargeXi < 0 ) {
      fHistMassXiMinus	       ->Fill( lInvMassXiMinus    );
          // MN Extra cut from Iouri, not sure that this is the right way to implement it
  //    if (TMath::Abs(lInvMassXiMinus-1.3217)>0.008)
	fHistMassOmegaMinus	       ->Fill( lInvMassOmegaMinus );
      if(lIsBachelorPion)	fHistMassWithCombPIDXiMinus    ->Fill( lInvMassXiMinus    );
      if(lIsBachelorKaon)	fHistMassWithCombPIDOmegaMinus ->Fill( lInvMassOmegaMinus );

      fHistDcaXiDaughtersvsInvMass->Fill(lDcaXiDaughters,lInvMassXiMinus);
      fHistDcaBachToPrimVertexvsInvMass->Fill(lDcaBachToPrimVertexXi,lInvMassXiMinus); 
      fHistXiCosineOfPointingAnglevsInvMass->Fill(lXiCosineOfPointingAngle,lInvMassXiMinus);
      fHistMassLambdaAsCascDghtervsInvMass->Fill(lInvMassLambdaAsCascDghter,lInvMassXiMinus);
      fHistDcaV0DaughtersXivsInvMass->Fill(lDcaV0DaughtersXi,lInvMassXiMinus);
      fHistDcaV0ToPrimVertexXivsInvMass->Fill(lDcaV0ToPrimVertexXi,lInvMassXiMinus);
    }
	
    if ( lChargeXi > 0 ) {
      fHistMassXiPlus		       ->Fill( lInvMassXiPlus     );
                 // MN Extra cut from Iouri, not sure that this is the right way to implement it 
    //  if (TMath::Abs(lInvMassXiPlus-1.3217)>0.008)
	fHistMassOmegaPlus	       ->Fill( lInvMassOmegaPlus  );
      if(lIsBachelorPion)	fHistMassWithCombPIDXiPlus     ->Fill( lInvMassXiPlus     );
      if(lIsBachelorKaon)	fHistMassWithCombPIDOmegaPlus  ->Fill( lInvMassOmegaPlus  );
    }
	
	
	// II.Fill.Step 6 : inv mass plots 2D, 3D
    if ( lChargeXi < 0 ) {
      f2dHistEffMassLambdaVsEffMassXiMinus->Fill( lInvMassLambdaAsCascDghter, lInvMassXiMinus ); 
      f2dHistEffMassXiVsEffMassOmegaMinus ->Fill( lInvMassXiMinus, lInvMassOmegaMinus );
      f2dHistXiRadiusVsEffMassXiMinus     ->Fill( lXiRadius, lInvMassXiMinus );
      f2dHistXiRadiusVsEffMassOmegaMinus  ->Fill( lXiRadius, lInvMassOmegaMinus );
//	f3dHistXiPtVsEffMassVsYXiMinus      ->Fill( lXiTransvMom, lInvMassXiMinus,    lRapXi    );
//	f3dHistXiPtVsEffMassVsYOmegaMinus   ->Fill( lXiTransvMom, lInvMassOmegaMinus, lRapOmega );
    } else {
      f2dHistEffMassLambdaVsEffMassXiPlus ->Fill( lInvMassLambdaAsCascDghter, lInvMassXiPlus );
      f2dHistEffMassXiVsEffMassOmegaPlus  ->Fill( lInvMassXiPlus, lInvMassOmegaPlus );
      f2dHistXiRadiusVsEffMassXiPlus      ->Fill( lXiRadius, lInvMassXiPlus);
      f2dHistXiRadiusVsEffMassOmegaPlus   ->Fill( lXiRadius, lInvMassOmegaPlus );
//	f3dHistXiPtVsEffMassVsYXiPlus       ->Fill( lXiTransvMom, lInvMassXiPlus,    lRapXi    );
//	f3dHistXiPtVsEffMassVsYOmegaPlus    ->Fill( lXiTransvMom, lInvMassOmegaPlus, lRapOmega );
    }
	
	// - Filling the AliCFContainers related to PID
	
    Double_t lContainerPIDVars[4] = {0.0};
	
	// Xi Minus		
    if ( lChargeXi < 0 ) {
      lContainerPIDVars[0] = lXiTransvMom       ;
      lContainerPIDVars[1] = lInvMassXiMinus    ;
      lContainerPIDVars[2] = lRapXi             ;
      lContainerPIDVars[3] = lSPDTrackletsMultiplicity;       
			
		// No PID
      fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 0); // No PID
		// TPC PID
      if ( lIsBachelorPionForTPC  )
        fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
      if ( lIsBachelorPionForTPC && 
	   lIsPosProtonForTPC     )
	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
      if ( lIsBachelorPionForTPC && 
	   lIsPosProtonForTPC    && 
	   lIsNegPionForTPC       )
	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
      if ( lIsBachelorPion        )
	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
      if ( lIsBachelorPion       && 
	   lIsPosInXiProton    )
	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
      if (lIsBachelorPion     && 
	  lIsPosInXiProton && 
	  lIsNegInXiPion    )
	fCFContCascadePIDXiMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
	
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
    // Xi Plus		
    if ( lChargeXi > 0 ) {
      lContainerPIDVars[0] = lXiTransvMom       ;
      lContainerPIDVars[1] = lInvMassXiPlus     ;
      lContainerPIDVars[2] = lRapXi             ;
      lContainerPIDVars[3] = lSPDTrackletsMultiplicity;       
			
		// No PID
      fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 0); // No PID
        	// TPC PID
      if ( lIsBachelorPionForTPC  )
	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
      if ( lIsBachelorPionForTPC && 
	   lIsNegProtonForTPC     )
	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
      if ( lIsBachelorPionForTPC && 
	   lIsNegProtonForTPC    && 
	   lIsPosPionForTPC       )
	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
      if ( lIsBachelorPion        )
	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
      if ( lIsBachelorPion       && 
          lIsNegInXiProton    )
	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
      if (lIsBachelorPion     && 
	  lIsNegInXiProton && 
	  lIsPosInXiPion    )
	fCFContCascadePIDXiPlus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
	
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Minus		
    if ( lChargeXi < 0 ) {
      lContainerPIDVars[0] = lXiTransvMom       ;
      lContainerPIDVars[1] = lInvMassOmegaMinus ;
      lContainerPIDVars[2] = lRapOmega          ;
      lContainerPIDVars[3] = lSPDTrackletsMultiplicity;      
			
		// No PID
      fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 0); // No PID
        	// TPC PID
      if ( lIsBachelorKaonForTPC  )
	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 1); // TPC PID / 4-#sigma cut on Bachelor track
		
      if ( lIsBachelorKaonForTPC && 
           lIsPosProtonForTPC     )
	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 2); // TPC PID / 4-#sigma cut on Bachelor+Baryon tracks
		
      if ( lIsBachelorKaonForTPC && 
	   lIsPosProtonForTPC    && 
	   lIsNegPionForTPC       )
	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 3); // TPC PID / 4-#sigma cut on Bachelor+Baryon+Meson tracks
		
		// Combined PID
      if ( lIsBachelorKaon        )
	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 4); // Comb. PID / Bachelor
		
      if ( lIsBachelorKaon       && 
	   lIsPosInOmegaProton    )
	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 5); // Comb. PID / Bachelor+Baryon
		
      if (lIsBachelorKaon     && 
	  lIsPosInOmegaProton && 
	  lIsNegInOmegaPion    )
	fCFContCascadePIDOmegaMinus->Fill(lContainerPIDVars, 6); // Comb. PID / Bachelor+Baryon+Meson
    }
	
    lContainerPIDVars[0] = 0.; lContainerPIDVars[1] = 0.; lContainerPIDVars[2] = 0.; lContainerPIDVars[3] = 0.;
	
	// Omega Plus		
    if ( lChargeXi > 0 ) {
      lContainerPIDVars[0] = lXiTransvMom       ;
      lContainerPIDVars[1] = lInvMassOmegaPlus ;
      lContainerPIDVars[2] = lRapOmega          ;
      lContainerPIDVars[3] = lSPDTrackletsMultiplicity; 
			
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
	
        	
        // II.Fill.Step 7 : filling the AliCFContainer (optimisation of topological selections)
     if (fUseCFContCascadeCuts) { 
       Double_t lContainerCutVars[20] = {0.0};
                        
       lContainerCutVars[0]  = lDcaXiDaughters;
       lContainerCutVars[1]  = lDcaBachToPrimVertexXi;
       lContainerCutVars[2]  = lXiCosineOfPointingAngle;
       lContainerCutVars[3]  = lXiRadius;
       lContainerCutVars[4]  = lInvMassLambdaAsCascDghter;
       lContainerCutVars[5]  = lDcaV0DaughtersXi;
       lContainerCutVars[6]  = lV0CosineOfPointingAngleXi;
       lContainerCutVars[7]  = lV0RadiusXi;
       lContainerCutVars[8]  = lDcaV0ToPrimVertexXi;	
       lContainerCutVars[9]  = lDcaPosToPrimVertexXi;
       lContainerCutVars[10] = lDcaNegToPrimVertexXi;
       
       lContainerCutVars[13] = lXiTransvMom;
        
       lContainerCutVars[16] = lBestPrimaryVtxPos[2];
       lContainerCutVars[17] = nTrackWithTPCrefitMultiplicity;  // FIXME : nTrackWithTPCrefitMultiplicity not checked for AOD ... = 0
       lContainerCutVars[18] = lSPDTrackletsMultiplicity;       
       lContainerCutVars[19] = lBachTPCClusters;                // FIXME : BachTPCClusters          is not available for AOD ... = -1
        
       if ( lChargeXi < 0 ) {
         lContainerCutVars[11] = lInvMassXiMinus;
         lContainerCutVars[12] = 1.63;
         lContainerCutVars[14] = lRapXi;
         lContainerCutVars[15] = -1.;
//         if ( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )
       fCFContCascadeCuts->Fill(lContainerCutVars,0); // for Xi-

         lContainerCutVars[11] = 1.26;
         lContainerCutVars[12] = lInvMassOmegaMinus;
         lContainerCutVars[14] = -1.;
         lContainerCutVars[15] = lRapOmega;
//         if( lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )
       fCFContCascadeCuts->Fill(lContainerCutVars,2); // for Omega-
       } else {
         lContainerCutVars[11] = lInvMassXiPlus;
         lContainerCutVars[12] = 1.63;
         lContainerCutVars[14] = lRapXi;
         lContainerCutVars[15] = -1.; 
//         if ( lIsBachelorPionForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )
       fCFContCascadeCuts->Fill(lContainerCutVars,1); // for Xi+

         lContainerCutVars[11] = 1.26;
         lContainerCutVars[12] = lInvMassOmegaPlus;
         lContainerCutVars[14] = -1.;
         lContainerCutVars[15] = lRapOmega;
//         if( lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )
       fCFContCascadeCuts->Fill(lContainerCutVars,3); // for Omega+
       }
     } 
                        
        // II.Fill.Step 8 :  angular correlations
/*        
        if( lChargeXi < 0 ){
                if( lIsBachelorPionForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )      DoAngularCorrelation("Xi-",    lInvMassXiMinus,    lArrTrackID, lTVect3MomXi, lEta);
                if( lIsBachelorKaonForTPC   && lIsPosProtonForTPC    && lIsNegPionForTPC )      DoAngularCorrelation("Omega-", lInvMassOmegaMinus, lArrTrackID, lTVect3MomXi, lEta);
        }
        else{
                if( lIsBachelorPionForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )      DoAngularCorrelation("Xi+",    lInvMassXiPlus,    lArrTrackID, lTVect3MomXi, lEta);
                if( lIsBachelorKaonForTPC   && lIsNegProtonForTPC    && lIsPosPionForTPC )      DoAngularCorrelation("Omega+", lInvMassOmegaPlus, lArrTrackID, lTVect3MomXi, lEta);
        }
        
*/	
   }// end of the Cascade loop (ESD or AOD)
    
  
  // Post output data.
 PostData(1, fListHistCascade);
}

/*
void AliAnalysisTaskCheckCascadePbPb::DoAngularCorrelation( const Char_t   *lCascType, 
							      Double_t  lInvMassCascade, 
							const Int_t    *lArrTrackID,
							      TVector3 &lTVect3MomXi, 
							      Double_t  lEtaXi ){
  // Perform the Delta(Phi)Delta(Eta) analysis 
  // by properly filling the THnSparseF 
	
	TString lStrCascType( lCascType );
	
	Double_t lCascPdgMass = 0.0;
	if( lStrCascType.Contains("Xi") )       lCascPdgMass = 1.3217;
	if( lStrCascType.Contains("Omega") )    lCascPdgMass = 1.6724;
	
	if( lInvMassCascade > lCascPdgMass + 0.010) return;
	if( lInvMassCascade < lCascPdgMass - 0.010) return;
	// Check the Xi- candidate is within the proper mass window m0 +- 10 MeV
	
	
	// 1st loop: check there is no track with a higher pt ...
	// = The cascade is meant to be a leading particle : Pt(Casc) > any track in the event
	for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
	{// Loop over all the tracks of the event
	
		AliVTrack *lCurrentTrck = dynamic_cast<AliVTrack*>( (InputEvent())->GetTrack( TrckIdx ) );
			if (!lCurrentTrck ) {
				AliWarning("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
				continue;
			}
		if(lTVect3MomXi.Pt() < lCurrentTrck->Pt() ) return;	
		// Room for improvement: //FIXME
		// 1. There is a given resolution on pt : maybe release the cut Pt(casc) < Pt(track)*90% ?
		// 2. Apply this cut only when DeltaPhi(casc, track) > 90 deg = when track is in the away-side ?
		// 3. Anti-splitting cut (like in Femto analysis) ?
			
	}// end control loop
	
	// 2nd loop: filling loop
	for(Int_t TrckIdx = 0; TrckIdx < (InputEvent())->GetNumberOfTracks() ; TrckIdx++ )
	{// Loop over all the tracks of the event
	
		AliVTrack *lCurrentTrck = dynamic_cast<AliVTrack*>( (InputEvent())->GetTrack( TrckIdx ) );
			if (!lCurrentTrck ) {
				AliWarning("ERROR Correl. Study : Could not retrieve a track while looping over the event tracks ...");
				continue;
			}
				
		// Room for improvement: //FIXME
		// 1. Loop only on primary tracks ?	
		// 2. Exclude the tracks that build the condisdered cascade = the bachelor + the V0 dghters
		//     This may bias the outcome, especially for low multplicity events.
		// Note : For ESD event, track ID == track index.
			if(lCurrentTrck->GetID() == lArrTrackID[0]) continue;
			if(lCurrentTrck->GetID() == lArrTrackID[1]) continue;
			if(lCurrentTrck->GetID() == lArrTrackID[2]) continue;
			
		TVector3 lTVect3MomTrck(lCurrentTrck->Px(), lCurrentTrck->Py(), lCurrentTrck->Pz() );
		
		// 2 hypotheses made here :
		//   - The Xi trajectory is a straight line,
		//   - The Xi doesn't loose any energy by crossing the first layer(s) of ITS, if ever;
		//      So, meaning hyp: vect p(Xi) at the emission = vect p(Xi) at the decay vertex
		//      By doing this, we introduce a systematic error on the cascade Phi ...
		// Room for improvement: take into account the curvature of the Xi trajectory //FIXME
	
   		Double_t lHnSpFillVar[5] = {0.};
		lHnSpFillVar[0] = lTVect3MomXi.DeltaPhi(lTVect3MomTrck) * 180.0/TMath::Pi(); // Delta phi(Casc,Track) (deg)
		if(lHnSpFillVar[0] < -50.0) lHnSpFillVar[0] += 360.0; 		
		lHnSpFillVar[1] = lEtaXi - lCurrentTrck->Eta(); 			   // Delta eta(Casc,Track)
		lHnSpFillVar[2] = lTVect3MomXi.Pt();					   // Pt_{Casc}
		lHnSpFillVar[3] = lCurrentTrck->Pt();					   // Pt_{any track}
		lHnSpFillVar[4] = lInvMassCascade;					   // Eff. Inv Mass (control var)
		
		if(      lStrCascType.Contains("Xi-") )      fHnSpAngularCorrXiMinus    ->Fill( lHnSpFillVar );
		else if( lStrCascType.Contains("Xi+") )      fHnSpAngularCorrXiPlus     ->Fill( lHnSpFillVar );
		else if( lStrCascType.Contains("Omega-") )   fHnSpAngularCorrOmegaMinus ->Fill( lHnSpFillVar );
		else if( lStrCascType.Contains("Omega+") )   fHnSpAngularCorrOmegaPlus  ->Fill( lHnSpFillVar );
	
	}// end - Loop over all the tracks in the event

}
*/
Int_t AliAnalysisTaskCheckCascadePbPb::DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent) {
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

    }// end loop over all event tracks
    return  nTrackWithTPCrefitMultiplicity;
}





//________________________________________________________________________
void AliAnalysisTaskCheckCascadePbPb::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  TList *cRetrievedList = 0x0;
         cRetrievedList = (TList*)GetOutputData(1);
	if(!cRetrievedList){
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: ouput data container list not available\n"); return;
	}

  fHistTrackMultiplicityForCentrEvt = dynamic_cast<TH1F*> (   cRetrievedList->FindObject("fHistTrackMultiplicityForCentrEvt") );
  if (!fHistTrackMultiplicityForCentrEvt) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: fHistTrackMultiplicityForCentrEvt not available\n"); return;
	}
 
  fHistCascadeMultiplicityForCentrEvt = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistCascadeMultiplicityForCentrEvt") );
	if (!fHistCascadeMultiplicityForCentrEvt) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: fHistCascadeMultiplicityForCentrEvt not available\n"); return;
	}
	
  fHistMassXiMinus    = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiMinus") );	
	if (!fHistMassXiMinus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: fHistMassXiMinus not available\n"); return;
	}
  fHistMassXiPlus     = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassXiPlus") );
	if (!fHistMassXiPlus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: fHistMassXiPlus not available\n"); return;
	}	
  fHistMassOmegaMinus = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaMinus") );
	if (!fHistMassOmegaMinus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: fHistMassOmegaMinus not available\n"); return;
	}
  fHistMassOmegaPlus  = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistMassOmegaPlus") );	
	if (!fHistMassOmegaPlus) {
		AliWarning("ERROR - AliAnalysisTaskCheckCascadePbPb: fHistMassOmegaPlus not available\n"); return;
	}
  
  TCanvas *canCheckCascade = new TCanvas("AliAnalysisTaskCheckCascadePbPb","CheckCascade overview",10,10,1010,660);
  canCheckCascade->Divide(2,2);
  
  canCheckCascade->cd(1);
  canCheckCascade->cd(1)->SetLogy();
  fHistTrackMultiplicityForCentrEvt->SetMarkerStyle(kFullStar);  
  fHistTrackMultiplicityForCentrEvt->GetXaxis()->SetLabelFont(42);
  fHistTrackMultiplicityForCentrEvt->GetYaxis()->SetLabelFont(42);
  fHistTrackMultiplicityForCentrEvt->SetTitleFont(42, "xy");
  fHistTrackMultiplicityForCentrEvt->GetXaxis()->SetTitleOffset(1.1);
  fHistTrackMultiplicityForCentrEvt->DrawCopy("H");
  
  canCheckCascade->cd(2);  
  canCheckCascade->cd(2)->SetLogy();
  fHistCascadeMultiplicityForCentrEvt->SetMarkerStyle(kOpenSquare);
  fHistCascadeMultiplicityForCentrEvt->GetXaxis()->SetLabelFont(42);
  fHistCascadeMultiplicityForCentrEvt->GetYaxis()->SetLabelFont(42);
  fHistCascadeMultiplicityForCentrEvt->SetTitleFont(42, "xy");
  fHistCascadeMultiplicityForCentrEvt->GetXaxis()->SetTitleOffset(1.1);
  fHistCascadeMultiplicityForCentrEvt->DrawCopy("E");
  
  canCheckCascade->cd(3);  
  fHistMassXiMinus ->SetMarkerStyle(kFullCircle);
  fHistMassXiMinus ->SetMarkerSize(0.5);
  fHistMassXiMinus ->GetXaxis()->SetLabelFont(42);
  fHistMassXiMinus ->GetYaxis()->SetLabelFont(42);
  fHistMassXiMinus ->SetTitleFont(42, "xy");
  fHistMassXiMinus ->GetXaxis()->SetTitleOffset(1.1);
  fHistMassXiMinus ->GetYaxis()->SetTitleOffset(1.3);
   // fHistMassXiMinus->Rebin(2);
  fHistMassXiMinus ->GetXaxis()->SetRangeUser(1.24, 1.42);
  fHistMassXiMinus ->DrawCopy("E");
  
  fHistMassXiPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassXiPlus ->SetMarkerColor(kRed+2);
  fHistMassXiPlus ->SetLineColor(kRed+2);
  fHistMassXiPlus ->SetMarkerSize(0.5);
  // fHistMassXiPlus ->Rebin(2);
  fHistMassXiPlus ->DrawCopy("ESAME");
  
  
  TLegend *legendeXi =new TLegend(0.67,0.34,0.97,0.54);
 		legendeXi->SetTextFont(42);
 		legendeXi->SetTextSize(0.05);
 		legendeXi->SetFillColor(kWhite);
 		legendeXi->AddEntry( fHistMassXiMinus,"#Xi^{-} candidates","lp");
 		legendeXi->AddEntry( fHistMassXiPlus,"#Xi^{+} candidates","lp");
 		legendeXi->Draw();
  
  
  canCheckCascade->cd(4);  
  fHistMassOmegaPlus ->SetMarkerStyle(kOpenCircle);
  fHistMassOmegaPlus ->SetMarkerColor(kRed+2);
  fHistMassOmegaPlus ->SetLineColor(kRed+2);
  fHistMassOmegaPlus ->SetMarkerSize(0.5);
  fHistMassOmegaPlus ->GetXaxis()->SetLabelFont(42);
  fHistMassOmegaPlus ->GetYaxis()->SetLabelFont(42);
  fHistMassOmegaPlus ->SetTitleFont(42, "xy");
  fHistMassOmegaPlus ->GetXaxis()->SetTitleOffset(1.1);
  fHistMassOmegaPlus ->GetYaxis()->SetTitleOffset(1.25);
  // fHistMassOmegaPlus ->Rebin(2);
  fHistMassOmegaPlus ->GetXaxis()->SetRangeUser(1.6, 1.84);
  fHistMassOmegaPlus ->DrawCopy("E");
  
  fHistMassOmegaMinus->SetMarkerStyle(kFullCircle);
  fHistMassOmegaMinus->SetMarkerSize(0.5);
  // fHistMassOmegaMinus->Rebin(2);
  fHistMassOmegaMinus->DrawCopy("ESAME");

  
   TLegend *legendeOmega = new TLegend(0.67,0.34,0.97,0.54);
 		legendeOmega->SetTextFont(42);
 		legendeOmega->SetTextSize(0.05);
 		legendeOmega->SetFillColor(kWhite);
 		legendeOmega->AddEntry( fHistMassOmegaMinus,"#Omega^{-} candidates","lp");
 		legendeOmega->AddEntry( fHistMassOmegaPlus,"#Omega^{+} candidates","lp");
 		legendeOmega->Draw();

}
