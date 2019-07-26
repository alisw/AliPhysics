/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author:  Alexander Borissov (PNU), Dmitri Peressounko (RRC KI)         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 * updated with Alexander Borissov, December 2011                         *
 * updated for rho^0, phi meson analysis Alexander Borissov, April 2012   *
 * updated for Sigma^0 analysis Alexander Borissov, Jihye Song, Dec. 2014 *   
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------//
// Class used to do analysis on conversion + calorimeter pairs                //
// A lot of code cut-and-pasted from AliV0Reader and GammaConversion classes  //
//                                                                            //
//----------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////////

// root
#include "TChain.h"
#include "TH3.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"


/* ABB-open-keep it or ? whjy not include???
class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;
*/

// analysis
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskSigma0.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliESDCaloCluster.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"
#include "AliEMCALGeometry.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliCaloParticle.h"
#include "AliTriggerAnalysis.h"
#include "AliMultiplicity.h"
#include "AliESDcascade.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliV0ReaderV1.h"
// #include "AliAnalysisTaskGammaConvV1.h"
#include "AliVEvent.h"
#include "AliKFConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliConversionCuts.h"
#include "TParticle.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskSigma0.h"

class Riostream;
class TFile;

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSigma0)


//______________________________________________________________________________
AliAnalysisTaskSigma0::AliAnalysisTaskSigma0():
AliAnalysisTaskSE(),
    fTreeList(NULL),
     fTreeV0(NULL), 
// fTreeGammaV0(NULL),  
  fPIDResponse(0), fESDtrackCuts(0), 
// fPPVsMultUtils(0), fUtils(0),
     fkSaveV0Tree      ( kFALSE ),
//---> Variables for fTreeV0
      fTreeVariableChi2V0(0),
      fTreeVariableDcaV0Daughters(0),
      fTreeVariableDcaV0ToPrimVertex(0),
      fTreeVariableDcaPosToPrimVertex(0),
      fTreeVariableDcaNegToPrimVertex(0),
      fTreeVariableV0CosineOfPointingAngle(0),
      fTreeVariableV0Radius(0),
      fTreeVariablePt(0),
      fTreeVariableRapK0Short(0),
      fTreeVariableRapLambda(0),
      fTreeVariableInvMassK0s(0),
      fTreeVariableInvMassLambda(0),
      fTreeVariableInvMassAntiLambda(0),
      fTreeVariableAlphaV0(0),
      fTreeVariablePtArmV0(0),
      fTreeVariableNegEta(0),
      fTreeVariablePosEta(0),
      fTreeVariableNSigmasPosProton(0),
      fTreeVariableNSigmasPosPion(0),
      fTreeVariableNSigmasNegProton(0),
      fTreeVariableNSigmasNegPion(0),
      fTreeVariableDistOverTotMom(0),
      fTreeVariableLeastNbrCrossedRows(0),
      fTreeVariableLeastRatioCrossedRowsOverFindable(0),

    fV0Reader(NULL),
    fESDEvent(NULL),
    fESDpid(NULL),
//  fESDtrackCuts(NULL),
    fStack(NULL),
    fOutputContainer(NULL),
    fReaderGammas(NULL),
    tTreeEvent(NULL),
    fCentr(0.),
    fRunPeriod(-1),
    fMinOpeningAngleGhostCut(0.),
    fPHOSgeom(0x0),
    fEMCALgeom(0x0),
    fBadDistCutPHOS(3.3),
    fBadDistCutEMCAL(6.),
    fTrackEvent(NULL),
    fConvEvent(NULL),
    fPHOSEvent(NULL),
    fEMCALEvent(NULL),
    fGenpi0Event(NULL),
    fnSigmaAboveElectronLine(5.),
    fnSigmaBelowElectronLine(-3.),
    fnSigmaAbovePionLine(0.),
    fpnSigmaAbovePionLine(1.),
    fprobCut(0.),
//-------------------------------------------
// V0 cut variables
//-------------------------------------------
    fmaxR(180.),
    fmaxZ(240.),
    fetaCut(1.4),
    fptCut(0.070),    // was 0.020
    fptMaxCut(1.500),
    fchi2CutConversion(30.),   // FredBock also 30 for pp@7 TeV
//-------------------------------------------
    fZvtx(0.),
    fCentrality(0.),
    fGetCent(0),
    fMultiplicity(0.),
    fConvR(0.),
    fPhi0Mass(1.019455),
    fKaonMass(0.49677),
    fRho0Mass(0.77549),
    fPionMass(0.13957),
    fPHOSGeo(0),
    fEventCounter(0),
    fT0multA(0),
    fT0multC(0),
    fV0multA(0),
    fV0multC(0),
    fSPDmultClust(0),
    fSPDmultTracl(0),
//    fTriggerAnalysis(new AliTriggerAnalysis),
    fnCINT1B(0),
    fnCINT1A(0),
    fnCINT1C(0),
    fnCINT1E(0),
    fEtaCuts (),
//-------------------------------------------
// Tree variables
//-------------------------------------------
    fLambdaMod(0),
    fLambdaMass(0),
    fLambdaPx(0),
    fLambdaPy(0),
    fLambdaPz(0),
    fLambdaArmPt(0),
    fLambdaArmAlpha(0),
    fLambdaEta(0),
    fLambdaCosPointingAngle(0),
    fLambdaDCAtoPVPos(0),
    fLambdaDCAtoPVNeg(0),
    fLambdaDCADaughters(0),
    fLambdaRadius(0),

    fGammaMass(0),
    fGammaPx(0),
    fGammaPy(0),
    fGammaPz(0),
    fGammaCosPointingAngle(0),
    fGammaDCADaughters(0),
    fGammaDCAtoPVNeg(0),
    fGammaDCAtoPVPos(0), 
 
    fGammaEta(0),
    fGammaArmPt(0),
    fGammaArmAlpha(0),
    fGammaZConv(0),
    fGammaChi2(0),
    fGammaRadius(0),
    fGammaDCAzToPrimVtx(0),
    fSigmaMass(0),
    fSigmaPt(0),
    fSigmaArmPt(0),
    fSigmaArmAlpha(0),

    fLambdaTPx(0),
    fLambdaTPy(0),
    fLambdaTPz(0),
    fGammaTPx(0),
    fGammaTPy(0),
    fGammaTPz(0),
    fSigmaTPx(0),
    fSigmaTPy(0),
    fSigmaTPz(0),
  


    fCentralityV0M(0),
	fInputEvent(0x0),
	fMCEvent(0x0),
        fMCStack(0x0),
        fIsMC(kFALSE)  



//-------------------------------------------

{
    // Default constructor
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    
    Int_t nBin=9 ;
    for(Int_t i=0;i<nBin;i++){
        fPHOSEvents[i]=0 ;
        fEMCALEvents[i]=0;
        fConvEvents[i]=0;
        fGenpi0Events[i]=0;
    }
    for(Int_t i=0; i<6; i++){
        fPHOSBadMap[i]= 0x0 ;
    }
}


//______________________________________________________________________________
AliAnalysisTaskSigma0::AliAnalysisTaskSigma0(const char* name):
AliAnalysisTaskSE(name),
    fTreeList(NULL),
fTreeV0(NULL), 
// fTreeGammaV0(NULL),  
fPIDResponse(0), fESDtrackCuts(0), 
// fPPVsMultUtils(0), fUtils(0),
     fkSaveV0Tree      ( kFALSE ),
//---> Variables for fTreeV0
      fTreeVariableChi2V0(0),
      fTreeVariableDcaV0Daughters(0),
      fTreeVariableDcaV0ToPrimVertex(0),
      fTreeVariableDcaPosToPrimVertex(0),
      fTreeVariableDcaNegToPrimVertex(0),
      fTreeVariableV0CosineOfPointingAngle(0),
      fTreeVariableV0Radius(0),
      fTreeVariablePt(0),
      fTreeVariableRapK0Short(0),
      fTreeVariableRapLambda(0),
      fTreeVariableInvMassK0s(0),
      fTreeVariableInvMassLambda(0),
      fTreeVariableInvMassAntiLambda(0),
      fTreeVariableAlphaV0(0),
      fTreeVariablePtArmV0(0),
      fTreeVariableNegEta(0),
      fTreeVariablePosEta(0),
      fTreeVariableNSigmasPosProton(0),
      fTreeVariableNSigmasPosPion(0),
      fTreeVariableNSigmasNegProton(0),
      fTreeVariableNSigmasNegPion(0),
      fTreeVariableDistOverTotMom(0),
      fTreeVariableLeastNbrCrossedRows(0),
      fTreeVariableLeastRatioCrossedRowsOverFindable(0),

/* /---> Variables for fTreeGammaV0
      fTreeGammaVarChi2V0(0),
      fTreeGammaVarDcaV0Daughters(0),
      fTreeGammaVarDcaV0ToPrimVertex(0),
      fTreeGammaVarDcaPosToPrimVertex(0),
      fTreeGammaVarDcaNegToPrimVertex(0),
      fTreeGammaVarV0CosineOfPointingAngle(0),
      fTreeGammaVarV0Radius(0),
      fTreeGammaVarPt(0),
      fTreeGammaVarRapK0Short(0),
      fTreeGammaVarRapLambda(0),
      fTreeGammaVarInvMassK0s(0),
      fTreeGammaVarInvMassLambda(0),
      fTreeGammaVarInvMassAntiLambda(0),
      fTreeGammaVarAlphaV0(0),
      fTreeGammaVarPtArmV0(0),
      fTreeGammaVarNegEta(0),
      fTreeGammaVarPosEta(0),
      fTreeGammaVarNSigmasPosProton(0),
      fTreeGammaVarNSigmasPosPion(0),
      fTreeGammaVarNSigmasNegProton(0),
      fTreeGammaVarNSigmasNegPion(0),
      fTreeGammaVarDistOverTotMom(0),
      fTreeGammaVarLeastNbrCrossedRows(0),
      fTreeGammaVarLeastRatioCrossedRowsOverFindable(0),
*/

/*      fTreeVariableCentV0M(0),
      fTreeVariableCentV0A(0),
      fTreeVariableCentV0C(0),
      fTreeVariableCentV0MEq(0),
      fTreeVariableCentV0AEq(0),
      fTreeVariableCentV0CEq(0),
      fTreeVariableCentV0B(0),
      fTreeVariableCentV0Apartial(0),
      fTreeVariableCentV0Cpartial(0),
      fTreeVariableCentV0S(0),
      fTreeVariableCentV0SB(0),
      fTreeVariableRefMultEta8(0),
      fTreeVariableRefMultEta5(0),
      fTreeVariableRunNumber(0), */

    fV0Reader(NULL),
    fESDEvent(NULL),
    fESDpid(NULL),
//  fESDtrackCuts(NULL),
    fStack(NULL),
    fOutputContainer(NULL),
    fReaderGammas(NULL),
    tTreeEvent(NULL),
    fCentr(0.),
    fRunPeriod(-1),
    fMinOpeningAngleGhostCut(0.),
    fPHOSgeom(0x0),
    fEMCALgeom(0x0),
    fBadDistCutPHOS(3.3),
    fBadDistCutEMCAL(6.),
    fTrackEvent(NULL),
    fConvEvent(NULL),
    fPHOSEvent(NULL),
    fEMCALEvent(NULL),
    fGenpi0Event(NULL),

    fnSigmaAboveElectronLine(5.),
    fnSigmaBelowElectronLine(-3.),
    fnSigmaAbovePionLine(0.),
    fpnSigmaAbovePionLine(1.),
    fprobCut(0.),
//-------------------------------------------
// V0 cut variables
//-------------------------------------------
    fmaxR(180.),
    fmaxZ(240.),
    fetaCut(1.4),
    fptCut(0.010),       // was 0.070
    fptMaxCut(1.00),     // was 1.50
    fchi2CutConversion(30.),
//-------------------------------------------
    fZvtx(0.),
    fCentrality(0.),
    fGetCent(0),
    fMultiplicity(0.),
    fConvR(0.),
    fPhi0Mass(1.019455),
    fKaonMass(0.49677),
    fRho0Mass(0.77549),
    fPionMass(0.13957),
    fPHOSGeo(0),
    fEventCounter(0),
    fT0multA(0),
    fT0multC(0),
    fV0multA(0),
    fV0multC(0),
    fSPDmultClust(0),
    fSPDmultTracl(0),
 //   fTriggerAnalysis(new AliTriggerAnalysis),
    fnCINT1B(0),
    fnCINT1A(0),
    fnCINT1C(0),
    fnCINT1E(0),
    fEtaCuts (),

//-------------------------------------------
// Tree variables
//-------------------------------------------
    fLambdaMod(0),
    fLambdaMass(0),
    fLambdaPx(0),
    fLambdaPy(0),
    fLambdaPz(0),
    fLambdaArmPt(0),
    fLambdaArmAlpha(0),
    fLambdaEta(0),
    fLambdaCosPointingAngle(0),
    fLambdaDCAtoPVPos(0),
    fLambdaDCAtoPVNeg(0),
    fLambdaDCADaughters(0),
    fLambdaRadius(0),

    fGammaMass(0),
    fGammaPx(0),
    fGammaPy(0),
    fGammaPz(0),
    fGammaCosPointingAngle(0),
    fGammaDCADaughters(0),
    fGammaDCAtoPVNeg(0),
    fGammaDCAtoPVPos(0), 
    fGammaEta(0),
    fGammaArmPt(0),
    fGammaArmAlpha(0),
    fGammaZConv(0),
    fGammaChi2(0),
    fGammaRadius(0),
    fGammaDCAzToPrimVtx(0),
    fSigmaMass(0),
    fSigmaPt(0),
    fSigmaArmPt(0),
    fSigmaArmAlpha(0),

    fLambdaTPx(0),
    fLambdaTPy(0),
    fLambdaTPz(0),
    fGammaTPx(0),
    fGammaTPy(0),
    fGammaTPz(0),
    fSigmaTPx(0),
    fSigmaTPy(0),
    fSigmaTPz(0),
 

    fCentralityV0M(0),
    fInputEvent(0x0),
    fMCEvent(0x0),
    fMCStack(0x0),
    fIsMC(0)

{
    // Common I/O in slot 0
    DefineInput (0, TChain::Class());
    //  DefineOutput(3, TTree::Class());
    //  DefineOutput(2, TTree::Class());
    
    
    // Your private output
    DefineOutput(1, TList::Class());
    
    // Initialize the PHOS geometry
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    fEtaCuts[0] = 1.0;
    fEtaCuts[1] = 0.8;
    fEtaCuts[2] = 0.7;
    fEtaCuts[3] = 0.6;
    fEtaCuts[4] = 0.5;
    fEtaCuts[5] = 0.4;



    //fEtaCuts[0] = 1.4; fEtaCuts[1] = 1.2;fEtaCuts[2] = 1.0;fEtaCuts[3] = 0.8;fEtaCuts[4] = 0.6;fEtaCuts[5] = 0.4;
    
    
    Int_t nBin=9 ;
    for(Int_t i=0;i<nBin;i++){
        fPHOSEvents[i]=0 ;
        fEMCALEvents[i]=0;
        fConvEvents[i]=0;
        fGenpi0Events[i]=0;
    }
    for(Int_t i=0; i<6; i++){
        fPHOSBadMap[i]= 0x0 ;
    }
    for(Int_t izvtx=0;izvtx<10; izvtx++){
        fPHOSEvents[izvtx]=new TList() ;
        fEMCALEvents[izvtx] =new TList();
        fConvEvents[izvtx] =new TList() ;
        fGenpi0Events[izvtx] =new TList() ;
    }

    DefineOutput(2, TTree::Class()); // V0 Tree ???

}



//_____________________________________________________
AliAnalysisTaskSigma0::~AliAnalysisTaskSigma0()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------
    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
        AliAnalysisManager::kProofAnalysis) {
        
        if(fOutputContainer){
            fOutputContainer->Clear() ;
            delete fOutputContainer ;
        }
        
	/*       	if (fTreeGammaV0) {
	  delete fTreeGammaV0;
	  fTreeGammaV0 = 0x0;
	  } */
 


	if(fTreeList){
	  fTreeList->Clear();
	  delete fTreeList;
	} 
     
        if(tTreeEvent){
            tTreeEvent->Clear();
            delete tTreeEvent;
        }
        
	if (fTreeV0) {
	  fTreeV0->Clear();
	  delete fTreeV0;
	}
 


        
        if(fPHOSgeom){
            delete fPHOSgeom ;
            fPHOSgeom=0x0 ;
        }
        
	if(fV0Reader){
	  delete fV0Reader;
	  fV0Reader = 0x0;
	}
        
        if(fReaderGammas){
            delete fReaderGammas;
            fReaderGammas = 0x0;
        }
        
        for(Int_t ivtx=0; ivtx<10; ivtx++){
            if(fPHOSEvents[ivtx]){
                delete fPHOSEvents[ivtx] ;
                fPHOSEvents[ivtx]=0x0 ;
            }
            if(fEMCALEvents[ivtx]){
                delete fEMCALEvents[ivtx] ;
                fEMCALEvents[ivtx]=0x0 ;
            }
            if(fConvEvents[ivtx]){
                delete fConvEvents[ivtx] ;
                fConvEvents[ivtx]=0x0 ;
            }
            if(fGenpi0Events[ivtx]){
                delete fGenpi0Events[ivtx] ;
                fGenpi0Events[ivtx]=0x0 ;
            }
        }
        
    }
}
//_____________________________________________________
void AliAnalysisTaskSigma0::Init()
{
    printf(" Initialization... \n");
}


/*  /________________________________________________________________________
void AliAnalysisTaskSigma0::ConnectInputData(Option_t *)
{
    
    TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
    if (!tree) {
        Printf("ERROR: Could not read chain from input slot 0");
    } else {
        
        AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
        
        if (esdH) {
            fESD = esdH->GetEvent();
            if(fESD) {
                fESDfriend = (AliESDfriend*)fESD->FindListObject("AliESDfriend");
                if (!fESDfriend){
                    AliError("No friend found");
                }
            }
        } else {
            Printf("ERROR: Could not get ESDInputHandler");
        }
        
    }
}
*/

//____________________________________________________________
void AliAnalysisTaskSigma0::UserCreateOutputObjects()
{
  //    printf("  UserCreateOutputObjects... \n");
    if(fDebug)gDirectory->Print() ;

    if(fOutputContainer != NULL){
        delete fOutputContainer;
        fOutputContainer = NULL;
    }
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
    
    fTreeList = new TList();
    fTreeList->SetOwner(kTRUE);
    fOutputContainer->Add(fTreeList);
    
    //Adding the histograms to the output container
    fOutputContainer->Add(new TH1F("hQAEvents","Events processed",10,-0.5,9.5)) ;
    fOutputContainer->Add(new TH1F("hQAEventsTrig","Events processed",21,-0.5,20.5)) ;

//  TH1I *hEventTrigger = new TH1I("hEventTrigger","Events processed",11,0.5,10.5);
//   hEventTrigger->GetXaxis()->SetBinLabel(1,"total");
//   fOutputContainer->Add(hEventTrigger);
    
//    fOutputContainer->Add(new TH1F("R2Conv","R2 of V0",400,0.,200)) ;
//    fOutputContainer->Add(new TH1F("R3Conv","R3 of V0",400,0.,200)) ;
//    fOutputContainer->Add(new TH1F("v0sum","v0sum",200,0,20000));
    
    
 
//------------------------------------------------
// setup local variables : number of bin, minimum, maximum
//------------------------------------------------
    char key[55] ;

    Int_t npt=500 ; // Number of Pt bins
    Int_t npt2=150 ; // Double_t ptmin2 = 0. ;
    //    Double_t ptmin=0.;
    Double_t ptmax=15. ;   // was 50.
    Double_t ptmax2=15. ;

    //    Double_t mSigmax=1.4; Double_t mSigmin=1.1;
    Double_t mLammax=1.20; Double_t mLammin=1.05;
    Int_t nLambins = 150;




    
    //    Int_t nBinMass = 70; // Number of mass bins
    //    Int_t nR=50 ; // Number of radius bin
    Double_t Rmax=500. ;

    
  
   
    Double_t mSigmax2=1.4;  Double_t mSigmin2=1.1;  Int_t nBinMass2 = 280 ;
    
    //------------------------------------------------
  
    printf(" INIT V0Reader ------------------------------- \n");
    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
    if(!fV0Reader){printf("START Error: No V0 Reader");return;} // GetV0Reader
    
    //    fV0Reader = 0;
    if(fV0Reader){
      cout << "found V0Reader ---------------" << endl;
	if((AliConvEventCuts*)fV0Reader->GetEventCuts()){
	   cout << "found" << endl;
	   if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms()){
	      cout << "found event histos -----------" << endl;
	      fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
	   }  
        }    
    }  
    if(fV0Reader) {
      if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts()) {
	if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms()) {
	  cout << "found Conversion histos -------------"  << endl;
	  fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	}
      }
    }
    printf("INIT V0Reader FINISHED ------------------------------- \n");

	//	if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
	//		if (fV0Reader->GetV0FindingEfficiencyHistograms())
	//		fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

    
    fOutputContainer->Add(new TH1F("hQAMult","Multiplicity",400,-10.,40000.)) ;
    /*    fOutputContainer->Add(new TH1F("hQACentr1","Centr1",400,-10.,400.)) ;
    fOutputContainer->Add(new TH1F("hQACentr2","Centr2",400,-10.,400.)) ;
    fOutputContainer->Add(new TH1F("hQACentr3","Centr3",400,-10.,400.)) ; */


    fOutputContainer->Add(new TH1F("hConvMult","Conv Multiplicity",51,-0.5,50.5)) ;
    fOutputContainer->Add(new TH1F("hRecGammas","V0 reconstr. gammas V0 reader",51,-0.5,50.5)) ;

    fOutputContainer->Add(new TH1F("hLamMult","Lam Multiplicity",51,-0.5,50.5)) ;
    fOutputContainer->Add(new TH1F("hPHOSMult","Lam Multiplicity",51,-0.5,50.5)) ;
    
    fOutputContainer->Add(new TH1F("hNlamEv","N Lam in ev",21,-0.5,20.5)) ;
    fOutputContainer->Add(new TH1F("hNalamEv","N aLam in ev",21,-0.5,50.5)) ;
    fOutputContainer->Add(new TH1F("hNlamminalamEv","N Lam aLam in ev",41,-20.5,20.5)) ;
    
  /*
    fOutputContainer->Add(new TH1F("hT0AmC","T0 A-C", 400, -2000, 2000 )) ;
    fOutputContainer->Add(new TH1F("hT0ApC","T0 A+C", 400, -2000, 2000 )) ;
    
    fOutputContainer->Add(new TH1F("h2T0AmC","2nd T0 A-C", 400, -2000, 2000 )) ;
    fOutputContainer->Add(new TH1F("h2T0ApC","2nd T0 A+C", 400, -2000, 2000 )) ;
 */


  fOutputContainer->Add(new TH2F("hQA_PHOS_mod1_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod2_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod3_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;

    
    fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
    fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
    fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity",150,0.,150.));
    
    //    fOutputContainer->Add(new TH2F("hdEdxTrack","6dEdx Sig of accepted tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hLamdEdxTrack","7dEdx of accepted tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hAntiLamdEdxTrack","8dEdx of accepted tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hTotProtdEdxTrack","9dEdx of accepted Prot tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hTotPiondEdxTrack","9dEdx of accepted Pion tracks",100,0.,10.,150,0.,150.)) ;
    

//------------------------------------------------
// Histograms for MC
//------------------------------------------------

    fOutputContainer->Add(new TH1F("hMCgenSig0","Primary Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenSig0Eta1","Primary Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenPSig0Eta1","Primary PSig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenASig0Eta1","Primary ASig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenPSig0Eta1Rap","Primary PSig0 gen Rapidity",100, -2., 2 )) ;
    fOutputContainer->Add(new TH1F("hMCgenASig0Eta1Rap","Primary ASig0 gen Rapidity",100, -2., 2 )) ;
    fOutputContainer->Add(new TH1F("hMCRposGen","Primary3 Sig0 gen Rapidity",40, 0., 1. )) ;
    fOutputContainer->Add(new TH1F("hMCgenrecSig0","Primary Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrecSig0Eta1","Primary Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenSig0Rap","Primary Sig0 gen Rapidity",100, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hMCgenrecSig0Rap","Primary Sig0 gen Rapidity",100, -10., 10 )) ;
    //    fOutputContainer->Add(new TH2F("hMCgenSig0RapEta","Primary Sig0 gen Rapidity Eta",30,-1.5, 1.5, 30,-1.5,1.5 )) ;
    // fOutputContainer->Add(new TH2F("hMCgenSig0RapPt","Primary Sig0 gen Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
    //   fOutputContainer->Add(new TH2F("hMCgenrecSig0RapPt","Primary Sig0 genrec Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
    //  fOutputContainer->Add(new TH2F("hMCgenrec2Sig0RapPt","Primary Sig0 genrec2 Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
    fOutputContainer->Add(new TH2F("hMCgenrec3Sig0RapPt","Primary Sig0 genrec3 Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
    fOutputContainer->Add(new TH2F("hMCgenrec3Sig0mvsPt_uncorr","PP Mass vs pt Mix",nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
    fOutputContainer->Add(new TH2F("hMCgenrec3Sig0mvsPt_corr","PP Mass vs pt Mix",nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));


    fOutputContainer->Add(new TH2F("hMCgenrec3MvsPtLamgenGamgen","CC Mass vs ptnArPod",nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
    fOutputContainer->Add(new TH2F("hMCgenrec3MvsPtLamgenGamrec","CC Mass vs ptnArPod",nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
    fOutputContainer->Add(new TH2F("hMCgenrec3MvsPtLamrecGamgen","CC Mass vs ptnArPod",nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
    fOutputContainer->Add(new TH2F("hMCgenrec3MvsPtLamrecGamrec","CC Mass vs ptnArPod",nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
    
    /* END new from 6mar14 */
    fOutputContainer->Add(new TH1F("hMCgenrec2Sig0","Primary2 Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrec2Sig0Eta1","Primary2 Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrec2Sig0Rap","Primary2 Sig0 gen Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3Sig0","Primary3 Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3Sig0Eta1","Primary3 Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3PSig0Eta1","Primary3 PSig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3ASig0Eta1","Primary3 ASig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3PSig0Eta1Rap","Primary3 PSig0 gen Rapidity",100, -2., 2 )) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3ASig0Eta1Rap","Primary3 ASig0 gen Rapidity",100, -2., 2 )) ;
    fOutputContainer->Add(new TH1F("hMCgenrec3Sig0Rap","Primary3 Sig0 gen Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH2F("hMCgenrec3Sig0dMrecvsdPt","GenRec AllSig0 dm vs dpt",40, -0.04, 0.04, 40,-0.4,0.4)) ;
    fOutputContainer->Add(new TH2F("hMCgenrec3Sig0dMvsdPt","GenRec AllSig0 dm vs dpt",40, -0.04, 0.04, 40,-0.4,0.4)) ;
    fOutputContainer->Add(new TH2F("hMCgenrec3PSig0dMvsdPt","GenRec PSig0 dm vs dpt",40, -0.04, 0.04, 40,-0.4,0.4)) ;
    fOutputContainer->Add(new TH2F("hMCgenrec3ASig0dMvsdPt","GenRec ASig0 dm vs dpt",40, -0.04, 0.04, 40,-0.4,0.4)) ;
    
    //    fOutputContainer->Add(new TH2F("hMCgenrec2LamdMvsdPt","GenRec Lam dm vs dpt",40, -0.04, 0.04,     40,-1.,1.)) ;
    // fOutputContainer->Add(new TH2F("hMCgenrec2ALamdMvsdPt","GenRec ALam dm vs dpt",40, -0.04, 0.04,   40,-1.,1.)) ;
    fOutputContainer->Add(new TH1F("hMCgenLam","Primary Lambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenPLam","Primary PLambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenALam","Primary ALambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenLamEta1","Primary Lambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenPLamEta1","Primary PLambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenALamEta1","Primary ALambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenLamSig0","Primary Lambda Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenLamSig0Eta1","Primary Lambda Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenLamRap","Primary Lam gen Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hMCgenPLamRap","Primary3 PLam gen Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hMCgenALamRap","Primary3 Alam gen Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hMCgenLamSig0Rap","Primary3 Sig0 gen Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hMCgenGamSig0","Primary Gamma Sig0 gen",npt,0.,ptmax/5.)) ;
    fOutputContainer->Add(new TH1F("hMCgenGamSig0Eta1","Primary Gamma Sig0 gen eta1",npt,0.,ptmax/5.)) ;
    fOutputContainer->Add(new TH1F("hMCgenGamSig0PHOS","Primary Gamma Sig0 gen in PHOS",npt,0.,ptmax/5.)) ;
    fOutputContainer->Add(new TH1F("hMCgenGamSig0PHOSEta1","Primary Gamma Sig0 gen in PHOS eta1",npt,0.,ptmax/5.)) ;
    fOutputContainer->Add(new TH1F("EGamPHOS","E Gamma in PHOS",npt,0.,ptmax/5.)) ;


    fOutputContainer->Add(new TH1F("hRecPLam","Primary Lambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecPLamEta1","Primary Lambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecALam","Primary Lambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecALamEta1","Primary Lambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec2PLam","Primary PLambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec2ALam","Primary ALambda gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec2PLamEta1","Primary PLambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec2ALamEta1","Primary ALambda gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecLamSig0","Primary Lambda Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecLamSig0Eta1","Primary Lambda Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecLamSig0Rap","Lam form Sig0 Rapidity",200, -10., 10 )) ;
    fOutputContainer->Add(new TH1F("hRecGamSig0","Primary Lambda Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRecGamSig0Eta1","Primary Lambda Sig0 gen eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec3LamSig0","Primary Lambda Sig0 gen3",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec3LamSig0Eta1","Primary Lambda Sig0 gen3 eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec3GamSig0","Primary Lambda Sig0 gen3",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec3GamSig0Eta1","Primary Lambda Sig0 gen3 eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec4GamSig0","Primary Lambda Sig0 gen3",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec4GamSig0Eta1","Primary Lambda Sig0 gen3 eta1",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec4LamSig0","Primary Lambda Sig0 gen3",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hRec4LamSig0Eta1","Primary Lambda Sig0 gen3 eta1",npt,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH1F("hMCSigma0Phot","Primary photons",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCSigma0GammaConv","Converted photons",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCSigma0GammaV0","Converted photons with V0",npt,0.,ptmax)) ;
    
    
    fOutputContainer->Add(new TH2F("hLamMvsPt10","10DT m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hALamMvsPt10","10DT m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hLamMvsPt11","10DT m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hALamMvsPt11","10DT m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH2F("hPhi1LamMvsPt11","10DTphi1 m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hPhi1ALamMvsPt11","10DTphi1 m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH2F("hPhi2LamMvsPt11","10DTphi1 m vs pt", 100, mLammin, mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hPhi2ALamMvsPt11","10DTphi1 m vs pt",100, mLammin, mLammax, 100,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH1F("hRLam2","R Lam", 400, 0,100. )) ;
    fOutputContainer->Add(new TH1F("hRALam2","R ALam",400,0,100. )) ;
    

    //    fOutputContainer->Add(new TH2F("hLamPtPvsM11","Lam ptPvsM",   50, 0, 5, 50,0.,5. )) ;
    // fOutputContainer->Add(new TH2F("hALamPtPvsM11","ALam ptPvsM",50, 0, 5, 50,0.,5. )) ;
    
    // fOutputContainer->Add(new TH2F("hLamMvsPt12","12DT m vs pt", nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
    // fOutputContainer->Add(new TH2F("hALamMvsPt12","12DT m vs pt",nLambins, mLammin,mLammax, 100,0.,ptmax)) ;
        
    //   fOutputContainer->Add(new TH2F("hLamALam_MvsPt","Lam-ALam M vs Pt",nLambins, mLammin*2, mLammax*2.5, 100,0., ptmax); 
    //  fOutputContainer->Add(new TH2F("hLamA_MvsPt","Lam-ALam M vs Pt",nLambins, mLammin, mLammax, 100,0., ptmax)) ; 
    // fOutputContainer->Add(new TH2F("hALamA_MvsPt","Lam-ALam M vs Pt",nLambins, mLammin, mLammax, 100,0., ptmax)) ; 


    fOutputContainer->Add(new TH2F("hArPodLam","MC Armenteros-Podolanski Lam ;#alpha;p_{t}",100,-1.0,1.0, 70,0,0.7) );
    fOutputContainer->Add(new TH2F("hArPodALam","MC Armenteros-Podolanski ALam ;#alpha;p_{t}",100,-1.0,1.0, 70,0,0.7) );
    //    fOutputContainer->Add(new TH2F("hArPodGConv","MC Armenteros-Podolanski Gamma Conv ;#alpha;p_{t}",100,-1.0,1.0, 70,0,0.7) );
    // fOutputContainer->Add(new TH1F("hCGamPtConv","Conv gamma Pt2",  npt*2, ptmin, 15. )) ;
    
    fOutputContainer->Add(new TH2F("hMCArPodSig0","Armenteros-Podolanski Sig0 ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    fOutputContainer->Add(new TH2F("hMCArPodSig02","Armenteros-Podolanski Sig02  ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    // fOutputContainer->Add(new TH2F("hMCArPodSig03","Armenteros-Podolanski NoSim Sig03  ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
        
    //    fOutputContainer->Add(new TH2F("hMCccMvsArpodPt","CC Mass vs ptnArPod",nBinMass, mSigmin,mSigmax ,npt,0.,ptmax));
    
    /*    fOutputContainer->Add(new TH2F("hMCV0MvsWidth","MC w vs m 0",   100,0.,2., 1,-0.1,0.1 )) ;
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth1","MC w vs m 1",   100,0.,2., 1,-0.1,0.1 )) ;
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth2","MC w vs m 2",100,0.,2., 1,-0.1,0.1 )) ;
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth3","MC w vs m 3",100,0.,2., 1,-0.1,0.1 )) ;
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth4","MC w vs m 4",100,0.,2., 1,-0.1,0.1 )) ;
 
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth5","MC w vs m 2",100,0.,2., 1,-0.1,0.1 )) ;
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth6","MC w vs m 3",100,0.,2., 1,-0.1,0.1 )) ;
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth7","MC w vs m 4",100,0.,2., 1,-0.1,0.1 )) ;
    */


    //Single photon spectrum  //Conversion
    fOutputContainer->Add(new TH1F("hSingleConvOnFly","Single photon spectrum",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hAllMvsWidth","DT w vs m",300,0.,TMath::Pi(),400,0.,1.)) ;
    fOutputContainer->Add(new TH1F("hCheckFB","hCheckFB",15,0,15));
    fOutputContainer->Add(new TH2F("hAMbeforeFB","hAMbeforeFB",200,-1,1,1000,0,1.));
    
    // followings are important information to get efficiency
    for(Int_t is=0;   is<=5 ; is++){
        
        sprintf(key,"hMCPSig0PtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCPSig0PtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
	//        sprintf(key,"hMCPSig0PrimPtGen%d",is) ;
	//        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
	//  sprintf(key,"hMCPSig0PrimPtRec%d",is) ;
	// fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCASig0PtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCASig0PtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
	// sprintf(key,"hMCASig0PrimPtGen%d",is) ;
	// fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
	// sprintf(key,"hMCASig0PrimPtRec%d",is) ;
	// fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCPLamPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCPLamPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
	//        sprintf(key,"hMCPLamPrimPtGen%d",is) ;
	//        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
	// sprintf(key,"hMCPLamPrimPtRec%d",is) ;
	// fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCALamPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCALamPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
	//        sprintf(key,"hMCALamPrimPtGen%d",is) ;
	//        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
	// sprintf(key,"hMCALamPrimPtRec%d",is) ;
	// fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        
	/*        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCgenrec3MvsPtLamgenGamRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,key, nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
        sprintf(key,"hMCgenrec3MvsPtLamrecGamGen%d",is);
        fOutputContainer->Add(new TH2F(key,key, nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
        sprintf(key,"hMCgenrec3MvsPtLamrecGamRec%d",is);
        fOutputContainer->Add(new TH2F(key,key, nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
	*/
        
    }
        
    // end of histos in  fOutputContainer
    fOutputContainer->SetName(GetName());
    
//------------------------------------------------
// fTree Branch definitions
//------------------------------------------------
  

//-----------BASIC-INFO---------------------------
    tTreeEvent = new TTree("tTreeEvent","Events");
/* 0*/    tTreeEvent->Branch("fLambdaMod",&fLambdaMod,"fLambdaMod/F");
/* 4*/    tTreeEvent->Branch("fLambdaMass",&fLambdaMass, "fLambdaMass/F");
/* 1*/    tTreeEvent->Branch("fLambdaPx",&fLambdaPx, "fLambdaPx/F");
/* 2*/    tTreeEvent->Branch("fLambdaPy",&fLambdaPy, "fLambdaPy/F");
/* 3*/    tTreeEvent->Branch("fLambdaPz",&fLambdaPz, "fLambdaPz/F");
/*10*/    tTreeEvent->Branch("fLambdaArmPt",&fLambdaArmPt,"fLambdaArmPt/F");
/*11*/    tTreeEvent->Branch("fLambdaArmAlpha",&fLambdaArmAlpha,"fLambdaArmAlpha/F");
/*12*/    tTreeEvent->Branch("fLambdaEta",&fLambdaEta,"fLambdaEta/F");
/* 5*/    tTreeEvent->Branch("fLambdaCosPointingAngle",&fLambdaCosPointingAngle,"fLambdaCosPointingAngle/F");
/* 8*/    tTreeEvent->Branch("fLambdaDCAtoPVPos",&fLambdaDCAtoPVPos,"fLambdaDCAtoPVPos/F");
/* 7*/    tTreeEvent->Branch("fLambdaDCAtoPVNeg",&fLambdaDCAtoPVNeg,"fLambdaDCAtoPVNeg/F");
/* 6*/    tTreeEvent->Branch("fLambdaDCADaughters",&fLambdaDCADaughters,"fLambdaDCADaughters/F");
/* 9*/    tTreeEvent->Branch("fLambdaRadius",&fLambdaRadius,"fLambdaRadius/F");

/*16*/    tTreeEvent->Branch("fGammaMass",&fGammaMass, "fGammaMass/F");    
/*13*/    tTreeEvent->Branch("fGammaPx",&fGammaPx, "fGammaPx/F");
/*14*/    tTreeEvent->Branch("fGammaPy",&fGammaPy, "fGammaPy/F");
/*15*/    tTreeEvent->Branch("fGammaPz",&fGammaPz, "fGammaPz/F");
/*17*/    tTreeEvent->Branch("fGammaCosPointingAngle",&fGammaCosPointingAngle,"fGammaCosPointingAngle/F");

/* 6*/    tTreeEvent->Branch("fGammaDCADaughters",&fGammaDCADaughters,"fGammaDCADaughters/F");
/* 7*/    tTreeEvent->Branch("fGammaDCAtoPVNeg",&fGammaDCAtoPVNeg,"fGammaDCAtoPVNeg/F");
          tTreeEvent->Branch("fGammaDCAtoPVPos",&fGammaDCAtoPVPos,"fGammaDCAtoPVPos/F");
/*21*/    tTreeEvent->Branch("fGammaEta",&fGammaEta,"fGammaEta/F");

/*19*/    tTreeEvent->Branch("fGammaArmPt",&fGammaArmPt,"fGammaArmPt/F");
/*20*/    tTreeEvent->Branch("fGammaArmAlpha",&fGammaArmAlpha,"fGammaArmAlpha/F");
/*23*/    tTreeEvent->Branch("fGammaZConv",&fGammaZConv,"fGammaZConv/F");
/*22*/    tTreeEvent->Branch("fGammaChi2",&fGammaChi2,"fGammaChi2/F");
/*18*/    tTreeEvent->Branch("fGammaRadius",&fGammaRadius,"fGammaRadius/F");
          tTreeEvent->Branch("fGammaDCAzToPrimVtx",&fGammaDCAzToPrimVtx,"fGammaDCAzToPrimVtx/F");

    
/*24*/    tTreeEvent->Branch("fSigmaMass",&fSigmaMass,"fSigmaMass/F");
/*25*/    tTreeEvent->Branch("fSigmaPt",&fSigmaPt, "fSigmaPt/F");
/*26*/    tTreeEvent->Branch("fSigmaArmPt",&fSigmaArmPt,"fSigmaArmPt/F");
/*27*/    tTreeEvent->Branch("fSigmaArmAlpha",&fSigmaArmAlpha,"fSigmaArmAlpha/F");

/**/    tTreeEvent->Branch("fLambdaTPx",&fLambdaTPx, "fLambdaTPx/F");
/**/    tTreeEvent->Branch("fLambdaTPy",&fLambdaTPy, "fLambdaTPy/F");
/**/    tTreeEvent->Branch("fLambdaTPz",&fLambdaTPz, "fLambdaTPz/F");
/**/    tTreeEvent->Branch("fGammaTPx",&fGammaTPx, "fGammaTPx/F");
/**/    tTreeEvent->Branch("fGammaTPy",&fGammaTPy, "fGammaTPy/F");
/**/    tTreeEvent->Branch("fGammaTPz",&fGammaTPz, "fGammaTPz/F");
/**/    tTreeEvent->Branch("fSigmaTPx",&fSigmaTPx, "fSigmaTPx/F");
/**/    tTreeEvent->Branch("fSigmaTPy",&fSigmaTPy, "fSigmaTPy/F");
/**/    tTreeEvent->Branch("fSigmaTPz",&fSigmaTPz, "fSigmaTPz/F");
/*28*/    tTreeEvent->Branch("fCentr",&fCentr,"fCentr/F");
    
    fTreeList->Add(tTreeEvent);
    


    //>>>>>>>>>  Create Basic V0 Output Tree
    fTreeV0 = new TTree( "fTreeV0", "V0 Candidates");
    //------------------------------------------------
    // fTreeV0 Branch definitions
    //------------------------------------------------

    //-----------BASIC-INFO---------------------------
    fTreeV0->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"fTreeVariableChi2V0/F");
    fTreeV0->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
    fTreeV0->Branch("fTreeVariableDcaV0ToPrimVertex",&fTreeVariableDcaV0ToPrimVertex,"fTreeVariableDcaV0ToPrimVertex/F");
    fTreeV0->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
    fTreeV0->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
    fTreeV0->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
    fTreeV0->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
    fTreeV0->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
    fTreeV0->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
    fTreeV0->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
    fTreeV0->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
    fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
    fTreeV0->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
    fTreeV0->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
    fTreeV0->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
    fTreeV0->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
    fTreeV0->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
    fTreeV0->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
    fTreeV0->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
    fTreeV0->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
    fTreeV0->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
    fTreeV0->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
    fTreeV0->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
    fTreeV0->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
    /* skip abb-27nov15 -----------MULTIPLICITY-INFO--------------------
    fTreeV0->Branch("fTreeVariableCentV0M",&fTreeVariableCentV0M,"fTreeVariableCentV0M/F");
    fTreeV0->Branch("fTreeVariableCentV0A",&fTreeVariableCentV0A,"fTreeVariableCentV0A/F");
    fTreeV0->Branch("fTreeVariableCentV0C",&fTreeVariableCentV0C,"fTreeVariableCentV0C/F");
    fTreeV0->Branch("fTreeVariableCentV0MEq",&fTreeVariableCentV0MEq,"fTreeVariableCentV0MEq/F");
    fTreeV0->Branch("fTreeVariableCentV0AEq",&fTreeVariableCentV0AEq,"fTreeVariableCentV0AEq/F");
    fTreeV0->Branch("fTreeVariableCentV0CEq",&fTreeVariableCentV0CEq,"fTreeVariableCentV0CEq/F");
    fTreeV0->Branch("fTreeVariableCentV0B",&fTreeVariableCentV0B,"fTreeVariableCentV0B/F");
    fTreeV0->Branch("fTreeVariableCentV0Apartial",&fTreeVariableCentV0Apartial,"fTreeVariableCentV0Apartial/F");
    fTreeV0->Branch("fTreeVariableCentV0Cpartial",&fTreeVariableCentV0Cpartial,"fTreeVariableCentV0Cpartial/F");
    fTreeV0->Branch("fTreeVariableCentV0S",&fTreeVariableCentV0S,"fTreeVariableCentV0S/F");
    fTreeV0->Branch("fTreeVariableCentV0SB",&fTreeVariableCentV0SB,"fTreeVariableCentV0SB/F");
    fTreeV0->Branch("fTreeVariableRefMultEta8",&fTreeVariableRefMultEta8,"fTreeVariableRefMultEta8/I");
    fTreeV0->Branch("fTreeVariableRefMultEta5",&fTreeVariableRefMultEta5,"fTreeVariableRefMultEta5/I");
    //Don't do this if not explicitly requested, takes up too much space
    if ( fkSaveExtendedRefMultInfo )
        fTreeV0->Branch("fTreeVariableRefMultDiffEta",fTreeVariableRefMultDiffEta,"fTreeVariableRefMultDiffEta[20]/I");
    fTreeV0->Branch("fTreeVariableRunNumber",&fTreeVariableRunNumber,"fTreeVariableRunNumber/I");
    */
    fTreeList->Add(fTreeV0);
 

   //>>>>>>>>>  Create Basic V0 Output Tree
   /*  fTreeGammaV0 = new TTree( "fTreeGammaV0", "V0 Candidates for Gamma");

    //-----------BASIC-INFO---------------------------
    fTreeGammaV0->Branch("fTreeGammaVarChi2V0",&fTreeGammaVarChi2V0,"fTreeGammaVarChi2V0/F");
    fTreeGammaV0->Branch("fTreeGammaVarDcaV0Daughters",&fTreeGammaVarDcaV0Daughters,"fTreeGammaVarDcaV0Daughters/F");
    fTreeGammaV0->Branch("fTreeGammaVarDcaV0ToPrimVertex",&fTreeGammaVarDcaV0ToPrimVertex,"fTreeGammaVarDcaV0ToPrimVertex/F");
    fTreeGammaV0->Branch("fTreeGammaVarDcaPosToPrimVertex",&fTreeGammaVarDcaPosToPrimVertex,"fTreeGammaVarDcaPosToPrimVertex/F");
    fTreeGammaV0->Branch("fTreeGammaVarDcaNegToPrimVertex",&fTreeGammaVarDcaNegToPrimVertex,"fTreeGammaVarDcaNegToPrimVertex/F");
    fTreeGammaV0->Branch("fTreeGammaVarV0Radius",&fTreeGammaVarV0Radius,"fTreeGammaVarV0Radius/F");
    fTreeGammaV0->Branch("fTreeGammaVarPt",&fTreeGammaVarPt,"fTreeGammaVarPt/F");
    fTreeGammaV0->Branch("fTreeGammaVarRapK0Short",&fTreeGammaVarRapK0Short,"fTreeGammaVarRapK0Short/F");
    fTreeGammaV0->Branch("fTreeGammaVarRapLambda",&fTreeGammaVarRapLambda,"fTreeGammaVarRapLambda/F");
    fTreeGammaV0->Branch("fTreeGammaVarInvMassK0s",&fTreeGammaVarInvMassK0s,"fTreeGammaVarInvMassK0s/F");
    fTreeGammaV0->Branch("fTreeGammaVarInvMassLambda",&fTreeGammaVarInvMassLambda,"fTreeGammaVarInvMassLambda/F");
    fTreeGammaV0->Branch("fTreeGammaVarInvMassAntiLambda",&fTreeGammaVarInvMassAntiLambda,"fTreeGammaVarInvMassAntiLambda/F");
    fTreeGammaV0->Branch("fTreeGammaVarV0CosineOfPointingAngle",&fTreeGammaVarV0CosineOfPointingAngle,"fTreeGammaVarV0CosineOfPointingAngle/F");
    fTreeGammaV0->Branch("fTreeGammaVarAlphaV0",&fTreeGammaVarAlphaV0,"fTreeGammaVarAlphaV0/F");
    fTreeGammaV0->Branch("fTreeGammaVarPtArmV0",&fTreeGammaVarPtArmV0,"fTreeGammaVarPtArmV0/F");
    fTreeGammaV0->Branch("fTreeGammaVarLeastNbrCrossedRows",&fTreeGammaVarLeastNbrCrossedRows,"fTreeGammaVarLeastNbrCrossedRows/I");
    fTreeGammaV0->Branch("fTreeGammaVarLeastRatioCrossedRowsOverFindable",&fTreeGammaVarLeastRatioCrossedRowsOverFindable,"fTreeGammaVarLeastRatioCrossedRowsOverFindable/F");
    fTreeGammaV0->Branch("fTreeGammaVarDistOverTotMom",&fTreeGammaVarDistOverTotMom,"fTreeGammaVarDistOverTotMom/F");
    fTreeGammaV0->Branch("fTreeGammaVarNSigmasPosProton",&fTreeGammaVarNSigmasPosProton,"fTreeGammaVarNSigmasPosProton/F");
    fTreeGammaV0->Branch("fTreeGammaVarNSigmasPosPion",&fTreeGammaVarNSigmasPosPion,"fTreeGammaVarNSigmasPosPion/F");
    fTreeGammaV0->Branch("fTreeGammaVarNSigmasNegProton",&fTreeGammaVarNSigmasNegProton,"fTreeGammaVarNSigmasNegProton/F");
    fTreeGammaV0->Branch("fTreeGammaVarNSigmasNegPion",&fTreeGammaVarNSigmasNegPion,"fTreeGammaVarNSigmasNegPion/F");
    fTreeGammaV0->Branch("fTreeGammaVarNegEta",&fTreeGammaVarNegEta,"fTreeGammaVarNegEta/F");
    fTreeGammaV0->Branch("fTreeGammaVarPosEta",&fTreeGammaVarPosEta,"fTreeGammaVarPosEta/F"); 
    fTreeList->Add(fTreeGammaV0);
   */

   //<<<<<------------------------------------

    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    /* Helper
    if(! fPPVsMultUtils ) {
        fPPVsMultUtils = new AliPPVsMultUtils();
    } 
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }
    */


    PostData(1, fOutputContainer);
    //    PostData(2, fTreeV0);

    //  PostData(2, fOutputContainer);
    //  PostData(3, fOutputContainer);
    
    
    //   PostData(2, tSigma0);
    //   PostData(3, tTreeEvent);
  
}

//======================================================
//_____________________________________________________
void AliAnalysisTaskSigma0::UserExec(Option_t */*option*/)
{
    // Called for each event
    
  //  printf("111 Execute analysis for current event... -1 \n");
    // Select conversion and calorimeter photons
    
    //Make sure that old events are added to the list of "previous" and new are cleaned
    //If event is bad lists should be empty
    //We either add current events to the list of "previous" or remove
    //Reject only those events where no photons anywhere
    //However, this will result in wrong absolute normalization.
    //Correct this if it is important for your analysiz

    const Int_t knEventsToKeep=100 ;
    if(fPHOSEvent && fEMCALEvent && fConvEvent && fGenpi0Event  ){
        
        
        if(   fPHOSEvent->GetEntriesFast()>0 || fEMCALEvent->GetEntriesFast()>0
           || fConvEvent->GetEntriesFast()>0 || fGenpi0Event->GetEntriesFast()>0    ){
            
   
            Int_t izvtx = (Int_t)((fZvtx+10.)/2.) ;
            if(izvtx<0)izvtx=0 ;
            if(izvtx>9)izvtx=9 ;
            
            TList * prevPHOS = fPHOSEvents[izvtx] ;
            TList * prevEMCAL = fEMCALEvents[izvtx] ;
            TList * prevConv = fConvEvents[izvtx] ;
            TList * prevLam = fGenpi0Events[izvtx] ;
            prevPHOS->AddFirst(fPHOSEvent) ;
            fPHOSEvent=0;
            prevEMCAL->AddFirst(fEMCALEvent) ;
            fEMCALEvent=0 ;
            prevConv->AddFirst(fConvEvent) ;
            fConvEvent=0 ;
            prevLam->AddFirst(fGenpi0Event) ;
            fGenpi0Event=0 ;
            
            if(prevPHOS->GetSize()>knEventsToKeep){//Remove redundant events
                TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
                prevPHOS->RemoveLast() ;
                delete tmp ;
            }
            if(prevEMCAL->GetSize()>knEventsToKeep){//Remove redundant events
                TClonesArray * tmp = static_cast<TClonesArray*>(prevEMCAL->Last()) ;
                prevEMCAL->RemoveLast() ;
                delete tmp ;
            }
            if(prevConv->GetSize()>knEventsToKeep){//Remove redundant events
                TClonesArray * tmp = static_cast<TClonesArray*>(prevConv->Last()) ;
                prevConv->RemoveLast() ;
                delete tmp ;
            }
            if(prevLam->GetSize()>knEventsToKeep){//Remove redundant events
                TClonesArray * tmp = static_cast<TClonesArray*>(prevLam->Last()) ;
                prevLam->RemoveLast() ;
                delete tmp ;
            }
        }
    }
    
    
    if(fPHOSEvent)
    fPHOSEvent->Clear() ;
    else
    fPHOSEvent = new TClonesArray("AliCaloParticle",10) ;
    if(fEMCALEvent)
    fEMCALEvent->Clear() ;
    else
    fEMCALEvent = new TClonesArray("AliCaloParticle",10) ;
    if(fConvEvent)
    fConvEvent->Clear() ;
    else
    fConvEvent = new TClonesArray("AliCaloParticle",10) ;
    if(fTrackEvent)
    fTrackEvent->Clear() ;
    else
    fTrackEvent = new TClonesArray("AliCaloParticle",10) ;
    
    if(fGenpi0Event)
    fGenpi0Event->Clear() ;
    else
    fGenpi0Event = new TClonesArray("AliCaloParticle",10) ;
    
    
    //First try to find Stack information.
    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
        if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent())
        fStack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
    }
   
    // Connect to the InputEvent
    // Appropriate for ESD analysis !
    
   // AliESDInputHandler *esdHandler=dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   AliESDInputHandler *esdHandler = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  
    
    if(!fESDpid){
        if( esdHandler && esdHandler->GetESDpid()){
            fESDpid=new AliESDpid(*(esdHandler->GetESDpid())) ;
        }
        else {
            fESDpid=new AliESDpid;
            Double_t alephParameters[5];
            if(fStack){// simulation
                alephParameters[0] = 2.15898e+00/50.;
                alephParameters[1] = 1.75295e+01;
                alephParameters[2] = 3.40030e-09;
                alephParameters[3] = 1.96178e+00;
                alephParameters[4] = 3.91720e+00;
                fESDpid->GetTOFResponse().SetTimeResolution(80.);
            }
            else{// data
                alephParameters[0] = 0.0283086;
                alephParameters[1] = 2.63394e+01;
                alephParameters[2] = 5.04114e-11;
                alephParameters[3] = 2.12543e+00;
                alephParameters[4] = 4.88663e+00;
                fESDpid->GetTOFResponse().SetTimeResolution(130.);
                fESDpid->GetTPCResponse().SetMip(47.9);
            }
            
            fESDpid->GetTPCResponse().SetBetheBlochParameters(
                                                              alephParameters[0],alephParameters[1],alephParameters[2],
                                                              alephParameters[3],alephParameters[4]);
            fESDpid->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
        }
    }
    
    if(!fESDtrackCuts){
        
        ///////////////////////////////////////////////
        // Track Cuts for ESD analysis
        fESDtrackCuts = new AliESDtrackCuts();
        
        //    fESDtrackCuts->SetPtRange(.15,1000); upto 18nov13
        
        fESDtrackCuts->SetPtRange(0.001,1000);
        
        fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
        fESDtrackCuts->SetMinNClustersTPC(70);
        fESDtrackCuts->SetRequireTPCRefit(kTRUE);
        
     }
    

    fESDEvent=(AliESDEvent*)InputEvent();
    //    fEventCounter = 0;
     
    // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
    if(fEventCounter == 0) {
        for(Int_t mod=0; mod<5; mod++) {
            if(! fESDEvent->GetPHOSMatrix(mod)) continue;
            fPHOSGeo->SetMisalMatrix( fESDEvent->GetPHOSMatrix(mod),mod) ;
            //    Printf("PHOS geo matrix %p for module # %d is set\n", fESDEvent->GetPHOSMatrix(mod), mod);
        }
    }
    
    FillHistogram("hQAEvents",0.) ; // total number of event
    FillHistogram("hQAEventsTrig",0.) ;

    //    PostData(1, fOutputContainer);    
    //    fEventCounter++;
    //    if( 1>0 ) return;
    
    TString trgClasses = fESDEvent->GetFiredTriggerClasses();
    if (trgClasses.Contains("FAST") && !trgClasses.Contains("ALL")) {
        AliWarning(Form("Skip event with triggers %s",trgClasses.Data() ));
    }
    FillHistogram("hQAEventsTrig",3.) ;
    
    
    // Event selection flags
    Bool_t eventVtxExist    = kFALSE;
    Bool_t eventPileup      = kFALSE;
    
    // Checks if we have a primary vertex  // Get primary vertices form ESD
    if      (fESDEvent->GetPrimaryVertexTracks()->GetNContributors()>0)
    eventVtxExist    = kTRUE;
    else if (fESDEvent->GetPrimaryVertexSPD()   ->GetNContributors()>0)
    eventVtxExist    = kTRUE;
    
    const AliESDVertex *esdVertex5 = fESDEvent->GetPrimaryVertex();
    Double_t vtx5[3];
    vtx5[0] = esdVertex5->GetX();
    vtx5[1] = esdVertex5->GetY();
    vtx5[2] = esdVertex5->GetZ();
    
       
    FillHistogram("hNvertexTracks",esdVertex5->GetNContributors());
    FillHistogram("hZvertex"      ,esdVertex5->GetZ());
    if( !  (TMath::Abs(esdVertex5->GetZ()) < 10. ) ) return;
   
    
    
    //Init geometry if not done yet
    InitGeometry();
  
    
    FillHistogram("hQAEventsTrig",11.) ;

    
    // Check for pileup and fill pileup histograms
    if (fESDEvent->IsPileupFromSPD()) {
        eventPileup = kTRUE;
        FillHistogram("hQAEventsTrig",12.) ;
        return;
    }
    

    
         
    FillHistogram("hQAEventsTrig",13.) ;
 //   FillHistogram("hQAEventsTrig",14) ;
    
    // Get Centrality information
    // fGetCent = fESDEvent->GetCentrality();

    /*
    AliCentrality *cent = fESDEvent->GetCentrality();
    if(cent)  fCentralityV0M = cent->GetCentralityPercentile("V0M");
    fCentr= fCentralityV0M ;
    //   printf(" fcent %f fcentv0 %f \n", fGetCent, fCentralityV0M );
    FillHistogram("hQAMult",fCentr) ;


    Float_t lCentrality = -100;
    Float_t centralityV0M = -100;
    Float_t centralityV0A = -100;
    
    
    //    fCentrality = fESDEvent->GetCentrality();
    centralityV0M = cent->GetCentralityPercentile("V0M");
    centralityV0A = cent->GetCentralityPercentile("V0A");
    lCentrality =  cent->GetCentralityPercentile("V0C");

    FillHistogram("hQACentr1",centralityV0M);
    FillHistogram("hQACentr2",centralityV0A);
    FillHistogram("hQACentr3",lCentrality);
    */


    //printf(" cent %f %f %f \n", centralityV0M, centralityV0A, lCentrality);

    //Calculate charged multiplicity
    Int_t trackCounter = 0;
    for (Int_t i=0;i<fESDEvent->GetNumberOfTracks();++i) {
        AliESDtrack *track = new AliESDtrack(*fESDEvent->GetTrack(i)) ;
        if( fESDtrackCuts->AcceptTrack(track)  ){
            trackCounter++;
        }
        delete track;
    }

    fMultiplicity = trackCounter;       
    FillHistogram("hTrackMult", fMultiplicity) ;
    FillHistogram("hQAEventsTrig",14.) ;
 

    //Fill MC histograms if MC is present
    ProcessMC();

    //    SelectV0 ();

    SelectLambda ();
    

    Int_t nLam=fGenpi0Event->GetEntriesFast() ;
    if( nLam < 1 ) return;    
    FillHistogram("hQAEventsTrig",15.) ;
    FillHistogram("hLamMult", nLam) ;

    //        SelectConvPhotons() ;
    SelectPhotonsFB();
    // SelectGammaV0();

    Int_t nConv=fConvEvent->GetEntriesFast() ;
    if( nConv>0) { 
      FillHistogram("hQAEventsTrig",16.) ;
      FillHistogram("hConvMult", nConv) ;
    }      

    SelectPHOSPhotons();
    Int_t nPHOS=fPHOSEvent->GetEntriesFast() ;
    if( nPHOS>0) {
      FillHistogram("hQAEventsTrig",17.) ;
      FillHistogram("hPHOSMult", nPHOS) ;
    }



    if( nConv < 1 && nPHOS<1  ) return;

    //   printf("========+++Event %d !!!  CONV-1 nConv %d nLam %d \n \n \n", fEventCounter, nConv, nLam );
    FillHistogram("hQAEventsTrig",18) ;

        
    
    PostData(1, fOutputContainer);
    
    fEventCounter++;
    
}
//______________________________________________________________________
void AliAnalysisTaskSigma0::InitGeometry()
{
    //If not done yet, create Geometry for PHOS and EMCAL
    //and read misalignment matrixes from ESD/AOD (AOD not implemented yet)
    //
    
    if(fPHOSgeom && fEMCALgeom){ //already initialized
        return ;
    }
    
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
    if(!esd )
    AliFatal("Can not read geometry matrixes from ESD/AOD: NO ESD") ;
    if(!fPHOSgeom){//reading PHOS matrixes
        fPHOSgeom = new AliPHOSGeoUtils("IHEP","");
        for(Int_t mod=0; mod<5; mod++){
            if  (esd){
                const TGeoHMatrix* m=esd->GetPHOSMatrix(mod) ;
                //  const TGeoHMatrix* m=esd->GetPHOSMisalMatrix(mod) ;
                if(m)
                fPHOSgeom->SetMisalMatrix(m, mod) ;
            }
        }
    }
}
//________________________________________________________________
void AliAnalysisTaskSigma0::SelectPHOSPhotons(){
    //SelectPHOSPhotons Loop over all CaloClusters
    
    const Double_t keMin= fptCut ;  //  0.070 ; // 0.2;    // 0.075 ;  0.1 ???
    const Double_t keMax= fptMaxCut; // 1.500 ; // 0.2;    // 0.075 ;  0.1 ???
    const Int_t knMin=1;   const Int_t knMax=4;

    //vertex
    Double_t vtx[3];
    vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
    vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
    vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();
    char key[55] ;
    Int_t n= fESDEvent->GetNumberOfCaloClusters() ;
    Int_t inPHOS = 0;
    for (Int_t i=0; i<n; i++) {
        AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(i);
        if(!clu->IsPHOS())
        continue ;
        TLorentzVector p ;
        clu ->GetMomentum(p ,vtx);
        if(p.Energy()<keMin)    continue ;
        if(p.Energy()>keMax)    continue ;
        if( !(clu->GetNCells()>=knMin && clu->GetNCells()<=knMax) )    continue ;
        
        Int_t iMod,iX,iZ ;
        Float_t xyz[3] = {0,0,0};
        clu->GetPosition(xyz);   //Global position in ALICE system
        TVector3 global(xyz) ;
        Int_t relid[4] ;
        if(!fPHOSgeom->GlobalPos2RelId(global,relid)){
            printf("PHOS_beyond: x=%f, y=%f, z=%f \n",xyz[0],xyz[1],xyz[2]) ;
            continue ;
        }
        iMod=relid[0] ;
        // abb 11apr13 for 2013 data - skip mod 2  if (iMod == 2) continue ;  // current for esd13cp2 - 2013c data
        iX=relid[2];
        iZ=relid[3] ;
        if(!IsGoodChannel("PHOS",iMod,iX,iZ))  continue ;
        
        //21oct13 for 10dp2
        if( iMod == 1 && ( ( iX == 20 && iZ == 11 ) ||  		       ( iX == 20 && iZ == 12 ) ||
                          ( iX == 21 && iZ == 11 ) ||		       ( iX == 48 && iZ == 54 ) ||
                          ( iX == 51 && iZ == 2 ) ||		       ( iX == 51 && iZ == 3 ) ||
                          ( iX == 52 && iZ == 11 ) ||		       ( iX == 52 && iZ == 12 ) ||
                          ( iX == 60 && iZ == 2 ) ||		       ( iX == 61 && iZ == 2 )     ) ) continue;
        Double_t phostof;
        phostof =  clu->GetTOF() * 1.e9;
   
        Float_t dXmodule[3] = {-2.30, -2.11, -1.53}; // X-shift in local system for module 1,2,3
        Float_t dZmodule[3] = {-0.40, +0.52, +0.80}; // Z-shift in local system for module 1,2,3
        TVector3 local;
        fPHOSgeom->Global2Local(local,global,iMod) ;
        fPHOSgeom->Local2Global(iMod,local.X()+dXmodule[iMod],local.Z()+dZmodule[iMod],global);
        for (Int_t ixyz=0; ixyz<3; ixyz++) xyz[ixyz]=global[ixyz] ;
        clu->SetPosition(xyz) ;
        
        // Ajust absotule calibration per module
        //    Double_t recalib[3] = {0.9942,0.9822,1.0072}; //Offi, is it for 2010 also???
        Double_t recalib[3] = {1.0, 1.0, 1.0};
        p *= recalib[iMod-1];
        //------ YK 04.03.2013
        if(inPHOS>=fPHOSEvent->GetSize())fPHOSEvent->Expand(2*inPHOS) ;
        AliCaloParticle * ph = new((*fPHOSEvent)[inPHOS]) AliCaloParticle(p) ;
        inPHOS++ ;
        Double_t dist = clu->GetEmcCpvDistance() ;
        //    FillHistogram("PHOSdist",dist,clu->E() ) ;
        // try for Yura 14feb14
        if(       clu->E() <  1 &&                  dist < 28.) continue;
        else if ( clu->E() >= 1 &&  clu->E() < 2 && dist < 16.) continue;
        else if ( clu->E() >= 2 &&                  dist < 12.) continue;
        
        ph->SetCPVBit(clu->GetEmcCpvDistance()>10.) ;
        ph->SetCPV2Bit(clu->GetEmcCpvDistance()>15.) ;
        ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
        //    Bool_t closeToBad=(clu->GetDistanceToBadChannel()>fBadDistCutPHOS) ;
        ph->SetModule(iMod) ;
        //Fill QA
	if(clu->E()> keMin  && iMod <= 3  ){
            sprintf(key,"hQA_PHOS_mod%d_soft",iMod) ;
            FillHistogram(key,iX-0.5, iZ-0.5,1.) ;
	}
	FillHistogram("EGamPHOS", clu->E() ) ;

    }
    //QA: number of clusters/event in run //QA: average cluster energy
} // END of SelectPHOSPhotons(){
//______________________________________________
//______________________________________________________
//_____________________________________________________________
//______________________________________________________________________
void AliAnalysisTaskSigma0::SelectV0(){
  //------------------------------------------------
  // Fill V0 Tree as needed
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
  
  Double_t fMinV0Pt = 0;
  Double_t fMaxV0Pt = 100;

  Int_t nv0s = 0;
  nv0s =  fESDEvent->GetNumberOfV0s();

  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
    {   // This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*) fESDEvent)->GetV0(iV0);
      if (!v0) continue;

      Double_t tDecayVertexV0[3];
      v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
      
      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
      Double_t lV0TotalMomentum = TMath::Sqrt( tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

      lPt = v0->Pt();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();
      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

      UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
      UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

      Double_t lMomPos[3];
      v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3];
      v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

      AliESDtrack *pTrack=((AliESDEvent*) fESDEvent)->GetTrack(lKeyPos);
      AliESDtrack *nTrack=((AliESDEvent*) fESDEvent)->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
	Printf("ERROR: Could not retreive one of the daughter track");
	continue;
      }

      //Daughter Eta for Eta selection, afterwards
      fTreeVariableNegEta = nTrack->Eta();
      fTreeVariablePosEta = pTrack->Eta();

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->GetSign() == nTrack->GetSign()) {
	continue;
      }

      //________________________________________________________________________
      // Track quality cuts
      Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
      Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
      fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
      if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
	fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
      
      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;


      if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;

      //GetKinkIndex condition
      if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

      //Findable clusters > 0 condition
      if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;

      //Compute ratio Crossed Rows / Findable clusters
      //Note: above test avoids division by zero!
      Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF() +1.e-10  ));
      Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF() +1.e-10  ));

      fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
      if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
	fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

      //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
      if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
      
      //End track Quality Cuts
      //________________________________________________________________________

    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = fESDEvent->GetPrimaryVertex();
    // const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    //    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

    Double_t  lMagneticField = fESDEvent->GetMagneticField();

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );

      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->GetChi2V0();
      //      printf(" lChi2V0 %f GetChi2V0 %f  \n",lChi2V0,  v0->GetChi2V0() );

      lDcaV0Daughters = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
      lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
      fTreeVariableV0CosineOfPointingAngle=lV0CosineOfPointingAngle;

      // Getting invariant mass infos directly from ESD
      v0->ChangeMassHypothesis(310);
      lInvMassK0s = v0->GetEffMass();
      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();
      lAlphaV0 = v0->AlphaV0();
      lPtArmV0 = v0->PtArmV0();

      fTreeVariablePt = v0->Pt();
      fTreeVariableChi2V0 = lChi2V0;
      fTreeVariableDcaV0ToPrimVertex = lDcaV0ToPrimVertex;
      fTreeVariableDcaV0Daughters = lDcaV0Daughters;
      fTreeVariableV0CosineOfPointingAngle = lV0CosineOfPointingAngle;
      fTreeVariableV0Radius = lV0Radius;
      fTreeVariableDcaPosToPrimVertex = lDcaPosToPrimVertex;
      fTreeVariableDcaNegToPrimVertex = lDcaNegToPrimVertex;
      fTreeVariableInvMassK0s = lInvMassK0s;
      fTreeVariableInvMassLambda = lInvMassLambda;
      fTreeVariableInvMassAntiLambda = lInvMassAntiLambda;
      fTreeVariableRapK0Short = lRapK0Short;
      fTreeVariableRapLambda = lRapLambda;
      fTreeVariableAlphaV0 = lAlphaV0;
      fTreeVariablePtArmV0 = lPtArmV0;

      //Official means of acquiring N-sigmas
      fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
      fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
      fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
      fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
      
      //This requires an Invariant Mass Hypothesis afterwards
      fTreeVariableDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
						);
      fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

      //      printf(" fTreeVariableDistOverTotMom %f \n", fTreeVariableDistOverTotMom );

      //      fTreeV0->Fill();

      /* ABB SKIP Copy Multiplicity information
      fTreeVariableCentV0M = fCentrality_V0M;
      fTreeVariableCentV0A = fCentrality_V0A;
      fTreeVariableCentV0C = fCentrality_V0C;
      fTreeVariableCentV0MEq = fCentrality_V0MEq;
      fTreeVariableCentV0AEq = fCentrality_V0AEq;
      fTreeVariableCentV0CEq = fCentrality_V0CEq;
      fTreeVariableCentV0B = fCentrality_V0B;
      fTreeVariableCentV0Apartial = fCentrality_V0Apartial;
      fTreeVariableCentV0Cpartial = fCentrality_V0Cpartial;
      fTreeVariableCentV0S = fCentrality_V0S;
      fTreeVariableCentV0SB = fCentrality_V0SB;
      fTreeVariableRefMultEta8 = fRefMultEta8;
      fTreeVariableRefMultEta5 = fRefMultEta5;
      fTreeVariableRunNumber = fRunNumber;
      for(Int_t i=0; i<20; i++) fTreeVariableRefMultDiffEta[i] = fRefMultDiffEta[i];   */

      //------------------------------------------------
      // Fill Tree!
      //------------------------------------------------
      
      // The conditionals are meant to decrease excessive
      // memory usage!

      //First Selection: Reject OnFly
      if( lOnFlyStatus == 0 ) {
	//Second Selection: rough 20-sigma band, parametric.
	//K0Short: Enough to parametrize peak broadening with linear function.
	Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
	Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
	//Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
	//[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
	Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + 
	  (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
	Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - 
	  (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
	//Do Selection
	if( (fTreeVariableInvMassLambda     < lUpperLimitLambda  && fTreeVariableInvMassLambda > lLowerLimitLambda     ) ||
	    (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda     ) ||
	    (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short    ) ) {
	  //Pre-selection in case this is AA...
	  if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && fkSaveV0Tree ) fTreeV0->Fill();
	}
      }
    }// This is the end of the V0 loop
  
  //------------------------------------------------
  // Fill V0 tree over.
  //------------------------------------------------
} // END of SelectV0

//______________________________________________________
//_____________________________________________________________
//______________________________________________________________________
void AliAnalysisTaskSigma0::SelectGammaV0(){
  //------------------------------------------------
  // Fill V0 Tree as needed
  //------------------------------------------------

  //Variable definition
  //  Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
  // Double_t lChi2V0 = 0;
  //   Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
  //  Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
  // Double_t lV0CosineOfPointingAngle = 0;
  //  Double_t lV0Radius = 0, lPt = 0;
  //  Double_t  lPt = 0;
  //  Double_t lRapK0Short = 0, lRapLambda = 0;
  // Double_t lInvMassK0s = 0, lInvMassLambda = 0 ; // , lInvMassAntiLambda = 0;
  //  Double_t lAlphaV0 = 0, lPtArmV0 = 0;
  
  //  Double_t fMinV0Pt = 0;
  //  Double_t fMaxV0Pt = 100;

  //  Int_t nv0s = 0;

  if(1>0 ) return;

}


  /* abb 11jan16
  nv0s =  fESDEvent->GetNumberOfV0s();

  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
    {   // This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*) fESDEvent)->GetV0(iV0);
      if (!v0) continue;

      Double_t tDecayVertexV0[3];
      v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
      
      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
      Double_t lV0TotalMomentum = TMath::Sqrt( tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

      lPt = v0->Pt();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();
      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

      UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
      UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

      Double_t lMomPos[3];
      v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3];
      v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

      AliESDtrack *pTrack=((AliESDEvent*) fESDEvent)->GetTrack(lKeyPos);
      AliESDtrack *nTrack=((AliESDEvent*) fESDEvent)->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
	Printf("ERROR: Could not retreive one of the daughter track");
	continue;
      }

      //Daughter Eta for Eta selection, afterwards
      fTreeGammaVarNegEta = nTrack->Eta();
      fTreeGammaVarPosEta = pTrack->Eta();

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->GetSign() == nTrack->GetSign()) {
	continue;
      }

      //________________________________________________________________________
      // Track quality cuts
      Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
      Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
      fTreeGammaVarLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
      if( lNegTrackCrossedRows < fTreeGammaVarLeastNbrCrossedRows )
	fTreeGammaVarLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
      
      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;


      if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;

      //GetKinkIndex condition
      if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

      //Findable clusters > 0 condition
      if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;

      //Compute ratio Crossed Rows / Findable clusters
      //Note: above test avoids division by zero!
      Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
      Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));

      fTreeGammaVarLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
      if( lNegTrackCrossedRowsOverFindable < fTreeGammaVarLeastRatioCrossedRowsOverFindable )
	fTreeGammaVarLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

      //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
      if ( fTreeGammaVarLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
      
      //End track Quality Cuts
      //________________________________________________________________________

    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = fESDEvent->GetPrimaryVertex();
    // const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    //    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

    Double_t  lMagneticField = fESDEvent->GetMagneticField();

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );

      v0->ChangeMassHypothesis(0);
      lInvMassK0s = v0->GetEffMass();
      if( lInvMassK0s > 2 ) continue;

      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->GetChi2V0();
      //      printf(" lChi2V0 %f GetChi2V0 %f  \n",lChi2V0,  v0->GetChi2V0() );

      lDcaV0Daughters = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
      lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
      fTreeGammaVarV0CosineOfPointingAngle=lV0CosineOfPointingAngle;

      // Getting invariant mass infos directly from ESD
      v0->ChangeMassHypothesis(0);
      lInvMassK0s = v0->GetEffMass();
      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();
      lAlphaV0 = v0->AlphaV0();
      lPtArmV0 = v0->PtArmV0();

      fTreeGammaVarPt = v0->Pt();
      fTreeGammaVarChi2V0 = lChi2V0;
      fTreeGammaVarDcaV0ToPrimVertex = lDcaV0ToPrimVertex;
      fTreeGammaVarDcaV0Daughters = lDcaV0Daughters;
      fTreeGammaVarV0CosineOfPointingAngle = lV0CosineOfPointingAngle;
      fTreeGammaVarV0Radius = lV0Radius;
      fTreeGammaVarDcaPosToPrimVertex = lDcaPosToPrimVertex;
      fTreeGammaVarDcaNegToPrimVertex = lDcaNegToPrimVertex;
      fTreeGammaVarInvMassK0s = lInvMassK0s;
      fTreeGammaVarInvMassLambda = lInvMassLambda;
      fTreeGammaVarInvMassAntiLambda = lInvMassAntiLambda;
      fTreeGammaVarRapK0Short = lRapK0Short;
      fTreeGammaVarRapLambda = lRapLambda;
      fTreeGammaVarAlphaV0 = lAlphaV0;
      fTreeGammaVarPtArmV0 = lPtArmV0;

      //Official means of acquiring N-sigmas
      fTreeGammaVarNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
      fTreeGammaVarNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
      fTreeGammaVarNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
      fTreeGammaVarNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
      
      //This requires an Invariant Mass Hypothesis afterwards
      fTreeGammaVarDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
						);
      fTreeGammaVarDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

      // Fill Tree!
      //------------------------------------------------
      
      // The conditionals are meant to decrease excessive
      // memory usage!

      //First Selection: Reject OnFly
      if( lOnFlyStatus == 0 ) {
	//Second Selection: rough 20-sigma band, parametric.
	//K0Short: Enough to parametrize peak broadening with linear function.
	Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeGammaVarPt;
	Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeGammaVarPt;
	//Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
	//[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
	Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeGammaVarPt + 
	  (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeGammaVarPt);
	Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeGammaVarPt - 
	  (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeGammaVarPt);
	//Do Selection
	if( (fTreeGammaVarInvMassLambda     < lUpperLimitLambda  && fTreeGammaVarInvMassLambda > lLowerLimitLambda     ) ||
	    (fTreeGammaVarInvMassAntiLambda < lUpperLimitLambda  && fTreeGammaVarInvMassAntiLambda > lLowerLimitLambda     ) ||
	    (fTreeGammaVarInvMassK0s        < lUpperLimitK0Short && fTreeGammaVarInvMassK0s        > lLowerLimitK0Short    ) ) {
	  //Pre-selection in case this is AA...
	  if ( TMath::Abs(fTreeGammaVarNegEta)<0.8 && TMath::Abs(fTreeGammaVarPosEta)<0.8 && fkSaveV0Tree ) fTreeV0->Fill();
	}
      }
    }// This is the end of the V0 loop
  
  //------------------------------------------------
  // Fill V0 tree over.
  //------------------------------------------------
  } // END of Select  GAMMA V0
  */
//______________________________________________
//______________________________________________________
//_____________________________________________________________
//______________________________________________________________________
void AliAnalysisTaskSigma0::SelectLambda(){     //Fill list of Lambdas
        
  Double_t minPnSigmaAbovePionLine = 1. ;
  Double_t maxPnSigmaAbovePionLine = 3. ;
  Double_t nSigmaAbovePionLine = 0 ;
  const Bool_t useImprovedVertex=kTRUE ;
    
  //No primary vertex in event => scip
  if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) {    return;  }
    
  Double_t  lMagneticField = fESDEvent->GetMagneticField();
    
  const AliESDVertex *esdVertex5 = fESDEvent->GetPrimaryVertex();
  Double_t vtx5[3];
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  Double_t  lPrimaryVtxPosition[3];
  lPrimaryVtxPosition[0] = vtx5[0]; lPrimaryVtxPosition[1] = vtx5[1]; lPrimaryVtxPosition[2] = vtx5[2];
    
  Int_t nV0=fESDEvent->GetNumberOfV0s() ;
  Int_t inLam=0 ;   Double_t ModLam = 0;
    
  Int_t NlamEv = 0;  Int_t NalamEv = 0; Int_t NlamminalamEv = 0;
  Double_t lampospxprev = -100 ;
  Double_t lamnegpxprev = -100 ;
  Double_t alampospxprev = -100 ;
  Double_t alamnegpxprev = -100 ;
  Double_t MassLam = 1.115683 ;
  Double_t MassLamLim = 0.050 ;
  Double_t MassPi = 0.13957018 ;
  Double_t MassPr = 0.938272046 ;
  Double_t mlamprev = 0;  Double_t malamprev = 0;
    
  //-old-upto-25nov15    Double_t ptminLam =  0.25;
  Double_t ptminLam =  0.0;  // as in Lam analysis
  Double_t MminLam = 1.080;       Double_t MmaxLam = 1.160;
  //    Double_t MminLam = 1.110;       Double_t MmaxLam = 1.122;
  
  for(Int_t iv0=0; iv0<nV0;iv0++){    // Main Loop over V0s
    AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;
        
    AliESDtrack * pos = fESDEvent->GetTrack(v0->GetPindex()) ;
    AliESDtrack * neg = fESDEvent->GetTrack(v0->GetNindex()) ;
    const AliExternalTrackParam * paramPos = v0->GetParamP() ;
    const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
    if(pos->GetSign() <0){//change tracks
      pos=neg ;
      neg=fESDEvent->GetTrack(v0->GetPindex()) ;
      paramPos=paramNeg ;
      paramNeg=v0->GetParamP() ;
    }
        
    //cuts ------------    //remove like sign pairs
    if(pos->GetSign() == neg->GetSign())     continue ;                
    if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) continue ;        
    if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) || !(neg->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
             
    Int_t    lOnFlyStatus        = 0;
    Bool_t    isOnFly        = 0;
    // VO's main characteristics to check the reconstruction cuts
    lOnFlyStatus       = v0->GetOnFlyStatus();
    isOnFly = v0->GetOnFlyStatus();

    // if( lOnFlyStatus == 0 ) continue;  // Used up to 14jan16

    if( !(lOnFlyStatus == 0) ) continue;  // reverse cut, check 14jan16

    //      lChi2V0            = v0->GetChi2V0();
    //   Double_t  lDcaV0Daughters    = v0->GetDcaV0Daughters();
    //      lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
    Double_t  lV0cosPointAngle   = v0->GetV0CosineOfPointingAngle(vtx5[0],vtx5[1], vtx5[2]);
        
    Double_t  lV0Position[3];
    v0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);
        
    Double_t  lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
      
        
    // DCA between daughter and Primary Vertex:
    Double_t  lDcaPosToPrimVertex=TMath::Abs(pos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField));
    Double_t lDcaNegToPrimVertex=TMath::Abs(pos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField));
    Double_t lDcaV0Daughters = v0->GetDcaV0Daughters();
    ///////////////values for cuts//////////////////////////////////
    if ((lDcaPosToPrimVertex < 0.1) || (lDcaNegToPrimVertex < 0.1) ||(lV0cosPointAngle < 0.998) 
     || (lV0Radius < 0.0) || (lV0Radius > 100. ) ) continue;
        
    //select V0 finder 
        
    Double_t v0x=0.,v0y=0.,v0z=0.;
    v0->GetXYZ(v0x,v0y,v0z) ;
    //    Double_t r=TMath::Sqrt(v0x*v0x + v0y*v0y) ;
    // printf ("Lambda r %f V0R %f \n", r, lV0Radius);
        
    Int_t LamRecFl = 0;
    // Lambda -> P+ pi-  ---------------
    AliKFParticle negKFKpim(*paramNeg,211);
    AliKFParticle posKFKprot(*paramPos,2212);
    AliKFParticle lamKF(negKFKpim,posKFKprot) ;
    lamKF.SetMassConstraint(MassLam, 0.2 );
        
    if(useImprovedVertex){
      AliKFVertex primaryVertexImproved(*(fESDEvent->GetPrimaryVertex()));
      //if Vtx do created
      if(primaryVertexImproved.GetNContributors()>1){
	primaryVertexImproved+= lamKF;
	lamKF.SetProductionVertex(primaryVertexImproved);
      }
    }
    Double_t mlam=0., widthlam=0., ptlam=0., etalam=0., philam=0 ;
    lamKF.GetMass(mlam, widthlam );
    TLorentzVector lamLV;
    lamLV.SetXYZM( lamKF.GetPx(), lamKF.GetPy(), lamKF.GetPz(),mlam ) ;  //Produces slightly better pi0 width
    ptlam= lamLV.Pt();
    etalam = abs( lamLV.Eta() ) ;
    philam = lamLV.Phi();
        
    // Anti-Lambda -> P- pi+  ---------------
    AliKFParticle negKFKapim(*paramNeg,2212);
    AliKFParticle posKFKaprot(*paramPos,211);
    AliKFParticle alamKF(negKFKapim,posKFKaprot) ;
    alamKF.SetMassConstraint(MassLam, 0.2 );
    //    alamKF.SetModule(-1);
    if(useImprovedVertex){
      AliKFVertex primaryVertexImproved(*(fESDEvent->GetPrimaryVertex()));
      //if Vtx do created
      if(primaryVertexImproved.GetNContributors()>1){
	primaryVertexImproved+= alamKF;
	alamKF.SetProductionVertex(primaryVertexImproved);
      }
    }
    Double_t malam=0., widthalam=0., ptalam=0., etaalam=0., phialam=0 ;
    alamKF.GetMass(malam, widthalam );
    TLorentzVector alamLV;
    alamLV.SetXYZM( alamKF.GetPx(), alamKF.GetPy(), alamKF.GetPz(), malam ) ;
    ptalam= alamLV.Pt();
    etaalam = abs( alamLV.Eta() );
    phialam = alamLV.Phi();
        
    //    FillHistogram("R2Conv", r ) ;
    
    Bool_t sigmaPionPos=kTRUE;    Bool_t sigmaPionNeg=kTRUE;
        
    if(pos->P()>minPnSigmaAbovePionLine && pos->P()<maxPnSigmaAbovePionLine ){
      if(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)<nSigmaAbovePionLine){
	sigmaPionPos=kFALSE;
	//  //          continue ;
      }
    }
    if(neg->P()>minPnSigmaAbovePionLine && neg->P()<maxPnSigmaAbovePionLine){
      if(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)<nSigmaAbovePionLine){
	sigmaPionNeg=kFALSE;
	// 	//          continue ;
      }
    }
        
    if( (mlam > MminLam &&  mlam < MmaxLam && ptlam > ptminLam )  ) {
                
      FillHistogram("hLamdEdxTrack", pos->GetP(), pos->GetTPCsignal()) ;
      FillHistogram("hPhi2LamMvsPt11",mlam, ptlam );
      FillHistogram("hRLam2", lV0Radius );
      // FillHistogram("hLamPtPvsM11",pos->Pt(), neg->Pt() );
      // printf("Lam-v0 m %f pt %f \n", mlam, ptlam );
	  
      FillHistogram("hRecPLam",ptlam );
      if (  etalam < fetaCut ) 	 FillHistogram("hRecPLamEta1",ptlam );
      //	 Int_t ilam = lamKF.GetLabel();
      //	 if (  etalam < fEtaCuts[1] )   printf("ptlam %f \n", ptlam );
    }
            
    if ( (malam > MminLam && malam < MmaxLam && ptalam> ptminLam )  ) {
      FillHistogram("hAntiLamdEdxTrack", neg->GetP(), neg->GetTPCsignal()) ;
      FillHistogram("hPhi2ALamMvsPt11",malam, ptalam );
      FillHistogram("hRALam2", lV0Radius );
      // FillHistogram("hALamPtPvsM11",pos->Pt(), neg->Pt() );
      //	 printf("ALam-v0 m %f pt %f \n", malam, ptalam );
      FillHistogram("hRecALam",ptalam );
      if ( etaalam  < fetaCut ) 	FillHistogram("hRecALamEta1",ptalam );
      // Int_t ialam = alamKF.GetLabel() ;
      //	 printf("ilam %d \n", ialam );
    }
            
    FillHistogram("hLamMvsPt10",mlam, ptlam );
    FillHistogram("hALamMvsPt10",malam, ptalam );
    // FillHistogram("R3Conv", r ) ;
    
    //    Double_t Rlam0 =   sqrt( v0->Xv()*v0->Xv() +  v0->Yv()*v0->Yv() + v0->Zv()*v0->Zv() );
    //    FillHistogram("hMClam0MassPt0All", v0->GetEffMass(),  v0->Pt() );
    // FillHistogram("hMClam0RPt0All", Rlam0,  v0->Pt() );
            
    
    ModLam = 0;
    if  ( (mlam > MminLam && mlam < MmaxLam && ptlam > ptminLam ) ||
	  (malam > MminLam && malam < MmaxLam && ptalam> ptminLam )
	  ) {                     // PLam and ALam fill into AliCaloParticle
      FillHistogram("hPhi1LamMvsPt11",mlam, ptlam );
      FillHistogram("hPhi1ALamMvsPt11",malam, ptalam );
	  
      Double_t posp[3]= { pos->Px(),  pos->Py(),  pos->Pz() };
      Double_t negp[3]= { neg->Px(),  neg->Py(),  neg->Pz() };
      Double_t moth[3]= { lamKF.GetPx(), lamKF.GetPy(), lamKF.GetPz() };
      Double_t arpod[2]= {0,0};
      GetArPod( posp, negp, moth, arpod );
      //	  printf("Arpod %f %f mlam %f malam %f \n", arpod[1],arpod[2],mlam, malam);
                
      if(inLam>=fGenpi0Event->GetSize())fGenpi0Event->Expand(inLam*2) ;
      
      if( (mlam > MminLam  && mlam < MmaxLam && ptlam > ptminLam )   &&
	  ( (arpod[1] > 0.2 && arpod[1] < 0.9) && (arpod[0] > 0.01 && arpod[0] < 0.17)) ) {
	
	AliCaloParticle * ph = new((*fGenpi0Event)[inLam]  ) AliCaloParticle(lamLV) ;
	    
	ModLam = 100.;
	ph->SetXt(ModLam);   // Xt will be renamed
	ph->SetIsolation(iv0);

	ph->SetPprx(pos->Px());
	ph->SetPpry(pos->Py());
	ph->SetPprz(pos->Pz());
	    
	ph->SetPpix(neg->Px());
	ph->SetPpiy(neg->Py());
	ph->SetPpiz(neg->Pz());
	    
	ph->SetV0CosPointingAngle(lV0cosPointAngle);
	ph->SetV0DCADaughters(lDcaV0Daughters);
	ph->SetV0IPNeg(lDcaNegToPrimVertex);
	ph->SetV0IPPos(lDcaPosToPrimVertex);
	ph->SetV0Radius(lV0Radius);
                    
	ph->SetArmPt(arpod[0] ) ;
	ph->SetArmAlpha(arpod[1]);
                    
                    
	NlamEv++;
	FillHistogram("hArPodLam",  arpod[1], arpod[0] ) ;
	lampospxprev = pos->Px(); lamnegpxprev = neg->Px() ; mlamprev = mlam;
	inLam++ ;
	    
      }
      if( (malam > MminLam && malam < MmaxLam && ptalam > ptminLam ) &&
	  ( (arpod[1] < -0.2 && arpod[1] > -0.9) && (arpod[0] > 0.01 && arpod[0] < 0.17) )  ) {
		  
	AliCaloParticle * ph = new((*fGenpi0Event)[inLam]  ) AliCaloParticle(alamLV) ;
	ModLam = -100.;
	ph->SetXt(ModLam);
	    
	ph->SetIsolation(iv0);
	
	ph->SetPpix(pos->Px());
	ph->SetPpiy(pos->Py());
	ph->SetPpiz(pos->Pz());
	    
	ph->SetPprx(neg->Px());
	ph->SetPpry(neg->Py());
	ph->SetPprz(neg->Pz());
	    
	ph->SetV0CosPointingAngle(lV0cosPointAngle);
	ph->SetV0DCADaughters(lDcaV0Daughters);
	ph->SetV0IPNeg(lDcaNegToPrimVertex);
	ph->SetV0IPPos(lDcaPosToPrimVertex);
	ph->SetV0Radius(lV0Radius);
	ph->SetV0Chi2(v0->GetChi2V0());
                    
	alampospxprev = pos->Px() ; alamnegpxprev = neg->Px() ; malamprev = malam;
		    
	inLam++ ;
	NalamEv++;
	FillHistogram("hArPodALam",  arpod[1], arpod[0] ) ;
	
      }

      if( ModLam == 100 ) {
	FillHistogram("hLamMvsPt11",mlam, ptlam );
	FillHistogram("hTotProtdEdxTrack", pos->GetP(), pos->GetTPCsignal()) ;
	FillHistogram("hTotPiondEdxTrack", neg->GetP(), neg->GetTPCsignal()) ;
      }
      else if ( ModLam == -100 ) {
	FillHistogram("hALamMvsPt11",malam, ptalam );
	FillHistogram("hTotPiondEdxTrack", pos->GetP(), pos->GetTPCsignal()) ;
	FillHistogram("hTotProtdEdxTrack", neg->GetP(), neg->GetTPCsignal()) ;
      }

      //>>>>>>>>>>>>>>>>>      
      if (ModLam != 0 ) {  // fill TreeV0 a la DDCh
	Double_t tDecayVertexV0[3];
	v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
      
	Double_t tV0mom[3];
	v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
	Double_t lV0TotalMomentum = TMath::Sqrt( tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

	lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

	//	Double_t lPt = v0->Pt();
	Double_t lRapK0Short = v0->RapK0Short();
	Double_t lRapLambda  = v0->RapLambda();
	//      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

	UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
	UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

	Double_t lMomPos[3];
	v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
	Double_t lMomNeg[3];
	v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

	AliESDtrack *pTrack=((AliESDEvent*) fESDEvent)->GetTrack(lKeyPos);
	AliESDtrack *nTrack=((AliESDEvent*) fESDEvent)->GetTrack(lKeyNeg);
	if (!pTrack || !nTrack) {
	  Printf("ERROR: Could not retreive one of the daughter track");
	  continue;
	}

	//Daughter Eta for Eta selection, afterwards
	fTreeVariableNegEta = nTrack->Eta();
	fTreeVariablePosEta = pTrack->Eta();

	// Filter like-sign V0 (next: add counter and distribution)
	//  if ( pTrack->GetSign() == nTrack->GetSign())	continue;

	//________________________________________________________________________
	// Track quality cuts
	Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
	Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
	fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
	if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
	  fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
      
	//CHECK !!! TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
	// if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
	// if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
	
	//      if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;

	//GetKinkIndex condition
	//      if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

	//Findable clusters > 0 condition
	//      if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;

	//Compute ratio Crossed Rows / Findable clusters
	//Note: above test avoids division by zero!
	Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
	Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));

	fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
	if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
	  fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

	//Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
	//      if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
	
	//End track Quality Cuts

	//classical Proton-proton like selection
	const AliESDVertex *lPrimaryBestESDVtx     = fESDEvent->GetPrimaryVertex();
	// const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
	//    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

	Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
	lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
	
	Double_t  lMagneticField = fESDEvent->GetMagneticField();
	
	lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );
	
	lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );

	//	lOnFlyStatus = v0->GetOnFlyStatus();
	Double_t lChi2V0 = v0->GetChi2V0();
	//      printf(" lChi2V0 %f GetChi2V0 %f  \n",lChi2V0,  v0->GetChi2V0() );

	//	lDcaV0Daughters = v0->GetDcaV0Daughters();
	Double_t lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
	Double_t lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
	fTreeVariableV0CosineOfPointingAngle=lV0CosineOfPointingAngle;

	// Getting invariant mass infos directly from ESD
	v0->ChangeMassHypothesis(310);
	Double_t lInvMassK0s = v0->GetEffMass();
	v0->ChangeMassHypothesis(3122);
	Double_t lInvMassLambda = v0->GetEffMass();
	v0->ChangeMassHypothesis(-3122);
      	Double_t lInvMassAntiLambda = v0->GetEffMass();
     	Double_t lAlphaV0 = v0->AlphaV0();
     	Double_t lPtArmV0 = v0->PtArmV0();

	fTreeVariablePt = v0->Pt();
	fTreeVariableChi2V0 = lChi2V0;
	fTreeVariableDcaV0ToPrimVertex = lDcaV0ToPrimVertex;
	fTreeVariableDcaV0Daughters = lDcaV0Daughters;
	fTreeVariableV0CosineOfPointingAngle = lV0CosineOfPointingAngle;
	fTreeVariableV0Radius = lV0Radius;
	fTreeVariableDcaPosToPrimVertex = lDcaPosToPrimVertex;
	fTreeVariableDcaNegToPrimVertex = lDcaNegToPrimVertex;
	fTreeVariableInvMassK0s = lInvMassK0s;
	fTreeVariableInvMassLambda = lInvMassLambda;
	fTreeVariableInvMassAntiLambda = lInvMassAntiLambda;
	fTreeVariableRapK0Short = lRapK0Short;
	fTreeVariableRapLambda = lRapLambda;
	fTreeVariableAlphaV0 = lAlphaV0;
	fTreeVariablePtArmV0 = lPtArmV0;

	//Official means of acquiring N-sigmas
	fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
	fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
	fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
	fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
      
	//This requires an Invariant Mass Hypothesis afterwards
	fTreeVariableDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
						);
	fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

	//      printf(" fTreeVariableDistOverTotMom %f \n", fTreeVariableDistOverTotMom );
	//      fTreeV0->Fill();

	// Fill Tree!
	// The conditionals are meant to decrease excessive memory usage!

	//First Selection: Reject OnFly
	if( lOnFlyStatus == 0 ) {
	  //Second Selection: rough 20-sigma band, parametric.
	  //K0Short: Enough to parametrize peak broadening with linear function.
	  Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
	  Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
	  //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
	  //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
	  Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + 
	    (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
	  Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - 
	    (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
	  //Do Selection
	  if( (fTreeVariableInvMassLambda     < lUpperLimitLambda  && fTreeVariableInvMassLambda > lLowerLimitLambda     ) ||
	      (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda     ) ||
	      (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short    ) ) {
	    //Pre-selection in case this is AA...
	    if ( 1<0 && TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8  ) fTreeV0->Fill(); // fill a lot!
	  }
	}
      }  // END fill TreeV0 a la DDCh
      
    } // end fill PLam or ALam into AliCaloParticle 

    FillHistogram("hNlamEv",  NlamEv );
    FillHistogram("hNalamEv",  NalamEv );
    NlamminalamEv = NlamEv - NalamEv;
    FillHistogram("hNlamminalamEv", NlamminalamEv  );
    //	printf("\n ModLam %f \n", ModLam );            
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    //Fill/check MC information for Lambda
    if( fStack && abs(ModLam)>0  ){                // trace back works for PHYTHIA-10f6a, not yet for PHOJET - 10f6 !!!
      //abb    if ( 1>0 ) {
      TParticle * negMC = fStack->Particle(TMath::Abs(neg->GetLabel()));
      TParticle * posMC = fStack->Particle(TMath::Abs(pos->GetLabel()));
	  
      //	    printf("MC Lambda mneg %f mpos %f ModLam %f ev %d evcount %d \n",
      // negMC->GetMass(),posMC->GetMass(),ModLam,fESDEvent->GetEventNumberInFile(),fEventCounter);
	  
      if(negMC && posMC){
	    
	Int_t lblMotherPosV0Dghter = negMC->GetFirstMother() ;
	Int_t lblMotherNegV0Dghter = posMC->GetFirstMother();
	    
	if( ! (lblMotherPosV0Dghter == lblMotherNegV0Dghter && lblMotherPosV0Dghter > -1) ) continue;
	    
	//printf("\n MC Lambda Mother neg %d pos %d \n", lblMotherPosV0Dghter , lblMotherNegV0Dghter );
	    
	//	      if(negMC->GetMother(0) != posMC->GetMother(0)) continue ;
	if( (( abs( negMC->GetMass() - 0.938272)< 0.01 && abs( posMC->GetMass() - 0.139570)<0.01   ) ||
	     ( abs( posMC->GetMass() - 0.938272)< 0.01 && abs( negMC->GetMass() - 0.139570)<0.01) ) ) {
                        
	  TParticle * v0Lam = fStack->Particle(negMC->GetMother(0) ) ;
                        
	  if( !v0Lam ) continue;
	  Double_t rapLam = Rapidity(  v0Lam->Pt(),  v0Lam->Pz(),  v0Lam->GetMass() );
                        
	  TLorentzVector Pigen;
	  TLorentzVector Prgen;
                        
	  TLorentzVector Pirec;
	  TLorentzVector Prrec;
	  if( ModLam > 0 ) {
	    Prrec.SetXYZM(pos->Px(),pos->Py(),pos->Pz(), MassPr ) ;
	    Pirec.SetXYZM(neg->Px(),neg->Py(),neg->Pz(), MassPi ) ;
                
	    Prgen.SetXYZM(posMC->Px(),posMC->Py(),posMC->Pz(), MassPr ) ;
	    Pigen.SetXYZM(negMC->Px(),negMC->Py(),negMC->Pz(), MassPi ) ;
	  }
	  else if(  ModLam < 0 ) {
	    Pirec.SetXYZM(pos->Px(),pos->Py(),pos->Pz(), MassPi ) ;
	    Prrec.SetXYZM(neg->Px(),neg->Py(),neg->Pz(), MassPr ) ;
			    
	    Pigen.SetXYZM(posMC->Px(),posMC->Py(),posMC->Pz(), MassPi ) ;
	    Prgen.SetXYZM(negMC->Px(),negMC->Py(),negMC->Pz(), MassPr ) ;
	  }
	      
	  //printf("recv0-lam pt %f Lamsim pt %f Alam pt %f evcount %d \n \n",v0Lam->Pt(),ptlam,ptalam,fEventCounter);
	      
	  if(  (v0Lam && abs( v0Lam->GetMass() - MassLam )< MassLamLim )
	       // && ( (abs( v0Lam->Pt() - ptlam) < 0.050  ) || (abs( v0Lam->Pt() - ptalam)< 0.050) )
	       ) {
	    if( ModLam > 0 ) {
	      FillHistogram("hMCgenrec2LamdMvsdPt", v0Lam->GetMass()-lamLV.M(),  v0Lam->Pt()- lamLV.Pt() );
	      FillHistogram("hRec2PLam", v0Lam->Pt() );
	      if( rapLam < fetaCut )   FillHistogram("hRec2PLamEta1", v0Lam->Pt() );
	    }
                            
	    else if ( ModLam<0 ) {
	      FillHistogram("hMCgenrec2ALamdMvsdPt",v0Lam->GetMass()-alamLV.M(),v0Lam->Pt()-alamLV.Pt());
	      FillHistogram("hRec2ALam", v0Lam->Pt() );
	      if( rapLam < fetaCut )   FillHistogram("hRec2ALamEta1", v0Lam->Pt() );
	    }
                                                        
	    TParticle * v0Sig = fStack->Particle(v0Lam->GetMother(0) );
	    //		  printf("mSig %f \n",  v0Sig->GetMass() );
                            
	    Int_t NDt =v0Sig->GetNDaughters();
	    if ( NDt != 2  ||  v0Lam->IsPrimary() ) continue;
	    
	    if( abs( v0Sig->GetPdgCode() ) == 3212  && NDt == 2)  {    // Sigma0
		  
	      if(  v0Sig && abs( v0Sig->GetMass() - 1.192642)< 0.01 ) {
		    
		Int_t Daughter1  = v0Sig->GetFirstDaughter();
		Int_t Daughter2  = v0Sig->GetLastDaughter();
                                    
		Double_t m=0. , mcalc=0. ;
		m =      v0Sig->GetMass();
		mcalc = v0Sig->GetCalcMass();
		// printf("S0 m %f mc %f \n", m, mcalc  );
		//	    printf(" m %f w %f cm %f PDG %d \n", m, width, mcalc, particle->GetPdgCode()  );
		    
		FillHistogram("hMCgenrecSig0",  v0Sig->Pt() );
		Double_t rap = Rapidity(  v0Sig->Pt(),  v0Sig->Pz(),  v0Sig->GetMass() );
		FillHistogram("hMCgenrecSig0Rap", rap );
		//		FillHistogram("hMCgenrecSig0RapPt", rap, v0Sig->Pt()  );
                                    
		//	   Float_t rapLam = Rapidity( v0Lam->Pt(), v0Lam->Pz(), v0Lam->GetMass() );  // 25nov doubled defined
		FillHistogram("hRecLamSig0Rap", rapLam );
		if( abs(rapLam) < fetaCut ) {
		  FillHistogram("hMCgenrecSig0Eta1",  v0Sig->Pt() );
		  FillHistogram("hRecLamSig0Eta1",  v0Lam->Pt() );
		  if(  ModLam > 0 ) FillHistogram("hRecLamSig0",  lamLV.Pt() );
		  else if ( ModLam < 0 )  FillHistogram("hRecLamSig0",  alamLV.Pt() );
		  
		  LamRecFl = 1;
		}
		// printf(" MC Lam 2 mSig %f \n",  v0Sig->GetMass() );
		TParticle* daught1 = (TParticle *)fStack->Particle( Daughter1 );
		TParticle* daught2 = (TParticle *)fStack->Particle( Daughter2 );
		//printf(" MC Lam 3 mSig %f d1-pdg %d d2-pdg %d \n",v0Sig->GetMass(),daught1->GetPdgCode(),daught2->GetPdgCode());
		    
		if ( abs( daught1->GetPdgCode() ) == 3122 && abs( daught2->GetPdgCode())  == 22 ) {
		  // FillHistogram("hMCrecSig0PtLamdaVSGamma0",daught1->Pt(), daught2->Pt()  );
		  FillHistogram("hRec4LamSig0",  daught1->Pt() );
		  // if(  abs( v0Sig->Eta() ) < 1 )  FillHistogram("hRec4LamSig0Eta1", daught1->Pt() );
		  if(  abs(rap) < fetaCut ) FillHistogram("hRec4LamSig0Eta1",  daught1->Pt() );
		  // abb-19may15- see histos up to here
		}
	      }
	    }   // generated Sig0 from recostructed Lam
	  }  // 1 more if to enshure that v0Lam as mother exists!!!
	}
      }
    }   // end fstack - excustion into MC
            

  }  // end of look over V0s

  // Look on Lam+ALam inv mass
  /* if ( NlamEv > 0 && NalamEv>0 ) {
    for(Int_t ilam1 = 0; ilam1 < inLam ; ilam1++){
      AliCaloParticle * Lam1 = static_cast<AliCaloParticle*>(fGenpi0Event->At(ilam1)) ;
      Double_t LambdaMod1 = Lam1->GetXt();
      for(Int_t ilam2 = ilam2+1; ilam2 <= inLam ; ilam2++){
	AliCaloParticle * Lam2 = static_cast<AliCaloParticle*>(fGenpi0Event->At(ilam2)) ;
	Double_t LambdaMod2 = Lam2->GetXt();
	if (LambdaMod1 == LambdaMod2 ) continue;  

	TLorentzVector LamALam= *Lam1 + *Lam2 ;
	Double_t  m=LamALam.M() ;
	Double_t pt=LamALam.Pt() ;
	FillHistogram("hLamALam_MvsPt",m,pt ) ;

       	if( LambdaMod1> 0 ) {
	  FillHistogram("hLamA_MvsPt", Lam1->M(), Lam1->Pt()  ) ;
 	  FillHistogram("hALamA_MvsPt", Lam2->M(),Lam2->Pt()  ) ;
	}
	else {
	  FillHistogram("hALamA_MvsPt", Lam1->M(),Lam1->Pt()  ) ;
 	  FillHistogram("hLamA_MvsPt", Lam2->M(),Lam2->Pt()  ) ;
	}  
      }
    }
  }  // end lool over LAM-ALAM 
  */

} // <<<<<<<<<<<<<<<<< End of SelectLambda, checked 8jan16

//________________________________________________________________________
void AliAnalysisTaskSigma0::SelectPhotonsFB()
{       
  Double_t MassPi = 0.13957018 ;
  Double_t MassPr = 0.938272046 ;
  Double_t MassSig0 = 1.192642 ;

  Double_t SigmaMassWindow = 0.050 ;

  Int_t inConv=0 ;    
  FillHistogram("hCheckFB",0.5);

  fV0Reader=0;
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");

  //   printf("\n \n  PhotonsFB NEW EVENT VOREADER %d \n", fEventCounter );

  if(fStack) fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = kFALSE;
  if(fStack)  fIsMC = kTRUE;

  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader // origianl - can not get it!!!

  //    if(( AliConversionCuts*) fV0Reader->GetConversionCuts() ) ;  //ABB NEW place 15may15
        
  fReaderGammas = fV0Reader->GetReconstructedGammas();    
  Int_t nRecGammas = fReaderGammas->GetEntriesFast();
  FillHistogram("hRecGammas",nRecGammas);


  //25sep16-abb  if(fV0Reader)  if( (AliConversionCuts*) fV0Reader->GetConversionCuts() );
  fV0Reader->GetConversionCuts();

  //    if(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms()) 
  //  fOutputContainer->Add(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
 
  FillHistogram("hCheckFB",1.5);
    
  //     printf("PhotonsFB2  VOREADER %d \n", fReaderGammas->GetEntriesFast() );    
  //    for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
     
  for(Int_t i = 0; i <= nRecGammas ; i++){
    //      printf("--PhotonsFB VOREADER after RecGammas\n");
    FillHistogram("hCheckFB",2.5);
        
    // AliKFConversionPhoton* PhotonCandidate = (AliKFConversionPhoton*) fReaderGammas->At(i);
    
    AliAODConversionPhoton* PhotonCandidate =dynamic_cast<AliAODConversionPhoton*>(fReaderGammas->At(i));

    if(!PhotonCandidate) continue;
    //      printf("--PhotonsFB VOREADER after PhotonCandidate \n");
    //      Double_t PAngle  =    GetCosineOfPointingAngle(PhotonCandidate, fESDEvent) ;
    //      printf("Pangle %f  \n", PAngle);

    TLorentzVector photFB;
    photFB.SetXYZM(PhotonCandidate->GetPx(),PhotonCandidate->GetPy(), PhotonCandidate->GetPz(),0.) ;
    
    // fGammaMass = 
    fGammaPx = PhotonCandidate->GetPx() ;
    fGammaPy = PhotonCandidate->GetPy() ;
    fGammaPz = PhotonCandidate->GetPz() ;
    // fGammaCosPointingAngle = 
    fGammaEta = PhotonCandidate->GetPhotonEta() ; //fEtaPhoton
    //  fGammaArmPt(0) = 
    //  fGammaArmAlpha(0),
    // fGammaZConv(0) = 

    //      fGammaChi2 =  PhotonCandidate->GetPhotonQuality(); // iCatPhoton
    fGammaRadius =  PhotonCandidate->GetConversionRadius();  // fRConvPhoton  
    fGammaDCAzToPrimVtx = PhotonCandidate->GetDCAzToPrimVtx();
    fGammaChi2 = PhotonCandidate->GetChi2perNDF();
    Double_t gammaPt =  PhotonCandidate->Pt();
    fGammaArmPt = gammaPt;

    // printf("pt %f Qual %d R %f DCAz %f \n",  PhotonCandidate->Pt(),  PhotonCandidate->GetPhotonQuality(),  PhotonCandidate->GetConversionRadius(), PhotonCandidate->GetDCAzToPrimVtx() ); 
    // printf("pt %f chi2 %d R %f DCAz %f \n",  gammaPt,  fGammaChi2, fGammaRadius,  fGammaDCAzToPrimVtx);

    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = fESDEvent->GetPrimaryVertex();

    // const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    //    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

    //    Float_t dca = 1000;
    //     PhotonCandidate->GetDistanceOfClossetApproachToPrimVtx(lBestPrimaryVtxPos,dca);

    /*    Double_t  lMagneticField = fESDEvent->GetMagneticField();
	  lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) );
						    lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
						    lBestPrimaryVtxPos[1],
						    lMagneticField) ); */

    //    printf(" lBestPrimaryVtxPos[1] %f dca %f \n", lBestPrimaryVtxPos[1], dca );

    FillHistogram("hCheckFB",3.5);        
    FillHistogram("hAMbeforeFB",PhotonCandidate->GetArmenterosAlpha(), PhotonCandidate->GetArmenterosQt());

    
    //ABB?	ProcessMCParticles();
    //	if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(PhotonCandidate , fStack,kTRUE )) ;
        
    if(inConv>=fConvEvent->GetSize())fConvEvent->Expand(inConv*2) ;
    new((*fConvEvent)[inConv]) AliCaloParticle(photFB) ;

    //      phg->SetV0CosPointingAngle(lV0cosPointAngle);
    /*      photFB->SetV0Radius( PhotonCandidate.GetConversionRadius();  
	    photFB->SetArmPt(  PhotonCandidate.GetArmenterosQt()  );
	    photFB->SetArmAlpha( PhotonCandidate.GetArmenterosAlpha()  );
	    photFB->SetV0Chi2(  PhotonCandidate.GetPhotonQuality()  );
    */
    //      phg->SetZConv(v0z);      

    inConv++ ;
    //	printf("---VOREADER in conv %d \n", inConv );      
    FillHistogram("hCheckFB",4.5);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Fill MC information for GammaConverted
    if(fStack){

      AliStack *fMCStack= fMCEvent->Stack();
	
      //25sep16-abb      if( (AliConversionCuts*) fV0Reader->GetConversionCuts() );
	
      fV0Reader->GetConversionCuts();

      //Int_t labelp=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelPositive())->GetLabel());
      //Int_t labeln=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelNegative())->GetLabel());
	
      Int_t labelp = PhotonCandidate->GetMCLabelPositive() ;
      Int_t labeln = PhotonCandidate->GetMCLabelNegative() ;
      
      if( abs( labelp - labeln) != 1 ) continue;
	
      FillHistogram("hCheckFB",5.5);

      //	printf("PhotonsFB:  fStack done label pos  %d neg %d \n",  labelp , labeln);	
     
      TParticle *negativeMC = fMCStack->Particle(labeln);
      TParticle *positiveMC = fMCStack->Particle(labelp);

      //		if(fPositiveMCParticle&&fNegativeMCParticle){
      //			fCurrentMotherKF->SetMCLabelPositive(labelp);
      //	fCurrentMotherKF->SetMCLabelNegative(labeln);		}
      
      /*  		Int_t MClabelPos = PhotonCandidate->GetMCLabelPositive() ;
			Int_t MClabelNeg = PhotonCandidate->GetMCLabelNegative() ;
			//abb wrong way, 18may15
			TParticle * negativeMC = fStack->Particle( MClabelNeg  );
			TParticle * positiveMC = fStack->Particle( MClabelPos );
      */

      //	    printf("M-neg %f  M-pos %f \n",  negativeMC->GetMass(), positiveMC->GetMass()  );

      if(negativeMC && positiveMC){
	if(negativeMC->GetMother(0) != positiveMC->GetMother(0))
	  continue ;
      }
	    
      if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
	continue;
      }
      if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
	continue;
      }
               
      //           AliKFParticle gamKF( PhotonCandidate, 5);
      TParticle * GamGen = fStack->Particle( negativeMC->GetMother(0)  );  // Generatred Gamma
      
      FillHistogram("hCheckFB",6.5);
      //      FillHistogram("hMCV0MvsWidth", GamGen->GetMass(), 0. ) ;
                
      if(negativeMC->GetUniqueID() != 5 || positiveMC->GetUniqueID() !=5) continue ; // id 5 is conversion
      // FillHistogram("hMCV0MvsWidth1",  GamGen->GetMass(), 0. ) ;

      //	    printf(" PDG gen %d  MClabelPos  %d  \n", GamGen->GetPdgCode(), MClabelPos );

      if(GamGen->GetPdgCode() == 22){    // indeed this Gamma reconstracted from Gamma
	      
	// FillHistogram("hMCV0MvsWidth2",  GamGen->GetMass(), 0. ) ;  
	//!!! add histo	      
	/* Yura 12.3.14 primary Gamma has no Mother ! */
	if( GamGen->IsPrimary() ) continue;

	FillHistogram("hCheckFB",7.5);
	
	// FillHistogram("hMCV0MvsWidth3",  GamGen->GetMass(), 0. ) ;  

	TParticle * v0Sig2 = fStack->Particle(GamGen->GetMother(0));
	      
	// FillHistogram("hMCV0MvsWidth5",  v0Sig2->GetMass(), 0. );

	//	      printf("=========================mSig %f \n",  v0Sig2->GetMass() );
	// abb-19may15 see a few histos above, but not below => take wrong pos-neg e+e- from MC  
	// if( abs(v0Sig2->GetMass() - 1.192642) < 0.05 ) {

	if(  abs(v0Sig2->GetPdgCode() ) == 3212) {
	  
	  // FillHistogram("hMCV0MvsWidth4",   GamGen->GetMass(), 0. ) ;  		
	  // FillHistogram("hMCV0MvsWidth6",  v0Sig2->GetMass(), 0. );
		
	  Int_t NDt1 =v0Sig2->GetNDaughters();
                
	  if( abs( v0Sig2->GetPdgCode() ) == 3212  && NDt1 == 2)  {    // Sigma0
	    // if ( 1>0 ) {

	    FillHistogram("hCheckFB",8.5);
	    
	    // FillHistogram("hMCV0MvsWidth7",  v0Sig2->GetMass(), 0. );

	    Int_t Daughter21  = v0Sig2->GetFirstDaughter();
	    Int_t Daughter22  = v0Sig2->GetLastDaughter();
	    Double_t mrec=0.,  mcalc=0. ;
	    mrec =     v0Sig2->GetMass();
	    mcalc = v0Sig2->GetCalcMass();
	    //	      	    printf("---- m %f w %f cm %f  \n", m, width, mcalc  );
	    FillHistogram("hMCgenrec2Sig0",  v0Sig2->Pt() );
	    Double_t rap = Rapidity(  v0Sig2->Pt(),  v0Sig2->Pz(),  v0Sig2->GetMass() );
	    FillHistogram("hMCgenrec2Sig0Rap", rap );
	    FillHistogram("hRecGamSig0",  GamGen->Pt() );
	    
	    FillHistogram("hMCgenrec2Sig0RapPt", rap,  v0Sig2->Pt() );
	    if( abs(rap) < fetaCut ) {
	      FillHistogram("hMCgenrec2Sig0Eta1",  v0Sig2->Pt() );
	      FillHistogram("hRecGamSig0Eta1",  GamGen->Pt() );
	    }
		  
	    TParticle* LamGen21 = (TParticle *)fStack->Particle( Daughter21 );
	    //	  Int_t labelLp21 = LamGen21->GetMCLabelPositive() ;
	    //		  Int_t labelLn21 = LamGen21->GetMCLabelNegative() ;	

	    TLorentzVector LGen21;
	    LGen21.SetXYZM(LamGen21->Px(),LamGen21->Py(),LamGen21->Pz(),LamGen21->GetCalcMass()) ;  //Produces slightly better pi0 width
	    TLorentzVector LGen21m;
	    LGen21m.SetXYZM(LamGen21->Px(),LamGen21->Py(),LamGen21->Pz(),LamGen21->GetMass()) ;  //Produces slightly better pi0 width
		  
                            
	    TParticle* GamGen22 = (TParticle *)fStack->Particle( Daughter22 );
	    TLorentzVector GGen22;
	    GGen22.SetXYZM(GamGen22->Px(),GamGen22->Py(),GamGen22->Pz(),GamGen22->GetCalcMass()) ;  //Produces slightly better pi0 width
	    TLorentzVector GGen22m;
	    GGen22m.SetXYZM(GamGen22->Px(),GamGen22->Py(),GamGen22->Pz(),GamGen22->GetMass()) ;  //Produces slightly better pi0 width
                            
	    //	  if ( (abs( LamGen21->GetPdgCode() ) == 3122 && abs( GamGen22->GetPdgCode())  == 22)
	    //old-skip		 &&  (inLam >0 && LamRecFl > 0)
	    //     ) {

	    FillHistogram("hCheckFB",9.5);		    
	    if ( 1>0 ) {
	      // Search here Sig0 - Loop over all Lambda reconstr in the event
	      Int_t nLam=fGenpi0Event->GetEntriesFast() ;
	      for(Int_t icon2=0; icon2<nLam; icon2++){ 
		
		FillHistogram("hCheckFB",10.5);	
		AliCaloParticle * LamRec= static_cast<AliCaloParticle*>(fGenpi0Event->At(icon2)) ;

		// ABB add 7jan16 >>>>>>>> To file recondtructed Lambda from Sigma0 decay through rec. V0
		Int_t iV0 = LamRec->GetIsolation();
		AliESDv0 * v0 = fESDEvent->GetV0(iV0) ;
		AliESDtrack * pos = fESDEvent->GetTrack(v0->GetPindex()) ;
		AliESDtrack * neg = fESDEvent->GetTrack(v0->GetNindex()) ;
		TParticle * negMC = fStack->Particle(TMath::Abs(neg->GetLabel()));
		TParticle * posMC = fStack->Particle(TMath::Abs(pos->GetLabel()));
		if( !( negMC && posMC) ) continue;
		//  Int_t lblMotherPosV0Dghter = negMC->GetFirstMother() ;
		//  Int_t lblMotherNegV0Dghter = posMC->GetFirstMother();
	    
		if(! (( abs( negMC->GetMass() - 0.938272)< 0.01 && abs( posMC->GetMass() - 0.139570)<0.01   ) ||
		      ( abs( posMC->GetMass() - 0.938272)< 0.01 && abs( negMC->GetMass() - 0.139570)<0.01) ) ) continue;

		// ABB 8jan16 - may be use better selection, but not through the masses?

		TParticle * v0Lam = fStack->Particle(negMC->GetMother(0) ) ;
		//		      Int_t labelLp = v0Lam->GetMCLabelPositive() ;
		// Int_t labelLn = v0Lam->GetMCLabelNegative() ;	
		//  if( abs( labelLp - labelLn) != 1 ) continue;	
		FillHistogram("hCheckFB",11.5);

		// check that LamGen21 which get from Sigma0 is identical to v0Lam, may be improved, was 0.010 8jan15 
		if(! (  (fabs(LamGen21->Px()-v0Lam->Px())< 0.0010) &&  (fabs(LamGen21->Px()-v0Lam->Px())<0.0010) 
			&&   (fabs(LamGen21->Px()-v0Lam->Px())<0.0010) ) ) continue ;
		
		// if( !( labelLp21 ==  labelLp && labelLn21 == labelLn ) ) continue; 
		// final chack that Lambda came from the same Sigma0 as gamma
		// end ABB add 7jan16 <<<<<<<<<<
		      
		FillHistogram("hCheckFB",12.5);

		TLorentzVector sig0rec=photFB  + *LamRec;
		      
		//  Double_t Mlamrec = LamRec->M(); //unused variable
		Double_t corr_factor = 1. ;   // MassLam / Mlamrec;
		TLorentzVector ph2corr;
		ph2corr	= corr_factor * *LamRec;
		TLorentzVector sig0= photFB + ph2corr ;   // with corrected M_Lam
		if( ! (abs( sig0.M() - MassSig0 ) < SigmaMassWindow ) ) continue;
		      

		TLorentzVector sig0LrecGgen =  GGen22 + *LamRec ;
		TLorentzVector sig0LgenGrec = photFB  + LGen21  ;
		TLorentzVector sig0LgenGgen =  GGen22 + LGen21  ;
		TLorentzVector sig0LgenGgenM =  GGen22m + LGen21m  ;

		fLambdaTPx = LGen21.Px();
		fLambdaTPy = LGen21.Py();
		fLambdaTPz = LGen21.Pz();
		fGammaTPx = GGen22.Px();
		fGammaTPy = GGen22.Py();
		fGammaTPz = GGen22.Pz();
		fSigmaTPx = sig0LgenGgenM.Px();
		fSigmaTPy = sig0LgenGgenM.Py();
		fSigmaTPz = sig0LgenGgenM.Pz();


		FillHistogram("hMCgenrec3MvsPtLamgenGamgenM",sig0LgenGgenM.M(),sig0LgenGgenM.Pt() );
		FillHistogram("hMCgenrec3MvsPtLamgelimnGamgen",sig0LgenGgen.M(),sig0LgenGgen.Pt() );
		FillHistogram("hMCgenrec3MvsPtLamgenGamrec",sig0LgenGrec.M(),sig0LgenGrec.Pt() );
		FillHistogram("hMCgenrec3MvsPtLamrecGamgen",sig0LrecGgen.M(), sig0LrecGgen.Pt() );
		FillHistogram("hMCgenrec3MvsPtLamrecGamrec",sig0rec.M(), sig0rec.Pt() );
                                    
		FillHistogram("hMCgenrec3Sig0dMrecvsdPt",v0Sig2->GetMass()-sig0rec.M(),v0Sig2->Pt()-sig0rec.Pt());
		FillHistogram("hMCgenrec3Sig0mvsPt_uncorr",sig0rec.M(),sig0rec.Pt());
		      
		//!!! one more if needed - that Mass Sig0 is correct - +- 30 MeV around 3jul14
		printf("MC---DATA---Sig0-CC  M %f PT %f  E %f nLam %d \n", sig0.M(), sig0.Pt(), sig0.E(), nLam );
		      
		//		FillHistogram("hMCrec2Sig0PtLamdaVSGamma0",LamGen21->Pt(), GamGen22->Pt()  );
		//FillHistogram("hMCgenrec3Sig0RapPt", rap,  v0Sig2->Pt() );
		// FillHistogram("hMCgenrec3Sig0mvsPt_corr",sig0.M(),sig0.Pt());
		      
                                    
		//	  Double_t ModLam =  LamRec->GetXt(); unused variable
		TLorentzVector Pirec;
		Pirec.SetXYZM(LamRec->GetPpix(),LamRec->GetPpiy(),LamRec->GetPpiz(), MassPi) ;
		TLorentzVector Prrec;
		Prrec.SetXYZM(LamRec->GetPprx(),LamRec->GetPpry(),LamRec->GetPprz(), MassPr) ;
                                    
		//		  printf("Pirec %f %f %f Prrec %f %f %f \n",LamRec->GetPpix(),LamRec->GetPpiy(),LamRec->GetPpiz(),
		//	 LamRec->GetPprx(),LamRec->GetPpry(),LamRec->GetPprz() );
                                    
		Int_t NDt2 = LamGen21->GetNDaughters();
		if(  NDt2 != 2 ) continue;
		      
		Int_t DLam21  = LamGen21->GetFirstDaughter();
		Int_t DLam22  = LamGen21->GetLastDaughter();
		//		    printf("Pdg1 %d 2 %d \n", DLam21, DLam22 );
		      
                                    
		TParticle* PPrgen = (TParticle *)fStack->Particle( DLam21 );
		TParticle* PPigen = (TParticle *)fStack->Particle( DLam22 );
		
		//		  TParticle* LamGen21 = (TParticle *)fStack->Particle( Daughter21 );
		      
		TLorentzVector Pigen;		  TLorentzVector Prgen;
		Prgen.SetXYZM(PPrgen->Px(),PPrgen->Py(),PPrgen->Pz(), MassPr ) ;
		Pigen.SetXYZM(PPigen->Px(),PPigen->Py(),PPigen->Pz(), MassPi ) ;
		      
		TLorentzVector LamPrrecPigen =  Prrec + Pigen  ;
		TLorentzVector LamPrgenPirec =  Prgen + Pirec ;
		TLorentzVector LamPrrecPirec =  Prrec + Pirec ;
		TLorentzVector LamPrgenPigen =  Prgen + Pigen  ;
		      
                                    
		if( abs(rap) < fetaCut ) {
			
		  FillHistogram("hMCgenrec3Sig0Eta1",  v0Sig2->Pt() );
		  if(     LamGen21->GetPdgCode() == 3122 ) {
		    FillHistogram("hMCgenrec3PSig0Eta1",  v0Sig2->Pt() );
		    FillHistogram("hMCgenrec3PSig0Eta1Rap",  rap );
		  }
		  else if(LamGen21->GetPdgCode() == -3122) {
		    FillHistogram("hMCgenrec3ASig0Eta1",  v0Sig2->Pt() );
		    FillHistogram("hMCgenrec3ASig0Eta1Rap",  rap );
		  }
                                        
		  char key[55] ;
		  for (Int_t irap = 0; irap < 6; irap++) {
		    if( abs(rap) < fEtaCuts[irap] ) {
		      
		      if(   LamGen21->GetPdgCode() == 3122 ) {
			sprintf(key,"hMCPSig0PtRec%d",irap) ;
			FillHistogram(key,  v0Sig2->Pt() );
			//		      if (v0Sig2->IsPrimary() ) {
			//		sprintf(key,"hMCPSig0PrimPtRec%d",irap) ;
			//	FillHistogram(key,  v0Sig2->Pt() );
			//}
		      }
		      else  if(   LamGen21->GetPdgCode() == -3122 ) {
			sprintf(key,"hMCASig0PtRec%d",irap) ;
			FillHistogram(key,  v0Sig2->Pt() );
			//if (v0Sig2->IsPrimary() ) {
			//	sprintf(key,"hMCASig0PrimPtRec%d",irap) ;
			//	FillHistogram(key,  v0Sig2->Pt() );
			//}
		      }
		    }
		  }
                        
		  FillHistogram("hRec3GamSig0Eta1",  GamGen->Pt() );
		  FillHistogram("hRec4GamSig0Eta1",   GamGen22->Pt() );
		  FillHistogram("hRec3LamSig0Eta1",  LamGen21->Pt() );
                        
		  FillHistogram("hMCgenrec3Sig0dMvsdPt",v0Sig2->GetMass()- sig0.M(),v0Sig2->Pt()-sig0.Pt());
		  if(LamGen21->GetPdgCode()==3122) FillHistogram("hMCgenrec3PSig0dMvsdPt",
								 v0Sig2->GetMass()- sig0.M(),v0Sig2->Pt()-sig0.Pt());
			else if(LamGen21->GetPdgCode()==-3122) FillHistogram("hMCgenrec3ASig0dMvsdPt",
									     v0Sig2->GetMass()-sig0.M(),v0Sig2->Pt()-sig0.Pt());
			
		  //printf("MC-RUN %d Ev %d pdf D1 %d D2 %d \n",
		  //fESDEvent->GetRunNumber(),fESDEvent->GetEventNumberInFile(),LamGen21->GetPdgCode(),GamGen22->GetPdgCode() );
		  //		  nLam = nLam/LamRecFl ;
                        
		  if(  abs(  LamGen21->Pt() -   LamRec->Pt() ) < 0.01   ) {
		    // or LamRecFl // Have to check that PROPER Lambda is reconstructed!!!
		    FillHistogram("hMCgenrec3Sig0",  v0Sig2->Pt() );
		    
		    //   Double_t rap = Rapidity(  v0Sig2->Pt(),  v0Sig2->Pz(),  v0Sig2->GetMass() );
		    FillHistogram("hMCgenrec3Sig0Rap", rap );
		    // FillHistogram("hMCrec3Sig0PtLamdaVSGamma0",LamGen21->Pt(), LamRec->Pt()  );
		    FillHistogram("hRec3GamSig0",  GamGen->Pt() );
		    FillHistogram("hRec4GamSig0",   GamGen22->Pt() );
		    FillHistogram("hRec3LamSig0",  LamGen21->Pt() );
		  }
		}
	      }
	    } // if abs(LamGen21==Lambda && GamGen22=gamma)
	  }
	} // end of mass == Sig0
      }
    }   // end of fStack - MC check
    
  }   // end of the loop of all conversion photons        
}
//_______________________________________________________________________________
void AliAnalysisTaskSigma0::SelectConvPhotons(){
    //Fill list of conversion photons, that is scan v0s and select photon-like
    
    Double_t MassPi = 0.13957018 ;
    Double_t MassPr = 0.938272046 ;
    Double_t MassSig0 = 1.192642 ;
   
    //set some constants
    const Double_t kCutSigmaMass=0.0001;  //Constraint on photon mass
    const Bool_t useImprovedVertex=kTRUE ; //Use verted with converted photon?
    //  const Double_t zrSlope = TMath::Tan(2*TMath::ATan(TMath::Exp(-fetaCut)));
    //  const Double_t zrSlope12 = TMath::Tan(2*TMath::ATan(TMath::Exp(-1.2)));
    const Double_t kzrSlope09 = TMath::Tan(2*TMath::ATan(TMath::Exp(-0.9)));
    const Double_t kzOffset = 7.;
    //upto v25  const Double_t kminR = 5.;
    const Double_t kminR = 0.01;
    // Double_t MassLam = 1.115683 ; //unused variable
   
    //  Bool_t fToUseCF(kFALSE);
    //  Bool_t fConvCFCont(kFALSE);
    
    //No primary vertex in event => scip
    if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0)   return;
    
    Double_t  lMagneticField = fESDEvent->GetMagneticField();
    
    const AliESDVertex *esdVertex5 = fESDEvent->GetPrimaryVertex();
    Double_t vtx5[3];
    vtx5[0] = esdVertex5->GetX();
    vtx5[1] = esdVertex5->GetY();
    vtx5[2] = esdVertex5->GetZ();
    Double_t  lPrimaryVtxPosition[3];
    lPrimaryVtxPosition[0] = vtx5[0]; lPrimaryVtxPosition[1] = vtx5[1]; lPrimaryVtxPosition[2] = vtx5[2];
    
    Int_t nV0=fESDEvent->GetNumberOfV0s() ;
    Int_t inConv=0 ;
    
    Double_t minPnSigmaAbovePionLine = 1. ;
    Double_t maxPnSigmaAbovePionLine = 3. ;
    Double_t nSigmaAbovePionLine = 0 ;
    //    Int_t NlamEv = 0;   Int_t NalamEv = 0;  Int_t NlamminalamEv = 0;
    
    for(Int_t iv0=0; iv0<nV0;iv0++){
        AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;
        
        AliESDtrack * pos = fESDEvent->GetTrack(v0->GetPindex()) ;
        AliESDtrack * neg = fESDEvent->GetTrack(v0->GetNindex()) ;
        const AliExternalTrackParam * paramPos = v0->GetParamP() ;
        const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
        if(pos->GetSign() <0){//change tracks
            pos=neg ;
            neg=fESDEvent->GetTrack(v0->GetPindex()) ;
            paramPos=paramNeg ;
            paramNeg=v0->GetParamP() ;
        }
        
        //cuts ------------    //remove like sign pairs
        if(pos->GetSign() == neg->GetSign())  continue ;
        
        if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) continue ;
        
        if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) || !(neg->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
        
        // Quality tracks cuts: fESDtrackCuts
        //     if ( !(fESDtrackCuts->IsSelected(pos)) || !(fESDtrackCuts->IsSelected(neg)) ) continue;
        
        Int_t    lOnFlyStatus        = 0;
        Bool_t    isOnFly        = 0;
        // VO's main characteristics to check the reconstruction cuts
        lOnFlyStatus       = v0->GetOnFlyStatus();
        isOnFly = v0->GetOnFlyStatus();
        if( lOnFlyStatus == 0 ) continue;
        
        //  lChi2V0            = v0->GetChi2V0();
        //  Double_t  lDcaV0Daughters    = v0->GetDcaV0Daughters();
        //  lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
        Double_t  lV0cosPointAngle   = v0->GetV0CosineOfPointingAngle(vtx5[0],vtx5[1], vtx5[2]);
        
        Double_t  lV0Position[3];
        v0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);
        
        Double_t  lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
        //      lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
        //		                   TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
        //		                   TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));
        //  Double_t lV0tDecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
        //				   TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2));
        // DCA between daughter and Primary Vertex:
        Double_t  lDcaPosToPrimVertex = TMath::Abs(pos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
        
        Double_t lDcaNegToPrimVertex = TMath::Abs(pos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
        
        if ((lDcaPosToPrimVertex < 0.1) || (lDcaNegToPrimVertex < 0.1) ||
            // (lDcaV0Daughters > 1.00) ||
            (lV0cosPointAngle < 0.998) || (lV0Radius < 0.0) || (lV0Radius > 100. ) ) continue;  // (lV0Radius < 0.5) was before 1.10.13
        
        Double_t v0x=0.,v0y=0.,v0z=0.;
        v0->GetXYZ(v0x,v0y,v0z) ;
        Double_t r=TMath::Sqrt(v0x*v0x + v0y*v0y) ;
        
        AliKFParticle negKF(*paramNeg,11);
        AliKFParticle posKF(*paramPos,-11);
        AliKFParticle photKF(negKF,posKF) ;
        // printf("st 1: px=%f, py=%f, pz=%f, E=%f \n",photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;
        photKF.SetMassConstraint(0,kCutSigmaMass);
        // printf("st 2: px=%f, py=%f, pz=%f, E=%f \n",photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;
        
        if(useImprovedVertex){
            AliKFVertex primaryVertexImproved(*(fESDEvent->GetPrimaryVertex()));
            primaryVertexImproved+=photKF;
            photKF.SetProductionVertex(primaryVertexImproved);
        }
        //     printf("st 3: px=%f, py=%f, pz=%f, E=%f \n",photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;
        Double_t m, width ;
        photKF.GetMass(m,width);
        
    //    FillHistogram("TR-0-mass",m) ;  // Not official
    //    FillHistogram("TR-0-width",width) ; // Not official
        
        
        TLorentzVector GamRec;
        //	AliCaloParticle GamRec;
        GamRec.SetXYZM(photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),0.) ;  //Produces slightly better pi0 width
        
        
        if ( ! v0->GetOnFlyStatus() ) continue;
        
        if(neg->GetNcls(1) <2 || pos->GetNcls(1) <2)     continue ;
        
        if(pos->GetSign() == neg->GetSign())      continue ;
        
    //    FillHistogram("TR-mass",m) ; // Not official
    //    FillHistogram("TR-width",width) ; // Not official
    //    FillHistogram("TR-radius",r) ; // Not official
        
        //  cut-4 TPC track quality cut - both tracks are refitted
        if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) ||
           !(neg->GetStatus() & AliESDtrack::kTPCrefit) )  continue;
        Bool_t isKink=kFALSE ;
        if( neg->GetKinkIndex(0) > 0 ||
           pos->GetKinkIndex(0) > 0) {
            isKink=kTRUE;
        }
        if( isKink ) continue;
        
        //  cut-6 First rough PID for e+e- pairs (epm) and \pi+ \pi- pairs (pipm), will be used later
        // ABB: Do we need it?
        Bool_t isepmpid=kTRUE ;
        Bool_t ispipmpid=kTRUE ;
        if( fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)<-4. ||
           fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)>5. ||
           fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)<-4. ||
           fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)>5. ){
            isepmpid=kFALSE ;
        }
        if(  fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)<-4. ||
           fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)>5. ||
           fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)<-4. ||
           fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)>5. ){
            ispipmpid=kFALSE ;
        }
        
        minPnSigmaAbovePionLine = 0.5 ;
        maxPnSigmaAbovePionLine = 100. ;
        nSigmaAbovePionLine = 0 ;
        
        Bool_t isepmdEdx=kTRUE ;
        if(pos->P()>minPnSigmaAbovePionLine && pos->P()<maxPnSigmaAbovePionLine ){
            if(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)<nSigmaAbovePionLine){
                isepmdEdx=kFALSE;
            }
        }
        if(neg->P()>minPnSigmaAbovePionLine && neg->P()<maxPnSigmaAbovePionLine){
            if(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)<nSigmaAbovePionLine){
                isepmdEdx=kFALSE;
            }
        }
        
        const Double_t kSigmaAroundLine=1. ;
        Bool_t isdEdx=isepmdEdx ;
        
        const Double_t kMinPPionRejection=0.5 ;
        if(neg->P()<kMinPPionRejection ){
            if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion))<kSigmaAroundLine){
                //printf("... dEdx 8 \n") ;
                isdEdx=kFALSE;
            }
        }
        if(pos->P()<kMinPPionRejection ){
            if( TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion))<kSigmaAroundLine){
                //printf("... dEdx 9 \n") ;
                isdEdx=kFALSE;
            }
        }
        Bool_t isRout=kTRUE ;
        Bool_t isRin=kFALSE ;
        if(   r>kminR && r<fmaxR  ){ // cuts on distance from collision point
            isRout = kTRUE;
        }
        if ( r < kminR )  isRin=kTRUE ;
        
        Bool_t ismaxZgam=kTRUE ;
        Bool_t isminZpi=kFALSE ;
        Double_t fminZ = 5.;
        if(TMath::Abs(v0z) > fmaxZ ){ // cuts out regions where we do not reconstruct
            ismaxZgam=kFALSE ;
        }
        if ( TMath::Abs(v0z) < fminZ ){
            isminZpi=kTRUE ;
        }
        
        // open
        Bool_t isStrictZ=kFALSE ;
        if((TMath::Abs(v0z)*kzrSlope09)-kzOffset < r )
        isStrictZ=kTRUE ;
        
        if(photKF.GetNDF()<=0)   continue;
        
        Double_t chi2V0 = photKF.GetChi2()/photKF.GetNDF();
        if(chi2V0 > fchi2CutConversion || chi2V0 <=0)      continue;
        Bool_t isStrictChi=kFALSE ;
        if(chi2V0 < 0.7*fchi2CutConversion && chi2V0 >0){
            isStrictChi=kTRUE;
        }
        
        const Double_t kWideEtaCut=9.9 ;
        if(TMath::Abs(GamRec.Eta())> kWideEtaCut){
            // printf("... ETA \n") ;
            continue;
        }
        if(TMath::Abs(paramPos->Eta())> kWideEtaCut ||
           TMath::Abs(paramNeg->Eta())> kWideEtaCut ){
            //printf("... ETA pls mns \n") ;
            continue ;
        }
        
        // cut-12 use 0.050 as CC
        // assume that it is relevant for both e+e- and \pi+\pi- pairs
        if(GamRec.Pt()<fptCut){
            // printf("... pt \n") ;
            continue;
        }
        if( GamRec.Pt() > fptMaxCut)  continue;
        
        Double_t w=PlanarityAngle(paramPos,paramNeg) ;
        FillHistogram("hAllMvsWidth",w,m) ;
        
        //Single photon spectrum
        //Double_t pt=GamRec.Pt() ;  // Not official
        if(isOnFly){
	  Double_t posp[3]= { posKF.Px(),  posKF.Py(),  posKF.Pz() };
	  Double_t negp[3]= { negKF.Px(),  negKF.Py(),  negKF.Pz() };
	  Double_t moth[3]= { GamRec.Px(), GamRec.Py(), GamRec.Pz() };
	  Double_t arpod[2]= {0,0};
	  GetArPod( posp, negp, moth, arpod );
            
	  if( !( abs(  arpod[1] ) < 0.8 && arpod[0] < 0.06 ) ) continue;
	  
	  //	  FillHistogram("hArPodGConv",  arpod[1], arpod[0] ) ;
	  // FillHistogram("hCGamPtConv", GamRec.Pt() ) ;
            
	  // printf(" CONVERTED PHOTON RUN %d Ev %d \n",  fESDEvent->GetRunNumber(),  fESDEvent->GetEventNumberInFile() );
	  // Fill TLIST for Sigma0
	  if(inConv>=fConvEvent->GetSize())fConvEvent->Expand(inConv*2) ;
	  AliCaloParticle *phg = new((*fConvEvent)[inConv]) AliCaloParticle(GamRec) ;
	  phg->SetV0CosPointingAngle(lV0cosPointAngle);
	  phg->SetV0Radius(lV0Radius);
	  phg->SetArmPt(arpod[0]);
	  phg->SetArmAlpha(arpod[1]);
	  phg->SetV0Chi2(chi2V0);
	  phg->SetZConv(v0z);
            
	  inConv++ ;

	  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  //Fill MC information for GammaConverted
	  if(fStack){
	    //abb    if ( 1>0 ) {
	    TParticle * negativeMC = fStack->Particle(TMath::Abs(neg->GetLabel()));
	    TParticle * positiveMC = fStack->Particle(TMath::Abs(pos->GetLabel()));
	    
	    if(negativeMC && positiveMC){
	      if(negativeMC->GetMother(0) != positiveMC->GetMother(0))
		continue ;
	    }
	    
	    if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
	      continue;
	    }
	    if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
	      continue;
	    }
                
	    TParticle * GamGen = fStack->Particle(negativeMC->GetMother(0));  // Generatred Gamma
                
	    //	    FillHistogram("hMCV0MvsWidth",width,m) ;
                
	    if(negativeMC->GetUniqueID() != 5 || positiveMC->GetUniqueID() !=5){ // id 5 is conversion
	      continue; }
	   

 
	    //	    FillHistogram("hMCV0MvsWidth1",width,m) ;

	    if(GamGen->GetPdgCode() == 22){    // indeed this Gamma reconstracted from Gamma
	      //     FillHistogram("hMCV0MvsWidth2",width,m) ;
	      //!!! add histo	      
	      /* Yura 12.3.14 primary Gamma has no Mother ! */
	     
	      if( GamGen->IsPrimary() ) continue;
	      TParticle * v0Sig2 = fStack->Particle(GamGen->GetMother(0));
	      
	      // 	      FillHistogram("hMCV0MvsWidth3",width,m) ;

	      //       	printf("=========================mSig %f \n",  v0Sig2->GetMass() );
	      if( abs(v0Sig2->GetMass() - 1.192642)<0.05 ) {

		// FillHistogram("hMCV0MvsWidth4",width,m) ;		

		Int_t NDt1 =v0Sig2->GetNDaughters();                
	       	if( abs( v0Sig2->GetPdgCode() ) == 3212  && NDt1 == 2)  {    // Sigma0
		//		if ( 1>0 ) {

		  Int_t Daughter21  = v0Sig2->GetFirstDaughter();
		  Int_t Daughter22  = v0Sig2->GetLastDaughter();
		  Double_t mrec=0.,  mcalc=0. ;
		  mrec =     v0Sig2->GetMass();
		  mcalc = v0Sig2->GetCalcMass();
		  //	      	    printf("---- m %f w %f cm %f  \n", m, width, mcalc  );
		  FillHistogram("hMCgenrec2Sig0",  v0Sig2->Pt() );
		  Double_t rap = Rapidity(  v0Sig2->Pt(),  v0Sig2->Pz(),  v0Sig2->GetMass() );
		  FillHistogram("hMCgenrec2Sig0Rap", rap );
		  FillHistogram("hRecGamSig0",  GamGen->Pt() );
		  
		  FillHistogram("hMCgenrec2Sig0RapPt", rap,  v0Sig2->Pt() );
		  if( abs(rap) < fetaCut ) {
		    FillHistogram("hMCgenrec2Sig0Eta1",  v0Sig2->Pt() );
		    FillHistogram("hRecGamSig0Eta1",  GamGen->Pt() );
		  }
		  
		  TParticle* LamGen21 = (TParticle *)fStack->Particle( Daughter21 );
		  TLorentzVector LGen21;
		  LGen21.SetXYZM(LamGen21->Px(),LamGen21->Py(),LamGen21->Pz(),LamGen21->GetCalcMass()) ;  //Produces slightly better pi0 width
		  TLorentzVector LGen21m;
		  LGen21m.SetXYZM(LamGen21->Px(),LamGen21->Py(),LamGen21->Pz(),LamGen21->GetMass()) ;  //Produces slightly better pi0 width
		  
                            
		  TParticle* GamGen22 = (TParticle *)fStack->Particle( Daughter22 );
		  TLorentzVector GGen22;
		  GGen22.SetXYZM(GamGen22->Px(),GamGen22->Py(),GamGen22->Pz(),GamGen22->GetCalcMass()) ;  //Produces slightly better pi0 width
		  TLorentzVector GGen22m;
		  GGen22m.SetXYZM(GamGen22->Px(),GamGen22->Py(),GamGen22->Pz(),GamGen22->GetMass()) ;  //Produces slightly better pi0 width
                            
		  //	  if ( (abs( LamGen21->GetPdgCode() ) == 3122 && abs( GamGen22->GetPdgCode())  == 22)
		       //old-skip		 &&  (inLam >0 && LamRecFl > 0)
		  //     ) {
		    
		  if ( 1>0 ) {
		    // Search here Sig0 - Loop over all Lambda reconstr in the event
		    Int_t nLam=fGenpi0Event->GetEntriesFast() ;
		    for(Int_t icon2=0; icon2<nLam; icon2++){
		      AliCaloParticle * LamRec= static_cast<AliCaloParticle*>(fGenpi0Event->At(icon2)) ;
		      
		      TLorentzVector sig0rec= GamRec + *LamRec;
		      
		      //  Double_t Mlamrec = LamRec->M(); //unused variable
		      Double_t corr_factor = 1. ;    // MassLam / Mlamrec;
		      TLorentzVector ph2corr;
		      ph2corr	= corr_factor * *LamRec;
		      TLorentzVector sig0= GamRec + ph2corr ;   // with corrected M_Lam
		      if( ! (abs( sig0.M() - MassSig0 ) < 0.025 ) ) continue;
		      
		      TLorentzVector sig0LrecGgen =  GGen22 + *LamRec ;
		      TLorentzVector sig0LgenGrec =  GamRec + LGen21  ;
		      TLorentzVector sig0LgenGgen =  GGen22 + LGen21  ;
		      TLorentzVector sig0LgenGgenM =  GGen22m + LGen21m  ;
		      FillHistogram("hMCgenrec3MvsPtLamgenGamgenM",sig0LgenGgenM.M(),sig0LgenGgenM.Pt() );
		      FillHistogram("hMCgenrec3MvsPtLamgenGamgen",sig0LgenGgen.M(),sig0LgenGgen.Pt() );
		      FillHistogram("hMCgenrec3MvsPtLamgenGamrec",sig0LgenGrec.M(),sig0LgenGrec.Pt() );
		      FillHistogram("hMCgenrec3MvsPtLamrecGamgen",sig0LrecGgen.M(), sig0LrecGgen.Pt() );
		      FillHistogram("hMCgenrec3MvsPtLamrecGamrec",sig0rec.M(), sig0rec.Pt() );
                                    
		      FillHistogram("hMCgenrec3Sig0dMrecvsdPt",v0Sig2->GetMass()-sig0rec.M(),v0Sig2->Pt()-sig0rec.Pt());
		      FillHistogram("hMCgenrec3Sig0mvsPt_uncorr",sig0rec.M(),sig0rec.Pt());
		      
		      //!!! one more if needed - that Mass Sig0 is correct - +- 30 MeV around 3jul14
		      // printf("MC---DATA---Sig0-CC  M %f PT %f  E %f nLam %d \n", sig0.M(), sig0.Pt(), sig0.E(), nLam );
		      
		      FillHistogram("hMCrec2Sig0PtLamdaVSGamma0",LamGen21->Pt(), GamGen22->Pt()  );
		      FillHistogram("hMCgenrec3Sig0RapPt", rap,  v0Sig2->Pt() );
		      FillHistogram("hMCgenrec3Sig0mvsPt_corr",sig0.M(),sig0.Pt());
		      
                                    
		      //	  Double_t ModLam =  LamRec->GetXt(); unused variable
                                    
		      TLorentzVector Pirec;
		      Pirec.SetXYZM(LamRec->GetPpix(),LamRec->GetPpiy(),LamRec->GetPpiz(), MassPi) ;
		      TLorentzVector Prrec;
		      Prrec.SetXYZM(LamRec->GetPprx(),LamRec->GetPpry(),LamRec->GetPprz(), MassPr) ;
		      
                                    
		      //		  printf("Pirec %f %f %f Prrec %f %f %f \n",LamRec->GetPpix(),LamRec->GetPpiy(),LamRec->GetPpiz(),
		      //	 LamRec->GetPprx(),LamRec->GetPpry(),LamRec->GetPprz() );
		      
                                    
		      Int_t NDt2 = LamGen21->GetNDaughters();
		      if(  NDt2 != 2 ) continue;
		      
		      Int_t DLam21  = LamGen21->GetFirstDaughter();
		      Int_t DLam22  = LamGen21->GetLastDaughter();
		      //		    printf("Pdg1 %d 2 %d \n", DLam21, DLam22 );
		      
                                    
		      TParticle* PPrgen = (TParticle *)fStack->Particle( DLam21 );
		      TParticle* PPigen = (TParticle *)fStack->Particle( DLam22 );
		      
		      //		  TParticle* LamGen21 = (TParticle *)fStack->Particle( Daughter21 );
		      
		      TLorentzVector Pigen;		  TLorentzVector Prgen;
		      Prgen.SetXYZM(PPrgen->Px(),PPrgen->Py(),PPrgen->Pz(), MassPr ) ;
		      Pigen.SetXYZM(PPigen->Px(),PPigen->Py(),PPigen->Pz(), MassPi ) ;
		      
		      TLorentzVector LamPrrecPigen =  Prrec + Pigen  ;
		      TLorentzVector LamPrgenPirec =  Prgen + Pirec ;
		      TLorentzVector LamPrrecPirec =  Prrec + Pirec ;
		      TLorentzVector LamPrgenPigen =  Prgen + Pigen  ;
		      
                                    
		      if( abs(rap) < fetaCut ) {
			
			FillHistogram("hMCgenrec3Sig0Eta1",  v0Sig2->Pt() );
			if(     LamGen21->GetPdgCode() == 3122 ) {
			  FillHistogram("hMCgenrec3PSig0Eta1",  v0Sig2->Pt() );
			  FillHistogram("hMCgenrec3PSig0Eta1Rap",  rap );
			}
			else if(LamGen21->GetPdgCode() == -3122) {
			  FillHistogram("hMCgenrec3ASig0Eta1",  v0Sig2->Pt() );
			  FillHistogram("hMCgenrec3ASig0Eta1Rap",  rap );
			}
                                        
			char key[55] ;
			for (Int_t irap = 0; irap < 6; irap++) {
			  if( abs(rap) < fEtaCuts[irap] ) {
			    
			    if(   LamGen21->GetPdgCode() == 3122 ) {
			      sprintf(key,"hMCPSig0PtRec%d",irap) ;
			      FillHistogram(key,  v0Sig2->Pt() );
			      //		      if (v0Sig2->IsPrimary() ) {
			      //		sprintf(key,"hMCPSig0PrimPtRec%d",irap) ;
			      //	FillHistogram(key,  v0Sig2->Pt() );
			      //}
			    }
			    else  if(   LamGen21->GetPdgCode() == -3122 ) {
			      sprintf(key,"hMCASig0PtRec%d",irap) ;
			      FillHistogram(key,  v0Sig2->Pt() );
			      //if (v0Sig2->IsPrimary() ) {
			      //	sprintf(key,"hMCASig0PrimPtRec%d",irap) ;
			      //	FillHistogram(key,  v0Sig2->Pt() );
			      //}
			    }
			  }
			}
                        
			FillHistogram("hRec3GamSig0Eta1",  GamGen->Pt() );
			FillHistogram("hRec4GamSig0Eta1",   GamGen22->Pt() );
			FillHistogram("hRec3LamSig0Eta1",  LamGen21->Pt() );
                        
			FillHistogram("hMCgenrec3Sig0dMvsdPt",v0Sig2->GetMass()- sig0.M(),v0Sig2->Pt()-sig0.Pt());
			if(LamGen21->GetPdgCode()==3122) FillHistogram("hMCgenrec3PSig0dMvsdPt",
								       v0Sig2->GetMass()- sig0.M(),v0Sig2->Pt()-sig0.Pt());
			else if(LamGen21->GetPdgCode()==-3122) FillHistogram("hMCgenrec3ASig0dMvsdPt",
									     v0Sig2->GetMass()-sig0.M(),v0Sig2->Pt()-sig0.Pt());
			
			//printf("MC-RUN %d Ev %d pdf D1 %d D2 %d \n",
			//fESDEvent->GetRunNumber(),fESDEvent->GetEventNumberInFile(),LamGen21->GetPdgCode(),GamGen22->GetPdgCode() );
			//		  nLam = nLam/LamRecFl ;
                        
			if(  abs(  LamGen21->Pt() -   LamRec->Pt() ) < 0.01   ) {
			  // or LamRecFl // Have to check that PROPER Lambda is reconstructed!!!
			  FillHistogram("hMCgenrec3Sig0",  v0Sig2->Pt() );
			  
			  //   Double_t rap = Rapidity(  v0Sig2->Pt(),  v0Sig2->Pz(),  v0Sig2->GetMass() );
			  FillHistogram("hMCgenrec3Sig0Rap", rap );
			  FillHistogram("hMCrec3Sig0PtLamdaVSGamma0",LamGen21->Pt(), LamRec->Pt()  );
			  FillHistogram("hRec3GamSig0",  GamGen->Pt() );
			  FillHistogram("hRec4GamSig0",   GamGen22->Pt() );
			  FillHistogram("hRec3LamSig0",  LamGen21->Pt() );
			}
		      }
		    }
		  } // if abs(LamGen21==Lambda && GamGen22=gamma)
		}
	      } // end of mass == Sig0
	    }
	  }   // end of fStack - MC check
        }    // end of if
    }// end of V0 loop
    
} // END of SelectConvPhoton

//_____________________________________________________
//___________________________________________________________________________
void AliAnalysisTaskSigma0::SelectTracks(){
    //Select pi+- tracks for search of Sigma+- 1385, PID have to be applied later

    // Int_t inTracks=0;
    Int_t inPHOS=0;
    Double_t mpi = 0.13957018;
    
    Int_t nTrk=fESDEvent->GetNumberOfTracks();
    for (Int_t iTracks = 0; iTracks < nTrk; iTracks++) {
        AliESDtrack* track = (AliESDtrack*)fESDEvent->GetTrack(iTracks);
        
        if(!track)      continue;
        if(  track->Pt() > 1.5 ) continue;
        
        Double_t trkTOF = 0;
        trkTOF = track->GetTOFsignal();
        
        //    Float_t trkRxy = 0;    Float_t trkRz = 0;
        //    track->GetImpactParameters(trkRxy, trkRz );
        //    trkTOF = track->GetDCA();
        
        //    Float_t tpcRxy = 0;    Float_t tpcRz = 0;
        //    track->GetImpactParameters(tpcRxy, tpcRz );
        //    trkTOF = track->GetDCA();
        
        Double_t trkM = 0;
        trkM = track->GetMass();
        
        if( !(trkM > 0.12 && trkM < 0.142) ) continue;
        //    if( trkM < 0.9 ) continue;
        
        if(  track->GetTPCsignal() < 30 ) continue;
        //    Double_t trkPID = 0;
        // trkPID = track->GetESDpid();
        
        //    Double_t r[10] = {0.};     track->GetESDpid(r);   trkPID = r[3];
        
        
        // to account the mass po pi+-
        TLorentzVector pi0rec;
        pi0rec.SetXYZM( track->Px(), track->Py(), track->Pz(), mpi ) ;  //Produces slightly better pi0 width
        
	//        FillHistogram("hdEdxTrack", track->GetP(), track->GetTPCsignal()) ;
        
        //    if(inTracks >= fTrackEvent->GetSize()) fTrackEvent->Expand(inTracks*2) ;
        //We assume here that track is massless
        //  new((*fTrackEvent)[inTracks]) AliCaloParticle(track->Px(),track->Py(),track->Pz(),track->P()) ;
        // inTracks++ ;
        
        
        if(inPHOS>=fPHOSEvent->GetSize())fPHOSEvent->Expand(2*inPHOS) ;
        //   AliCaloParticle * ph = new((*fPHOSEvent)[inPHOS]) AliCaloParticle(pi0rec); //unused variable
        //    AliCaloParticle * ph = new((*fPHOSEvent)[inPHOS]) AliCaloParticle(track->Px(),track->Py(),track->Pz(),track->P() ) ;
        inPHOS++ ;
    }  // end loop over tracks
} // end selectTracks


//_________________________________________________________________________
//___________________________________________________________________________
void AliAnalysisTaskSigma0::ProcessMC(){
    
    //fill histograms for efficiensy etc. calculation
    
    if(!fStack) return ;
    //    printf("ProcessMC \n");
    
    const Double_t kRcut = 1. ; //cut for primary particles
    Double_t vtx[3];
    vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
    vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
    vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();
    
    Int_t Daughter1  = 0;  Int_t Daughter2  = 0;
    
    //    if(  fStack->GetNtrack() < 2 ) return;
    
    //---------First Lambda -----------------------------
    for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
        TParticle* particle = (TParticle *)fStack->Particle(iTracks);
        
        //    Int_t mPdgSign   = particle->GetPdgCode();
        //    Double_t charge  = particle->GetPDG()->Charge();
        //     Int_t mPdg       = TMath::Abs(mPdgSign);
        //      Int_t mStatus    = particle->GetStatusCode() ;
        //  Int_t iParent    = particle->GetFirstMother() ;
        //  Int_t iParent2   = particle->GetSecondMother() ;
        Int_t NDt = particle->GetNDaughters();
        if( NDt == 2 ){
            Daughter1  = particle->GetFirstDaughter();
            Daughter2  = particle->GetLastDaughter();
        }
        // Double_t mEta    = particle->Eta();    // Double_t mPhi    = particle->Phi();    //  Double_t mPt = particle->Pt();
        //  Double_t mE      = particle->Energy();    //     Double_t mM      = particle->GetMass();  // Double_t mP = particle->P();
        
        // if(  abs( particle->GetPdgCode()) == 3324 && NDt == 2)  {
	//        printf("MC----itr %d status %d PID=%d Parent %d Daut1=%d Daut2=%d ", iTracks,mStatus,mPdgSign,iParent,Daughter1,Daughter2);
        // printf("MC E %f Pt=%f M=%f  Pdg %d \n", mE,  mPt, mM, mPdg );
	//     printf("MC  itr %d \n", iTracks );     //    }
        
	if(  abs( particle->GetPdgCode() ) == 3122  && NDt == 2)  { 

	  //	  printf("   all Lambda generated \n" );
            
            Float_t ptLam =  particle->Pt();
            Float_t pdgLam =  particle->GetPdgCode();
            
            FillHistogram("hMCgenLam",  particle->Pt() );
            if(  particle->GetPdgCode() == 3122 )    FillHistogram("hMCgenPLam",  particle->Pt() );
            else if( particle->GetPdgCode() == -3122 )    FillHistogram("hMCgenALam",  particle->Pt() );
            
            Float_t rapLam = Rapidity(particle->Pt(), particle->Pz(), particle->GetMass());
            FillHistogram("hMCgenLamRap",  rapLam );
            if(  particle->GetPdgCode() == 3122 )    FillHistogram("hMCgenPLamRap",  rapLam );
            else if( particle->GetPdgCode() == -3122 )    FillHistogram("hMCgenALamRap",  rapLam );
            
            char key[55] ;
            for (Int_t irap = 0; irap < 6; irap++) {
                if( abs(rapLam) < fEtaCuts[irap] ) {
                    if( pdgLam>0 ) {
                        sprintf(key,"hMCPLamPtGen%d",irap) ;
                        FillHistogram(key,  ptLam );
			//                       if( particle->IsPrimary() ) {
			//   sprintf(key,"hMCPLamPrimPtGen%d",irap) ;
                        //    FillHistogram(key,  ptLam );
                        // }
                    }
                    else  if( pdgLam<0 ) {
                        sprintf(key,"hMCALamPtGen%d",irap) ;
                        FillHistogram(key,  ptLam );
			//  if( particle->IsPrimary() ) {
                        //    sprintf(key,"hMCALamPrimPtGen%d",irap) ;
                        //    FillHistogram(key,  ptLam );
			// }
                    }
                }
            }

	    // 25sep16-abb            
	    //            if( abs( rapLam < fetaCut ) ) {
	    if( abs( rapLam ) < fetaCut  ) {
                FillHistogram("hMCgenLamEta1",  particle->Pt() );
                if(  particle->GetPdgCode() == 3122 )    FillHistogram("hMCgenPLamEta1",  particle->Pt() );
                else if( particle->GetPdgCode() == -3122 )    FillHistogram("hMCgenALamEta1",  particle->Pt() );
            }
        }
	
        
        // Sigma0 generated
        if( abs( particle->GetPdgCode() ) == 3212  && NDt == 2)  {    

	  //	  printf(" Sigma0 generated \n");
     

	  Double_t m=0., mcalc=0. ;
            m =      particle->GetMass();
            mcalc = particle->GetCalcMass();
            Float_t rap = Rapidity(particle->Pt(), particle->Pz(), particle->GetMass());
            Float_t ptSig0 = particle->Pt();
            
            if( ptSig0 < 1 ) continue; // can reconstruct Sigma0 with pT>1 GeV/c ONLY!
            
            Float_t pdgSig0 = particle->GetPdgCode();
            
	    //	               printf(" m %f m-calc %f PDG %d \n", m,  mcalc, particle->GetPdgCode() == 3212 );
            FillHistogram("hMCgenSig0",  particle->Pt() );
            FillHistogram("hMCgenSig0Rap",  rap );
	    //            FillHistogram("hMCgenSig0RapEta",  rap,  particle->Eta() );
            // FillHistogram("hMCgenSig0RapPt", rap,  particle->Pt() );
            
            char key[55] ;
            for (Int_t irap = 0; irap < 6; irap++) {
                if( abs(rap) < fEtaCuts[irap] ) {
                    if( pdgSig0>0 ) {
                        sprintf(key,"hMCPSig0PtGen%d",irap) ;
                        FillHistogram(key,  ptSig0 );
                        
			//                        if( particle->IsPrimary() ) {
                        //    sprintf(key,"hMCPSig0PrimPtGen%d",irap) ;
                        //    FillHistogram(key,  ptSig0 );
			// }
                    }
                    else  if( pdgSig0<0 ) {
                        sprintf(key,"hMCASig0PtGen%d",irap) ;
                        FillHistogram(key,  ptSig0 );
			// if( particle->IsPrimary() ) {
                        //    sprintf(key,"hMCASig0PrimPtGen%d",irap) ;
                        //    FillHistogram(key,  ptSig0 );
                        // }
                    }
                }
            }
            //      if( abs( particle->Eta() < 1 ) ) FillHistogram("hMCgenSig0Eta1",  particle->Pt() );
            if( abs(rap) < fetaCut  ) {
                //	Int_t FlSig0eta1 = 1;
                
                Double_t pos[3] =  { particle->Vx(),  particle->Vy(),  particle->Vz() };
                //	Double_t cneg[3] =  { ph2c->X(),  ph2c->Y(),  ph2c->Z() };
                const AliESDVertex *esdVertex5 = fESDEvent->GetPrimaryVertex();
                //	Double_t vtx0[3] = {0,0,0}; // don't rely on ESD vertex, assume (0,0,0)
                Double_t vtx5[3];
                vtx5[0] = esdVertex5->GetX();
                vtx5[1] = esdVertex5->GetY();
                vtx5[2] = esdVertex5->GetZ();
                Double_t xvpos =  vtx5[0] - pos[0]/pos[2]*( pos[2] - vtx5[2] ) ;
                Double_t yvpos =  vtx5[1] - pos[1]/pos[2]*( pos[2] - vtx5[2] ) ;
                //	Double_t xvneg =  vtx5[0] - neg[0]/neg[2]*( cneg[2] - vtx5[2] ) ;
                //	Double_t yvneg =  vtx5[1] - neg[1]/neg[2]*( cneg[2] - vtx5[2] ) ;
                Double_t Rpos = sqrt ( (xvpos - vtx5[0])*(xvpos - vtx5[0]) + (yvpos - vtx5[1])*(yvpos - vtx5[1]) );
                //	Double_t Rneg = sqrt ( (xvneg - vtx5[0])*(xvneg - vtx5[0]) + (yvneg - vtx5[1])*(yvneg - vtx5[1]) );
                
                //	printf("VertXYZ %f %f %f Rpos %f \n \n", vtx5[0], vtx5[1], vtx5[2], Rpos );
                //	printf("SIG0 %f %f %f \n \n\n",  particle->Vx(),  particle->Vy(),  particle->Vz() );
                
                if( Rpos < 1 ){
                    FillHistogram("hMCgenSig0Eta1",  particle->Pt() );
                    if(       particle->GetPdgCode() == 3212 ) {
                        FillHistogram("hMCgenPSig0Eta1",  particle->Pt() );
                        FillHistogram("hMCgenPSig0Eta1Rap",  rap );
                    }
                    else if(  particle->GetPdgCode() == -3212) {
                        FillHistogram("hMCgenASig0Eta1",  particle->Pt() );
                        FillHistogram("hMCgenASig0Eta1Rap",  rap );
                    }
                    FillHistogram("hMCRposGen", Rpos );
                }
            }
            
	    //            FillHistogram("hMCsig0MassPt0All", m,  particle->Pt() );
            // FillHistogram("hMCsig0RPt0All", particle->R(), particle->Pt()  );
            
            TParticle* daught1 = (TParticle *)fStack->Particle( Daughter1 );
            TParticle* daught2 = (TParticle *)fStack->Particle( Daughter2 );
            
	    //            printf("Sig0 E %f P=%f Pt=%f M=%f \n", mE, mP ,mPt, mM );
            if ( !(abs( daught1->GetPdgCode() ) == 3122 && abs( daught2->GetPdgCode())  == 22) ) { //23apr14-open
                Double_t pos2[3]= { daught1->Px(),  daught1->Py(),  daught1->Pz() };
                Double_t neg2[3]= { daught2->Px(),  daught2->Py(),  daught2->Pz() };
                Double_t moth2[3]= { particle->Px(), particle->Py(),  particle->Pz() };
                Double_t arpod[2]= {0,0};
                
                GetArPod( pos2, neg2, moth2, arpod );
		/*               if( (abs(arpod[1]) > 0.6 && arpod[0]< 0.12) ) {
                    FillHistogram("hMCccMvsArpodPt",m,  particle->Pt() ) ;
                    FillHistogram("hMCArPodSig03",  arpod[1], arpod[0] ) ;
		    } */
            }
            
            if ( abs( daught1->GetPdgCode() ) == 3122 && abs( daught2->GetPdgCode())  == 22 ) {
                FillHistogram("hMCsig0PtLamdaVSGamma0",daught1->Pt(), daught2->Pt()  );
                
                Double_t pos2[3]= { daught1->Px(),  daught1->Py(),  daught1->Pz() };
                Double_t neg2[3]= { daught2->Px(),  daught2->Py(),  daught2->Pz() };
                Double_t moth2[3]= { particle->Px(), particle->Py(),  particle->Pz() };
                Double_t arpod[2]= {0,0};
                
                FillHistogram("hMCgenLamSig0", daught1->Pt() );
                Float_t rapLam = Rapidity(particle->Pt(), particle->Pz(), particle->GetMass());
                FillHistogram("hMCgenLamSig0Rap", rapLam );
                //	if( abs( particle->Eta() < 1 ) ) FillHistogram("hMCgenLamSig0Eta1",  daught1->Pt() );
                
                if( abs(rap) < fetaCut ) FillHistogram("hMCgenLamSig0Eta1",  daught1->Pt() );
                
                FillHistogram("hMCgenGamSig0", daught2->Pt() );

		// 25sep16
		//               	if( abs(  daught2->Eta() < 0.12 ) ) {
		if( abs(  daught2->Eta() ) < 0.12  ) {
		  FillHistogram("hMCgenGamSig0PHOS",  daught2->Pt() );
		  if( (daught2->Phi() > 220.*3.1415492/360. && daught2->Phi() < 320.*3.1415492/360.) ) { 
		  FillHistogram("hMCgenGamSig0PHOSEta1",  daught2->Pt() );
		  }
		}

                if( abs(rap) < fetaCut ) FillHistogram("hMCgenGamSig0Eta1",  daught2->Pt() );
                
                GetArPod( pos2, neg2, moth2, arpod );
                FillHistogram("hMCArPodSig0",  arpod[1], arpod[0] ) ;
                if( daught2->Energy() > 0.1 )  FillHistogram("hMCArPodSig02",  arpod[1], arpod[0] ) ;
                
		//                printf("Sig0->Gamma E %f Pt=%f E/pt=%f \n",daught2->Energy(), daught2->Pt(), daught2->Energy()/daught2->Pt() );
                
                
            }
        }        
    }
    
    //------------- now photons ----------------

    if( 1 ) return;

    for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
        TParticle* particle = (TParticle *)fStack->Particle(iTracks);
        if(particle->GetPdgCode() != 22)
        continue ;
        
        if(particle->R() >kRcut)
        continue ;
        
        if(TMath::Abs(particle->Eta())>0.9)
        continue ;
        
        Double_t pt = particle->Pt() ;
        //Total number of pi0 with creation radius <1 cm
        FillHistogram("hMCSigma0Phot",pt) ;
        
        //     Int_t mod ;
        //     Double_t x,z ;
	// Bool_t hitPHOS = 0; //  fPHOSgeom->ImpactOnEmc(particle, mod, z,x) ;
	//   Bool_t hitEMCAL= 0; //  fEMCALgeom->Impact(particle) ;
        
        //Photons in PHOS/EMCAL acceptance
        //      if(hitPHOS)
        //        FillHistogram("hMCSigma0GammaPHOSacc",pt) ;
        //      if(hitEMCAL)
        //        FillHistogram("hMCSigma0GammaEMCALacc",pt) ;
        
        //number of photons converted
        Bool_t converted = kFALSE ;
        if(particle->GetNDaughters()==2){
            TParticle * e1=fStack->Particle(particle->GetFirstDaughter()) ;
            TParticle * e2=fStack->Particle(particle->GetLastDaughter()) ;
            if(TMath::Abs(e1->GetPdgCode())==11 && TMath::Abs(e2->GetPdgCode())==11){ //conversion
                if(e1->R()<180.)
                converted = kTRUE ;
            }
        }
        if(converted)
        FillHistogram("hMCSigma0GammaConv",pt) ;
        
        printf(" Converted photons with V0 \n");
        TLorentzVector pConv ;
        Bool_t foundV0=kFALSE ;
        for(Int_t iv0=0; iv0<fESDEvent->GetNumberOfV0s();iv0++){
            AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;
            
            TParticle * negativeMC = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(v0->GetNindex())->GetLabel()));
            TParticle * positiveMC = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(v0->GetPindex())->GetLabel()));
            
            if(negativeMC && positiveMC){
                if(negativeMC->GetMother(0) != positiveMC->GetMother(0))
                continue ;
            }
            
            if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
                continue;
            }
            if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
                continue;
            }
            
            TParticle * GamGen = fStack->Particle(negativeMC->GetMother(0));
            Bool_t same = (GamGen == particle) ;
            TParticle * tmp = GamGen ;
            while(!same && tmp->GetFirstMother()>=0){
                tmp = fStack->Particle(tmp->GetFirstMother());
                same = (tmp == particle) ;
            }
            if(same){
                foundV0 = kTRUE ;
                const AliExternalTrackParam * paramPos = v0->GetParamP() ;
                const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
                AliKFParticle negKF(*paramNeg,11);
                AliKFParticle posKF(*paramPos,-11);
                pConv.SetXYZM(negKF.Px()+posKF.Px(),negKF.Py()+posKF.Py(),negKF.Pz()+negKF.Pz(),0.) ;
                break ;
            }
        }
	/*        if(foundV0){
            FillHistogram("hMCSigma0GammaV0",pt) ;
            //       FillHistogram("hMCSigma0GammaV0_devsE",
            //     (particle->Energy()-pConv.E())/particle->Energy(),particle->Energy()) ;
        }
	*/
        
        //Registered in PHOS/EMCAL
	/*        Bool_t cluInPHOS = kFALSE, cluInEMCAL=kFALSE ;
        TLorentzVector pCalo ;
        for (Int_t i=0; i<fESDEvent->GetNumberOfCaloClusters(); i++) {
            AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(i);
            Int_t iprim = clu->GetLabel() ; //# of particle hit PHOS/EMCAL
            Bool_t matched = kFALSE ;
            while(iprim>=0 ) {
                if(iprim==iTracks){
                    matched=kTRUE ;
                    break ;
                }
                else{
                    iprim=fStack->Particle(iprim)->GetFirstMother() ;
                }
            }
            if(!matched)
            continue ;
            if(clu->IsPHOS() && hitPHOS){
                cluInPHOS=kTRUE ;
                clu->GetMomentum(pCalo ,vtx);
                break ;
            }
            if(!clu->IsPHOS() && hitEMCAL){
                cluInEMCAL=kTRUE ;
                clu->GetMomentum(pCalo ,vtx);
                break ;
            }
        }
	*/
        
        //      if(cluInPHOS){
        //        FillHistogram("hMCSigma0Gamma_PHOSclu",pt) ;
        //        FillHistogram("hMCSigma0Gamma_PHOSclu_recE",pCalo.E()) ;
        
        //        FillHistogram("hMCSigma0Gamma_PHOSclu_devsE",
        // 		     (particle->Energy()-pCalo.E())/particle->Energy(),particle->Energy()) ;
        //      }
        //      if(cluInEMCAL){
        //        FillHistogram("hMCSigma0Gamma_EMCALclu",pt) ;
        //        FillHistogram("hMCSigma0Gamma_EMCALclu_recE",pCalo.E()) ;
        //        FillHistogram("hMCSigma0Gamma_EMCALclu_devsE",
        // 		     (particle->Energy()-pCalo.E())/particle->Energy(),particle->Energy()) ;
        //      }
    }
}


//_____________________________________________________________________________
void AliAnalysisTaskSigma0::FillHistogram(const char * key,Double_t x)const{
    //FillHistogram
    TH1F * tmp = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
    
    if(!tmp){
        AliInfo(Form("can not find histogram <%s> ",key)) ;
        return ;
    }
    tmp->Fill(x) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma0::FillHistogram(const char * key,Double_t x,Double_t y)const{
    //FillHistogram
    TObject * tmp = fOutputContainer->FindObject(key) ;
    if(!tmp){
        AliInfo(Form("can not find histogram <%s> ",key)) ;
        return ;
    }
    if(tmp->IsA() == TClass::GetClass("TH1F")){
        ((TH1F*)tmp)->Fill(x,y) ;
        return ;
    }
    if(tmp->IsA() == TClass::GetClass("TH2F")){
        ((TH2F*)tmp)->Fill(x,y) ;
        return ;
    }
    AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskSigma0::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
    //Fills 1D histograms with key
    TObject * tmp = fOutputContainer->FindObject(key) ;
    if(!tmp){
        AliInfo(Form("can not find histogram <%s> ",key)) ;
        return ;
    }
    if(tmp->IsA() == TClass::GetClass("TH2F")){
        ((TH2F*)tmp)->Fill(x,y,z) ;
        return ;
    }
    if(tmp->IsA() == TClass::GetClass("TH3F")){
        ((TH3F*)tmp)->Fill(x,y,z) ;
        return ;
    }
}
//______________________________________________________________________________
Double_t AliAnalysisTaskSigma0
::PlanarityAngle(const AliExternalTrackParam * paramPos,const AliExternalTrackParam * paramNeg)const {
    //calculate angle between e+e- plain and perpendicular to MF
    //We need sign of MagField to calculate orienation
    
    TVector3 u(paramPos->Px()+paramNeg->Px(),paramPos->Py()+paramNeg->Py(),paramPos->Pz()+paramNeg->Pz()) ;
    u.Unit() ;
    TVector3 vPos(paramPos->Px(),paramPos->Py(),paramPos->Pz()) ;
    TVector3 vNeg(paramNeg->Px(),paramNeg->Py(),paramNeg->Pz()) ;
    TVector3 v=vPos.Cross(vNeg) ;
    TVector3 w = u.Cross(v);
    TVector3 z(0,0,1.);
    TVector3 ua=u.Cross(z);
    Double_t wa = w.Angle(ua);
    Double_t mfield=fESDEvent->GetMagneticField() ;
    if(mfield>0.)
    return wa;      //normal field
    else
    return TMath::Pi()-wa ; //reverse field
    
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskSigma0::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz){
    //Check if this channel belogs to the good ones
    
    if(strcmp(det,"PHOS")==0){
        if(mod>5 || mod<1){
            AliError(Form("No bad map for PHOS module %d ",mod)) ;
            return kTRUE ;
        }
        if(!fPHOSBadMap[mod]){ //Bad map not set
            return kTRUE ;
        }
        if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0) {
            //      printf(".,");
            return kFALSE ;
        }
        else
        return kTRUE ;
    }
    
    /*  else{
     if(strcmp(det,"EMCAL")==0){
     if(mod>9 || mod<0){
     AliError(Form("No bad map for EMCAL module %d ",mod)) ;
     return kTRUE ;
     }
     if(!fEMCALBadMap[mod]){ //Bad map not set
     return kTRUE ;
     }
     if(fEMCALBadMap[mod]->GetBinContent(ix,iz)>0)
     return kFALSE ;
     else
     return kTRUE ;
     }
     else{
     AliError(Form("Can not find bad channels for detector %s ",det)) ;
     }
     }
     */
    
    return kTRUE ;
}
// //______________________________________________________________________________
// void AliAnalysisTaskSigma0::GetArmenterosQtAlfa(AliKFParticle* positiveKFParticle, AliKFParticle * negativeKFParticle, AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] ){
//   //see header file for documentation

//   TVector3 momentumVectorPositiveKF(positiveKFParticle->GetPx(),positiveKFParticle->GetPy(),positiveKFParticle->GetPz());
//   TVector3 momentumVectorNegativeKF(negativeKFParticle->GetPx(),negativeKFParticle->GetPy(),negativeKFParticle->GetPz());
//   TVector3 vecV0(gammaKFCandidate->GetPx(),gammaKFCandidate->GetPy(),gammaKFCandidate->GetPz());

//   Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
//   Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));

//   Float_t alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
//     ((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;


//   Float_t qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);

//   armenterosQtAlfa[0]=qt;
//   armenterosQtAlfa[1]=alfa;

// }

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigma0::TestLambda(Double_t pt,Double_t l1,Double_t l2){
    
    Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
    Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
    Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
    Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
    Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
    Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
    0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
    0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
    return (R2<2.5*2.5) ;    
}


//_____________________________________________________
void AliAnalysisTaskSigma0::FillCorr(TLorentzVector * trig, const Int_t itype  )
{
  //OLD     
  //  Int_t hnum = -1;
  if( 1>0 ) return;
}

//______________________________________________________________________________
void AliAnalysisTaskSigma0::GetArPod( Double_t pos[3], Double_t neg[3], Double_t moth[3],  Double_t arpod[2] ){
    
    //see header file for documentation
    
    TVector3 momentumVectorPositiveKF(pos[0],pos[1],pos[2]);
    TVector3 momentumVectorNegativeKF(neg[0],neg[1],neg[2]);
    TVector3 vecV0(moth[0],moth[1],moth[2]);
    
    Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
    Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));
    
    Float_t alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
    ((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
    
    Float_t qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);
    
    arpod[0]=qt;
    arpod[1]=alfa;
    
}

//____________________________________________________________________
Double_t  AliAnalysisTaskSigma0::Rapidity(Double_t pt, Double_t pz, Double_t m)
{
    //
    // calculates rapidity keeping the sign in case E == pz
    //
    
    Double_t energy = TMath::Sqrt(pt*pt+pz*pz+m*m);
    if (energy != TMath::Abs(pz))
    return 0.5*TMath::Log((energy+pz)/(energy-pz));
    
    Printf("W- mt=0");
    return TMath::Sign(1.e30,pz);
}

