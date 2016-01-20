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
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliCaloParticle.h"
//#include "AliTriggerAnalysis.h"
#include "AliMultiplicity.h"
#include "AliESDcascade.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliV0ReaderV1.h"
#include "AliVEvent.h"
#include "AliKFConversionPhoton.h"
#include "AliConversionCuts.h"
#include "TParticle.h"
#include "AliCentrality.h"

class Riostream;
class TFile;

ClassImp(AliAnalysisTaskSigma0)


//______________________________________________________________________________
AliAnalysisTaskSigma0::AliAnalysisTaskSigma0():
AliAnalysisTaskSE(),
    fV0Reader(NULL),
    fESDEvent(NULL),
    fESDpid(NULL),
    fESDtrackCuts(NULL),
    fStack(NULL),
    fOutputContainer(NULL),
    fReaderGammas(NULL),
    fTreeList(NULL),
    fSigma0(NULL),
    fTreeEvent(NULL),
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
    fLeadingTrack(-1),
    fLeadingPHOS(-1),
    fLeadingEMCAL(-1),
    fLeadingConv(-1),
    fAbsLeading(-1),
    fLeadingGenpi0(-1),
    fELeadingTrack(0.),
    fELeadingPHOS(0.),
    fELeadingEMCAL(0.),
    fELeadingConv(0.),
    fELeadingGenpi0(0.),
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
    fGammaEta(0),
    fGammaArmPt(0),
    fGammaArmAlpha(0),
    fGammaZConv(0),
    fGammaChi2(0),
    fGammaRadius(0),

    fSigmaMass(0),
    fSigmaPt(0),
    fSigmaArmPt(0),
    fSigmaArmAlpha(0),
    fCentralityV0M(0)
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
    fV0Reader(NULL),
    fESDEvent(NULL),
    fESDpid(NULL),
    fESDtrackCuts(NULL),
    fStack(NULL),
    fOutputContainer(NULL),
    fReaderGammas(NULL),
    fTreeList(NULL),
    fSigma0(NULL),
    fTreeEvent(NULL),
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
    fLeadingTrack(-1),
    fLeadingPHOS(-1),
    fLeadingEMCAL(-1),
    fLeadingConv(-1),
    fAbsLeading(-1),
    fLeadingGenpi0(-1),
    fELeadingTrack(0.),
    fELeadingPHOS(0.),
    fELeadingEMCAL(0.),
    fELeadingConv(0.),
    fELeadingGenpi0(0.),
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
    fptCut(0.070),
    fptMaxCut(1.500),
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
    fGammaEta(0),
    fGammaArmPt(0),
    fGammaArmAlpha(0),
    fGammaZConv(0),
    fGammaChi2(0),
    fGammaRadius(0),

    fSigmaMass(0),
    fSigmaPt(0),
    fSigmaArmPt(0),
    fSigmaArmAlpha(0),
    fCentralityV0M(0)
{
    // Common I/O in slot 0
    DefineInput (0, TChain::Class());
    //  DefineOutput(3, TTree::Class());
    //  DefineOutput(2, TTree::Class());
    
    
    // Your private output
    DefineOutput(1, TList::Class());
    
    // Initialize the PHOS geometry
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
    fEtaCuts[0] = 1.4;
    fEtaCuts[1] = 1.2;
    fEtaCuts[2] = 1.0;
    fEtaCuts[3] = 0.8;
    fEtaCuts[4] = 0.6;
    fEtaCuts[5] = 0.4;
    
    
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
        
        if(fTreeList){
            fTreeList->Clear();
            delete fTreeList;
        }
     
        if(fTreeEvent){
            fTreeEvent->Clear();
            delete fTreeEvent;
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
    fOutputContainer->Add(new TH1F("hQAEvents","Events processed",1,0.,1.)) ;
    fOutputContainer->Add(new TH1F("hQAEventsTrig","Events processed",21,-0.25,20.25)) ;

//  TH1I *hEventTrigger = new TH1I("hEventTrigger","Events processed",11,0.5,10.5);
//   hEventTrigger->GetXaxis()->SetBinLabel(1,"total");
//   fOutputContainer->Add(hEventTrigger);
    
    fOutputContainer->Add(new TH1F("R2Conv","R2 of V0",400,0.,200)) ;
    fOutputContainer->Add(new TH1F("R3Conv","R3 of V0",400,0.,200)) ;
    fOutputContainer->Add(new TH1F("v0sum","v0sum",200,0,20000));
    
    
 
//------------------------------------------------
// setup local variables : number of bin, minimum, maximum
//------------------------------------------------
    char key[55] ;

    Int_t npt=500 ; // Number of Pt bins
    Int_t npt2=150 ; // Double_t ptmin2 = 0. ;
    Double_t ptmin=0.;
    Double_t ptmax=15. ;   // was 50.
    Double_t ptmax2=15. ;
    
    Int_t nBinMass = 70; // Number of mass bins
    Int_t nR=50 ; // Number of radius bin
    Double_t Rmax=500. ;

    
  
   
    Double_t mSigmax2=1.4;  Double_t mSigmin2=1.1;  Int_t nBinMass2 = 280 ;
    
    //------------------------------------------------

    for(Int_t isol=0;   isol<1 ; isol++){       //  isol<12; isol++){   -offi
        
        sprintf(key,"hMCsig0MassPt%dAll",isol) ;
        fOutputContainer->Add(new TH2F(key,"Mass vs pt Sigma0",nBinMass,1.,1.7,npt,0.,ptmax)) ;
  
        sprintf(key,"hMCsig0RPt%dAll",isol) ;
        fOutputContainer->Add(new TH2F(key,"R vs pt Sigma0",npt,0.,Rmax, nR,0.,ptmax )) ;
        
        sprintf(key,"hMCsig0PtLamdaVSGamma%d",isol) ;
        fOutputContainer->Add(new TH2F(key,"pT L  vs g from Sigma0",nR,0., 5, nR,0., 2.5 )) ;
                
        sprintf(key,"hMCrecSig0PtLamdaVSGamma%d",isol) ;
        fOutputContainer->Add(new TH2F(key,"pT L  vs g from Sigma0rec",nR,0., 5, nR,0., 2.5 )) ;
        
        sprintf(key,"hMCrec2Sig0PtLamdaVSGamma%d",isol) ;
        fOutputContainer->Add(new TH2F(key,"pT L  vs g from Sigma0rec2",nR,0., 5, nR,0., 2.5 )) ;
        
        sprintf(key,"hMCrec3Sig0PtLamdaVSGamma%d",isol) ;
        fOutputContainer->Add(new TH2F(key,"pT L  vs g from Sigma0rec2",nR,0., 5, nR,0., 5. )) ;
        
        sprintf(key,"hMClam0MassPt%dAll",isol) ;
        fOutputContainer->Add(new TH2F(key,"Mass vs pt Lambda0",nBinMass,1.,1.7,npt,0.,ptmax)) ;
        
        sprintf(key,"hMClam0RPt%dAll",isol) ;
        fOutputContainer->Add(new TH2F(key,"Radius vs pt Lambda0",npt,0.,Rmax,  nR,0.,ptmax)) ;
    }
    
    fOutputContainer->Add(new TH1F("hQAMult","Multiplicity",400,-10.,40000.)) ;
    fOutputContainer->Add(new TH1F("hQACentr1","Centr1",400,-10.,400.)) ;
    fOutputContainer->Add(new TH1F("hQACentr2","Centr2",400,-10.,400.)) ;
    fOutputContainer->Add(new TH1F("hQACentr3","Centr3",400,-10.,400.)) ;


    fOutputContainer->Add(new TH1F("hConvMult","Conv Multiplicity",51,-0.5,50.5)) ;
    fOutputContainer->Add(new TH1F("hLamMult","Lam Multiplicity",51,-0.5,50.5)) ;
    
    fOutputContainer->Add(new TH1F("hNlamEv","N Lam in ev",21,-0.5,20.5)) ;
    fOutputContainer->Add(new TH1F("hNalamEv","N aLam in ev",21,-0.5,50.5)) ;
    fOutputContainer->Add(new TH1F("hNlamminalamEv","N Lam aLam in ev",41,-20.5,20.5)) ;
    
  /*
    fOutputContainer->Add(new TH1F("hT0AmC","T0 A-C", 400, -2000, 2000 )) ;
    fOutputContainer->Add(new TH1F("hT0ApC","T0 A+C", 400, -2000, 2000 )) ;
    
    fOutputContainer->Add(new TH1F("h2T0AmC","2nd T0 A-C", 400, -2000, 2000 )) ;
    fOutputContainer->Add(new TH1F("h2T0ApC","2nd T0 A+C", 400, -2000, 2000 )) ;
 */
    
    fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
    fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
    fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity",150,0.,150.));
    
    fOutputContainer->Add(new TH2F("hdEdxTrack","6dEdx Sig of accepted tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hLamdEdxTrack","7dEdx of accepted tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hAntiLamdEdxTrack","8dEdx of accepted tracks",100,0.,10.,150,0.,150.)) ;
    fOutputContainer->Add(new TH2F("hTotLamdEdxTrack","9dEdx of accepted tracks",100,0.,10.,150,0.,150.)) ;
    
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
    fOutputContainer->Add(new TH2F("hMCgenSig0RapEta","Primary Sig0 gen Rapidity Eta",30,-1.5, 1.5, 30,-1.5,1.5 )) ;
    fOutputContainer->Add(new TH2F("hMCgenSig0RapPt","Primary Sig0 gen Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
    fOutputContainer->Add(new TH2F("hMCgenrecSig0RapPt","Primary Sig0 genrec Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
    fOutputContainer->Add(new TH2F("hMCgenrec2Sig0RapPt","Primary Sig0 genrec2 Rapidity Pt",30,-1.5,1.5, 50,0.,15. )) ;
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
    
    fOutputContainer->Add(new TH2F("hMCgenrec2LamdMvsdPt","GenRec Lam dm vs dpt",40, -0.04, 0.04,     40,-1.,1.)) ;
    fOutputContainer->Add(new TH2F("hMCgenrec2ALamdMvsdPt","GenRec ALam dm vs dpt",40, -0.04, 0.04,   40,-1.,1.)) ;
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
    fOutputContainer->Add(new TH1F("hMCgenGamSig0","Primary Lambda Sig0 gen",npt,0.,ptmax)) ;
    fOutputContainer->Add(new TH1F("hMCgenGamSig0Eta1","Primary Lambda Sig0 gen eta1",npt,0.,ptmax)) ;
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
    
    Double_t mSigmax=1.4; Double_t mSigmin=1.1;
    Double_t mLammax=1.20; Double_t mLammin=1.05;
    
    
    fOutputContainer->Add(new TH2F("hLamMvsPt10","10DT m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hALamMvsPt10","10DT m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hLamMvsPt11","10DT m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hALamMvsPt11","10DT m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH2F("hPhi1LamMvsPt11","10DTphi1 m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hPhi1ALamMvsPt11","10DTphi1 m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH2F("hPhi2LamMvsPt11","10DTphi1 m vs pt",100, 1.11,1.13, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hPhi2ALamMvsPt11","10DTphi1 m vs pt",100, 1.11,1.13, 100,0.,ptmax)) ;
    
    fOutputContainer->Add(new TH1F("hRLam2","R Lam", 400, 0,100. )) ;
    fOutputContainer->Add(new TH1F("hRALam2","R ALam",400,0,100. )) ;
    
    fOutputContainer->Add(new TH2F("hLamPtPvsM11","Lam ptPvsM",   50, 0, 5, 50,0.,5. )) ;
    fOutputContainer->Add(new TH2F("hALamPtPvsM11","ALam ptPvsM",50, 0, 5, 50,0.,5. )) ;
    
    fOutputContainer->Add(new TH2F("hLamMvsPt12","12DT m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
    fOutputContainer->Add(new TH2F("hALamMvsPt12","12DT m vs pt",200, mLammin,mLammax, 100,0.,ptmax)) ;
        
    fOutputContainer->Add(new TH2F("hArPodLam","MC Armenteros-Podolanski Lam ;#alpha;p_{t}",100,-1.0,1.0, 70,0,0.7) );
    fOutputContainer->Add(new TH2F("hArPodALam","MC Armenteros-Podolanski ALam ;#alpha;p_{t}",100,-1.0,1.0, 70,0,0.7) );
    fOutputContainer->Add(new TH2F("hArPodGConv","MC Armenteros-Podolanski Gamma Conv ;#alpha;p_{t}",100,-1.0,1.0, 70,0,0.7) );
    fOutputContainer->Add(new TH1F("hCGamPtConv","Conv gamma Pt2",  npt*2, ptmin, 15. )) ;
    
    fOutputContainer->Add(new TH2F("hMCArPodSig0","Armenteros-Podolanski Sig0 ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    fOutputContainer->Add(new TH2F("hMCArPodSig02","Armenteros-Podolanski Sig02  ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
    fOutputContainer->Add(new TH2F("hMCArPodSig03","Armenteros-Podolanski NoSim Sig03  ;#alpha;p_{t}",100,-1.0,1.0, 80,0,0.8) );
        
    fOutputContainer->Add(new TH2F("hMCccMvsArpodPt","CC Mass vs ptnArPod",nBinMass, mSigmin,mSigmax ,npt,0.,ptmax));
    
    fOutputContainer->Add(new TH2F("hMCV0MvsWidth","MC w vs m",300,0.,TMath::Pi(),400,0.,1.)) ;
 
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
        
        sprintf(key,"hMCPSig0PrimPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCPSig0PrimPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCASig0PtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCASig0PtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCASig0PrimPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCASig0PrimPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCPLamPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCPLamPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCPLamPrimPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCPLamPrimPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCALamPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCALamPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        sprintf(key,"hMCALamPrimPtGen%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCALamPrimPtRec%d",is) ;
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        
        
        fOutputContainer->Add(new TH1F(key,key, npt,0.,ptmax)) ;
        sprintf(key,"hMCgenrec3MvsPtLamgenGamRec%d",is) ;
        fOutputContainer->Add(new TH2F(key,key, nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
        sprintf(key,"hMCgenrec3MvsPtLamrecGamGen%d",is);
        fOutputContainer->Add(new TH2F(key,key, nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
        sprintf(key,"hMCgenrec3MvsPtLamrecGamRec%d",is);
        fOutputContainer->Add(new TH2F(key,key, nBinMass2, mSigmin2,mSigmax2 ,npt2,0.,ptmax2));
        
    }
        
    // end of histos in  fOutputContainer
    fOutputContainer->SetName(GetName());
    
//------------------------------------------------
// fTree Branch definitions
//------------------------------------------------
  
//-----------BASIC-INFO---------------------------
    

    fTreeEvent = new TTree("fTreeEvent","Events");
/* 0*/    fTreeEvent->Branch("fLambdaMod",&fLambdaMod,"fLambdaMod/D");
/* 1*/    fTreeEvent->Branch("fLambdaPx",&fLambdaPx, "fLambdaPx/D");
/* 2*/    fTreeEvent->Branch("fLambdaPy",&fLambdaPy, "fLambdaPy/D");
/* 3*/    fTreeEvent->Branch("fLambdaPz",&fLambdaPz, "fLambdaPz/D");
/* 4*/    fTreeEvent->Branch("fLambdaMass",&fLambdaMass, "fLambdaMass/D");
/* 5*/    fTreeEvent->Branch("fLambdaCosPointingAngle",&fLambdaCosPointingAngle,"fLambdaCosPointingAngle/D");
/* 6*/    fTreeEvent->Branch("fLambdaDCADaughters",&fLambdaDCADaughters,"fLambdaDCADaughters/D");
/* 7*/    fTreeEvent->Branch("fLambdaDCAtoPVNeg",&fLambdaDCAtoPVNeg,"fLambdaDCAtoPVNeg/D");
/* 8*/    fTreeEvent->Branch("fLambdaDCAtoPVPos",&fLambdaDCAtoPVPos,"fLambdaDCAtoPVPos/D");
/* 9*/    fTreeEvent->Branch("fLambdaRadius",&fLambdaRadius,"fLambdaRadius/D");
/*10*/    fTreeEvent->Branch("fLambdaArmPt",&fLambdaArmPt,"fLambdaArmPt/D");
/*11*/    fTreeEvent->Branch("fLambdaArmAlpha",&fLambdaArmAlpha,"fLambdaArmAlpha/D");
/*12*/    fTreeEvent->Branch("fLambdaEta",&fLambdaEta,"fLambdaEta/D");
    
/*13*/    fTreeEvent->Branch("fGammaPx",&fGammaPx, "fGammaPx/D");
/*14*/    fTreeEvent->Branch("fGammaPy",&fGammaPy, "fGammaPy/D");
/*15*/    fTreeEvent->Branch("fGammaPz",&fGammaPz, "fGammaPz/D");
/*16*/    fTreeEvent->Branch("fGammaMass",&fGammaMass, "fGammaMass/D");
/*17*/    fTreeEvent->Branch("fGammaCosPointingAngle",&fGammaCosPointingAngle,"fGammaCosPointingAngle/D");
/*18*/    fTreeEvent->Branch("fGammaRadius",&fGammaRadius,"fGammaRadius/D");
/*19*/    fTreeEvent->Branch("fGammaArmPt",&fGammaArmPt,"fGammaArmPt/D");
/*20*/    fTreeEvent->Branch("fGammaArmAlpha",&fGammaArmAlpha,"fGammaArmAlpha/D");
/*21*/    fTreeEvent->Branch("fGammaEta",&fGammaEta,"fGammaEta/D");
/*22*/    fTreeEvent->Branch("fGammaChi2",&fGammaChi2,"fGammaChi2/D");
/*23*/    fTreeEvent->Branch("fGammaZConv",&fGammaZConv,"fGammaZConv/D");
    
/*24*/    fTreeEvent->Branch("fSigmaMass",&fSigmaMass,"fSigmaMass/D");
/*25*/    fTreeEvent->Branch("fSigmaPt",&fSigmaPt, "fSigmaPt/D");
/*26*/    fTreeEvent->Branch("fSigmaArmPt",&fSigmaArmPt,"fSigmaArmPt/D");
/*27*/    fTreeEvent->Branch("fSigmaArmAlpha",&fSigmaArmAlpha,"fSigmaArmAlpha/D");
/*28*/    fTreeEvent->Branch("fCentr",&fCentr,"fCentr/D");
    
    fTreeList->Add(fTreeEvent);
    
    PostData(1, fOutputContainer);
    //  PostData(2, fOutputContainer);
    //  PostData(3, fOutputContainer);
    
    
    //   PostData(2, fSigma0);
    //   PostData(3, fTreeEvent);
  
}

//======================================================
//_____________________________________________________
void AliAnalysisTaskSigma0::UserExec(Option_t */*option*/)
{
    // Called for each event
    
      printf(" Execute analysis for current event... -1 \n");
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
    
    fLeadingTrack=-1;
    fLeadingPHOS=-1;
    fLeadingEMCAL=-1;
    fLeadingConv=-1;
    fAbsLeading=-1;
    fELeadingTrack=0.;
    fELeadingPHOS=0.;
    fELeadingEMCAL=0.;
    fELeadingConv=0.;
    
    //First try to find Stack information.
    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
        if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent())
        fStack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
    }
   
    // Connect to the InputEvent
    // Appropriate for ESD analysis !
    
    AliESDInputHandler *esdHandler=dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
  
    
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
    

    
    
    //Fill MC histograms if MC is present
    fAbsLeading = -100;
    ProcessMC();
    
    FillHistogram("hQAEventsTrig",13.) ;
 //   FillHistogram("hQAEventsTrig",14) ;
    
    // Get Centrality information
    // fGetCent = fESDEvent->GetCentrality();
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

    printf(" cent %f %f %f \n", centralityV0M, centralityV0A, lCentrality);

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
    FillHistogram("hQAEventsTrig",8.) ;
    
    
  
//------------------------------------------------
// Select tracks
//------------------------------------------------

    //   SelectTracks() ;

    SelectConvPhotons() ;
    Int_t nConv=fConvEvent->GetEntriesFast() ;
    if( nConv < 1 ) return;
  
    SelectLambda ();
    Int_t nLam=fGenpi0Event->GetEntriesFast() ;
    if( nLam < 1 ) return;

    //    printf("CONV nConv %d nLam %d \n", nConv, nLam );

   //  SelectPHOSPhotons() ;
    //  SelectPhotonsFB();
//------------------------------------------------
    
    

 
    FillHistogram("hConvMult", nConv) ;
    FillHistogram("hLamMult", nLam) ;
    FillHistogram("hQAEventsTrig",15) ;
    
    
    
    
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
    fELeadingPHOS=0.;  fLeadingPHOS = -1 ;
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
        if(fELeadingPHOS<p.Pt()){
            fLeadingPHOS = inPHOS ;
            fELeadingPHOS=p.Pt();
        }
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
    }
    //QA: number of clusters/event in run //QA: average cluster energy
} // END of SelectPHOSPhotons(){

//______________________________________________
//______________________________________________________
//_____________________________________________________________
//______________________________________________________________________
void AliAnalysisTaskSigma0::SelectLambda(){
    //Fill list of Lambdas
    
    fELeadingConv=0.;
    fLeadingConv = -1 ;
    
    
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
    
    Int_t NlamEv = 0;  Int_t NalamEv = 0;
    Double_t lampospxprev = -100 ;
    Double_t lamnegpxprev = -100 ;
    Double_t alampospxprev = -100 ;
    Double_t alamnegpxprev = -100 ;
    Double_t MassLam = 1.115683 ;
    Double_t MassLamLim = 0.050 ;
    Double_t MassPi = 0.13957018 ;
    Double_t MassPr = 0.938272046 ;
    Double_t mlamprev = 0;  Double_t malamprev = 0;
    
    Double_t ptminLam =  0.25;
    //    Double_t MminLam = 1.080;       Double_t MmaxLam = 1.150;
    Double_t MminLam = 1.110;       Double_t MmaxLam = 1.122;
  
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
        if(pos->GetSign() == neg->GetSign())     continue ;
        
        
        if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) continue ;
        
        if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) || !(neg->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
     
        
        Int_t    lOnFlyStatus        = 0;
        Bool_t    isOnFly        = 0;
        // VO's main characteristics to check the reconstruction cuts
        lOnFlyStatus       = v0->GetOnFlyStatus();
        isOnFly = v0->GetOnFlyStatus();
        if( lOnFlyStatus == 0 ) continue;
        
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
        if ((lDcaPosToPrimVertex < 0.1) || (lDcaNegToPrimVertex < 0.1) ||(lV0cosPointAngle < 0.998) || (lV0Radius < 0.0) || (lV0Radius > 1000) ) continue;
        
        //select V0 finder  OPEN
        
        Double_t v0x=0.,v0y=0.,v0z=0.;
        v0->GetXYZ(v0x,v0y,v0z) ;
        Double_t r=TMath::Sqrt(v0x*v0x + v0y*v0y) ;
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
        
         FillHistogram("R2Conv", r ) ;
    
        
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
        
	if( (mlam > MminLam && mlam < MmaxLam && ptlam > ptminLam )  ) {
                
	  FillHistogram("hLamdEdxTrack", pos->GetP(), pos->GetTPCsignal()) ;
	  FillHistogram("hPhi2LamMvsPt11",mlam, ptlam );
	  FillHistogram("hRLam2", lV0Radius );
	  FillHistogram("hLamPtPvsM11",pos->Pt(), neg->Pt() );
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
	  FillHistogram("hALamPtPvsM11",pos->Pt(), neg->Pt() );
	  //	 printf("ALam-v0 m %f pt %f \n", malam, ptalam );
	  FillHistogram("hRecALam",ptalam );
	  if ( etaalam  < fetaCut ) 	FillHistogram("hRecALamEta1",ptalam );
	  // Int_t ialam = alamKF.GetLabel() ;
	  //	 printf("ilam %d \n", ialam );
	}
            
	FillHistogram("hLamMvsPt10",mlam, ptlam );
	FillHistogram("hALamMvsPt10",malam, ptalam );
	FillHistogram("R3Conv", r ) ;
            
	Double_t Rlam0 =   sqrt( v0->Xv()*v0->Xv() +  v0->Yv()*v0->Yv() + v0->Zv()*v0->Zv() );
	FillHistogram("hMClam0MassPt0All", v0->GetEffMass(),  v0->Pt() );
	FillHistogram("hMClam0RPt0All", Rlam0,  v0->Pt() );
            
            
	ModLam = 0;
	if  ( (mlam > MminLam && mlam < MmaxLam && ptlam > ptminLam ) ||
	      (malam > MminLam && malam < MmaxLam && ptalam> ptminLam )
	      ) {
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
	} // end include  Lambda and  anti-ALambda
	//	printf("\n ModLam %f \n", ModLam );
            
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
	//Fill MC information for Lambda
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
                        
                        
	      //	     printf("Pdg1 %d 2 %d \n", negMC->GetPdgCode(), posMC->GetPdgCode() );
	      //	if(  !fStack->Particle(negMC->GetMother(0)) || !fStack->Particle(posMC->GetMother(0)) ) continue;
                        
	      TParticle * v0Lam = fStack->Particle(negMC->GetMother(0) ) ;
                        
	      if( !v0Lam ) continue;
	      Double_t rapLam = Rapidity(  v0Lam->Pt(),  v0Lam->Pz(),  v0Lam->GetMass() );
                        
	      //???-23jul14
	      if( v0Lam->IsPrimary() ) {              //have NO plots in  10f6a-PYTHiA, but have ones in 10f6 - PHOJET, 23jul14
		if( ModLam > 0 ) FillHistogram("hLamMvsPt12", v0Lam->GetMass(), v0Lam->Pt() );
		else if ( ModLam < 0 ) FillHistogram("hALamMvsPt12",  v0Lam->GetMass(), v0Lam->Pt() );
		continue; //???-23jul14  // Goes through in mc10f6a, but breaks in mc10f6, WHY? - open
	      } // /???-23jul14 }
	      
                        
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
		    FillHistogram("hMCgenrecSig0RapPt", rap, v0Sig->Pt()  );
                                    
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
		      FillHistogram("hMCrecSig0PtLamdaVSGamma0",daught1->Pt(), daught2->Pt()  );
		      FillHistogram("hRec4LamSig0",  daught1->Pt() );
		      // if(  abs( v0Sig->Eta() ) < 1 )  FillHistogram("hRec4LamSig0Eta1", daught1->Pt() );
		      if(  abs(rap) < fetaCut ) FillHistogram("hRec4LamSig0Eta1",  daught1->Pt() );
		    }
		  }
		}   // generated Sig0 from recostructed Lam
	      }  // 1 more if to enshure that v0Lam as mother exists!!!
	    }
	  }
	}   // end fstack - excustion into MC
            
	if ( r < 10 ) {
	  
	  FillHistogram("hLamMvsPt11",mlam, ptlam );
	  FillHistogram("hALamMvsPt11",malam, ptalam );
	  
	  if (  pos->GetP() > 0.5 &&  neg->GetP() > 0.5 ) {
	    FillHistogram("hTotLamdEdxTrack", pos->GetP(), pos->GetTPCsignal()) ;
	    FillHistogram("hTotLamdEdxTrack", neg->GetP(), neg->GetTPCsignal()) ;
	  }
	}
        // FULL end including MC for Lambdas  ----------------------
    }
}

//________________________________________________________________________
void AliAnalysisTaskSigma0::SelectPhotonsFB()
{
    
   
    Int_t inConv=0 ;
    
    FillHistogram("hCheckFB",0.5);
    
    
    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
    if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader // origianl
    
    
    fReaderGammas = fV0Reader->GetReconstructedGammas();
    
    if(fV0Reader)
    if((AliConversionCuts*)fV0Reader->GetConversionCuts())
    //		if(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
    //	 	fOutputContainer->Add(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
    //   printf("FB 2 \n");
    FillHistogram("hCheckFB",1.5);
    
    //    if(1>0 ) return;
    //    printf("VOREADER %d \n", fReaderGammas->GetEntriesFast() );
    
    for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
        
        FillHistogram("hCheckFB",2.5);
        
        AliKFConversionPhoton* PhotonCandidate = (AliKFConversionPhoton*) fReaderGammas->At(i);
        
        if(!PhotonCandidate) continue;
        
        TLorentzVector photFB;
        photFB.SetXYZM(PhotonCandidate->GetPx(),PhotonCandidate->GetPy(), PhotonCandidate->GetPz(),0.) ;
        
        FillHistogram("hCheckFB",3.5);
        
        FillHistogram("hAMbeforeFB",PhotonCandidate->GetArmenterosAlpha(), PhotonCandidate->GetArmenterosQt());
        
        if(inConv>=fConvEvent->GetSize())fConvEvent->Expand(inConv*2) ;
        new((*fConvEvent)[inConv]) AliCaloParticle(photFB) ;
        inConv++ ;
        
        FillHistogram("hCheckFB",4.5);
   }        
}
//_______________________________________________________________________________
void AliAnalysisTaskSigma0::SelectConvPhotons(){
    //Fill list of conversion photons, that is scan v0s and select photon-like
    fELeadingConv=0.;
    fLeadingConv = -1 ;
    
    Double_t MassPi = 0.13957018 ;
    Double_t MassPr = 0.938272046 ;
    
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
    Double_t MassSig0 = 1.192642 ;
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
    Int_t NlamEv = 0;   Int_t NalamEv = 0;  Int_t NlamminalamEv = 0;
    
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
            (lV0cosPointAngle < 0.998) || (lV0Radius < 0.0) || (lV0Radius > 1000) ) continue;  // (lV0Radius < 0.5) was before 1.10.13
        
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
	  
	  FillHistogram("hArPodGConv",  arpod[1], arpod[0] ) ;
	  FillHistogram("hCGamPtConv", GamRec.Pt() ) ;
            
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
                
                
                
	    if(negativeMC->GetUniqueID() != 5 || positiveMC->GetUniqueID() !=5){ // id 5 is conversion
	      continue; }
	    
	    if(GamGen->GetPdgCode() == 22){    // indeed this Gamma reconstracted from Gamma
	      FillHistogram("hMCV0MvsWidth",w,m) ;
	      //!!! add histo
	      
	      /* Yura 12.3.14 primary Gamma has no Mother ! */
	      if( GamGen->IsPrimary() ) continue;
	      TParticle * v0Sig2 = fStack->Particle(GamGen->GetMother(0));
	      
	      //       	printf("=========================mSig %f \n",  v0Sig2->GetMass() );
	      if( abs(v0Sig2->GetMass() - 1.192642)<0.01 ) {
		
		Int_t NDt1 =v0Sig2->GetNDaughters();
                
		if( abs( v0Sig2->GetPdgCode() ) == 3212  && NDt1 == 2)  {    // Sigma0
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
                            
		  if ( (abs( LamGen21->GetPdgCode() ) == 3122 && abs( GamGen22->GetPdgCode())  == 22)
		       //		 &&  (inLam >0 && LamRecFl > 0)
		       ) {
		    
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
		      //		  TLorentzVector sig0LgenGgenM =  GGen22m + LGen21m  ;
		      //		  FillHistogram("hMCgenrec3MvsPtLamgenGamgenM",sig0LgenGgenM.M(),sig0LgenGgenM.Pt() );
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
			      if (v0Sig2->IsPrimary() ) {
				sprintf(key,"hMCPSig0PrimPtRec%d",irap) ;
				FillHistogram(key,  v0Sig2->Pt() );
			      }
			    }
			    else  if(   LamGen21->GetPdgCode() == -3122 ) {
			      sprintf(key,"hMCASig0PtRec%d",irap) ;
			      FillHistogram(key,  v0Sig2->Pt() );
			      if (v0Sig2->IsPrimary() ) {
				sprintf(key,"hMCASig0PrimPtRec%d",irap) ;
				FillHistogram(key,  v0Sig2->Pt() );
			      }
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
    
    FillHistogram("hNlamEv",  NlamEv );
    FillHistogram("hNalamEv",  NalamEv );
    NlamminalamEv = NlamEv - NalamEv;
    FillHistogram("hNlamminalamEv", NlamminalamEv  );
} // END of SelectConvPhoton

//_____________________________________________________
//___________________________________________________________________________
void AliAnalysisTaskSigma0::SelectTracks(){
    //Select pi+- tracks for search of Sigma+- 1385, PID have to be applied later
    fELeadingTrack=0.;
    fLeadingTrack = -1 ;
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
        
        FillHistogram("hdEdxTrack", track->GetP(), track->GetTPCsignal()) ;
        
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
    //   printf("ProcessMC \n");
    
    const Double_t kRcut = 1. ; //cut for primary particles
    Double_t vtx[3];
    vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
    vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
    vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();
    
    Int_t Daughter1  = 0;  Int_t Daughter2  = 0;
    
    if(  fStack->GetNtrack() < 2 ) return;
    
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
        // printf("MC----itr %d status %d PID=%d Parent %d Daut1=%d Daut2=%d ", iTracks,mStatus,mPdgSign,iParent,Daughter1,Daughter2);
        // printf("MC E %f Pt=%f M=%f  Pdg %d \n", mE,  mPt, mM, mPdg );
        // printf("MC   Pdg %d  itr %d \n", mPdgSign, iTracks );     //    }
        
        if(  abs( particle->GetPdgCode() ) == 3122  && NDt == 2)  { // all Lambda generated
            
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
                        if( particle->IsPrimary() ) {
                            sprintf(key,"hMCPLamPrimPtGen%d",irap) ;
                            FillHistogram(key,  ptLam );
                        }
                    }
                    else  if( pdgLam<0 ) {
                        sprintf(key,"hMCALamPtGen%d",irap) ;
                        FillHistogram(key,  ptLam );
                        if( particle->IsPrimary() ) {
                            sprintf(key,"hMCALamPrimPtGen%d",irap) ;
                            FillHistogram(key,  ptLam );
                        }
                    }
                }
            }
            
            if( TMath::Abs(rapLam) < fetaCut ) {
                FillHistogram("hMCgenLamEta1",  particle->Pt() );
                if(  particle->GetPdgCode() == 3122 )    FillHistogram("hMCgenPLamEta1",  particle->Pt() );
                else if( particle->GetPdgCode() == -3122 )    FillHistogram("hMCgenALamEta1",  particle->Pt() );
            }
        }
        
        // Sigma0 generated
        if( abs( particle->GetPdgCode() ) == 3212  && NDt == 2)  {    // Sigma0
            Double_t m=0., mcalc=0. ;
            m =      particle->GetMass();
            mcalc = particle->GetCalcMass();
            Float_t rap = Rapidity(particle->Pt(), particle->Pz(), particle->GetMass());
            Float_t ptSig0 = particle->Pt();
            
            if( ptSig0 < 1 ) continue; // can reconstruct Sigma0 with pT>1 GeV/c ONLY!
            
            Float_t pdgSig0 = particle->GetPdgCode();
            
            //    printf(" m %f w %f cm %f PDG %d \n", m, width, mcalc, particle->GetPdgCode() == 3212 );
            FillHistogram("hMCgenSig0",  particle->Pt() );
            FillHistogram("hMCgenSig0Rap",  rap );
            FillHistogram("hMCgenSig0RapEta",  rap,  particle->Eta() );
            FillHistogram("hMCgenSig0RapPt", rap,  particle->Pt() );
            
            char key[55] ;
            for (Int_t irap = 0; irap < 6; irap++) {
                if( abs(rap) < fEtaCuts[irap] ) {
                    if( pdgSig0>0 ) {
                        sprintf(key,"hMCPSig0PtGen%d",irap) ;
                        FillHistogram(key,  ptSig0 );
                        
                        if( particle->IsPrimary() ) {
                            sprintf(key,"hMCPSig0PrimPtGen%d",irap) ;
                            FillHistogram(key,  ptSig0 );
                        }
                    }
                    else  if( pdgSig0<0 ) {
                        sprintf(key,"hMCASig0PtGen%d",irap) ;
                        FillHistogram(key,  ptSig0 );
                        if( particle->IsPrimary() ) {
                            sprintf(key,"hMCASig0PrimPtGen%d",irap) ;
                            FillHistogram(key,  ptSig0 );
                        }
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
            
            FillHistogram("hMCsig0MassPt0All", m,  particle->Pt() );
            FillHistogram("hMCsig0RPt0All", particle->R(), particle->Pt()  );
            
            TParticle* daught1 = (TParticle *)fStack->Particle( Daughter1 );
            TParticle* daught2 = (TParticle *)fStack->Particle( Daughter2 );
            
            //      printf("Sig0 E %f P=%f Pt=%f M=%f \n", mE, mP ,mPt, mM );
            if ( !(abs( daught1->GetPdgCode() ) == 3122 && abs( daught2->GetPdgCode())  == 22) ) { //23apr14-open
                Double_t pos2[3]= { daught1->Px(),  daught1->Py(),  daught1->Pz() };
                Double_t neg2[3]= { daught2->Px(),  daught2->Py(),  daught2->Pz() };
                Double_t moth2[3]= { particle->Px(), particle->Py(),  particle->Pz() };
                Double_t arpod[2]= {0,0};
                
                GetArPod( pos2, neg2, moth2, arpod );
                if( (abs(arpod[1]) > 0.6 && arpod[0]< 0.12) ) {
                    FillHistogram("hMCccMvsArpodPt",m,  particle->Pt() ) ;
                    FillHistogram("hMCArPodSig03",  arpod[1], arpod[0] ) ;
                }
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
                //	if( abs( particle->Eta() < 1 ) ) FillHistogram("hMCgenGamSig0Eta1",  daught2->Pt() );
                if( abs(rap) < fetaCut ) FillHistogram("hMCgenGamSig0Eta1",  daught2->Pt() );
                
                GetArPod( pos2, neg2, moth2, arpod );
                FillHistogram("hMCArPodSig0",  arpod[1], arpod[0] ) ;
                if( daught2->Energy() > 0.1 )  FillHistogram("hMCArPodSig02",  arpod[1], arpod[0] ) ;
                
                //	printf("Sig0->Gamma E %f Pt=%f E/pt=%f \n",daught2->Energy(), daught2->Pt(), daught2->Energy()/daught2->Pt() );
                
                fAbsLeading = 100;
                
            }
        }        
    }
    
    //------------- now photons ----------------
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
        Bool_t hitPHOS = 0; //  fPHOSgeom->ImpactOnEmc(particle, mod, z,x) ;
        Bool_t hitEMCAL= 0; //  fEMCALgeom->Impact(particle) ;
        
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
        
        //Converted photons with V0
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
        if(foundV0){
            FillHistogram("hMCSigma0GammaV0",pt) ;
            //       FillHistogram("hMCSigma0GammaV0_devsE",
            //     (particle->Energy()-pConv.E())/particle->Energy(),particle->Energy()) ;
        }
        
        //Registered in PHOS/EMCAL
        Bool_t cluInPHOS = kFALSE, cluInEMCAL=kFALSE ;
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

//_________________________________________________________________________________
void AliAnalysisTaskSigma0::REvalIsolation(TLorentzVector * ph,const Int_t isolation ){
    // Int_t AliAnalysisTaskSigma0::EvalIsolation( AliCaloParticle * ph){
    
    //  printf("  Check if this particle is isolated ITYPE %d  \n",isolation);
    //We use several cone radii and epsilons looking on charged particles and EMCAL ones
    //   const Double_t kConeR1=0.2 ;   const Double_t kConeR2=0.3 ;   const Double_t kConeR3=0.4 ;
    //   const Double_t ptthresh=0.5 ;
    
    Double_t  phiTrig = ph->Phi();
    Double_t  etaTrig = ph->Eta();
    Double_t  MggTrig = ph->M();
    Double_t  PtTrig = ph->Pt();
    char key[55] ;
    
    
    Int_t nbin=9 ;
    Double_t Rcone[]={0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.55 } ;
    
    
    //  Int_t nbin=8 ;
    //  Double_t Rcone[]={0.05, 0.1, 0.2, 0.4, 0.6 , 0.8, 1.0 , 1.2, 1.4, 1.6, 1.25 } ;
    //  Double_t Rcone[]={0.05, 0.1, 0.2, 0.4, 0.6 , 0.8, 1.0 , 1.52 } ;
    
    Double_t FSumPtcone[nbin];    Double_t BSumPtcone[nbin];   Double_t BxEcone[nbin];  Double_t SSumPtcone[nbin]; 
    Double_t TSumPtcone = 0;    Double_t TVECSumPtcone = 0;  
    
    Int_t FNcone[nbin];    Int_t BNcone[nbin];    Int_t SNcone[nbin]; 
    //  Double_t xecone[nbin];  Double_t multall[nbin];
    
    for(Int_t i=0 ; i< nbin; i++){
        FSumPtcone[i] = 0;  BSumPtcone[i] = 0;  BxEcone[i] = 0; SSumPtcone[i] = 0;  
        FNcone[i] =0;  BNcone[i]=0 ;   SNcone[i]=0; 
    }
    
    Int_t nstart = 0;
    Int_t nfin = nbin;
    
    //  Int_t ConeSum =0;   Int_t ConeProp = 0;   Int_t SumPt = 0;
    Int_t ntracks=fTrackEvent->GetEntriesFast() ;
    
    for(Int_t itr=0; itr<ntracks; itr++){
        AliCaloParticle * track = (AliCaloParticle*)fTrackEvent->At(itr) ;
        if(track==ph) //this particle
        continue ;     
        
        nstart = 0;
        nfin = nbin;
        
        if ( isolation == -555 ){
            nstart = 9;
            nfin = 10;
        } 
        
        Double_t deleta = etaTrig - track->Eta() ;
        Double_t deletaneg = -etaTrig - track->Eta() ;
        
        Double_t delphi = phiTrig - track->Phi() ;      
        if(delphi <= -TMath::PiOver2()) delphi+=TMath::TwoPi();
        if(delphi > 3*TMath::PiOver2()) delphi-=TMath::TwoPi();
        
        for(Int_t icone= nstart; icone<nfin; icone++){
            Double_t i2 = icone ; 
            Double_t a2 = 2 ;  
            Double_t kind2 = pow( a2, i2) ;
            Int_t kind3 = float( kind2 );
            
            if(  ( isolation>=0 && (isolation&kind3) == kind3 ) || isolation== -555  ){
                
                //	if( MggTrig > 0.075 &&  MggTrig < 0.255 )  FillCorr( ph, 20+icone ) ;
                //	if( MggTrig > 0.11 &&  MggTrig < 0.16 ) {   // 12jul13
                
                if( MggTrig > 0.115 &&  MggTrig < 0.155 ) {
                    sprintf(key,"PP_PhiCor_peak_%d",icone) ;
                    FillHistogram(key, delphi ,PtTrig) ;
                    
                    sprintf(key,"PP_r_BxE_peak_%d",icone) ; 
                    if( PtTrig >= 6 && delphi> 0.66*TMath::Pi() && delphi<1.33*TMath::Pi() ) FillHistogram(key,track->Pt()*TMath::Cos(delphi)/PtTrig, PtTrig);
                    
                    //	  FillCorr(ph, icone ) ;
                    //	  if (isolation == -555 ) printf("peak icone %d m %d \n", icone, m);
                    
                }
                //12jul13	else if( (MggTrig > 0.08 &&  MggTrig < 0.11) ||  (MggTrig > 0.16 &&  MggTrig < 0.21) ) { 
                else if( (MggTrig > 0.075 &&  MggTrig < 0.115) ||  (MggTrig > 0.155 &&  MggTrig < 0.255) ) {
                    sprintf(key,"PP_PhiCor_side_%d",icone) ;
                    FillHistogram(key, delphi ,PtTrig) ;
                    
                    sprintf(key,"PP_r_BxE_side_%d",icone) ; 
                    if( PtTrig >= 6 && delphi> 0.66*TMath::Pi() && delphi<1.33*TMath::Pi() ) FillHistogram(key,track->Pt()*TMath::Cos(delphi)/PtTrig,PtTrig);
                    
                    //	  FillCorr(ph, icone+10 ) ;	  
                    //  if (isolation == -555 ) printf("side icone %d m2 %d \n", icone, m2 );
                    
                }
            }
        }
        
        Double_t delphineg = delphi + TMath::Pi() ;
        if(delphineg <= -TMath::PiOver2()) delphineg+=TMath::TwoPi();
        if(delphineg > 3*TMath::PiOver2()) delphineg-=TMath::TwoPi();
        
        Double_t delphiside = delphi + TMath::PiOver2() ;
        if(delphiside <= -TMath::PiOver2()) delphiside+=TMath::TwoPi();
        if(delphiside > 3*TMath::PiOver2()) delphiside-=TMath::TwoPi();
        
        
        Double_t dr    = TMath::Sqrt(deleta * deleta + delphi * delphi);
        Double_t drneg = TMath::Sqrt(deletaneg * deletaneg + delphineg * delphineg );
        //    Double_t drbg1 = TMath::Sqrt(deleta * deleta + delphineg * delphineg);
        // Double_t drbg2 = TMath::Sqrt(deletaneg * deletaneg + delphi * delphi);
        
        TSumPtcone += track->Pt(); 
        TVECSumPtcone += track->Pt()*TMath::Cos(delphi) ; 
        
        
        if ( isolation == -555 ){
            nstart = 9;
            nfin = 10;
            Int_t iicone = 10;
            if(  fabs(delphi) < TMath::TwoPi()/6 ) {
                FSumPtcone[iicone] +=  track->Pt(); 	 
                FNcone[iicone] += 1;
            }	
            else if(  fabs(delphi) > TMath::TwoPi()/6 &&  fabs(delphi) < TMath::TwoPi()/3 ) {  
                SSumPtcone[iicone] +=  track->Pt(); 	
                SNcone[iicone] += 1;
            }
            else if(  fabs(delphi) > TMath::TwoPi()/3 ) {
                BSumPtcone[iicone] +=   track->Pt(); 	
                // BxEcone[iicone] +=   fabs( track->Pt() *TMath::Cos(delphi) )/PtTrig ;
                BNcone[iicone] += 1; 
            }
        } 
        
        for(Int_t icone= nstart; icone< nfin; icone++){
            if( dr < Rcone[icone] ) {
                FSumPtcone[icone] +=  track->Pt(); 	 
                FNcone[icone] += 1;
            }
            if( drneg < 0.4 ) {
                BSumPtcone[icone] +=   track->Pt(); 	
                //	BxEcone[icone] +=   fabs( track->Pt()*TMath::Cos(delphi) )/PtTrig ;
                BNcone[icone] += 1; 
            } //opposide BG
            
            if( fabs(delphi) > TMath::TwoPi()/6 && fabs(delphi) < TMath::TwoPi()/3 ) { //was ( fabs(delphiside) < 0.4 ) {
                SSumPtcone[icone] +=  track->Pt(); 	
                SNcone[icone] += 1;
            } // side BG
            
            //            if (isolation == -555 ){
            //	      printf("  COnes Mtr %f ITYPE %d nst %d nfin %d \n", MggTrig, isolation, nstart, nfin);
            //	      printf("  co   FSumPtcone %f  BSumPtcone %f SSumPtcone %f  TSumPtcone %f icone %d Rcone %f \n",
            //		     FSumPtcone[icone], BSumPtcone[icone], SSumPtcone[icone], TSumPtcone, icone,  Rcone[icone] );  }
            
        }	
    }    // end of track loop 
    
    nstart = 0;
    nfin = nbin;
    
    if( isolation == -555) { 
        nstart = 9;
        nfin = 10;
    }
    
    
    for(Int_t icone= nstart; icone<nfin; icone++){
        Double_t i2 = icone ; 
        Double_t a2 = 2 ;  
        Double_t kind2 = pow( a2, i2) ;
        Int_t kind3 = float( kind2 );
        
        if(  ( isolation>=0 && (isolation&kind3) == kind3 ) || isolation== -555  ){
            
            //	printf("  Check2 Mtr %f ITYPE %d nst %d nfin %d \n", MggTrig, isolation, nstart, nfin);
            //	printf("  chk2   FSumPtcone %f  BSumPtcone %f SSumPtcone %f  TSumPtcone %f icone %d Rcone %f \n",
            //     FSumPtcone[icone], BSumPtcone[icone], SSumPtcone[icone], TSumPtcone, icone,  Rcone[icone] );
            
            
            if( MggTrig > 0.115 &&  MggTrig < 0.155 ) {
                //	if( MggTrig > 0.11 &&  MggTrig < 0.17 ) {
                // sprintf(key,"PP_R_FsumPt_peak_%d",icone) ;
                // FillHistogram(key, FSumPtcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_BsumPt_peak_%d",icone) ;
                // FillHistogram(key, BSumPtcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_SsumPt_peak_%d",icone) ;
                // FillHistogram(key, SSumPtcone[icone] ,PtTrig) ;
                
                sprintf(key,"PP_R_TsumPt_peak_%d",icone) ;
                FillHistogram(key, TSumPtcone, PtTrig) ;
                
                sprintf(key,"PP_R_TVECsumPt_peak_%d",icone) ;
                FillHistogram(key, TVECSumPtcone, PtTrig) ;
                
                //  sprintf(key,"PP_R_FMult_peak_%d",icone) ;
                // FillHistogram(key, FNcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_BMult_peak_%d",icone) ;
                // FillHistogram(key, BNcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_SMult_peak_%d",icone) ;
                // FillHistogram(key, SNcone[icone] ,PtTrig) ;
                
                sprintf(key,"PP_R_TMult_peak_%d",icone) ;
                FillHistogram(key, fMultiplicity ,PtTrig) ;	       
                
                sprintf(key,"PP_r_pi0_BPt_peak_%d",icone) ;
                if(  BSumPtcone[icone] > 0 && PtTrig > 0 )	FillHistogram(key, BSumPtcone[icone]/PtTrig, PtTrig  ) ;
                
                
                sprintf(key,"PP_r_pi0_SPt_peak_%d",icone) ;
                if(  SSumPtcone[icone] > 0 && PtTrig > 0 )	FillHistogram(key, SSumPtcone[icone]/PtTrig, PtTrig )  ;
                
                sprintf(key,"PP_r_pi0_TPt_peak_%d",icone) ;
                if(  TSumPtcone > 0  && PtTrig > 0 )	FillHistogram(key, TSumPtcone/PtTrig, PtTrig  ) ;
                
                //	  sprintf(key,"PP_pi0_T0multAC_peak_%d",icone) ;
                //	  FillHistogram(key, fT0multA, fT0multC   ) ;
                
                //	  sprintf(key,"PP_pi0_V0multAC_peak_%d",icone) ;
                //	  FillHistogram(key, fV0multA, fV0multC   ) ;
                
                
                //	  sprintf(key,"hT0AmC_peak_%d",icone) ;
                //	  FillHistogram(key, fT0multA-fT0multC );
                
                //	  sprintf(key,"hT0ApC_peak_%d",icone) ;
                //	  FillHistogram(key, fT0multA+fT0multC );
                
                //	  sprintf(key,"PP_pi0_SPDclustracl_peak_%d",icone) ;
                //	  FillHistogram(key, fSPDmultClust,  fSPDmultTracl  ) ;
                
                //	  sprintf(key,"PP_pi0_T0multAC_peak_%d",icone) ;
                //	  FillHistogram(key, fT0multA+ fT0multC,  fT0multA-fT0multC   ) ;
                
                
            } 
            
            
            
            else if( (MggTrig > 0.075 &&  MggTrig < 0.115) ||  (MggTrig > 0.155 &&  MggTrig < 0.255) ) {
                
                //	else if( (MggTrig > 0.08 &&  MggTrig < 0.11) ||  (MggTrig > 0.18 &&  MggTrig < 0.25 ) )  {
                // sprintf(key,"PP_R_FsumPt_side_%d",icone) ;
                //  FillHistogram(key, FSumPtcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_BsumPt_side_%d",icone) ;
                //  FillHistogram(key, BSumPtcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_SsumPt_side_%d",icone) ;
                //  FillHistogram(key, SSumPtcone[icone] ,PtTrig) ;
                
                sprintf(key,"PP_R_TsumPt_side_%d",icone) ;
                FillHistogram(key, TSumPtcone, PtTrig) ;
                
                sprintf(key,"PP_R_TVECsumPt_side_%d",icone) ;
                FillHistogram(key, TVECSumPtcone, PtTrig) ;
                
                // sprintf(key,"PP_R_FMult_side_%d",icone) ;
                // FillHistogram(key, FNcone[icone] ,PtTrig) ;
                
                //	  sprintf(key,"PP_R_BMult_side_%d",icone) ;
                // FillHistogram(key, BNcone[icone] ,PtTrig) ;
                
                // sprintf(key,"PP_R_SMult_side_%d",icone) ;
                // FillHistogram(key, SNcone[icone] ,PtTrig) ;
                
                sprintf(key,"PP_R_TMult_side_%d",icone) ;
                FillHistogram(key, fMultiplicity ,PtTrig) ;
                
                sprintf(key,"PP_r_pi0_BPt_side_%d",icone) ;
                if(  BSumPtcone[icone] > 0 && PtTrig > 0 )	FillHistogram(key, BSumPtcone[icone]/PtTrig, PtTrig  ) ;
                
                sprintf(key,"PP_r_pi0_SPt_side_%d",icone) ;
                if(  SSumPtcone[icone] > 0 && PtTrig > 0 )	FillHistogram(key, SSumPtcone[icone]/PtTrig, PtTrig  ) ;
                
                sprintf(key,"PP_r_pi0_TPt_side_%d",icone) ;
                if(  TSumPtcone > 0 && PtTrig > 0 )	FillHistogram(key, TSumPtcone/PtTrig, PtTrig  ) ;
                
                //	  sprintf(key,"PP_pi0_T0multAC_side_%d",icone) ;
                //	  FillHistogram(key, fT0multA+fT0multC,  fT0multA-fT0multC  ) ;
                
                //	  sprintf(key,"hT0AmC_side_%d",icone) ;
                //	  FillHistogram(key, fT0multA-fT0multC );
                
                //	  sprintf(key,"hT0ApC_side_%d",icone) ;
                //	  FillHistogram(key, fT0multA+fT0multC );
                
                //	  sprintf(key,"PP_pi0_V0multAC_side_%d",icone) ;
                //	  FillHistogram(key, fV0multA, fV0multC   ) ;
                
                //	  sprintf(key,"PP_pi0_SPDclustracl_side_%d",icone) ;
                //	  FillHistogram(key, fSPDmultClust,  fSPDmultTracl  ) ;
                
            }
        }
    }  // 2nd loop over the cones  
}

//_________________________________________________________________________________
Int_t AliAnalysisTaskSigma0::EvalIsolation(TLorentzVector * ph ){
    
    //  printf("                                          Check if this particle is isolated \n");
    //We use several cone radii and epsilons looking on charged particles and EMCAL ones
    
    Int_t nbin=9 ;
    Double_t Rcone[]={0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 , 0.8, 1.55 } ;
    Double_t SumPtcone[nbin];  Int_t IsolFlag[nbin]; 
    //  Double_t xecone[nbin];  Double_t multall[nbin];
    const Double_t kEpsilon=0.1 ;   // was 0.05 0.1
    //const Double_t kEpsilon2=0.05 ;
    //   const Double_t ptthresh=0.5 ;
    
    for(Int_t i=0; i<nbin; i++){  SumPtcone[i] = 0;  }
    
    Int_t isolation = 0;
    
    if(!ph) return 0 ;
    
    //Sum of energies in cones, tracks and clusters in EMCAL
    //   Double_t eCone1 = 0;   Double_t eCone2 = 0;   Double_t eCone3 = 0;   
    // Int_t iCone1 = 1;   Int_t iCone2 = 1;   Int_t iCone3 = 1;
    
    Double_t  phiTrig = ph->Phi();
    Double_t  etaTrig = ph->Eta();
    
    //   Double_t  MggTrig = ph->M();
    
    //   printf(" MggTrig %d \n",MggTrig );
    // Int_t ConeSum =0;   Int_t ConeProp = 0;   Int_t SumPt = 0;
    
    Int_t n=fTrackEvent->GetEntriesFast() ;
    
    for(Int_t itr=0; itr<n; itr++){
        AliCaloParticle * track = (AliCaloParticle*)fTrackEvent->At(itr) ;
        if(track==ph) //this particle
        continue ;     
        
        Double_t deleta = etaTrig - track->Eta() ;
        Double_t delphi = phiTrig - track->Phi() ;      
        //      while(delphi<-TMath::Pi()) delphi+=TMath::TwoPi() ; // while(delphi>TMath::Pi()) delphi-=TMath::TwoPi() ;
        if(delphi <= -TMath::PiOver2()) delphi+=TMath::TwoPi();
        if(delphi > 3*TMath::PiOver2()) delphi-=TMath::TwoPi();
        
        Double_t dr    = TMath::Sqrt(deleta * deleta + delphi * delphi);
        
        for(Int_t icone=0; icone<nbin; icone++){
            if( dr < Rcone[icone] ) {
                SumPtcone[icone] += track->Pt() ;
                //	 VSumPtcone[icone] += track->Pt() *TMath::Cos(delphi) ;	// Ncone[icone] += 1;
            }   // 
        }	// cones loop    
    }   // track loop
    isolation = 0; 
    for(Int_t icone=0; icone<nbin; icone++){
        IsolFlag[icone] = 0; 
        
        if( SumPtcone[icone] < kEpsilon* ph->Pt() )  IsolFlag[icone] = 1;
        
        //    if( SumPtcone[icone] <= 0.5 )  IsolFlag[icone] = 1;
        
        //    printf(" sumtr %f eps*pt %f isolfl %d \n", SumPtcone[icone] , kEpsilon* ph->Pt(), IsolFlag[icone] );    
        Double_t i2 = icone ; 
        Double_t a2 = 2 ;  
        Double_t kind2 = pow( a2, i2) ;
        Int_t kind3 = float( kind2 );
        isolation = isolation + int(  kind3*IsolFlag[icone] );
        //    printf("isolation %d  icone %d  \n", isolation,  icone );
        
    }
    return isolation ;		    
}

//_________________________________________________________________________________
Int_t AliAnalysisTaskSigma0::EvalMCIsolation(TParticle * ph){
    
    // Check if this particle is isolated. 
    //We use several cone radii and epsilons.
    //As well we look at charged particles and EMCAL ones
    
    const Double_t kConeR1=0.2 ;
    const Double_t kConeR2=0.3 ;
    const Double_t kConeR3=0.4 ;
    
    const Double_t kEpsilon1=0.1 ;
    const Double_t kEpsilon2=0.05 ;
    
    if(!ph) return 0 ;
    
    //   Int_t isolation = 15;
    
    //Sum of energies in cones, tracks and clusters in EMCAL
    Double_t eCone1 = 0;
    Double_t eCone2 = 0;
    Double_t eCone3 = 0;
    Double_t eCone1EMCAL = 0;
    Double_t eCone2EMCAL = 0;
    Double_t eCone3EMCAL = 0;
    
    Double_t  phiTrig = ph->Phi();
    Double_t  etaTrig = ph->Eta();
    
    //   Int_t n=fTrackEvent->GetEntriesFast() ;
    
    //   for(Int_t itr=0; itr<n; itr++){
    //     AliCaloParticle * track = (AliCaloParticle*)fTrackEvent->At(itr) ;
    
    //   printf(".MCisol pt-trig %f Y %f R %f PdG %d \n", ph->Pt(),ph->Y(), ph->R(), ph->GetPdgCode() );
    
    for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
        TParticle* particle = (TParticle *)fStack->Particle(iTracks) ;
        
        //    Int_t itr=iTracks; 
        //printf("...MCisol pt-trig %f Y %f R %f PdG %d itr %d \n"
        //    ,particle->Pt(),particle->Y(),particle->R(),particle->GetPdgCode(),iTracks);    
        
        if( particle  ==ph) //this particle
        continue ;     
        
        Int_t pdg = fabs(particle->GetPdgCode() ) ;
        
        // whi pi0 included - no reason -> skip it!
        // if( !(pdg==11 || pdg == 13 || pdg == 211 || pdg == 111 || pdg == 321 || pdg == 2212 ||  pdg == 22)) continue;     
        if( ! (pdg==11 || pdg == 13 || pdg == 211 ||  pdg == 321 || pdg == 2212 ||  pdg == 22  ) ) continue;    
        //     e+-       mu            pi+-           K+-         p                gamma
        
        //     Bool_t hitEMCAL = fEMCALgeom->Impact( particle ) ;             
        //     if( !hitEMCAL  ) continue ;
        //   printf(" hit %d pt-part %f pdg %d \n", hitEMCAL, particle->Pt(), pdg );
        
        if( TMath::Abs( particle->Eta() ) > 0.9 ) continue;
        
        Double_t deleta = etaTrig - particle->Eta() ;
        Double_t delphi = phiTrig - particle->Phi() ;      
        //      while(delphi<-TMath::Pi()) delphi+=TMath::TwoPi() ;
        //      while(delphi>TMath::Pi()) delphi-=TMath::TwoPi() ;
        if(delphi <= -TMath::PiOver2()) delphi+=TMath::TwoPi();
        if(delphi > 3*TMath::PiOver2()) delphi-=TMath::TwoPi();
        
        Double_t dr    = TMath::Sqrt(deleta * deleta + delphi * delphi);
        
        if ( pdg != 22  ) {
            if(dr<kConeR3){
                eCone3+= particle->Pt() ;
                if(dr<kConeR2){
                    eCone2+= particle->Pt() ;
                    if(dr<kConeR1){
                        eCone1+= particle->Pt() ;
                    }
                }
            }	
        }
        
        else if ( pdg == 22 ) {
            if(dr<kConeR3){
                eCone3EMCAL += particle->Pt() ;
                if(dr<kConeR2){
                    eCone2EMCAL += particle->Pt() ;
                    if(dr<kConeR1){
                        eCone1EMCAL += particle->Pt() ;
                    }
                }
            }	
        }
    }        
    
    //Fill QA histgams
    Double_t ptTrig=ph->Pt() ;
    
    Int_t iEMCone1E1 = 0;    Int_t iEMCone2E1 = 0;    Int_t iEMCone3E1 = 0;    
    Int_t iEMCone1E2 = 0;    Int_t iEMCone2E2 = 0;    Int_t iEMCone3E2 = 0;
    
    
    //Fill Bits
    Int_t iCone1E1 = (kEpsilon1*ptTrig > eCone1) ;
    Int_t iCone2E1 = (kEpsilon1*ptTrig > eCone2) ;
    Int_t iCone3E1 = (kEpsilon1*ptTrig > eCone3) ;
    
    Int_t iCone1E2 = (kEpsilon2*ptTrig > eCone1) ;
    Int_t iCone2E2 = (kEpsilon2*ptTrig > eCone2) ;
    Int_t iCone3E2 = (kEpsilon2*ptTrig > eCone3) ;
    
    Int_t mcisolation=   iCone1E1+  2*iCone2E1   +4*iCone3E1+
    8*iEMCone1E1+16*iEMCone2E1+32*iEMCone3E1+
    64*iCone1E2 +128*iCone2E2 +256*iCone3E2+
    512*iEMCone1E2+1024*iEMCone2E2+2048*iEMCone3E2;
    
    return mcisolation ;		    
}


//_____________________________________________________
void AliAnalysisTaskSigma0::FillCorr(TLorentzVector * trig, const Int_t itype  )
{
    //  Double_t xesum = 0; //  Double_t xesumpos = 0;   Double_t xesumneg = 0; 
    //  Double_t yesum = 0;
    
    Int_t n=fTrackEvent->GetEntriesFast() ;
    Int_t itype2 = itype*10;
    Double_t ptTrig = trig->Pt();
    
    if(  ptTrig <= 1 ) return;
    char key0[155];   char key2[155]; char key3[155];   char key5[55]; 
    char key1[155];  char key4[155]; 
    
    Int_t hnum = -1;
    
    if( (ptTrig > 5 && ptTrig <= 8  )   )  hnum  = itype2 ;
    else if( (ptTrig > 8 && ptTrig <= 10  )  )  hnum  = itype2 + 1 ;
    else if( (ptTrig >10 && ptTrig <= 16  )  )  hnum  = itype2 + 2 ;
    else if( (ptTrig >16 && ptTrig <= 50  )  )  hnum  = itype2 + 3 ;
    
    
    if( hnum == -1 ) return;
    
    
    sprintf(key5,"Multiplicity-%d",hnum);
    FillHistogram(key5,fMultiplicity);
    
    sprintf(key3,"hmgg-%d",hnum);
    FillHistogram(key3, trig->M() );      
    
    for(Int_t i=0; i<n;i++){
        
        AliCaloParticle * partn = static_cast<AliCaloParticle*>(fTrackEvent->At(i)) ;
        if(trig==partn) //for the case of track trigger
        continue ;
        
        Double_t dphi=trig->Phi()-partn->Phi() ;
        // while(dphi<-TMath::PiOver2())dphi+=TMath::TwoPi(); //while(dphi>3*TMath::PiOver2())dphi-=TMath::TwoPi() ;
        
        if(dphi <= -TMath::PiOver2()) dphi+=TMath::TwoPi();
        if(dphi > 3*TMath::PiOver2()) dphi-=TMath::TwoPi();
        
        // if( dphi > 2/3*TMath::Pi() &&  dphi < 5/3*TMath::Pi() ) {
        
        Double_t xe=-partn->Pt()*TMath::Cos(dphi)/trig->Pt() ;     
        Double_t pttr = partn->Pt();
        
        Double_t cosdelphi = TMath::Cos(dphi) ;
        if( cosdelphi == 0 ) cosdelphi = 1.e-6;
        //    Double_t rat= -pttr/ptTrig * cosdelphi / fabs( cosdelphi )  ;
        
        sprintf(key0,"hdelphi-%d",hnum); 
        FillHistogram(key0,dphi); 
        
        if ( dphi > 0.66*TMath::Pi() &&  dphi < 1.33*TMath::Pi() && pttr > 0.2   ){
            sprintf(key2,"hxe-%d",hnum); 
            FillHistogram(key2,xe); 
            
            sprintf(key1,"pTtrig-%d",hnum) ;
            FillHistogram(key1, ptTrig ) ;
            
            sprintf(key4,"phi-trk-%d",hnum) ;
            FillHistogram(key4,  pttr ) ;
            
        } // if of backward cone
        
    }  // loop over tracks
    
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

