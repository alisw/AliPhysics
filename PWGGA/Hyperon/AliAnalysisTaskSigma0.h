
#ifndef ALIANALYSISTASKHYPERON_H
#define ALIANALYSISTASKHYPERON_H



//#ifndef SIGMA0TEST_H
//#define SIGMA0TEST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//---------------------------------------------
// Class used to prepeare lists of photons in calorimeters,
// converted photon etc. and fill few QA histograms
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "AliV0ReaderV1.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "AliKFConversionPhoton.h"

class AliESDInputHandler;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class TList;
class TLorentzVector;
class TParticle ;

class AliCFContainer ;
class AliStack;
class AliESDpid ;
class AliESDtrackCuts ;
class AliEMCALGeometry ;
class AliPHOSGeoUtils ;
class AliExternalTrackParam ;
class AliKFParticle ;
class AliCaloParticle ;
class AliPHOSGeometry;
class AliTriggerAnalysis;
class AliMultiplicity;
class TH2I ;
class TTree;


#include "TH2I.h"

class AliAnalysisTaskSigma0 : public AliAnalysisTaskSE
{
    
public:
    AliAnalysisTaskSigma0();
    AliAnalysisTaskSigma0(const char* name);
    virtual ~AliAnalysisTaskSigma0() ;// virtual destructor
    
    // Implementation of interface methods
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t * /*option*/){}
    void SetPHOSBadMap(Int_t mod,TH2I * h)
    {
        if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
        fPHOSBadMap[mod]=new TH2I(*h) ;
    }
    
    //PID setters
    void SetDEdxCuts(Double_t sEUp=5., Double_t sEDn=-3., Double_t sPiUp=0., Double_t sPiDn=1.){
        fnSigmaAboveElectronLine= sEUp; fnSigmaBelowElectronLine=sEDn;
        fnSigmaAbovePionLine=sPiUp; fpnSigmaAbovePionLine=sPiDn;}
    void SetConvProbCut(Double_t prob=0.){ fprobCut=prob;}
    void SetConvMaxRCut(Double_t maxR=180.){fmaxR=maxR ;}
    void SetConvMaxZCut(Double_t maxZ=240.){fmaxZ=maxZ ;}
    void SetConvMaxEtaCut(Double_t eta=0.9){fetaCut=eta ;}
    void SetConvMinPtCut(Double_t minPt=0.02){fptCut=minPt ;}
    void SetConvMaxPtCut(Double_t maxPt=1.500){fptMaxCut=maxPt ;}
    
    void SetConvMaxChi2Cut(Double_t chi2=30.){fchi2CutConversion=chi2 ;}
    
    void GetArPod( Double_t pos[3], Double_t neg[3], Double_t moth[3],  Double_t arpod[2] );
    
    Double_t  Rapidity(Double_t pt, Double_t pz, Double_t m) ;
    
protected:
    void InitGeometry() ; //Create PHOS/EMCAL geometry
    void ProcessMC();
    
    void SelectTracks() ; //collects V0s in event
    void SelectConvPhotons() ; //collects V0s in event
    void SelectPhotonsFB() ;
    void SelectLambda() ;
    void SelectV0() ;
  void SelectGammaV0() ;
    void SelectMCLambda() ;
    void SelectPHOSPhotons(); //collects PHOS photons in event
    void SelectEMCALPhotons(); //collects EMCAL photons in event
    void FindLeaders(); //Mark leading particles
    
    void FillCorr(TLorentzVector * trig, const Int_t itype = 100) ;
     
    AliMCEvent   *fMCEvent;
    
    Bool_t IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz) ; //apply bad map
    
    void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
    
    
    
    Double_t PlanarityAngle(const AliExternalTrackParam * pos, const AliExternalTrackParam * neg)const ;
    void GetArmenterosQtAlfa(AliKFParticle* positiveKFParticle, AliKFParticle * negativeKFParticle,
                             AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] ) ;
    //Checks dispersion of PHOS clusters
    Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;
    
private:
    AliAnalysisTaskSigma0(const AliAnalysisTaskSigma0&); // Not implemented
    AliAnalysisTaskSigma0& operator=(const AliAnalysisTaskSigma0&); // Not implemented
    
    enum{
        kCaloPIDdisp = BIT(14),
        kCaloPIDtof  = BIT(15),
        kCaloPIDneutral= BIT(16)
    };
    enum{
        kConvOnFly= BIT(14),
        kConvKink = BIT(15),
        kConvdEdx = BIT(16),
        kConvProb = BIT(17),
        kConvR    = BIT(18),
        kConvZR   = BIT(19),
        kConvNDF  = BIT(20),
        kConvEta  = BIT(21),
        kConvPlan = BIT(22)
    };
    

    
protected:
  		
    
    //  AliV0ReaderV1 *fV0Reader;
    
    AliESDEvent* fESDEvent;    //!pointer to the ESDEvent
    AliESDpid * fESDpid ;      //class for Track PID calculation
  
    TTree *fTreeV0;
    AliPIDResponse  *   fPIDResponse;


    AliESDtrackCuts * fESDtrackCuts; //class for charged multiplicity estimation
    AliStack * fStack;         //! pointer to the MC particle stack
    TList * fOutputContainer;  //final histogram container
    TClonesArray *fReaderGammas;
    
    TList * fTreeList;
    //    TTree *fTreeV0;

    //    TTree *tSigma0;
    TTree *tTreeEvent;
    
    Double_t fCentr ;
    
   
    Int_t fRunPeriod ;
    
    Double_t fMinOpeningAngleGhostCut; // minimum angle cut
    
    AliPHOSGeoUtils  *fPHOSgeom;      //!PHOS geometry
    AliEMCALGeometry *fEMCALgeom;     //!EMCAL geometry
    Double_t fBadDistCutPHOS ; //Cut on distance to bad channel
    Double_t fBadDistCutEMCAL ; //Cut on distance to bad channel
    
    void SetCentrality(Double_t centr=0.){fCentrality=centr ; }   //Centrality of the event
    void SetMultiplicity(Double_t mult=0.){fMultiplicity=mult ; }
    void SetConvR(Double_t conv=0.){fConvR=conv ; }   //R of conversion
    
    
    //Containers for storing previous events
    // 10 bins for vtx class
    TList * fPHOSEvents[10] ;      //Container for PHOS photons
    TList * fEMCALEvents[10] ;     //Container for EMCAL photons
    TList * fConvEvents[10] ;      //Container for conversion photons
    //  TList * fLamEvents[10] ;      //Container for conversion photons
    TList * fGenpi0Events[10] ;      //Container for generated pi0
    
    
    //Current event
    TClonesArray * fTrackEvent ;   //tracks in the current event
    TClonesArray * fConvEvent ;    //Conversion photons in current event
    TClonesArray * fPHOSEvent ;    //PHOS photons in current event
    TClonesArray * fEMCALEvent ;   //EMCAL  photons in current event
    TClonesArray * fGenpi0Event ;   //EMCAL  photons in current event
    
    
    Double_t fnSigmaAboveElectronLine; //fnSigmaAboveElectronLine
    Double_t fnSigmaBelowElectronLine; //fnSigmaBelowElectronLine
    Double_t fnSigmaAbovePionLine;     //fnSigmaAbovePionLine
    Double_t fpnSigmaAbovePionLine;    //fpnSigmaAbovePionLine
    Double_t fprobCut;                 //fprobCut
    Double_t fmaxR ;                   //fmaxR
    Double_t fmaxZ ;                   //fmaxZ
    Double_t fetaCut ;                 //fetaCut
    Double_t fptCut ;                  //fptCut
    Double_t fptMaxCut ;
    Double_t fchi2CutConversion ;      //fchi2CutConversion
    Double_t fZvtx ;                   //Z vertex in currecnt event
    
    Double_t fCentrality ;    //Centrality of the event
    AliCentrality *fGetCent;
    Double_t fMultiplicity ;    //Centrality of the event
    Double_t fConvR ;         //Centrality of the event
    
    Double_t fPhi0Mass;
    Double_t fKaonMass;
    Double_t fRho0Mass;
    Double_t fPionMass;
    AliPHOSGeometry  *fPHOSGeo;  // PHOS geometry
    Int_t fEventCounter;         // number of analyzed events
    Double_t  fT0multA;
    Double_t  fT0multC;
    Double_t  fV0multA;
    Double_t  fV0multC;
    Double_t  fSPDmultClust;
    Double_t  fSPDmultTracl;
    Double_t fInputEvent;
    Float_t  fCentralityV0M;
    //  Double_t fMCEvent;
    Double_t fMCStack;
    Bool_t fIsMC;  
   
    
    
    //    AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation
    
    Int_t fnCINT1B;           // Number of CINT1B triggers
    Int_t fnCINT1A;           // Number of CINT1A triggers
    Int_t fnCINT1C;           // Number of CINT1C triggers
    Int_t fnCINT1E;           // Number of CINT1E triggers
    Double_t  fEtaCuts[15];
    
    
    
    
    //  Int_t fGammaV0s[1000] ;         //correspondence between final conv photon and V0
    //  Int_t fpipmV0s[1000] ;         //correspondence between pi+ pi- pair and V0
    // Int_t fGammaPHOS[1000] ;        //correspondence between final conv photon and V0
    // Int_t fGammaEMCAL[1000] ;       //correspondence between final conv photon and V0
    
    
    //================= Variables for Tree ==================//
    

    Float_t fLambdaMod;
     Float_t fLambdaMass;
     Float_t fLambdaPx;
     Float_t fLambdaPy;
     Float_t fLambdaPz;
     Float_t fLambdaArmPt;
     Float_t fLambdaArmAlpha;
     Float_t fLambdaEta;
     Float_t fLambdaCosPointingAngle;
     Float_t fLambdaDCAtoPVPos;
    Float_t fLambdaDCAtoPVNeg;
     Float_t fLambdaDCADaughters;
     Float_t fLambdaRadius;
    
    
     Float_t fGammaMass;
     Float_t fGammaPx;
     Float_t fGammaPy;
     Float_t fGammaPz;
    Float_t fGammaCosPointingAngle;
    Float_t fGammaDCADaughters;
  Float_t fGammaRadius;   
 Float_t fGammaDCAtoPVNeg;
    Float_t fGammaDCAtoPVPos;
    Float_t fGammaDCAzToPrimVtx;
     Float_t fGammaEta;
     Float_t fGammaArmPt;
     Float_t fGammaArmAlpha;
     Float_t fGammaZConv;
     Float_t fGammaChi2;
   
    
     Float_t fSigmaMass;
     Float_t fSigmaPt;
     Float_t fSigmaArmPt;
     Float_t fSigmaArmAlpha;
     //    Float_t  fCentralityV0M;

    Float_t   fLambdaTPx;
    Float_t   fLambdaTPy;
    Float_t   fLambdaTPz;
    Float_t   fGammaTPx;
    Float_t    fGammaTPy;
    Float_t    fGammaTPz;
    Float_t   fSigmaTPx;
    Float_t  fSigmaTPy;
    Float_t  fSigmaTPz;
 
    
    Float_t  fkSaveV0Tree;

 Float_t  fTreeVariableChi2V0;
    Float_t    fTreeVariableDcaV0Daughters;
    Float_t    fTreeVariableDcaV0ToPrimVertex;
    Float_t    fTreeVariableDcaPosToPrimVertex;
    Float_t    fTreeVariableDcaNegToPrimVertex;
    Float_t    fTreeVariableV0CosineOfPointingAngle;
    Float_t    fTreeVariableV0Radius;
    Float_t    fTreeVariablePt;
    Float_t     fTreeVariableRapK0Short;
    Float_t     fTreeVariableRapLambda;
     Float_t     fTreeVariableInvMassK0s;
     Float_t     fTreeVariableInvMassLambda;
     Float_t     fTreeVariableInvMassAntiLambda;
     Float_t     fTreeVariableAlphaV0;
     Float_t     fTreeVariablePtArmV0;
     Float_t     fTreeVariableNegEta;
     Float_t     fTreeVariablePosEta;
     Float_t     fTreeVariableNSigmasPosProton;
     Float_t     fTreeVariableNSigmasPosPion;
     Float_t      fTreeVariableNSigmasNegProton;
     Float_t       fTreeVariableNSigmasNegPion;
     Float_t      fTreeVariableDistOverTotMom;
     Float_t      fTreeVariableLeastNbrCrossedRows;
     Float_t     fTreeVariableLeastRatioCrossedRowsOverFindable;

 AliV0ReaderV1 *fV0Reader;

  
private:
    
    //  TH2I * fPHOSBadMap[6] ;   //Container for PHOS bad channels map
    TH2I * fEMCALBadMap[10] ; //Container for EMCAL Bad channels map  
    
    TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map
    
    ClassDef(AliAnalysisTaskSigma0, 3); // Analysis task for conversion + calorimeters
};


#endif //ALIANALYSISTASKHYPERON_H

// #endif //SIGMA0TEST_H
