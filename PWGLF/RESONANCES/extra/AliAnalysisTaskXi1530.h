#ifndef ALIANALYSISTASKXi1530_H
#define ALIANALYSISTASKXi1530_H
//
// Class AliAnalysisTaskXi1530
//
// AliAnalysisTaskXi1530
//  author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//        , Beomkyu  KIM (kimb@cern.ch)
//

#include <deque>
#include <THnSparse.h>
#include <AliAnalysisTaskSE.h>
class AliAnalysisTask;
class AliESDtrackCuts;
class AliESDEvent;
class AliAODEvent;
class AliStack;
class AliPIDResponse;
class AliPIDCombined;
class THistManager;

class AliAnalysisTaskXi1530RunTable {
public:
    enum {kPP,kPA,kAA,kUnknownCollType};
    AliAnalysisTaskXi1530RunTable();
    AliAnalysisTaskXi1530RunTable(Int_t runnumber);
    ~AliAnalysisTaskXi1530RunTable();
    
    Bool_t IsPP(){
        return fCollisionType==kPP;
    }
    Bool_t IsPA(){
        return fCollisionType==kPA;
    }
    Bool_t IsAA(){
        return fCollisionType==kAA;
    }
private:
    Int_t  fCollisionType=kPP; //! Is proton-proton collisions?
};

class AliAnalysisTaskXi1530 : public AliAnalysisTaskSE {
public:
    enum {kSD=0,kDD,kND,kCD,kAllProc,
        kXiStarCode          = 3324, // Xi(1530)^0 MC code
        kXiCode              = 3312, // Xi- MC code
        kLambdaCode          = 3122, // Lambda MC code
        kProtonCode          = 2212, // Proton+ MC code
        kPionCode            = 211}; // Pion+ MC code
    //PN = unlike sign, PP and NN are like signs
    
    AliAnalysisTaskXi1530();
    AliAnalysisTaskXi1530(const char *name, const char *option);
    AliAnalysisTaskXi1530(const AliAnalysisTaskXi1530& ap);
    AliAnalysisTaskXi1530& operator =(const AliAnalysisTaskXi1530& ap);
    ~AliAnalysisTaskXi1530();
    
    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t *);
    virtual void    Terminate(Option_t *);
    
    void SetOption   (char * option)    {fOption = option;}
    void SetFilterBit(UInt_t filterbit) {fFilterBit = filterbit;}
    void SetMixing   (Bool_t setmixing) {fsetmixing = setmixing;}
    void SetIsAA     (Bool_t isaa)      {IsAA = isaa;}
    void SetIsMC     (Bool_t ismc)      {IsMC = ismc;}
    void SetnMix     (Int_t nMix)       {fnMix = nMix;}
    void SetHighMult (Bool_t highmult)  {IsHighMult = highmult;}
    
    // Set Functions for the cut study & Systematic study
    void SetTPCNsigXi1530PionCut (Double_t nXi1530PionCut)  {fTPCNsigXi1530PionCut = nXi1530PionCut;}
    void SetTPCNsigLambdaProtonCut (Double_t nLambdaProtonCut)  {fTPCNsigLambdaProtonCut = nLambdaProtonCut;}
    void SetTPCNsigLambdaPionCut (Double_t nLambdaPionCut)  {fTPCNsigLambdaPionCut = nLambdaPionCut;}
    void SetTPCNsigBachelorPionCut (Double_t nBachelorPionCut)  {fTPCNsigBachelorPionCut = nBachelorPionCut;}
    void SetXi1530PionEtaCut (Double_t nXi1530PionEtaCut)  {fXi1530PionEtaCut = nXi1530PionEtaCut;}
    void SetXiEtaCut (Double_t nXiEtaCut)  {fXiEtaCut = nXiEtaCut;}
    void SetXi1530PionZVertexCut (Double_t nXi1530PionZVertexCut)  {fXi1530PionZVertexCut = nXi1530PionZVertexCut;}
    void SetDCADist_LambdaDaughtersCut (Double_t nDCADist_LambdaDaughtersCut)  {fDCADist_LambdaDaughtersCut = nDCADist_LambdaDaughtersCut;}
    void SetDCADist_XiDaughtersCut (Double_t nDCADist_XiDaughtersCut)  {fDCADist_XiDaughtersCut = nDCADist_XiDaughtersCut;}
    
    void SetDCADist_Lambda_PVCut (Double_t nDCADist_Lambda_PVCut)  {fDCADist_Lambda_PVCut = nDCADist_Lambda_PVCut;}
    
    void SetV0CosineOfPointingAngleCut (Double_t nV0CosineOfPointingAngleCut)  {fV0CosineOfPointingAngleCut = nV0CosineOfPointingAngleCut;}
    void SetCascadeCosineOfPointingAngleCut (Double_t nCascadeCosineOfPointingAngleCut)  {fCascadeCosineOfPointingAngleCut = nCascadeCosineOfPointingAngleCut;}
    void SetXiMassWindowCut (Double_t nXiMassWindowCut)  {fXiMassWindowCut = nXiMassWindowCut;}
    void SetXi1530RapidityCut (Double_t nXi1530RapidityCut)  {fXi1530RapidityCut = nXi1530RapidityCut;}
    
    
    
    Bool_t  GoodTracksSelection();
    Bool_t  GoodCascadeSelection();
    void FillTracks();
    
    Double_t GetMultiplicty(AliVEvent *fEvt);
    Bool_t SelectVertex2015pp(AliESDEvent *esd, Bool_t checkSPDres, Bool_t requireSPDandTrk, Bool_t checkProximity);
    Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
    Bool_t IsMCEventTrueINEL0();
    Bool_t IsTrueXi1530(AliESDcascade* Xi, AliVTrack* pion);
    void FillMCinput(AliStack* fMCStack);
    void FillTrackToEventPool();
    
    TAxis AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax);
    TAxis AxisVar( TString name, std::vector<Double_t> bin );
    TAxis AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0);
    TAxis AxisStr( TString name, std::vector<TString> bin );
    THnSparse * CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt="");
    THnSparse * CreateTHnSparse(TString name, TString title, TString templ, Option_t * opt="");
    Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w=1.);
    Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w=1.);
    
private:
    typedef std::vector< AliVTrack* > tracklist;
    typedef std::deque< tracklist >  eventpool;
    typedef std::vector< std::vector<eventpool> > mixingpool;
    
    TString                         fOption;
    
    AliESDtrackCuts*                fTrackCuts=nullptr; //!
    AliESDtrackCuts*                fTrackCuts2=nullptr; //!
    AliVEvent*                      fEvt=nullptr; //!
    UInt_t                          fFilterBit;
    AliAnalysisTaskXi1530RunTable*  fRunTable=nullptr; //!
    
    Double_t                        fCent=-1;
    Double_t                        ftrackmult=-1;
    Double_t                        fZ=-30;
    std::vector < UInt_t >          goodtrackindices; //!
    std::vector < UInt_t >          goodcascadeindices; //!
    
    AliPIDResponse                 *fPIDResponse=nullptr; //!
    AliPIDCombined                 *fPIDCombined=nullptr; //!
    //Histograms below are main
    
    
    mixingpool                      fEMpool; //!
    TAxis                           binCent; //!
    TAxis                           binZ; //!
    Int_t                           fnMix = 10;
    Int_t                           centbin = -1 ;
    Int_t                           zbin = -1 ;
    
    Double_t                        fTPCNsigXi1530PionCut   = 3.0;
    Double_t                        fTPCNsigLambdaProtonCut = 3.0;
    Double_t                        fTPCNsigLambdaPionCut   = 3.0;
    Double_t                        fTPCNsigBachelorPionCut = 3.0;
    Double_t                        fXi1530PionEtaCut       = 0.8;
    Double_t                        fXiEtaCut               = 0.8;
    Double_t                        fXi1530PionZVertexCut   = 2.0;
    Double_t                        fDCADist_LambdaDaughtersCut = 1.6;
    Double_t                        fDCADist_XiDaughtersCut     = 1.6;
    Double_t                        fDCADist_Lambda_PVCut       = 0.07;
    Double_t                        fV0CosineOfPointingAngleCut       = 0.97;
    Double_t                        fCascadeCosineOfPointingAngleCut       = 0.07;
    Double_t                        fXiMassWindowCut       = 0.007;
    Double_t                        fXi1530RapidityCut       = 0.5;
    
    Bool_t                          fsetmixing = kTRUE;
    Bool_t                          IsAA=kFALSE;
    Bool_t                          IsMC=kFALSE;
    Bool_t                          IsPS = kFALSE;
    Bool_t                          IsINEL0Rec = kFALSE;
    Bool_t                          IsINEL0True = kFALSE;
    Bool_t                          IsHighMult = kFALSE;
    THistManager*                   fHistos=nullptr; //!
    TClonesArray*                   fMCArray=nullptr; //!
    AliStack*                       fMCStack=nullptr; //!
    Int_t                           fNTracks = 0;
    Int_t                           fNCascade = 0;
    Double_t                        PVx = 999;
    Double_t                        PVy = 999;
    Double_t                        PVz = 999;
    Double_t                        bField = 999;
    ClassDef(AliAnalysisTaskXi1530, 5);
    //1: Frist version
    //2: Add Track cut2 for the Xi daughter particles
    //3: Add FillMixingPool function
    //4: Add Cut parameters to header and add "Set" fuction for cut study&Systematic study
    //5: include AliAnalysisTaskSE.h to avoid compile problem.
};

#endif
