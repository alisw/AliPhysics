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
class THnSparse;
class AliAnalysisTask;
class AliAnalysisTaskSE;
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
    Double_t                        fptcut = 0.15;
    Double_t                        fetacut = 0.8;
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
    ClassDef(AliAnalysisTaskXi1530, 3);
    //1: Frist version
    //2: Add Track cut2 for the Xi daughter particles
    //3: Add FillMixingPool function
};

#endif
