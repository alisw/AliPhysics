#ifndef ALIANLYSISTASKLMEECOCKTAILMC_cxx
#define ALIANLYSISTASKLMEECOCKTAILMC_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

class AliAnalysisTaskLMeeCocktailMC : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskLMeeCocktailMC();
    AliAnalysisTaskLMeeCocktailMC(const char *name);
    virtual ~AliAnalysisTaskLMeeCocktailMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();

    // Configuration functions
    void SetCollisionSystem(Int_t collisionSystem){fcollisionSystem = collisionSystem;}
    void SetWriteTTree(Bool_t WriteTTree){fWriteTTree = WriteTTree;}
    void SetMaxEta(Float_t maxEta = 0.8){fMaxEta = maxEta;}
    void SetMinPt(Float_t MinPt = 0.2){fMinPt = MinPt;}
    void SetMaxPt(Float_t MaxPt = 8.0){fMaxPt = MaxPt;}
    void SetNBinsMee(Float_t val = 1200) {fNBinsMee = val;}
    void SetMinMee(Float_t val = 0.0) {fMinMee = val;}
    void SetMaxMee(Float_t val = 6.) {fMaxMee = val;}
    void SetNBinsPtee(Float_t val = 400) {fNBinsPtee = val;}
    void SetMinPtee(Float_t val = 0.0) {fMinPtee = val;}
    void SetMaxPtee(Float_t val = 10.0) {fMaxPtee = val;}
    void SetResolType(Int_t ResolType = 2){fResolType = ResolType;}
    void SetALTweight(Int_t ALTweightType = 1){fALTweightType = ALTweightType;}
    void SetResFileName(TString name);
    void SetResFileLocal(Bool_t localres) {fLocalRes = localres; }
    void SetEffFileName(TString name);
    void SetMinOpAng(Float_t MinOpAng = 0.050){fMinOpAng = MinOpAng;}
    void SetRapidityCut(Bool_t DoRapidityCut) {fDoRapidityCut = DoRapidityCut; }
    void SetMaxRapPair(Float_t maxRap = 0.8){fMaxRap = maxRap;}

    TH1F *GetEffHisto() {return fhwEffpT;}


    // For resolution smearing (from Theos LightFlavorGenerator)
    TObjArray       *fArr;
    TObjArray       *fArrResoPt;
    TObjArray       *fArrResoEta;
    TObjArray       *fArrResoPhi_Pos;
    TObjArray       *fArrResoPhi_Neg;
    //TLorentzVector ApplyResolution(TLorentzVector vec);
    TLorentzVector ApplyResolution(TLorentzVector vec, Char_t ech, Int_t Run);
    Double_t  PhiV(TLorentzVector l1, TLorentzVector l2);

  protected:
    AliVEvent*            fInputEvent;                // current event
    AliMCEvent*           fMCEvent;                   // corresponding MC event
    TList*                fOutputContainer;           // Output container

    Int_t*                fParticleList;              // array with particle Pdg values
    TString*              fParticleListNames;         // array with particle names

    const Int_t nInputParticles = 17;

    // Event histograms
    TH1F*                 fHistNEvents;               // number of events histo

    // Multiplicity weight histos (input):
    TH1F* fhwEffpT;
    TH1F* fhwMultpT;
    TH1F* fhwMultmT;
    TH1F* fhwMultpT2;
    TH1F* fhwMultmT2;

    // DCA input templates:
    TH1F** fh_DCAtemplates;
    const Int_t nbDCAtemplate = 6;

    //VPH histogram (pT) and function (mass)
    TF1* ffVPHpT;
    TH1F* fhKW;

    // output histograms
    // before smearing+acceptance cuts:
    TH1F** fmee_orig;
    TH2F** fpteevsmee_orig;
    TH1F** fmotherpT_orig;
    TH1F** fphi_orig;
    TH1F** frap_orig;
    TH1F** fmee_orig_wALT;
    TH2F** fpteevsmee_orig_wALT;
    TH1F** fmotherpT_orig_wALT;
    // after smearing + acceptance cuts
    TH1F** fmee;
    TH2F** fpteevsmee;
    TH1F** fphi;
    TH1F** frap;
    TH2F* fDCAeevsmee;
    TH2F* fDCAeevsptee;
    TH1F** fmee_wALT;
    TH2F** fpteevsmee_wALT;
    // LS, ULS histos
    TH2F* fULS_orig;
    TH2F* fLSpp_orig;
    TH2F* fLSmm_orig;
    TH2F* fULS;
    TH2F* fLSpp;
    TH2F* fLSmm;

    //TTree:
    Float_t fd1DCA;
    Float_t fd2DCA;
    Float_t fpairDCA;
    Float_t fd1origpt;
    Float_t fd1origp;
    Float_t fd1origeta;
    Float_t fd1origphi;
    Float_t fd2origpt;
    Float_t fd2origp;
    Float_t fd2origeta;
    Float_t fd2origphi;
    Float_t fd1pt;
    Float_t fd1p;
    Float_t fd1eta;
    Float_t fd1phi;
    Float_t fd2pt;
    Float_t fd2p;
    Float_t fd2eta;
    Float_t fd2phi;
    Float_t feeorigpt;
    Float_t feeorigp;
    Float_t feeorigm;
    Float_t feeorigeta;
    Float_t feeorigphi;
    Float_t feeorigphiv;
    Float_t feept;
    Float_t feemt;
    Float_t feep;
    Float_t feem;
    Float_t feeeta;
    Float_t feephi;
    Float_t feephiv;
    Float_t fmotherpt;
    Float_t fmothermt;
    Float_t fmotherp;
    Float_t fmotherm;
    Float_t fmothereta;
    Float_t fmotherphi;
    Int_t fID;
    UInt_t fdectyp;
    Int_t fdau3pdg;
    Double_t fweight;
    Double_t fwEffpT;
    Double_t fwMultpT;
    Double_t fwMultmT;
    Double_t fwMultpT2;
    Double_t fwMultmT2;
    Bool_t fpass;

    TString     fFileName;    // Name of the input file (resolution)
    TFile*      fFile;        //! Pointer to input file
    TString     fFileNameDCA;    // Name of the input file (DCA)
    TFile*      fFileDCA;        //! Pointer to input file
    TString     fFileNameEff;    // Name of the input file (Eff)
    TString     fFileNameEffLocal;    // Name of the input file (Eff)
    TFile*      fFileEff;        //! Pointer to input file
    TString     fFileNameWM;    // Name of the input file (weight multiplicity)
    TFile*      fFileWM;        //! Pointer to input file
    TString     fFileNameVPH;    // Name of the input file (VPH)
    TFile*      fFileVPH;        //! Pointer to input file
    TString     fResolDataSetName; //Specify multiplicity class and data set for Run 2 data
    Bool_t      fLocalRes;         // Local resolution file

    //tree
    TTree*               teeTTree;

    Int_t                 fIsMC;                      // MC flag
    Float_t              fMaxEta;                     // Max single electron Eta
    Float_t              fMinPt;                      // Min single electron Pt
    Float_t              fMaxPt;                      // Max single electron Pt
    Float_t              fNBinsMee;                   // number of bins in mass
    Float_t              fMinMee;                     // minimum mass
    Float_t              fMaxMee;                     // max mass
    Float_t              fNBinsPtee;                  // number of bins in pair pt
    Float_t              fMinPtee;                    // minimum pair pt
    Float_t              fMaxPtee;                    // minimum pair pt
    Float_t              fMinOpAng;                   // Min opening angle for pairs in mrad
    Bool_t               fWriteTTree;
    Int_t              fcollisionSystem;
    Int_t                fResolType;
    Int_t               fALTweightType;
    Bool_t               fDoRapidityCut;              // Possible to cut on the rapidity of the dielectron pair
    Float_t              fMaxRap;                     // Max rapidity cut

  private:
    AliAnalysisTaskLMeeCocktailMC(const AliAnalysisTaskLMeeCocktailMC&); // Prevent copy-construction
    AliAnalysisTaskLMeeCocktailMC &operator=(const AliAnalysisTaskLMeeCocktailMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskLMeeCocktailMC, 5);
};

#endif
