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
    void SetMaxY(Float_t maxy){fMaxY = maxy;}
    void SetMinPt(Float_t MinPt){fMinPt = MinPt;}
    void SetCollisionSystem(Int_t collisionSystem){fcollisionSystem = collisionSystem;}
    void SetWriteTTree(Bool_t WriteTTree){fWriteTTree = WriteTTree;}

    // For resolution smearing (from Theos LightFlavorGenerator)
    TObjArray       *fArr;
    TLorentzVector ApplyResolution(TLorentzVector vec);
    Double_t  PhiV(TLorentzVector l1, TLorentzVector l2);

  protected:
    AliVEvent*            fInputEvent;                // current event
    AliMCEvent*           fMCEvent;                   // corresponding MC event
    AliStack*             fMCStack;                   // stack belonging to MC event
    TList*                fOutputContainer;           // Output container
    
    Int_t*                fParticleList;              // array with particle Pdg values
    TString*              fParticleListNames;         // array with particle names
    
    // histograms events
    TH1F*                 fHistNEvents;               // number of events histo

    //TTree:
    TH1F** fmee;
    TH1F** fphi;
    TH1F** frap;
    TH2F** fpteevsmee;
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
    Float_t feep;
    Float_t feem;
    Float_t feeeta;
    Float_t feephi;
    Float_t feephiv;
    Float_t fmotherpt;
    Float_t fmotherp;
    Float_t fmotherm;
    Float_t fmothereta;
    Float_t fmotherphi;
    Int_t fID;
    Double_t fweight;
    Bool_t fpass;

    TString     fFileName;    // Name of the input file
    TFile*      fFile;        //! Pointer to input file
    //tree
    TTree*               teeTTree; 


    Int_t                 fIsMC;                      // MC flag
    Float_t              fMaxY;                      // Max y
    Float_t              fMinPt;
    Bool_t               fWriteTTree;
    Int_t              fcollisionSystem;
    
  private:
    AliAnalysisTaskLMeeCocktailMC(const AliAnalysisTaskLMeeCocktailMC&); // Prevent copy-construction
    AliAnalysisTaskLMeeCocktailMC &operator=(const AliAnalysisTaskLMeeCocktailMC&); // Prevent assignment

    ClassDef(AliAnalysisTaskLMeeCocktailMC, 1);
};

#endif
