//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 21 13:10:29 2012 by ROOT version 5.33/03
// from TTree fTree/V0Candidates
// found on file: DavidsV0MC_offline.root
//////////////////////////////////////////////////////////

#ifndef AliV0CutVariations_h
#define AliV0CutVariations_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class TH1F;
class TH2F;
class TH3F;
class TH1D;


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AliV0CutVariations : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           fTreeVariablePrimaryStatus;
   Int_t           fTreeVariablePrimaryStatusMother;
   Float_t         fTreeVariableChi2V0;
   Float_t         fTreeVariableDcaV0Daughters;
   Float_t         fTreeVariableDcaPosToPrimVertex;
   Float_t         fTreeVariableDcaNegToPrimVertex;
   Float_t         fTreeVariableV0Radius;
   Float_t         fTreeVariablePt;
   Float_t         fTreeVariablePtMC;
   Float_t         fTreeVariableRapK0Short;
   Float_t         fTreeVariableRapLambda;
   Float_t         fTreeVariableRapMC;
   Float_t         fTreeVariableInvMassK0s;
   Float_t         fTreeVariableInvMassLambda;
   Float_t         fTreeVariableInvMassAntiLambda;
   Float_t         fTreeVariableAlphaV0;
   Float_t         fTreeVariablePtArmV0;
   Float_t         fTreeVariableNegTransvMomentum;
   Float_t         fTreeVariablePosTransvMomentum;
   Float_t         fTreeVariableNegTransvMomentumMC;
   Float_t         fTreeVariablePosTransvMomentumMC;
   Int_t           fTreeVariableLeastNbrCrossedRows;
   Float_t         fTreeVariableLeastRatioCrossedRowsOverFindable;
   Int_t           fTreeVariablePID;
   Int_t           fTreeVariablePIDPositive;
   Int_t           fTreeVariablePIDNegative;
   Int_t           fTreeVariablePIDMother;
   Float_t         fTreeVariablePtXiMother;
   Float_t         fTreeVariableV0CosineOfPointingAngle;
   Int_t           fTreeVariableMultiplicity;
   Float_t         fTreeVariableDistOverTotMom;
   Float_t         fTreeVariableNSigmasPosProton;
   Float_t         fTreeVariableNSigmasPosPion;
   Float_t         fTreeVariableNSigmasNegProton;
   Float_t         fTreeVariableNSigmasNegPion;
   Float_t         fTreeVariableNegEta;
   Float_t         fTreeVariablePosEta;
   Float_t         fTreeVariableV0CreationRadius;
   Int_t           fTreeVariableIndexStatus;
   Int_t           fTreeVariableIndexStatusMother;
   Bool_t          fTreeVariableIsNonInjected;

   // List of branches
   TBranch        *b_fTreeVariablePrimaryStatus;   //!
   TBranch        *b_fTreeVariablePrimaryStatusMother;   //!
   TBranch        *b_Chi2V0;   //!
   TBranch        *b_fTreeVariableDcaV0Daughters;   //!
   TBranch        *b_fTreeVariableDcaPosToPrimVertex;   //!
   TBranch        *b_fTreeVariableDcaNegToPrimVertex;   //!
   TBranch        *b_fTreeVariableV0Radius;   //!
   TBranch        *b_fTreeVariablePt;   //!
   TBranch        *b_fTreeVariablePtMC;   //!
   TBranch        *b_fTreeVariableRapK0Short;   //!
   TBranch        *b_fTreeVariableRapLambda;   //!
   TBranch        *b_fTreeVariableRapMC;   //!
   TBranch        *b_fTreeVariableInvMassK0s;   //!
   TBranch        *b_fTreeVariableInvMassLambda;   //!
   TBranch        *b_fTreeVariableInvMassAntiLambda;   //!
   TBranch        *b_fTreeVariableAlphaV0;   //!
   TBranch        *b_fTreeVariablePtArmV0;   //!
   TBranch        *b_fTreeVariableNegTransvMomentum;   //!
   TBranch        *b_fTreeVariablePosTransvMomentum;   //!
   TBranch        *b_fTreeVariableNegTransvMomentumMC;   //!
   TBranch        *b_fTreeVariablePosTransvMomentumMC;   //!
   TBranch        *b_fTreeVariableLeastNbrCrossedRows;   //!
   TBranch        *b_fTreeVariableLeastRatioCrossedRowsOverFindable;   //!
   TBranch        *b_fTreeVariablePID;   //!
   TBranch        *b_fTreeVariablePIDPositive;   //!
   TBranch        *b_fTreeVariablePIDNegative;   //!
   TBranch        *b_fTreeVariablePIDMother;   //!
   TBranch        *b_fTreeVariablePtMother;   //!
   TBranch        *b_fTreeVariableV0CosineOfPointingAngle;   //!
   TBranch        *b_fTreeVariableMultiplicity;   //!
   TBranch        *b_fTreeVariableDistOverTotMom;   //!
   TBranch        *b_fTreeVariableNSigmasPosProton;   //!
   TBranch        *b_fTreeVariableNSigmasPosPion;   //!
   TBranch        *b_fTreeVariableNSigmasNegProton;   //!
   TBranch        *b_fTreeVariableNSigmasNegPion;   //!
   TBranch        *b_fTreeVariableNegEta;   //!
   TBranch        *b_fTreeVariablePosEta;   //!
   TBranch        *b_fTreeVariableV0CreationRadius;   //!
   TBranch        *b_fTreeVariableIndexStatus;   //!
   TBranch        *b_fTreeVariableIndexStatusMother;   //!
   TBranch        *b_fTreeVariableIsNonInjected;   //!

  AliV0CutVariations(TTree * /*tree*/ =0);
   virtual ~AliV0CutVariations() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;}
  void SetCentrality(Double_t min,Double_t max) {fCMin=min; fCMax=max;}
  void SetSelectNonInjected(Bool_t is=kTRUE) {fSelectNonInjected=is;}
  Bool_t AcceptV0();
  Bool_t AcceptTracks();
  Bool_t AcceptPID(Int_t code);

private:
  Bool_t fIsMC;         // MC flag
  Bool_t fSelectNonInjected;// non-injected flag
  Double_t fCMin;       // Min centrality
  Double_t fCMax;       // Max centrality
  Double_t fCPA;        // cos(PA) threshold
  Double_t fDCA;        // threshold for the DCA between V0 daughters
  Double_t fTPCcr;      // threshold for the number of crossed TPC pad rows
  Double_t fTPCcrfd;    // threshold for the ratio of TPC crossed/findable rows
  Double_t fDCApv;      // threshold for the DCA wrt the primary vertex

  //TList       *fOutput; //! The list of histograms

  TH1F *fMult;       //! Track multiplicity
  TH2F* fdEdx;       //! dEdx
  TH2F* fdEdxPid;    //! dEdx with PID

  TH1F *fCosPA;      //! cos(PA)
  TH1F *fDtrDCA;     //! DCA between V0 daughters
  TH1F *fTPCrows;    //! number of crossed TPC pad rows
  TH1F *fTPCratio;   //! ratio of TPC crossed/findable rows
  TH1F *fPrimDCA;    //! DCA wrt the primary vertex

  TH2F* fK0sM;       //! Mass for K0s
  TH2F* fK0sSi;      //! Side-band subtracted LvsP  for K0s 
  TH2F* fK0sMC;      //! LvsP for the K0s from the Monte Carlo stack 
  TH2F* fK0sAs;      //! LvsP for the K0s associated with the Monte Carlo 


  TH2F* fLambdaM;    //! Mass for Lambdas
  TH2F* fLambdaSi;   //! Side-band subtrated LvsP for Lambda
  TH2F* fLambdaMC;   //! LvsP for Lambdas from the Monte Carlo stack
  TH2F* fLambdaAs;   //! LvsP for Lambdas associated with the Monte Carlo

  TH1D* fLambdaEff;  //! Efficiency for Lambda  
  TH1D* fLambdaPt;   //! Pt spectrum for Lambda

  TH2F* fLambdaBarM;  //! Mass for anti-Lambdas
  TH2F* fLambdaBarSi; //! Side-band subtrated LvsP for anti-Lambda
  TH2F* fLambdaBarMC; //! LvsP for anti-Lambdas from the Monte Carlo stack
  TH2F* fLambdaBarAs; //! LvsP for anti-Lambdas associated with the Monte Carlo

  TH1D* fLambdaBarEff;  //! Efficiency for anti-Lambda  
  TH1D* fLambdaBarPt;   //! Pt spectrum for anti-Lambda

  TH3F* fLambdaFromXi;//! LvsPvsPxi for Lambdas from Xis associated with MC 
  TH2F* fXiM;         //! Mass for Xis
  TH1F* fXiSiP;       //! Side-band subtracted Pt for reconstructed Xi

  TH3F* fLambdaBarFromXiBar;//! LvsPvsPxi for anti-Lambdas from anti-Xis associated with MC 
  TH2F* fXiBarM;         //! Mass for anti-Xis
  TH1F* fXiBarSiP;       //! Side-band subtracted Pt for reconstructed anti-Xi

   ClassDef(AliV0CutVariations,1);
};

#endif

#ifdef AliV0CutVariations_cxx
void AliV0CutVariations::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   if (fIsMC) {
   fChain->SetBranchAddress("fTreeVariablePrimaryStatus", &fTreeVariablePrimaryStatus, &b_fTreeVariablePrimaryStatus);
   fChain->SetBranchAddress("fTreeVariablePrimaryStatusMother", &fTreeVariablePrimaryStatusMother, &b_fTreeVariablePrimaryStatusMother);
   fChain->SetBranchAddress("fTreeVariablePtMC", &fTreeVariablePtMC, &b_fTreeVariablePtMC);
   fChain->SetBranchAddress("fTreeVariableRapMC", &fTreeVariableRapMC, &b_fTreeVariableRapMC);
   fChain->SetBranchAddress("fTreeVariableNegTransvMomentum", &fTreeVariableNegTransvMomentum, &b_fTreeVariableNegTransvMomentum);
   fChain->SetBranchAddress("fTreeVariablePosTransvMomentum", &fTreeVariablePosTransvMomentum, &b_fTreeVariablePosTransvMomentum);
   fChain->SetBranchAddress("fTreeVariableNegTransvMomentumMC", &fTreeVariableNegTransvMomentumMC, &b_fTreeVariableNegTransvMomentumMC);
   fChain->SetBranchAddress("fTreeVariablePosTransvMomentumMC", &fTreeVariablePosTransvMomentumMC, &b_fTreeVariablePosTransvMomentumMC);
   fChain->SetBranchAddress("fTreeVariablePID", &fTreeVariablePID, &b_fTreeVariablePID);
   fChain->SetBranchAddress("fTreeVariablePIDPositive", &fTreeVariablePIDPositive, &b_fTreeVariablePIDPositive);
   fChain->SetBranchAddress("fTreeVariablePIDNegative", &fTreeVariablePIDNegative, &b_fTreeVariablePIDNegative);
   fChain->SetBranchAddress("fTreeVariablePIDMother", &fTreeVariablePIDMother, &b_fTreeVariablePIDMother);
   fChain->SetBranchAddress("fTreeVariablePtXiMother", &fTreeVariablePtXiMother, &b_fTreeVariablePtMother);
   fChain->SetBranchAddress("fTreeVariableV0CreationRadius", &fTreeVariableV0CreationRadius, &b_fTreeVariableV0CreationRadius);
   fChain->SetBranchAddress("fTreeVariableIndexStatus", &fTreeVariableIndexStatus, &b_fTreeVariableIndexStatus);
   fChain->SetBranchAddress("fTreeVariableIndexStatusMother", &fTreeVariableIndexStatusMother, &b_fTreeVariableIndexStatusMother);
   fChain->SetBranchAddress("fTreeVariableIsNonInjected", &fTreeVariableIsNonInjected, &b_fTreeVariableIsNonInjected);
   }

   fChain->SetBranchAddress("fTreeVariableChi2V0", &fTreeVariableChi2V0, &b_Chi2V0);
   fChain->SetBranchAddress("fTreeVariableDcaV0Daughters", &fTreeVariableDcaV0Daughters, &b_fTreeVariableDcaV0Daughters);
   fChain->SetBranchAddress("fTreeVariableDcaPosToPrimVertex", &fTreeVariableDcaPosToPrimVertex, &b_fTreeVariableDcaPosToPrimVertex);
   fChain->SetBranchAddress("fTreeVariableDcaNegToPrimVertex", &fTreeVariableDcaNegToPrimVertex, &b_fTreeVariableDcaNegToPrimVertex);
   fChain->SetBranchAddress("fTreeVariableV0Radius", &fTreeVariableV0Radius, &b_fTreeVariableV0Radius);
   fChain->SetBranchAddress("fTreeVariablePt", &fTreeVariablePt, &b_fTreeVariablePt);
   fChain->SetBranchAddress("fTreeVariableRapK0Short", &fTreeVariableRapK0Short, &b_fTreeVariableRapK0Short);
   fChain->SetBranchAddress("fTreeVariableRapLambda", &fTreeVariableRapLambda, &b_fTreeVariableRapLambda);
   fChain->SetBranchAddress("fTreeVariableInvMassK0s", &fTreeVariableInvMassK0s, &b_fTreeVariableInvMassK0s);
   fChain->SetBranchAddress("fTreeVariableInvMassLambda", &fTreeVariableInvMassLambda, &b_fTreeVariableInvMassLambda);
   fChain->SetBranchAddress("fTreeVariableInvMassAntiLambda", &fTreeVariableInvMassAntiLambda, &b_fTreeVariableInvMassAntiLambda);
   fChain->SetBranchAddress("fTreeVariableAlphaV0", &fTreeVariableAlphaV0, &b_fTreeVariableAlphaV0);
   fChain->SetBranchAddress("fTreeVariablePtArmV0", &fTreeVariablePtArmV0, &b_fTreeVariablePtArmV0);
   fChain->SetBranchAddress("fTreeVariableLeastNbrCrossedRows", &fTreeVariableLeastNbrCrossedRows, &b_fTreeVariableLeastNbrCrossedRows);
   fChain->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable", &fTreeVariableLeastRatioCrossedRowsOverFindable, &b_fTreeVariableLeastRatioCrossedRowsOverFindable);
   fChain->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle", &fTreeVariableV0CosineOfPointingAngle, &b_fTreeVariableV0CosineOfPointingAngle);
   fChain->SetBranchAddress("fTreeVariableMultiplicity", &fTreeVariableMultiplicity, &b_fTreeVariableMultiplicity);
   fChain->SetBranchAddress("fTreeVariableDistOverTotMom", &fTreeVariableDistOverTotMom, &b_fTreeVariableDistOverTotMom);
   fChain->SetBranchAddress("fTreeVariableNSigmasPosProton", &fTreeVariableNSigmasPosProton, &b_fTreeVariableNSigmasPosProton);
   fChain->SetBranchAddress("fTreeVariableNSigmasPosPion", &fTreeVariableNSigmasPosPion, &b_fTreeVariableNSigmasPosPion);
   fChain->SetBranchAddress("fTreeVariableNSigmasNegProton", &fTreeVariableNSigmasNegProton, &b_fTreeVariableNSigmasNegProton);
   fChain->SetBranchAddress("fTreeVariableNSigmasNegPion", &fTreeVariableNSigmasNegPion, &b_fTreeVariableNSigmasNegPion);
   fChain->SetBranchAddress("fTreeVariableNegEta", &fTreeVariableNegEta, &b_fTreeVariableNegEta);
   fChain->SetBranchAddress("fTreeVariablePosEta", &fTreeVariablePosEta, &b_fTreeVariablePosEta);
}

Bool_t AliV0CutVariations::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef AliV0CutVariations_cxx
