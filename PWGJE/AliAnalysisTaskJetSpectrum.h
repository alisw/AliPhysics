#ifndef ALIANALYSISTASKJETSPECTRUM_H
#define ALIANALYSISTASKJETSPECTRUM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************

#include "AliAnalysisTaskSE.h"
#include  "THnSparse.h" // cannot forward declare ThnSparseF

////////////////
class AliJetHeader;
class AliESDEvent;
class AliAODEvent;
class AliAODJet;
class AliGenPythiaEventHeader;

class TList;
class TChain;
class TH2F;
class TH3F;
class TProfile;



class AliAnalysisTaskJetSpectrum : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetSpectrum();
    AliAnalysisTaskJetSpectrum(const char* name);
    virtual ~AliAnalysisTaskJetSpectrum() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual Bool_t Notify();


    virtual void SetExternalWeight(Float_t f){fExternalWeight = f;}
    virtual void SetUseExternalWeightOnly(Bool_t b){fUseExternalWeightOnly = b;}
    virtual void SetAODInput(Bool_t b){fUseAODInput = b;}
    virtual void SetLimitGenJetEta(Bool_t b){fLimitGenJetEta = b;}
    virtual void SetRecEtaWindow(Float_t f){fRecEtaWindow = f;}
    virtual void SetAnalysisType(Int_t i){fAnalysisType = i;}
    virtual void SetBranchGen(const char* c){fBranchGen = c;}
    virtual void SetBranchRec(const char* c){fBranchRec = c;}

    // Helper
    //

    enum {kAnaMC =  0x1};
    enum {kMaxJets =  4};
    enum {kMaxCorrelation =  3};

 private:

    AliAnalysisTaskJetSpectrum(const AliAnalysisTaskJetSpectrum&);
    AliAnalysisTaskJetSpectrum& operator=(const AliAnalysisTaskJetSpectrum&);
    
    static const Float_t fgkJetNpartCut[kMaxCorrelation]; // n Part cut


    AliJetHeader *fJetHeaderRec; // store jet header for rec
    AliJetHeader *fJetHeaderGen; // store jet header for rec
    AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD

    TString       fBranchRec;  // AOD branch name for reconstructed
    TString       fBranchGen;  // AOD brnach for genereated

    Bool_t        fUseAODInput;           // use AOD input
    Bool_t        fUseExternalWeightOnly; // use only external weight
    Bool_t        fLimitGenJetEta;        // Limit the eta of the generated jets
    Int_t         fAnalysisType;          // Analysis type 
    Float_t       fExternalWeight;        // external weight
    Float_t       fRecEtaWindow;          // eta window used for corraltion plots between rec and gen 

    TProfile*     fh1Xsec;   // pythia cross section and trials
    TH1F*         fh1Trials; // trials are added
    TH1F*         fh1PtHard;  // Pt har of the event...       
    TH1F*         fh1PtHardNoW;  // Pt har of the event...       
    TH1F*         fh1PtHardTrials;  // Number of trials 
    TH1F*         fh1NGenJets; // Number of gen jets
    TH1F*         fh1NRecJets; // Number of  rev jets
    TH1F*         fh1E[kMaxJets];       // Jet Energy       
    TH1F*         fh1PtRecIn[kMaxJets];       // Jet pt for all      
    TH1F*         fh1PtRecOut[kMaxJets];      // Jet pt with corellated generated jet    
    TH1F*         fh1PtGenIn[kMaxJets];       // Detection efficiency for given p_T.gen
    TH1F*         fh1PtGenOut[kMaxJets];      // gen pT of found jets

    TH2F*         fh2PtFGen[kMaxJets];   // correlation betwen genreated and found  jet pT
    TH2F*         fh2PhiFGen[kMaxJets];  // correlation betwen genreated and found  jet phi
    TH2F*         fh2EtaFGen[kMaxJets];  // correlation betwen genreated and found  jet eta
    TH2F*         fh2Frag[kMaxJets];    // fragmentation function
    TH2F*         fh2FragLn[kMaxJets];  // fragmetation in xi
    TH2F*         fh2PtGenDeltaPhi[kMaxJets];  // difference between generated and found  jet phi
    TH2F*         fh2PtGenDeltaEta[kMaxJets];  // difference between generated and found  jet eta
    TH2F*         fh2PtRecDeltaR[kMaxJets];     // distance of rec jet i to j > i p_T,j
    TH2F*         fh2PtGenDeltaR[kMaxJets];     // distance of jet i to j > i vs p_T,j

    TH3F*         fh3PtRecGenHard[kMaxJets];      // correlation beetwen pt hard, rec and gen                             
    TH3F*         fh3PtRecGenHardNoW[kMaxJets];  // correlation beetwen pt hard, rec and gen no weight                                              
    TH3F*         fh3RecEtaPhiPt[kMaxJets]; // correlation between eta phi and rec pt                                              
    TH3F*         fh3RecEtaPhiPtNoGen[kMaxJets]; // correlation between eta phi and rec pt of jets without a partner                    
    TH3F*         fh3GenEtaPhiPtNoFound[kMaxJets]; // correlation between eta phi and gen pt of jets without a partner
    TH3F*         fh3GenEtaPhiPt[kMaxJets]; // correlation between eta phi and gen pt of jets without a partner                              
    // ========= Multiplicity dependence ======

    // ==========TODO , flavaor dependence ========                           
    // ============================================                       
                                     

    // ============= TODO , phi dependence ========
    // ============================================                            
  
    TList *fHistList; // Output list
    
    ///////// For 2 dimensional unfolding //////////////////
    TH1F*         fh1JetMultiplicity;   // JetMultiplicity
    TH2F*         fh2ERecZRec;          // rec e and z
    TH2F*         fh2EGenZGen;          // gen e and z
    TH2F*         fh2Efficiency;        // 2 dimensional efficiency
    TH3F*         fh3EGenERecN;         // correlation rec energy, gen energy N particles 
    THnSparseF*   fhnCorrelation[kMaxCorrelation]; // 4d Correllation for unfolding
    ///////////////////////////////////////////////////////


    ClassDef(AliAnalysisTaskJetSpectrum, 2) // Analysis task for standard jet analysis
};
 
#endif
