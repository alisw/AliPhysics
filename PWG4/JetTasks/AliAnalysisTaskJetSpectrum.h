#ifndef ALIANALYSISTASKJETSPECTRUM_H
#define ALIANALYSISTASKJETSPECTRUM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
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
    virtual void SetAnalysisType(Int_t i){fAnalysisType = i;}
    virtual void SetBranchGen(char* c){fBranchGen = c;}
    virtual void SetBranchRec(char* c){fBranchRec = c;}

    // Helper
    static void GetClosestJets(AliAODJet *genJets,Int_t &nGenJets,
			       AliAODJet *recJets,Int_t &nRecJets,
			       Int_t *iGenIndex,Int_t *iRecIndex,Int_t iDebug = 0);

  //

    enum {kAnaMC =  0x1};
    enum {kMaxJets =  5};
    enum {kMaxCorrelation =  3};

 private:

    AliAnalysisTaskJetSpectrum(const AliAnalysisTaskJetSpectrum&);
    AliAnalysisTaskJetSpectrum& operator=(const AliAnalysisTaskJetSpectrum&);
    
    static const Float_t fgkJetNpartCut[kMaxCorrelation];


    AliJetHeader *fJetHeaderRec;
    AliJetHeader *fJetHeaderGen;
    AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD

    TString       fBranchRec;  // AOD branch name for reconstructed
    TString       fConfigRec;  // Name of the Config file 
    TString       fBranchGen;  // AOD brnach for genereated
    TString       fConfigGen;  // Name of the Config file (if any)

    Bool_t        fUseAODInput;
    Bool_t        fUseExternalWeightOnly;
    Bool_t        fLimitGenJetEta;
    Int_t         fAnalysisType;
    Float_t       fExternalWeight;    

    TProfile*     fh1Xsec;   // pythia cross section and trials
    TH1F*         fh1Trials; // trials are added
    TH1F*         fh1PtHard;  // Pt har of the event...       
    TH1F*         fh1PtHard_NoW;  // Pt har of the event...       
    TH1F*         fh1PtHard_Trials;  // Number of trials 
    TH1F*         fh1NGenJets;
    TH1F*         fh1NRecJets;
    TH1F*         fh1E[kMaxJets];       // Jet Energy       
    TH1F*         fh1PtRecIn[kMaxJets];       // Jet pt for all      
    TH1F*         fh1PtRecOut[kMaxJets];      // Jet pt with corellated generated jet    
    TH1F*         fh1PtGenIn[kMaxJets];       // Detection efficiency for given p_T.gen
    TH1F*         fh1PtGenOut[kMaxJets];      //


    
    TH2F*         fh2PtFGen[kMaxJets];  //
    TH2F*         fh2PhiFGen[kMaxJets];  //
    TH2F*         fh2EtaFGen[kMaxJets];  //
    TH2F*         fh2Frag[kMaxJets];    // fragmentation function
    TH2F*         fh2FragLn[kMaxJets];  //
    TH2F*         fh2PtGenDeltaPhi[kMaxJets];  
    TH2F*         fh2PtGenDeltaEta[kMaxJets];  

    TH3F*         fh3PtRecGenHard[kMaxJets];  //                              
    TH3F*         fh3PtRecGenHard_NoW[kMaxJets];  //                  
    TH3F*         fh3RecEtaPhiPt[kMaxJets]; // 
    TH3F*         fh3RecEtaPhiPt_NoGen[kMaxJets]; // 
    TH3F*         fh3GenEtaPhiPt_NoFound[kMaxJets]; //                    
    TH3F*         fh3GenEtaPhiPt[kMaxJets]; //                                                                                                        
    // ========= Multiplicity dependence ======

    // ==========TODO , flavaor dependence ========                           
    // ============================================                       
                                     

    // ============= TODO , phi dependence ========
    // ============================================                            
  
    TList *fHistList; // Output list
    
    ///////// For 2 dimensional unfolding //////////////////
    TH1F*         fh1JetMultiplicity;   
    TH2F*         fh2ERecZRec;
    TH2F*         fh2EGenZGen;
    TH2F*         fh2Efficiency;
    TH3F*         fh3EGenERecN;
    THnSparseF*   fhnCorrelation[kMaxCorrelation];
    ///////////////////////////////////////////////////////


    ClassDef(AliAnalysisTaskJetSpectrum, 1) // Analysis task for standard jet analysis
};
 
#endif
