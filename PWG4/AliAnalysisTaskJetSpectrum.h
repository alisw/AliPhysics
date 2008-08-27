#ifndef ALIANALYSISTASKJETSPECTRUM_H
#define ALIANALYSISTASKJETSPECTRUM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
class AliJetFinder;
class AliESDEvent;
class AliAODEvent;
class AliAODJet;
class AliGenPythiaEventHeader;

class TList;
class TChain;
class TH2F;
class TH3F;

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

    virtual void SetExternalWeight(Float_t f){fExternalWeight = f;}
    virtual void SetUseExternalWeightOnly(Bool_t b){fUseExternalWeightOnly = b;}
    virtual void SetAODInput(Bool_t b){fUseAODInput = b;}
    virtual void SetAnalysisType(Int_t b){fAnalysisType = b;}
    virtual void SetBranchGen(char* c){fBranchGen = c;}
    virtual void SetBranchRec(char* c){fBranchRec = c;}

    // Helper
    static void GetClosestJets(AliAODJet *genJets,Int_t &nGenJets,
			       AliAODJet *recJets,Int_t &nRecJets,
			       Int_t *iGenIndex,Int_t *iRecIndex,Int_t iDebug = 0);

  //

    enum {kAnaMC =  0x1};

 private:

    AliAnalysisTaskJetSpectrum(const AliAnalysisTaskJetSpectrum&);
    AliAnalysisTaskJetSpectrum& operator=(const AliAnalysisTaskJetSpectrum&);
    

    enum {kMaxJets =  5};

    AliJetFinder *fJetFinderRec;
    AliJetFinder *fJetFinderGen;
    AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD

    TString       fBranchRec;  // AOD branch name for reconstructed
    TString       fConfigRec;  // Name of the Config file 
    TString       fBranchGen;  // AOD brnach for genereated
    TString       fConfigGen;  // Name of the Config file (if any)

    Bool_t        fUseAODInput;
    Bool_t        fUseExternalWeightOnly;
    Int_t         fAnalysisType;
    Float_t       fExternalWeight;    

    TH1F*         fh1PtHard;  // Pt har of the event...       
    TH1F*         fh1PtHard_NoW;  // Pt har of the event...       
    TH1F*         fh1PtHard_Trials;  // Number of trials 
    TH1F*         fh1PtHard_Trials_NoW;  // Number of trials 
    TH1F*         fh1NGenJets;
    TH1F*         fh1NRecJets;
    TH1F*         fh1E[kMaxJets];       // Jet Energy       
    TH1F*         fh1PtRecIn[kMaxJets];       // Jet pt for all      
    TH1F*         fh1PtRecOut[kMaxJets];      // Jet pt with corellated generated jet    
    TH1F*         fh1PtGenIn[kMaxJets];       // Detection efficiency for given p_T.gen
    TH1F*         fh1PtGenOut[kMaxJets];      //


    
    TH2F*         fh2PtFGen[kMaxJets];  //
    TH2F*         fh2Frag[kMaxJets];    // fragmentation function
    TH2F*         fh2FragLn[kMaxJets];  //

    TH3F*         fh3PtRecGenHard[kMaxJets];  //                              
    TH3F*         fh3PtRecGenHard_NoW[kMaxJets];  //                  
    TH3F*         fh3RecEtaPhiPt[kMaxJets]; // 
    TH3F*         fh3RecEtaPhiPt_NoGen[kMaxJets]; // 
    TH3F*         fh3RecEtaPhiPt_NoFound[kMaxJets]; //                    
    TH3F*         fh3MCEtaPhiPt[kMaxJets]; //                                                                                                        
    // ========= Multiplicity dependence ======

    // ==========TODO , flavaor dependence ========                           
    // ============================================                       
                                     

    // ============= TODO , phi dependence ========
    // ============================================                            
  
    TList *fHistList; // Output list


    ClassDef(AliAnalysisTaskJetSpectrum, 1) // Analysis task for standard jet analysis
};
 
#endif
