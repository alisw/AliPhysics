#ifndef ALIANALYSISTASKJFSYSTEMATICS_H
#define ALIANALYSISTASKJFSYSTEMATICS_H
 
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



class AliAnalysisTaskJFSystematics : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJFSystematics();
    AliAnalysisTaskJFSystematics(const char* name);
    virtual ~AliAnalysisTaskJFSystematics() {;}
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
    virtual void SetAnalysisType(UInt_t i){fAnalysisType = i;}
    virtual void SetBranchGen(const char* c){fBranchGen = c;}
    virtual void SetBranchRec(const char* c){fBranchRec = c;}

    // ========= TODO Multiplicity dependence ======
    // ========= TODO z-dependence? ======
    // ========= TODO flavor dependence ========                           
    // ============================================                       
    enum {kSysJetOrder = 1, kSysTypes};
 
    enum {kMaxJets =  6}; // limit to 6 jets...

 private:

    AliAnalysisTaskJFSystematics(const AliAnalysisTaskJFSystematics&);
    AliAnalysisTaskJFSystematics& operator=(const AliAnalysisTaskJFSystematics&);


    static const Int_t fgkSysBins[kSysTypes];
    static const char* fgkSysName[kSysTypes];

    AliJetHeader *fJetHeaderRec;
    AliJetHeader *fJetHeaderGen;
    AliAODEvent  *fAOD; // where we take the jets from can be input or output AOD

    TString       fBranchRec;  // AOD branch name for reconstructed
    TString       fBranchGen;  // AOD brnach for genereated

    Bool_t        fUseAODInput;           // use AOD input
    Bool_t        fUseExternalWeightOnly; // use only external weight
    Bool_t        fLimitGenJetEta;        // Limit the eta of the generated jets
    UInt_t        fAnalysisType;          // Analysis type 
    Float_t       fExternalWeight;        // external weight
    Float_t       fRecEtaWindow;          // eta window used for corraltion plots between rec and gen 
    Float_t       fAvgTrials;             // average number of trials from pyxsec.root or pysec_hists.root in case trials are not avaiable from the MC Header
    // Event histograms
    TProfile*     fh1Xsec;    // pythia cross section and trials
    TH1F*         fh1Trials;  // trials are added
    TH1F*         fh1PtHard;  // Pt har of the even
    TH1F*         fh1PtHardNoW;  // Pt hard of the event without trials
    TH1F*         fh1PtHardTrials;  // Number of trials 
    TH1F*         fh1NGenJets;      // number of generated jets
    TH1F*         fh1NRecJets;      // number of reconstructed jets

    TH1F*         fh1PtRecIn;       // Jet pt for all      
    TH1F*         fh1PtRecOut;      // Jet pt with corellated generated jet    
    TH1F*         fh1PtGenIn;       // Detection efficiency for given p_T.gen
    TH1F*         fh1PtGenOut;      // gen pT of found jets

    TH2F*         fh2PtFGen;           // correlation betwen genreated and found  jet pT
    TH2F*         fh2PhiFGen;          // correlation betwen genreated and found  jet phi
    TH2F*         fh2EtaFGen;          // correlation betwen genreated and found  jet eta
    TH2F*         fh2PtGenDeltaPhi;   // difference between generated and found  jet phi
    TH2F*         fh2PtGenDeltaEta;   // difference between generated and found  jet eta


    TH3F*         fh3RecInEtaPhiPt;         // correlation between eta phi and rec pt                           
    TH3F*         fh3RecOutEtaPhiPt;        // correlation between eta phi and rec pt of jets with a partner   
    TH3F*         fh3GenInEtaPhiPt;         // correlation between eta phi and gen pt 
    TH3F*         fh3GenOutEtaPhiPt;        // correlation between eta phi and gen pt of jets with a partner       

    THnSparseF*     fhnCorrelation;          // correlation can be used for unfolding
                                     
  
    TList *fHistList; // Output list
    

    ClassDef(AliAnalysisTaskJFSystematics, 2) // Analysis task for standard jet analysis
};
 
#endif
