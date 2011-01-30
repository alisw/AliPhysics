#ifndef ALIANALYSISTASKJETSPECTRUM2_H
#define ALIANALYSISTASKJETSPECTRUM2_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************

#include  "AliAnalysisTaskSE.h"
#include  "THnSparse.h" // cannot forward declare ThnSparseF

////////////////
class AliJetHeader;
class AliESDEvent;
class AliAODEvent;
class AliAODExtension;
class AliAODJet;
class AliAODJetEventBackground;
class AliGenPythiaEventHeader;
class AliCFManager;

class TList;
class TChain;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TSTring;


class AliAnalysisTaskJetSpectrum2 : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJetSpectrum2();
    AliAnalysisTaskJetSpectrum2(const char* name);
    virtual ~AliAnalysisTaskJetSpectrum2() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() { Init(); }
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual Bool_t Notify();

    virtual void SetUseGlobalSelection(Bool_t b){fUseGlobalSelection = b;}
    virtual void SetEventClass(Int_t i){fEventClass = i;}
    virtual void SetExternalWeight(Float_t f){fExternalWeight = f;}
    virtual void SetUseExternalWeightOnly(Bool_t b){fUseExternalWeightOnly = b;}
    virtual void SetAODJetInput(Bool_t b){fUseAODJetInput = b;}
    virtual void SetAODTrackInput(Bool_t b){fUseAODTrackInput = b;}
    virtual void SetAODMCInput(Bool_t b){fUseAODMCInput = b;}
    virtual void SetLimitGenJetEta(Bool_t b){fLimitGenJetEta = b;}
    virtual void SetBkgSubtraction(Bool_t b){fBkgSubtraction = b;}
    virtual void SetFillCorrBkg(Int_t i){fFillCorrBkg = i;}
    virtual void SetJetEtaWindow(Float_t f){fJetRecEtaWindow = f;}
    virtual void SetTrackEtaWindow(Float_t f){fTrackRecEtaWindow = f;}
    virtual void SetNMatchJets(Short_t f){fNMatchJets = f;}
    virtual void SetMinJetPt(Float_t f){fMinJetPt = f;}
    virtual void SetMinTrackPt(Float_t f){fMinTrackPt = f;}
    virtual void SetDeltaPhiWindow(Float_t f){fDeltaPhiWindow = f;}
    virtual void SetAnalysisType(Int_t i){fAnalysisType = i;}
    virtual void SetBranchGen(const char* c){fBranchGen = c;}
    virtual void SetBranchRec(const char* c){fBranchRec = c;}
    virtual void SetBranchBkg(const char* c){fBranchBkg = c;}  
    virtual void SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
    virtual void SetTrackTypeRec(Int_t i){fTrackTypeRec = i;}
    virtual void SetFilterMask(UInt_t i){fFilterMask = i;}
    virtual void SetEventSelectionMask(UInt_t i){fEventSelectionMask = i;}
    virtual void SetNonStdFile(char* c){fNonStdFile = c;} 


    // Helper
    //

    // we have different cases
    // AOD reading -> MC from AOD
    // ESD reading -> MC from Kinematics
    // this has to match with our selection of input events
    enum {kTrackUndef = 0, kTrackAOD, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};
    enum {kAnaMC =  0x1, kAnaMCESD = 0x2};
    enum {kMaxJets = 4};
    enum {kJetRec = 0, kJetGen, kJetTypes}; //
    enum {kMaxCorrelation =  3};
    
    // 
    // 0 all jets
    // 1 all jet in eta window
    // 2 all jets with partner
    // 3 all jets in eta window with partner
    // 4 all jets with partner in eta window
    enum {kStep0 = 0, kStep1, kStep2, kStep3, kStep4,kMaxStep};


 private:

    AliAnalysisTaskJetSpectrum2(const AliAnalysisTaskJetSpectrum2&);
    AliAnalysisTaskJetSpectrum2& operator=(const AliAnalysisTaskJetSpectrum2&);

    void    MakeJetContainer();
    Int_t   GetListOfTracks(TList *list,Int_t type);
    void    FillTrackHistos(TList &particlesList,int iType);

    Int_t   GetListOfJets(TList *list,TClonesArray* jarray,Int_t type);
    void    FillJetHistos(TList &jetsList,TList &particlesList,Int_t iType);

    void    FillMatchHistos(TList &recJetsList,TList &genJetsList);

    Bool_t  JetSelected(AliAODJet *jet);
    AliJetHeader *fJetHeaderRec;//! The jet header that can be fetched from the userinfo
    AliJetHeader *fJetHeaderGen;//! The jet header that can fetched from the userinfo
    AliAODEvent  *fAODIn; //! where we take the jets from 
    AliAODEvent  *fAODOut; //! where we take the jets from 
    AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
    THnSparseF   *fhnJetContainer[kMaxStep*2];   //! like particle container in corrfw with different steps need AliCFContainer with Scale(), and clone() to do the same
    THnSparseF   *fhnCorrelation;           //! response matrix for unfolding 
    THnSparseF   *fhnCorrelationPhiZRec;       //! response matrix for unfolding in max Z rec bins

    TF1          *f1PtScale;                //! correction function to correct to the average true jet energy depending on p_T,rec

    TString       fBranchRec;  // AOD branch name for reconstructed
    TString       fBranchGen;  // AOD brnach for genereated
    TString       fBranchBkg;  //AOD branch for background 
    TString       fNonStdFile; // name of delta aod file to catch the extension

    Bool_t        fUseAODJetInput;        // take jet from input AOD not from ouptu AOD
    Bool_t        fUseAODTrackInput;      // take track from input AOD not from ouptu AOD
    Bool_t        fUseAODMCInput;         // take MC from input AOD not from ouptu AOD
    Bool_t        fUseGlobalSelection;    // Limit the eta of the generated jets
    Bool_t        fUseExternalWeightOnly; // use only external weight
    Bool_t        fLimitGenJetEta;        // Limit the eta of the generated jets
    Bool_t        fBkgSubtraction;        // flag for bckg subtraction
    Short_t       fNMatchJets;            // number of leading jets considered from the list
    Int_t         fFillCorrBkg;           // flag for filling bckg response matrix
    UInt_t        fFilterMask;            // filter bit for slecected tracks
    UInt_t        fEventSelectionMask;    // Selection information used to filter events
    Int_t         fAnalysisType;          // Analysis type 
    Int_t         fTrackTypeRec;          // type of tracks used for FF 
    Int_t         fTrackTypeGen;          // type of tracks used for FF 
    Int_t         fEventClass;            // event class to be looked at for this instance of the task
    Float_t       fAvgTrials;             // Average nimber of trials
    Float_t       fExternalWeight;        // external weight
    Float_t       fJetRecEtaWindow;       // eta window for rec jets
    Float_t       fTrackRecEtaWindow;     // eta window for rec tracks
    Float_t       fMinJetPt;              // limits the jet p_T in addition to what already is done in the jet finder, this is important for jet matching for JF with lo threshold
    Float_t       fMinTrackPt;            // limits the track p_T 
    Float_t       fDeltaPhiWindow;        // minium angle between dijets


    TProfile*     fh1Xsec;   //! pythia cross section and trials
    TH1F*         fh1Trials; //! trials are added
    TH1F*         fh1PtHard;  //! Pt har of the event...       
    TH1F*         fh1PtHardNoW;  //! Pt har of the event without weigt       
    TH1F*         fh1PtHardTrials;  //! Number of trials 
    TH1F*         fh1ZVtx;          //! z-vtx distribution
    TH1F*         fh1NGenJets;      //! nr of gen jets
    TH1F*         fh1NRecJets;      //! nr of rec jets
    TH1F*         fh1PtTrackRec;    //! track pt
    TH1F*         fh1SumPtTrackRec; //! sum over all track pT    
    TH1F*         fh1SumPtTrackAreaRec; //! sum over all track pT    
    TH1F*         fh1TmpRho;            //! just temporary histo for calculation


    TH1F*         fh1PtRecIn[kMaxJets];  //! Jet pt for all this info is also in the THNsparse      
    TH1F*         fh1PtGenIn[kMaxJets];  //! Jet pt with corellated generated jet    

    TH1F*         fh1PtJetsRecIn;  //! Jet pt for all jets
    TH1F*         fh1PtJetsGenIn;  //! Jet pt for all jets
    TH1F*         fh1PtJetsLeadingRecIn;  //! Jet pt for all jets
    TH1F*         fh1PtTracksRecIn;  //! track pt for all tracks
    TH1F*         fh1PtTracksLeadingRecIn;  //! track pt for all tracks
    TH1F*         fh1PtTracksGenIn;  //! track pt for all tracks
    

    TH2F*         fh2NRecJetsPt;            //! Number of found jets above threshold
    TH2F*         fh2NRecTracksPt;          //! Number of found tracks above threshold
    TH2F*         fh2NGenJetsPt;            //! Number of found jets above threshold
    TH2F*         fh2NGenTracksPt;          //! Number of found tracks above threshold
    TH2F*         fh2PtFGen;                //! found vs generated 
    TH2F*         fh2RelPtFGen;             //! relative difference between generated and found 
    TH2F*         fh2RhoPtRec[kMaxJets];    //! jet shape variable rho
    TH2F*         fh2PsiPtRec[kMaxJets];    //! jet shape variable psi
    TH2F*         fh2RhoPtGen[kMaxJets];    //! 
    TH2F*         fh2PsiPtGen[kMaxJets];    //!
    TH2F*         fh2FragRec[kMaxJets];     //! fragmentation function
    TH2F*         fh2FragLnRec[kMaxJets];   //! fragmetation in xi
    TH2F*         fh2FragGen[kMaxJets];     //! fragmentation function
    TH2F*         fh2FragLnGen[kMaxJets];   //! fragmetation in xi



    // Jet histos second go

    TH1F*         fh1NJets[kJetTypes];      //! nr of gen jets
    TH1F*         fh1SumPtTrack[kJetTypes]; //! sum over all track pT    

    TH1F*         fh1PtIn[kJetTypes][kMaxJets];  //! Jet pt  
    TH1F*         fh1PtJetsIn[kJetTypes];  //! Jet pt for all jets
    TH1F*         fh1PtTracksIn[kJetTypes];  //! track pt for all tracks
    TH1F*         fh1PtTracksLeadingIn[kJetTypes];  //! track pt for all tracks
    
    TH2F*         fh2NJetsPt[kJetTypes];    //! Number of found jets above threshold
    TH2F*         fh2NTracksPt[kJetTypes];  //! Number of found jets above threshold
    TH2F*         fh2LeadingJetPtJetPhi[kJetTypes];     //! Phi distribution of accepted leading jets 
    TH2F*         fh2LeadingTrackPtTrackPhi[kJetTypes]; //! phi distribution of accepted leading tracks
    TH2F*         fh2RhoPt[kJetTypes][kMaxJets];     //! jet shape variable rho
    TH2F*         fh2PsiPt[kJetTypes][kMaxJets];     //! jet shape variable psi
    TH2F*         fh2PhiPt[kJetTypes][kMaxJets];       //! phi distribution correlation of jets      
    TH2F*         fh2EtaPt[kJetTypes][kMaxJets];       //! phi distribution correlation of jets      
    TH2F*         fh2PhiEta[kJetTypes][kMaxJets];      //! eta phi distribution of jet      



    TH1F*   fh1DijetMinv[kJetTypes];            //! dijet inv mass
    TH2F*   fh2DijetDeltaPhiPt[kJetTypes];      //! dijet delta phi vs pt
    TH2F*   fh2DijetAsymPt[kJetTypes];          //! dijet asym vs pt after delta phi cut
    TH2F*   fh2DijetPt2vsPt1[kJetTypes];        //! dijet pt2 vs pt1
    TH2F*   fh2DijetDifvsSum[kJetTypes];        //! dijet dif vs sum

    //background histos

    TH1F*         fh1Bkg1; //! background estimate, all jets
    TH1F*         fh1Bkg2; //! background estimate, wo 2 hardest jet 
    TH1F*         fh1Bkg3; //! background estimate, random jets
    TH1F*         fh1Sigma1; //! background fluctuations, all jets
    TH1F*         fh1Sigma2; //! background fluctuations, wo hardest jet
    TH1F*         fh1Sigma3; //! background fluctuations,random jets
    TH1F*         fh1Area1; //! average background jet area, all jets 
    TH1F*         fh1Area2; //! average background jet area, wo 2 hardest jet 
    TH1F*         fh1Area3; //! average background jet area, random jets
    TH1F*         fh1Ptjet; //! rec jet spectrum
    TH1F*         fh1Ptjetsub1;//! subtracted jet spectrum (Bkg1) 
    TH1F*         fh1Ptjetsub2; //! subtracted jet spectrum (Bkg2)
    TH1F*         fh1Ptjetsub3; //! subtracted jet spectrum (Bkg3)
    TH1F*         fh1Ptjethardest; //! rec hardest jet spectrum
    TH1F*         fh1Ptjetsubhardest1;//! subtracted hardest jet spectrum (Bkg1)
    TH1F*         fh1Ptjetsubhardest2;//! subtracted hardest jet spectrum (Bkg2)
    TH1F*         fh1Ptjetsubhardest3;//! subtracted hardest jet spectrum (Bkg3)
    TH2F*         fh2Rhovspthardest1;//! rho vs hardest subtracted jet pt (Bkg1)
    TH2F*         fh2Rhovspthardest2;//! rho vs hardest subtracted jet pt (Bkg2)
    TH2F*         fh2Rhovspthardest3;//! rho vs hardest subtracted jet pt (Bkg3)
    TH2F*         fh2Errorvspthardest1;//! relative error vs hardest subtracted jet pt (Bkg1)
    TH2F*         fh2Errorvspthardest2;//! relative error vs hardest subtracted jet pt (Bkg2)
    TH2F*         fh2Errorvspthardest3;//! relative error vs hardest subtracted jet pt (Bkg3)    
    TList *fHistList;                  //! Output list
   

    ClassDef(AliAnalysisTaskJetSpectrum2, 11) // Analysis task for standard jet analysis
};
 
#endif
