#ifndef AliAnalysisTaskPPJetSpectra_cxx
#define AliAnalysisTaskPPJetSpectra_cxx

class THnSparse;
class TNtuple;
class TH2F;
class AliAODEvent;
class AliESDEvent;
class AliAODExtension;
class AliAODJet;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPPJetSpectra : public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskPPJetSpectra();
    AliAnalysisTaskPPJetSpectra(const char* name);
    virtual ~AliAnalysisTaskPPJetSpectra() {}

    virtual void        Init();
    virtual Bool_t      Notify();
    virtual void        UserCreateOutputObjects();
    virtual void        UserExec(Option_t *option);
    virtual void        Terminate(Option_t *);

    void                SetVertexCuts(Int_t, Double_t, Double_t);
    void                SetTrackFilter(UInt_t i) {nTrackFilter = i;}
    void                SetEventSelectionMask(UInt_t i) {fEvtSelectionMask = i;}
    void                SetEventClass(Float_t i) {fEventClass = i;}
    void                SetTrackCuts(Double_t,Double_t,Double_t);
    void                SetJetCuts(Double_t, Double_t, Double_t);

    void                SetRecJetBranch(TString);
    void                SetGenJetBranch(TString);
    void                SetRecBckgBranch(TString);
    void                SetGenBckgBranch(TString);
    void                SetTrackType(Int_t i) {fTrackType = i;}
    void                SetParticleType(Int_t i) {fParticleType = i;}
    void                SetNonStdFile(TString s) {fNonStdFile = s;}
    void				SetDoUEanalysis(Bool_t i) {kDoUEanalysis = i;}

    void                UseMC(Bool_t i) {fUseMC = i;}

    enum                {kNone = 0, kAOD, kAODMC, kAODMC2};
  private:

    Bool_t              EventSelection(Double_t*);
    Int_t               GetListOfTracks(Int_t, TList*);
    Int_t               GetListOfJets(TClonesArray*, TList*, Bool_t);
    void                FillJetContainer(TList*, THnSparseF*);
    Double_t            GetUE(TList*, TList*, Double_t,THnSparseF*);
    Double_t            GetBckgUE(TList*, Double_t,Bool_t, THnSparseF*);
    Int_t               CorrectForUE(TList*,Double_t,TList*,THnSparseF*);
    void                MatchJets(Bool_t,TList*,TList*,Float_t, THnSparseF*);
    Int_t	        CheckPtBin(Double_t);
    void		DoUEAnalysis(TList*, Double_t, Double_t);

    TList*              fOutputList;

    AliESDEvent*        fESD;
    AliAODEvent*        fAOD;
    AliAODEvent*        fAODIn;
    AliAODEvent*        fAODOut;
    AliAODExtension*    fAODExt;
    TString             fNonStdFile;
    Int_t               fDebug;
    Bool_t              fUseMC;

    UInt_t              fEvtSelectionMask;
    Float_t             fEventClass;
    Int_t               nVtxContCut;
    Double_t            fVtxZcut;
    Double_t            fVtxRcut;

    UInt_t              nTrackFilter;
    Double_t            trackPtMin;
    Double_t            trackPtMax;
    Double_t            trackEtaAbsMax;

    Double_t            jetPtMin;
    Double_t            jetEtaCut;
    Double_t            jetZmax;

    THnSparseF*         fhnEvent;

    THnSparseF*         fhnTracks;
    THnSparseF*         fhnMC;
    THnSparseF*         fhnMC2;

    THnSparseF*         fhnRecJetsNoCut;
    THnSparseF*         fhnGenJetsNoCut;
    THnSparseF*         fhnRecJetsCut;
    THnSparseF*         fhnGenJetsCut;
    THnSparseF*         fhnRecBckg;
    THnSparseF*         fhnGenBckg;

    THnSparseF*         fhnRecJetsTrackUEcor;
    THnSparseF*         fhnGenJetsTrackUEcor;
    THnSparseF*         fhnRecJetsBckgUEcor;
    THnSparseF*         fhnGenJetsBckgUEcor;

    THnSparseF*         fhnTrackUE;
    THnSparseF*         fhnParticleUE;
    THnSparseF*         fhnBckgRecUE;
    THnSparseF*         fhnBckgGenUE;

    TString             fRecJetBranch;
    TString             fGenJetBranch;
    TString             fRecBckgBranch;
    TString             fGenBckgBranch;

    Double_t            fRecJetR;
    Double_t            fGenJetR;
    Double_t            fRecBckgR;
    Double_t            fGenBckgR;

    Int_t               fTrackType;
    Int_t               fParticleType;

    THnSparseF*		fhnMatching;
    THnSparseF*		fhnTrackCorrMatching;
    THnSparseF*		fhnBckgCorrMatching;

    THnSparseF*		fhnTrackUEanal;

    Bool_t 		kDoUEanalysis;
    Int_t               fRejectPileUp;

    // dummies for assignment op and copy ctr
    AliAnalysisTaskPPJetSpectra(const AliAnalysisTaskPPJetSpectra&);
    AliAnalysisTaskPPJetSpectra &operator=(const AliAnalysisTaskPPJetSpectra&);

    ClassDef(AliAnalysisTaskPPJetSpectra, 2);
};

#endif
