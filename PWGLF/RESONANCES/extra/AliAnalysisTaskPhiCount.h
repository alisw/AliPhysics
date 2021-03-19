// TODO LIST
// TODO: You're all set!

#ifndef AliAnalysisTaskPhiCount_H
#define AliAnalysisTaskPhiCount_H

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;
class AliAODTrack;

class AliAnalysisTaskPhiCount : public AliAnalysisTaskSE
{
    //
    public:
    //
    // Constructors
    //
                                AliAnalysisTaskPhiCount     ( );
                                AliAnalysisTaskPhiCount     ( const char *      name );
    virtual                    ~AliAnalysisTaskPhiCount     ( );
    //
    // Required implementations
    //
    virtual void                UserCreateOutputObjects     ( );
    virtual void                UserExec                    ( Option_t*         option );
    virtual void                Terminate                   ( Option_t*         option );
    //
    // Not implemented
    //
                                AliAnalysisTaskPhiCount     ( const AliAnalysisTaskPhiCount& );
    AliAnalysisTaskPhiCount&    operator =                  ( const AliAnalysisTaskPhiCount& );
    //
    // Setters & Getters
    //
    void                        SetMCFlag                   ( Bool_t    MCFlag )        { kMCbool = MCFlag; };
    void                        SetPhiFlag                  ( Bool_t    PhiFlag )       { kPhibool = PhiFlag; };
    void                        SetKaonFlag                 ( Bool_t    KaonFlag )      { kKaonbool = KaonFlag; };
    void                        SetFilterBit                ( Int_t     FilterBit )     { kFilterBit = FilterBit; };
    void                        SetVertexCut                ( Float_t   VertexCut )     { kVertexCut = VertexCut; };
    void                        SetkSgTPC_Alone             ( Float_t   kTPCSigma )     { kSgTPC_Alone = kTPCSigma; };
    void                        SetkSgTPC_TOFVt             ( Float_t   kTPCSigma )     { kSgTPC_TOFVt = kTPCSigma; };
    void                        SetkSgTOF_Veto              ( Float_t   kTOFSigma )     { kSgTOF_Veto = kTOFSigma; };
    void                        SetfRunName                 ( TString   kRunName )      { fRunName = kRunName; };
    Bool_t                      GetMCFlag                   ( )                         { return kMCbool; };
    Bool_t                      GetPhiFlag                  ( )                         { return kPhibool; };
    Bool_t                      GetKaonFlag                 ( )                         { return kKaonbool; };
    Int_t                       GetFilterBit                ( )                         { return kFilterBit; };
    Float_t                     GetVertexCut                ( )                         { return kVertexCut; };
    AliAODVertex               *fGetPrimaryVertex           ( )                         { return fPrimaryVertex; };
    //
    private:
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> General
    //
    void                        fPostData                   ( );
    void                        fFillTrees                  ( );
    void                        fSetZero                    ( );
    Double_t                   *fGetDCA                     ( );
    //
    Bool_t                      kMCbool;                    //  MC Flag
    Bool_t                      kPhibool;                   //  Phi tree Flag
    Bool_t                      kKaonbool;                  //  Kaon tree Flag
    Int_t                       kFilterBit;                 //  Filterbit
    Float_t                     kVertexCut;                 //  VertexCut
    TString                     fRunName;                   //  MultiRun name
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> QC & Selection
    //
    //>->   Event Selection
    //
    Bool_t                      fIsEventCandidate           ( );
    Bool_t                      fIsEventMultiplicityAvailable ( );
    Bool_t                      fIsEventPileUp              ( );
    void                        fSetEventMask               ( Int_t iMaskBit );
    void                        fSetTrueEventMask           ( Int_t iMaskBit );
    bool                        fCheckMask                  ( Int_t iMaskBit );
    bool                        fCheckTrueMask              ( Int_t iMaskBit );
    void                        fStoreTruePhi               ( Int_t iMaskBit );
    bool                        fCheckINELgt0               ( AliAODMCParticle* );
    void                        fFillEventEnumerate         ( Int_t iIndex );
    //
    AliAODEvent                *fAOD;                       //! input event AOD Format
    AliESDEvent                *fESD;                       //! input event ESD Format
    AliMCEvent                 *fMCD;                       //! input event MC
    TClonesArray               *AODMCTrackArray;            //! MC Tracks Array
    AliPIDResponse             *fPIDResponse;               //! PID Response obj
    AliAODVertex               *fPrimaryVertex;             //! AOD Vertex obj
    //
    //>->->     Vertex Selection
    //
    Bool_t                      fIsPrimaryVertexCandidate   ( );
    //
    //>->   Track Selection
    //
    Bool_t                      fIsTrackCandidate           ( );
    Bool_t                      fAssignTrack                ( );
    //
    AliAODTrack                *fCurrent_Track;             //! Track under scrutiny
    Int_t                       fCurrent_Track_Charge;      //! Track Charge
    Float_t                     fCurrent_Track_Momentum;    //! Track momentum
    Float_t                     fCurrent_Track_TransMom;    //! Track transverse momentum
    Float_t                     fCurrent_Track_Eta;         //! Track Eta
    Float_t                     fCurrent_Track_Phi;         //! Track Phi
    Float_t                     fCurrent_Track_DCAXY;       //! Track DCA xy
    Float_t                     fCurrent_Track_DCAZ;        //! Track DCA z
    //
    //>->->     Track QC
    //
    void                        fQC_TRK                     ( );
    void                        fQC_TRK_Kaons               ( );
    //
    //
    //>->   PID Selection
    //
    //
    //>->->     PID QC
    //
    void                        fQC_PID                     ( );
    void                        fQC_PID_Kaons               ( );
    void                        fQC_PID_Sel_Kaons           ( );
    Double_t                    fTOFBeta                    ( );
    //
    Bool_t                      fIsTPCAvailable;            //! TPC availabilty flag
    Bool_t                      fIsTOFAvailable;            //! TOF availabilty flag
    Float_t                     fBetaFromTOFSignal;         //! Particle beta from TOF signal
    Float_t                     fTPCSignal;                 //! Particle dE/dX in TPC
    Float_t                     kSgTPC_Alone;               // TPC Alone Sigma limit
    Float_t                     kSgTPC_TOFVt;               // TPC Alone Sigma limit
    Float_t                     kSgTOF_Veto;                // TPC Alone Sigma limit
    //
    //>->->     PID Kaons QC
    //
    Bool_t                      fIsKaonCandidate            ( );
    //
    //>->   QC Output
    //
    TList                      *fQCOutputList;              //! QC output list
    //
    //>->->     Event
    //
    TH1D                       *fQC_Event_Enumerate;        //! Analysis output list
    TH1F                       *fQC_Event_Vertex_Fll;       //! Analysis output list
    TH1F                       *fQC_Event_Vertex_Cut;       //! Analysis output list
    TH1F                       *fQC_Event_Enum_Mult;        //! Event vs Multiplicity
    //
    //>->->     Tracks
    //
    TH1F                       *fQC_Tracks_P_Momentum;      //! Analysis output list
    TH1F                       *fQC_Tracks_P_TransMom;      //! Analysis output list
    TH1F                       *fQC_Tracks_P_Eta;           //! Analysis output list
    TH1F                       *fQC_Tracks_P_Phi;           //! Analysis output list
    TH1F                       *fQC_Tracks_M_Momentum;      //! Analysis output list
    TH1F                       *fQC_Tracks_M_TransMom;      //! Analysis output list
    TH1F                       *fQC_Tracks_M_Eta;           //! Analysis output list
    TH1F                       *fQC_Tracks_M_Phi;           //! Analysis output list
    TH2F                       *fQC_Tracks_DCAXY_P;         //! Analysis output list
    TH2F                       *fQC_Tracks_DCAZ_P;          //! Analysis output list
    TH2F                       *fQC_Tracks_DCAXY_PT;        //! Analysis output list
    TH2F                       *fQC_Tracks_DCAZ_PT;         //! Analysis output list
    //
    //>->->->       Tracks: Kaons
    //
    TH1F                       *fQC_Kaons_P_Momentum;       //! Analysis output list
    TH1F                       *fQC_Kaons_P_TransMom;       //! Analysis output list
    TH1F                       *fQC_Kaons_P_Eta;            //! Analysis output list
    TH1F                       *fQC_Kaons_P_Phi;            //! Analysis output list
    TH1F                       *fQC_Kaons_M_Momentum;       //! Analysis output list
    TH1F                       *fQC_Kaons_M_TransMom;       //! Analysis output list
    TH1F                       *fQC_Kaons_M_Eta;            //! Analysis output list
    TH1F                       *fQC_Kaons_M_Phi;            //! Analysis output list
    TH2F                       *fQC_Kaon2_SigmaTPC_VETO_P;    //! gg
    TH2F                       *fQC_Kaon2_SigmaTPC_VETO_PT;    //! gg
    TH2F                       *fQC_Kaon2_SigmaTPC_P;    //! gg
    TH2F                       *fQC_Kaon2_SigmaTPC_PT;    //! gg
    TH2F                       *fQC_Kaon2_SigmaTOF_P;    //! gg
    TH2F                       *fQC_Kaon2_SigmaTOF_PT;    //! gg
    TH2F                       *fQC_Kaons_DCAXY_P;          //! Analysis output list
    TH2F                       *fQC_Kaons_DCAZ_P;           //! Analysis output list
    TH2F                       *fQC_Kaons_DCAXY_PT;         //! Analysis output list
    TH2F                       *fQC_Kaons_DCAZ_PT;          //! Analysis output list
    TH2F                       *fQC_Kaons_P_TPCSignal_P;      //! ee
    TH2F                       *fQC_Kaons_P_TOFSignal_P;      //! ee
    TH2F                       *fQC_Kaons_M_TPCSignal_P;      //! ee
    TH2F                       *fQC_Kaons_M_TOFSignal_P;      //! ee
    TH2F                       *fQC_Kaons_P_TPCSignal_PT;      //! ee
    TH2F                       *fQC_Kaons_P_TOFSignal_PT;      //! ee
    TH2F                       *fQC_Kaons_M_TPCSignal_PT;      //! ee
    TH2F                       *fQC_Kaons_M_TOFSignal_PT;      //! ee
    //
    //>->->     PID
    //
    TH2F                       *fQC_PID_SignalTPC_P;        //! Analysis output list
    TH2F                       *fQC_PID_SignalTOF_P;        //! Analysis output list
    TH2F                       *fQC_PID_SignalTPC_PT;       //! Analysis output list
    TH2F                       *fQC_PID_SignalTOF_PT;       //! Analysis output list
    //
    //>->->->       PID: Kaons
    //
    TH2F                       *fQC_Kaons_SigmaTPC_P;       //! Analysis output list
    TH2F                       *fQC_Kaons_SigmaTOF_P;       //! Analysis output list
    TH2F                       *fQC_Kaons_SigmaTPC_PT;      //! Analysis output list
    TH2F                       *fQC_Kaons_SigmaTOF_PT;      //! Analysis output list
    TH2F                       *fQC_Kaons_SigmaTOF_TPC;     //! Analysis output list
    TH2F                       *fQC_Kaon2_SigmaTOF_TPC;     //! edede
    TH2F                       *fQC_Kaons_SignalTPC_P;      //! Analysis output list
    TH2F                       *fQC_Kaons_SignalTOF_P;      //! Analysis output list
    TH2F                       *fQC_Kaons_SignalTPC_PT;     //! Analysis output list
    TH2F                       *fQC_Kaons_SignalTOF_PT;     //! Analysis output list
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> Output
    //
    // Event Variables
    //
    //>->   General Utilities
    Float_t                     fMultiplicity;              //! Event Multiplicity
    Int_t                       fCurrentRun;                //! Current Run Number
    Int_t                       fKaonLabels     [1024];     //! Kaon Labels
    Int_t                       fnPhiRec;                   //! Recordable Phi Number
    AliAODMCParticle*           fPhiRecParticles[1024];     //! Recordable Phi Labels
    //
    //>->->->   Data Event Mask
    UChar_t                     fEventMask;                 //! Event Mask
    Bool_t                      fIsINELgt0;                 //! Check the event is INEL > 0
    //
    //>->->->   MC Event Mask
    UChar_t                     fTrueEventMask;             //! True Event Mask
    Bool_t                      fIsTrueINELgt0;             //! Check the event is True INEL > 0
    //
    // Trees
    //
    TTree                      *fPhiCandidate;              //! output tree for Signal
    //
    Bool_t                      fIsPhiCandidate             ( TLorentzVector fPhi );
    //
    UChar_t                     fnPhi;                      //! Number of Phis produced found
    Float_t                     fInvMass        [1024];     //! Invariant Mass
    Float_t                     fPhiPx          [1024];     //! Phi Px
    Float_t                     fPhiPy          [1024];     //! Phi Py
    Float_t                     fPhiPz          [1024];     //! Phi Pz
    UChar_t                     fiKaon          [1024];     //! iKaon
    UChar_t                     fjKaon          [1024];     //! jKaon
    UChar_t                     fNature         [1024];     //! Nature
    //
    TTree                      *fKaonCandidate;             //! output tree for Signal
    UChar_t                     fnKaon;                     //! Number of Phis produced found
    Float_t                     fKaonPx         [1024];     //! Kaon Px
    Float_t                     fKaonPy         [1024];     //! Kaon Py
    Float_t                     fKaonPz         [1024];     //! Kaon Pz
    Char_t                      fCharge         [1024];     //! Kaon Charge
    Char_t                      fTOFSigma       [1024];     //! PID TOF
    Char_t                      fTPCSigma       [1024];     //! PID TPC
    //
    TTree                      *fPhiEfficiency;             //! output tree for MC Truth
    UChar_t                     fnPhiTru;                   //! Number of Phis produced found
    Float_t                     fPhiTruPx       [1024];     //! Phi Px
    Float_t                     fPhiTruPy       [1024];     //! Phi Py
    Float_t                     fPhiTruPz       [1024];     //! Phi Pz
    UChar_t                     fSelection      [1024];     //! Selection integer
    //
    TTree                      *fKaonEfficiency;            //! output tree for MC Truth
    UChar_t                     fnKaonTru;                  //! Number of Kaons produced found
    //
    // List
    TList                      *fAnalysisOutputList;        //! Analysis output list
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> TODO: TBD WIP
    //
    Bool_t                    fIsPhiGen                   ( AliAODMCParticle* particle );
    Bool_t                    fIsPhiRec                   ( AliAODMCParticle* particle );
    Bool_t                    fIsPhi                      ( AliAODMCParticle* particle );
    Bool_t                    fIsCandidateTruPhi          ( AliAODMCParticle* piKaon, AliAODMCParticle* pjKaon );
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> Class Definition
    //
    ClassDef(AliAnalysisTaskPhiCount, 1);
};

#endif
