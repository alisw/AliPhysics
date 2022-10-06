// TODO LIST
// TODO: Add Invariant Mass to the Efficiency ( Migrate from canidate MC tree ? Add rec px, py, pz ? ) to check the Invariant Mass dependence of the Phi meson ( Rec (Gen IM) / Gen (Gen IM) ) not ( Rec (Rec IM) / Gen (Gen IM) )

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
    void                        SetMCFlag                   ( Bool_t    MCFlag )            { kMCbool = MCFlag; };
    void                        SetPhiFlag                  ( Bool_t    PhiFlag )           { kPhibool = PhiFlag; };
    void                        SetKaonFlag                 ( Bool_t    KaonFlag )          { kKaonbool = KaonFlag; };
    void                        SetSPCompute                ( Bool_t    SPFlag )            { kComputeSpherocity = SPFlag; };
    void                        SetSPWighted                ( Bool_t    SPWeightFlag )      { kSpherocityPTWeight = SPWeightFlag; };
    void                        SetRTCompute                ( Bool_t    RTFlag )            { kComputeRT = RTFlag; };
    void                        SetFilterBit                ( Int_t     FilterBit )         { kFilterBit = FilterBit; };
    void                        SetVertexCut                ( Float_t   VertexCut )         { kVertexCut = VertexCut; };
    void                        SetDCAzCut                  ( Float_t   DCAzCut )           { kDCAzCut = DCAzCut; };
    void                        SetMinTPCclusters           ( Float_t   MinTPCclusters )    { kMinTPCclusters = MinTPCclusters; };
    void                        SetChi2TPCcluster           ( Float_t   Chi2TPCcluster )    { kChi2TPCcluster = Chi2TPCcluster; };
    void                        SetChi2ITScluster           ( Float_t   Chi2ITScluster )    { kChi2ITScluster = Chi2ITScluster; };
    void                        SetNSigmaPtDepXYDCA         ( Float_t   NSigmaPtDepXYDCA )  { kNSigmaPtDepXYDCA = NSigmaPtDepXYDCA; };
    void                        SetChi2TPCGlobal            ( Float_t   Chi2TPCGlobal )     { kChi2TPCGlobal = Chi2TPCGlobal; };
    void                        SetTPCClsOverFndbl          ( Float_t   TPCClsOverFndbl )   { kTPCClsOverFndbl = TPCClsOverFndbl; };
    void                        SetkSgTPC_Alone             ( Float_t   TPCSigma )          { kSgTPC_Alone = TPCSigma; };
    void                        SetkSgTPC_TOFVt             ( Float_t   TPCSigma )          { kSgTPC_TOFVt = TPCSigma; };
    void                        SetkSgTOF_Veto              ( Float_t   TOFSigma )          { kSgTOF_Veto = TOFSigma; };
    void                        SetfRunName                 ( TString   RunName )           { fRunName = RunName; };
    void                        SetkTriggerMask             ( ULong64_t TriggerMask )       { kTriggerMask = TriggerMask; };
    void                        AddkTriggerMask             ( ULong64_t TriggerMask )       { kTriggerMask += TriggerMask; };
    void                        fSetEventCountLabels        ( TH1D *    fEvCount );
    Bool_t                      GetMCFlag                   ( )                             { return kMCbool; };
    Bool_t                      GetPhiFlag                  ( )                             { return kPhibool; };
    Bool_t                      GetKaonFlag                 ( )                             { return kKaonbool; };
    Int_t                       GetFilterBit                ( )                             { return kFilterBit; };
    Float_t                     GetVertexCut                ( )                             { return kVertexCut; };
    AliAODVertex               *fGetPrimaryVertex           ( )                             { return fPrimaryVertex; };
    //
    private:
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> General
    //
    void                        fPostData                   ( );
    void                        fFillTrees                  ( );
    void                        fSetZero                    ( );
    //
    Bool_t                      kMCbool;                    //  MC Flag
    Bool_t                      kComputeSpherocity;         //  Spherocity Flag
    Bool_t                      kSpherocityPTWeight;        //  Spherocity PT Weighted Flag
    Bool_t                      kComputeRT;                 //  RT Flag
    Bool_t                      kPhibool;                   //  Phi tree Flag
    Bool_t                      kKaonbool;                  //  Kaon tree Flag
    Float_t                     kVertexCut;                 //  VertexCut
    Int_t                       kFilterBit;                 //  Filterbit
    Float_t                     kDCAzCut;                   //  DCAz Cut
    Float_t                     kNSigmaPtDepXYDCA;          //  XY-DCA Sigmas Cut in pT
    Float_t                     kMinTPCclusters;            //  Minimum TPC Clusters
    Float_t                     kChi2TPCcluster;            //  Maximum Chi2 per TPC cluster
    Float_t                     kChi2TPCGlobal;             //  Maximum Chi2 per TPC Global
    Float_t                     kChi2ITScluster;            //  Maximum Chi2 per ITS cluster
    Float_t                     kTPCClsOverFndbl;           //  TPC Crossed Rows over findable
    TString                     fRunName;                   //  MultiRun name
    ULong64_t                   kTriggerMask;               //  TriggerMask
    //
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> QC & Selection
    //
    //>->   Event Selection
    //
    Bool_t                      fIsEventCandidate           ( );
    Bool_t                      fIsEventMultiplicityAvailable ( );
    Bool_t                      fIsEventPileUp              ( );
    bool                        fIsSPDClusterVsTrackletBG   ( );
    void                        fSetEventMask               ( Int_t iMaskBit );
    void                        fSetTrueEventMask           ( Int_t iMaskBit );
    bool                        fCheckMask                  ( Int_t iMaskBit );
    bool                        fCheckTrueMask              ( Int_t iMaskBit );
    void                        fStoreTruePhi               ( Int_t iMaskBit );
    bool                        fCheckINELgt0               ( AliAODMCParticle* );
    bool                        fCheckINELgt0               ( );
    void                        fFillEventEnumerate         ( Int_t iIndex );
    void                        fFillEventEnumerate         ( TString iBinName );
    void                        uCalculateSpherocity        ( );
    void                        uCalculateRT                ( );
    void                        uGuessCollisionSystem       ( );
    Float_t                     uTrackLengthInActiveTPC     ( AliExternalTrackParam* fCurrentTrackExternalParameters, Double_t deltaY, Double_t deltaZ );
    //
    AliAODEvent                *fAOD;                       //! input event AOD Format
    AliESDEvent                *fESD;                       //! input event ESD Format
    AliMCEvent                 *fMCD;                       //! input event MC
    TClonesArray               *AODMCTrackArray;            //! MC Tracks Array
    AliPIDResponse             *fPIDResponse;               //! PID Response obj
    AliAODVertex               *fPrimaryVertex;             //! AOD Vertex obj
    AliPPVsMultUtils           *fMultSelection;             //! Multiplicity Utility
    Bool_t                      fIs_p_p;                    //! Collision system flag for pp
    Bool_t                      fIs_p_Pb;                   //! Collision system flag for pPb
    Bool_t                      fIs_Pb_Pb;                  //! Collision system flag for PbPb
    //
    //>->->     Vertex Selection
    //
    Bool_t                      fIsPrimaryVertexCandidate   ( );
    //
    //>->   Track Selection
    //
    Bool_t                      fIsTrackCandidate           ( );
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
    Float_t                     kSgTPC_TOFVt;               // TPC TOF Veto Sigma limit
    Float_t                     kSgTOF_Veto;                // TOF Veto Sigma limit
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
    TH1D                       *fQC_Event_Enum_FLL;         //! Event Selection Counter
    TH1D                       *fQC_Event_Enum_V0M;         //! Event Multiplicity Counter w/ V0M
    TH1D                       *fQC_Event_Enum_TRK;         //! Event Multiplicity Counter w/ SPD Tracklets
    TH2D                       *fQC_Event_Enum_V0T;         //! Event Multiplicity Counter w/ SPD Tracklets
    TH1F                       *fQC_Event_Vertex_Fll;       //! Event Vertex Position
    TH1F                       *fQC_Event_Vertex_Cut;       //! Event Vertex Position w/ Cuts
    TH1F                       *fQC_Event_Spherocity;       //! Event Spherocity
    TH1F                       *fQC_Event_FullRT;           //! Event R Transverse
    TH1F                       *fQC_Event_FullRF;           //! Event R Forward
    TH1F                       *fQC_Event_FullRB;           //! Event R Backward
    //
    //>->->     Tracks
    //
    TH2F                       *fQC_Tracks_Momentum;        //! Acc. tracks Momentum
    TH2F                       *fQC_Tracks_TMomentum;       //! Acc. tracks Transverse Momentum
    TH2F                       *fQC_Tracks_Eta;             //! Acc. tracks Eta
    TH2F                       *fQC_Tracks_Phi;             //! Acc. tracks Phi
    TH2F                       *fQC_Tracks_V0M;             //! Acc. tracks vs V0M Multiplicity
    TH2F                       *fQC_Tracks_TRK;             //! Acc. tracks vs SPD Tracklets Multiplicity
    TH3F                       *fQC_Tracks_DCAXY_P;         //! Acc. tracks XY-DCA in Momentum
    TH3F                       *fQC_Tracks_DCAXY_PT;        //! Acc. tracks XY-DCA in Transverse Momentum
    TH3F                       *fQC_Tracks_DCAZ_P;          //! Acc. tracks Z-DCA in Momentum
    TH3F                       *fQC_Tracks_DCAZ_PT;         //! Acc. tracks Z-DCA in Transverse Momentum
    TH3F                       *fQC_Tracks_TOF_P;           //! Acc. tracks TOF Signal in Momentum
    TH3F                       *fQC_Tracks_TOF_PT;          //! Acc. tracks TOF Signal in Transverse Momentum
    TH3F                       *fQC_Tracks_TPC_P;           //! Acc. tracks TPC Signal in Momentum
    TH3F                       *fQC_Tracks_TPC_PT;          //! Acc. tracks TPC Signal in Transverse Momentum
    TH1F                       *fQC_Tracks_TPC_CLS;         //! Acc. tracks TPC Clusters
    TH1F                       *fQC_Tracks_TPC_FRC;         //! Acc. tracks TPC Clusters over findable
    TH1F                       *fQC_Tracks_TPC_CHI;         //! Acc. tracks TPC Chi^2 over Clusters
    TH1F                       *fQC_Tracks_TPC_CNS;         //! Acc. tracks TPC Global Constrained
    TH1F                       *fQC_Tracks_ITS_CHI;         //! Acc. tracks ITS Chi^2 over Clusters
    //
    //>->->->       Tracks: Kaons
    //
    //
    TH2F                       *fQC_Kaons_Momentum;         //! Kaons Momentum
    TH2F                       *fQC_Kaons_TMomentum;        //! Kaons Transverse Momentum
    TH2F                       *fQC_Kaons_Eta;              //! Kaons Eta
    TH2F                       *fQC_Kaons_Phi;              //! Kaons Phi
    TH2F                       *fQC_Kaons_V0M;              //! Kaons vs V0M Multiplicity
    TH2F                       *fQC_Kaons_TRK;              //! Kaons vs SPD Tracklets Multiplicity
    TH3F                       *fQC_Kaons_DCAXY_P;          //! Kaons XY-DCA in Momentum
    TH3F                       *fQC_Kaons_DCAXY_PT;         //! Kaons XY-DCA in Transverse Momentum
    TH3F                       *fQC_Kaons_DCAZ_P;           //! Kaons Z-DCA in Momentum
    TH3F                       *fQC_Kaons_DCAZ_PT;          //! Kaons Z-DCA in Transverse Momentum
    TH3F                       *fQC_Kaons_TOF_P;            //! Kaons TOF Signal in Momentum
    TH3F                       *fQC_Kaons_TOF_PT;           //! Kaons TOF Signal in Transverse Momentum
    TH3F                       *fQC_Kaons_TPC_P;            //! Kaons TPC Signal in Momentum
    TH3F                       *fQC_Kaons_TPC_PT;           //! Kaons TPC Signal in Transverse Momentum
    //
    //>->->     PID
    //
    TH3F                       *fQC_PID_TOF_Kaons_P;        //! n#sigma_{TOF}^{kaons} in Momentum
    TH3F                       *fQC_PID_TOF_Kaons_PT;       //! n#sigma_{TOF}^{kaons} in Transverse Momentum
    TH3F                       *fQC_PID_TPC_Kaons_P;        //! n#sigma_{TPC}^{kaons} in Momentum
    TH3F                       *fQC_PID_TPC_Kaons_PT;       //! n#sigma_{TPC}^{kaons} in Transverse Momentum
    TH3F                       *fQC_PID_TPC_TOF_Kaons_PT;   //! n#sigma_{TPC}^{kaons} vs n#sigma_{TOF}^{kaons}
    //
    //>->->->       PID: Kaons
    //
    TH3F                       *fQC_PID_TOF_NSig_SEL_Kaons_P;   //! Analysis output list
    TH3F                       *fQC_PID_TOF_NSig_SEL_Kaons_PT;  //! Analysis output list
    TH3F                       *fQC_PID_TPC_NSig_SEL_Kaons_P;   //! Analysis output list
    TH3F                       *fQC_PID_TPC_NSig_SEL_Kaons_PT;  //! Analysis output list
    TH3F                       *fQC_PID_TOF_Sgnl_SEL_Kaons_P;   //! Analysis output list
    TH3F                       *fQC_PID_TOF_Sgnl_SEL_Kaons_PT;  //! edede
    TH3F                       *fQC_PID_TPC_Sgnl_SEL_Kaons_P;   //! Analysis output list
    TH3F                       *fQC_PID_TPC_Sgnl_SEL_Kaons_PT;  //! Analysis output list
    //
    //>->->     GENERAL
    //
    TH2F                       *fQC_Phi_InvMass_Rec;            //! Test
    TH2F                       *fQC_Phi_InvMass_Gen;            //! Test
    TH2F                       *fQC_Phi_InvMass_Eff;            //! Test
    //
    
    //>->->->->->->->->->->->->->->->->->->->->->->->->->-> Output
    //
    // Event Variables
    //
    //>->   General Utilities
    Float_t                     fCurrent_SPH;               //! Event Spherocity
    Float_t                     fCurrent_V0M;               //! Event Multiplicity
    Float_t                     fCurrent_TRK;               //! Event Multiplicity
    Float_t                     fCurrent_RT;                //! Event RTransverse
    Int_t                       fCurrent_Run;               //! Current Run Number
    Int_t                       fKaonLabels     [2048];     //! Kaon Labels
    Int_t                       fnPhiRec;                   //! Recordable Phi Number
    AliAODMCParticle*           fPhiRecParticles[2048];     //! Recordable Phi Labels
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
    Int_t                       fnPhi;                      //! Number of Phis produced found
    Float_t                     fInvMass        [2048];     //! Invariant Mass
    Float_t                     fPhiPx          [2048];     //! Phi Px
    Float_t                     fPhiPy          [2048];     //! Phi Py
    Float_t                     fPhiPz          [2048];     //! Phi Pz
    Int_t                       fiKaon          [2048];     //! iKaon
    Int_t                       fjKaon          [2048];     //! jKaon
    Float_t                     fTrueInvMass    [2048];     //! True Invariant Mass
    //
    TTree                      *fKaonCandidate;             //! output tree for Signal
    Int_t                       fnKaon;                     //! Number of Phis produced found
    Float_t                     fKaonPx         [2048];     //! Kaon Px
    Float_t                     fKaonPy         [2048];     //! Kaon Py
    Float_t                     fKaonPz         [2048];     //! Kaon Pz
    Char_t                      fCharge         [2048];     //! Kaon Charge
    Char_t                      fTOFSigma       [2048];     //! PID TOF
    Char_t                      fTPCSigma       [2048];     //! PID TPC
    //
    TTree                      *fPhiEfficiency;             //! output tree for MC Truth
    Int_t                       fnPhiTru;                   //! Number of Phis produced found
    Float_t                     fPhiTruPx       [2048];     //! Phi Px
    Float_t                     fPhiTruPy       [2048];     //! Phi Py
    Float_t                     fPhiTruPz       [2048];     //! Phi Pz
    UChar_t                     fSelection      [2048];     //! Selection integer
    //
    TTree                      *fKaonEfficiency;            //! output tree for MC Truth
    Int_t                     fnKaonTru;                  //! Number of Kaons produced found
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
    enum    TRU_EvMask {
        kTRU_NOTRG      =   BIT(0),
        kTRU_INELGT0    =   BIT(1),
        kTRU_NOSPDVTX   =   BIT(2),
        kTRU_TRKSPDMM   =   BIT(3),
        kTRU_VTXCUT     =   BIT(4),
        kTRU_HAST10VTX  =   BIT(5),
    };
    enum    DAT_EvMask {
        kDAT_INELGT0    =   BIT(0)
    };
    //
    ClassDef(AliAnalysisTaskPhiCount, 1);
};

#endif
