// TODO LIST
// TODO: You're all set!

#ifndef AliAnalysisTaskPhiCount_H
#define AliAnalysisTaskPhiCount_H

#include "AliAnalysisTaskSE.h"
#include "AliPPVsMultUtils.h"

class AliPIDResponse;
class AliAODTrack;

class AliAnalysisTaskPhiCount : public AliAnalysisTaskSE
{
    public:
    // Constructors
                                AliAnalysisTaskPhiCount     ( );
                                AliAnalysisTaskPhiCount     ( const char *      name );
        virtual                 ~AliAnalysisTaskPhiCount    ( );

    // Required implementations
        virtual void            UserCreateOutputObjects     ( );
        virtual void            UserExec                    ( Option_t*         option );
        virtual void            Terminate                   ( Option_t*         option );
    
    // Setters
        void                    fSetMCFlag                  ( Bool_t MCFlag )           { kMCbool = MCFlag; }
        void                    fSetPhiFlag                 ( Bool_t PhiFlag )          { kPhibool = PhiFlag; }
        void                    fSetKaonFlag                ( Bool_t KaonFlag )         { kKaonbool = KaonFlag; }
    
    private:
        void                    fSetZero                    ();
        void                    fPostData                   ( Bool_t            fEventEfficiency, Int_t nPhi, Int_t nPhiTru,Int_t nKaon, Int_t nKaonTru);
        bool                    fIsPrimaryVertexCandidate   ( AliAODEvent*      event );
        bool                    fIsTrackCandidate           ( AliAODTrack *     track );
        bool                    fIsKaonCandidate            ( AliAODTrack *     track );
        bool                    fSetKaonPID                 ( AliAODTrack *     track );
        bool                    fIsPhiCandidate             ( TLorentzVector    fPhi );
        bool                    fIsPhiGen                   ( AliAODMCParticle* particle );
        bool                    fIsPhiRec                   ( AliAODMCParticle* particle );
        bool                    fIsPhi                      ( AliAODMCParticle* particle );
        bool                    fIsKaonTruPhi               ( AliAODMCParticle* piKaon, AliAODMCParticle* pjKaon );
        void                    fFillPIDHist                ( AliAODTrack *     track, Int_t iIndex );
        void                    fFillVtxHist                ( Int_t iIndex );
        AliAODVertex*           fGetPrimaryVertex           ( AliAODEvent*      event ) { return fIsPrimaryVertexCandidate(event) ? fPrimaryVertex : nullptr; };
    
        AliAODEvent*            fAOD;               //! input event AOD Format
        AliESDEvent*            fESD;               //! input event ESD Format
        AliMCEvent*             fMCD;               //! input event MC

        TClonesArray*           AODMCTrackArray;    //! MC Tracks Array
        AliPPVsMultUtils*       fMultUtil;          //! ww
        
        Bool_t                  kMCbool;            // MC Flag
        Bool_t                  kPhibool;           // Phi tree Flag
        Bool_t                  kKaonbool;          // Kaon tree Flag
        
        // Tree
        TTree*                  fKaonCandidate;     //! output tree for Signal
        TTree*                  fPhiCandidate;      //! output tree for Signal
        TTree*                  fKaonEfficiency;    //! output tree for MC Truth
        TTree*                  fPhiEfficiency;     //! output tree for MC Truth

        // Event Variables
        Float_t                 fMultiplicity;      //! Event Multiplicity
        Float_t                 fMultiplicit2;      //! Event Multiplicity
        Float_t                 fMultiplicit3;      //! Event Multiplicity
        Int_t                   fKaonLabels [1024]; //! Kaon Labels
         
        // Tree Variables ( PhiCandidate )
        UChar_t                 fnPhi;              //! Number of Phis produced found
        Float_t                 fInvMass    [1024]; //! Invariant Mass
        Float_t                 fPhiPx      [1024]; //! Phi Px
        Float_t                 fPhiPy      [1024]; //! Phi Py
        Float_t                 fPhiPz      [1024]; //! Phi Pz
        UChar_t                 fiKaon      [1024]; //! iKaon
        UChar_t                 fjKaon      [1024]; //! jKaon
    
        // Tree Variables ( KaonCandidate )
        UChar_t                 fnKaon;             //! Number of Phis produced found
        Float_t                 fKaonPx     [1024]; //! Kaon Px
        Float_t                 fKaonPy     [1024]; //! Kaon Py
        Float_t                 fKaonPz     [1024]; //! Kaon Pz
        Char_t                  fCharge     [1024]; //! Kaon Charge
        Char_t                  fTOFSigma   [1024]; //! PID TOF
        Char_t                  fTPCSigma   [1024]; //! PID TPC
         
        // Tree Variables ( PhiEfficiency )
        UChar_t                 fnPhiTru;           //! Number of Phis produced found
        Float_t                 fPhiTruPx   [1024]; //! Phi Px
        Float_t                 fPhiTruPy   [1024]; //! Phi Py
        Float_t                 fPhiTruPz   [1024]; //! Phi Pz
        UChar_t                 fSelection  [1024]; //! Selection integer
        
        // List
        TList*                  fAnalysisOutputList;//! Analysis output list
        TList*                  fQCOutputList;      //! Analysis output list
        
        // List Variables
        TH1F*                   fHistEvntEff;       //! histogram of Event Efficiency
        TH1F*                   fHistVertex0;       //! histogram of Vertex
        TH2F*                   fHistTPCPID0;       //! histogram of TPC PID
        TH2F*                   fHistTOFPID0;       //! histogram of TOF PID
        TH1F*                   fHistVertex1;       //! histogram of Vertex
        TH2F*                   fHistTPCPID1;       //! histogram of TPC PID
        TH2F*                   fHistTOFPID1;       //! histogram of TOF PID
        TH2F*                   fHistTPCPID2;       //! histogram of TPC PID
        TH2F*                   fHistTOFPID2;       //! histogram of TOF PID
        TH2F*                   fHistTOFPID3;       //! histogram of TOF PID Sigma
        TH2F*                   fHistTPCPID3;       //! histogram of TPC PID Sigma

        AliAnalysisTaskPhiCount   (const AliAnalysisTaskPhiCount&); // not implemented
        AliAnalysisTaskPhiCount&  operator = (const AliAnalysisTaskPhiCount&); // not implemented
    
        AliPIDResponse*         fPIDResponse;       //! PID Response obj
    
        AliAODVertex*           fPrimaryVertex;     //! AOD Vertex obj
    
        ClassDef(AliAnalysisTaskPhiCount, 1);
};

#endif
