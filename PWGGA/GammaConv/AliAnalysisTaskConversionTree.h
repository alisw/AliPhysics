#ifndef AliAnalysisConversionTree_cxx
#define AliAnalysisConversionTree_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include <vector>
#include "AliV0ReaderV1.h"
#include "AliConversionPhotonCuts.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskConversionTree : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskConversionTree();
    AliAnalysisTaskConversionTree(const char *name);
    virtual ~AliAnalysisTaskConversionTree();

    virtual void   UserCreateOutputObjects  ();
    virtual Bool_t Notify                   ();
    virtual void   UserExec                 ( Option_t *option );
    virtual void   Terminate                ( Option_t * );

    void SetV0Reader                        ( AliV0ReaderV1 *v0Reader )                   { fV0Reader=v0Reader                  ; }
    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetConversionCuts                  ( AliConversionPhotonCuts* conversionCuts, 
                                              Bool_t IsHeavyIon )                         {
                                                                                            fConversionCuts=conversionCuts      ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }
    void SetEventCuts                       ( AliConvEventCuts* conversionCuts,
                                              Bool_t IsHeavyIon )                         {
                                                                                            fEventCuts=conversionCuts           ;
                                                                                            fIsHeavyIon = IsHeavyIon            ;
                                                                                          }

    void FillType                           ( Double_t fillTree,
                                              Bool_t fillHistorams)                       {
                                                                                            ffillTree = fillTree                ;
                                                                                          }
    void SetIsMC                            ( Bool_t isMC )                               { fIsMC = isMC                        ; }

  private:
        
    AliAnalysisTaskConversionTree     ( const AliAnalysisTaskConversionTree& ); // Prevent copy-construction
    AliAnalysisTaskConversionTree &operator=( const AliAnalysisTaskConversionTree& ); // Prevent assignment

    void ProcessQATree              ( AliAODConversionPhoton *gamma );
    void ProcessQA                  ( AliAODConversionPhoton *gamma );
    void RelabelAODPhotonCandidates ( Bool_t mode );
    void ProcessTrueQAESD           ( AliAODConversionPhoton *TruePhotonCandidate, 
                                      AliESDtrack *elec, 
                                      AliESDtrack *posi );
    void ProcessTrueQAAOD           ( AliAODConversionPhoton *TruePhotonCandidate, 
                                      AliAODTrack *elec,
                                      AliAODTrack *posi );
    UInt_t IsTruePhotonESD          ( AliAODConversionPhoton *TruePhotonCandidate );
    UInt_t IsTruePhotonAOD          ( AliAODConversionPhoton *TruePhotonCandidate );
    UInt_t GetTrueMotherInfoESD     ( AliAODConversionPhoton *TruePhotonCandidate );
    UInt_t GetTrueMotherInfoAOD     ( AliAODConversionPhoton *TruePhotonCandidate );
    void ResetBuffer();
        
    AliV0ReaderV1*              fV0Reader;                  //
    TString                     fV0ReaderName;
    TClonesArray*               fConversionGammas;          //
    AliConversionPhotonCuts*    fConversionCuts;            // Cuts used by the V0Reader
    AliConvEventCuts*           fEventCuts;                 // Cuts used by the V0Reader
    AliVEvent*                  fInputEvent;                //
    AliMCEvent*                 fMCEvent;                   //
    TTree*                      fAnalysisTree;                    //
    Bool_t                      fIsHeavyIon;                //
    Double_t                    ffillTree;                  //
    TList*                      fOutputList;                //
    Int_t                       fBuffer_NConversionCandidates; //
    Float_t*                    fBuffer_ConversionCandidate_E; //
    Float_t*                    fBuffer_ConversionCandidate_Px; //
    Float_t*                    fBuffer_ConversionCandidate_Py; //
    Float_t*                    fBuffer_ConversionCandidate_Pz; //
    Float_t*                    fBuffer_ConversionCandidate_Qt; //
    Float_t*                    fBuffer_ConversionCandidate_Alpha; //
    Float_t*                    fBuffer_ConversionCandidate_PsiPair; //
    Float_t*                    fBuffer_ConversionCandidate_Chi2; //
    Float_t*                    fBuffer_ConversionCandidate_CosPA; //
    Float_t*                    fBuffer_ConversionCandidate_Eta; //
    Float_t*                    fBuffer_ConversionCandidate_Phi; //
    Float_t*                    fBuffer_ConversionCandidate_ConvPointX; //
    Float_t*                    fBuffer_ConversionCandidate_ConvPointY; //
    Float_t*                    fBuffer_ConversionCandidate_ConvPointZ; //
    UInt_t*                     fBuffer_ConversionCandidate_MC_Type;                      //
    Int_t*                      fBuffer_ConversionCandidate_MC_Mother_ID;
    Int_t*                      fBuffer_ConversionCandidate_MC_Mother_PDG;
    UInt_t*                     fBuffer_ConversionCandidate_MC_Mother_Type;
    Float_t*                    fBuffer_ConversionCandidate_MC_Mother_TruePt;
    Int_t*                      fBuffer_ConversionCandidate_MC_GrandMother_PDG;
    Bool_t                      fIsMC;                      //
    Int_t                       fnGammaCandidates;          //
    Int_t*                      fMCStackPos;                //[fnGammaCandidates]
    Int_t*                      fMCStackNeg;                //[fnGammaCandidates]
    
    ClassDef(AliAnalysisTaskConversionTree, 1);
};
const Int_t kMaxConvCandidates = 200;

#endif
