#ifndef AliAnalysisConversionQA_cxx
#define AliAnalysisConversionQA_cxx

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


class AliAnalysisTaskConversionQA : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskConversionQA();
    AliAnalysisTaskConversionQA(const char *name);
    virtual ~AliAnalysisTaskConversionQA();

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

    void FillType                           ( Bool_t fillTree, 
                                              Bool_t fillHistorams)                       {
                                                                                            ffillTree = fillTree                ;
                                                                                            ffillHistograms = fillHistorams     ;
                                                                                          }
    void SetIsMC                            ( Bool_t isMC )                               { fIsMC = isMC                        ; }

  private:
        
    AliAnalysisTaskConversionQA     ( const AliAnalysisTaskConversionQA& ); // Prevent copy-construction
    AliAnalysisTaskConversionQA &operator=( const AliAnalysisTaskConversionQA& ); // Prevent assignment

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
    void CountTracks                ();
    void SetLogBinningXTH2          ( TH2* histoRebin );
        
    AliV0ReaderV1*              fV0Reader;                  //
    TString                     fV0ReaderName;
    TClonesArray*               fConversionGammas;          //
    AliConversionPhotonCuts*    fConversionCuts;            // Cuts used by the V0Reader
    AliConvEventCuts*           fEventCuts;                 // Cuts used by the V0Reader
    AliVEvent*                  fInputEvent;                //
    Int_t                       fNumberOfESDTracks;         //
    AliMCEvent*                 fMCEvent;                   //
    TTree*                      fTreeQA;                    //
    Bool_t                      fIsHeavyIon;                //
    Bool_t                      ffillTree;                  //
    Bool_t                      ffillHistograms;            //
    TList*                      fOutputList;                //
    TList*                      fTreeList;                  //
    TList*                      fESDList;                   //
    TH1F*                       hVertexZ;                   //
    TH1I*                       hNGoodESDTracks;            //
    TH1I*                       hNV0Tracks;                 //
    TH1I*                       hNContributorsVertex;       //
    TH2F*                       hITSClusterPhi;             //
    TH1F*                       hGammaPt;                   //
    TH1F*                       hGammaPhi;                  //
    TH1F*                       hGammaPhi_Pos;              //
    TH1F*                       hGammaPhi_Neg;              //
    TH1F*                       hGammaEta;                  //
    TH1F*                       hGammaChi2perNDF;           //
    TH1F*                       hGammaPsiPair;              //
    TH2F*                       hGammaArmenteros;           //
    TH1F*                       hGammaCosinePointingAngle;  //
    TH1F*                       hGammaInvMass;              //
    TH2F*                       hElecPt;                    //
    TH2F*                       hElecEta;                   //
    TH2F*                       hElecPhi;                   //
    TH1F*                       hElecNfindableClsTPC;       //
    TH1F*                       hPosiNfindableClsTPC;       //
    TH1F*                       hElecClsTPC;                //
    TH1F*                       hPosiClsTPC;                //
    TH2F*                       hElectrondEdxP;             //
    TH2F*                       hElectronITSdEdxP;          //
    TH2F*                       hElectronTOFP;              //
    TH2F*                       hElectronNSigmadEdxP;       //
    TH2F*                       hElectronNSigmadEdxEta;     //
    TH2F*                       hElectronNSigmaPiondEdxP;   //
    TH2F*                       hElectronNSigmaITSP;        //
    TH2F*                       hElectronNSigmaTOFP;        //
    TH2F*                       hPositrondEdxP;             //
    TH2F*                       hPositronITSdEdxP;          //
    TH2F*                       hPositronTOFP;              //
    TH2F*                       hPositronNSigmadEdxP;       //
    TH2F*                       hPositronNSigmadEdxEta;     //
    TH2F*                       hPositronNSigmaPiondEdxP;   //
    TH2F*                       hPositronNSigmaITSP;        //
    TH2F*                       hPositronNSigmaTOFP;        //
    TH2F*                       hInvMassPair;               //
    Float_t                     fGammaPt;                   //
    Float_t                     fGammaTheta;                //
    Float_t                     fGammaChi2NDF;              //
    TVectorF                    fGammaPhotonProp;           //
    TVectorF                    fGammaConvCoord;            //
    TVectorF                    fDaughterProp;              //
    UInt_t                      fKind;                      //
    Bool_t                      fIsMC;                      //
    Int_t                       fnGammaCandidates;          //
    Int_t*                      fMCStackPos;                //[fnGammaCandidates]
    Int_t*                      fMCStackNeg;                //[fnGammaCandidates]
    
    ClassDef(AliAnalysisTaskConversionQA, 7);
};

#endif

