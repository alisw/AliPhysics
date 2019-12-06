#ifndef AliAnalysisQA_cxx
#define AliAnalysisQA_cxx

#include "AliAnalysisTaskSE.h"
#include "AliConversionPhotonBase.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "AliV0ReaderV1.h"
#include "AliConversionPhotonCuts.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "TClonesArray.h"



using namespace std;


class AliAnalysisTaskQA : public AliAnalysisTaskSE{

  public:

    AliAnalysisTaskQA();
    AliAnalysisTaskQA(const char *name);
    virtual ~AliAnalysisTaskQA();

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

    /* void FillType                           ( Double_t fillTree, Bool_t fillHistorams)    { ffillTree = fillTree                ; */
    /*                                                                                         ffillHistograms = fillHistorams     ;} */
    void FillType                           (Bool_t fillHistorams)                        { ffillHistograms = fillHistorams     ;}
                                            
    void SetIsMC                            ( Bool_t isMC )                               { fIsMC = isMC                        ; }

    // TStrings with active or inactive branches
    void SetTreeActiveBranch(TString b)   {fActiveBranches+=b+";";}
    void SetTreeInactiveBranch(TString b) {fInactiveBranches+=b+";";}
    void SetWriteVariableTree(Bool_t a)   {fWriteVariableTree = a;}


    
  private:
        
    AliAnalysisTaskQA     ( const AliAnalysisTaskQA& ); // Prevent copy-construction
    AliAnalysisTaskQA &operator=( const AliAnalysisTaskQA& ); // Prevent assignment

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
    //    void CountTracks                ();
    void SetLogBinningXTH2          ( TH2* histoRebin );
        
    AliV0ReaderV1*              fV0Reader;                  //
    TString                     fV0ReaderName;
    TClonesArray*               fConversionGammas;          //
    AliConversionPhotonCuts*    fConversionCuts;            // Cuts used by the V0Reader
    AliConvEventCuts*           fEventCuts;                 // Cuts used by the V0Reader
    AliVEvent*                  fInputEvent;                //
    Int_t                       fNumberOfESDTracks;         //
    AliMCEvent*                 fMCEvent;                   //
  //    TTree*                      fTreeQA;                    //
    Bool_t                      fIsHeavyIon;                //
    //    Double_t                    ffillTree;                  //
    Bool_t                      ffillHistograms;            //
    TList*                      fOutputList;                //
    TList*                      fESDList;                   //
    TH1F*                       hCentralityV0A;             // centrality distribution for selected events
    TH1I*                       hBunch;                     // centrality distribution for selected events
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

    TTree*                      fTree;                      //           

    //Event properties
                                Float_t fCentralityV0M;//
                                Float_t fCentralityV0A;//
                                Float_t fCentralityV0C;//
                                Int_t   fRunNumber;  
                                Float_t fVertexZ;//
                                Int_t   fBunch;  //
				Int_t   fGoodESDTracks;//  

    //Gamma candiates
                                Float_t ftheta;//   eta = -TMath::Log(TMath::Tan(theta/2.))
				Float_t fpt;//      p = pt * TMath::CosH(eta)
                                Float_t fphi;//
                                Float_t fchi2;//
                                Float_t fqt;//
                                Float_t falpha;//
                                Float_t fpsipair;//
                                Float_t fcosPA;//
                                Float_t fInvMass;//
                                Float_t fX;//
                                Float_t fY;//
                                Float_t fZ;//
                                Float_t fR;//
                                Int_t   fQual;//
                                Float_t fDCAz;//
                                Float_t fDCAr;//

    //Track param
                                Float_t fele_theta;         //ele_eta = -TMath::Log(TMath::Tan(ele_theta/2.))
				Float_t fele_pt;            //ele_p = ele_pt * TMath::CosH(ele_eta)
                                Float_t fele_phi;           //
				Float_t fele_nSigmaTPC;     //
                                Float_t fele_nSigmaTPCpion; //
				Float_t fele_nSigmaTOF;     //
				Float_t fele_nSigmaITS;     //			
                                Float_t fele_TPCsignal;     //
                                Float_t fele_TOFsignal;     //
				Float_t fele_ITSsignal;     //
                                Float_t fele_Cls;           //
				Float_t fele_NfindableCls;  //
				Bool_t  fele_SPD1;          //
				Bool_t  fele_SPD2;          //
				Bool_t  fele_SDD1;          //
				Bool_t  fele_SDD2;          //
				Bool_t  fele_SSD1;          //
				Bool_t  fele_SSD2;          //
				
                                Float_t fpos_theta;         //pos_eta = -TMath::Log(TMath::Tan(pos_theta/2.))
                                Float_t fpos_pt;            //pos_p = pos_pt * TMath::CosH(pos_eta)
                                Float_t fpos_phi;           //
				Float_t fpos_nSigmaTPC;     //
                                Float_t fpos_nSigmaTPCpion; //
				Float_t fpos_nSigmaTOF;     //
				Float_t fpos_nSigmaITS;     //
                                Float_t fpos_TPCsignal;     //
                                Float_t fpos_TOFsignal;     //
				Float_t fpos_ITSsignal;     //
                                Float_t fpos_Cls;           //
				Float_t fpos_NfindableCls;  //
				Bool_t  fpos_SPD1;          //
				Bool_t  fpos_SPD2;          //
				Bool_t  fpos_SDD1;          //
				Bool_t  fpos_SDD2;          //
				Bool_t  fpos_SSD1;          //
				Bool_t  fpos_SSD2;          //


    UInt_t                      fKind;                      //
    Bool_t                      fIsMC;                      //
    Int_t                       fnGammaCandidates;          //
    Int_t*                      fMCStackPos;                //[fnGammaCandidates]
    Int_t*                      fMCStackNeg;                //[fnGammaCandidates]

    TString fActiveBranches;    // list of active output tree branches
    TString fInactiveBranches;  // list of inactive output tree branches
    Bool_t  fWriteVariableTree; // flag to decide whether to write the candidate variables on a tree variables

    
    ClassDef(AliAnalysisTaskQA, 1);
};

#endif

