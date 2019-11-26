#ifndef AliAnalysisQA_cxx
#define AliAnalysisQA_cxx

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
    TTree*                      fTreeQA;                    //
    Bool_t                      fIsHeavyIon;                //
    //    Double_t                    ffillTree;                  //
    Bool_t                      ffillHistograms;            //
    TList*                      fOutputList;                //
    TList*                      fESDList;                   //
    TH1F*                       hCentrality;                // centrality distribution for selected events
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
    Float_t                     fGammaPt;                   //
    Float_t                     fGammaTheta;                //
    Float_t                     fGammaChi2NDF;              //
    TVectorF                    fGammaPhotonProp;           //
    TVectorF                    fGammaConvCoord;            //
    TVectorF                    fDaughterProp;              //


    TTree*                      fTree;                //           

    //Event properties
                                Float_t Centrality;//
                                Float_t VertexZ;//
                                Int_t   Bunch;  //
				Int_t   GoodESDTracks;//  

    //Gamma candiates
                                Float_t theta;//   eta = -TMath::Log(TMath::Tan(theta/2.))
				Float_t pt;//      p = pt * TMath::CosH(eta)
                                Float_t phi;//
                                Float_t chi2;//
                                Float_t qt;//
                                Float_t alpha;//
                                Float_t psipair;//
                                Float_t cosPA;//
                                Float_t InvMass;//
                                Float_t X;//
                                Float_t Y;//
                                Float_t Z;//
                                Float_t R;//
                                Int_t   Qual;//
                                Float_t DCAz;//
                                Float_t DCAr;//

    //Track param
                                Float_t ele_theta;//   ele_eta = -TMath::Log(TMath::Tan(ele_theta/2.))
				Float_t ele_pt;//      ele_p = ele_pt * TMath::CosH(ele_eta)
                                Float_t ele_phi;//
				Float_t ele_nSigmaTPC;//
                                Float_t ele_nSigmaTPCpion;//
				Float_t ele_nSigmaTOF;//
				Float_t ele_nSigmaITS;//			
                                Float_t ele_TPCsignal;//
                                Float_t ele_TOFsignal;//
				Float_t ele_ITSsignal;//
                                Float_t ele_Cls;//
				Float_t ele_NfindableCls;//
				Bool_t  ele_SPD1;           //
				Bool_t  ele_SPD2;           //
				Bool_t  ele_SDD1;           //
				Bool_t  ele_SDD2;           //
				Bool_t  ele_SSD1;           //
				Bool_t  ele_SSD2;           //
				
                                Float_t pos_theta;         //pos_eta = -TMath::Log(TMath::Tan(pos_theta/2.))
                                Float_t pos_pt;            //pos_p = pos_pt * TMath::CosH(pos_eta)
                                Float_t pos_phi;           //
				Float_t pos_nSigmaTPC;     //
                                Float_t pos_nSigmaTPCpion; //
				Float_t pos_nSigmaTOF;     //
				Float_t pos_nSigmaITS;     //
                                Float_t pos_TPCsignal;     //
                                Float_t pos_TOFsignal;     //
				Float_t pos_ITSsignal;     //
                                Float_t pos_Cls;           //
				Float_t pos_NfindableCls;  //
				Bool_t  pos_SPD1;           //
				Bool_t  pos_SPD2;           //
				Bool_t  pos_SDD1;           //
				Bool_t  pos_SDD2;           //
				Bool_t  pos_SSD1;           //
				Bool_t  pos_SSD2;           //




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

