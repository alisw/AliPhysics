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
#include "AliConversionCuts.h"
#include "TList.h"
#include "AliStack.h"
#include "TClonesArray.h"


using namespace std;


class AliAnalysisTaskConversionQA : public AliAnalysisTaskSE{

 public:
   
   AliAnalysisTaskConversionQA();
   AliAnalysisTaskConversionQA(const char *name);
   virtual ~AliAnalysisTaskConversionQA();

   virtual void   UserCreateOutputObjects();
   virtual Bool_t Notify();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);

   void SetV0Reader(AliV0ReaderV1 *v0Reader){fV0Reader=v0Reader;}
   void SetConversionCuts(AliConversionCuts* conversionCuts,Bool_t IsHeavyIon ){
      fConversionCuts=conversionCuts;
      fIsHeavyIon = IsHeavyIon;
   }
   void FillType(Bool_t fillTree, Bool_t fillHistorams){
      ffillTree = fillTree;
      ffillHistograms = fillHistorams;
   }
   void SetIsMC(Bool_t isMC){fIsMC = isMC;}
   
 private:
    
   AliAnalysisTaskConversionQA(const AliAnalysisTaskConversionQA&); // Prevent copy-construction
   AliAnalysisTaskConversionQA &operator=(const AliAnalysisTaskConversionQA&); // Prevent assignment

   void ProcessQATree(AliAODConversionPhoton *gamma);
   void ProcessQA(AliAODConversionPhoton *gamma);
   void RelabelAODPhotonCandidates(Bool_t mode);
   void ProcessTrueQAESD(AliAODConversionPhoton *TruePhotonCandidate, AliESDtrack *elec, AliESDtrack *posi);
   void ProcessTrueQAAOD(AliAODConversionPhoton *TruePhotonCandidate, AliAODTrack *elec, AliAODTrack *posi);
   UInt_t IsTruePhotonESD(AliAODConversionPhoton *TruePhotonCandidate);
   UInt_t IsTruePhotonAOD(AliAODConversionPhoton *TruePhotonCandidate);
   void CountTracks();
   
	
   AliV0ReaderV1 *fV0Reader;    
   TClonesArray *fConversionGammas;
   AliConversionCuts *fConversionCuts; // Cuts used by the V0Reader
   AliVEvent *fInputEvent;
   Int_t fNumberOfESDTracks;
   AliMCEvent *fMCEvent;
   AliStack *fMCStack;
   TTree *fTreeQA;
   Bool_t fIsHeavyIon;
   Bool_t ffillTree;
   Bool_t ffillHistograms;
   TList *fOutputList;
   TList *fTreeList;
   TList *fESDList;
   TH1F *hVertexZ;
   TH1I *hNGoodESDTracks;
   TH1I *hNV0Tracks;
   TH1I *hNContributorsVertex;
   TH2F *hITSClusterPhi;
   TH1F *hGammaPt;
   TH1F *hGammaPhi;
   TH1F *hGammaEta;
   TH1F *hGammaChi2perNDF;
   TH1F *hGammaPsiPair;
   TH1F *hGammaQt;
   TH1F *hGammaCosinePointingAngle;
   TH2F *hGammaXY;
   TH2F *hGammaZR;
   TH2F *hElecPt;
   TH2F *hElecEta;
   TH2F *hElecPhi;
   TH1F *hElecNfindableClsTPC;
   TH1F *hPosiNfindableClsTPC;
   TH2F *hElecAsymP;
   TList *fTrueList;
   TH2F *hTrueResolutionR;
   TH2F *hTrueResolutionZ;
   TH2F *hTrueResolutionPhi;
   TH1F *hTrueGammaPt;
   TH1F *hTrueGammaPhi;
   TH1F *hTrueGammaEta;
   TH1F *hTrueGammaMass;
   TH1F *hTrueGammaChi2perNDF;
   TH1F *hTrueGammaPsiPair;
   TH1F *hTrueGammaQt;
   TH1F *hTrueGammaCosinePointingAngle;
   TH2F *hTrueGammaXY;
   TH2F *hTrueGammaZR;
   TH2F *hTrueElecPt;
   TH2F *hTrueElecEta;
   TH2F *hTrueElecPhi;
   TH1F *hTrueElecNfindableClsTPC;
   TH1F *hTruePosiNfindableClsTPC;
   TH2F *hTrueElecAsymP;
   Float_t fGammaPt;
   Float_t fGammaTheta;
   Float_t fGammaChi2NDF;
   Float_t fGammaPhotonProp[5];
   Float_t fGammaConvCoord[5];
   Float_t fDaughterProp[20];
   UInt_t fKind;
   Bool_t fIsMC;
   Int_t fnGammaCandidates;
   Int_t *fMCStackPos;     //[fnGammaCandidates]
   Int_t *fMCStackNeg;     //[fnGammaCandidates]
   ClassDef(AliAnalysisTaskConversionQA, 3);
};

#endif

