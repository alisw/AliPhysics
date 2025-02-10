#ifndef AliAnalysisTask_pd_CreateTrees_PairsOnly_H
#define AliAnalysisTask_pd_CreateTrees_PairsOnly_H

#include "AliAnalysisTaskSE.h"
#include "TObject.h"
#include "AliEventCuts.h"

class AliAODEvent;
class AliAODInputHandler;







class AliAnalysisTask_pd_CreateTrees_PairsOnly : public AliAnalysisTaskSE
{

  public:

    AliAnalysisTask_pd_CreateTrees_PairsOnly();
    AliAnalysisTask_pd_CreateTrees_PairsOnly(const char *name,Int_t CollisionSystem, Bool_t UseOpenCuts, Bool_t isMC, Bool_t SavePairsOnly);
    AliAnalysisTask_pd_CreateTrees_PairsOnly& operator = (const AliAnalysisTask_pd_CreateTrees_PairsOnly &task);
    AliAnalysisTask_pd_CreateTrees_PairsOnly(const AliAnalysisTask_pd_CreateTrees_PairsOnly &task);
    virtual ~AliAnalysisTask_pd_CreateTrees_PairsOnly();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    Double_t CalculateBetaTOF(AliAODTrack &track); 
    Double_t CalculateMassSquareTOF(AliAODTrack &track);
    Double_t CalculateSigmaMassSquareTOF(Double_t pT, Double_t massSq, Int_t ParticleSpecies, Int_t RunNumber);
    Bool_t CheckProtonCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, Int_t ParticleSpecies, Int_t RunNumber);
    Bool_t CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, Int_t ParticleSpecies, Int_t RunNumber);
    Double_t CalculateSigmadEdxTPC(AliAODTrack &Track, Int_t ParticleSpecies);
    Double_t CalculateSigmadEdxITS(AliAODTrack &Track, Int_t ParticleSpecies, Int_t RunNumber);


  private:

    AliAODEvent		  *fAODEvent;
    AliAODInputHandler	  *fAODHandler; 
    AliAODHeader	  *fHeader;
    AliPIDResponse	  *fPIDResponse;

    Int_t   fCollisionSystem;
    Bool_t  fUseOpenCuts;
    Bool_t  fIsMC;
    Bool_t  fSavePairsOnly;
    AliTimeRangeCut fTimeRangeCut;

    TTree *fSaveTree_Proton;
    Float_t   fProton_px;
    Float_t   fProton_py;
    Float_t   fProton_pz;
    Float_t   fProton_px_Generated;
    Float_t   fProton_py_Generated;
    Float_t   fProton_pz_Generated;
    Float_t   fProton_pTPC;
    Float_t   fProton_Eta;
    Float_t   fProton_Phi;
    Float_t   fProton_TPC_Chi2;
    Float_t   fProton_TPC_dEdx;
    Float_t   fProton_TPC_dEdx_Sigma;
    Float_t   fProton_TOF_Mass2;
    Float_t   fProton_TOF_Mass2_Sigma;
    Float_t   fProton_ITS_dEdx;
    Float_t   fProton_ITS_dEdx_Sigma;
    Float_t   fProton_DCAxy;
    Float_t   fProton_DCAz;
    Float_t   fProton_Event_Centrality;
    Float_t   fProton_Event_PrimaryVertexX;
    Float_t   fProton_Event_PrimaryVertexY;
    Float_t   fProton_Event_PrimaryVertexZ;
    Bool_t    fProton_Event_BField;
    UShort_t  fProton_TPC_nCrossedRows;
    UShort_t  fProton_TPC_nFindableCluster;
    UShort_t  fProton_TPC_nCluster;
    UShort_t  fProton_ITS_nCluster;
    Int_t     fProton_PDG;
    Int_t     fProton_MotherPDG;
    Bool_t    fProton_FilterBit;
    Int_t    fProton_Event_Multiplicity;
    ULong64_t fProton_Event_Identifier;
    Int_t     fProton_Event_RunNumber;
    Bool_t    fProton_ITS_Layer0;
    Bool_t    fProton_ITS_Layer1;
    Bool_t    fProton_ITS_Layer2;
    Bool_t    fProton_ITS_Layer3;
    Bool_t    fProton_ITS_Layer4;
    Bool_t    fProton_ITS_Layer5;
    Bool_t    fProton_Event_IsFirstParticle;
    UInt_t    fProton_Event_TimeStamp;
    UInt_t  fProton_Event_RandomCrossCheckNumber;
    Int_t     fProton_ID;
    Float_t     fProton_x;
    Float_t     fProton_y;
    Float_t     fProton_z;


    TTree *fSaveTree_Deuteron;
    Float_t   fDeuteron_px;
    Float_t   fDeuteron_py;
    Float_t   fDeuteron_pz;
    Float_t   fDeuteron_px_Generated;
    Float_t   fDeuteron_py_Generated;
    Float_t   fDeuteron_pz_Generated;
    Float_t   fDeuteron_pTPC;
    Float_t   fDeuteron_Eta;
    Float_t   fDeuteron_Phi;
    Float_t   fDeuteron_TPC_Chi2;
    Float_t   fDeuteron_TPC_dEdx;
    Float_t   fDeuteron_TPC_dEdx_Sigma;
    Float_t   fDeuteron_TOF_Mass2;
    Float_t   fDeuteron_TOF_Mass2_Sigma;
    Float_t   fDeuteron_ITS_dEdx;
    Float_t   fDeuteron_ITS_dEdx_Sigma;
    Float_t   fDeuteron_DCAxy;
    Float_t   fDeuteron_DCAz;
    Float_t   fDeuteron_Event_Centrality;
    Float_t   fDeuteron_Event_PrimaryVertexX;
    Float_t   fDeuteron_Event_PrimaryVertexY;
    Float_t   fDeuteron_Event_PrimaryVertexZ;
    Bool_t    fDeuteron_Event_BField;
    UShort_t  fDeuteron_TPC_nCrossedRows;
    UShort_t  fDeuteron_TPC_nFindableCluster;
    UShort_t  fDeuteron_TPC_nCluster;
    UShort_t  fDeuteron_ITS_nCluster;
    Int_t     fDeuteron_PDG;
    Int_t     fDeuteron_MotherPDG;
    Bool_t    fDeuteron_FilterBit;
    Int_t    fDeuteron_Event_Multiplicity;
    ULong64_t fDeuteron_Event_Identifier;
    Int_t     fDeuteron_Event_RunNumber;
    Bool_t    fDeuteron_ITS_Layer0;
    Bool_t    fDeuteron_ITS_Layer1;
    Bool_t    fDeuteron_ITS_Layer2;
    Bool_t    fDeuteron_ITS_Layer3;
    Bool_t    fDeuteron_ITS_Layer4;
    Bool_t    fDeuteron_ITS_Layer5;
    Bool_t    fDeuteron_Event_IsFirstParticle;
    UInt_t    fDeuteron_Event_TimeStamp;
    UInt_t  fDeuteron_Event_RandomCrossCheckNumber;
    Int_t     fDeuteron_ID;
    Float_t     fDeuteron_x;
    Float_t     fDeuteron_y;
    Float_t     fDeuteron_z;


    TTree *fSaveTree_AntiProton;
    Float_t   fAntiProton_px;
    Float_t   fAntiProton_py;
    Float_t   fAntiProton_pz;
    Float_t   fAntiProton_px_Generated;
    Float_t   fAntiProton_py_Generated;
    Float_t   fAntiProton_pz_Generated;
    Float_t   fAntiProton_pTPC;
    Float_t   fAntiProton_Eta;
    Float_t   fAntiProton_Phi;
    Float_t   fAntiProton_TPC_Chi2;
    Float_t   fAntiProton_TPC_dEdx;
    Float_t   fAntiProton_TPC_dEdx_Sigma;
    Float_t   fAntiProton_TOF_Mass2;
    Float_t   fAntiProton_TOF_Mass2_Sigma;
    Float_t   fAntiProton_ITS_dEdx;
    Float_t   fAntiProton_ITS_dEdx_Sigma;
    Float_t   fAntiProton_DCAxy;
    Float_t   fAntiProton_DCAz;
    Float_t   fAntiProton_Event_Centrality;
    Float_t   fAntiProton_Event_PrimaryVertexX;
    Float_t   fAntiProton_Event_PrimaryVertexY;
    Float_t   fAntiProton_Event_PrimaryVertexZ;
    Bool_t    fAntiProton_Event_BField;
    UShort_t  fAntiProton_TPC_nCrossedRows;
    UShort_t  fAntiProton_TPC_nFindableCluster;
    UShort_t  fAntiProton_TPC_nCluster;
    UShort_t  fAntiProton_ITS_nCluster;
    Int_t     fAntiProton_PDG;
    Int_t     fAntiProton_MotherPDG;
    Bool_t    fAntiProton_FilterBit;
    Int_t    fAntiProton_Event_Multiplicity;
    ULong64_t fAntiProton_Event_Identifier;
    Int_t     fAntiProton_Event_RunNumber;
    Bool_t    fAntiProton_ITS_Layer0;
    Bool_t    fAntiProton_ITS_Layer1;
    Bool_t    fAntiProton_ITS_Layer2;
    Bool_t    fAntiProton_ITS_Layer3;
    Bool_t    fAntiProton_ITS_Layer4;
    Bool_t    fAntiProton_ITS_Layer5;
    Bool_t    fAntiProton_Event_IsFirstParticle;
    UInt_t    fAntiProton_Event_TimeStamp;
    UInt_t  fAntiProton_Event_RandomCrossCheckNumber;
    Int_t     fAntiProton_ID;
    Float_t     fAntiProton_x;
    Float_t     fAntiProton_y;
    Float_t     fAntiProton_z;


    TTree *fSaveTree_AntiDeuteron;
    Float_t   fAntiDeuteron_px;
    Float_t   fAntiDeuteron_py;
    Float_t   fAntiDeuteron_pz;
    Float_t   fAntiDeuteron_px_Generated;
    Float_t   fAntiDeuteron_py_Generated;
    Float_t   fAntiDeuteron_pz_Generated;
    Float_t   fAntiDeuteron_pTPC;
    Float_t   fAntiDeuteron_Eta;
    Float_t   fAntiDeuteron_Phi;
    Float_t   fAntiDeuteron_TPC_Chi2;
    Float_t   fAntiDeuteron_TPC_dEdx;
    Float_t   fAntiDeuteron_TPC_dEdx_Sigma;
    Float_t   fAntiDeuteron_TOF_Mass2;
    Float_t   fAntiDeuteron_TOF_Mass2_Sigma;
    Float_t   fAntiDeuteron_ITS_dEdx;
    Float_t   fAntiDeuteron_ITS_dEdx_Sigma;
    Float_t   fAntiDeuteron_DCAxy;
    Float_t   fAntiDeuteron_DCAz;
    Float_t   fAntiDeuteron_Event_Centrality;
    Float_t   fAntiDeuteron_Event_PrimaryVertexX;
    Float_t   fAntiDeuteron_Event_PrimaryVertexY;
    Float_t   fAntiDeuteron_Event_PrimaryVertexZ;
    Bool_t    fAntiDeuteron_Event_BField;
    UShort_t  fAntiDeuteron_TPC_nCrossedRows;
    UShort_t  fAntiDeuteron_TPC_nFindableCluster;
    UShort_t  fAntiDeuteron_TPC_nCluster;
    UShort_t  fAntiDeuteron_ITS_nCluster;
    Int_t     fAntiDeuteron_PDG;
    Int_t     fAntiDeuteron_MotherPDG;
    Bool_t    fAntiDeuteron_FilterBit;
    Int_t    fAntiDeuteron_Event_Multiplicity;
    ULong64_t fAntiDeuteron_Event_Identifier;
    Int_t     fAntiDeuteron_Event_RunNumber;
    Bool_t    fAntiDeuteron_ITS_Layer0;
    Bool_t    fAntiDeuteron_ITS_Layer1;
    Bool_t    fAntiDeuteron_ITS_Layer2;
    Bool_t    fAntiDeuteron_ITS_Layer3;
    Bool_t    fAntiDeuteron_ITS_Layer4;
    Bool_t    fAntiDeuteron_ITS_Layer5;
    Bool_t    fAntiDeuteron_Event_IsFirstParticle;
    UInt_t    fAntiDeuteron_Event_TimeStamp;
    UInt_t  fAntiDeuteron_Event_RandomCrossCheckNumber;
    Int_t     fAntiDeuteron_ID;
    Float_t     fAntiDeuteron_x;
    Float_t     fAntiDeuteron_y;
    Float_t     fAntiDeuteron_z;


    ClassDef(AliAnalysisTask_pd_CreateTrees_PairsOnly,1);

};






#endif
