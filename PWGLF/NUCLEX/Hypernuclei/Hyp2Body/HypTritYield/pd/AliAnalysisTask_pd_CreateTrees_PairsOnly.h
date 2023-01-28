#ifndef AliAnalysisTask_pd_CreateTrees_PairsOnly_H
#define AliAnalysisTask_pd_CreateTrees_PairsOnly_H

#include "AliAnalysisTaskSE.h"
#include "TObject.h"

class AliAODEvent;
class AliAODInputHandler;







class AliAnalysisTask_pd_CreateTrees_PairsOnly : public AliAnalysisTaskSE
{

  public:

    AliAnalysisTask_pd_CreateTrees_PairsOnly();
    AliAnalysisTask_pd_CreateTrees_PairsOnly(const char *name,int CollisionSystem);
    AliAnalysisTask_pd_CreateTrees_PairsOnly& operator = (const AliAnalysisTask_pd_CreateTrees_PairsOnly &task);
    AliAnalysisTask_pd_CreateTrees_PairsOnly(const AliAnalysisTask_pd_CreateTrees_PairsOnly &task);
    virtual ~AliAnalysisTask_pd_CreateTrees_PairsOnly();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    double CalculateBetaTOF(AliAODTrack &track); 
    double CalculateMassSquareTOF(AliAODTrack &track);
    double CalculateDeuteronSigmaMassSquareTOF(double pT, double massSq, bool isMatter);
    bool CheckProtonCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter);
    bool CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter);


  private:

    AliAODEvent		  *fAODEvent;
    AliAODInputHandler	  *fAODHandler; 
    AliAODHeader	  *fHeader;
    AliPIDResponse	  *fPIDResponse;

    int			  fCollisionSystem;

    TTree     *fSaveTree_Proton;
    float     fProton_px;
    float     fProton_py;
    float     fProton_pz;
    float     fProton_pTPC;
    float     fProton_Eta;
    float     fProton_Phi;
    float     fProton_TPC_Chi2overNDF;
    float     fProton_TPC_dEdx;
    float     fProton_TPC_dEdx_Sigma;
    float     fProton_TOF_Beta;
    float     fProton_TOF_Beta_Sigma;
    float     fProton_TOF_Mass2;
    float     fProton_TOF_Mass2_Sigma;
    float     fProton_ITS_dEdx;
    float     fProton_ITS_dEdx_Sigma;
    float     fProton_DCAxy;
    float     fProton_DCAz;
    float     fProton_Event_Centrality;
    float     fProton_Event_PrimaryVertexZ;
    float     fProton_Event_BField;
    unsigned short    fProton_TPC_nCrossedRows;
    unsigned short    fProton_TPC_nSharedCluster;
    unsigned short    fProton_TPC_nClusterFindable;
    unsigned short    fProton_TPC_nCluster;
    unsigned short    fProton_ITS_nCluster;
    unsigned short    fProton_nParticlesPerEvent;
    unsigned int      fProton_ID;
    unsigned int      fProton_Event_Identifier;


    TTree     *fSaveTree_Deuteron;
    float     fDeuteron_px;
    float     fDeuteron_py;
    float     fDeuteron_pz;
    float     fDeuteron_pTPC;
    float     fDeuteron_Eta;
    float     fDeuteron_Phi;
    float     fDeuteron_TPC_Chi2overNDF;
    float     fDeuteron_TPC_dEdx;
    float     fDeuteron_TPC_dEdx_Sigma;
    float     fDeuteron_TOF_Beta;
    float     fDeuteron_TOF_Beta_Sigma;
    float     fDeuteron_TOF_Mass2;
    float     fDeuteron_TOF_Mass2_Sigma;
    float     fDeuteron_ITS_dEdx;
    float     fDeuteron_ITS_dEdx_Sigma;
    float     fDeuteron_DCAxy;
    float     fDeuteron_DCAz;
    float     fDeuteron_Event_Centrality;
    float     fDeuteron_Event_PrimaryVertexZ;
    float     fDeuteron_Event_BField;
    unsigned short    fDeuteron_TPC_nCrossedRows;
    unsigned short    fDeuteron_TPC_nSharedCluster;
    unsigned short    fDeuteron_TPC_nClusterFindable;
    unsigned short    fDeuteron_TPC_nCluster;
    unsigned short    fDeuteron_ITS_nCluster;
    unsigned short    fDeuteron_nParticlesPerEvent;
    unsigned int      fDeuteron_ID;
    unsigned int      fDeuteron_Event_Identifier;


    TTree     *fSaveTree_AntiProton;
    float     fAntiProton_px;
    float     fAntiProton_py;
    float     fAntiProton_pz;
    float     fAntiProton_pTPC;
    float     fAntiProton_Eta;
    float     fAntiProton_Phi;
    float     fAntiProton_TPC_Chi2overNDF;
    float     fAntiProton_TPC_dEdx;
    float     fAntiProton_TPC_dEdx_Sigma;
    float     fAntiProton_TOF_Beta;
    float     fAntiProton_TOF_Beta_Sigma;
    float     fAntiProton_TOF_Mass2;
    float     fAntiProton_TOF_Mass2_Sigma;
    float     fAntiProton_ITS_dEdx;
    float     fAntiProton_ITS_dEdx_Sigma;
    float     fAntiProton_DCAxy;
    float     fAntiProton_DCAz;
    float     fAntiProton_Event_Centrality;
    float     fAntiProton_Event_PrimaryVertexZ;
    float     fAntiProton_Event_BField;
    unsigned short    fAntiProton_TPC_nCrossedRows;
    unsigned short    fAntiProton_TPC_nSharedCluster;
    unsigned short    fAntiProton_TPC_nClusterFindable;
    unsigned short    fAntiProton_TPC_nCluster;
    unsigned short    fAntiProton_ITS_nCluster;
    unsigned short    fAntiProton_nParticlesPerEvent;
    unsigned int      fAntiProton_ID;
    unsigned int      fAntiProton_Event_Identifier;


    TTree     *fSaveTree_AntiDeuteron;
    float     fAntiDeuteron_px;
    float     fAntiDeuteron_py;
    float     fAntiDeuteron_pz;
    float     fAntiDeuteron_pTPC;
    float     fAntiDeuteron_Eta;
    float     fAntiDeuteron_Phi;
    float     fAntiDeuteron_TPC_Chi2overNDF;
    float     fAntiDeuteron_TPC_dEdx;
    float     fAntiDeuteron_TPC_dEdx_Sigma;
    float     fAntiDeuteron_TOF_Beta;
    float     fAntiDeuteron_TOF_Beta_Sigma;
    float     fAntiDeuteron_TOF_Mass2;
    float     fAntiDeuteron_TOF_Mass2_Sigma;
    float     fAntiDeuteron_ITS_dEdx;
    float     fAntiDeuteron_ITS_dEdx_Sigma;
    float     fAntiDeuteron_DCAxy;
    float     fAntiDeuteron_DCAz;
    float     fAntiDeuteron_Event_Centrality;
    float     fAntiDeuteron_Event_PrimaryVertexZ;
    float     fAntiDeuteron_Event_BField;
    unsigned short    fAntiDeuteron_TPC_nCrossedRows;
    unsigned short    fAntiDeuteron_TPC_nSharedCluster;
    unsigned short    fAntiDeuteron_TPC_nClusterFindable;
    unsigned short    fAntiDeuteron_TPC_nCluster;
    unsigned short    fAntiDeuteron_ITS_nCluster;
    unsigned short    fAntiDeuteron_nParticlesPerEvent;
    unsigned int      fAntiDeuteron_ID;
    unsigned int      fAntiDeuteron_Event_Identifier;









    ClassDef(AliAnalysisTask_pd_CreateTrees_PairsOnly,1);

};






#endif
