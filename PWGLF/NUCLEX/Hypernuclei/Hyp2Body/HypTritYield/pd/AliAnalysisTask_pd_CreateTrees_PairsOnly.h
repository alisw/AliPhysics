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
    AliAnalysisTask_pd_CreateTrees_PairsOnly(const char *name,int CollisionSystem, bool UseOpenCuts);
    AliAnalysisTask_pd_CreateTrees_PairsOnly& operator = (const AliAnalysisTask_pd_CreateTrees_PairsOnly &task);
    AliAnalysisTask_pd_CreateTrees_PairsOnly(const AliAnalysisTask_pd_CreateTrees_PairsOnly &task);
    virtual ~AliAnalysisTask_pd_CreateTrees_PairsOnly();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    double CalculateBetaTOF(AliAODTrack &track); 
    double CalculateMassSquareTOF(AliAODTrack &track);
    double CalculateSigmaMassSquareTOF(double pT, double massSq, int ParticleSpecies, int RunNumber);
    bool CheckProtonCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, bool UseOpenCuts);
    bool CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, bool UseOpenCuts);
    double CalculateSigmadEdxITS(AliAODTrack &Track, int ParticleSpecies, int RunNumber);


  private:

    AliAODEvent		  *fAODEvent;
    AliAODInputHandler	  *fAODHandler; 
    AliAODHeader	  *fHeader;
    AliPIDResponse	  *fPIDResponse;

    int			  fCollisionSystem;
    bool		  fUseOpenCuts;

    TTree     *fSaveTree_Proton;
    float     fProton_px;
    float     fProton_py;
    float     fProton_pz;
    float     fProton_pTPC;
    float     fProton_Eta;
    float     fProton_Phi;
    float     fProton_TPC_Chi2;
    float     fProton_TPC_dEdx;
    float     fProton_TPC_dEdx_nSigma;
    float     fProton_TOF_Mass2;
    float     fProton_TOF_Mass2_nSigma;
    float     fProton_ITS_dEdx;
    float     fProton_ITS_dEdx_nSigma;
    float     fProton_DCAxy;
    float     fProton_DCAz;
    float     fProton_Event_Centrality;
    float     fProton_Event_PrimaryVertexZ;
    float     fProton_Event_BField;
    unsigned short    fProton_TPC_nCrossedRows;
    unsigned short    fProton_TPC_nSharedCluster;
    unsigned short    fProton_TPC_nFindableCluster;
    unsigned short    fProton_TPC_nCluster;
    unsigned short    fProton_ITS_nCluster;
    unsigned int      fProton_ID;
    unsigned int      fProton_Event_Multiplicity;
    unsigned long     fProton_Event_Identifier;


    TTree     *fSaveTree_Deuteron;
    float     fDeuteron_px;
    float     fDeuteron_py;
    float     fDeuteron_pz;
    float     fDeuteron_pTPC;
    float     fDeuteron_Eta;
    float     fDeuteron_Phi;
    float     fDeuteron_TPC_Chi2;
    float     fDeuteron_TPC_dEdx;
    float     fDeuteron_TPC_dEdx_nSigma;
    float     fDeuteron_TOF_Mass2;
    float     fDeuteron_TOF_Mass2_nSigma;
    float     fDeuteron_ITS_dEdx;
    float     fDeuteron_ITS_dEdx_nSigma;
    float     fDeuteron_DCAxy;
    float     fDeuteron_DCAz;
    float     fDeuteron_Event_Centrality;
    float     fDeuteron_Event_PrimaryVertexZ;
    float     fDeuteron_Event_BField;
    unsigned short    fDeuteron_TPC_nCrossedRows;
    unsigned short    fDeuteron_TPC_nSharedCluster;
    unsigned short    fDeuteron_TPC_nFindableCluster;
    unsigned short    fDeuteron_TPC_nCluster;
    unsigned short    fDeuteron_ITS_nCluster;
    unsigned int      fDeuteron_ID;
    unsigned int      fDeuteron_Event_Multiplicity;
    unsigned long     fDeuteron_Event_Identifier;


    TTree     *fSaveTree_AntiProton;
    float     fAntiProton_px;
    float     fAntiProton_py;
    float     fAntiProton_pz;
    float     fAntiProton_pTPC;
    float     fAntiProton_Eta;
    float     fAntiProton_Phi;
    float     fAntiProton_TPC_Chi2;
    float     fAntiProton_TPC_dEdx;
    float     fAntiProton_TPC_dEdx_nSigma;
    float     fAntiProton_TOF_Mass2;
    float     fAntiProton_TOF_Mass2_nSigma;
    float     fAntiProton_ITS_dEdx;
    float     fAntiProton_ITS_dEdx_nSigma;
    float     fAntiProton_DCAxy;
    float     fAntiProton_DCAz;
    float     fAntiProton_Event_Centrality;
    float     fAntiProton_Event_PrimaryVertexZ;
    float     fAntiProton_Event_BField;
    unsigned short    fAntiProton_TPC_nCrossedRows;
    unsigned short    fAntiProton_TPC_nSharedCluster;
    unsigned short    fAntiProton_TPC_nFindableCluster;
    unsigned short    fAntiProton_TPC_nCluster;
    unsigned short    fAntiProton_ITS_nCluster;
    unsigned int      fAntiProton_ID;
    unsigned int      fAntiProton_Event_Multiplicity;
    unsigned long     fAntiProton_Event_Identifier;


    TTree     *fSaveTree_AntiDeuteron;
    float     fAntiDeuteron_px;
    float     fAntiDeuteron_py;
    float     fAntiDeuteron_pz;
    float     fAntiDeuteron_pTPC;
    float     fAntiDeuteron_Eta;
    float     fAntiDeuteron_Phi;
    float     fAntiDeuteron_TPC_Chi2;
    float     fAntiDeuteron_TPC_dEdx;
    float     fAntiDeuteron_TPC_dEdx_nSigma;
    float     fAntiDeuteron_TOF_Mass2;
    float     fAntiDeuteron_TOF_Mass2_nSigma;
    float     fAntiDeuteron_ITS_dEdx;
    float     fAntiDeuteron_ITS_dEdx_nSigma;
    float     fAntiDeuteron_DCAxy;
    float     fAntiDeuteron_DCAz;
    float     fAntiDeuteron_Event_Centrality;
    float     fAntiDeuteron_Event_PrimaryVertexZ;
    float     fAntiDeuteron_Event_BField;
    unsigned short    fAntiDeuteron_TPC_nCrossedRows;
    unsigned short    fAntiDeuteron_TPC_nSharedCluster;
    unsigned short    fAntiDeuteron_TPC_nFindableCluster;
    unsigned short    fAntiDeuteron_TPC_nCluster;
    unsigned short    fAntiDeuteron_ITS_nCluster;
    unsigned int      fAntiDeuteron_ID;
    unsigned int      fAntiDeuteron_Event_Multiplicity;
    unsigned long     fAntiDeuteron_Event_Identifier;

    TList     *fHistoList;
    TH2F      *h_Proton_TOF_m2_NoTOFcut;
    TH2F      *h_Deuteron_TOF_m2_NoTOFcut;
    TH2F      *h_AntiProton_TOF_m2_NoTOFcut;
    TH2F      *h_AntiDeuteron_TOF_m2_NoTOFcut;
    TH2F      *h_Proton_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_Deuteron_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_AntiProton_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut;






    ClassDef(AliAnalysisTask_pd_CreateTrees_PairsOnly,1);

};






#endif
