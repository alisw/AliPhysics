#ifndef AliAnalysisTask_Ld_CreateTrees_PairsOnly_H
#define AliAnalysisTask_Ld_CreateTrees_PairsOnly_H

#include "AliAnalysisTaskSE.h"
#include "TObject.h"

class AliAODEvent;
class AliAODInputHandler;







class AliAnalysisTask_Ld_CreateTrees_PairsOnly : public AliAnalysisTaskSE
{

  public:

    AliAnalysisTask_Ld_CreateTrees_PairsOnly();
    AliAnalysisTask_Ld_CreateTrees_PairsOnly(const char *name,int CollisionSystem, bool UseOpenCuts, bool isMC, bool SavePairsOnly);
    AliAnalysisTask_Ld_CreateTrees_PairsOnly& operator = (const AliAnalysisTask_Ld_CreateTrees_PairsOnly &task);
    AliAnalysisTask_Ld_CreateTrees_PairsOnly(const AliAnalysisTask_Ld_CreateTrees_PairsOnly &task);
    virtual ~AliAnalysisTask_Ld_CreateTrees_PairsOnly();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    double CalculateBetaTOF(AliAODTrack &track); 
    double CalculateMassSquareTOF(AliAODTrack &track);
    double CalculateSigmaMassSquareTOF(double pT, double massSq, int ParticleSpecies, int RunNumber);
    bool CheckProtonCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber);
    bool CheckPionCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, double Pion_DCAxy_min, double Pion_DCAz_min);
    bool CheckLambdaCuts(AliAODv0 &v0,double PrimaryVertexPos[3], AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, double MassInvariant, double WrongMassInvariant);
    bool CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber);
    double CalculateSigmadEdxITS(AliAODTrack &Track, int ParticleSpecies, int RunNumber);
    double CalculateInvariantMassLambda(double Momentum1[3], double Momentum2[3], bool isMatter);


  private:

    AliAODEvent		  *fAODEvent;
    AliAODInputHandler	  *fAODHandler; 
    AliAODHeader	  *fHeader;
    AliPIDResponse	  *fPIDResponse;

    int			  fCollisionSystem;
    bool		  fUseOpenCuts;
    bool		  fIsMC;
    bool		  fSavePairsOnly;

    TTree *fSaveTree_Lambda;
    float fLambda_px;
    float fLambda_py;
    float fLambda_pz;
    float fLambda_px_Generated;
    float fLambda_py_Generated;
    float fLambda_pz_Generated;
    float fLambda_Eta;
    float fLambda_Phi;
    float fLambda_TransverseRadius;
    float fLambda_CosinePointingAngle;
    float fLambda_DCAv0ToPrimaryVertex;
    float fLambda_DCAv0Daughters;
    float fLambda_Alpha;
    float fLambda_qT;
    float fLambda_DecayLength;
    float fLambda_OpenAngle;
    int	  fLambda_PDG_Daughter1;
    int	  fLambda_PDG_Daughter2;
    int	  fLambda_PDG_v01;
    int	  fLambda_PDG_v02;
    int	  fLambda_PDG_Mother1;
    int	  fLambda_PDG_Mother2;
    bool  fLambda_SameV0;
    float fLambda_Event_Centrality;
    float fLambda_Event_PrimaryVertexZ;
    float fLambda_Event_BField;
    unsigned int fLambda_Event_Multiplicity;
    unsigned long fLambda_Event_Identifier;
    bool fLambda_Event_IsFirstParticle;
    unsigned short fLambda_Event_RandomCrossCheckNumber;

    float fLambda_Daughter_Proton_px;
    float fLambda_Daughter_Proton_py;
    float fLambda_Daughter_Proton_pz;
    float fLambda_Daughter_Proton_px_Generated;
    float fLambda_Daughter_Proton_py_Generated;
    float fLambda_Daughter_Proton_pz_Generated;
    float fLambda_Daughter_Proton_px_DecayVertex;
    float fLambda_Daughter_Proton_py_DecayVertex;
    float fLambda_Daughter_Proton_pz_DecayVertex;
    float fLambda_Daughter_Proton_pTPC;
    float fLambda_Daughter_Proton_Eta;
    float fLambda_Daughter_Proton_Phi;
    float fLambda_Daughter_Proton_TPC_Chi2;
    float fLambda_Daughter_Proton_TPC_dEdx;
    float fLambda_Daughter_Proton_TPC_dEdx_nSigma;
    float fLambda_Daughter_Proton_TOF_Mass2;
    float fLambda_Daughter_Proton_TOF_Mass2_nSigma;
    float fLambda_Daughter_Proton_ITS_dEdx;
    float fLambda_Daughter_Proton_ITS_dEdx_nSigma;
    float fLambda_Daughter_Proton_DCAxy;
    float fLambda_Daughter_Proton_DCAz;
    
    unsigned short fLambda_Daughter_Proton_TPC_nCrossedRows;
    unsigned short fLambda_Daughter_Proton_TPC_nSharedCluster;
    unsigned short fLambda_Daughter_Proton_TPC_nFindableCluster;
    unsigned short fLambda_Daughter_Proton_TPC_nCluster;
    unsigned short fLambda_Daughter_Proton_ITS_nCluster;
    unsigned int fLambda_Daughter_Proton_ID;

    float fLambda_Daughter_AntiPion_px;
    float fLambda_Daughter_AntiPion_py;
    float fLambda_Daughter_AntiPion_pz;
    float fLambda_Daughter_AntiPion_px_Generated;
    float fLambda_Daughter_AntiPion_py_Generated;
    float fLambda_Daughter_AntiPion_pz_Generated;
    float fLambda_Daughter_AntiPion_px_DecayVertex;
    float fLambda_Daughter_AntiPion_py_DecayVertex;
    float fLambda_Daughter_AntiPion_pz_DecayVertex;
    float fLambda_Daughter_AntiPion_pTPC;
    float fLambda_Daughter_AntiPion_Eta;
    float fLambda_Daughter_AntiPion_Phi;
    float fLambda_Daughter_AntiPion_TPC_Chi2;
    float fLambda_Daughter_AntiPion_TPC_dEdx;
    float fLambda_Daughter_AntiPion_TPC_dEdx_nSigma;
    float fLambda_Daughter_AntiPion_TOF_Mass2;
    float fLambda_Daughter_AntiPion_TOF_Mass2_nSigma;
    float fLambda_Daughter_AntiPion_ITS_dEdx;
    float fLambda_Daughter_AntiPion_ITS_dEdx_nSigma;
    float fLambda_Daughter_AntiPion_DCAxy;
    float fLambda_Daughter_AntiPion_DCAz;

    unsigned short fLambda_Daughter_AntiPion_TPC_nCrossedRows;
    unsigned short fLambda_Daughter_AntiPion_TPC_nSharedCluster;
    unsigned short fLambda_Daughter_AntiPion_TPC_nFindableCluster;
    unsigned short fLambda_Daughter_AntiPion_TPC_nCluster;
    unsigned short fLambda_Daughter_AntiPion_ITS_nCluster;
    unsigned int fLambda_Daughter_AntiPion_ID;


    TTree     *fSaveTree_Deuteron;
    float     fDeuteron_px;
    float     fDeuteron_py;
    float     fDeuteron_pz;
    float     fDeuteron_px_Generated;
    float     fDeuteron_py_Generated;
    float     fDeuteron_pz_Generated;
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
    int		      fDeuteron_PDG;
    int		      fDeuteron_MotherPDG;
    bool	      fDeuteron_FilterBit;
    unsigned int      fDeuteron_ID;
    unsigned int      fDeuteron_Event_Multiplicity;
    unsigned long     fDeuteron_Event_Identifier;
    bool	      fDeuteron_ITS_Layer0;
    bool	      fDeuteron_ITS_Layer1;
    bool	      fDeuteron_ITS_Layer2;
    bool	      fDeuteron_ITS_Layer3;
    bool	      fDeuteron_ITS_Layer4;
    bool	      fDeuteron_ITS_Layer5;
    bool fDeuteron_Event_IsFirstParticle;
    unsigned short fDeuteron_Event_RandomCrossCheckNumber;

    TTree *fSaveTree_AntiLambda;
    float fAntiLambda_px;
    float fAntiLambda_py;
    float fAntiLambda_pz;
    float fAntiLambda_px_Generated;
    float fAntiLambda_py_Generated;
    float fAntiLambda_pz_Generated;
    float fAntiLambda_Eta;
    float fAntiLambda_Phi;
    float fAntiLambda_TransverseRadius;
    float fAntiLambda_CosinePointingAngle;
    float fAntiLambda_DCAv0ToPrimaryVertex;
    float fAntiLambda_DCAv0Daughters;
    float fAntiLambda_Alpha;
    float fAntiLambda_qT;
    float fAntiLambda_DecayLength;
    float fAntiLambda_OpenAngle;
    int	  fAntiLambda_PDG_Daughter1;
    int	  fAntiLambda_PDG_Daughter2;
    int	  fAntiLambda_PDG_v01;
    int	  fAntiLambda_PDG_v02;
    int	  fAntiLambda_PDG_Mother1;
    int	  fAntiLambda_PDG_Mother2;
    bool  fAntiLambda_SameV0;
    float fAntiLambda_Event_Centrality;
    float fAntiLambda_Event_PrimaryVertexZ;
    float fAntiLambda_Event_BField;
    unsigned int fAntiLambda_Event_Multiplicity;
    unsigned long fAntiLambda_Event_Identifier;
    bool fAntiLambda_Event_IsFirstParticle;
    unsigned short fAntiLambda_Event_RandomCrossCheckNumber;

    float fAntiLambda_Daughter_AntiProton_px;
    float fAntiLambda_Daughter_AntiProton_py;
    float fAntiLambda_Daughter_AntiProton_pz;
    float fAntiLambda_Daughter_AntiProton_px_Generated;
    float fAntiLambda_Daughter_AntiProton_py_Generated;
    float fAntiLambda_Daughter_AntiProton_pz_Generated;
    float fAntiLambda_Daughter_AntiProton_px_DecayVertex;
    float fAntiLambda_Daughter_AntiProton_py_DecayVertex;
    float fAntiLambda_Daughter_AntiProton_pz_DecayVertex;
    float fAntiLambda_Daughter_AntiProton_pTPC;
    float fAntiLambda_Daughter_AntiProton_Eta;
    float fAntiLambda_Daughter_AntiProton_Phi;
    float fAntiLambda_Daughter_AntiProton_TPC_Chi2;
    float fAntiLambda_Daughter_AntiProton_TPC_dEdx;
    float fAntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma;
    float fAntiLambda_Daughter_AntiProton_TOF_Mass2;
    float fAntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma;
    float fAntiLambda_Daughter_AntiProton_ITS_dEdx;
    float fAntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma;
    float fAntiLambda_Daughter_AntiProton_DCAxy;
    float fAntiLambda_Daughter_AntiProton_DCAz;

    unsigned short fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows;
    unsigned short fAntiLambda_Daughter_AntiProton_TPC_nSharedCluster;
    unsigned short fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster;
    unsigned short fAntiLambda_Daughter_AntiProton_TPC_nCluster;
    unsigned short fAntiLambda_Daughter_AntiProton_ITS_nCluster;
    unsigned int fAntiLambda_Daughter_AntiProton_ID;

    float fAntiLambda_Daughter_Pion_px;
    float fAntiLambda_Daughter_Pion_py;
    float fAntiLambda_Daughter_Pion_pz;
    float fAntiLambda_Daughter_Pion_px_Generated;
    float fAntiLambda_Daughter_Pion_py_Generated;
    float fAntiLambda_Daughter_Pion_pz_Generated;
    float fAntiLambda_Daughter_Pion_px_DecayVertex;
    float fAntiLambda_Daughter_Pion_py_DecayVertex;
    float fAntiLambda_Daughter_Pion_pz_DecayVertex;
    float fAntiLambda_Daughter_Pion_pTPC;
    float fAntiLambda_Daughter_Pion_Eta;
    float fAntiLambda_Daughter_Pion_Phi;
    float fAntiLambda_Daughter_Pion_TPC_Chi2;
    float fAntiLambda_Daughter_Pion_TPC_dEdx;
    float fAntiLambda_Daughter_Pion_TPC_dEdx_nSigma;
    float fAntiLambda_Daughter_Pion_TOF_Mass2;
    float fAntiLambda_Daughter_Pion_TOF_Mass2_nSigma;
    float fAntiLambda_Daughter_Pion_ITS_dEdx;
    float fAntiLambda_Daughter_Pion_ITS_dEdx_nSigma;
    float fAntiLambda_Daughter_Pion_DCAxy;
    float fAntiLambda_Daughter_Pion_DCAz;

    unsigned short fAntiLambda_Daughter_Pion_TPC_nCrossedRows;
    unsigned short fAntiLambda_Daughter_Pion_TPC_nSharedCluster;
    unsigned short fAntiLambda_Daughter_Pion_TPC_nFindableCluster;
    unsigned short fAntiLambda_Daughter_Pion_TPC_nCluster;
    unsigned short fAntiLambda_Daughter_Pion_ITS_nCluster;
    unsigned int fAntiLambda_Daughter_Pion_ID;



    TTree     *fSaveTree_AntiDeuteron;
    float     fAntiDeuteron_px;
    float     fAntiDeuteron_py;
    float     fAntiDeuteron_pz;
    float     fAntiDeuteron_px_Generated;
    float     fAntiDeuteron_py_Generated;
    float     fAntiDeuteron_pz_Generated;
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
    int		      fAntiDeuteron_PDG;
    int		      fAntiDeuteron_MotherPDG;
    bool	      fAntiDeuteron_FilterBit;
    unsigned int      fAntiDeuteron_ID;
    unsigned int      fAntiDeuteron_Event_Multiplicity;
    unsigned long     fAntiDeuteron_Event_Identifier;
    bool	      fAntiDeuteron_ITS_Layer0;
    bool	      fAntiDeuteron_ITS_Layer1;
    bool	      fAntiDeuteron_ITS_Layer2;
    bool	      fAntiDeuteron_ITS_Layer3;
    bool	      fAntiDeuteron_ITS_Layer4;
    bool	      fAntiDeuteron_ITS_Layer5;
    bool fAntiDeuteron_Event_IsFirstParticle;
    unsigned short fAntiDeuteron_Event_RandomCrossCheckNumber;

    TList     *fHistoList;
    TH2F      *h_Proton_TOF_m2_NoTOFcut;
    TH2F      *h_Deuteron_TOF_m2_NoTOFcut;
    TH2F      *h_AntiProton_TOF_m2_NoTOFcut;
    TH2F      *h_AntiDeuteron_TOF_m2_NoTOFcut;
    TH2F      *h_Proton_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_Deuteron_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_AntiProton_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_AntiDeuteron_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_Pion_TOF_m2_NoTOFcut;
    TH2F      *h_AntiPion_TOF_m2_NoTOFcut;
    TH2F      *h_Pion_ITS_dEdx_NoTOFcutNoITScut;
    TH2F      *h_AntiPion_ITS_dEdx_NoTOFcutNoITScut;





    ClassDef(AliAnalysisTask_Ld_CreateTrees_PairsOnly,1);

};






#endif
