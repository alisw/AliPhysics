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
    bool CheckLambdaCuts(AliAODv0 &v0,double PrimaryVertexPos[3], AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber, double MassInvLambda, double MassInvWrongLambda, double MassInvKaonShort, double MassInvPhoton);
    bool CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter, int RunNumber);
    double CalculateSigmadEdxITS(AliAODTrack &Track, int ParticleSpecies, int RunNumber);
    double CalculateInvariantMassLambda(double Momentum1[3], double Momentum2[3], int WhichMassHypothesis);


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
    Float_t fLambda_px;
    Float_t fLambda_py;
    Float_t fLambda_pz;
    Float_t fLambda_px_Generated;
    Float_t fLambda_py_Generated;
    Float_t fLambda_pz_Generated;
    Float_t fLambda_Eta;
    Float_t fLambda_Phi;
    Float_t fLambda_TransverseRadius;
    Float_t fLambda_CosinePointingAngle;
    Float_t fLambda_DCAv0ToPrimaryVertex;
    Float_t fLambda_DCAv0Daughters;
    Float_t fLambda_Alpha;
    Float_t fLambda_qT;
    Float_t fLambda_DecayLength;
    Float_t fLambda_OpenAngle;
    Int_t	  fLambda_PDG_Daughter1;
    Int_t	  fLambda_PDG_Daughter2;
    Int_t	  fLambda_PDG_v01;
    Int_t	  fLambda_PDG_v02;
    Int_t	  fLambda_PDG_Mother1;
    Int_t	  fLambda_PDG_Mother2;
    Bool_t  fLambda_SameV0;
    Float_t fLambda_Event_Centrality;
    Float_t fLambda_Event_PrimaryVertexZ;
    Bool_t fLambda_Event_BField;
    UInt_t fLambda_Event_Multiplicity;
    ULong64_t fLambda_Event_Identifier;
    Bool_t fLambda_Event_IsFirstParticle;

    Float_t fLambda_Daughter_Proton_px;
    Float_t fLambda_Daughter_Proton_py;
    Float_t fLambda_Daughter_Proton_pz;
    Float_t fLambda_Daughter_Proton_px_Generated;
    Float_t fLambda_Daughter_Proton_py_Generated;
    Float_t fLambda_Daughter_Proton_pz_Generated;
    Float_t fLambda_Daughter_Proton_px_DecayVertex;
    Float_t fLambda_Daughter_Proton_py_DecayVertex;
    Float_t fLambda_Daughter_Proton_pz_DecayVertex;
    Float_t fLambda_Daughter_Proton_pTPC;
    Float_t fLambda_Daughter_Proton_Eta;
    Float_t fLambda_Daughter_Proton_Phi;
    Float_t fLambda_Daughter_Proton_TPC_Chi2;
    Float_t fLambda_Daughter_Proton_TPC_dEdx;
    Float_t fLambda_Daughter_Proton_TPC_dEdx_nSigma;
    Float_t fLambda_Daughter_Proton_TOF_Mass2;
    Float_t fLambda_Daughter_Proton_TOF_Mass2_nSigma;
    Float_t fLambda_Daughter_Proton_ITS_dEdx;
    Float_t fLambda_Daughter_Proton_ITS_dEdx_nSigma;
    Float_t fLambda_Daughter_Proton_DCAxy;
    Float_t fLambda_Daughter_Proton_DCAz;
    
    UShort_t fLambda_Daughter_Proton_TPC_nCrossedRows;
    UShort_t fLambda_Daughter_Proton_TPC_nSharedCluster;
    UShort_t fLambda_Daughter_Proton_TPC_nFindableCluster;
    UShort_t fLambda_Daughter_Proton_TPC_nCluster;
    UShort_t fLambda_Daughter_Proton_ITS_nCluster;

    Float_t fLambda_Daughter_AntiPion_px;
    Float_t fLambda_Daughter_AntiPion_py;
    Float_t fLambda_Daughter_AntiPion_pz;
    Float_t fLambda_Daughter_AntiPion_px_Generated;
    Float_t fLambda_Daughter_AntiPion_py_Generated;
    Float_t fLambda_Daughter_AntiPion_pz_Generated;
    Float_t fLambda_Daughter_AntiPion_px_DecayVertex;
    Float_t fLambda_Daughter_AntiPion_py_DecayVertex;
    Float_t fLambda_Daughter_AntiPion_pz_DecayVertex;
    Float_t fLambda_Daughter_AntiPion_pTPC;
    Float_t fLambda_Daughter_AntiPion_Eta;
    Float_t fLambda_Daughter_AntiPion_Phi;
    Float_t fLambda_Daughter_AntiPion_TPC_Chi2;
    Float_t fLambda_Daughter_AntiPion_TPC_dEdx;
    Float_t fLambda_Daughter_AntiPion_TPC_dEdx_nSigma;
    Float_t fLambda_Daughter_AntiPion_TOF_Mass2;
    Float_t fLambda_Daughter_AntiPion_TOF_Mass2_nSigma;
    Float_t fLambda_Daughter_AntiPion_ITS_dEdx;
    Float_t fLambda_Daughter_AntiPion_ITS_dEdx_nSigma;
    Float_t fLambda_Daughter_AntiPion_DCAxy;
    Float_t fLambda_Daughter_AntiPion_DCAz;

    UShort_t fLambda_Daughter_AntiPion_TPC_nCrossedRows;
    UShort_t fLambda_Daughter_AntiPion_TPC_nSharedCluster;
    UShort_t fLambda_Daughter_AntiPion_TPC_nFindableCluster;
    UShort_t fLambda_Daughter_AntiPion_TPC_nCluster;
    UShort_t fLambda_Daughter_AntiPion_ITS_nCluster;


    TTree     *fSaveTree_Deuteron;
    Float_t     fDeuteron_px;
    Float_t     fDeuteron_py;
    Float_t     fDeuteron_pz;
    Float_t     fDeuteron_px_Generated;
    Float_t     fDeuteron_py_Generated;
    Float_t     fDeuteron_pz_Generated;
    Float_t     fDeuteron_pTPC;
    Float_t     fDeuteron_Eta;
    Float_t     fDeuteron_Phi;
    Float_t     fDeuteron_TPC_Chi2;
    Float_t     fDeuteron_TPC_dEdx;
    Float_t     fDeuteron_TPC_dEdx_nSigma;
    Float_t     fDeuteron_TOF_Mass2;
    Float_t     fDeuteron_TOF_Mass2_nSigma;
    Float_t     fDeuteron_ITS_dEdx;
    Float_t     fDeuteron_ITS_dEdx_nSigma;
    Float_t     fDeuteron_DCAxy;
    Float_t     fDeuteron_DCAz;
    Float_t     fDeuteron_Event_Centrality;
    Float_t     fDeuteron_Event_PrimaryVertexZ;
    Bool_t     fDeuteron_Event_BField;
    UShort_t    fDeuteron_TPC_nCrossedRows;
    UShort_t    fDeuteron_TPC_nSharedCluster;
    UShort_t    fDeuteron_TPC_nFindableCluster;
    UShort_t    fDeuteron_TPC_nCluster;
    UShort_t    fDeuteron_ITS_nCluster;
    Int_t		      fDeuteron_PDG;
    Int_t		      fDeuteron_MotherPDG;
    Bool_t	      fDeuteron_FilterBit;
    UInt_t      fDeuteron_Event_Multiplicity;
    ULong64_t     fDeuteron_Event_Identifier;
    Bool_t	      fDeuteron_ITS_Layer0;
    Bool_t	      fDeuteron_ITS_Layer1;
    Bool_t	      fDeuteron_ITS_Layer2;
    Bool_t	      fDeuteron_ITS_Layer3;
    Bool_t	      fDeuteron_ITS_Layer4;
    Bool_t	      fDeuteron_ITS_Layer5;
    Bool_t fDeuteron_Event_IsFirstParticle;

    TTree *fSaveTree_AntiLambda;
    Float_t fAntiLambda_px;
    Float_t fAntiLambda_py;
    Float_t fAntiLambda_pz;
    Float_t fAntiLambda_px_Generated;
    Float_t fAntiLambda_py_Generated;
    Float_t fAntiLambda_pz_Generated;
    Float_t fAntiLambda_Eta;
    Float_t fAntiLambda_Phi;
    Float_t fAntiLambda_TransverseRadius;
    Float_t fAntiLambda_CosinePointingAngle;
    Float_t fAntiLambda_DCAv0ToPrimaryVertex;
    Float_t fAntiLambda_DCAv0Daughters;
    Float_t fAntiLambda_Alpha;
    Float_t fAntiLambda_qT;
    Float_t fAntiLambda_DecayLength;
    Float_t fAntiLambda_OpenAngle;
    Int_t	  fAntiLambda_PDG_Daughter1;
    Int_t	  fAntiLambda_PDG_Daughter2;
    Int_t	  fAntiLambda_PDG_v01;
    Int_t	  fAntiLambda_PDG_v02;
    Int_t	  fAntiLambda_PDG_Mother1;
    Int_t	  fAntiLambda_PDG_Mother2;
    Bool_t  fAntiLambda_SameV0;
    Float_t fAntiLambda_Event_Centrality;
    Float_t fAntiLambda_Event_PrimaryVertexZ;
    Bool_t fAntiLambda_Event_BField;
    UInt_t fAntiLambda_Event_Multiplicity;
    ULong64_t fAntiLambda_Event_Identifier;
    Bool_t fAntiLambda_Event_IsFirstParticle;

    Float_t fAntiLambda_Daughter_AntiProton_px;
    Float_t fAntiLambda_Daughter_AntiProton_py;
    Float_t fAntiLambda_Daughter_AntiProton_pz;
    Float_t fAntiLambda_Daughter_AntiProton_px_Generated;
    Float_t fAntiLambda_Daughter_AntiProton_py_Generated;
    Float_t fAntiLambda_Daughter_AntiProton_pz_Generated;
    Float_t fAntiLambda_Daughter_AntiProton_px_DecayVertex;
    Float_t fAntiLambda_Daughter_AntiProton_py_DecayVertex;
    Float_t fAntiLambda_Daughter_AntiProton_pz_DecayVertex;
    Float_t fAntiLambda_Daughter_AntiProton_pTPC;
    Float_t fAntiLambda_Daughter_AntiProton_Eta;
    Float_t fAntiLambda_Daughter_AntiProton_Phi;
    Float_t fAntiLambda_Daughter_AntiProton_TPC_Chi2;
    Float_t fAntiLambda_Daughter_AntiProton_TPC_dEdx;
    Float_t fAntiLambda_Daughter_AntiProton_TPC_dEdx_nSigma;
    Float_t fAntiLambda_Daughter_AntiProton_TOF_Mass2;
    Float_t fAntiLambda_Daughter_AntiProton_TOF_Mass2_nSigma;
    Float_t fAntiLambda_Daughter_AntiProton_ITS_dEdx;
    Float_t fAntiLambda_Daughter_AntiProton_ITS_dEdx_nSigma;
    Float_t fAntiLambda_Daughter_AntiProton_DCAxy;
    Float_t fAntiLambda_Daughter_AntiProton_DCAz;

    UShort_t fAntiLambda_Daughter_AntiProton_TPC_nCrossedRows;
    UShort_t fAntiLambda_Daughter_AntiProton_TPC_nSharedCluster;
    UShort_t fAntiLambda_Daughter_AntiProton_TPC_nFindableCluster;
    UShort_t fAntiLambda_Daughter_AntiProton_TPC_nCluster;
    UShort_t fAntiLambda_Daughter_AntiProton_ITS_nCluster;

    Float_t fAntiLambda_Daughter_Pion_px;
    Float_t fAntiLambda_Daughter_Pion_py;
    Float_t fAntiLambda_Daughter_Pion_pz;
    Float_t fAntiLambda_Daughter_Pion_px_Generated;
    Float_t fAntiLambda_Daughter_Pion_py_Generated;
    Float_t fAntiLambda_Daughter_Pion_pz_Generated;
    Float_t fAntiLambda_Daughter_Pion_px_DecayVertex;
    Float_t fAntiLambda_Daughter_Pion_py_DecayVertex;
    Float_t fAntiLambda_Daughter_Pion_pz_DecayVertex;
    Float_t fAntiLambda_Daughter_Pion_pTPC;
    Float_t fAntiLambda_Daughter_Pion_Eta;
    Float_t fAntiLambda_Daughter_Pion_Phi;
    Float_t fAntiLambda_Daughter_Pion_TPC_Chi2;
    Float_t fAntiLambda_Daughter_Pion_TPC_dEdx;
    Float_t fAntiLambda_Daughter_Pion_TPC_dEdx_nSigma;
    Float_t fAntiLambda_Daughter_Pion_TOF_Mass2;
    Float_t fAntiLambda_Daughter_Pion_TOF_Mass2_nSigma;
    Float_t fAntiLambda_Daughter_Pion_ITS_dEdx;
    Float_t fAntiLambda_Daughter_Pion_ITS_dEdx_nSigma;
    Float_t fAntiLambda_Daughter_Pion_DCAxy;
    Float_t fAntiLambda_Daughter_Pion_DCAz;

    UShort_t fAntiLambda_Daughter_Pion_TPC_nCrossedRows;
    UShort_t fAntiLambda_Daughter_Pion_TPC_nSharedCluster;
    UShort_t fAntiLambda_Daughter_Pion_TPC_nFindableCluster;
    UShort_t fAntiLambda_Daughter_Pion_TPC_nCluster;
    UShort_t fAntiLambda_Daughter_Pion_ITS_nCluster;



    TTree     *fSaveTree_AntiDeuteron;
    Float_t     fAntiDeuteron_px;
    Float_t     fAntiDeuteron_py;
    Float_t     fAntiDeuteron_pz;
    Float_t     fAntiDeuteron_px_Generated;
    Float_t     fAntiDeuteron_py_Generated;
    Float_t     fAntiDeuteron_pz_Generated;
    Float_t     fAntiDeuteron_pTPC;
    Float_t     fAntiDeuteron_Eta;
    Float_t     fAntiDeuteron_Phi;
    Float_t     fAntiDeuteron_TPC_Chi2;
    Float_t     fAntiDeuteron_TPC_dEdx;
    Float_t     fAntiDeuteron_TPC_dEdx_nSigma;
    Float_t     fAntiDeuteron_TOF_Mass2;
    Float_t     fAntiDeuteron_TOF_Mass2_nSigma;
    Float_t     fAntiDeuteron_ITS_dEdx;
    Float_t     fAntiDeuteron_ITS_dEdx_nSigma;
    Float_t     fAntiDeuteron_DCAxy;
    Float_t     fAntiDeuteron_DCAz;
    Float_t     fAntiDeuteron_Event_Centrality;
    Float_t     fAntiDeuteron_Event_PrimaryVertexZ;
    Bool_t     fAntiDeuteron_Event_BField;
    UShort_t    fAntiDeuteron_TPC_nCrossedRows;
    UShort_t    fAntiDeuteron_TPC_nSharedCluster;
    UShort_t    fAntiDeuteron_TPC_nFindableCluster;
    UShort_t    fAntiDeuteron_TPC_nCluster;
    UShort_t    fAntiDeuteron_ITS_nCluster;
    Int_t		      fAntiDeuteron_PDG;
    Int_t		      fAntiDeuteron_MotherPDG;
    Bool_t	      fAntiDeuteron_FilterBit;
    UInt_t      fAntiDeuteron_Event_Multiplicity;
    ULong64_t     fAntiDeuteron_Event_Identifier;
    Bool_t	      fAntiDeuteron_ITS_Layer0;
    Bool_t	      fAntiDeuteron_ITS_Layer1;
    Bool_t	      fAntiDeuteron_ITS_Layer2;
    Bool_t	      fAntiDeuteron_ITS_Layer3;
    Bool_t	      fAntiDeuteron_ITS_Layer4;
    Bool_t	      fAntiDeuteron_ITS_Layer5;
    Bool_t fAntiDeuteron_Event_IsFirstParticle;

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
