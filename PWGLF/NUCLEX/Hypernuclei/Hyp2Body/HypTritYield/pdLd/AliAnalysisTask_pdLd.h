
#ifndef AliAnalysisTask_pdLd_H
#define AliAnalysisTask_pdLd_H

#include "AliAnalysisTaskSE.h"
#include "TObject.h"

class AliAODEvent;
class AliAODInputHandler;
class AliEventPoolManager;






class AliAnalysisTask_pdLd : public AliAnalysisTaskSE
{

  public:

    AliAnalysisTask_pdLd();
    AliAnalysisTask_pdLd(const char *name,int CollisionSystem);
    AliAnalysisTask_pdLd& operator = (const AliAnalysisTask_pdLd &task);
    AliAnalysisTask_pdLd(const AliAnalysisTask_pdLd &task);
    virtual ~AliAnalysisTask_pdLd();

    void UserCreateOutputObjects();
    void UserExec(Option_t *);
    void Terminate(Option_t *);
    double CalculateBetaTOF(AliAODTrack &track); 
    double CalculateMassSquareTOF(AliAODTrack &track);
    double CalculateRelativeMomentum(TLorentzVector &Pair, TLorentzVector &Particle1, TLorentzVector &Particle2);
    bool ClosePairRejection(double bfield, AliAODTrack &track1, AliAODTrack &track2);
    double RecalculatePhi(double bfield, AliAODTrack &track, double radius);
    double CalculateDeuteronSigmaMassSquareTOF(double pT, double massSq, bool isMatter);
    bool CheckProtonCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter);
    bool CheckDeuteronCuts(AliAODTrack &Track,AliPIDResponse &fPIDResponse, bool isMatter);


  private:

    AliAODEvent		  *fAODEvent;
    AliAODInputHandler	  *fAODHandler; 
    AliAODHeader	  *fHeader;
    AliPIDResponse	  *fPIDResponse;
    AliEventPoolManager	  *fEventPoolManager;
    AliStack		  *fStack;
    TList		  *fHistList_Event;
    TList		  *fHistList_Proton;
    TList		  *fHistList_Deuteron;
    TList		  *fHistList_ProtonDeuteron;
    TList		  *fHistList_AntiProton;
    TList		  *fHistList_AntiDeuteron;
    TList		  *fHistList_AntiProtonAntiDeuteron;

    int			  fCollisionSystem;

    // histograms for event
    TH1F      *fHist_Event_CutCounter;
    TH1F      *fHist_Event_PrimVertexZ;
    TH1F      *fHist_Event_Centrality;

    // histograms for protons
    TH1F      *fHist_Proton_CutCounter;
    TH1F      *fHist_Proton_pT;
    TH1F      *fHist_Proton_p;
    TH1F      *fHist_Proton_pTPC;
    TH1F      *fHist_Proton_Eta;
    TH1F      *fHist_Proton_Phi;
    TH1F      *fHist_Proton_TPC_nCluster;
    TH1F      *fHist_Proton_TPC_CrossedRows;
    TH1F      *fHist_Proton_TPC_RatioRowsCluster;
    TH1F      *fHist_Proton_TPC_SharedCluster;
    TH1F      *fHist_Proton_TPC_Chi2overNDF;
    TH1F      *fHist_Proton_ITS_nCluster;
    TH1F      *fHist_Proton_TOF_Signal;

    TH2F      *fHist_Proton_DCAxy;
    TH2F      *fHist_Proton_DCAz;
    TH2F      *fHist_Proton_TPC_nSigma_pT;
    TH2F      *fHist_Proton_TPC_nSigma_p;
    TH2F      *fHist_Proton_TOF_nSigma_pT;
    TH2F      *fHist_Proton_TOF_nSigma_p;
    TH2F      *fHist_Proton_ITS_nSigma_pT;
    TH2F      *fHist_Proton_ITS_nSigma_p;
    TH2F      *fHist_Proton_TPC_dEdx_pT;
    TH2F      *fHist_Proton_TPC_dEdx_p;
    TH2F      *fHist_Proton_TOF_Beta_pT;
    TH2F      *fHist_Proton_TOF_Beta_p;
    TH2F      *fHist_Proton_TOF_MassSquare_pT;
    TH2F      *fHist_Proton_TOF_MassSquare_p;
    TH2F      *fHist_Proton_ITS_dEdx_pT;
    TH2F      *fHist_Proton_ITS_dEdx_p;

    // histograms for deuterons
    TH1F      *fHist_Deuteron_CutCounter;
    TH1F      *fHist_Deuteron_pT;
    TH1F      *fHist_Deuteron_p;
    TH1F      *fHist_Deuteron_pTPC;
    TH1F      *fHist_Deuteron_Eta;
    TH1F      *fHist_Deuteron_Phi;
    TH1F      *fHist_Deuteron_TPC_nCluster;
    TH1F      *fHist_Deuteron_TPC_CrossedRows;
    TH1F      *fHist_Deuteron_TPC_RatioRowsCluster;
    TH1F      *fHist_Deuteron_TPC_SharedCluster;
    TH1F      *fHist_Deuteron_TPC_Chi2overNDF;
    TH1F      *fHist_Deuteron_ITS_nCluster;
    TH1F      *fHist_Deuteron_TOF_Signal;

    TH2F      *fHist_Deuteron_DCAxy;
    TH2F      *fHist_Deuteron_DCAz;
    TH2F      *fHist_Deuteron_TPC_nSigma_pT;
    TH2F      *fHist_Deuteron_TPC_nSigma_p;
    TH2F      *fHist_Deuteron_TOF_nSigma_pT;
    TH2F      *fHist_Deuteron_TOF_nSigma_p;
    TH2F      *fHist_Deuteron_TOF_nSigmaMassSq_pT;
    TH2F      *fHist_Deuteron_TOF_nSigmaMassSq_p;
    TH2F      *fHist_Deuteron_ITS_nSigma_pT;
    TH2F      *fHist_Deuteron_ITS_nSigma_p;
    TH2F      *fHist_Deuteron_TPC_dEdx_pT;
    TH2F      *fHist_Deuteron_TPC_dEdx_p;
    TH2F      *fHist_Deuteron_TOF_Beta_pT;
    TH2F      *fHist_Deuteron_TOF_Beta_p;
    TH2F      *fHist_Deuteron_TOF_MassSquare_pT;
    TH2F      *fHist_Deuteron_TOF_MassSquare_p;
    TH2F      *fHist_Deuteron_ITS_dEdx_pT;
    TH2F      *fHist_Deuteron_ITS_dEdx_p;
    
    // histograms for p-d pairs
    TH1F      *fHist_ProtonDeuteron_CutCounter;
    TH1F      *fHist_ProtonDeuteron_SED;
    TH1F      *fHist_ProtonDeuteron_MED;
    TH1F      *fHist_ProtonDeuteron_RPD;
    TH1F      *fHist_ProtonDeuteron_PairsPerEvent;
    TH1F      *fHist_ProtonDeuteron_EventsForMixing;
    TH2F      *fHist_ProtonDeuteron_AngleOfPairs;
    TH2F      *fHist_ProtonDeuteron_PairMultiplicity;
    TH2F      *fHist_ProtonDeuteron_pT; 
    TH2F      *fHist_ProtonDeuteron_Eta; 
    TH2F      *fHist_ProtonDeuteron_Centrality; 
    TH2F      *fHist_ProtonDeuteron_VertexZ; 
    TH2F      *fHist_ProtonDeuteron_UsedEventsInPool;


    // histograms for antiprotons
    TH1F      *fHist_AntiProton_CutCounter;
    TH1F      *fHist_AntiProton_pT;
    TH1F      *fHist_AntiProton_p;
    TH1F      *fHist_AntiProton_pTPC;
    TH1F      *fHist_AntiProton_Eta;
    TH1F      *fHist_AntiProton_Phi;
    TH1F      *fHist_AntiProton_TPC_nCluster;
    TH1F      *fHist_AntiProton_TPC_CrossedRows;
    TH1F      *fHist_AntiProton_TPC_RatioRowsCluster;
    TH1F      *fHist_AntiProton_TPC_SharedCluster;
    TH1F      *fHist_AntiProton_TPC_Chi2overNDF;
    TH1F      *fHist_AntiProton_ITS_nCluster;
    TH1F      *fHist_AntiProton_TOF_Signal;

    TH2F      *fHist_AntiProton_DCAxy;
    TH2F      *fHist_AntiProton_DCAz;
    TH2F      *fHist_AntiProton_TPC_nSigma_pT;
    TH2F      *fHist_AntiProton_TPC_nSigma_p;
    TH2F      *fHist_AntiProton_TOF_nSigma_pT;
    TH2F      *fHist_AntiProton_TOF_nSigma_p;
    TH2F      *fHist_AntiProton_ITS_nSigma_pT;
    TH2F      *fHist_AntiProton_ITS_nSigma_p;
    TH2F      *fHist_AntiProton_TPC_dEdx_pT;
    TH2F      *fHist_AntiProton_TPC_dEdx_p;
    TH2F      *fHist_AntiProton_TOF_Beta_pT;
    TH2F      *fHist_AntiProton_TOF_Beta_p;
    TH2F      *fHist_AntiProton_TOF_MassSquare_pT;
    TH2F      *fHist_AntiProton_TOF_MassSquare_p;
    TH2F      *fHist_AntiProton_ITS_dEdx_pT;
    TH2F      *fHist_AntiProton_ITS_dEdx_p;

    // histograms for deuterons
    TH1F      *fHist_AntiDeuteron_CutCounter;
    TH1F      *fHist_AntiDeuteron_pT;
    TH1F      *fHist_AntiDeuteron_p;
    TH1F      *fHist_AntiDeuteron_pTPC;
    TH1F      *fHist_AntiDeuteron_Eta;
    TH1F      *fHist_AntiDeuteron_Phi;
    TH1F      *fHist_AntiDeuteron_TPC_nCluster;
    TH1F      *fHist_AntiDeuteron_TPC_CrossedRows;
    TH1F      *fHist_AntiDeuteron_TPC_RatioRowsCluster;
    TH1F      *fHist_AntiDeuteron_TPC_SharedCluster;
    TH1F      *fHist_AntiDeuteron_TPC_Chi2overNDF;
    TH1F      *fHist_AntiDeuteron_ITS_nCluster;
    TH1F      *fHist_AntiDeuteron_TOF_Signal;

    TH2F      *fHist_AntiDeuteron_DCAxy;
    TH2F      *fHist_AntiDeuteron_DCAz;
    TH2F      *fHist_AntiDeuteron_TPC_nSigma_pT;
    TH2F      *fHist_AntiDeuteron_TPC_nSigma_p;
    TH2F      *fHist_AntiDeuteron_TOF_nSigma_pT;
    TH2F      *fHist_AntiDeuteron_TOF_nSigma_p;
    TH2F      *fHist_AntiDeuteron_TOF_nSigmaMassSq_pT;
    TH2F      *fHist_AntiDeuteron_TOF_nSigmaMassSq_p;
    TH2F      *fHist_AntiDeuteron_ITS_nSigma_pT;
    TH2F      *fHist_AntiDeuteron_ITS_nSigma_p;
    TH2F      *fHist_AntiDeuteron_TPC_dEdx_pT;
    TH2F      *fHist_AntiDeuteron_TPC_dEdx_p;
    TH2F      *fHist_AntiDeuteron_TOF_Beta_pT;
    TH2F      *fHist_AntiDeuteron_TOF_Beta_p;
    TH2F      *fHist_AntiDeuteron_TOF_MassSquare_pT;
    TH2F      *fHist_AntiDeuteron_TOF_MassSquare_p;
    TH2F      *fHist_AntiDeuteron_ITS_dEdx_pT;
    TH2F      *fHist_AntiDeuteron_ITS_dEdx_p;
    
    // histograms for ap-ad pairs
    TH1F      *fHist_AntiProtonAntiDeuteron_CutCounter;
    TH1F      *fHist_AntiProtonAntiDeuteron_SED;
    TH1F      *fHist_AntiProtonAntiDeuteron_MED;
    TH1F      *fHist_AntiProtonAntiDeuteron_RPD;
    TH1F      *fHist_AntiProtonAntiDeuteron_PairsPerEvent;
    TH1F      *fHist_AntiProtonAntiDeuteron_EventsForMixing;
    TH2F      *fHist_AntiProtonAntiDeuteron_AngleOfPairs;
    TH2F      *fHist_AntiProtonAntiDeuteron_PairMultiplicity;
    TH2F      *fHist_AntiProtonAntiDeuteron_pT; 
    TH2F      *fHist_AntiProtonAntiDeuteron_Eta; 
    TH2F      *fHist_AntiProtonAntiDeuteron_Centrality; 
    TH2F      *fHist_AntiProtonAntiDeuteron_VertexZ; 
    TH2F      *fHist_AntiProtonAntiDeuteron_UsedEventsInPool;


    std::vector<int>	*ProtonTrackArray;
    std::vector<int>	*DeuteronTrackArray;
    std::vector<int>	*Lambdav0Array;
    std::vector<int>	*AntiProtonTrackArray;
    std::vector<int>	*AntiDeuteronTrackArray;
    std::vector<int>	*AntiLambdav0Array;



    ClassDef(AliAnalysisTask_pdLd,1);

};







class AliAODTrackTiny : public TObject
{

  public:

    AliAODTrackTiny() : x(-999),y(-999),z(-999){};
    virtual ~AliAODTrackTiny(){};

    void InitFromTrack(const AliAODTrack *Track, double Centrality, double PrimaryVertexZ){
      if(!Track)
      {
	std::cout << "WARNING: Input track for AliAODTrackTiny does not exist" << std::endl;
	return;
      }

      x = Track->Px();
      y = Track->Py();
      z = Track->Pz();
      centrality = Centrality;
      vertex = PrimaryVertexZ;

    }

    double Px() const { return x;}
    double Py() const { return y;}
    double Pz() const { return z;}
    double GetPrimaryVertexZ() const { return vertex;}
    double GetCentrality() const { return centrality;}

  private:

    double x;
    double y;
    double z;
    double vertex;
    double centrality;


  ClassDef(AliAODTrackTiny,1);

};

#endif
