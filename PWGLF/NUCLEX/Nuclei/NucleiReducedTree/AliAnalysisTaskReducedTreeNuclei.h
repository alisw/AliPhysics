#ifndef AliAnalysisTaskReducedTreeNuclei_cxx
#define AliAnalysisTaskReducedTreeNuclei_cxx

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliMultSelection;
class AliPIDResponse;
class AliAODEvent;
class AliAODTrack;
class AliAnalysisUtils;
class TList;
class TTree;
class TH1F;

#include "AliEventCuts.h"

class AliAnalysisTaskReducedTreeNuclei : public AliAnalysisTaskSE {
   
public:
   AliAnalysisTaskReducedTreeNuclei();
   AliAnalysisTaskReducedTreeNuclei(const char *name);
   virtual ~AliAnalysisTaskReducedTreeNuclei();
   
   virtual void   UserCreateOutputObjects();
   virtual void   UserExec (Option_t *option);
   
   void FillTritonTree(Bool_t fillTri){ fFillTri=fillTri; }
   void FillHypTritonTree(Bool_t fillHypTri){ fFillHypTri=fillHypTri; }
   
   

   Bool_t   GetInputEvent ();
   Bool_t   PassedBasicTrackQualityCuts (AliAODTrack *track);
   Bool_t   IsHeliumCandidate           (AliAODTrack *track);
   Bool_t   IsTritonCandidate           (AliAODTrack *track);
   Bool_t   IsHyperTritonCandidate      (AliAODTrack *track1,AliAODTrack *track2);
   Double_t    GetDCAxy                    (AliAODTrack *track);
   Double_t    GetDCAz                     (AliAODTrack *track);
   
   virtual void   Terminate(Option_t *);
   
private:
   AliAODEvent    *fAODevent;//!
   AliPIDResponse *fPIDResponse;//!
   AliEventCuts   fAODeventCuts;// Event cuts
   AliAnalysisUtils *fUtils;//!

   AliAODVZERO* fAODVZERO;//!
   // globle varibles
   Bool_t fFillTri;
   Bool_t fFillHypTri;

   TList          *fQAList;//!
   
   // Event Selection Tree
//    TTree *TreeEventSelection;//!

   TList          *fOutputList;//!
   // Event histograms
   TH2D *histoEventSelection; //!
   TH2D *histoEventMultiplicity; //!

   TH1D *histoEventMultV0M_MB; //!
   TH1D *histoEventMultV0M_HM; //!
   
   //Reduced Trees
   TTree *reducedTree_Helium;//!
   TTree *reducedTree_Triton;//!
   TTree *reducedTree_HyperTriton;//!
   
   //Variables (Helium)
   // trigger 
   bool isTrigINT7;//
   bool isTrigHighMult;//
//    Long64_t triggerMask;
   Int_t magFieldSign;//
   
   Int_t SelectionStep;//
   //check for more estimators, e.g. SPD, TPC track multiplicity ...
   Double_t multPercentile_V0M;//
   Double_t multPercentile_V0A;//
   Double_t multPercentile_V0C;//
   Double_t multPercentile_OnlineV0M;//
   Double_t multPercentile_OnlineV0A;//
   Double_t multPercentile_OnlineV0C;//
   Double_t multPercentile_ADM;//
   Double_t multPercentile_ADA;//
   Double_t multPercentile_ADC;//
   Double_t multPercentile_SPDClusters;//
   Double_t multPercentile_SPDTracklets;//
   Double_t multPercentile_RefMult05;//
   Double_t multPercentile_RefMult08;//
   Double_t multPercentile_CL1;//
   Double_t multPercentile_ZNA;//
   
   Int_t Ntrk_V0M;//
   Int_t Ntrk_V0A;//
   Int_t Ntrk_V0C;//
   Int_t Ntrk_OnlineV0M;//
   Int_t Ntrk_OnlineV0A;//
   Int_t Ntrk_OnlineV0C;//
   Int_t Ntrk_ADM;//
   Int_t Ntrk_ADA;//
   Int_t Ntrk_ADC;//
   Int_t Ntrk_SPDClusters;//
   Int_t Ntrk_SPDTracklets;//
   Int_t Ntrk_RefMult05;//
   Int_t Ntrk_RefMult08;//
   
   Int_t nVertexContributors;//
   Double_t xVertex;//
   Double_t yVertex;//
   Double_t zVertex;//
   
   
   Double_t px;//
   Double_t py;//
   Double_t pz;//
   
   Double_t TPCmomentum;//
   Double_t TRDmomentum;//
   Double_t integratedLength;//
   Double_t timeOfFlight;//
   Double_t beta;//
   Double_t gamma;//
   Double_t mass;//
   
   Int_t trackID;//
   Double_t eta;//
   Double_t phi;//
   Double_t theta;//
   Double_t y;//
   Int_t q;//
   Double_t dcaxy;//
   Double_t dcaz;//
   
   Int_t nTPC_Clusters;//
   Int_t nTRD_Clusters;//
   Int_t nITS_Clusters;//
   Int_t nTPC_FindableClusters;//
   Int_t nTPC_CrossedRows;//
   Int_t nTPC_Clusters_dEdx;//
   bool HasPointOnITSLayer0;//
   bool HasPointOnITSLayer1;//
   bool HasPointOnITSLayer2;//
   bool HasPointOnITSLayer3;//
   bool HasPointOnITSLayer4;//
   bool HasPointOnITSLayer5;//
   bool HasSharedPointOnITSLayer0;//
   bool HasSharedPointOnITSLayer1;//
   bool HasSharedPointOnITSLayer2;//
   bool HasSharedPointOnITSLayer3;//
   bool HasSharedPointOnITSLayer4;//
   bool HasSharedPointOnITSLayer5;//
   Double_t chi2_TPC;//check
   Double_t chi2_NDF;//check
   Double_t chi2_ITS;//check
   
   //ITSsignal, TPCsignal,TOFsignal,TRDsignal,HMPID signal
   Double_t ITSsignal;
   Double_t TPCsignal;
   Double_t TOFsignal;
   Double_t TRDsignal;
   Double_t HMPIDsignal;
   
   Double_t nSigmaITS_He4;//
   Double_t nSigmaTPC_He4;//
   Double_t nSigmaTOF_He4;//
   Double_t nSigmaTRD_He4;//
   Double_t nSigmaHMPID_He4;//
   
   Double_t nSigmaITS_He3;//
   Double_t nSigmaTPC_He3;//
   Double_t nSigmaTOF_He3;//
   Double_t nSigmaTRD_He3;//
   Double_t nSigmaHMPID_He3;//
   
   Double_t nSigmaITS_Trit;//
   Double_t nSigmaTPC_Trit;//
   Double_t nSigmaTOF_Trit;//
   Double_t nSigmaTRD_Trit;//
   Double_t nSigmaHMPID_Trit;//
   
   Double_t nSigmaITS_Deut;//
   Double_t nSigmaTPC_Deut;//
   Double_t nSigmaTOF_Deut;//
   Double_t nSigmaTRD_Deut;//
   Double_t nSigmaHMPID_Deut;//
   
   Double_t nSigmaITS_Prot;//
   Double_t nSigmaTPC_Prot;//
   Double_t nSigmaTOF_Prot;//
   Double_t nSigmaTRD_Prot;//
   Double_t nSigmaHMPID_Prot;//
   
   Double_t nSigmaITS_Pion;//
   Double_t nSigmaTPC_Pion;//
   Double_t nSigmaTOF_Pion;//
   Double_t nSigmaTRD_Pion;//
   Double_t nSigmaHMPID_Pion;//
   
   Double_t nSigmaITS_Kaon;//
   Double_t nSigmaTPC_Kaon;//
   Double_t nSigmaTOF_Kaon;//
   Double_t nSigmaTRD_Kaon;//
   Double_t nSigmaHMPID_Kaon;//
   
   Double_t nSigmaITS_Elec;//
   Double_t nSigmaTPC_Elec;//
   Double_t nSigmaTOF_Elec;//
   Double_t nSigmaTRD_Elec;//
   Double_t nSigmaHMPID_Elec;//
   
   //Variables for HyperTriton - First Daughter
   Double_t px_Daughter1;//
   Double_t py_Daughter1;//
   Double_t pz_Daughter1;//
   
   Double_t TPCmomentum_Daughter1;//
   
   Int_t trackID_Daughter1;//
   Double_t eta_Daughter1;//
   Double_t phi_Daughter1;//
   Double_t theta_Daughter1;//
   Double_t y_Daughter1;//
   Int_t q_Daughter1;//
   Double_t dcaxy_Daughter1;//
   Double_t dcaz_Daughter1;//
   
   Int_t nTPC_Clusters_Daughter1;//
   Int_t nTPC_FindableClusters_Daughter1;//
   Int_t nTPC_CrossedRows_Daughter1;//
   Int_t nTPC_Clusters_dEdx_Daughter1;//
   Int_t nTRD_Clusters_Daughter1;
   Int_t nITS_Clusters_Daughter1;//
   bool HasPointOnITSLayer0_Daughter1;//
   bool HasPointOnITSLayer1_Daughter1;//
   bool HasPointOnITSLayer2_Daughter1;//
   bool HasPointOnITSLayer3_Daughter1;//
   bool HasPointOnITSLayer4_Daughter1;//
   bool HasPointOnITSLayer5_Daughter1;//
   bool HasSharedPointOnITSLayer0_Daughter1;//
   bool HasSharedPointOnITSLayer1_Daughter1;//
   bool HasSharedPointOnITSLayer2_Daughter1;//
   bool HasSharedPointOnITSLayer3_Daughter1;//
   bool HasSharedPointOnITSLayer4_Daughter1;//
   bool HasSharedPointOnITSLayer5_Daughter1;//
   
   Double_t chi2_TPC_Daughter1;//
   Double_t chi2_NDF_Daughter1;//
   Double_t chi2_ITS_Daughter1;//
   
   Double_t nSigmaITS_He3_Daughter1;//
   Double_t nSigmaTPC_He3_Daughter1;//
   Double_t nSigmaTOF_He3_Daughter1;//
   
   Double_t nSigmaITS_Pion_Daughter1;//
   Double_t nSigmaTPC_Pion_Daughter1;//
   Double_t nSigmaTOF_Pion_Daughter1;//
   
   Double_t nSigmaITS_Trit_Daughter1;//
   Double_t nSigmaTPC_Trit_Daughter1;//
   Double_t nSigmaTOF_Trit_Daughter1;//
   
   
   //Variables for HyperTriton - Second Daughter
   Double_t px_Daughter2;//
   Double_t py_Daughter2;//
   Double_t pz_Daughter2;//
   
   Double_t TPCmomentum_Daughter2;//
   
   Int_t trackID_Daughter2;//
   Double_t eta_Daughter2;//
   Double_t phi_Daughter2;//
   Double_t theta_Daughter2;//
   Double_t y_Daughter2;//
   Int_t q_Daughter2;//
   Double_t dcaxy_Daughter2;//
   Double_t dcaz_Daughter2;//
   
   Int_t nTPC_Clusters_Daughter2;//
   Int_t nTPC_FindableClusters_Daughter2;//
   Int_t nTPC_CrossedRows_Daughter2;//
   Int_t nTPC_Clusters_dEdx_Daughter2;//
   Int_t nTRD_Clusters_Daughter2;
   Int_t nITS_Clusters_Daughter2;//
   bool HasPointOnITSLayer0_Daughter2;//
   bool HasPointOnITSLayer1_Daughter2;//
   bool HasPointOnITSLayer2_Daughter2;//
   bool HasPointOnITSLayer3_Daughter2;//
   bool HasPointOnITSLayer4_Daughter2;//
   bool HasPointOnITSLayer5_Daughter2;//
   bool HasSharedPointOnITSLayer0_Daughter2;//
   bool HasSharedPointOnITSLayer1_Daughter2;//
   bool HasSharedPointOnITSLayer2_Daughter2;//
   bool HasSharedPointOnITSLayer3_Daughter2;//
   bool HasSharedPointOnITSLayer4_Daughter2;//
   bool HasSharedPointOnITSLayer5_Daughter2;//
   
   Double_t chi2_TPC_Daughter2;//
   Double_t chi2_NDF_Daughter2;//
   Double_t chi2_ITS_Daughter2;//
   
   Double_t nSigmaITS_He3_Daughter2;//
   Double_t nSigmaTPC_He3_Daughter2;//
   Double_t nSigmaTOF_He3_Daughter2;//
   
   Double_t nSigmaITS_Pion_Daughter2;//
   Double_t nSigmaTPC_Pion_Daughter2;//
   Double_t nSigmaTOF_Pion_Daughter2;//
   
   Double_t nSigmaITS_Trit_Daughter2;//
   Double_t nSigmaTPC_Trit_Daughter2;//
   Double_t nSigmaTOF_Trit_Daughter2;//
   
   //Pair Variables
   Double_t cosPointingAngle;//
   Double_t dcaV0Daughters;//
   Double_t radius;//
   Double_t DecayLength;
   Double_t massV0;//
   Double_t alphaV0;//
   Double_t qtV0;//
   
   
   AliAnalysisTaskReducedTreeNuclei(const AliAnalysisTaskReducedTreeNuclei&);
   AliAnalysisTaskReducedTreeNuclei& operator=(const AliAnalysisTaskReducedTreeNuclei&);
   
  ClassDef(AliAnalysisTaskReducedTreeNuclei, 4);
};
#endif
