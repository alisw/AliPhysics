#ifndef AliAnalysisTaskCMEv2A_cxx
#define AliAnalysisTaskCMEv2A_cxx

//class TList;
//class TH1F;
//class TH2F;
//class TProfile;
//class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCMEv2A : public AliAnalysisTaskSE
{
  
 public:
  
  AliAnalysisTaskCMEv2A();
  AliAnalysisTaskCMEv2A(const char *name);
  virtual ~AliAnalysisTaskCMEv2A();
  virtual void UserCreateOutputObjects(); // makes TList and histograms for output
  virtual void UserExec(Option_t *option); // event processor
  
  void SetParameters(); // set parameters used in code, e.g. TPC efficiency corrections
  
  float calc3(float xn, float yn, float x2n, float y2n, float M); // 3 particle corr
  float calc4(float xn, float yn, float x2n, float y2n, float M); // 4 particle corr
  
  int debug; // debug level controls amount of output statements
  int GetDebug(){return debug;}
  void SetDebug(const int x){debug = x;}
  
  bool doMC; // doMC determines whether to look for MC event header
  bool GetDoMC(){return doMC;}
  void SetDoMC(const bool x){doMC = x;}
  
  ULong64_t trigger; // trigger selection
  ULong64_t GetTrigger(){return trigger;}
  void SetTrigger(const ULong64_t x){trigger = x;}
  
  bool dopupcut; // dopupcut determines whether to apply pileup cuts
  bool GetDoPileupCut(){return dopupcut;}
  void SetDoPileupCut(const bool x){dopupcut = x;}
  
  bool doeffcorr; // doeffcorr determines whether to apply efficiency corrections
  bool GetDoEffCorrection(){return doeffcorr;}
  void SetDoEffCorrection(const bool x){doeffcorr = x;}
  
  int centhandle; // cenrality definition selection
  int GetCentHandle(){return centhandle;}
  void SetCentHandle(const int x){centhandle = x;}
  
  int fbit; // AOD filter bit selection
  int GetFilterBit(){return fbit;}
  void SetFilterBit(const int x){fbit = x;}
  
  float zvtxcut; // z-vertex selection for collision
  float GetZvtxCut(){return zvtxcut;}
  void SetZvtxCut(const float x){zvtxcut = x;}
  
  float centcut; // centrality restriction for V0M and TRK
  float GetCentCut(){return centcut;}
  void SetCentCut(const float x){centcut = x;}
  
  int nclscut; // ncls cut for all tracks
  int GetNclsCut(){return nclscut;}
  void SetNclsCut(const float x){nclscut = x;}
  
  float dcacutz; // dcaz cut for all tracks
  float GetDCAcutZ(){return dcacutz;}
  void SetDCAcutZ(const float x){dcacutz = x;}
  
  float dcacutxy; // dcaxy cut for all tracks
  float GetDCAcutXY(){return dcacutxy;}
  void SetDCAcutXY(const float x){dcacutxy = x;}
  
  bool dodcacuts; // dodcacuts determines whether to apply the DCA cuts, which will cause errors depending on the filter bit
  bool GetDoDCAcuts(){return dodcacuts;}
  void SetDoDCAcuts(const bool x){dodcacuts = x;}
  
  float outeta; // outer eta value for all Q-vector components
  float GetOutEta(){return outeta;}
  void SetOutEta(const float x){outeta = x;}
  
  float ineta; // inner eta value for gap Q-vector components
  float GetInEta(){return ineta;}
  void SetInEta(const float x){ineta = x;}
  
  float excleta; // inner eta value for gap Q-vector components
  float GetExclEta(){return excleta;}
  void SetExclEta(const float x){excleta = x;}
  
  float ptmin; // minimum pt for Q-vector components
  float GetPtMin(){return ptmin;}
  void SetPtMin(const float x){ptmin = x;}
  
  float ptmax; // maximum pt for Q-vector components
  float GetPtMax(){return ptmax;}
  void SetPtMax(const float x){ptmax = x;}
  
  bool doacuts; // doacuts determines whether to apply all the same tracking cuts for charge asymmetry
  bool GetDoAcuts(){return doacuts;}
  void SetDoAcuts(const bool x){doacuts = x;}
  
  float nspid; // number of sigma for pid cuts
  int GetSigmaPID(){return nspid;}
  void SetSigmaPID(const int x){nspid = x;}
  
  int cbinlo; // lower centrality bin for histogram array
  int GetCentBinLow(){return cbinlo;}
  void SetCentBinLow(const int x){cbinlo = x;}
  
  int cbinhi; // higher centrality bin for histogram array
  int GetCentBinHigh(){return cbinhi;}
  void SetCentBinHigh(const int x){cbinhi = x;}
  
  bool donested; // donested determines whether to look for MC event header
  bool GetDoNested(){return donested;}
  void SetDoNested(const bool x){donested = x;}
  
  bool dopaircut; // dopaircut determines whether to look for MC event header
  bool GetDoPairCut(){return dopaircut;}
  void SetDoPairCut(const bool x){dopaircut = x;}
  
  float centlo; // lower centrality value for differential cumulants
  float GetCentDCLow(){return centlo;}
  void SetCentDCLow(const float x){centlo = x;}
  
  float centhi; // higher centrality bin for differential cumulants
  float GetCentDCHigh(){return centhi;}
  void SetCentDCHigh(const float x){centhi = x;}
  
 private:
  
  float effTPC[50];//!
  float effHYB[50];//!
  
  float MCEF1_fb1[100];//!
  float MCEF1_fb128[100];//!
  float MCEF1_fb272[100];//!
  float MCEF1_dca[100];//!
  float MCEF3_fb1[100];//!
  float MCEF3_fb128[100];//!
  float MCEF3_fb272[100];//!
  float MCEF3_dca[100];//!
  
  TList       *fOutputList;//!
  TH1F        *fHistDCACOVerrors;//!
  TH1F        *fHistPt;//!
  TH1F        *fHistPhi;//!
  TH1F        *fHistEta;//!
  TH1F        *fHistCharge;//! 
  TH1F        *fHistTPCncls;    //!
  TH1F        *fHistDedx;       //!
  TH1F        *fHistDCAxy;      //!
  TH1F        *fHistDCAz;       //!
  TH1F        *fHistDCAxyAfter; //!
  TH1F        *fHistDCAzAfter;  //!
  TH1F        *fHistPosPt;         //!
  TH1F        *fHistPosPhi;        //!
  TH1F        *fHistPosEta;        //!
  TH1F        *fHistNegPt;         //!
  TH1F        *fHistNegPhi;        //!
  TH1F        *fHistNegEta;        //!
  TH1F        *fHistPtPion;        //!
  TH1F        *fHistPtPionH;       //!
  TH1F        *fHistPtPionL;       //!
  TH1F        *fHistNsigmaPion;    //!
  TH1F        *fHistPtProt;        //!
  TH1F        *fHistPtProtH;       //!
  TH1F        *fHistPtProtL;       //!
  TH1F        *fHistNsigmaProt;    //!
  TH1F        *fHistPtMC;       //!
  TH1F        *fHistPhiMC;      //!
  TH1F        *fHistEtaMC;      //!
  TH1F        *fHistChargeMC;   //! 
  TH1F        *fHistTPCnclsMC;  //!
  TH1F        *fHistDedxMC;     //!
  TH1F        *fHistDCAxyMC;    //!
  TH1F        *fHistDCAzMC;     //!
  TH1F        *fHistPlaneV0h2;  //!
  TH1F        *fHistPlaneV0Ah2; //!
  TH1F        *fHistPlaneV0Ch2; //!
  TH1F        *fHistPlaneV0ACDCh2;//!
  TH1F        *fHistZVtx;       //!
  TH1F        *fHistZVtxD;      //!
  TH1F        *fHistZVtxMC;     //!
  TH1F        *fHistZVtxDiffMC; //!
  TH1F        *fHistCentTRK;    //!
  TH1F        *fHistCentV0M;    //!
  TH1F        *fHistCentDIFF;   //!
  TH1F        *fHistCentDIAG;   //!
  TH1F        *fHistCtrkDIAG;   //!
  TH1F        *fHistVtxRDIAG;   //!
  TH1F        *fHistVtxSDIAG;   //!
  TH1F        *fHistVtxRDIBG;   //!
  TH1F        *fHistVtxSDIBG;   //!
  TH1F        *fHistCentTRKAVEkMB;//!
  TH1F        *fHistCentV0MAVEkMB;//!
  TH1F        *fHistCentTRKAVEkCentral;//!
  TH1F        *fHistCentV0MAVEkCentral;//!
  TH1F        *fHistCentTRKAVEkSemiCentral;//!
  TH1F        *fHistCentV0MAVEkSemiCentral;//!
  TH1F        *fHistCentTRKAVEkA3;//!
  TH1F        *fHistCentV0MAVEkA3;//!
  TH1F        *fHistCentTRKAVEkSel;//!
  TH1F        *fHistCentV0MAVEkSel;//!
  
  TProfile *h_MCq22_cent;//!
  TProfile *h_MCq22_centP;//!
  TProfile *h_MCq22_centN;//!
  TProfile *h_MCAq22_cent;//!
  TProfile *h_MCAq22_centP;//!
  TProfile *h_MCAq22_centN;//!
  TProfile *h_MCq22gap0_cent;//!
  TProfile *h_MCq22gap0_centP;//!
  TProfile *h_MCq22gap0_centN;//!
  TProfile *h_MCq22gap1_cent;//!
  TProfile *h_MCq22gap1_centP;//!
  TProfile *h_MCq22gap1_centN;//!
  
  TProfile *hMC_AT_X_deta;//!
  TProfile *hMC_AT_X_detaP;//!
  TProfile *hMC_AT_X_detaN;//!
  TProfile *hMC_diffq22_X_deta;//!
  TProfile *hMC_diffq22_X_detaP;//!
  TProfile *hMC_diffq22_X_detaN;//!
  TProfile *hMC_ATdiffq22_X_deta;//!
  TProfile *hMC_ATdiffq22_X_detaP;//!
  TProfile *hMC_ATdiffq22_X_detaN;//!
  TProfile *hMC_diffq32_X_deta;//!
  TProfile *hMC_diffq32_X_detaP;//!
  TProfile *hMC_diffq32_X_detaN;//!
  TProfile *hMC_ATdiffq32_X_deta;//!
  TProfile *hMC_ATdiffq32_X_detaP;//!
  TProfile *hMC_ATdiffq32_X_detaN;//!
  TProfile *hMC_diffq42_X_deta;//!
  TProfile *hMC_diffq42_X_detaP;//!
  TProfile *hMC_diffq42_X_detaN;//!
  TProfile *hMC_ATdiffq42_X_deta;//!
  TProfile *hMC_ATdiffq42_X_detaP;//!
  TProfile *hMC_ATdiffq42_X_detaN;//!
  
  TProfile *fProfMeanChargePt;//!
  TProfile *fProfMeanChargeEta;//!
  
  TProfile *h_X22_cent;//!
  TProfile *h_X22_centP;//!
  TProfile *h_X22_centN;//!
  TProfile *h_Y22_cent;//!
  TProfile *h_Y22_centP;//!
  TProfile *h_Y22_centN;//!
  TProfile *h_X22_lo_cent;//!
  TProfile *h_X22_lo_centP;//!
  TProfile *h_X22_lo_centN;//!
  TProfile *h_Y22_lo_cent;//!
  TProfile *h_Y22_lo_centP;//!
  TProfile *h_Y22_lo_centN;//!
  TProfile *h_X22_li_cent;//!
  TProfile *h_X22_li_centP;//!
  TProfile *h_X22_li_centN;//!
  TProfile *h_Y22_li_cent;//!
  TProfile *h_Y22_li_centP;//!
  TProfile *h_Y22_li_centN;//!
  TProfile *h_X22_ri_cent;//!
  TProfile *h_X22_ri_centP;//!
  TProfile *h_X22_ri_centN;//!
  TProfile *h_Y22_ri_cent;//!
  TProfile *h_Y22_ri_centP;//!
  TProfile *h_Y22_ri_centN;//!
  TProfile *h_X22_ro_cent;//!
  TProfile *h_X22_ro_centP;//!
  TProfile *h_X22_ro_centN;//!
  TProfile *h_Y22_ro_cent;//!
  TProfile *h_Y22_ro_centP;//!
  TProfile *h_Y22_ro_centN;//!
  
  TProfile *h_q22_cent;//!
  TProfile *h_q22_centP;//!
  TProfile *h_q22_centN;//!
  TProfile *h_q22_centPP;//!
  TProfile *h_q22_centNN;//!
  TProfile *h_q23_cent;//!
  TProfile *h_q23_centP;//!
  TProfile *h_q23_centN;//!
  TProfile *h_q24_cent;//!
  TProfile *h_q24_centP;//!
  TProfile *h_q24_centN;//!
  TProfile *h_q22gap0_cent;//!
  TProfile *h_q22gap0_centP;//!
  TProfile *h_q22gap0_centN;//!
  TProfile *h_q22gap1_cent;//!
  TProfile *h_q22gap1_centP;//!
  TProfile *h_q22gap1_centN;//!
  
  TProfile *h_q32_cent;//!
  TProfile *h_q32_centP;//!
  TProfile *h_q32_centN;//!
  TProfile *h_q32_centPP;//!
  TProfile *h_q32_centNN;//!
  TProfile *h_q32gap0_cent;//!
  TProfile *h_q32gap0_centP;//!
  TProfile *h_q32gap0_centN;//!
  TProfile *h_q32gap1_cent;//!
  TProfile *h_q32gap1_centP;//!
  TProfile *h_q32gap1_centN;//!
  
  TProfile *h_q42_cent;//!
  TProfile *h_q42_centP;//!
  TProfile *h_q42_centN;//!
  TProfile *h_q42_centPP;//!
  TProfile *h_q42_centNN;//!
  TProfile *h_q42gap0_cent;//!
  TProfile *h_q42gap0_centP;//!
  TProfile *h_q42gap0_centN;//!
  TProfile *h_q42gap1_cent;//!
  TProfile *h_q42gap1_centP;//!
  TProfile *h_q42gap1_centN;//!
  
  
  
  TProfile *h_q22qasym_cent[10];//!
  TProfile *h_q22qasymP_cent[10];//!
  TProfile *h_q22qasymN_cent[10];//!
  TProfile *h_q22qasym_gap0_cent[10];//!
  TProfile *h_q22qasymP_gap0_cent[10];//!
  TProfile *h_q22qasymN_gap0_cent[10];//!
  TProfile *h_q22qasym_gap1_cent[10];//!
  TProfile *h_q22qasymP_gap1_cent[10];//!
  TProfile *h_q22qasymN_gap1_cent[10];//!
  TProfile *h_q24qasym_cent[10];//!
  TProfile *h_q24qasymP_cent[10];//!
  TProfile *h_q24qasymN_cent[10];//!
  TProfile *h_q32qasym_cent[10];//!
  TProfile *h_q32qasymP_cent[10];//!
  TProfile *h_q32qasymN_cent[10];//!
  TProfile *h_q32qasym_gap0_cent[10];//!
  TProfile *h_q32qasymP_gap0_cent[10];//!
  TProfile *h_q32qasymN_gap0_cent[10];//!
  TProfile *h_q32qasym_gap1_cent[10];//!
  TProfile *h_q32qasymP_gap1_cent[10];//!
  TProfile *h_q32qasymN_gap1_cent[10];//!
  TProfile *h_q34qasym_cent[10];//!
  TProfile *h_q34qasymP_cent[10];//!
  TProfile *h_q34qasymN_cent[10];//!
  TProfile *h_q42qasym_cent[10];//!
  TProfile *h_q42qasymP_cent[10];//!
  TProfile *h_q42qasymN_cent[10];//!
  TProfile *h_q42qasym_gap0_cent[10];//!
  TProfile *h_q42qasymP_gap0_cent[10];//!
  TProfile *h_q42qasymN_gap0_cent[10];//!
  TProfile *h_q42qasym_gap1_cent[10];//!
  TProfile *h_q42qasymP_gap1_cent[10];//!
  TProfile *h_q42qasymN_gap1_cent[10];//!
  TProfile *h_q44qasym_cent[10];//!
  TProfile *h_q44qasymP_cent[10];//!
  TProfile *h_q44qasymN_cent[10];//!
  
  TH1F *h_eta_pos_F1;//!
  TH1F *h_eta_pos_F3;//!
  TH1F *h_eta_neg_F1;//!
  TH1F *h_eta_neg_F3;//!
  TH1F *h_cut_eta_pos_F1;//!
  TH1F *h_cut_eta_pos_F3;//!
  TH1F *h_cut_eta_neg_F1;//!
  TH1F *h_cut_eta_neg_F3;//!
  
  TProfile *h_AT_eta;//!
  TProfile *h_AT_etaF1;//!
  TProfile *h_AT_etaF3;//!
  TProfile *h_AT_etaMC;//!
  TH2F *h2_AT_eta;//!
  TH2F *h2_AT_etaF1;//!
  TH2F *h2_AT_etaF3;//!
  TH2F *h2_AT_etaMC;//!

  TProfile *h_AT_cut_eta;//!
  TProfile *h_AT_cut_etaF1;//!
  TProfile *h_AT_cut_etaF3;//!
  TProfile *h_AT_cut_etaMC;//!
  TH2F *h2_AT_cut_eta;//!
  TH2F *h2_AT_cut_etaF1;//!
  TH2F *h2_AT_cut_etaF3;//!
  TH2F *h2_AT_cut_etaMC;//!

  TProfile *h_A_cent;//!
  TProfile *h_rA_cent;//!
  TProfile *h_r1A_cent;//!
  TProfile *h_r2A_cent;//!
  TProfile *h_A_centF1;//!
  TProfile *h_A_centF3;//!
  TProfile *h_A_centMC;//!
  TH2F *h2_A_cent;//!
  TH2F *h2_rA_cent;//!
  TH2F *h2_r1A_cent;//!
  TH2F *h2_r2A_cent;//!
  TH2F *h2_A_centF1;//!
  TH2F *h2_A_centF3;//!
  TH2F *h2_A_centMC;//!


  TProfile *h_Aq22_cent;//!
  TProfile *h_Aq22_centP;//!
  TProfile *h_Aq22_centN;//!
  TProfile *h_Aq22_centPP;//!
  TProfile *h_Aq22_centNN;//!
  TProfile *h_Aq22_SL_cent;//!
  TProfile *h_Aq22_SL_centP;//!
  TProfile *h_Aq22_SL_centN;//!
  TProfile *h_Aq22_SLL_cent;//!
  TProfile *h_Aq22_SLL_centP;//!
  TProfile *h_Aq22_SLL_centN;//!
  TProfile *h_Aq22_SLR_cent;//!
  TProfile *h_Aq22_SLR_centP;//!
  TProfile *h_Aq22_SLR_centN;//!
  TProfile *h_Aq32_cent;//!
  TProfile *h_Aq32_centP;//!
  TProfile *h_Aq32_centN;//!
  TProfile *h_Aq32_centPP;//!
  TProfile *h_Aq32_centNN;//!
  TProfile *h_Aq32_SL_cent;//!
  TProfile *h_Aq32_SL_centP;//!
  TProfile *h_Aq32_SL_centN;//!
  TProfile *h_Aq32_SLL_cent;//!
  TProfile *h_Aq32_SLL_centP;//!
  TProfile *h_Aq32_SLL_centN;//!
  TProfile *h_Aq32_SLR_cent;//!
  TProfile *h_Aq32_SLR_centP;//!
  TProfile *h_Aq32_SLR_centN;//!
  TProfile *h_Aq42_cent;//!
  TProfile *h_Aq42_centP;//!
  TProfile *h_Aq42_centN;//!
  TProfile *h_Aq42_centPP;//!
  TProfile *h_Aq42_centNN;//!
  TProfile *h_Aq42_SL_cent;//!
  TProfile *h_Aq42_SL_centP;//!
  TProfile *h_Aq42_SL_centN;//!
  TProfile *h_Aq42_SLL_cent;//!
  TProfile *h_Aq42_SLL_centP;//!
  TProfile *h_Aq42_SLL_centN;//!
  TProfile *h_Aq42_SLR_cent;//!
  TProfile *h_Aq42_SLR_centP;//!
  TProfile *h_Aq42_SLR_centN;//!

  TProfile *h_diffq22_pt;//!
  TProfile *h_diffq22_ptP;//!
  TProfile *h_diffq22_ptN;//!
  TProfile *h_diffAq22_pt;//!
  TProfile *h_diffAq22_ptP;//!
  TProfile *h_diffAq22_ptN;//!
  TProfile *h_diffq22_eta;//!
  TProfile *h_diffq22_etaP;//!
  TProfile *h_diffq22_etaN;//!
  TProfile *h_diffAq22_eta;//!
  TProfile *h_diffAq22_etaP;//!
  TProfile *h_diffAq22_etaN;//!
  TProfile *h_diffq32_pt;//!
  TProfile *h_diffq32_ptP;//!
  TProfile *h_diffq32_ptN;//!
  TProfile *h_diffAq32_pt;//!
  TProfile *h_diffAq32_ptP;//!
  TProfile *h_diffAq32_ptN;//!
  TProfile *h_diffq32_eta;//!
  TProfile *h_diffq32_etaP;//!
  TProfile *h_diffq32_etaN;//!
  TProfile *h_diffAq32_eta;//!
  TProfile *h_diffAq32_etaP;//!
  TProfile *h_diffAq32_etaN;//!
  TProfile *h_diffq42_pt;//!
  TProfile *h_diffq42_ptP;//!
  TProfile *h_diffq42_ptN;//!
  TProfile *h_diffAq42_pt;//!
  TProfile *h_diffAq42_ptP;//!
  TProfile *h_diffAq42_ptN;//!
  TProfile *h_diffq42_eta;//!
  TProfile *h_diffq42_etaP;//!
  TProfile *h_diffq42_etaN;//!
  TProfile *h_diffAq42_eta;//!
  TProfile *h_diffAq42_etaP;//!
  TProfile *h_diffAq42_etaN;//!


  TProfile *h_AT_X_deta;//!
  TProfile *h_AT_X_detaP;//!
  TProfile *h_AT_X_detaN;//!
  TProfile *h_diffq22_X_deta;//!
  TProfile *h_diffq22_X_detaP;//!
  TProfile *h_diffq22_X_detaN;//!
  TProfile *h_ATdiffq22_X_deta;//!
  TProfile *h_ATdiffq22_X_detaP;//!
  TProfile *h_ATdiffq22_X_detaN;//!
  TProfile *h_diffq32_X_deta;//!
  TProfile *h_diffq32_X_detaP;//!
  TProfile *h_diffq32_X_detaN;//!
  TProfile *h_ATdiffq32_X_deta;//!
  TProfile *h_ATdiffq32_X_detaP;//!
  TProfile *h_ATdiffq32_X_detaN;//!
  TProfile *h_diffq42_X_deta;//!
  TProfile *h_diffq42_X_detaP;//!
  TProfile *h_diffq42_X_detaN;//!
  TProfile *h_ATdiffq42_X_deta;//!
  TProfile *h_ATdiffq42_X_detaP;//!
  TProfile *h_ATdiffq42_X_detaN;//!

  TProfile *h_AT_X_deta_F1;//!
  TProfile *h_AT_X_detaP_F1;//!
  TProfile *h_AT_X_detaN_F1;//!
  TProfile *h_diffq22_X_deta_F1;//!
  TProfile *h_diffq22_X_detaP_F1;//!
  TProfile *h_diffq22_X_detaN_F1;//!
  TProfile *h_ATdiffq22_X_deta_F1;//!
  TProfile *h_ATdiffq22_X_detaP_F1;//!
  TProfile *h_ATdiffq22_X_detaN_F1;//!
  TProfile *h_diffq32_X_deta_F1;//!
  TProfile *h_diffq32_X_detaP_F1;//!
  TProfile *h_diffq32_X_detaN_F1;//!
  TProfile *h_ATdiffq32_X_deta_F1;//!
  TProfile *h_ATdiffq32_X_detaP_F1;//!
  TProfile *h_ATdiffq32_X_detaN_F1;//!
  TProfile *h_diffq42_X_deta_F1;//!
  TProfile *h_diffq42_X_detaP_F1;//!
  TProfile *h_diffq42_X_detaN_F1;//!
  TProfile *h_ATdiffq42_X_deta_F1;//!
  TProfile *h_ATdiffq42_X_detaP_F1;//!
  TProfile *h_ATdiffq42_X_detaN_F1;//!

  TProfile *h_AT_X_deta_F3;//!
  TProfile *h_AT_X_detaP_F3;//!
  TProfile *h_AT_X_detaN_F3;//!
  TProfile *h_diffq22_X_deta_F3;//!
  TProfile *h_diffq22_X_detaP_F3;//!
  TProfile *h_diffq22_X_detaN_F3;//!
  TProfile *h_ATdiffq22_X_deta_F3;//!
  TProfile *h_ATdiffq22_X_detaP_F3;//!
  TProfile *h_ATdiffq22_X_detaN_F3;//!
  TProfile *h_diffq32_X_deta_F3;//!
  TProfile *h_diffq32_X_detaP_F3;//!
  TProfile *h_diffq32_X_detaN_F3;//!
  TProfile *h_ATdiffq32_X_deta_F3;//!
  TProfile *h_ATdiffq32_X_detaP_F3;//!
  TProfile *h_ATdiffq32_X_detaN_F3;//!
  TProfile *h_diffq42_X_deta_F3;//!
  TProfile *h_diffq42_X_detaP_F3;//!
  TProfile *h_diffq42_X_detaN_F3;//!
  TProfile *h_ATdiffq42_X_deta_F3;//!
  TProfile *h_ATdiffq42_X_detaP_F3;//!
  TProfile *h_ATdiffq42_X_detaN_F3;//!

  TProfile *h_ATXdiffq22_X_detaP;//!
  TProfile *h_ATXdiffq22_X_detaN;//!
  TProfile *h_ATYdiffq22_X_detaP;//!
  TProfile *h_ATYdiffq22_X_detaN;//!
  TProfile *h_ATZdiffq22_X_detaP;//!
  TProfile *h_ATZdiffq22_X_detaN;//!

  TProfile *h_ATXdiffq32_X_detaP;//!
  TProfile *h_ATXdiffq32_X_detaN;//!
  TProfile *h_ATYdiffq32_X_detaP;//!
  TProfile *h_ATYdiffq32_X_detaN;//!
  TProfile *h_ATZdiffq32_X_detaP;//!
  TProfile *h_ATZdiffq32_X_detaN;//!

  TProfile *h_ATXdiffq42_X_detaP;//!
  TProfile *h_ATXdiffq42_X_detaN;//!
  TProfile *h_ATYdiffq42_X_detaP;//!
  TProfile *h_ATYdiffq42_X_detaN;//!
  TProfile *h_ATZdiffq42_X_detaP;//!
  TProfile *h_ATZdiffq42_X_detaN;//!

  TProfile *h_AT_S_deta;//!
  TProfile *h_AT_S_detaP;//!
  TProfile *h_AT_S_detaN;//!
  TProfile *h_diffq22_S_deta;//!
  TProfile *h_diffq22_S_detaP;//!
  TProfile *h_diffq22_S_detaN;//!
  TProfile *h_ATdiffq22_S_deta;//!
  TProfile *h_ATdiffq22_S_detaP;//!
  TProfile *h_ATdiffq22_S_detaN;//!
  TProfile *h_diffq32_S_deta;//!
  TProfile *h_diffq32_S_detaP;//!
  TProfile *h_diffq32_S_detaN;//!
  TProfile *h_ATdiffq32_S_deta;//!
  TProfile *h_ATdiffq32_S_detaP;//!
  TProfile *h_ATdiffq32_S_detaN;//!
  TProfile *h_diffq42_S_deta;//!
  TProfile *h_diffq42_S_detaP;//!
  TProfile *h_diffq42_S_detaN;//!
  TProfile *h_ATdiffq42_S_deta;//!
  TProfile *h_ATdiffq42_S_detaP;//!
  TProfile *h_ATdiffq42_S_detaN;//!


  TProfile *h_AT_X_dpt;//!
  TProfile *h_AT_X_dptP;//!
  TProfile *h_AT_X_dptN;//!
  TProfile *h_diffq22_X_dpt;//!
  TProfile *h_diffq22_X_dptP;//!
  TProfile *h_diffq22_X_dptN;//!
  TProfile *h_ATdiffq22_X_dpt;//!
  TProfile *h_ATdiffq22_X_dptP;//!
  TProfile *h_ATdiffq22_X_dptN;//!
  TProfile *h_diffq32_X_dpt;//!
  TProfile *h_diffq32_X_dptP;//!
  TProfile *h_diffq32_X_dptN;//!
  TProfile *h_ATdiffq32_X_dpt;//!
  TProfile *h_ATdiffq32_X_dptP;//!
  TProfile *h_ATdiffq32_X_dptN;//!
  TProfile *h_diffq42_X_dpt;//!
  TProfile *h_diffq42_X_dptP;//!
  TProfile *h_diffq42_X_dptN;//!
  TProfile *h_ATdiffq42_X_dpt;//!
  TProfile *h_ATdiffq42_X_dptP;//!
  TProfile *h_ATdiffq42_X_dptN;//!



  TProfile *h_a1q22_cent;//!
  TProfile *h_a1q22_centP;//!
  TProfile *h_a1q22_centN;//!

  TProfile *h_a2q22_cent;//!
  TProfile *h_a2q22_centP;//!
  TProfile *h_a2q22_centN;//!

  TProfile *h_rAq22_X1_cent;//!
  TProfile *h_rAq22_X1_centP;//!
  TProfile *h_rAq22_X1_centN;//!
  TProfile *h_rAq22_X2_cent;//!
  TProfile *h_rAq22_X2_centP;//!
  TProfile *h_rAq22_X2_centN;//!
  TProfile *h_rAq22_X3_cent;//!
  TProfile *h_rAq22_X3_centP;//!
  TProfile *h_rAq22_X3_centN;//!
  TProfile *h_rAq22_X4_cent;//!
  TProfile *h_rAq22_X4_centP;//!
  TProfile *h_rAq22_X4_centN;//!


   
  AliAnalysisTaskCMEv2A(const AliAnalysisTaskCMEv2A&); // copy constructor not implemented
  AliAnalysisTaskCMEv2A& operator=(const AliAnalysisTaskCMEv2A&); // assignment operator not implemented
  
  ClassDef(AliAnalysisTaskCMEv2A, 1); // Add this class as ROOT class (inherit from TObject)
};

#endif
