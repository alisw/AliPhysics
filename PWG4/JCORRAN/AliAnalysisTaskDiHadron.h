#ifndef AliAnalysisTaskDiHadron_cxx
#define AliAnalysisTaskDiHadron_cxx


class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "AliAnalysisTask.h"


class AliAnalysisTaskDiHadron : public AliAnalysisTask{
  //Default constructor
 public:


  AliAnalysisTaskDiHadron(const char *name="AliAnalysisTaskDiHadron");

  virtual ~AliAnalysisTaskDiHadron() {}
  //virtual ~AliAnalysisTaskDiHadron();

 

  //AliAnalysisTask Functions
 
  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *);
 

  void SetEfficiencies(Float_t EffFitPt, TF1 *FitLow, TF1 *FitHigh, Int_t NFitLowParam, Int_t NFitHighParam, Float_t *FitLowParam, Float_t *FitHighParam);
  void SetBins(Int_t nBinPhi, Int_t nBinEta, Int_t nBinPhiEtaPhi, Int_t nBinPhiEtaEta, Int_t nBinPhi3, Int_t nBinEta3, Float_t dPhiMin, Float_t dPhiMax, Int_t NTPtBins, Int_t NMixBins, Int_t NCentBins,Int_t NAPtBins, Int_t NAPt3Bins, Int_t NVertexBins, Int_t NXEBin,Float_t *PtTrigArray, Float_t *PtAssocArray, Float_t *PtAssoc3Array1, Float_t *PtAssoc3Array2, Int_t *CentArrayMin, Int_t *CentArrayMax, Float_t *XEArray);
  void SetOptions(Int_t fEfficiencyCorr, Int_t DEBUG,Int_t fMCHistos);
  void SetCuts(Int_t MinClutersTPC, Float_t MinClusterRatio, Float_t MaxTPCchi2, Int_t MinClustersITS, Float_t EtaCut, Float_t TrigEtaCut, Float_t NearPhiCut, Float_t XECut, Float_t MaxDCA, Float_t MaxDCAXY, Float_t MaxDCAZ, Int_t DCA2D, Int_t TPCRefit, Int_t ITSRefit, Int_t SPDCut, Float_t MinPtAssoc, Float_t MaxPtAssoc, Float_t VzCut, Int_t NIDs, char **TrigIDArray);
  Int_t CheckVertex(AliESDEvent *rESD);
  Int_t CheckTrigger(AliESDEvent *rESD);
  Int_t TrackCuts(AliESDEvent *rESD, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks);
  Int_t TrackCutsMC(AliMCEvent *rMC, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks);

  
 private:

//Maximum numbers for creating the arrays
  enum{NUMBER_OF_PT_BINS=20, 
       NUMBER_OF_CENT_BINS=10, 
       NUMBER_OF_APT_BINS=50,
       NUMBER_OF_APT3_BINS=50,
       NUMBER_OF_TRIGGER_IDS=1,
       NUMBER_OF_VERTEX_BINS=20,
       NUMBER_OF_XE_BINS=20,
       NUMBER_OF_EVENTS_TO_MIX=100};
  // virtual void SetParameters();
  AliESDEvent *fESD;
  AliMCEvent *fMC;
  TList *fOutput;
  
   //Cuts
  
    Int_t fMinClustersTPC;
    Float_t fMinClusterRatio;
    Float_t fMaxTPCchi2;
    Int_t fMinClustersITS;
    Float_t fEtaCut;
    Float_t fTrigEtaCut;
    Float_t fNearPhiCut;
    Float_t fXECut;
    Float_t fMaxDCA;
    Float_t fMaxDCAXY;
    Float_t fMaxDCAZ;
    Int_t fDCA2D;
    Int_t fTPCRefit;
    Int_t fITSRefit;
    Int_t fSPDCut;
    Float_t fMinPtAssoc;
    Float_t fMaxPtAssoc;
    Float_t fVzCut;
    Int_t fEfficiencyCorr;//Toggle correcting of efficiencies when filling histograms
    Int_t DEBUG;
    //Binning
    Int_t fnBinPhi;
    Int_t fnBinEta;
    Int_t fnBinPhiEtaPhi;
    Int_t fnBinPhiEtaEta;
    Int_t fnBinPhi3;
    Int_t fnBinEta3;
    Float_t fPi;
    Float_t fdPhiMin;
    Float_t fdPhiMax;
    //Parameters for settings
    Int_t fNTPtBins;
    Int_t fNMix;
    Int_t fNCentBins;
    Int_t fNAPtBins;
    Int_t fNAPt3Bins;
    Int_t fNVertexBins;
    Int_t fNXEBins;
    Int_t fNIDs;
    Float_t fEffFitPt;
    Int_t fNFitLowParam;
    Int_t fNFitHighParam;
    Int_t fMCHistos;
    //TF1 *fFitLow;
    //  TF1 *fFitHigh;
    Int_t fNFitLow;
    Int_t fNFitHigh;

    TF1 *fFitLow;
    TF1 *fFitHigh;
    Float_t *fFitLowParam; 
    Float_t *fFitHighParam;
    Float_t *fPtTrigArray;
    Float_t *fPtAssocArray;
    Float_t *fPtAssoc3Array1;
    Float_t *fPtAssoc3Array2; 
    Int_t *fCentArrayMin;
    Int_t *fCentArrayMax;
    Float_t *fXEArray;
    char **fTrigIDArray;

    /*
    //Float_t *FitLowParam;
    //   Float_t *FitHighParam;
    //    Float_t fPtTrigArray[(NUMBER_OF_PT_BINS+1)];
    Float_t fPtAssocArray[(NUMBER_OF_APT_BINS+1)];
    Float_t fPtAssoc3Array1[NUMBER_OF_APT3_BINS];
    Float_t fPtAssoc3Array2[NUMBER_OF_APT3_BINS];
    Int_t fCentArrayMin[NUMBER_OF_CENT_BINS];
    Int_t fCentArrayMax[NUMBER_OF_CENT_BINS];
    Float_t fXEArray[(NUMBER_OF_XE_BINS+1)];
    char *fTrigIDArray[NUMBER_OF_TRIGGER_IDS];
  */
 
  Float_t fVertexArray[(NUMBER_OF_VERTEX_BINS+1)];
  
  //Histograms
  TH1F *fHistPt[NUMBER_OF_CENT_BINS][2];
  TH1F *fHistPtEff[NUMBER_OF_CENT_BINS][2];
  TH1F *fHistPtTrig[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH1F *fHistMult[2];
  TH1F *fHistMultTrig[NUMBER_OF_PT_BINS][2];

  TH2F *fHistPhi[NUMBER_OF_CENT_BINS][2];
  TH2F *fHistPhiTrig[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH2F *fHistDeltaPhi[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];
  TH2F *fHistDeltaPhiMix[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];

  TH2F *fHistPhiPt[NUMBER_OF_CENT_BINS][2];
  TH2F *fHistPhiTrigPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH2F *fHistDeltaPhiPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];
  TH2F *fHistDeltaPhiMixPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];

  TH2F *fHistEta[NUMBER_OF_CENT_BINS][2];
  TH2F *fHistEtaTrig[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];

  TH2F *fHistEtaPt[NUMBER_OF_CENT_BINS][2];
  TH2F *fHistEtaTrigPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];

  TH2F *fHistDeltaEtaN[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];
  TH2F *fHistDeltaEtaNMix[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];

  TH2F *fHistDeltaEtaNPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];
  TH2F *fHistDeltaEtaNMixPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];

  TH2F *fHistDeltaEtaA[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];
  TH2F *fHistDeltaEtaAMix[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];

  TH2F *fHistDeltaEtaAPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];
  TH2F *fHistDeltaEtaAMixPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][3][2];

  TH1F *fHistNEvents[NUMBER_OF_CENT_BINS][2];
  TH1F *fHistNTrigger[NUMBER_OF_CENT_BINS][2];
  TH1F *fHistNTriggerPt[NUMBER_OF_CENT_BINS][2];
  TH1F *fHistNMix[NUMBER_OF_CENT_BINS][2];

  TH3F *fHistPhiEta[NUMBER_OF_CENT_BINS][2];
  TH3F *fHistPhiEtaTrig[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH3F *fHistDeltaPhiEta[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH3F *fHistDeltaPhiEtaMix[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];

  TH3F *fHistPhiEtaPt[NUMBER_OF_CENT_BINS][2];
  TH3F *fHistPhiEtaTrigPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH3F *fHistDeltaPhiEtaPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH3F *fHistDeltaPhiEtaMixPt[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];

  TH1F *fHistXEN[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH1F *fHistXENMix[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH1F *fHistXEA[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];
  TH1F *fHistXEAMix[NUMBER_OF_PT_BINS][NUMBER_OF_CENT_BINS][2];

  //Three-Particle Histograms
  TH2F *fHistDeltaPhiPhi[NUMBER_OF_PT_BINS][NUMBER_OF_APT3_BINS][NUMBER_OF_CENT_BINS][4][2];
  TH2F *fHistDeltaPhiPhiMix[NUMBER_OF_PT_BINS][NUMBER_OF_APT3_BINS][NUMBER_OF_CENT_BINS][4][2];
  TH2F *fHistDeltaPhiPhiSS[NUMBER_OF_PT_BINS][NUMBER_OF_APT3_BINS][NUMBER_OF_CENT_BINS][4][2];
  TH2F *fHistDeltaEtaEta[NUMBER_OF_PT_BINS][NUMBER_OF_APT3_BINS][NUMBER_OF_CENT_BINS][4][2];
  TH2F *fHistDeltaEtaEtaMix[NUMBER_OF_PT_BINS][NUMBER_OF_APT3_BINS][NUMBER_OF_CENT_BINS][4][2];
  TH2F *fHistDeltaEtaEtaSS[NUMBER_OF_PT_BINS][NUMBER_OF_APT3_BINS][NUMBER_OF_CENT_BINS][4][2];


 

  //Arrays for event mixing
  Float_t *fMPt[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
  Float_t *fMPhi[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2]; 
  Int_t fMixTrack[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
  Float_t *fMEta[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
  Short_t *fMCharge[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
   Float_t *fMEff[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
  Int_t fMixPointer[NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
  Int_t fMixEnd[NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2]; 
  //additional ones to speed up the 3-particle correlations
  Short_t *fMPtAssoc3[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2][10];//using a bit extra memory fixing the pt bins to 10, but makes the arrays much easier
  Short_t *fMNPtAssoc3[NUMBER_OF_EVENTS_TO_MIX][NUMBER_OF_CENT_BINS][NUMBER_OF_VERTEX_BINS][2];
  

  //Arrays for Main Events
  Float_t *tPhi;
  Float_t *tEta;
  Float_t *tPt;
  Short_t *tCharge;
  Float_t *tEff;
  Int_t **tPtAssoc3;
  Int_t *tNPtAssoc3;

  AliAnalysisTaskDiHadron(const AliAnalysisTaskDiHadron&);//not implimented
  AliAnalysisTaskDiHadron& operator=(const AliAnalysisTaskDiHadron&);//not implimnted

  ClassDef(AliAnalysisTaskDiHadron,1);
  
};
#endif
  
  
