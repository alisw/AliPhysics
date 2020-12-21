#ifndef AliAnalysisTaskPHOSTrigPi0_cxx
#define AliAnalysisTaskPHOSTrigPi0_cxx

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class THnSparse;
class TGraphAsymmErrors;
class TRandom3;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliESDEvent ;
class AliESDHeader ;
class AliAODInputHandler;
class AliPHOSCalibData;
class AliESDtrack ;
class AliESDCaloCluster ;
class AliStack;
//class AliCaloPhoton_syano;
class AliVCaloTrigger;
class AliAODCaloTrigger;
class AliVCluster;
class AliPHOSEmcCalibData;
class AliCDBStorage;
class AliPIDResponse;

#include "TArrayD.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPHOSTrigPi0 : public AliAnalysisTaskSE {
 
 public:
  
  AliAnalysisTaskPHOSTrigPi0(const char *name = "AliAnalysisTaskPHOSTrigPi0");
  virtual ~AliAnalysisTaskPHOSTrigPi0();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
  void SetMC(Bool_t isMC){fIsMC = isMC;}
  void SetModuleCondition(Int_t mod, Bool_t isGood)         {fIsGoodMod[mod]=isGood;}
  void SetTRUCondition(Int_t mod, Int_t tru, Bool_t isGood) {fIsGoodTRU[mod][tru]=isGood;}
  void SetAnalysisTriggerClass(TString trigName)            {fIsFiredTrig = trigName;}
  void SetClusterCutParameters(Double_t time, Int_t cell, Double_t ene, Double_t disp) {fTOFcut=time; fMinNCell=cell; fMinEene=ene; fDispCut=disp;}

  void SetTOFCalibration(Int_t mod, TH2F* TOFmapLG,  TH2F* TOFmapHG){
    fTimingCalibMapLG[mod] = TOFmapLG;
    fTimingCalibMapHG[mod] = TOFmapHG;
  }
  void SetBadCell(Int_t mod, TH2I* hmap){
    fBadMap[mod] = (TH2I*)hmap->Clone();
    
  }
  void SetBadCell(Int_t mod, TH2F* hmap){
    fBadMap[mod] = (TH2I*)hmap->Clone();
  }
  void SetBadTile(Int_t mod, TH2I* hmap){
    fBadTile[mod] = (TH2I*)hmap->Clone();
  }
  void SetSPDTrackletMultiCorrections(TH1F* unit, TH1F* gap1, TH1F* gap2, Double_t fit_unit, Double_t fit_gap1, Double_t fit_gap2){
    fSPDMultiCorrUnit = (TH1F*)unit->Clone("fSPDMultiCorrUnit");
    fSPDMultiCorrGap1 = (TH1F*)gap1->Clone("fSPDMultiCorrGap1");
    fSPDMultiCorrGap2 = (TH1F*)gap2->Clone("fSPDMultiCorrGap2");
    fRefSPDMultiUnit = fit_unit;
    fRefSPDMultiGap1 = fit_gap1;
    fRefSPDMultiGap2 = fit_gap2;
  }
  void SetV0MultiCorrections(TH1F* fV0AC,Double_t fit_v0ac){
    fV0ACMultiCorr = fV0AC;
    fRefV0ACMulti = fit_v0ac;
  }
  void SetMeanMultiForMultiBin(Double_t unitrap, Double_t gaprap1, Double_t gaprap2, Double_t v0ac){
     fMeanSPDUnitRap = unitrap;
     fMeanSPDGapRap1 = gaprap1;
     fMeanSPDGapRap2 = gaprap2;
     fMeanV0AC       = v0ac;
  }
  void SetMCTriggerThreshold(Int_t Threshold){
    fTriggerThreshold = Threshold;
  }

  
 protected:
  AliAnalysisTaskPHOSTrigPi0(const AliAnalysisTaskPHOSTrigPi0&); // not implemented
  AliAnalysisTaskPHOSTrigPi0& operator=(const AliAnalysisTaskPHOSTrigPi0&); // not implemented
  
  void   EventFlagInit();
  void   SetPeriod();
  void   SetGeometry();
  void   SetCalibration();
  void   GetMCStack();
  void   ConnectInputData();
  void   SetVertex();
  void   DirectPhotonAnalysis(AliAODEvent* AOD);
  void   NeutralPion();
  void   MaxEnergyCellPos(AliAODCaloCells *cells, AliAODCluster* clu, Int_t& maxId);
  void   CalculatedClusterParameters(AliAODCaloCluster* clu, AliAODCaloCells* cells, Double_t calib,
				     Double_t& CalcM20, Double_t& CalcM02, Double_t& time1, Double_t& time2, Double_t& time3,Double_t& time4, Double_t& time5, Double_t& time6);
  void   SetClusterFlag(Int_t PDGcode,Int_t clustCharge,Bool_t &isPhoton,Bool_t &isElectron,Bool_t &isProton,Bool_t &isAntiProton,
			Bool_t &isNeutron,Bool_t &isAntiNeutron,Bool_t &isConversion,Bool_t &isChargedParticle,Bool_t &isNeutralParticle, Bool_t &isChargedPion);
  void   SetMultiplicity();
  void   SetTriggerInfo();
  void   ProcessMC();
  void   UpdataPHOSClusterPool();

  Int_t  GetTRUNum(Int_t cellX, Int_t cellZ);
  Int_t  FindPrimary(AliAODCluster*clu,  Bool_t& sure);
  Int_t  GetPDGCode(Int_t iLabel, Int_t &clustCharge, Int_t &momLabel, Bool_t &isRadius1cmCut);

  Double_t AddEnergyCalib(Int_t mod,Double_t ene);

  Bool_t SetPhysicsSelection(AliAODEvent* AOD);  
  Bool_t IsGoodEvent(AliAODEvent* AOD);
  Bool_t IsGoodCluster(Double_t ene, Int_t ncell, Int_t mod, Int_t ix, Int_t iz);
  Bool_t IsTrackCluster(Int_t nTrack, Double_t dx, Double_t dz, Double_t ene);
  Bool_t IsFireTile(Int_t *trig_relid, Int_t *cluster_relid);
  Bool_t IsTrigger(AliAODCaloTrigger *trigData, AliAODCaloCells *cells, Int_t maxAbsId, Double_t cluE, Double_t tof);
  Bool_t IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz);
  Bool_t IsGoodTile(const char * det, Int_t mod, Int_t ix, Int_t iz);
  
  TString GetTrackPID(AliAODTrack *track);

 protected:
  TList        *fOutputContainer[3];
  TList        *fPHOSClusterList[10][4];
  
  AliAODEvent         *fAOD[10];
  AliPHOSGeometry     *fPHOSGeo;
  AliPHOSEmcCalibData *fCalibDataEmc;
  AliPHOSCalibData    *fPHOSCalibData;
  AliAODCaloTrigger   *fTrigData;
  
  AliCDBStorage       *fCDBstorage;
  
  AliPIDResponse      *fPIDResponse;

  TRandom3            *fRandom;
  TVector3             fVertexVector;
  TClonesArray        *AODMCTrackArray;

  TF1                 *fEneSmearMeanMC;
  TF1                 *fEneSmearSigmaMC;
  TF1                 *fHistTrackCutX;
  TF1                 *fHistTrackCutZ;

  TH1F                *fSPDMultiCorrUnit;
  TH1F                *fSPDMultiCorrGap1;
  TH1F                *fSPDMultiCorrGap2;
  TH1F                *fV0ACMultiCorr;

  TH2F                *fTimingCalibMapLG[3];
  TH2F                *fTimingCalibMapHG[3];

  Int_t fRunNumber;
  Int_t fMinNCell;
  Int_t fMultiBinSPDUnitRap;
  Int_t fMultiBinSPDGapRap1;
  Int_t fMultiBinSPDGapRap2;
  Int_t fMultiBinV0AC;
  Int_t fZvtxBin;
  Int_t fTriggerThreshold;

  
  Double_t fTOFcut;
  Double_t fDispCut;
  Double_t fMinEene;
  Double_t fRefSPDMultiUnit;
  Double_t fRefSPDMultiGap1;
  Double_t fRefSPDMultiGap2;
  Double_t fRefV0ACMulti;
  Double_t fMCTracksUnit;
  Double_t fMCTracksGap1;
  Double_t fMCTracksGap2;
  Double_t fMCTracksV0AC;
  Double_t fMeanSPDUnitRap;
  Double_t fMeanSPDGapRap1;
  Double_t fMeanSPDGapRap2;
  Double_t fMeanV0AC;

  TString fPeriod;
  TString fIsFiredTrig;

  Bool_t fIsMC;
  Bool_t fIsPileUp;
  Bool_t fIsVtxOut10cm;
  Bool_t fIsVtxNoCont;
  Bool_t f0PH0Event;
  Bool_t f0PH0Event_HighEnergy;
  Bool_t fTOFcut0PH0Event;
  Bool_t fTOFcut0PH0Event_HighEnergy;
  Bool_t fIsGoodMod[3];
  Bool_t fIsGoodTRU[3][8];
  
  TH2I* fBadMap[3];
  TH2I* fBadTile[3];

  TH1F* fHistClustEne[4][4][4];
  TH1F* fHistClustEneMIP[4][4][4];
  TH2F* fHistClustEp[4][4][4];
  TH2F* fHistClustM02[4][4][4];
  TH2F* fHistClustM02M20[4][4][4];
  TH2F* fHistClustM20[4][4][4];
  TH2F* fHistClustDx[4][4][4];
  TH2F* fHistClustDz[4][4][4];
  TH2F* fHistClustDr[4][4][4];
  TH2F* fHistClustNcells[4][4][4];

  TH1F* fHistTrueClustEne[10][4][4][4];
  TH1F* fHistTrueClustEneMIP[10][4][4][4];
  TH2F* fHistTrueClustEp[10][4][4][4];
  TH2F* fHistTrueClustM02[10][4][4][4];
  TH2F* fHistTrueClustM02M20[10][4][4][4];
  TH2F* fHistTrueClustM20[10][4][4][4];
  TH2F* fHistTrueClustDx[10][4][4][4];
  TH2F* fHistTrueClustDz[10][4][4][4];
  TH2F* fHistTrueClustDr[10][4][4][4];
  TH2F* fHistTrueClustNcells[10][4][4][4];

  TH1F* fHist0PH0ClustEne[4][4][4];
  TH1F* fHist0PH0ClustEneMIP[4][4][4];
  TH2F* fHist0PH0ClustEp[4][4][4];
  TH2F* fHist0PH0ClustM02[4][4][4];
  TH2F* fHist0PH0ClustM02M20[4][4][4];
  TH2F* fHist0PH0ClustM20[4][4][4];
  TH2F* fHist0PH0ClustDx[4][4][4];
  TH2F* fHist0PH0ClustDz[4][4][4];
  TH2F* fHist0PH0ClustDr[4][4][4];
  TH2F* fHist0PH0ClustNcells[4][4][4];

  TH1F* fHistTrue0PH0ClustEne[10][4][4][4];
  TH1F* fHistTrue0PH0ClustEneMIP[10][4][4][4];
  TH2F* fHistTrue0PH0ClustEp[10][4][4][4];
  TH2F* fHistTrue0PH0ClustM02[10][4][4][4];
  TH2F* fHistTrue0PH0ClustM02M20[10][4][4][4];
  TH2F* fHistTrue0PH0ClustM20[10][4][4][4];
  TH2F* fHistTrue0PH0ClustDx[10][4][4][4];
  TH2F* fHistTrue0PH0ClustDz[10][4][4][4];
  TH2F* fHistTrue0PH0ClustDr[10][4][4][4];
  TH2F* fHistTrue0PH0ClustNcells[10][4][4][4];
  
  TH2F* fHistMassTwoGammas[2][2][3][4][4];
  TH2F* fHistMixMassTwoGammas[2][2][3][4][4];
  TH2F* fHistTrueMassTwoGammas[2][2][4][5][4];

  TH2F* fHistMassTwoGammasMultiUnitRap[2][2][3][10];
  TH2F* fHistMassTwoGammasMultiGapRap1[2][2][3][10];
  TH2F* fHistMassTwoGammasMultiGapRap2[2][2][3][10];
  TH2F* fHistMassTwoGammasMultiV0AC[2][2][3][10];

  TH2F* fHistMixMassTwoGammasMultiUnitRap[2][2][3][10];
  TH2F* fHistMixMassTwoGammasMultiGapRap1[2][2][3][10];
  TH2F* fHistMixMassTwoGammasMultiGapRap2[2][2][3][10];
  TH2F* fHistMixMassTwoGammasMultiV0AC[2][2][3][10];
  
  TH2F* fHistClustTOF[4];
  TH2F* fHistClustTOFv1[4];
  TH2F* fHistClustTOFv2[4];
  TH2F* fHistClustTOFv3[4];
  TH2F* fHistClustTOFv4[4];
  TH2F* fHistClustTOFv5[4];
  TH2F* fHistClustTOFv6[4];

  TH1F* fHistClustEneAll;
  TH1F* fHistClustEneAllTOFcut1;
  TH1F* fHistClustEneAllTOFcut2;
  TH1F* fHistClustEneAllTOFcut3;
  TH1F* fHistClustEneAllTOFcut4;
  TH1F* fHistClustEneNcellCut;
  TH1F* fHistClustEneNcellCutTOFcut1;
  TH1F* fHistClustEneNcellCutTOFcut2;
  TH1F* fHistClustEneNcellCutTOFcut3;
  TH1F* fHistClustEneNcellCutTOFcut4;
  TH1F* fHistGoodTRUClustEneAll;
  TH1F* fHistGoodTRUClustEneAllTOFcut1;
  TH1F* fHistGoodTRUClustEneAllTOFcut2;
  TH1F* fHistGoodTRUClustEneAllTOFcut3;
  TH1F* fHistGoodTRUClustEneAllTOFcut4;
  TH1F* fHistGoodTRUClustEneNcellCut;
  TH1F* fHistGoodTRUClustEneNcellCutTOFcut1;
  TH1F* fHistGoodTRUClustEneNcellCutTOFcut2;
  TH1F* fHistGoodTRUClustEneNcellCutTOFcut3;
  TH1F* fHistGoodTRUClustEneNcellCutTOFcut4;
  TH1F* fHist0PH0ClustEneAll;
  TH1F* fHist0PH0ClustEneAllTOFcut1;
  TH1F* fHist0PH0ClustEneAllTOFcut2;
  TH1F* fHist0PH0ClustEneAllTOFcut3;
  TH1F* fHist0PH0ClustEneAllTOFcut4;
  TH1F* fHist0PH0ClustEneNcellCut;
  TH1F* fHist0PH0ClustEneNcellCutTOFcut1;
  TH1F* fHist0PH0ClustEneNcellCutTOFcut2;
  TH1F* fHist0PH0ClustEneNcellCutTOFcut3;
  TH1F* fHist0PH0ClustEneNcellCutTOFcut4;
  TH1F* fHistGoodTRU0PH0ClustEneAll;
  TH1F* fHistGoodTRU0PH0ClustEneAllTOFcut1;
  TH1F* fHistGoodTRU0PH0ClustEneAllTOFcut2;
  TH1F* fHistGoodTRU0PH0ClustEneAllTOFcut3;
  TH1F* fHistGoodTRU0PH0ClustEneAllTOFcut4;
  TH1F* fHistGoodTRU0PH0ClustEneNcellCut;
  TH1F* fHistGoodTRU0PH0ClustEneNcellCutTOFcut1;
  TH1F* fHistGoodTRU0PH0ClustEneNcellCutTOFcut2;
  TH1F* fHistGoodTRU0PH0ClustEneNcellCutTOFcut3;
  TH1F* fHistGoodTRU0PH0ClustEneNcellCutTOFcut4;

  TH2F* fHistCellTimeHG;
  TH2F* fHistCellTimeLG;
  TH2F* fHistCellTimeHG_HE;
  TH2F* fHistCellTimeLG_HE;

  TH2F* fHist0PH0CellTimeHG;
  TH2F* fHist0PH0CellTimeLG;
  TH2F* fHist0PH0CellTimeHG_HE;
  TH2F* fHist0PH0CellTimeLG_HE;
  
  TH2F* fHistClustOccupMap[3];
  TH2F* fHist0PH0ClustOccupMap[3];

  TH1F* fHistTotalEvents;
  TH1F* fHistAnalyzedEvents;
  TH1F* fHistPUEvents;
  TH1F* fHistPhysSelectionEvents;
  TH1F* fHistAnalyzed0PH0Events;

  TH2F* fHistSPDTrackletsUnitRap;
  TH2F* fHistSPDTrackletsGapRap1;
  TH2F* fHistSPDTrackletsGapRap2;
  TH2F* fHistCorrectedSPDTrackletsUnitRap;
  TH2F* fHistCorrectedSPDTrackletsGapRap1;
  TH2F* fHistCorrectedSPDTrackletsGapRap2;

  TH2F* fHistCorrelationMCTrackSPDTrackletsUnitRap;
  TH2F* fHistCorrelationMCTrackSPDTrackletsGapRap1;
  TH2F* fHistCorrelationMCTrackSPDTrackletsGapRap2;
  TH2F* fHistCorrelationMCTrackCorrectedSPDTrackletsUnitRap;
  TH2F* fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap1;
  TH2F* fHistCorrelationMCTrackCorrectedSPDTrackletsGapRap2;
  TH2F* fHistCorrelationMCTrackV0AC;
  TH2F* fHistCorrelationMCTrackCorrectedV0AC;
  TH2F* fHistCorrelationSPDTrackletsV0AC;
  TH2F* fHistCorrelationCorrectedSPDTrackletsCorrectedV0AC;
    
  TH2F* fHistV0AMulti;
  TH2F* fHistV0CMulti;
  TH2F* fHistV0ACMulti;
  TH2F* fHistCorrectedV0ACMulti;

  TH2F* fHistV0Timing;
  TH2F* fHistV0TimingPhysSel;

  TH1F* fHistZvertexPosition; 

  TH2F* fHistdEdxElesigma;
  
  TH1F* fHistdEventClassUnitRap;
  TH1F* fHistdEventClassGapRap1;
  TH1F* fHistdEventClassGapRap2;
  TH1F* fHistdEventClassV0AC;
  
  TH1F* fHistdEventClassPHIUnitRap;
  TH1F* fHistdEventClassPHIGapRap1;
  TH1F* fHistdEventClassPHIGapRap2;
  TH1F* fHistdEventClassPHIV0AC;
  
  TH1F* fHistSPDMultiUnitRapV0ACEventClass[10];
  TH1F* fHistSPDMultiGapRap1V0ACEventClass[10];
  TH1F* fHistSPDMultiGapRap2V0ACEventClass[10];
  
  TH1F* fHistV0ACMultiSPDUnitRapEventClass[10];
  TH1F* fHistV0ACMultiSPDGapRap1EventClass[10];
  TH1F* fHistV0ACMultiSPDGapRap2EventClass[10];

  TH1F* fHistPi0UnitRapV0AC[10];
  TH1F* fHistPi0AcceptV0AC[10] ;
  TH1F* fHistPi0UnitRapSPDUnitRap[10];
  TH1F* fHistPi0AcceptSPDUnitRap[10] ;
  TH1F* fHistPi0UnitRapSPDGapRap1[10];
  TH1F* fHistPi0AcceptSPDGapRap1[10] ;
  TH1F* fHistPi0UnitRapSPDGapRap2[10];
  TH1F* fHistPi0AcceptSPDGapRap2[10] ;

  TH1F* fHistGammaUnitRap;
  TH1F* fHistGammaPi0UnitRap;
  TH1F* fHistGammaEtaUnitRap;
  TH1F* fHistGammaOmegaUnitRap;
  TH1F* fHistGammaEtaPrimeUnitRap;
  TH1F* fHistGammaPhiUnitRap;
  TH1F* fHistGammaRhoUnitRap;
  TH1F* fHistGammaOtherUnitRap;
  TH1F* fHistGammaAccept;
  TH1F* fHistGammaPi0Accept;
  TH1F* fHistGammaEtaAccept;
  TH1F* fHistGammaOmegaAccept;
  TH1F* fHistGammaEtaPrimeAccept;
  TH1F* fHistGammaPhiAccept;
  TH1F* fHistGammaRhoAccept;
  TH1F* fHistGammaOtherAccept;

  TH1F* fHistRadiusPi0;
  TH2F* fHistRadius2DPi0;
  TH2F* fHistRadius2DRecPi0;

  TH2F* fHistPHOSAcceptance;

  ClassDef(AliAnalysisTaskPHOSTrigPi0, 3); // PHOS analysis task
};

#endif
