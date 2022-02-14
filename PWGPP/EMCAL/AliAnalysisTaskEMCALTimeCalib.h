#ifndef ALIANALYSISTASKEMCALTIMECALIB_H
#define ALIANALYSISTASKEMCALTIMECALIB_H

//_________________________________________________________________________
/// \class AliAnalysisTaskEMCALTimeCalib
/// \ingroup EMCALPerformance 
/// \brief Task to work on Time Calibration for EMCal/DCal.
///
/// Derived from "Exercice with a task to work on T0 from TOF or T0" by Hugues Delagrange (SUBATECH)"
/// History:
/// TimeTaskMB2_v2     wrt TimeTaskMB add histos time_vs_energy
/// TimeTaskMB2_v2b    (memory optimization.. t_vs_E splitted in hi and low res  part)
/// TimeTaskMB2_v2m    memory leak patch: sterilized ComputeT0TOF.. temporarly it is not called (to be understood and fixed
/// TimeTaskMB2_v3m   fix "dÃ¨calage" in the binning 
/// TimeTaskMB2_v4a   Implement single cell histos+ correction from average computed on run 128503 
/// TimeTaskMB2_v4b   Implement BC offsetCorrectio + histo to look at partial RCU
/// TimeTaskMB2_v4d   histoaming convention ""=raw  "Corr"= cell/cell correction "BCCorr"= cell/cell+BC 
/// TimeTaskMB2_v4e   check of BCnb vs time vs SM new histos
/// TimeTaskMB2_v4f   new histos for time for each BC 
/// TimeTaskMB2_v5   new histos for time for each BC calc average for each BC -> new ref  
///
/// 2015.06.03 Extended to DCal, added setters and getters
/// 2015.06.11 Added AODs and lego train
/// 2015.06.29 Geometry removed from Nofify to UserCreateOutputObjects
///            Added second step of calibration with time cut
/// 2015.07.08 Added geometry check in Notify and set once, 
///            added T0 time from TOF histos
///            added reference file protection
/// 2015.07.14 corrected geometry pointer in terminate
///            added reference histograms for low gain
///            modified loading of reference file
///            modified default parameters in AddTask
/// 2015.07.17 additional selection criteria for clusters with low gain cell 
/// 2015.10.09 Modification to be consistent with OADB. Info about calibration constant
///            for cell with absId=0 is kept in 'underflow bin' 
/// 2015.11.11 Modification to calibrate run by run by additional L1 phase
/// 2015.11.19 Added histogram settings
/// 2016.01.18 L1 phases read once at the beginning of lego train
/// 2016.01.23 revert L1 phase and 100ns shift for runs before LHC15n muon calo pass1
/// 2016.01.28 correction for L1 phase shift for runs before LHC15n muon calo pass1
/// 2016.02.02 added flag to fill heavy histograms
/// 2016.02.03 added bad channel map
/// 2016.02.08 added control histograms for low gain separatelly
/// 2017.06.16 added correct triggers + calibratin on most ene cell in cluster
///
/// \author Hugues Delagrange+, SUBATECH
/// \author Marie Germain <marie.germain@subatech.in2p3.fr>, SUBATECH
/// \author Adam Matyja <adam.tomasz.matyja@ifj.edu.pl>, INP PAN Cracow
/// \date Jun 3, 2015
//_________________________________________________________________________

class TH1F;
class TH1D;
class TH2D;
class TH1C;

//class AliESDEvent;
//class AliESDCaloCluster;
//class AliAODCaloCluster;
class AliVCluster;
//class AliAODEvent;
class AliVEvent;
//class AliTOFT0maker;

#include "AliAnalysisTaskSE.h"
#include <fstream>

class AliAnalysisTaskEMCALTimeCalib : public AliAnalysisTaskSE 
{
 public:

   enum { kNSM = 20, kNBCmask = 4 };

   AliAnalysisTaskEMCALTimeCalib() : AliAnalysisTaskSE(),
    fPARvec(),
    fCurrentPARs(),
    fCurrentPARIndex(0),
    fIsPARRun(0),
    fRunNumber(-1),
    fOutputList(0),
    fgeom(0),
    fGeometryName(0),
    fMinClusterEnergy(0),
    fMaxClusterEnergy(0),
    fMinNcells(0),
    fMaxNcells(0),
    fMinLambda0(0),
    fMaxLambda0(0),
    fMinLambda0LG(0),
    fMaxLambda0LG(0),
    fMaxRtrack(0),
    fMinCellEnergy(0),
    fReferenceFileName(0),
    fReferenceRunByRunFileName(0),
    fPileupFromSPD(kFALSE),
    fMinTime(0),
    fMaxTime(0),
    fMostEneCellOnly(kFALSE),
    fRawTimeNbins (0),
    fRawTimeMin   (0),
    fRawTimeMax   (0),
    fPassTimeNbins(0),
    fPassTimeMin  (0),
    fPassTimeMax  (0),
    fEnergyNbins  (0),
    fEnergyMin(0),
    fEnergyMax(0),
    fEnergyLGNbins(0),
    fEnergyLGMin(0),
    fEnergyLGMax(0),
    fFineNbins(0),
    fFineTmin(0),
    fFineTmax(0),
    fL1PhaseList(0),
    fBadReco(kFALSE),
    fFillHeavyHisto(kFALSE),
    fOneHistAllBCs(kFALSE),
    fTimeECorrection(kFALSE),
    fEMCALTimeEShiftCorrection(0),
    fEMCALRecalibrationFactors(NULL),
    fBadChannelMapArray(0),
    fBadChannelMapSet(kFALSE),
    fSetBadChannelMapSource(0),
    fBadChannelFileName(0),
    fhcalcEvtTime(0),
    fhEvtTimeHeader(0),
    fhEvtTimeDiff(0),
    fhEventType(0),
    fhTOFT0vsEventNumber(0),
    fhTcellvsTOFT0(0),
    fhTcellvsTOFT0HD(0),
    fhTcellvsSM(0),
    fhEneVsAbsIdHG(0),
    fhEneVsAbsIdLG(0),
    fhTimeVsBC(0),
    fhTimeSumSq(),
    fhTimeEnt(),
    fhTimeSum(),
    fhTimeSumSqAllBCs(0x0),
    fhTimeEntAllBCs(0x0),
    fhTimeSumAllBCs(0x0),
    fhTimeLGSumSq(),
    fhTimeLGEnt(),
    fhTimeLGSum(),
    fhTimeLGSumSqAllBCs(0x0),
    fhTimeLGEntAllBCs(0x0),
    fhTimeLGSumAllBCs(0x0),
    fhAllAverageBC(),
    fhAllAverageLGBC(),
    fhAllAverageAllBCs(0x0),
    fhAllAverageLGAllBCs(0x0),
    fhRefRuns(0),
    fhTimeDsup(),
    fhTimeDsupBC(),
    fhTimeDsupLG(),
    fhTimeDsupLGBC(),
    fhRawTimeVsIdBC(),
    fhRawTimeSumBC(),
    fhRawTimeEntriesBC(),
    fhRawTimeSumSqBC(),
    fhRawTimeVsIdLGBC(),
    fhRawTimeSumLGBC(),
    fhRawTimeEntriesLGBC(),
    fhRawTimeSumSqLGBC(),
    fhRawTimePARs(),
    fhRawTimeLGPARs(),
    fhRawCorrTimeVsIdBC(),
    fhRawCorrTimeVsIdLGBC(),
    fhTimeVsIdBC(),
    fhTimeVsIdLGBC(),
    fhTimeVsIdAllBCs(0x0),
    fhTimeVsIdLGAllBCs(0x0)
    { ; }
  
  AliAnalysisTaskEMCALTimeCalib(const char *name);
  virtual ~AliAnalysisTaskEMCALTimeCalib() { ; }
  
  // struct for storing PAR info
  struct PARInfo {
      Int_t runNumber;
      Int_t numPARs;
      std::vector<ULong64_t> PARGlobalBCs;
      PARInfo() : runNumber(0), numPARs(0), PARGlobalBCs() {}
  };

  //  virtual void   LocalInit();
  //virtual Bool_t Notify();
  virtual void   NotifyRun();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // Getters and setters
  Int_t    GetRunNumber()         { return fRunNumber         ; }
  Double_t GetMinClusterEnergy()  { return fMinClusterEnergy  ; }
  Double_t GetMaxClusterEnergy()  { return fMaxClusterEnergy  ; }
  Int_t    GetMinNcells()         { return fMinNcells         ; }
  Int_t    GetMaxNcells()         { return fMaxNcells         ; }
  Double_t GetMinLambda0()        { return fMinLambda0        ; }
  Double_t GetMaxLambda0()        { return fMaxLambda0        ; }
  Double_t GetMinLambda0LG()      { return fMinLambda0LG      ; }
  Double_t GetMaxLambda0LG()      { return fMaxLambda0LG      ; }
  Double_t GetMaxRtrack()         { return fMaxRtrack         ; }
  Double_t GetMinCellEnergy()     { return fMinCellEnergy     ; }
  TString  GetReferenceFileName() { return fReferenceFileName ; }
  TString  GetReferenceRunByRunFileName() { return fReferenceRunByRunFileName ; }
  TString  GetGeometryName()      { return fGeometryName      ; }
  Double_t GetMinTime()           { return fMinTime           ; }
  Double_t GetMaxTime()           { return fMaxTime           ; }
  TString  GetBadChannelFileName() { return fBadChannelFileName ; }

  void SetRunNumber        (Int_t    v) { fRunNumber        = v ; }
  void SetMinClusterEnergy (Double_t v) { fMinClusterEnergy = v ; }
  void SetMaxClusterEnergy (Double_t v) { fMaxClusterEnergy = v ; }
  void SetMinNcells        (Int_t    v) { fMinNcells        = v ; }  
  void SetMaxNcells        (Int_t    v) { fMaxNcells        = v ; }   
  void SetMinLambda0       (Double_t v) { fMinLambda0       = v ; }	   
  void SetMaxLambda0       (Double_t v) { fMaxLambda0       = v ; }	   
  void SetMinLambda0LG     (Double_t v) { fMinLambda0LG     = v ; }	   
  void SetMaxLambda0LG     (Double_t v) { fMaxLambda0LG     = v ; }	   
  void SetMaxRtrack        (Double_t v) { fMaxRtrack        = v ; }	   
  void SetMinCellEnergy    (Double_t v) { fMinCellEnergy    = v ; }	   
  void SetReferenceFileName(TString  v) { fReferenceFileName= v ; }
  void SetReferenceRunByRunFileName(TString  v) { fReferenceRunByRunFileName= v ; }
  void SetGeometryName     (TString  v) { fGeometryName     = v ; }
  void SetMinTime          (Double_t v) { fMinTime          = v ; }	   
  void SetMaxTime          (Double_t v) { fMaxTime          = v ; }	
  void SetBadChannelFileName(TString  v) { fBadChannelFileName = v ; }

  //histogram settings
  void SetRawTimeHisto (Int_t nbins,Double_t lower,Double_t upper) { 
    fRawTimeNbins = nbins;
    fRawTimeMin = lower ;
    fRawTimeMax = upper ;
  }
  void SetPassTimeHisto (Int_t nbins,Double_t lower,Double_t upper) { 
    fPassTimeNbins = nbins;
    fPassTimeMin = lower ;
    fPassTimeMax = upper ;
  }
  void SetEnergyHistoHG (Int_t nbins,Double_t lower,Double_t upper) { 
    fEnergyNbins = nbins;
    fEnergyMin = lower ;
    fEnergyMax = upper ;
  }
  void SetEnergyHistoLG (Int_t nbins,Double_t lower,Double_t upper) { 
    fEnergyLGNbins = nbins;
    fEnergyLGMin = lower ;
    fEnergyLGMax = upper ;
  }

  void SetFineT0Histo (Int_t nbins,Double_t lower,Double_t upper) { 
    fFineNbins = nbins;
    fFineTmin = lower ;
    fFineTmax = upper ;
  }


  // Switches
  void SwitchOnPileupFromSPD()  { fPileupFromSPD = kTRUE ; }
  void SwitchOffPileupFromSPD() { fPileupFromSPD = kFALSE ; }

  void SwitchOnMostEneCellOnly()  { fMostEneCellOnly = kTRUE ; }
  void SwitchOffMostEneCellOnly() { fMostEneCellOnly = kFALSE ; }
  
  void SwitchOnBadReco()  { fBadReco = kTRUE ; }
  void SwitchOffBadReco() { fBadReco = kFALSE ; }

  void SwithOnFillHeavyHisto()  { fFillHeavyHisto = kTRUE ; }
  void SwithOffFillHeavyHisto() { fFillHeavyHisto = kFALSE ; }

  void SetBadChannelMapSource(Int_t v) { fSetBadChannelMapSource = v ; }
  Int_t GetBadChannelMapSource() { return fSetBadChannelMapSource ; }

  void SetDefaultCuts();
  void LoadReferenceHistos(); //loaded once per period to the memory

  void LoadReferenceRunByRunHistos(); //loaded once to the memory at the beginning, phases for all runs 
  void SetL1PhaseReferenceForGivenRun();//set refernce L1phase per run 

  void SetL1PhaseReferencePAR();//set reference L1phase for specific PAR in run
  void SetPARInfo(TString PARfilename);//for given run, load in PAR info from file

  void LoadBadChannelMap(); //load bad channel map, main
  void LoadBadChannelMapFile(); //load bad channel map from file
  void LoadBadChannelMapOADB();//load bad channel map from OADB
  Int_t GetEMCALChannelStatus(Int_t iSM , Int_t iCol, Int_t iRow) const { 
    if(fBadChannelMapArray) return (Int_t) ((TH2I*)fBadChannelMapArray->At(iSM))->GetBinContent(iCol,iRow); 
    else return 0;}//Channel is ok by default
  Int_t GetEMCALChannelStatus(Int_t absId) const {
    if(fBadChannelMapArray) return (Int_t) ((TH2I*)fBadChannelMapArray->At(0))->GetBinContent(absId+1);
    else return 0;}//Channel is ok by default

  static void ProduceCalibConsts(TString inputFile="time186319testWOL0.root",TString outputFile="Reference.root",Bool_t isFinal=kFALSE, Bool_t oneHistoAllBCs=kFALSE, Bool_t isPAR=kFALSE, Bool_t doFit=kFALSE);
  static void ProduceOffsetForSMsV2(Int_t runNumber,TString inputFile="Reference.root",TString outputFile="ReferenceSM.root",Bool_t offset100=kTRUE, Bool_t justL1phase=kTRUE,TString PARfilename="");

  void SwithOnFillOneHistAllBCs()  { fOneHistAllBCs = kTRUE ; }
  void SwitchOnTimeECorrection()  { fTimeECorrection = kTRUE  ; }
  void SwitchOffTimeECorrection()  { fTimeECorrection = kFALSE  ; }

  void CorrectCellTimeVsE(Float_t energy, Float_t & celltime, Bool_t isHighGain) const;
  Double_t GetLowGainSlewing(Double_t energy) const;

  private:
  
  // variables and functions needed for PAR handling
  std::vector<PARInfo> fPARvec; ///< vector of PAR info for all runs
  PARInfo fCurrentPARs;         //! Par Info for current Run Number
  Int_t fCurrentPARIndex;       //! Which PAR the currnt event is after
  Bool_t fIsPARRun;             //! Does current run have PAR info? 
  void GetPARInfoForRunNumber(Int_t runnum);

  
  virtual void PrepareTOFT0maker();
  Bool_t SetEMCalGeometry();
  Bool_t AcceptCluster(AliVCluster* clus);
  Bool_t CheckCellRCU(Int_t nSupMod,Int_t icol,Int_t irow);
  Bool_t IsLowGainCellInCluster(AliVCluster* clus);

  Int_t InitEDepTimeCalibration();
  Int_t InitRecalib();
  Float_t GetEMCALChannelRecalibrationFactor(Int_t iSM , Int_t iCol, Int_t iRow) const;
  void SetEMCALChannelRecalibrationFactors(Int_t iSM , const TH2F* h);

  TString GetPass();

  // data members
  Int_t          fRunNumber ; //!<! run number
  
//  /// pointer to get T0 from TOF
//  AliTOFT0maker *fTOFmaker;   //->
  
  /// pointer to output list
  TList         *fOutputList; //->
  
  /// pointer to EMCal geometry
  AliEMCALGeometry *fgeom;              ///< EMCAL geometry
  TString        fGeometryName ;        ///< geometry name
   
  // setable variables for cuts
  
  Double_t       fMinClusterEnergy ;    ///< minimum cluster energy
  Double_t       fMaxClusterEnergy ;    ///< maximum cluster energy
  
  Int_t          fMinNcells ;           ///< minimum number of cells in cluster
  Int_t          fMaxNcells ;           ///< maximum number of cells in cluster
  
  Double_t       fMinLambda0 ;          ///< minimum cluster lambda0
  Double_t       fMaxLambda0 ;          ///< maximum cluster lambda0
  Double_t       fMinLambda0LG ;        ///< minimum cluster lambda0 Low Gain
  Double_t       fMaxLambda0LG ;        ///< maximum cluster lambda0 Low Gain
  
  Double_t       fMaxRtrack ;           ///< maximum cluster track distance
  
  Double_t       fMinCellEnergy;        ///< minimum cell energy

  TString        fReferenceFileName ;   //!<! name of reference file (for one period)
  TString        fReferenceRunByRunFileName ;   ///< name of reference file (run-by-run)

  Bool_t         fPileupFromSPD ;       ///< flag to set PileupFromSPD

  Double_t       fMinTime ;             ///< minimum cluster time after correction
  Double_t       fMaxTime ;             ///< maximum cluster time after correction

  Bool_t         fMostEneCellOnly ;     ///< flag to use calibration on most energetic cell in cluster only
  
  //histogram settings
  Int_t          fRawTimeNbins ;        ///< number of bins of histo with raw time
  Double_t       fRawTimeMin   ;        ///< lower range of histo with raw time
  Double_t       fRawTimeMax   ;        ///< upper range of histo with raw time
  Int_t          fPassTimeNbins;        ///< number of bins of histo with time in passX
  Double_t       fPassTimeMin  ;        ///< lower range of histo with time in passX
  Double_t       fPassTimeMax  ;        ///< upper range of histo with time in passX
  Int_t          fEnergyNbins  ;        ///< number of bins of histo with energy HG
  Double_t       fEnergyMin    ;        ///< lower range of histo with energy	 HG
  Double_t       fEnergyMax    ;        ///< upper range of histo with energy    HG
  Int_t          fEnergyLGNbins;        ///< number of bins of histo with energy LG
  Double_t       fEnergyLGMin  ;        ///< lower range of histo with energy	 LG
  Double_t       fEnergyLGMax  ;        ///< upper range of histo with energy    LG
  Int_t          fFineNbins    ;        ///< number of bins of histo with T0 time
  Double_t       fFineTmin     ;        ///< lower range of histo with T0 time
  Double_t       fFineTmax     ;        ///< upper range of histo with T0 time

  TObjArray     *fL1PhaseList;          ///< array with phases for set of runs 
  Bool_t         fBadReco;              ///< flag to apply 100ns shift and L1 shift

  Bool_t         fFillHeavyHisto;       ///< flag to fill heavy histograms

  Bool_t	 fOneHistAllBCs;		///< flag to use one histogram for all the BCs instead of four
  Bool_t   fTimeECorrection;  ///< Switch on or off the energy dependent time recalibration

  TSpline3*  fEMCALTimeEShiftCorrection;  ///< Spline to correct energy dependent time shift for high gain cells
  TObjArray* fEMCALRecalibrationFactors;  ///< Array of histograms with map of recalibration factors, EMCAL

  // bad channel map
  TObjArray     *fBadChannelMapArray;   ///< bad channel map array
  Bool_t         fBadChannelMapSet;     ///< flag whether bad channel map is set
  Int_t          fSetBadChannelMapSource;///< switch to load BC map 0-no BC,1-OADB,2-file
  TString        fBadChannelFileName ;  ///< name of file with bad channels

  // histograms
  TH1F          *fhcalcEvtTime;         //!<! spectrum calcolot0[0]
  TH1F          *fhEvtTimeHeader;       //!<! spectrum time from header
  TH1F          *fhEvtTimeDiff;         //!<! spectrum time difference
  TH1F          *fhEventType;           //!<! event type
  TH1F          *fhTOFT0vsEventNumber;  //!<! TOF T0 evolution as a function of time 
  TH2F          *fhTcellvsTOFT0;        //!<! time of cell vs TOF T0 time
  TH2F          *fhTcellvsTOFT0HD;      //!<! time of cell vs TOF T0 time for higher energy threshold
  TH2F          *fhTcellvsSM;           //!<! cell time vs SM
  TH2F          *fhEneVsAbsIdHG;        //!<! energy of each cell for high gain cells with strange time
  TH2F          *fhEneVsAbsIdLG;        //!<! energy of each cell for low gain cells with strange time
  TH2F          *fhTimeVsBC;            //!<!cell time vs BC

  // histos for storing the time values per cell for further averaging;
  TH1F		*fhTimeSumSq  [kNBCmask]; //!<!  4
  TH1F		*fhTimeEnt    [kNBCmask]; //!<!  4
  TH1F		*fhTimeSum    [kNBCmask]; //!<!  4
  TH1F		*fhTimeSumSqAllBCs  	; //!
  TH1F		*fhTimeEntAllBCs    	; //!
  TH1F		*fhTimeSumAllBCs    	; //!
  TH1F		*fhTimeLGSumSq[kNBCmask]; //!<!  4
  TH1F		*fhTimeLGEnt  [kNBCmask]; //!<!  4
  TH1F		*fhTimeLGSum  [kNBCmask]; //!<!  4
  TH1F		*fhTimeLGSumSqAllBCs  	; //!
  TH1F		*fhTimeLGEntAllBCs    	; //!
  TH1F		*fhTimeLGSumAllBCs    	; //!

  // histos with reference values after the first iteration  
  TH1F		*fhAllAverageBC   [kNBCmask]; ///> 4 BCmask High gain
  TH1F		*fhAllAverageLGBC [kNBCmask]; ///> 4 BCmask Low gain

  TH1S		*fhAllAverageAllBCs	; ///> High gain
  TH1S		*fhAllAverageLGAllBCs   ; ///> Low gain

  // histo with reference values run-by-run after the first iteration 
  TH1C		*fhRefRuns; ///< 20 entries per run: nSM

  // control histos
  TH2F		*fhTimeDsup  [kNSM];            //!<! 20 SM high gain
  TH2F		*fhTimeDsupBC[kNSM][kNBCmask];  //!<! 20 x 4 high gain
  TH2F		*fhTimeDsupLG  [kNSM];            //!<! 20 SM low gain
  TH2F		*fhTimeDsupLGBC[kNSM][kNBCmask];  //!<! 20 x 4 low gain

  //main histos for raw time
  TH2F          *fhRawTimeVsIdBC     [kNBCmask]; //!<! 4 BCmask HG
  TH1F          *fhRawTimeSumBC      [kNBCmask]; //!<! 4 BCmask HG
  TH1F          *fhRawTimeEntriesBC  [kNBCmask]; //!<! 4 BCmask HG
  TH1F          *fhRawTimeSumSqBC    [kNBCmask]; //!<! 4 BCmask HG
  TH2F          *fhRawTimeVsIdLGBC   [kNBCmask]; //!<! 4 BCmask LG
  TH1F          *fhRawTimeSumLGBC    [kNBCmask]; //!<! 4 BCmask LG
  TH1F          *fhRawTimeEntriesLGBC[kNBCmask]; //!<! 4 BCmask LG
  TH1F          *fhRawTimeSumSqLGBC  [kNBCmask]; //!<! 4 BCmask LG

  //histos for correction of Raw Time with PAR
  std::vector<std::vector<TH2F*>> fhRawTimePARs;//!<!
  std::vector<std::vector<TH2F*>> fhRawTimeLGPARs;//!<!

  //histos for raw time after wrong reconstruction correction (100ns and L1 phase)
  TH2F          *fhRawCorrTimeVsIdBC  [kNBCmask]; //!<! 4 BCmask HG
  TH2F          *fhRawCorrTimeVsIdLGBC[kNBCmask]; //!<! 4 BCmask LG

  //histos for raw time after wrong reconstruction correction (100ns and L1 phase) and new L1 phase
  TH2F          *fhTimeVsIdBC  [kNBCmask]; //!<! 4 BCmask HG
  TH2F          *fhTimeVsIdLGBC[kNBCmask]; //!<! 4 BCmask LG
  TH2F          *fhTimeVsIdAllBCs        ; //! HG
  TH2F          *fhTimeVsIdLGAllBCs      ; //! LG

  /// Copy constructor not implemented.
  AliAnalysisTaskEMCALTimeCalib(           const AliAnalysisTaskEMCALTimeCalib&);
  
  /// Assignment operator not implemented.
  AliAnalysisTaskEMCALTimeCalib& operator=(const AliAnalysisTaskEMCALTimeCalib&); 
  
/// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALTimeCalib, 7) ;
/// \endcond
};

#endif
