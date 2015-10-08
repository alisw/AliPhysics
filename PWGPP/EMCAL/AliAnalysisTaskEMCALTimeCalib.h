#ifndef AliAnalysisTaskEMCALTimeCalib_h
#define AliAnalysisTaskEMCALTimeCalib_h

//_________________________________________________________________________
/// \class AliAnalysisTaskEMCALTimeCalib
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
///
/// \author Hugues Delagrange+, SUBATECH
/// \author Marie Germain <marie.germain@subatech.in2p3.fr>, SUBATECH
/// \author Adam Matyja <adam.tomasz.matyja@ifj.edu.pl>, INP PAN Cracow
/// \date Jun 3, 2015
//_________________________________________________________________________

class TH1F;
class TH1D;
class TH2D;

//class AliESDEvent;
//class AliESDCaloCluster;
//class AliAODCaloCluster;
class AliVCluster;
//class AliAODEvent;
class AliVEvent;
class AliTOFT0maker;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALTimeCalib : public AliAnalysisTaskSE 
{
 public:

  enum { kNSM = 20, kNBCmask = 4 };

   AliAnalysisTaskEMCALTimeCalib() : AliAnalysisTaskSE(),
    fRunNumber(-1),
    fTOFmaker(0),
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
    fReferenceFileName(),
    fPileupFromSPD(kFALSE),
    fMinTime(0),
    fMaxTime(0),
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
    fhTimeLGSumSq(),
    fhTimeLGEnt(),
    fhTimeLGSum(),
    fhAllAverageBC(),
    fhAllAverageLGBC(),
    fhTimeDsup(),
    fhTimeDsupBC(),
    fhRawTimeVsIdBC(),
    fhRawTimeSumBC(),
    fhRawTimeEntriesBC(),
    fhRawTimeSumSqBC(),
    fhRawTimeVsIdLGBC(),
    fhRawTimeSumLGBC(),
    fhRawTimeEntriesLGBC(),
    fhRawTimeSumSqLGBC()
    { ; }
  
  AliAnalysisTaskEMCALTimeCalib(const char *name);
  virtual ~AliAnalysisTaskEMCALTimeCalib() { ; }
  
  //  virtual void   LocalInit();
  virtual Bool_t Notify();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // Getters and setters
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
  TString  GetGeometryName()      { return fGeometryName      ; }
  Double_t GetMinTime()           { return fMinTime           ; }
  Double_t GetMaxTime()           { return fMaxTime           ; }

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
  void SetGeometryName     (TString  v) { fGeometryName     = v ; }
  void SetMinTime          (Double_t v) { fMinTime          = v ; }	   
  void SetMaxTime          (Double_t v) { fMaxTime          = v ; }	

  // Switches
  void SwitchOnPileupFromSPD()  { fPileupFromSPD = kTRUE ; }
  void SwitchOffPileupFromSPD() { fPileupFromSPD = kFALSE ; }


  void SetDefaultCuts();
  void LoadReferenceHistos();
  static void ProduceCalibConsts(TString inputFile="time186319testWOL0.root",TString outputFile="Reference.root",Bool_t isFinal=kFALSE);

 private:
  
  virtual void PrepareTOFT0maker();
  Bool_t SetEMCalGeometry();
  Bool_t AcceptCluster(AliVCluster* clus);
  Bool_t CheckCellRCU(Int_t nSupMod,Int_t icol,Int_t irow);
  Bool_t IsLowGainCellInCluster(AliVCluster* clus);

  // data members
  Int_t          fRunNumber ; //!<! run number
  
  /// pointer to get T0 from TOF
  AliTOFT0maker *fTOFmaker;   //->
  
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

  TString        fReferenceFileName ;   //!<! name of reference file

  Bool_t         fPileupFromSPD ;       ///< flag to set PileupFromSPD

  Double_t       fMinTime ;             ///< minimum cluster time after correction
  Double_t       fMaxTime ;             ///< maximum cluster time after correction


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
  TH1F		*fhTimeLGSumSq[kNBCmask]; //!<!  4
  TH1F		*fhTimeLGEnt  [kNBCmask]; //!<!  4
  TH1F		*fhTimeLGSum  [kNBCmask]; //!<!  4

  // histos with reference values after the first iteration  
  TH1F		*fhAllAverageBC   [kNBCmask]; ///> 4 BCmask High gain
  TH1F		*fhAllAverageLGBC [kNBCmask]; ///> 4 BCmask Low gain

  // control histos
  TH2F		*fhTimeDsup  [kNSM];            //!<! 20 SM
  TH2F		*fhTimeDsupBC[kNSM][kNBCmask];  //!<! 20 x 4

  //main histos for raw time
  TH2F          *fhRawTimeVsIdBC     [kNBCmask]; //!<! 4 BCmask HG
  TH1F          *fhRawTimeSumBC      [kNBCmask]; //!<! 4 BCmask HG
  TH1F          *fhRawTimeEntriesBC  [kNBCmask]; //!<! 4 BCmask HG
  TH1F          *fhRawTimeSumSqBC    [kNBCmask]; //!<! 4 BCmask HG
  TH2F          *fhRawTimeVsIdLGBC   [kNBCmask]; //!<! 4 BCmask LG
  TH1F          *fhRawTimeSumLGBC    [kNBCmask]; //!<! 4 BCmask LG
  TH1F          *fhRawTimeEntriesLGBC[kNBCmask]; //!<! 4 BCmask LG
  TH1F          *fhRawTimeSumSqLGBC  [kNBCmask]; //!<! 4 BCmask LG

  /// Copy constructor not implemented.
  AliAnalysisTaskEMCALTimeCalib(           const AliAnalysisTaskEMCALTimeCalib&);
  
  /// Assignment operator not implemented.
  AliAnalysisTaskEMCALTimeCalib& operator=(const AliAnalysisTaskEMCALTimeCalib&); 
  
/// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALTimeCalib, 2) ;
/// \endcond
};

#endif
