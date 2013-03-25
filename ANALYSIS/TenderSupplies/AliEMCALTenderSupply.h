#ifndef ALIEMCALTENDERSUPPLY_H
#define ALIEMCALTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  EMCAL tender, apply corrections to EMCAl clusters                 //
//  and do track matching.                                            //
//  Author: Deepa Thomas (Utrecht University)                         //
//  Later mods/rewrite: Jiri Kral (University of Jyvaskyla)           //
//  S. Aiola, C. Loizides : Make it work for AODs                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliTenderSupply.h"

class TTree;
class TClonesArray;

class AliVCluster;
class AliEMCALRecoUtils;
class AliEMCALGeometry;
class TGeoHMatrix;
class TTree;
class TFile;
class TString;
class AliEMCALClusterizer;
class AliEMCALAfterBurnerUF;
class AliEMCALRecParam;
class AliAnalysisTaskSE;
class AliVEvent;

class AliEMCALTenderSupply: public AliTenderSupply {
  
public:
  AliEMCALTenderSupply();
  AliEMCALTenderSupply(const char *name, const AliTender *tender=NULL);
  AliEMCALTenderSupply(const char *name, AliAnalysisTaskSE *task);
  virtual ~AliEMCALTenderSupply();

  enum NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3,kBeamTest=4,kBeamTestCorrected=5};
  enum MisalignSettings{kdefault=0,kSurveybyS=1,kSurveybyM=2};

  virtual void Init();
  virtual void ProcessEvent();

  void     SetTask(AliAnalysisTaskSE *task)               { fTask = task                     ;}
  void     SetDefaults();

  void     SetEMCALGeometryName(const char *name)         { fEMCALGeoName = name             ;}
  TString  EMCALGeometryName()                      const { return fEMCALGeoName             ;}

  void     SetDebugLevel(Int_t level)                     { fDebugLevel=level                ;}

  void     SetBasePath(const Char_t *basePath)            { fBasePath = basePath             ;}
 
  void     SetConfigFileName(const char *name)            { fConfigName = name               ;}

  void     SetNonLinearityFunction(Int_t fun)             { fNonLinearFunc = fun             ;}
  Int_t    GetNonLinearityFunction() const                { return fNonLinearFunc            ;}
  
  void     SetNonLinearityThreshold(Int_t threshold)      { fNonLinearThreshold = threshold  ;} //only for Alexei's non linearity correction
  Int_t    GetNonLinearityThreshold()               const { return fNonLinearThreshold       ;}

  void     SwitchOnNonLinearityCorrection()               { fDoNonLinearity = kTRUE          ;}
  void     SwitchOffNonLinearityCorrection()              { fDoNonLinearity = kFALSE         ;}

  void     SwitchOnReCalibrateCluster()                   { fReCalibCluster = kTRUE          ;}
  void     SwitchOffReCalibrateCluster()                  { fReCalibCluster = kFALSE         ;}

  void     SwitchOnRecalculateClusPos()                   { fRecalClusPos = kTRUE            ;}
  void     SwitchOffRecalculateClusPos()                  { fRecalClusPos = kFALSE           ;}

  void     SetMisalignmentMatrixSurvey(Int_t misalignSurvey) { fMisalignSurvey = misalignSurvey ;}
  Int_t    GetMisalignmentMatrixSurvey() const               { return fMisalignSurvey           ;}    

  void     SwitchOnCellFiducialRegion()                   { fFiducial = kTRUE                ;}
  void     SwitchOffCellFiducialRegion()                  { fFiducial = kFALSE               ;}

  void     SetNumberOfCellsFromEMCALBorder(Int_t n)       { fNCellsFromEMCALBorder = n       ;}
  Int_t    GetNumberOfCellsFromEMCALBorder()        const { return fNCellsFromEMCALBorder    ;}

  void     SwitchOnRecalDistBadChannel()                  { fRecalDistToBadChannels = kTRUE  ;}
  void     SwitchOffRecalDistBadChannel()                 { fRecalDistToBadChannels = kFALSE ;}

  void     SwitchOnRecalShowerShape()                     { fRecalShowerShape = kTRUE        ;}
  void     SwitchOffRecalShowerShape()                    { fRecalShowerShape = kFALSE       ;}

  Float_t  GetRCut()                                const { return fRcut                     ;}
  void     SetRCut(Float_t rcut)                          { fRcut = rcut                     ;}

  Double_t GetMass()                                const { return fMass                     ;}
  void     SetMass(Double_t mass)                         { fMass = mass                     ;}

  Double_t GetStep()                                const { return fStep                     ;}
  void     SetStep(Double_t step)                         { fStep = step                     ;}

  Double_t GetEtaCut()                              const { return fEtacut                   ;}
  void     SetEtaCut(Double_t eta)                        { fEtacut = eta                    ;}

  Double_t GetPhiCut()                              const { return fPhicut                   ;}
  void     SetPhiCut(Double_t phi)                        { fPhicut = phi                    ;}

  Float_t  GetExoticCellFraction()                  const { return fExoticCellFraction       ;}
  void     SetExoticCellFraction(Float_t f)               { fExoticCellFraction = f          ;}

  Float_t  GetExoticCellDiffTime()                  const { return fExoticCellDiffTime       ;}
  void     SetExoticCellDiffTime(Float_t f)               { fExoticCellDiffTime = f          ;}

  Float_t  GetExoticCellMinAmplitude()              const { return fExoticCellMinAmplitude   ;}
  void     SetExoticCellMinAmplitude(Float_t f)           { fExoticCellMinAmplitude = f      ;}

  void     SwitchOnReclustering()                         { fReClusterize = kTRUE            ;}
  void     SwitchOffReclustering()                        { fReClusterize = kFALSE           ;}

  void     SwitchOnCutEtaPhiSum()                         { fCutEtaPhiSum=kTRUE;      fCutEtaPhiSeparate=kFALSE ;}
  void     SwitchOnCutEtaPhiSeparate()                    { fCutEtaPhiSeparate=kTRUE; fCutEtaPhiSum=kFALSE      ;}
  
  void     SwitchOnLoadOwnGeometryMatrices()              { fLoadGeomMatrices = kTRUE        ;}
  void     SwitchOffLoadOwnGeometryMatrices()             { fLoadGeomMatrices = kFALSE       ;}
  void     SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fEMCALMatrix[i]    = m           ;}
  
  AliEMCALRecParam   *GetRecParam() const                 { return fRecParam                 ;} 
  void                SetRecParam(AliEMCALRecParam *p)    { fRecParam = p                    ;}
 
  AliEMCALRecoUtils  *GetRecoUtils() const                { return fEMCALRecoUtils           ;}
  void                SetRecoUtils(AliEMCALRecoUtils *r)  { fEMCALRecoUtils = r              ;}
  
  void     SwitchOnUpdateCell()                           { fUpdateCell = kTRUE              ;} 
  void     SwitchOffUpdateCell()                          { fUpdateCell = kFALSE             ;}  

  void     SwitchOnBadCellRemove()                        { fBadCellRemove = kTRUE           ;} 
  void     SwitchOffBadCellRemove()                       { fBadCellRemove = kFALSE          ;}  

  void     SwitchOnClusterBadChannelCheck()               { fClusterBadChannelCheck = kTRUE  ;} 
  void     SwitchOffClusterBadChannelCheck()              { fClusterBadChannelCheck = kFALSE ;}  

  void     SwitchOnExoticCellRemove()                     { fRejectExoticCells = kTRUE       ;} 
  void     SwitchOffExoticCellRemove()                    { fRejectExoticCells = kFALSE      ;}  

  void     SwitchOnClusterExoticChannelCheck()            { fRejectExoticClusters = kTRUE    ;} 
  void     SwitchOffClusterExoticChannelCheck()           { fRejectExoticClusters = kFALSE   ;}  

  void     SwitchOnCalibrateEnergy()                      { fCalibrateEnergy = kTRUE         ;} 
  void     SwitchOffCalibrateEnergy()                     { fCalibrateEnergy = kFALSE        ;}  

  void     SwitchOnCalibrateTime()                        { fCalibrateTime = kTRUE           ;} 
  void     SwitchOffCalibrateTime()                       { fCalibrateTime = kFALSE          ;}  

  void     SwitchOnUpdateCellOnly()                       { fDoUpdateOnly = kTRUE            ;} 
  void     SwitchOffUpdateCellOnly()                      { fDoUpdateOnly = kFALSE           ;}  
 
  void     SwitchOnTrackMatch()                           { fDoTrackMatch = kTRUE            ;}
  void     SwitchOffTrackMatch()                          { fDoTrackMatch = kFALSE           ;}

  void     SetPass(const char* p)                         { fFilepass = p; fGetPassFromFileName = kFALSE; }
  void     UnsetPass()                                    { fFilepass = ""; fGetPassFromFileName = kTRUE; }
  
  //MC label methods
  
  void     RemapMCLabelForAODs(Int_t &label);
  void     SwitchOnRemapMCLabelForAODs()                  { fRemapMCLabelForAODs  = kTRUE   ; }
  void     SwitchOffRemapMCLabelForAODs()                 { fRemapMCLabelForAODs  = kFALSE  ; }
  
  void     SetClustersMCLabelFromOriginalClusters() ;  
  void     SwitchOnUseClusterMCLabelForCell(Int_t opt = 2) { fSetCellMCLabelFromCluster = opt ; }
  void     SwitchOffUseClusterMCLabelForCell()             { fSetCellMCLabelFromCluster = 0   ; }

  
private:

  AliVEvent* GetEvent();
  TString    GetBeamType();
  Bool_t     RunChanged() const;
  Int_t      InitBadChannels();
  Bool_t     InitClusterization();
  Int_t      InitRecParam();
  Bool_t     InitMisalignMatrix();
  Int_t      InitRecalib();
  Int_t      InitRunDepRecalib();
  Int_t      InitTimeCalibration();
  void       Clusterize();
  void       FillDigitsArray();
  void       GetPass();
  void       RecPoints2Clusters(TClonesArray *clus);
  void       RecalibrateCells();
  void       UpdateCells();
  void       UpdateClusters();

  AliAnalysisTaskSE     *fTask;                   // analysis task
  Int_t                  fRun;                    // current run number
  AliEMCALGeometry      *fEMCALGeo;               // EMCAL geometry
  TString                fEMCALGeoName;           // name of geometry to use.
  AliEMCALRecoUtils     *fEMCALRecoUtils;         // pointer to EMCAL utilities for clusterization
  TString                fConfigName;             // name of analysis configuration file
  Int_t                  fDebugLevel;             // debug level
  Int_t                  fNonLinearFunc;          // non linearity function 
  Int_t                  fNonLinearThreshold;     // non linearity threshold value for kBeamTesh non linearity function   
  Bool_t                 fReCalibCluster;         // switch for Recalibrate clusters
  Bool_t                 fUpdateCell;             // flag cell update
  Bool_t                 fCalibrateEnergy;        // flag cell energy clibration
  Bool_t                 fCalibrateTime;          // flag cell time clSibration
  Bool_t                 fDoNonLinearity;         // nNon linearity correction flag
  Bool_t                 fBadCellRemove;          // zero bad cells
  Bool_t                 fRejectExoticCells;      // reject exotic cells
  Bool_t                 fRejectExoticClusters;   // recect clusters with exotic cells
  Bool_t                 fClusterBadChannelCheck; // check clusters for bad channels
  TGeoHMatrix           *fEMCALMatrix[12];        // geometry matrices with misalignments
  Bool_t                 fRecalClusPos;           // switch for applying missalignment
  Bool_t                 fFiducial;               // switch for checking cells in the fiducial region
  Int_t                  fNCellsFromEMCALBorder;  // number os cells from EMCAL border  
  Bool_t                 fRecalDistToBadChannels; // switch for recalculation cluster position from bad channel    
  Bool_t                 fRecalShowerShape;       // switch for recalculation of the shower shape
  TTree                 *fInputTree;              //!input data tree
  TFile                 *fInputFile;              //!input data file 
  Bool_t                 fGetPassFromFileName;    // get fFilepass from file name
  TString                fFilepass;               // input data pass number
  Double_t               fMass;                   // mass for track matching
  Double_t               fStep;                   // step size during track matching
  Bool_t                 fCutEtaPhiSum;           // swicth to apply residual cut together
  Bool_t                 fCutEtaPhiSeparate;      // swicth to apply residual cut separately
  Float_t                fRcut;                   // residual cut for track matching  
  Float_t                fEtacut;                 // eta cut for track matching  
  Float_t                fPhicut;                 // phi cut for track matching  
  TString                fBasePath;               // base folder path to get root files 
  Bool_t                 fReClusterize;           // switch for reclustering
  AliEMCALClusterizer   *fClusterizer;            //!clusterizer 
  Bool_t                 fGeomMatrixSet;          // set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices;       // matrices set from configuration, not get from geometry.root or from ESDs/AODs
  AliEMCALRecParam      *fRecParam;               // reconstruction parameters container
  Bool_t                 fDoTrackMatch;           // do track matching
  Bool_t                 fDoUpdateOnly;           // do only update of cells
  AliEMCALAfterBurnerUF *fUnfolder;               //!unfolding procedure
  TClonesArray          *fDigitsArr;              //!digits array
  TObjArray             *fClusterArr;             //!recpoints array
  Int_t                  fMisalignSurvey;         // misalignment matrix survey  
  Float_t                fExoticCellFraction;     // good cell if fraction < 1-ecross/ecell
  Float_t                fExoticCellDiffTime;     // if time of candidate to exotic and close cell is too different (in ns), it must be noisy, set amp to 0
  Float_t                fExoticCellMinAmplitude; // check for exotic only if amplitud is larger than this value

  // MC labels
  Int_t                  fOrgClusterCellId[12672]; // Array ID of cluster to wich the cell belongs in unmodified clusters
 
  Int_t                  fSetCellMCLabelFromCluster; // Use cluster MC label as cell label:
                                                     // 0 - get the MC label stored in cells
                                                     // 1 - assign to the cell the MC lable of the cluster
                                                     // 2 - get the original clusters id, add all the MC labels
  TClonesArray *         fTempClusterArr      ;      //! Temporal clusters array to recover the MC of original cluster if fSetCellMCLabelFromCluster=2
  Bool_t                 fRemapMCLabelForAODs ;      // Remap AOD cells MC label
  
  AliEMCALTenderSupply(            const AliEMCALTenderSupply&c);
  AliEMCALTenderSupply& operator= (const AliEMCALTenderSupply&c);
  
  ClassDef(AliEMCALTenderSupply, 15); // EMCAL tender task
};
#endif
