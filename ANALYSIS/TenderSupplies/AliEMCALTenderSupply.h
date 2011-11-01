#ifndef ALIEMCALTENDERSUPPLY_H
#define ALIEMCALTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  EMCAL tender, apply corrections to EMCAl clusters                 //
//  and do track matching.                                            //
//  Author: Deepa Thomas (Utrecht University)                         //
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

#include "AliEMCALRecParam.h"

class AliEMCALTenderSupply: public AliTenderSupply {
  
public:
  AliEMCALTenderSupply();
  AliEMCALTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliEMCALTenderSupply();

  enum NonlinearityFunctions{kPi0MC=0,kPi0GammaGamma=1,kPi0GammaConversion=2,kNoCorrection=3,kBeamTest=4,kBeamTestCorrected=5};
  enum MisalignSettngs{kdefault=0,kSurveybyS=1, kSurveybyM=2};

  virtual void Init();
  virtual void ProcessEvent();
  
  void     SetEMCALGeometryName(TString name)             { fEMCALGeoName = name             ;}
  TString  EMCALGeometryName()                      const { return fEMCALGeoName             ;}

  void     SetDebugLevel(Int_t level)                     { fDebugLevel=level                ;}

  void     SetBasePath(const Char_t *basePath)            { fBasePath = basePath             ;}
 
  void     SetConfigFileName(TString name)                { fConfigName = name               ;}

  void     SetNonLinearityFunction(Int_t fun)             { fNonLinearFunc = fun             ;}
  Int_t    GetNonLinearityFunction() const                { return fNonLinearFunc            ;}
  
  void     SetNonLinearityThreshold(Int_t threshold)      { fNonLinearThreshold = threshold  ;} //only for Alexei's non linearity correction
  Int_t    GetNonLinearityThreshold()               const { return fNonLinearThreshold       ;}

  void     SwitchOnReCalibrateCluster()                   { fReCalibCluster = kTRUE          ;}
  void     SwitchOffReCalibrateCluster()                  { fReCalibCluster = kFALSE         ;}

  void     SwitchOnRecalculateClusPos()                   { fRecalClusPos = kTRUE            ;}
  void     SwitchOffRecalculateClusPos()                  { fRecalClusPos = kFALSE           ;}

  void 	   SetMisalignmentMatrixSurvey(Int_t misalignSurvey) { fMisalignSurvey = misalignSurvey ;}
  Int_t	   GetMisalignmentMatrixSurvey() const               { return fMisalignSurvey           ;}	  

  void     SwitchOnCellFiducialRegion()                   { fFiducial = kTRUE                ;}
  void     SwitchOffCellFiducialRegion()                  { fFiducial = kFALSE               ;}

  void     SetNumberOfcellsFromEMCALBorder(Int_t n)       { fNCellsFromEMCALBorder = n       ;}
  Int_t    GetNumberOfcellsFromEMCALBorder()        const { return fNCellsFromEMCALBorder    ;}

  void     SwitchOnRecalDistBadChannel()                  { fRecalDistToBadChannels = kTRUE  ;}
  void     SwitchOffRecalDistBadChannel()                 { fRecalDistToBadChannels = kFALSE ;}

  Float_t  GetRCut()                                const { return fRcut                     ;}
  void     SetRCut(Float_t rcut)                          { fRcut = rcut                     ;}

  Double_t GetMass()                                const { return fMass                     ;}
  void     SetMass(Double_t mass)                         { fMass = mass                     ;}

  Double_t GetStep()                                const { return fStep                     ;}
  void     SetStep(Double_t step)                         { fStep = step                     ;}

  Double_t GetEtaCut()                              const { return fEtacut                   ;}
  void	   SetEtaCut(Double_t eta)                        { fEtacut = eta                    ;}

  Double_t GetPhiCut() 	                            const { return fPhicut                   ;}
  void	   SetPhiCut(Double_t phi)                        { fPhicut = phi                    ;}

  void     SwitchOnReclustering()                         { fReClusterize = kTRUE            ;}
  void     SwitchOffReclustering()                        { fReClusterize = kFALSE           ;}

  void     SwitchOnCutEtaPhiSum()                         { fCutEtaPhiSum=kTRUE;      fCutEtaPhiSeparate=kFALSE ;}
  void     SwitchOnCutEtaPhiSeparate()                    { fCutEtaPhiSeparate=kTRUE; fCutEtaPhiSum=kFALSE      ;}
  
  void     SwitchOnLoadOwnGeometryMatrices()              { fLoadGeomMatrices = kTRUE        ;}
  void     SwitchOffLoadOwnGeometryMatrices()             { fLoadGeomMatrices = kFALSE       ;}
  void     SetGeometryMatrixInSM(TGeoHMatrix* m, Int_t i) { fEMCALMatrix[i]    = m           ;}
  
  AliEMCALRecParam   *GetRecParam() const                 { return fRecParam                 ;} 
 
  AliEMCALRecoUtils  *GetRecoUtils() const                { return fEMCALRecoUtils           ;}

  void     SetOCDBPath(const char *path)                  { fOCDBpath = path                 ;}

  //Will update cell list by removing bad channels and recalibration + reclusterize	
  void     SwitchOnUpdateCell()                           { fUpdateCell = kTRUE              ;} 
  void     SwitchOffUpdateCell()                          { fUpdateCell = kFALSE             ;}	
 
private:

  Int_t    InitBadChannels();

  Bool_t   InitClusterization();

  void     InitRecParam();

  Bool_t   InitMisalignMatrix();

  Int_t    InitRecalib();

  void     Clusterize();

  void     FillDigitsArray();

  void     GetPass();

  void     RecPoints2Clusters(TClonesArray *clus);

  void     RecalibrateCells();

  void     UpdateCells();

  void     UpdateClusters();

  AliEMCALGeometry      *fEMCALGeo;               //! EMCAL geometry
  TString                fEMCALGeoName;           //  name of geometry to use.
  AliEMCALRecoUtils     *fEMCALRecoUtils;         //  pointer to EMCAL utilities for clusterization
  TString                fConfigName;             //  name of analysis configuration file
  Int_t                  fDebugLevel;             //  debug level
  Int_t                  fNonLinearFunc;          //  non linearity function 
  Int_t                  fNonLinearThreshold;     //  non linearity threshold value for kBeamTesh non linearity function   
  Bool_t                 fReCalibCluster;         //  switch for Recalibrate clusters
  Bool_t                 fUpdateCell;             //  Flag cell update
  TGeoHMatrix           *fEMCALMatrix[10];        //  geometry matrices with misalignments
  Bool_t                 fRecalClusPos;           //  switch for applying missalignment
  Bool_t                 fFiducial;               //  switch for checking cells in the fiducial region
  Int_t                  fNCellsFromEMCALBorder;  //  number os cells from EMCAL border	
  Bool_t                 fRecalDistToBadChannels; //  switch for recalculation cluster position from bad channel	  
  TTree                 *fInputTree;              //! input data tree
  TFile                 *fInputFile;              //! input data file 
  TString                fFilepass;               //! input data pass number
  Double_t               fMass;                   //  mass for track matching
  Double_t               fStep;                   //  step size during track matching
  Bool_t                 fCutEtaPhiSum;           //  swicth to apply residual cut together
  Bool_t                 fCutEtaPhiSeparate;      //  swicth to apply residual cut separately
  Float_t                fRcut;                   //  residual cut for track matching  
  Float_t                fEtacut;                 //  eta cut for track matching  
  Float_t                fPhicut;                 //  phi cut for track matching  
  TString                fBasePath;               //  base folder path to get root files 
  Bool_t                 fReClusterize;           //  switch for reclustering
  AliEMCALClusterizer   *fClusterizer;            //! clusterizer 
  Bool_t                 fGeomMatrixSet;          //  set geometry matrices only once, for the first event.         
  Bool_t                 fLoadGeomMatrices;       //  matrices set from configuration, not get from geometry.root or from ESDs/AODs
  AliEMCALRecParam      *fRecParam;               //  reconstruction parameters container
  Bool_t                 fRecParamSet;            //  Flag if rec param already set
  TString                fOCDBpath;               //  path with OCDB location
  AliEMCALAfterBurnerUF *fUnfolder;               //! unfolding procedure
  TClonesArray          *fDigitsArr;              //! digits array
  TObjArray             *fClusterArr;             //! recpoints array
  Int_t                  fMisalignSurvey;         //! misalignment matrix survey	

  AliEMCALTenderSupply(const AliEMCALTenderSupply&c);
  AliEMCALTenderSupply& operator= (const AliEMCALTenderSupply&c);
  
  ClassDef(AliEMCALTenderSupply, 7); // EMCAL tender task
};

#endif

