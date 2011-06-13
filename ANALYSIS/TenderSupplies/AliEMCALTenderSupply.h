#ifndef ALIEMCALTENDERSUPPLY_H
#define ALIEMCALTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  EMCAL tender, apply corrections to EMCAl clusters                 //
//  and do track matching                                             //
//  Author : Deepa Thomas (Utrecht University)                        //                      
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>

class AliVCluster;
class AliEMCALRecoUtils;
class AliEMCALGeometry;
class TGeoHMatrix;
class TTree;
class TFile;
class TString;

class AliEMCALTenderSupply: public AliTenderSupply {
  
public:
  AliEMCALTenderSupply();
  AliEMCALTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliEMCALTenderSupply();

  virtual void   Init();
  virtual void   ProcessEvent();
  
  
  void SetEMCALGeometryName(TString name)  { fEMCALGeoName = name  ;}
  TString EMCALGeometryName()       const  { return fEMCALGeoName  ;}

  void  SetDebugLevel(Int_t level)         {fDebugLevel=level      ;}

  void  SetBasePath(const Char_t *basePath) { fBasePath = basePath; }                 // mfasel: add function to set a path where to find the root files

  void  SetNonLinearityFunction(Int_t fun) { fNonLinearFunc = fun  ;}
  Int_t GetNonLinearityFunction() const    { return fNonLinearFunc ;}
  
  void  SetNonLinearityThreshold(Int_t threshold) { fNonLinearThreshold = threshold ;} //only for Alexei's non linearity correction
  Int_t GetNonLinearityThreshold() const          { return fNonLinearThreshold      ;}

  void  SwitchOnReCalibrateCluster()       { fReCalibCluster = kTRUE  ;}
  void  SwitchOffReCalibrateCluster()      { fReCalibCluster = kFALSE ;}

  void  SwitchOnReCalibrateCell()          { fReCalibCell = kTRUE  ;}
  void  SwitchOffReCalibrateCell()         { fReCalibCell = kFALSE ;}

  void  SwitchOnRecalculateClusPos()       { fRecalClusPos = kTRUE  ;}
  void  SwitchOffRecalculateClusPos()      { fRecalClusPos = kFALSE ;}

  void  SwitchOnCellFiducialRegion()       { fFiducial = kTRUE  ;}
  void  SwitchOffCellFiducialRegion()      { fFiducial = kFALSE ;}

  void  SetNumberOfcellsFromEMCALBorder(Int_t n) { fNCellsFromEMCALBorder = n    ;}
  Int_t GetNumberOfcellsFromEMCALBorder() const  { return fNCellsFromEMCALBorder ;}

  void  SwitchOnRecalDistBadChannel()      { fRecalDistToBadChannels = kTRUE  ;}
  void  SwitchOffRecalDistBadChannel()     { fRecalDistToBadChannels = kFALSE ;}

  Float_t  GetRCut()       const { return fRcut ;}
  void     SetRCut(Float_t Rcut) { fRcut = Rcut ;}

  Double_t GetMass()       const { return fMass ;}
  void     SetMass(Double_t mass){ fMass=mass   ;}

  Double_t GetStep()       const { return fStep ;}
  void     SetStep(Double_t step){ fStep=step   ;}

 
private:

  Bool_t  InitBadChannels();
  Bool_t  InitRecalib();
  Bool_t  InitMisalignMatrix();
  
  void SetClusterMatchedToTrack(AliESDEvent *event);
  void SetTracksMatchedToCluster(AliESDEvent *event);

  void RecalibrateCells();	 		//Recalibrate cells 

  AliEMCALGeometry   *fEMCALGeo;               //! EMCAL geometry
  TString       fEMCALGeoName;           // Name of geometry to use.
  
  AliEMCALRecoUtils  *fEMCALRecoUtils;         // Pointer to EMCAL utilities for clusterization
  TString       fConfigName;             // Name of analysis configuration file
  
  Int_t         fDebugLevel;             // debug level

  Int_t         fNonLinearFunc;          // Non linearity function 
  Int_t         fNonLinearThreshold;     // Non linearity threshold value for kBeamTesh non linearity function   

  Bool_t        fReCalibCluster; 	       // switch for Recalibrate clusters
  Bool_t        fReCalibCell; 	       // switch for Recalibrate cell

  TGeoHMatrix  *fEMCALMatrix[10];        // Geometry matrices with alignments
  Bool_t        fRecalClusPos; 	       // switch for applying missalignment

  Bool_t        fFiducial; 	       // switch for checking cells in the fiducial region
  Int_t         fNCellsFromEMCALBorder;  // number os cells from EMCAL border	
  Bool_t        fRecalDistToBadChannels; // switch for recalculation cluster position from bad channel	  

  TTree        *fInputTree;              // input data tree
  TFile        *fInputFile;              // input data file 
  TString       fFilepass;               // input data pass number

  Double_t      fMass;                   // mass for track matching
  Double_t      fStep;                   // step size during track matching
  Float_t       fRcut;                   // residual cut for track matching  
  
  const Char_t* fBasePath;                   // mfasel: Base Folder path to get root files 

  AliEMCALTenderSupply(const AliEMCALTenderSupply&c);
  AliEMCALTenderSupply& operator= (const AliEMCALTenderSupply&c);
  
  ClassDef(AliEMCALTenderSupply, 3); // TPC tender task
};


#endif

