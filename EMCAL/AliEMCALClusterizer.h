#ifndef ALIEMCALCLUSTERIZER_H
#define ALIEMCALCLUSTERIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
                            
/* $Id$ */

//_________________________________________________________________________
//  Base class for the clusterization algorithm (pure abstract)
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
// --- ROOT system ---
#include "AliLog.h"
#include "TObject.h" 
class TTree;

// --- Standard library ---

// --- AliRoot header files ---
class AliEMCALGeometry ;
class AliEMCALCalibData ;
class AliCaloCalibPedestal ;

class AliEMCALClusterizer : public TObject {

public:

  AliEMCALClusterizer() ;        // default ctor
  virtual ~AliEMCALClusterizer() ; // dtorEM
  AliEMCALClusterizer(AliEMCALGeometry* geometry);
  AliEMCALClusterizer(AliEMCALGeometry* geometry, AliEMCALCalibData * calib, AliCaloCalibPedestal * pedestal);

  virtual void    Digits2Clusters(Option_t *option) = 0;

  virtual Float_t Calibrate(const Float_t amp, const Float_t time, const Int_t cellId) ;  // Tranforms Amp to energy 
  virtual void    Init() ;
  virtual void    InitParameters() ; //{ AliInfo("Overload this method."); }

  //Get/Set reconstruction parameters
  virtual void  GetCalibrationParameters(void) ;
  virtual void  GetCaloCalibPedestal(void) ;
  virtual void  SetCalibrationParameters(AliEMCALCalibData * calib)   { fCalibData = calib ; }
  virtual void  SetCaloCalibPedestal(AliCaloCalibPedestal  * caloped) { fCaloPed   = caloped ; }
  
  virtual Float_t GetTimeMin()           const { return fTimeMin ; }
  virtual Float_t GetTimeMax()           const { return fTimeMax ; }
  virtual Float_t GetTimeCut()           const { return fTimeCut ; }
  virtual void    GetNumberOfClustersFound(int numb )const { numb = fNumberOfECAClusters ;} 
  virtual Float_t GetECAClusteringThreshold()        const { return fECAClusteringThreshold;}  
  virtual Float_t GetECALocalMaxCut()                const { return fECALocMaxCut;} 
  virtual Float_t GetECALogWeight()                  const { return fECAW0;}
  virtual Float_t GetMinECut()                       const { return fMinECut;}

  virtual void SetTimeMin(Float_t t)					        { fTimeMin = t ;}
  virtual void SetTimeMax(Float_t t)					        { fTimeMax = t ;}
  virtual void SetTimeCut(Float_t t)					        { fTimeCut = t ;}
  virtual void SetECAClusteringThreshold(Float_t th)  { fECAClusteringThreshold = th ; }
  virtual void SetMinECut(Float_t mine)               { fMinECut = mine; }
  virtual void SetECALocalMaxCut(Float_t cut)         { fECALocMaxCut = cut ; }
  virtual void SetECALogWeight(Float_t w)             { fECAW0 = w ; }
  virtual void SetUnfolding(Bool_t toUnfold = kTRUE ) {fToUnfold = toUnfold ;}  

  virtual void SetInput(TTree *digitsTree);
  virtual void SetOutput(TTree *clustersTree);
  
  virtual void Print(Option_t * option)const ;
  virtual void PrintRecPoints(Option_t * option);
  virtual void PrintRecoInfo();                        //*MENU*
  
  virtual const char * Version() const {Warning("Version", "Not Defined") ; return 0 ; } 

protected:

  virtual void MakeClusters() = 0;

  TClonesArray *fDigitsArr; // Array with EMCAL digits
  TTree *fTreeR;            // Tree with output clusters
  TObjArray    *fRecPoints; // Array with EMCAL clusters
  
  AliEMCALGeometry     * fGeom;       //! pointer to geometry for utilities
  AliEMCALCalibData    * fCalibData ; //! Calibration database if aval
  AliCaloCalibPedestal * fCaloPed   ; //! Tower status map if aval
  
  Float_t fADCchannelECA ;           // width of one ADC channel for EC section (GeV)
  Float_t fADCpedestalECA ;          // pedestal of ADC for EC section (GeV) 

  Float_t fTimeMin ;                 // Minimum time of physical signal in a cell/digit
  Float_t fTimeMax ;                 // Maximum time of physical signal in a cell/digit
  Float_t fTimeCut ;                 // Maximum time difference between the digits inside EMC cluster

  Bool_t  fDefaultInit;              //! Says if the task was created by defaut ctor (only parameters are initialized)
  Bool_t  fToUnfold ;                // To perform unfolding 
  Int_t   fNumberOfECAClusters ;     // number of clusters found in EC section
  
  Float_t fECAClusteringThreshold ;  // minimum energy to seed a EC digit in a cluster
  Float_t fECALocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fECAW0 ;                   // logarithmic weight for the cluster center of gravity calculation
  Float_t fMinECut;                  // Minimum energy for a digit to be a member of a cluster
  
private:
  AliEMCALClusterizer(const AliEMCALClusterizer &); //copy ctor
  AliEMCALClusterizer & operator = (const AliEMCALClusterizer &);
  
  
  
  ClassDef(AliEMCALClusterizer,2)  // Clusterization algorithm class 
} ;

#endif // AliEMCALCLUSTERIZER_H
