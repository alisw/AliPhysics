#ifndef ALIPHOSCLUSTERIZERV1_H
#define ALIPHOSCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log: AliPHOSClusterizerv1.h,v $
 * Revision 1.54  2007/10/01 20:24:08  kharlov
 * Memory leaks fixed
 *
 * Revision 1.53  2007/08/28 12:55:07  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.52  2007/08/07 14:16:00  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.51  2007/04/11 11:55:45  policheh
 * SetDistancesToBadChannels() added.
 *
 * Revision 1.50  2007/03/28 19:18:15  kharlov
 * RecPoints recalculation in TSM removed
 *
 * Revision 1.49  2007/03/06 06:51:27  kharlov
 * Calculation of cluster properties dep. on vertex posponed to TrackSegmentMaker
 *
 * Revision 1.48  2006/08/30 16:12:52  kharlov
 * Reconstruction of raw data from beam test 2006 (B.Polichtchouk)
 *
 * Revision 1.47  2006/08/25 16:56:30  kharlov
 * Compliance with Effective C++
 *
 * Revision 1.46  2006/08/01 12:20:17  cvetan
 * 1. Adding a possibility to read and reconstruct an old rcu formatted raw data. This is controlled by an option of AliReconstruction and AliPHOSReconstructor. 2. In case of raw data processing (without galice.root) create the default AliPHOSGeometry object. Most likely this should be moved to the CDB
 *
 * Revision 1.45  2006/04/29 20:26:46  hristov
 * Separate EMC and CPV calibration (Yu.Kharlov)
 *
 * Revision 1.44  2005/09/02 14:32:07  kharlov
 * Calibration of raw data
 *
 * Revision 1.43  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  Implementation version 1 of the clusterization algorithm                     
//  Performs clusterization (collects neighbouring active cells) and 
//  unfolding of the clusters with several local maxima.  
//  results are stored in TreeR#, branches PHOSEmcRP (EMC recPoints),
//  PHOSCpvRP (CPV RecPoints) and AliPHOSClusterizer
//
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
class TClonesArray ;
class TVector3 ;
// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSClusterizer.h"
class AliPHOSEmcRecPoint ; 
class AliPHOSDigit ;
class AliPHOSDigitizer ;
class AliPHOSGeometry ;

class AliPHOSClusterizerv1 : public AliPHOSClusterizer {
  
public:
  
  AliPHOSClusterizerv1() ;
  AliPHOSClusterizerv1(AliPHOSGeometry *geom);
  virtual ~AliPHOSClusterizerv1()  ;

  void    InitParameters() ;
  
  virtual Int_t   AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)const ; 
                               // Checks if digits are in neighbour cells 

  virtual void    GetNumberOfClustersFound(int * numb )const{  numb[0] = fNumberOfEmcClusters ; 
                                                               numb[1] = fNumberOfCpvClusters ; }

  virtual Float_t GetEmcClusteringThreshold()const{ return fEmcClusteringThreshold;}
  virtual Float_t GetEmcLocalMaxCut()const        { return fEmcLocMaxCut;} 
  virtual Float_t GetEmcLogWeight()const          { return fW0;}  
  virtual Float_t GetEmcTimeGate() const          { return fEmcTimeGate ; }
  virtual Float_t GetCpvClusteringThreshold()const{ return fCpvClusteringThreshold;  } 
  virtual Float_t GetCpvLocalMaxCut()const        { return fCpvLocMaxCut;} 
  virtual Float_t GetCpvLogWeight()const          { return fW0CPV;}  
  virtual Float_t GetEcoreRadius()const           { return fEcoreRadius;}  
  //  virtual const char *  GetRecPointsBranch() const{ return GetName() ;}

  virtual void    Digits2Clusters(Option_t *option);

  void Print(const Option_t * = "")const ;

  virtual void SetEmcClusteringThreshold(Float_t cluth)  { fEmcClusteringThreshold = cluth ; }
  virtual void SetEmcLocalMaxCut(Float_t cut)            { fEmcLocMaxCut = cut ; }
  virtual void SetEmcLogWeight(Float_t w)                { fW0 = w ; }
  virtual void SetEmcTimeGate(Float_t gate)              { fEmcTimeGate = gate ;}
  virtual void SetCpvClusteringThreshold(Float_t cluth)  { fCpvClusteringThreshold = cluth ; }
  virtual void SetCpvLocalMaxCut(Float_t cut)            { fCpvLocMaxCut = cut ; }
  virtual void SetCpvLogWeight(Float_t w)                { fW0CPV = w ; }
  virtual void SetUnfolding(Bool_t toUnfold = kTRUE )    { fToUnfold = toUnfold ;}
  virtual void SetCoreRadius(Float_t coreRadius)         { fEcoreRadius = coreRadius ;}
  //Switch to "on flyght" mode, without writing to TreeR and file  
  void SetWriting(Bool_t toWrite = kFALSE){fWrite = toWrite;} 
  static Double_t ShowerShape(Double_t x, Double_t z) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
  static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passed to MINUIT
  //  void Unload() ; 
  virtual const char * Version() const { return "clu-v1"; }  

protected:

  void           WriteRecPoints() ;
  virtual void   MakeClusters( ) ;            
  virtual Bool_t IsInEmc (AliPHOSDigit * digit)const ;     // Tells if id digit is in EMC
  virtual Bool_t IsInCpv (AliPHOSDigit * digit)const ;     // Tells if id digit is in CPV
  void           CleanDigits(TClonesArray * digits) ;
  void           SetDistancesToBadChannels();
  virtual Float_t Calibrate(Float_t amp, Int_t absId) const ;  // Tranforms ADC counts to energy   
  virtual Float_t CalibrateT(Float_t amp, Int_t absId) const ;  //Tranforms Sample counts to sec.
   
private:
  AliPHOSClusterizerv1(const AliPHOSClusterizerv1 & clu) ;
  AliPHOSClusterizerv1 & operator = (const AliPHOSClusterizerv1 & obj);

  Bool_t  FindFit(AliPHOSEmcRecPoint * emcRP, AliPHOSDigit ** MaxAt, Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit
  void    Init() ;

  virtual void   MakeUnfolding() ;
  void           UnfoldCluster(AliPHOSEmcRecPoint * iniEmc,Int_t Nmax, 
		       AliPHOSDigit ** maxAt,Float_t * maxAtEnergy ) ; //Unfolds cluster using TMinuit package
  void           PrintRecPoints(Option_t * option) ;

private:

  Bool_t  fDefaultInit;              //! Says if the task was created by defaut ctor (only parameters are initialized)
  Int_t   fEmcCrystals ;             // number of EMC cristals in PHOS

  Bool_t  fToUnfold ;                // To perform unfolding 
  Bool_t  fWrite ;                   // Write RecPoints to TreeR  
  Bool_t  fDigitsUsed[53760];        //Mark digits as already used in cluster (EMC:5*56*64 ; CPV: 5*56*128)
 
  Int_t   fNumberOfEmcClusters ;     // number of EMC clusters found
  Int_t   fNumberOfCpvClusters ;     // number of CPV clusters found
 
  Float_t fEmcClusteringThreshold ;  // minimum energy to start EMC cluster
  Float_t fCpvClusteringThreshold ;  // minimum energy to start CPV cluster
  Float_t fEmcLocMaxCut ;            // minimum energy difference to distinguish local maxima in a cluster
  Float_t fW0 ;                      // logarithmic weight for the cluster center of gravity calculation
  Float_t fCpvLocMaxCut ;            // minimum energy difference to distinguish local maxima in a CPV cluster
  Float_t fW0CPV ;                   // logarithmic weight for the CPV cluster center of gravity calculation
  Float_t fEmcTimeGate ;             // Maximum time difference between the digits in ont EMC cluster
  Float_t fEcoreRadius ;             // Radius within which the core energy is calculated, in cm

  ClassDef(AliPHOSClusterizerv1,7)   // Clusterizer implementation version 1

};

#endif // AliPHOSCLUSTERIZERV1_H
