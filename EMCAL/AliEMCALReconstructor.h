#ifndef ALIEMCALRECONSTRUCTOR_H
#define ALIEMCALRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Wrapping class for reconstruction
//*--
//*-- Author: Yves Schutz (SUBATECH) 
//*--         Dmitri Peressounko (SUBATECH & Kurchatov Institute)
// Reconstruction class. Redesigned from the old AliReconstructionner class and 
// derived from STEER/AliReconstructor. 
// 

// --- ROOT system ---
class TClonesArray;
class TTree;

// --- AliRoot header files ---
#include "AliReconstructor.h" 
#include "AliEMCALTracker.h" 
#include "AliEMCALRecParam.h"

class AliEMCALDigitizer ;
class AliEMCALClusterizer ;
class AliEMCALSDigitizer ;
class AliEMCALRecParam;
class AliESDEvent ;
class AliRawReader ;
class AliEMCALRawUtils;
class AliEMCALGeometry;
class AliEMCALCalibData ;
class AliCaloCalibPedestal ;
class AliEMCALTriggerElectronics;
class AliEMCALTriggerData;

class AliEMCALReconstructor : public AliReconstructor {
  
 public:
  
  AliEMCALReconstructor() ; //ctor            
  
  virtual ~AliEMCALReconstructor() ; //dtor
  
  virtual  void  Init() {;}
  virtual  void  InitClusterizer() const;
  
  using AliReconstructor::FillESD;
  virtual void   FillESD(TTree* digitsTree, TTree* clustersTree, AliESDEvent* esd) const;
  AliTracker*    CreateTracker () const {return new AliEMCALTracker;} 
  using AliReconstructor::Reconstruct;
  virtual void   Reconstruct(TTree* digitsTree, TTree* clustersTree) const ;
  
  virtual Bool_t HasDigitConversion() const {return kTRUE;};
  virtual void   ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;
  
  static void    SetRecParam(AliEMCALRecParam * recParam){ fgkRecParam = recParam;}
  
  void           ReadDigitsArrayFromTree(TTree *digitsTree) const;
  
  static const AliEMCALRecParam* GetRecParam() { 
    return dynamic_cast<const AliEMCALRecParam*>(AliReconstructor::GetRecoParam(6)); }
  
  static TClonesArray* GetDigitsArr() {return fgDigitsArr;}
  
  void           FillMisalMatrixes(AliESDEvent* esd)const ;

  //New class used to sort the matched tracks
  class  AliEMCALMatch : public TObject
  {
  public:
    AliEMCALMatch();
    AliEMCALMatch(const AliEMCALMatch& copy);
    virtual ~AliEMCALMatch() { }
    //----------------------------------------------------------------------------
    Int_t     Compare(const TObject *obj) const;
    Bool_t    IsSortable() const {return kTRUE;}
    Double_t  GetDistance() const {return fDistance;}
    Double_t  GetdEta() const {return fdEta;}
    Double_t  GetdPhi() const {return fdPhi;}
    Int_t     GetIndexT() const {return fIndexT;}
    void      SetIndexT(Int_t itr) {fIndexT=itr;}
    void      SetDistance(Double_t dist) {fDistance=dist;}
    void      SetdEta(Double_t dEta) {fdEta=dEta;}
    void      SetdPhi(Double_t dPhi) {fdPhi=dPhi;}
  private:
    Int_t      fIndexT;      // track index in 'fTracks' array
    Double_t   fDistance;    // track - cluster distance
    Double_t   fdEta;        // track - cluster residual in eta
    Double_t   fdPhi;        // track - cluster residual in phi
  };
  Bool_t CalculateResidual(AliESDtrack *track, AliESDCaloCluster *cluster, Float_t &dEta, Float_t &dPhi) const;
  
 private:
  
  AliEMCALReconstructor(const AliEMCALReconstructor &); //Not implemented
  AliEMCALReconstructor & operator = (const AliEMCALReconstructor &); //Not implemented
  
  AliEMCALGeometry           * fGeom;             // pointer to the EMCAL geometry
  static AliEMCALClusterizer * fgClusterizer;     // clusterizer
  static AliEMCALRawUtils    * fgRawUtils;        // raw utilities class 
  
  //Temporal arrays with clusters, digits, triggers, to be reused per event
  static TClonesArray        * fgDigitsArr;       //-> Array with EMCAL digits
  static TObjArray           * fgClustersArr;     //-> Array with EMCAL clusters
  static TClonesArray        * fgTriggerDigits;   //-> Array with EMCAL trigger digits
  
  //OCDB
  static const AliEMCALRecParam* fgkRecParam;     // reconstruction parameters for EMCAL
  AliEMCALCalibData          * fCalibData   ;     //! Calibration database if aval
  AliCaloCalibPedestal       * fPedestalData ;    //! Tower status database if aval
  
  //Trigger specific
  static AliEMCALTriggerElectronics* fgTriggerProcessor; // Trigger preprocessor  
  AliEMCALTriggerData        * fTriggerData;      // Trigger parameters data container

  //Track matching
  TList                      * fMatches;          //! collection of matches between tracks and clusters
  
  ClassDef(AliEMCALReconstructor,12)  // Reconstruction algorithm class (Base Class)
}; 

#endif // ALIEMCALRECONSTRUCTOR_H

