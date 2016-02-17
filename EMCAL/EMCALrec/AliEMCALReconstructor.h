#ifndef ALIEMCALRECONSTRUCTOR_H
#define ALIEMCALRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALReconstructor
/// \brief Wrapping class for EMCal reconstruction
///
/// Wrapping class for reconstruction.
///
/// It steers
///   * Raw data fitting to create digits
///   * Clusterization from digits
///   * ESD filling from clusters
///
/// Trigger data handling is also dealt here.
///
/// \author Yves Schutz (SUBATECH), originally 
/// \author Dmitri Peressounko (SUBATECH & Kurchatov Institute), PHOS copy
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-Grenoble
/// \author plus others that have put their hands, list them ...
/// 

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
class AliEMCALCalibTime ;
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

  //
  // New class used to sort the matched tracks
  //
  class  AliEMCALMatch : public TObject
  {
  public:
    AliEMCALMatch();
    AliEMCALMatch(const AliEMCALMatch& copy);
    AliEMCALMatch& operator = (const AliEMCALMatch& source) ;
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
    
    Int_t      fIndexT;      ///< Track index in 'fTracks' array
    Double_t   fDistance;    ///< Track - cluster distance
    Double_t   fdEta;        ///< Track - cluster residual in eta
    Double_t   fdPhi;        ///< Track - cluster residual in phi
  };
  
  Bool_t CalculateResidual(AliESDtrack *track, AliESDCaloCluster *cluster, Float_t &dEta, Float_t &dPhi) const;
  
 private:
  
  AliEMCALReconstructor              (const AliEMCALReconstructor &); /// Not implemented
  AliEMCALReconstructor & operator = (const AliEMCALReconstructor &); /// Not implemented
  
  AliEMCALGeometry           * fGeom;             ///< Access to the geometry
  static AliEMCALClusterizer * fgClusterizer;     ///< Access to the clusterization tools
  static AliEMCALRawUtils    * fgRawUtils;        ///< Access to raw fitting tools 
    
  /// Temporary array with digits, to be reused in the event.
  static TClonesArray        * fgDigitsArr;       //-> 
  
  /// Temporary array with clusters, to be reused in the event.
  static TObjArray           * fgClustersArr;     //->
  
  /// Temporary array with trigger digits, to be reused in the event.
  static TClonesArray        * fgTriggerDigits;   //->
  
  // OCDB
  static const AliEMCALRecParam* fgkRecParam;     ///<  Access to OCDB reconstruction parameters
  AliEMCALCalibData          * fCalibData   ;     //!<! Access to OCDB energy calibration database if available
  AliEMCALCalibTime          * fCalibTime   ;     //!<! Access to OCDB time calibration database if available
  AliCaloCalibPedestal       * fPedestalData;     //!<! Access to OCDB tower status database if available
  
  //Trigger specific
  static AliEMCALTriggerElectronics *fgTriggerProcessor; ///< Trigger preprocessor  
  static TClonesArray               *fgTriggerData;      ///< Trigger parameters data container

  //Track matching
  TList                      * fMatches;          //!<! Collection of matches between tracks and clusters
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALReconstructor,14)  // Reconstruction algorithm class (Base Class)
  /// \endcond

}; 

#endif // ALIEMCALRECONSTRUCTOR_H

