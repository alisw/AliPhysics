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

// --- ROOT system ---
class TList;
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

// --- Standard library ---



class AliEMCALReconstructor : public AliReconstructor {

public:

  AliEMCALReconstructor() ; //ctor            
  AliEMCALReconstructor(const AliEMCALReconstructor & rec);
   
  virtual ~AliEMCALReconstructor() ; //dtor

  virtual  void Init();
  Bool_t       Debug() const { return fDebug ; }

  using AliReconstructor::FillESD;
  virtual void FillESD(TTree* digitsTree, TTree* clustersTree, 
		       AliESDEvent* esd) const;
  AliTracker*  CreateTracker () const 
  {return new AliEMCALTracker;} 
  using AliReconstructor::Reconstruct;
  virtual void Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual Bool_t             HasDigitConversion() const {return kTRUE;};
  virtual void               ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;
  
  
  AliEMCALReconstructor & operator = (const AliEMCALReconstructor & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  
  static void   SetRecParam(AliEMCALRecParam * recParam){ fgkRecParam = recParam;}

  void   ReadDigitsArrayFromTree(TTree *digitsTree) const;

  TList *GetList() {return fList;}

  static const AliEMCALRecParam* GetRecParam() { 
    return dynamic_cast<const AliEMCALRecParam*>(AliReconstructor::GetRecoParam(6)); }

  static TClonesArray* GetDigitsArr() {return fgDigitsArr;}

  void FillMisalMatrixes(AliESDEvent* esd)const ;

private:
  
  Bool_t fDebug; //! verbosity controller

  TList *fList;  //! List of hists (only for trigger now)
  AliEMCALGeometry         *fGeom;           // pointer to the EMCAL geometry

  static AliEMCALClusterizer* fgClusterizer; // clusterizer
  static const AliEMCALRecParam*   fgkRecParam; // reconstruction
						// parameters for EMCAL
  static AliEMCALRawUtils*   fgRawUtils;  // raw utilities class -
					  // only need one per reco
  static TClonesArray*       fgDigitsArr; // Array with EMCAL digits
  AliEMCALCalibData * fCalibData  ;   //! Calibration database if aval

  ClassDef(AliEMCALReconstructor,7)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIEMCALRECONSTRUCTOR_H
