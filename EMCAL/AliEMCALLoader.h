#ifndef ALIEMCALLOADER_H
#define ALIEMCALLOADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  The AliEMCALLoader gets the TClonesArray and TObjArray for reading
//  Hits, Dgits, SDigits and RecPoints. Filling is managed in the GetEvent()
//  method. The objects are retrived from  the corresponding folders.  
//
//  It also provides acces methods to the calibration and simulation OCDB parameters 
//

//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//    


// --- ROOT system ---
#include "TClonesArray.h"
#include "TFolder.h"  
#include "TTree.h"
class TString ;
class TParticle ;

// --- AliRoot header files ---
#include "AliLoader.h"
#include "AliEMCALCalibData.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALSimParam.h"
#include "AliEMCALRecParam.h"

class AliLoader ;
class AliEMCAL ; 
class AliEMCALDigit ;
class AliEMCALSDigit ;
class AliEMCALRecPoint ; 

class AliEMCALLoader : public AliLoader {
  
 public:

  AliEMCALLoader();
  AliEMCALLoader(const Char_t *detname,const Char_t *eventfoldername); 
  AliEMCALLoader(const Char_t *name,TFolder *topfolder);
  
  virtual ~AliEMCALLoader() ; 

  virtual Int_t GetEvent();  // Overload to fill TClonesArray

  //Clean arrays methods
  virtual void    CleanHits() const {GetHitsDataLoader()->Clean();}       
  virtual void    CleanSDigits() const {GetSDigitsDataLoader()->Clean();}  
  virtual void    CleanDigits() const {GetDigitsDataLoader()->Clean();}  
  virtual void    CleanRecPoints() const {GetRecPointsDataLoader()->Clean();}  
  
  // Initialize arrays methods
  void MakeSDigitsArray() ;
  void MakeDigitsArray() ;
  void MakeRecPointsArray() ;
  
  // ************    TClonesArrays Access functions
  
  TClonesArray*  SDigits() {return (TClonesArray*)GetDetectorData(fgkECASDigitsBranchName);} //const { return fSDigits;}
  const AliEMCALDigit*  SDigit(Int_t index)  {
    if (SDigits())return (const AliEMCALDigit*) SDigits()->At(index);
    return 0x0; 
  }
  
  TClonesArray*   Digits() {return (TClonesArray*)GetDetectorData(fgkECADigitsBranchName);}//const { return fDigits;}
  const AliEMCALDigit *  Digit(Int_t index)  {
    if (Digits()) return (const AliEMCALDigit*) Digits()->At(index);
    return 0x0; 
  }
  
  TObjArray * RecPoints()  {return (TObjArray*)GetDetectorData(fgkECARecPointsBranchName);}//const { return fRecPoints;}
  const AliEMCALRecPoint * RecPoint(Int_t index)  {
    if (RecPoints())return (const AliEMCALRecPoint*) RecPoints()->At(index);
    return 0x0; 
  }
  
  void   SetDebug(Int_t level) {fDebug = level;} // Set debug level
  
  //OCDB access methods
  
  void  SetCalibData(AliEMCALCalibData* calibda)  { fgCalibData = calibda; }
  AliEMCALCalibData * CalibData();              // to get the calibration CDB object
  
  void  SetPedestalData(AliCaloCalibPedestal* caloped)  { fgCaloPed = caloped; }
  AliCaloCalibPedestal* PedestalData();         // to get the pedestal CDB object
  
  void  SetSimParam(AliEMCALSimParam* simparam)  { fgSimParam = simparam; }
  AliEMCALSimParam* SimulationParameters();     // to get the simulation parameter CDB object
  
  void  SetRecParam(AliEMCALRecParam* recparam)  { fgRecParam = recparam; }
  AliEMCALRecParam* ReconstructionParameters(Int_t eventType); // to get the reconstruction parameter CDB object

 private:
  
  // assignement operator requested by coding convention, but not needed
  AliEMCALLoader(const AliEMCALLoader &); //Not implemented
  const AliEMCALLoader & operator = (const AliEMCALLoader &); //Not implemented
  
  static const TString fgkECASDigitsBranchName;   //! Name of branch with ECA SDigits
  static const TString fgkECADigitsBranchName;    //! Name of branch with ECA Digits
  static const TString fgkECARecPointsBranchName; //! Name of branch with ECA Reconstructed Points
  
  Int_t  fDebug ;             // Debug level
	
  static AliEMCALCalibData    * fgCalibData;  //  calibration data 
  static AliCaloCalibPedestal * fgCaloPed;    //  dead map
  static AliEMCALSimParam     * fgSimParam;   //  sim param 
  static AliEMCALRecParam     * fgRecParam;   //  rec param 

  ClassDef(AliEMCALLoader,7)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 
    
};

#endif // AliEMCALLOADER_H
