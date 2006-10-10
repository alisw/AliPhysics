#ifndef ALIEMCALLOADER_H
#define ALIEMCALLOADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A singleton that returns various objects 
//  Should be used on the analysis stage to avoid confusing between different
//  branches of reconstruction tree: e.g. reading RecPoints and TS made from 
//  another set of RecPoints.
// 
//  The objects are retrived from folders.  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//    


// --- ROOT system ---
#include "TClonesArray.h"
#include "TFolder.h"  
#include "TTree.h"
class TString ;
class TParticle ;
class TTask ;

// --- AliRoot header files ---
#include "AliLoader.h"
#include "AliEMCALCalibData.h"

class AliLoader ;
class AliEMCAL ; 
class AliEMCALHit ;
class AliEMCALDigit ;
class AliEMCALSDigit ;
class AliEMCALRecPoint ; 

class AliEMCALLoader : public AliLoader {
  
 public:

  AliEMCALLoader();
  AliEMCALLoader(const AliEMCALLoader & obj);
  AliEMCALLoader(const Char_t *detname,const Char_t *eventfoldername); 
  
  virtual ~AliEMCALLoader() ; 

  // assignement operator requested by coding convention, but not needed
  const AliEMCALLoader & operator = (const AliEMCALLoader & ) {return *this;}

  virtual Int_t GetEvent();  // Overload to fill TClonesArray

  virtual void    CleanHits() const 
    { if (fHits) fHits->Clear(); AliLoader::CleanHits(); }
  virtual void    CleanSDigits() const
    { if (fSDigits) fSDigits->Clear(); AliLoader::CleanSDigits(); }
  virtual void    CleanDigits() const
    { if (fDigits) fDigits->Clear(); AliLoader::CleanDigits(); }
  virtual void    CleanRecPoints() const
    { if (fRecPoints) fRecPoints->Clear(); AliLoader::CleanRecPoints(); }

  // This does not work due to const
  /*
  virtual void   MakeHitsContainer() const { AliLoader::MakeHitsContainer(); TreeH()->Branch(fDetectorName,"TClonesArray",&fHits); }
  virtual void   MakeSDigitsContainer() const { AliLoader::MakeSDigitsContainer(); TreeS()->SetBranchAddress(fDetectorName,&fSDigits); }
  virtual void   MakeDigitsContainer() const { AliLoader::MakeDigitsContainer(); TreeD()->SetBranchAddress(fDetectorName,&fDigits); }
  virtual void   MakeRecPointsContainer() const { AliLoader::MakeRecPointsContainer(); TreeR()->SetBranchAddress(fgkECARecPointsBranchName,&fRecPoints); }
  */

  // ************    TClonesArrays Access functions

  TClonesArray*  Hits(void) { return fHits;}

  const AliEMCALHit*    Hit(Int_t index) {
    if (fHits)
      return (const AliEMCALHit*) fHits->At(index);
    return 0x0; 
  }

  TClonesArray*  SDigits()  { return fSDigits;}
  const AliEMCALDigit*  SDigit(Int_t index)  {
    if (fSDigits)
      return (const AliEMCALDigit*) fSDigits->At(index);
    return 0x0; 
  }

  TClonesArray*   Digits()  { return fDigits;}
  const AliEMCALDigit *  Digit(Int_t index)  {
    if (fDigits)
      return (const AliEMCALDigit*) fDigits->At(index);
    return 0x0; 
  }

  TObjArray * RecPoints()  { return fRecPoints;}
  const AliEMCALRecPoint * RecPoint(Int_t index)  {
    if (fRecPoints)
      return (const AliEMCALRecPoint*) fRecPoints->At(index);
    return 0x0; 
  }

  void   SetDebug(Int_t level) {fDebug = level;} // Set debug level

  //Calibration

  Int_t CalibrateRaw (Double_t energy, Int_t module, Int_t column, Int_t row);//take real calibration coefficients
  
  void  SetCalibData(AliEMCALCalibData* calibda)  { fgCalibData = calibda; }
  AliEMCALCalibData * CalibData(); // to get the calibration CDB object

private:
 
  static const TString fgkECARecPointsBranchName; //! Name of branch with ECA Reconstructed Points

  Int_t  fDebug ;             // Debug level

  // All data are stored in TTrees on file. 
  // These TCLonesArrays are temporary storage for reading or writing
  // (connected to TTrees with SetBranchAddress)
  TClonesArray     *fHits;         //! TClonesArray of hits (for tree reading)
  TClonesArray     *fDigits;       //! TClonesArray of digits (for tree reading)
  TClonesArray     *fSDigits;      //! TClonesArray of sdigits (for tree reading)
  TObjArray        *fRecPoints;    //! TClonesArray of recpoints (for tree reading)   
  
  static AliEMCALCalibData * fgCalibData;  //  calibration data 

  ClassDef(AliEMCALLoader,0)  // Algorithm class that provides methods to retrieve objects from a list knowing the index 
   
};

#endif // AliEMCALLOADER_H
