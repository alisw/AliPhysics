#ifndef ALIMUONTRIGGERLUT_H
#define ALIMUONTRIGGERLUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup calib
/// \class AliMUONTriggerLut
/// \brief MUON trigger look up table class 
///
//  Author: Philippe Crochet

#include <TNamed.h>

class TH3;
class TMap;

//----------------------------------------------
class AliMUONTriggerLut : public TNamed 
{
 public: 
  AliMUONTriggerLut();    // constructor
  virtual ~AliMUONTriggerLut();   // destructor

  Int_t Compare(const TObject* object) const;
  
  void GetLutOutput(Int_t circuit, Int_t xstrip, Int_t idev, Int_t ystrip, 
		    Int_t lutLpt[2], Int_t lutHpt[2]) const;

  void ReadFromFile(const char* filename);
  
  void SetContent(const char* hname, Int_t icirc, UChar_t istripX, 
                  UChar_t idev, Short_t value); 

  void SetLutCode(const UChar_t lutCode);

  void PrintLutCode();
  
 private:
  
    /// Not implemented copy constructor
  AliMUONTriggerLut (const AliMUONTriggerLut& AliMUONTriggerLut);
  /// Not implemented assignment operator
  AliMUONTriggerLut& operator=(const AliMUONTriggerLut& AliMUONTriggerLut); 

  void Add(TH3* h);

  Int_t Compare(TH3* h1, TH3* h2) const;
  
  Int_t GetMask(Int_t ystrip) const;

  void RegisterHistos();

private:
  TH3 *fLptPlus; ///< 3-d histogram with 234x32x31 bins Low pt Plus  
  TH3 *fLptMinu; ///< 3-d histogram with 234x32x31 bins Low pt Minus
  TH3 *fLptUnde; ///< 3-d histogram with 234x32x31 bins Low pt Undefined
  TH3 *fHptPlus; ///< 3-d histogram with 234x32x31 bins High pt Plus
  TH3 *fHptMinu; ///< 3-d histogram with 234x32x31 bins High pt Minus 
  TH3 *fHptUnde; ///< 3-d histogram with 234x32x31 bins High pt Undefined 
  TH3 *fAptPlus; ///< 3-d histogram with 234x32x31 bins All pt Plus  
  TH3 *fAptMinu; ///< 3-d histogram with 234x32x31 bins All pt Minus  
  TH3 *fAptUnde; ///< 3-d histogram with 234x32x31 bins All pt Undefined    

  TMap* fMap; //!<! from name to histo
  
  ClassDef(AliMUONTriggerLut,2) // Trigger Look up Table class

};
#endif






