#ifndef ALIRICHDIGITIZER_H
#define ALIRICHDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliDigitizer.h"

class AliRunDigitizer;
class AliRICHSDigit;
class AliHitMap;
class AliRICHhit;

class AliRICHDigitizer : public AliDigitizer 
{
public:
    enum {kBgTag = -1};
    
           AliRICHDigitizer();
           AliRICHDigitizer(AliRunDigitizer * manager);
  virtual ~AliRICHDigitizer();
        
  Bool_t Exists(const AliRICHSDigit *p);          //Compare pad hits
  void   Update(AliRICHSDigit *sdigit);           //Update a pad hit
  void   CreateNew(AliRICHSDigit *sdigit);        //Create a new hit
  Bool_t Init();                                  //virtual  Initialize merging and digitisation
  void   Exec(Option_t* option=0);                //virtual Do the main work
  Int_t  GetDebug() {return fDebug;}              //get debug level
  void   SetDebug(Int_t level){fDebug = level;}   //set debug level    
  
  AliRICHSDigit*   FirstPad(AliRICHhit *hit, TClonesArray *clusters);
  AliRICHSDigit*   NextPad(TClonesArray *clusters);
protected:
  TClonesArray *fHits;            //! List of hits for one track only
  TClonesArray *fSDigits;         //! List of clusters for one track only
  AliHitMap **fHitMap;            //! pointer to array of pointers to hitmaps
  Int_t fNch;                     //! chamber nr (loop variable)
  Int_t fTrack;                   //! track nr (loop variable)
  TObjArray *fTDList;             //! list of AliRICHTransientDigit
  Int_t fCounter;                 //! nr. of AliRICHTransientDigit
  Bool_t fSignal;                 //! kTRUE if signal file is processed
  Int_t fMask;                    //! mask dependent on input file
  Int_t fDigits[5];               //! array with digits
  Int_t fDebug;                   //! debug level
  
  ClassDef(AliRICHDigitizer,1)
};    
#endif
