#ifndef ALIPHOSRAW2DIGITS_H
#define ALIPHOSRAW2DIGITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
//  Base Class for PHOS     
//                  
/*-- Author: Maxim Volkov (RRC KI)
              Dmitri Peressounko (RRC KI & SUBATECH)
              Yuri Kharlov (IHEP & SUBATECH)     */

// --- ROOT system ---
#include "TTask.h"
class TClonesArray ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSGeometry ;
class AliPHOSBeamTestEvent ;
class AliPHOSConTableDB ;

class AliPHOSRaw2Digits : public TTask {

public:
  AliPHOSRaw2Digits() ;          // ctor
  AliPHOSRaw2Digits(const char * inputFileName) ;         
  AliPHOSRaw2Digits(AliPHOSRaw2Digits & r2d) ;          // cpy ctor
  virtual ~AliPHOSRaw2Digits() ; // dtor

  void Exec(const Option_t * = "") ;

  void SetBeamEnergy(Float_t energy){fBeamEnergy = energy ;}
  void SetInputFile(TString inname="Run_1234.fz"){fInName=inname ; }
  void SetDebugLevel(Int_t idebug=1){fDebug=idebug ;}

  //Set position of the target in the given run.
  //The reference system is following
  //Z axis along beam direction, from target to prototype (0-surface of prototype)
  //X axis along columns of prototype (0-center of prototype)
  //Y axis along raws of prototype    (0-center of prototype)
  void SetTargetPosition(Double_t * pos)
    {for(Int_t i=0;i<3;i++)fTarget[i]=pos[i] ;}
  void SetConTableDB(AliPHOSConTableDB * ctdb){fctdb = ctdb ;}
  void SetMaxEventsPerFile(Int_t nev=20000){fMaxPerFile = nev ;}
  void Print(const Option_t * = "")const ;
  AliPHOSRaw2Digits & operator = ( AliPHOSRaw2Digits & /*r2d*/ ) { return *this ; } 
  
private:
  Bool_t StartRootFiles(void) const ;
  Bool_t CloseRootFiles(void) ;
  Bool_t ProcessRawFile() ;
  void Swab4(void *from, void *to, size_t nwords) const ;
  void Swab2(void *from, void *to, size_t nwords) const ;
  Bool_t Init() ;
  void WriteDigits(void) ;

  TClonesArray * fDigits ;             //!list of final digits
  AliPHOSBeamTestEvent * fPHOSHeader ; //!PHOSBeamTest header 
  AliPHOSConTableDB * fctdb ;          //!
  Double_t fTarget[3] ;                //!Position of the target
  TFile * fHeaderFile ;                //!galice.root file
  TFile * fDigitsFile ;                //!file with digits
  Float_t fBeamEnergy ;    //BeamEnergy 
  Int_t   fMaxPerFile ;    //!Maximal number  of events per root file
  Int_t   fEvent ;         //Event number
  Int_t   fStatus ;        //status of input file: OK, not found etc.
  TString fInName ;        // FileName of the input file
  Bool_t  fDebug ;         //!
  Bool_t  fIsInitialized ; //!
 
  UInt_t  fMK1 ;     //!ZEBRA markers
  UInt_t  fMK2 ;     //!ZEBRA markers
  UInt_t  fMK3 ;     //!ZEBRA markers
  UInt_t  fMK4 ;     //!ZEBRA markers
  UInt_t  fCKW ;     //!ZEBRA markers

  ClassDef(AliPHOSRaw2Digits,1)  // description 

};

#endif // AliPHOSRAW2DIGITS_H
