#ifndef ALIPHOSALIENFILE_H
#define ALIPHOSALIENFILE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Description of logical filename in AliEn catalogue
// check here : /afs/cern.ch/user/p/peters/public/README.ALIEN 
//*--                  
//*-- Author: Yves Schutz (CERN)

// --- ROOT system ---
#include "TObject.h" 
#include "TString.h" 
#include "TGrid.h" 

// --- AliRoot header files ---

class AliPHOSAliEnFile : public TObject {

 public:
  
  AliPHOSAliEnFile() ; 
  AliPHOSAliEnFile(AliPHOSAliEnFile & lfn) : TObject(lfn) {
    lfn.Copy(*this) ; 
  } 
  virtual ~AliPHOSAliEnFile(void) ; 
  virtual void Copy(AliPHOSAliEnFile & lfn) ;
  
  void ListEvents() const ; 
  void ListRuns() const ; 
  TString GetRootDir() const { return fRoot ; }
  TString GetLFN() const ; 
  void Help() ; 
  
  Bool_t SetYearProd(TString year, TString prod) ; 
  Bool_t SetVers(TString vers) ; 
  Bool_t SetType(TString type) ; 

  Bool_t SetPath(TString year, TString prod, TString vers, TString type) ; 

  Bool_t SetRun(Int_t run) ; 
  Bool_t SetEvt(Int_t evt) ; 

  TString Pwd() const { return fPath ; }

  AliPHOSAliEnFile & operator = (const AliPHOSAliEnFile & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  
private:
  
  TGrid * fGrid ; //! connection to alien data catalogue 
  TString fRoot ; //! root directory
  TString fYear ; //! year of the DC 
  TString fProd ; //! production id 
  TString fVers ; //! aliroot tag version
  TString fType ; //! event type
  TString fRun  ; //! run number
  TString fEvt  ; //! event number
  TString fPath ; //! the lfn is fRoot/fYear/fProd/fVers/fType/fRun/fEvt
 

  ClassDef(AliPHOSAliEnFile,1)  
    
    };

#endif // AliPHOSALIENFILE_H
 
