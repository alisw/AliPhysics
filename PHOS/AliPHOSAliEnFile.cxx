/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */

//_________________________________________________________________________
// To navigate in the AliEn catalogue (very elementary)
// check here : /afs/cern.ch/user/p/peters/public/README.ALIEN                   
//-- Author: Yves Schutz (CERN)

// --- ROOT system ---
#include "TObjString.h"   
#include "TAlienResult.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSAliEnFile.h" 

ClassImp(AliPHOSAliEnFile) ;

//____________________________________________________________________________
AliPHOSAliEnFile::AliPHOSAliEnFile()
{
  // default ctor; Doing initialisation ; 
  
  fGrid = TGrid::Connect("alien://aliendb1.cern.ch:15000/?direct") ; 
  if ( !fGrid ) 
    Fatal("ctor", "Cannot connect to alien://aliendb1.cern.ch:15000/?direct") ; 
    
  fRoot = "/alice/production/aliprod" ; 
  if ( !fGrid->OpenDir(fRoot) ) 
    Fatal("ctor", "Cannot find directory %s ", fRoot.Data() ) ; 
  
  fYear = "" ; 
  fProd = "" ; 
  fVers = "" ; 
  fType = "" ; 
  fRun  = "" ; 
  fEvt  = "" ; 
  
  fPath += fRoot ; 

}

//____________________________________________________________________________
AliPHOSAliEnFile::~AliPHOSAliEnFile()
{
}

//____________________________________________________________________________
TString AliPHOSAliEnFile::GetLFN() const
{
  TString fileName(Pwd()) ;
  fileName += "galice.root" ; 
  if ( !fGrid->GetAccessPath(fileName) ) { 
    Warning("GetLFN", "file %s does not exist", fileName.Data()) ; 
    fileName = "" ; 
  }
  else 
    fileName.Prepend("alien://") ; 
  return fileName ; 
}

//____________________________________________________________________________
void AliPHOSAliEnFile::Copy(AliPHOSAliEnFile & lfn) 
{
  //Copy method used by the Copy ctor 
  fRoot = lfn.fRoot ; 
  fYear = lfn.fYear ; 
  fProd = lfn.fProd ; 
  fVers = lfn.fVers ; 
  fType = lfn.fType ; 
  fRun  = lfn.fRun ; 
  fEvt  = lfn.fEvt ; 
  TObject::Copy(lfn) ; 
}

//____________________________________________________________________________
void AliPHOSAliEnFile::Help()
{
  // Prints information on available lfn's
  
  Info("Info", "  ") ; 

}

//____________________________________________________________________________
void AliPHOSAliEnFile::ListEvents() const
{
  // list the available events for the current path and run selected

  char path[80] ; 
  sprintf(path, "%s/%s-%s/%s/%s/%s", fRoot.Data(), fYear.Data(), fProd.Data(), fVers.Data(), fType.Data(), fRun.Data()) ; 
  Info("ListEvents", "Searching %s", path) ; 
  Grid_ResultHandle_t gr = fGrid->Find(path, "galice.root") ; 
  TAlienResult ar(gr) ;
  ar.Print() ; 
}

//____________________________________________________________________________
void AliPHOSAliEnFile::ListRuns() const
{
  // list the available runs for the current path selected

  char path[80] ; 
  sprintf(path, "%s/%s-%s/%s/%s", fRoot.Data(), fYear.Data(), fProd.Data(), fVers.Data(), fType.Data()) ; 
  Info("ListEvents", "Searching %s", path) ; 
  Grid_ResultHandle_t gr = fGrid->OpenDir(path) ; 
  TAlienResult ar(gr) ;
  ar.Print() ; 
}

//____________________________________________________________________________
Bool_t AliPHOSAliEnFile::SetYearProd(TString year, TString prod)
{
  // set the year and verifies if the directory exists
  Bool_t rv = kFALSE ;
  char tempo[80] ; 
  sprintf(tempo, "/%s-%s", year.Data(), prod.Data()) ; 

  TString path(fRoot) ; 
  path += tempo ; 
  if ( !fGrid->OpenDir(path) )  
    Error("ctor", "Cannot find directory %s ", path.Data() ) ; 
  else {
    rv = kTRUE ; 
    fYear = year ;  
    fProd = prod ; 
    fPath = path ; 
  }
  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSAliEnFile::SetVers(TString vers) 
{
  // set the year and verifies if the directory exists
  Bool_t rv = kFALSE ;
  char tempo[80] ; 
  sprintf(tempo, "/%s-%s/%s", fYear.Data(), fProd.Data(), vers.Data()) ; 
  fVers = tempo ;  

  TString path(fRoot) ; 
  path += tempo ; 
 if ( !fGrid->OpenDir(path) )  
    Error("ctor", "Cannot find directory %s ", path.Data() ) ; 
  else {
    rv = kTRUE ; 
    fVers = vers ; 
    fPath = path ;
  }
  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSAliEnFile::SetType(TString type) 
{
  // set the year and verifies if the directory exists
  Bool_t rv = kFALSE ;
  char tempo[80] ; 
  sprintf(tempo, "/%s-%s/%s/%s", fYear.Data(), fProd.Data(), fVers.Data(), type.Data()) ; 
 
  TString path(fRoot) ;
  path += tempo ; 
  if ( !fGrid->OpenDir(path) )  
    Error("ctor", "Cannot find directory %s ", path.Data() ) ; 
  else {
    rv = kTRUE ; 
    fType = type ;  
    fPath = path ; 
  }
  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSAliEnFile::SetPath(TString year, TString prod, TString vers, TString type) 
{
  // set the year and verifies if the directory exists
  Bool_t rv = kFALSE ; 
  char tempo[80] ; 
  sprintf(tempo, "/%s-%s/%s/%s", year.Data(), prod.Data(), vers.Data(), type.Data()) ; 
    
  TString path(fRoot) ; 
  path += tempo ; 
  if ( !fGrid->OpenDir(path) )  
    Error("ctor", "Cannot find directory %s ", path.Data() ) ; 
  else {
    rv = kTRUE ; 
    fPath = path ; 
    fYear += year ; 
    fProd += prod ; 
    fVers += vers ; 
    fType += type ;
  }
  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSAliEnFile::SetRun(Int_t run) 
{
  // set the year and verifies if the directory exists
  Bool_t rv = kFALSE ;

  TString zero("00000") ; 
  TString srun ; 
  srun += run ; 
  Int_t nzero = zero.Length() - srun.Length() ; 
  Int_t index ;
  for (index = 0 ; index < nzero ; index++) 
    srun.Prepend("0") ; 

  char tempo[80] ; 
  sprintf(tempo, "/%s-%s/%s/%s/%s", fYear.Data(), fProd.Data(), fVers.Data(), fType.Data(), srun.Data()) ; 

  TString path(fRoot) ; 
  path += tempo ;
  if ( !fGrid->OpenDir(path) )  
    Error("ctor", "Cannot find directory %s ", path.Data() ) ; 
  else {
    rv = kTRUE ; 
    fRun = srun ; 
    fPath = path ; 
  }
  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSAliEnFile::SetEvt(Int_t evt) 
{
  // set the year and verifies if the directory exists
  Bool_t rv = kFALSE ; 

  TString zero("00000") ; 
  TString sevt ; 
  sevt += evt ; 
  Int_t nzero = zero.Length() - sevt.Length() ; 
  Int_t index ;
  for (index = 0 ; index < nzero ; index++) 
    sevt.Prepend("0") ;
 
  char tempo[80] ; 
  sprintf(tempo, "/%s-%s/%s/%s/%s/%s/", fYear.Data(), fProd.Data(), fVers.Data(), fType.Data(), fRun.Data(), sevt.Data()) ; 
  TString path(fRoot) ; 
  path += tempo ;
  if ( !fGrid->OpenDir(path) )  
    Error("ctor", "Cannot find directory %s ", path.Data() ) ; 
  else {
    rv = kTRUE ; 
    fEvt = sevt ; 
    fPath = path ; 
  }
  return rv ; 
}


