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
// Class provides interface between classes 
// consuming/producing calibration data, 
// such as AliPHOSClusterizer and AliPHOSCalibrator
// from one side and database (now just a root file, 
// later  - AliEn etc.)
// Use case: 
//  First Manager should be configured so that is reads data 
//  from a given source. e.g. from a root file:
//   AliPHOSCalibrManager * m = AliPHOSCalibrManager::GetInstance("MyCalibrationRootFile.root") ;
//  After that one can real/write data to this source:
//  AliPHOSCalibrationData mydata("Gains","v1") ;
//  ......
//  m->Write(&mydata) ;
//  AliPHOSCalibrationData myanotherdata("Gains","v2") ;
//  m->ReadFromRoot(&myanotherdata) ;
//
//*-- Author :  D.Peressounko (RRC KI & SUBATECH) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"

// --- Standard library ---
#include <fstream> 

// --- AliRoot header files ---
#include "AliPHOSConTableDB.h"
#include "AliPHOSCalibrManager.h"
#include "AliPHOSCalibrationData.h"
ClassImp(AliPHOSCalibrManager)


AliPHOSCalibrManager * AliPHOSCalibrManager::fgCaMa = 0 ;
//____________________________________________________________________________ 
AliPHOSCalibrManager::AliPHOSCalibrManager():TNamed() 
{
  // default ctor: not to be used
  fctdb = 0 ;
  fFileName="" ;
  Fatal("Default constructor","Should not use") ;
}
//____________________________________________________________________________ 
AliPHOSCalibrManager::AliPHOSCalibrManager(const char* filename ):
  TNamed("AliPHOSCalibrManager",filename){
  fFileName = filename ;	  
  fctdb = 0 ;
} 
//____________________________________________________________________________ 
  AliPHOSCalibrManager::~AliPHOSCalibrManager()
{
  //dtor
  TFile * f = gROOT->GetFile(fFileName) ;
  if(f && f->IsOpen())
    f->Close() ;
  if(fctdb)
    delete fctdb ;
  fctdb = 0 ;
}
//____________________________________________________________________________ 
AliPHOSCalibrManager * AliPHOSCalibrManager::GetInstance(void)
{ 
  // gets the instance of the unique object
 return fgCaMa ;
}
//____________________________________________________________________________ 
AliPHOSCalibrManager * AliPHOSCalibrManager::GetInstance(const char* filename )
{
  // Opens source (now only root file) and
  // returns the instance of the unique object

  if(!fgCaMa)
    fgCaMa = new AliPHOSCalibrManager(filename) ;
  else{
   if(strcmp(filename,fgCaMa->fFileName.Data())){
   //Close previous file, construct new Manager
     delete fgCaMa ;
     fgCaMa = new AliPHOSCalibrManager(filename) ;
   }
  }
  return fgCaMa ;	 
}
//____________________________________________________________________________
void AliPHOSCalibrManager::ReadFromASCII(AliPHOSCalibrationData &data,const char * filename)
{
  //reads calibration parameters from ascii file

  if(!fctdb){
    Error("ReadCalibrationParameters", "Specify Connections Table Database first") ;
    return ;
  }
  ifstream file(filename) ; 
  for(Int_t i = 1; i<=64; i++){
    Float_t inp ;
    file >> inp  ;      
    data.SetData(fctdb->Raw2AbsId(i),static_cast<Float_t>(inp));
  }
  file.close();   
}
//____________________________________________________________________________
void AliPHOSCalibrManager::ReadFromRoot(AliPHOSCalibrationData & data,Int_t run)
{
  //reads calibration parameters from root file

  //construct name
  TString searchname(data.GetCategory()) ;
  searchname+='_' ;
  searchname+=data.GetVersion() ;
  searchname+='_' ;
  TFile * file = gROOT->GetFile(fFileName) ;
  if(!file || !file->IsOpen())
    file = TFile::Open(fFileName) ;
  if(!file->IsOpen()){
    Error("ReadFromRoot","Can not open file %s\n",fFileName.Data() ) ;
    return ;
  }
  TList * list = file->GetListOfKeys() ;
  TIter next(list) ;
  TKey * key ;
  while((key=((TKey*)next()))){
    TString kname(key->GetName()) ;
    if(kname.BeginsWith(searchname)){
      TString sbegin = kname(searchname.Sizeof()-1,
		             kname.Last('_')-searchname.Sizeof()+1) ;
      TString send   = kname(kname.Last('_')+1,kname.Sizeof()) ;
      Int_t begin,end ;
      if(sscanf(sbegin.Data(),"%d",&begin) && 
	 sscanf(send.Data(),"%d",&end)){
	if(begin<=run && end >=run){
	  data=(*dynamic_cast<AliPHOSCalibrationData*>(file->Get(key->GetName()) ));
	  return ;
	}
      }
    }
  }
  Error("ReadFromRoot","Can not find key %s for run %d in file %s \n",searchname.Data(),run,fFileName.Data()) ;
}
//____________________________________________________________________________
void AliPHOSCalibrManager::WriteData(AliPHOSCalibrationData * data)
{
  //Writes data
  TFile * file = gROOT->GetFile(fFileName) ;
  if(!file || !file->IsOpen()){
    file = TFile::Open(fFileName,"UPDATE") ;
    if(!file->IsOpen()){
      Error("WriteData","Can not open file %s\n",fFileName.Data() ) ;
      return ;
     }
  }
  file->cd() ;
  TString kname(data->GetCategory()) ;
  kname+='_';
  kname+=data->GetVersion() ;
  kname+='_';
  Int_t begin,end ;
  data->GetValidityRange(begin,end) ;
  kname+=begin ;
  kname+='_' ;
  kname+=end ;
  data->Write(kname,TObject::kOverwrite) ;
  
}
//____________________________________________________________________________
AliPHOSCalibrManager& AliPHOSCalibrManager::operator=(AliPHOSCalibrManager const & cdb)
{
  //overloads = operator
  fFileName = cdb.fFileName;
  return *this ;
}
