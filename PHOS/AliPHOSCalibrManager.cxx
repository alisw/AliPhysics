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
#include "AliLog.h"  
#include "AliPHOSConTableDB.h"
#include "AliPHOSCalibrManager.h"
#include "AliPHOSCalibrationData.h"
ClassImp(AliPHOSCalibrManager)

  
  AliPHOSCalibrManager * AliPHOSCalibrManager::fgCaMa = 0 ;
//____________________________________________________________________________ 
AliPHOSCalibrManager::AliPHOSCalibrManager():TNamed() 
{
  // default ctor: not to be used
  fInputKind = 0 ;
  fFileName= "" ;
  fctdb = 0 ;
  AliFatal(Form("Should not be used")) ;
}
//____________________________________________________________________________ 
AliPHOSCalibrManager::AliPHOSCalibrManager(const char* filename,const char * kind ):
  TNamed("AliPHOSCalibrManager",filename){
  if(strcmp(kind,"root")==0){
    fInputKind = 0 ;
  }
  else{
    fInputKind = 1 ;
  }
  fFileName = filename;
  fctdb = 0 ;
} 
//____________________________________________________________________________ 
AliPHOSCalibrManager::~AliPHOSCalibrManager()
{
  //dtor
  if(fInputKind==0){
    TFile * f = gROOT->GetFile(fFileName) ;
    if(f && f->IsOpen())
      f->Close() ;
  }
}
//____________________________________________________________________________ 
AliPHOSCalibrManager * AliPHOSCalibrManager::GetInstance(void)
{ 
  // gets the instance of the unique object
  return fgCaMa ;
}
//____________________________________________________________________________ 
AliPHOSCalibrManager * AliPHOSCalibrManager::GetInstance(const char* filename,const char * kind)
{
  // Opens source (now only root file) and
  // returns the instance of the unique object
  
  if(!fgCaMa)
    fgCaMa = new AliPHOSCalibrManager(filename,kind) ;
  else{
    if(strcmp(filename,fgCaMa->fFileName.Data())){
      //Close previous file, construct new Manager
      delete fgCaMa ;
      fgCaMa = new AliPHOSCalibrManager(filename,kind) ;
    }
  }
  return fgCaMa ;	 
}
//____________________________________________________________________________
void AliPHOSCalibrManager::GetParameters(AliPHOSCalibrationData &data){ 
  if(fInputKind == 0){
    ReadFromRoot(data) ;
  }
  else{
    ReadFromASCII(data) ;
  }
}
//____________________________________________________________________________
void AliPHOSCalibrManager::ReadFromASCII(AliPHOSCalibrationData &data){
  // We read pedestals and gains from *.dat file with following format:
  //	0	0	0	0	37.09	1972.	// next nmodrows*nmodcols*ncryrows*ncrycols lines
  //	0	0	0	1	28.53	2072.	// contains <RR CC r c ped peak>
  //	0	0	0	2	30.93	1938.	//
  // where module is an array of 8*8 crystals and RR and CC are module raw and column position 

  if(!fctdb){
    AliError(Form("Specify Connections Table Database first")) ;
    return ;
  }
  
  FILE * file = fopen(fFileName, "r");
  if (!file) {
    Error("ReadFromASCII", "could not open file %s", fFileName.Data());
    return;
  }
  
  Bool_t isPed = kFALSE ;
  if(strcmp(data.GetCategory(),"Pedestals")==0){
    isPed = 1 ;
  }
  else{
    if(strcmp(data.GetCategory(),"Gains")){
      Error("ReadFromASCII","Unknown data category: %s",data.GetCategory()) ;
      return ;
    }
  }
  
  Int_t modRaw,modCol,raw,col;
  Float_t ped,pik;
  Int_t nread = 0 ;
  while(fscanf(file,"%d %d %d %d %f %f",&modRaw,&modCol,&raw,&col,&ped,&pik)==6){
    //Calculate plain crystal position:
    Int_t rawPosition = (modRaw*8+raw)*fctdb->GetNColumns()+modCol*8+col ;
    Int_t absId =  fctdb->Raw2AbsId(rawPosition);
    if(isPed){
      data.SetData(absId,ped) ;
    }
    else{
      if(pik!=0.)
	data.SetData(absId,1./pik);
      else
	data.SetData(absId,0.);
    }
    nread++ ;
  }    
  if(nread != fctdb->GetNColumns()*fctdb->GetNRaws()){
    Error("ReadFromASCII","Read %d parameters instead of %d\n",nread,fctdb->GetNColumns()*fctdb->GetNRaws()) ;
  }
  fclose(file) ;
  
}
//____________________________________________________________________________
void AliPHOSCalibrManager::ReadFromRoot(AliPHOSCalibrationData & data)
{
  //reads calibration parameters from root file
  
  //construct name
  TString searchname(data.GetCategory()) ;
  searchname+='_' ;
  searchname+=data.GetVersion() ;
  TFile * file = gROOT->GetFile(fFileName) ;
  if(!file || !file->IsOpen())
    file = TFile::Open(fFileName) ;
  if(!file->IsOpen()){
    AliError(Form("Can not open file %s\n",fFileName.Data() )) ;
    return ;
  }

  AliPHOSCalibrationData * tmp = dynamic_cast<AliPHOSCalibrationData*>(file->Get(searchname) ) ;
  if(!tmp)
    Error("ReadFromRoot","there is no data %s in file %s",fFileName.Data(),fFileName.Data()) ;
  else
    data=(*tmp);
  

//   TList * list = file->GetListOfKeys() ;
//   TIter next(list) ;
//   TKey * key ;
//   while((key=((TKey*)next()))){
//     TString kname(key->GetName()) ;
//     if(kname.BeginsWith(searchname)){
//       TString sbegin = kname(searchname.Sizeof()-1,
// 		             kname.Last('_')-searchname.Sizeof()+1) ;
//       TString send   = kname(kname.Last('_')+1,kname.Sizeof()) ;
//       Int_t begin,end ;
//       if(sscanf(sbegin.Data(),"%d",&begin) && 
// 	 sscanf(send.Data(),"%d",&end)){
// 	if(begin<=run && end >=run){
// 	  data=(*dynamic_cast<AliPHOSCalibrationData*>(file->Get(key->GetName()) ));
// 	  return ;
// 	}
//       }
//     }
//   }
//   Error("ReadFromRoot","Can not find key %s for run %d in file %s \n",searchname.Data(),run,fFileName.Data()) ;
}
//____________________________________________________________________________
void AliPHOSCalibrManager::WriteData(AliPHOSCalibrationData &data)
{
  //Writes data
  TFile * file = gROOT->GetFile(fFileName) ;
  if(!file || !file->IsOpen()){
    file = TFile::Open(fFileName,"UPDATE") ;
    if(!file->IsOpen()){
      AliError(Form("Can not open file %s\n",fFileName.Data() )) ;
      return ;
     }
  }
  file->cd() ;
  TString kname(data.GetCategory()) ;
  kname+='_';
  kname+=data.GetVersion() ;
//   kname+='_';
//   Int_t begin,end ;
//   data->GetValidityRange(begin,end) ;
//   kname+=begin ;
//   kname+='_' ;
//   kname+=end ;
  data.Write(kname,TObject::kOverwrite) ;
  
}
//____________________________________________________________________________
AliPHOSCalibrManager& AliPHOSCalibrManager::operator=(AliPHOSCalibrManager const & cdb)
{
  //overloads = operator
  fFileName = cdb.fFileName;
  return *this ;
}
