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
// Short description  
//
//*-- Author :  D.Peressounko (RRC KI & SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TFile.h"
// --- Standard library ---
#include <fstream>

// --- AliRoot header files ---
#include "AliPHOSConTableDB.h"
#include "AliPHOSCalibrationDB.h"
#include "AliPHOSGeometry.h"
ClassImp(AliPHOSCalibrationDB)


//____________________________________________________________________________ 
  AliPHOSCalibrationDB::AliPHOSCalibrationDB():TNamed() 
{
  fPedestals = 0 ;
  fSlopes = 0;
  fNChannels = 0 ;
  fctdb = 0 ;
}
//____________________________________________________________________________ 
AliPHOSCalibrationDB::AliPHOSCalibrationDB(const char * database):
  TNamed("AliPHOSCalibrationDB",database){
  //Creates the containers: we prepare places for all channels in PHOS

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance("IHEP","") ;
  fNChannels = geom->GetNModules()*geom->GetNPhi()*geom->GetNZ()+
               geom->GetNumberOfCPVPadsZ()*geom->GetNumberOfCPVPadsPhi()*geom->GetNModules() ;
  //Note, that to avoid shifting AbsId each time we access data, we do not use
  //first slot and index runs from 1 to fNChannels inclusively
  fPedestals = new TArrayF(fNChannels+1) ;
  fSlopes =  new TArrayF(fNChannels+1) ;
  fctdb = 0 ;
} 
//____________________________________________________________________________ 
  AliPHOSCalibrationDB::~AliPHOSCalibrationDB()
{
  if(fPedestals)
    delete fPedestals ;
  if(fSlopes)
    delete fSlopes ;
}
//____________________________________________________________________________ 
Float_t AliPHOSCalibrationDB::Calibrate(Int_t amp, Int_t absId)const 
{
  //if absID is known, return calibrated energy, else - zero
  if(absId <=0 || absId > fNChannels) 
    return 0.0000001 ;

  Float_t ret = (amp - fPedestals->At(absId))*fSlopes->At(absId); 
  if(ret > 0)
    return ret ;
  else
    return 0.0000001 ; //Should not be zero - to avoid FPE
}
//____________________________________________________________________________ 
void  AliPHOSCalibrationDB::SetAll(Float_t pedestal, Float_t slope){
  //Set all calibration parameters to the same value
  if(fPedestals)
    for(Int_t i=0;i<fNChannels;i++){
      fPedestals->AddAt(pedestal,i) ;
      fSlopes->AddAt(slope,i);
    }
  else
    Warning("SetAll", "Please, create me with non-default constructor!") ;
}
//____________________________________________________________________________
void AliPHOSCalibrationDB::ReadCalibrationParameters(const char * filename, Option_t* opt){
  //reads calibration parameters from ascii file

  if(strcmp(opt,"gains")==0){  //read gains
    if(!fctdb){
      Error("ReadCalibrationParameters", "Specify Connections Table Database first") ;
      return ;
    }
    ifstream gainfile(filename) ; 
    for(Int_t i = 1; i<=64; i++){
      Float_t slope ;
      gainfile >> slope  ;      
      fSlopes->AddAt(slope,fctdb->Raw2AbsId(i));
    }
    gainfile.close();   
  }
  else
    if(strstr(opt,"pedest")){  //read pedestals
      if(!fctdb){
	Error("ReadCalibrationParameters", "Specify Connections Table Database first") ;
	return ;
      }
      ifstream pfile(filename) ; 
      for(Int_t i = 1; i<=64; i++){
	Float_t pedestal ;
	pfile >> pedestal  ;      
	fPedestals->AddAt(pedestal,fctdb->Raw2AbsId(i));
      }
      pfile.close();
   }
   else{
     Warning("ReadCalibrationParameters", "Available options are\n `gains' : to read gains\n `pedestals : to read pedestals ") ;
   }
}
//____________________________________________________________________________
AliPHOSCalibrationDB& AliPHOSCalibrationDB::operator=(AliPHOSCalibrationDB const & cdb)
{
  //
  fNChannels = cdb.fNChannels;
  fFileName = cdb.fFileName;
  fPedestals=cdb.fPedestals ;
  fSlopes = cdb.fSlopes;
  return *this ;
}
