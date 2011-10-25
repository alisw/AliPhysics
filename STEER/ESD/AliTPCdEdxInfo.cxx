#include "AliTPCdEdxInfo.h"

//##################################################################
//
// Simple class to store TPC dE/dx info for different pad regions.
//
// Origin: Marian Ivanov, Alexander Kalweit 
//
//##################################################################



ClassImp(AliTPCdEdxInfo)

AliTPCdEdxInfo::AliTPCdEdxInfo():
  TObject(),
  fTPCsignalRegion(),
  fTPCsignalNRegion(),
  fTPCsignalNRowRegion()
{
  // Default constructor
  for (Int_t i=0;i<3; i++){
    fTPCsignalRegion[i]=0;
    fTPCsignalNRegion[i]=0;
    fTPCsignalNRowRegion[i]=0;
  }
  fTPCsignalRegion[3]=0;
 
}


AliTPCdEdxInfo::AliTPCdEdxInfo(const AliTPCdEdxInfo& source):
  TObject(),
  fTPCsignalRegion(),
  fTPCsignalNRegion(),
  fTPCsignalNRowRegion()
{
  //
  // copy constructor
  //
  Double32_t signal[4]; Char_t ncl[3]; Char_t nrows[3];
  source.GetTPCSignalRegionInfo(signal, ncl, nrows);
  for (Int_t i=0;i<3; i++){
    fTPCsignalRegion[i]=signal[i];
    fTPCsignalNRegion[i]=ncl[i];
    fTPCsignalNRowRegion[i]=nrows[i];
  }
  fTPCsignalRegion[3]=signal[3];
 
}



void  AliTPCdEdxInfo::GetTPCSignalRegionInfo(Double32_t signal[4], Char_t ncl[3], Char_t nrows[3]) const {
  //
  // Get the TPC dEdx variables per region
  //
  // Double32_t  fTPCsignalRegion[4]; // TPC dEdx signal in 4 different regions - 0 - IROC, 1- OROC medium, 2 - OROC long, 3- OROC all, (default truncation used)  
  // Char_t      fTPCsignalNRegion[3]; // number of clusters above threshold used in the dEdx calculation
  // Char_t      fTPCsignalNRowRegion[3]; // number of crosed rows used in the dEdx calculation - signal below threshold included
  //
  for (Int_t i=0; i<3; i++){
    signal[i]=fTPCsignalRegion[i];
    ncl[i]=fTPCsignalNRegion[i];
    nrows[i]=fTPCsignalNRowRegion[i];
  }
  signal[3]=fTPCsignalRegion[3];
  return; 
}


void  AliTPCdEdxInfo::SetTPCSignalRegionInfo(Double32_t signal[4], Char_t ncl[3], Char_t nrows[3]){
  //
  // Set the TPC dEdx variables per region
  //
  // Double32_t  fTPCsignalRegion[4]; // TPC dEdx signal in 4 different regions - 0 - IROC, 1- OROC medium, 2 - OROC long, 3- OROC all, (default truncation used)  
  // Char_t      fTPCsignalNRegion[3]; // number of clusters above threshold used in the dEdx calculation
  // Char_t      fTPCsignalNRowRegion[3]; // number of crosed rows used in the dEdx calculation - signal below threshold included
  //
  for (Int_t i=0;i<3; i++){
    fTPCsignalRegion[i]=signal[i];
    fTPCsignalNRegion[i]=ncl[i];
    fTPCsignalNRowRegion[i]=nrows[i];
  }
  fTPCsignalRegion[3]=signal[3];
  return;
}
