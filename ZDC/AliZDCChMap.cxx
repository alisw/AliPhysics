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

///////////////////////////////////////////////////////////////
//                                                           //
// Class for ZDC calibration -> ADC channels mapping         //
// author: Chiara Oppedisano			             //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliZDCChMap.h"

ClassImp(AliZDCChMap)

//________________________________________________________________
AliZDCChMap::AliZDCChMap():
TNamed()
{
  Reset();
}

//________________________________________________________________
AliZDCChMap::AliZDCChMap(const char* name):
TNamed()
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  int const kNModules = 10;
  int const kNChannels = 48;
  int const kNScChannels = 32;
  for(Int_t i=0; i<kNModules; i++){
    for(Int_t j=0; j<3; j++) fModuleMap[i][j] = 0;
  }
  for(Int_t i=0; i<kNChannels; i++){
    fADCModule[i] = -1;
    fADCChannel[i] = -1;
    fDetector[i] = -1;
    fSector[i] = -1;
    fADCSignalCode[i] = -1;
  }
  for(Int_t i=0; i<kNScChannels; i++){
    fScalerChannel[i] = -1;
    fScDetector[i] = -1;
    fScSector[i] = -1;
    fScSignalCode[i] = -1;
    //
    fTDCChannel[i] = -1;
    fTDCSignalCode[i] = -1;
  }
  
  
}

//________________________________________________________________
AliZDCChMap::AliZDCChMap(const AliZDCChMap& calibda) :
  TNamed(calibda)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  int const kNModules = 10;
  int const kNChannels = 48;
  int const kNScChannels = 32;
  for(Int_t i=0; i<kNModules; i++){
     for(Int_t j=0; j<3; j++) fModuleMap[i][j] = calibda.GetModuleMap(i,j);
  }
  for(int t=0; t<kNChannels; t++){
     fADCModule[t]  = calibda.GetADCModule(t);
     fADCChannel[t] = calibda.GetADCChannel(t);
     fDetector[t]   = calibda.GetDetector(t);
     fSector[t]     = calibda.GetSector(t);
     fADCSignalCode[t]  = calibda.GetADCSignalCode(t);
     if(t<kNScChannels){
       fScalerChannel[t] = calibda.GetScChannel(t);
       fScDetector[t]    = calibda.GetScDetector(t);
       fScSector[t]      = calibda.GetScSector(t);
       fScSignalCode[t]  = calibda.GetScSignalCode(t);
       //
       fTDCChannel[t] = calibda.GetTDCChannel(t);
       fTDCSignalCode[t] = calibda.GetTDCChannel(t);
     }
  }
}

//________________________________________________________________
AliZDCChMap &AliZDCChMap::operator =(const AliZDCChMap& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  int const kNModules = 10;
  int const kNChannels = 48;
  int const kNScChannels = 32;
  for(Int_t i=0; i<kNModules; i++){
     for(Int_t j=0; j<3; j++) fModuleMap[i][j] = calibda.GetModuleMap(i,j);
  }
  for(int t=0; t<kNChannels; t++){
     fADCModule[t]  = calibda.GetADCModule(t);
     fADCChannel[t] = calibda.GetADCChannel(t);
     fDetector[t]   = calibda.GetDetector(t);
     fSector[t]     = calibda.GetSector(t);
     fADCSignalCode[t]  = calibda.GetADCSignalCode(t);
     if(t<kNScChannels){
       fScalerChannel[t] = calibda.GetScChannel(t);
       fScDetector[t]    = calibda.GetScDetector(t);
       fScSector[t]      = calibda.GetScSector(t);
       fScSignalCode[t]  = calibda.GetScSignalCode(t);
       //
       fTDCChannel[t] = calibda.GetTDCChannel(t);
       fTDCSignalCode[t] = calibda.GetTDCChannel(t);
     }
  }

  return *this;
}

//________________________________________________________________
AliZDCChMap::~AliZDCChMap()
{
}

//________________________________________________________________
void AliZDCChMap::Reset()
{
  // Reset
  memset(fADCModule,0,48*sizeof(Int_t));
  memset(fADCChannel,0,48*sizeof(Int_t));
  memset(fDetector,0,48*sizeof(Int_t));
  memset(fSector,0,48*sizeof(Int_t));
  memset(fADCSignalCode,0,48*sizeof(Int_t));
  memset(fScalerChannel,0,32*sizeof(Int_t));
  memset(fScDetector,0,32*sizeof(Int_t));
  memset(fScSector,0,32*sizeof(Int_t));
  memset(fScSignalCode,0,32*sizeof(Int_t));
  memset(fTDCChannel,0,32*sizeof(Int_t));
  memset(fTDCSignalCode,0,32*sizeof(Int_t));
}                                                                                       


//________________________________________________________________
void  AliZDCChMap::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n\t ******************* AliZDCChMap object *******************\n\n");
   for(Int_t i=0; i<10; i++){
     printf("  ******** GEO %d mod. type %d #ch. %d\n",
      fModuleMap[i][0],fModuleMap[i][1],fModuleMap[i][2]);     
   } 
   printf("\n");
   for(Int_t i=0; i<48; i++) 
     printf(" ADC - mod. %d ch. %d signal %d -> detector %d sector %d\n",
      fADCModule[i], fADCChannel[i], fADCSignalCode[i], fDetector[i], fSector[i]);
   printf("\n");
   for(Int_t i=0; i<32; i++)
     if(fScalerChannel[i]!=-1)
       printf(" SCALER - ch. %d signal %d\n",
        fScalerChannel[i], fScSignalCode[i]);
   printf("\n");
   for(Int_t i=0; i<32; i++) 
     if(fTDCChannel[i]!=-1)
       printf(" TDC - ch. %d signal %d\n",
        fTDCChannel[i], fTDCSignalCode[i]);
   printf("\n\t **********************************************************\n\n");
 
} 
