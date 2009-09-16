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
  for(Int_t i=0; i<48; i++){
    fADCModule[i] = 0;
    fADCChannel[i] = 0;
    fDetector[i] = 0;
    fSector[i] = 0;
    fScalerChannel[i] = 0;
    fScDetector[i] = 0;
    fScSector[i] = 0;
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
  for(int t=0; t<48; t++){
     fADCModule[t]  = calibda.GetADCModule(t);
     fADCChannel[t] = calibda.GetADCChannel(t);
     fDetector[t]   = calibda.GetDetector(t);
     fSector[t]     = calibda.GetSector(t);
     if(t<32){
       fScalerChannel[t] = calibda.GetScChannel(t);
       fScDetector[t]    = calibda.GetScDetector(t);
       fScSector[t]      = calibda.GetScSector(t);
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
  for(int t=0; t<48; t++){
     fADCModule[t]  = calibda.GetADCModule(t);
     fADCChannel[t] = calibda.GetADCChannel(t);
     fDetector[t]   = calibda.GetDetector(t);
     fSector[t]     = calibda.GetSector(t);
     if(t<32){
       fScalerChannel[t] = calibda.GetScChannel(t);
       fScDetector[t]    = calibda.GetScDetector(t);
       fScSector[t]      = calibda.GetScSector(t);
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
  memset(fScalerChannel,0,32*sizeof(Int_t));
  memset(fScDetector,0,32*sizeof(Int_t));
  memset(fScSector,0,32*sizeof(Int_t));
}                                                                                       


//________________________________________________________________
void  AliZDCChMap::Print(Option_t *) const
{
   // Printing of calibration object
   printf("\n\n\t ******************* AliZDCChMap object *******************\n\n");
   for(Int_t i=0; i<48; i++) 
     printf(" ADC - mod. %d ch. %d -> detector %d sector %d\n",
      fADCModule[i], fADCChannel[i],fDetector[i], fSector[i]);
   for(Int_t i=0; i<32; i++) 
     printf(" SCALER - ch. %d -> detector %d sector %d\n",
      fScalerChannel[i],fScDetector[i], fScSector[i]);
   printf("\n\n\t **********************************************************\n\n");
 
} 
