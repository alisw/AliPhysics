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


//-------------------------------------------------------------------------
//  This class is used to store data in events where scalers are read
//   in SOFTWARE_TRIGGER events, no ADC datum filled in such events
//  Author: Chiara Oppedisano	Chiara.Oppedisano@to.infn.it	
//-------------------------------------------------------------------------

#include "AliESDZDCScalers.h"

ClassImp(AliESDZDCScalers)

//______________________________________________________________________________
AliESDZDCScalers::AliESDZDCScalers() :
  TObject()
{
  // constructor
  for(Int_t i=0; i<2; i++)  fScalerDown[i]=0;
  for(Int_t i=0; i<8; i++)  fScalerUp[i]=0;
  for(Int_t i=0; i<32; i++) fVMEScaler[i]=0;

}//______________________________________________________________________________
AliESDZDCScalers::AliESDZDCScalers(const AliESDZDCScalers& zdc) :
  TObject(zdc)
{
  // copy constructor
  for(Int_t i=0; i<32; i++){
     if(i<2) fScalerDown[i] = zdc.fScalerDown[i];
     if(i<8) fScalerUp[i] = zdc.fScalerUp[i];
     fVMEScaler[i] = zdc.fVMEScaler[i];
  }
}

//______________________________________________________________________________
AliESDZDCScalers& AliESDZDCScalers::operator=(const AliESDZDCScalers&zdc)
{
  // assigment operator
  if(this!=&zdc) {
    TObject::operator=(zdc);
    for(Int_t i=0; i<32; i++){
       if(i<2) fScalerDown[i] = zdc.fScalerDown[i];
       if(i<8) fScalerUp[i] = zdc.fScalerUp[i];
       fVMEScaler[i] = zdc.fVMEScaler[i];
    }
  } 
  return *this;
}

//______________________________________________________________________________
void AliESDZDCScalers::Copy(TObject &obj) const {
  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDZDCScalers *robj = dynamic_cast<AliESDZDCScalers*>(&obj);
  if(!robj)return; // not an AliESDZDCScalers
  *robj = *this;

}


//______________________________________________________________________________
void AliESDZDCScalers::Reset()
{
  // reset all data members
  for(Int_t i=0; i<32; i++){
     if(i<2) fScalerDown[i] = 0;
     if(i<8) fScalerUp[i] = 0;
     fVMEScaler[i] = 0;
  }
}

//______________________________________________________________________________
void AliESDZDCScalers::Print(const Option_t *) const
{
  //  Print ESD scaler events for the ZDC
  printf("\n");
  for(Int_t i=0; i<2; i++) printf("\tfScalerDown[%d] = %d \n",i,fScalerDown[i]);
  for(Int_t i=0; i<8; i++) printf("\tfScalerUp[%d] = %d \n",i,fScalerUp[i]);
  for(Int_t i=0; i<32; i++) printf("\tfVMEScaler[%d] = %d \n",i,fVMEScaler[i]);
  printf("\n");
}
