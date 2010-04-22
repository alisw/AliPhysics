/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notifce   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id: AliJRunHeader.cxx,v 1.1 2008/02/04 13:28:47 rak Exp $

////////////////////////////////////////////////////
//
//	\file AliJRunHeader.cxx
//	\brief
//	\author J. Rak, D.J.Kim, F.Krizek(University of Jyvaskyla)
//	\email: djkim@cc.jyu.fi
//	\version $Revision: 1.1 $
//	\date $Date: 2008/02/04 13:28:47 $
//
// Class encapsulation aliroot run header information
////////////////////////////////////////////////////
#include "AliJRunHeader.h"

ClassImp(AliJRunHeader)

//----------------------------
AliJRunHeader::AliJRunHeader():
  TNamed(),
  fRunNumber(-1),
  fL3MagnetPolarity(0),
  fMagneticFieldL3(0),
  fActiveTriggersAlice(),
  fSizeOfTableJCorran(0),
  fActiveTriggersJCorran()
{                     //constructor
 
  for(Int_t i=0;i<kRangeTriggerTableAlice;i++){
    fActiveTriggersAlice.Add(new TObjString("EMPTY"));
    fActiveTriggersJCorran.Add(new TObjString("EMPTY"));
  }

  fActiveTriggersAlice.SetOwner(kTRUE);
  fActiveTriggersJCorran.SetOwner(kTRUE);

  SetName("AliRunHeader");
  SetTitle("AliRunHeader");

}
//________________________________________________________________________

AliJRunHeader::AliJRunHeader(const AliJRunHeader& ap):
  TNamed(ap),
  fRunNumber(ap.fRunNumber),
  fL3MagnetPolarity(ap.fL3MagnetPolarity),
  fMagneticFieldL3(ap.fMagneticFieldL3),  
  fActiveTriggersAlice(ap.fActiveTriggersAlice),       
  fSizeOfTableJCorran(ap.fSizeOfTableJCorran),
  fActiveTriggersJCorran(ap.fActiveTriggersJCorran)
{
  //cpy ctor
}
//________________________________________________________________________

void AliJRunHeader::SetActiveTriggersAlice(TString *triggers){
  // fill aliroot trigger table
  for(Int_t t=0;t<kRangeTriggerTableAlice;t++){
    ((TObjString*) (fActiveTriggersAlice.At(t)))->SetString(triggers[t].Data());
  }

  fActiveTriggersAlice.SetOwner(kTRUE);
}
//________________________________________________________________________

void AliJRunHeader::SetActiveTriggersJCorran(TString *triggers, Int_t range){
  //fill jcorran trigger table
  for(Int_t t=0;t<range;t++){
    ((TObjString*) (fActiveTriggersJCorran.At(t)))->SetString(triggers[t].Data());
  }
  fActiveTriggersJCorran.SetOwner(kTRUE);
  fSizeOfTableJCorran = range;
}
//________________________________________________________________________

Int_t AliJRunHeader::GetActiveTriggerBitAlice(TString TriggerName){
  //get trigger bit corresponding to trigger name
  Int_t tbit=-1;
  for(Int_t t=0;t<kRangeTriggerTableAlice;t++){
    if(TriggerName.Contains(((TObjString*)(fActiveTriggersAlice.At(t)))->GetString())){
      tbit = t;
      break;
    }
  }
  return tbit;
}
//________________________________________________________________________

void AliJRunHeader::PrintOut(){
  //print object
  cout<<"RUN "<<fRunNumber<<endl;
  cout<<"Alice trigger table:"<<endl;
  for(Int_t i=0;i<kRangeTriggerTableAlice;i++){
    cout<<i<<"   "<<GetActiveTriggerAlice(i)<<endl;
  }
  cout<<"============================="<<endl;
  cout<<"JCorran trigger table:"<<endl;
  for(Int_t i=0;i<fSizeOfTableJCorran;i++){
    cout<<i<<"   "<<GetActiveTriggerJCorran(i)<<endl;
  }
  cout<<"============================="<<endl;
  cout<<"Magnet polarity "<<fL3MagnetPolarity<<endl;
  cout<<"B "<<fMagneticFieldL3<<endl;
  cout<<"============================="<<endl;

}
//_______________________________________________________________________
AliJRunHeader& AliJRunHeader::operator=(const  AliJRunHeader& header){
  //operator =
  if(this != &header) {
    TNamed::operator=(header);
    fRunNumber        = header.fRunNumber;
    fL3MagnetPolarity = header.fL3MagnetPolarity;
    fMagneticFieldL3  = header.fMagneticFieldL3;
    fActiveTriggersAlice   = header.fActiveTriggersAlice;
    fSizeOfTableJCorran    = header.fSizeOfTableJCorran;
    fActiveTriggersJCorran = header.fActiveTriggersJCorran;
  }
  return *this;
}


