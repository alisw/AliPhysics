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

//-----------------------------------------------------------------
//                AliGRPDCS class
//   This class deals with the DCS related info of the GRP
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//    Modified: Ernesto.Lopez.Torres@cern.ch  CEADEN-CERN
//-----------------------------------------------------------------
#include "Riostream.h"
#include "TError.h"

#include "AliGRPDCS.h"
#include "AliDCSValue.h"
#include "AliLog.h"

class TObjString;

#include <TObjArray.h>

ClassImp(AliGRPDCS)

//_______________________________________________________________
AliGRPDCS::AliGRPDCS():
  TObject(), fStartTime(0), fStopTime(0), fDCSArray(new TObjArray())
{
  // default constructor
  
}

//_______________________________________________________________
AliGRPDCS::AliGRPDCS(TObjArray *dcsArray, UInt_t fStart, UInt_t fStop):
  TObject(), fStartTime(fStart), fStopTime(fStop), fDCSArray(new TObjArray()) {
  // constructor
  fDCSArray = dcsArray;
}

//___________________________________________________________________________
AliGRPDCS::AliGRPDCS(const AliGRPDCS& grpDcs):
  TObject(grpDcs), 
  fStartTime(grpDcs.fStartTime), fStopTime(grpDcs.fStopTime), 
  fDCSArray(grpDcs.fDCSArray)
{
  //copy constructor

  //if (grpDcs.fDCSArray) fDCSArray = new TObjArray();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessDCS(Int_t iType)
{
  // process the dcs dps
  if(!fDCSArray->GetEntries()) {
    Error("AliGRPDCS::ProcessDCS","No entries in array!!!");
    return 0;
  }

  // 1:Bool - 2:Char - 3:Int - 4:UInt - 5:Float
  TString fDCSDataPointValue;
  switch (iType) {
    case 1: {
      fDCSDataPointValue += ProcessBoolean();
      break;
    }
    case 2: {
      fDCSDataPointValue += ProcessChar();
      break;
    }
    case 3: {
      fDCSDataPointValue += ProcessInt();
      break;
    }
    case 4: {
      fDCSDataPointValue += ProcessUInt();
      break;
    }
    case 5: {
      fDCSDataPointValue += ProcessFloat();
      break;
    }
//  case 6: {
//    fDCSDataPointValue += ProcessString();
//    break;
//  }
    default: {
      Error("AliGRPDCS::ProcessDCS", "Unknown type!!!" );
      break;
    }
  }//switch

  //cout<<fDCSDataPointValue.Data()<<endl;
  return fDCSDataPointValue.Data();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessBoolean()
{
  Bool_t fDCSBool = kTRUE;
  Bool_t previousBool = kTRUE;

  AliDCSValue *v = 0x0;

  for(Int_t iCount = 0; iCount < fDCSArray->GetEntries(); iCount++) {
    v = (AliDCSValue *)fDCSArray->At(iCount);
    if ((v->GetTimeStamp() < fStartTime) || (v->GetTimeStamp() > fStopTime)) {
      AliError(Form("DCS values for the parameter outside the queried interval"));
      continue;
    }
    fDCSBool = v->GetBool();
    if (iCount > 0) {
      if (fDCSBool != previousBool)
	AliError(Form("DCS values for the parameter changed from %d to %d within the queried interval", (Int_t)previousBool, (Int_t)fDCSBool));
    }
    previousBool = fDCSBool;
  }

  TString fDCSDataPointValue = (fDCSBool)? "1" : "0";

  return fDCSDataPointValue.Data();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessInt()
{
  Float_t fFDCSArraySum = 0.0, fFDCSArrayMean = 0.0;
  Int_t iCounts = 0;
  for(Int_t i = 0; i < fDCSArray->GetEntries(); i++) {
    AliDCSValue *v = (AliDCSValue *)fDCSArray->At(i);
    if((v->GetTimeStamp() >= fStartTime) &&(v->GetTimeStamp() <= fStopTime)) {
      fFDCSArraySum += v->GetInt();
      iCounts += 1;
    }
  }
  if(iCounts != 0) fFDCSArrayMean = fFDCSArraySum/iCounts;
  else fFDCSArrayMean = -10.;
 
  TString fDCSDataPointValue; fDCSDataPointValue += fFDCSArrayMean;

  return fDCSDataPointValue.Data();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessUInt()
{
  Float_t fFDCSArraySum = 0.0, fFDCSArrayMean = 0.0;
  Int_t iCounts = 0;
  for(Int_t i = 0; i < fDCSArray->GetEntries(); i++) {
    AliDCSValue *v = (AliDCSValue *)fDCSArray->At(i);
    if((v->GetTimeStamp() >= fStartTime) &&(v->GetTimeStamp() <= fStopTime)) {
      fFDCSArraySum += v->GetUInt();
      iCounts += 1;
    }
  }
  if(iCounts != 0) fFDCSArrayMean = fFDCSArraySum/iCounts;
  else fFDCSArrayMean = -10.;
   
  TString fDCSDataPointValue; fDCSDataPointValue += fFDCSArrayMean;

  return fDCSDataPointValue.Data();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessFloat()
{
  Float_t fFDCSArraySum = 0.0, fFDCSArrayMean = 0.0;
  Int_t iCounts = 0;
  //printf("Entries: %d\n",fDCSArray->GetEntries());
  //printf("Start: %d - Stop: %d\n",fStartTime,fStopTime);
  for(Int_t i = 0; i < fDCSArray->GetEntries(); i++) {
    AliDCSValue *v = (AliDCSValue *)fDCSArray->At(i);
    if((v->GetTimeStamp() >= fStartTime) &&(v->GetTimeStamp() <= fStopTime)) {
      fFDCSArraySum += v->GetFloat();
      iCounts += 1;
    }
  }
  if(iCounts != 0) fFDCSArrayMean = fFDCSArraySum/iCounts;
  else fFDCSArrayMean = -10.;
   
  TString fDCSDataPointValue; fDCSDataPointValue += fFDCSArrayMean;

  return fDCSDataPointValue.Data();
}

//_______________________________________________________________
const char* AliGRPDCS::ProcessChar()
{
  TString fDCSString;

  AliDCSValue *v = 0x0;

  //printf("Entries: %d\n",fDCSArray->GetEntries());
  //printf("Start: %d - Stop: %d\n",fStartTime,fStopTime);
  
  for(Int_t iCount = 0; iCount < fDCSArray->GetEntries(); iCount++) {
    v = (AliDCSValue *)fDCSArray->At(iCount);
    if ((v->GetTimeStamp() >= fStartTime) && (v->GetTimeStamp() <= fStopTime)) 
      AliError(Form("DCS values for the parameter changed within the queried interval"));
    if (v->GetTimeStamp() > fStopTime) continue;
    fDCSString = v->GetChar();
  }
  
  TString fDCSDataPointValue = fDCSString;
  
  return fDCSDataPointValue.Data();
}

/*
//_______________________________________________________________
const char* AliGRPDCS::ProcessString()
{
  TString fDCSString;

  AliDCSValue *v = 0x0;

  //printf("Entries: %d\n",fDCSArray->GetEntries());
  //printf("Start: %d - Stop: %d\n",fStartTime,fStopTime);
  
  for(Int_t iCount = 0; iCount < fDCSArray->GetEntries(); iCount++) {
    v = (TObjString *)fDCSArray->At(iCount);
    if ((v->GetTimeStamp() >= fStartTime) && (v->GetTimeStamp() <= fStopTime)) 
      AliError(Form("DCS values for the parameter changed within the queried interval"));
    if (v->GetTimeStamp() > fStopTime) continue;
    fDCSString = v->GetChar(); // FIXME  should be something like v->GetString()
  }
  
  TString fDCSDataPointValue = fDCSString;
  //cout<<"Returned String: "<<fDCSDataPointValue.Data()<<endl;

  return fDCSDataPointValue.Data();
}
*/
