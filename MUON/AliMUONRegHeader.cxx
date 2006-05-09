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

#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"



/// 
/// Regional structure for trigger raw data.
/// Each Reg structure contains 16 (at most) local card structure.
/// The structure includes the information of the Reg. boards and
/// regional inputs
/// 

ClassImp(AliMUONRegHeader)
 
 const Int_t  AliMUONRegHeader::fgkHeaderLength = 3;
 const Int_t  AliMUONRegHeader::fgkScalerLength = 8;
 const UInt_t AliMUONRegHeader::fgkEndOfReg     = 0xBEEFFACE;

//___________________________________________
AliMUONRegHeader::AliMUONRegHeader()
  :  TObject(),
     fWord(0),
     fL0(0), 
     fClk(0),
     fHold(0)

{
  //
  // ctor
  //
  fInput[0] = fInput[1] = 0;

  for (Int_t i = 0; i < 8; i++)
    fScaler[i] = 0;

  fLocalArray  = new TClonesArray("AliMUONLocalStruct",16);

}

//___________________________________________
AliMUONRegHeader::~AliMUONRegHeader()
{
  // 
  // dtor
  //
  fLocalArray->Delete();
  delete fLocalArray;
}

//___________________________________________
AliMUONRegHeader::AliMUONRegHeader(const AliMUONRegHeader& event)
  :  TObject(event)
{
  //
  // copy ctor
  //
  fWord = event.fWord;
  fL0   = event.fL0;
  fClk  = event.fClk;
  fHold = event.fHold;

  fInput[0] = event.fInput[0];
  fInput[1] = event.fInput[1];

  for (Int_t i = 0; i < 8; i++)
    fScaler[i] = event.fScaler[i];

  fLocalArray = new TClonesArray("AliMUONLocalStruct", 16);
  for (Int_t index = 0; index < (event.fLocalArray)->GetEntriesFast(); index++) {
    {new ((*fLocalArray)[fLocalArray->GetEntriesFast()]) 
        AliMUONLocalStruct(*(AliMUONLocalStruct*)(event.fLocalArray)->At(index));}
  }
}

//___________________________________________
AliMUONRegHeader& AliMUONRegHeader::operator=(const AliMUONRegHeader& event)
{
  // 
  // assignment operator
  //

  if (this == &event) return *this;

  fWord = event.fWord;
  fL0   = event.fL0;
  fClk  = event.fClk;
  fHold = event.fHold;

  fInput[0] = event.fInput[0];
  fInput[1] = event.fInput[1];

  for (Int_t i = 0; i < 8; i++)
    fScaler[i] = event.fScaler[i];

  fLocalArray = new TClonesArray("AliMUONLocalStruct", 16);
  for (Int_t index = 0; index < (event.fLocalArray)->GetEntriesFast(); index++) {
    {new ((*fLocalArray)[fLocalArray->GetEntriesFast()]) 
        AliMUONLocalStruct(*(AliMUONLocalStruct*)(event.fLocalArray)->At(index));}
  }

  return *this;
}

//___________________________________________
void AliMUONRegHeader::SetScalersNumbers()
{
  // set numbers for scaler events for Regional header
  // since this is provided by the experiment
  // put dummy numbers to check the monitoring
  
  fL0   = 1000;
  fClk  = 10000;
  fHold = 100; 
  
  for (Int_t i = 0; i < 8; i++)
    fScaler[i] = i;  
}

//___________________________________________
void AliMUONRegHeader::Clear(Option_t* )
{
  // Clear TClones arrays
  // instead of deleting
  //
  fLocalArray->Clear("C");
 
}
