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

#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"

/// \class AliMUONDarcHeader
/// Darc structure for trigger raw data.
/// Each DDL contains one Darc structure
/// The structure includes the information of the Darc boards
/// the Global board input and the global board output
/// The structure containes the information of the 8 (at most) 
/// regional structures.
///
/// \author Christian Finck

/// \cond CLASSIMP
ClassImp(AliMUONDarcHeader)
/// \endcond

 const Int_t AliMUONDarcHeader::fgkDarcHeaderLength   =  1;
 const Int_t AliMUONDarcHeader::fgkGlobalHeaderLength =  5;
 const Int_t AliMUONDarcHeader::fgkDarcScalerLength   =  8;
 const Int_t AliMUONDarcHeader::fgkGlobalScalerLength =  10;

 const UInt_t AliMUONDarcHeader::fgkEndOfDarc   = 0xDEADFACE;
 const UInt_t AliMUONDarcHeader::fgkEndOfGlobal = 0xDEADBEEF;

//___________________________________________
AliMUONDarcHeader::AliMUONDarcHeader()
  :  TObject(),
     fWord(0),
     fGlobalOutput(0),

     fGlobalL0(0), 
     fGlobalClk(0),
     fGlobalHold(0),      
     fGlobalSpare(0),     

     fDarcL0R(0),
     fDarcL1P(0),
     fDarcL1S(0),
     fDarcL2A(0),
     fDarcL2R(0),
     fDarcClk(0),
     fDarcHold(0),
     fDarcSpare(0),
     fRegHeaderArray(new TClonesArray("AliMUONRegHeader",8))
  

{
  /// ctor
  
  for (Int_t i = 0; i < 4; i++)
    fGlobalInput[i] = 0;

  for (Int_t i = 0; i < 6; i++)
    fGlobalScaler[i] = 0;

}

//___________________________________________
AliMUONDarcHeader::AliMUONDarcHeader(const AliMUONDarcHeader& event)
  :  TObject(event),
     fWord(event.fWord),
     fGlobalOutput(event.fGlobalOutput),
     fGlobalL0(event.fGlobalL0),
     fGlobalClk(event.fGlobalClk),
     fGlobalHold(event.fGlobalHold),   
     fGlobalSpare(event.fGlobalSpare),

     fDarcL0R(event.fDarcL0R),
     fDarcL1P(event.fDarcL1P),
     fDarcL1S(event.fDarcL1S),
     fDarcL2A(event.fDarcL2A),
     fDarcL2R(event.fDarcL2R),
     fDarcClk(event.fDarcClk),
     fDarcHold(event.fDarcHold),
     fDarcSpare(event.fDarcSpare),
     fRegHeaderArray(new TClonesArray("AliMUONRegHeader", 8))

{
  ///
  /// copy ctor
  ///
 
 for (Int_t i = 0; i < 4; i++)
    fGlobalInput[i] = event.fGlobalInput[i];

  for (Int_t i = 0; i < 6; i++)
    fGlobalScaler[i] = event.fGlobalScaler[i];

  for (Int_t index = 0; index < (event.fRegHeaderArray)->GetEntriesFast(); index++) {
    new ((*fRegHeaderArray)[fRegHeaderArray->GetEntriesFast()]) 
        AliMUONRegHeader(*(AliMUONRegHeader*)(event.fRegHeaderArray)->At(index));
  }
}

//___________________________________________
AliMUONDarcHeader& AliMUONDarcHeader::operator=(const AliMUONDarcHeader& event)
{
  /// 
  /// assignment operator
  ///
  if (this == &event) return *this;

  fWord         = event.fWord;
  fGlobalOutput = event.fGlobalOutput;
  fGlobalL0     = event.fGlobalL0;
  fGlobalClk    = event.fGlobalClk;
  fGlobalHold   = event.fGlobalHold;   
  fGlobalSpare  = event.fGlobalSpare;

  fDarcL0R   = event.fDarcL0R;
  fDarcL1P   = event.fDarcL1P;
  fDarcL1S   = event.fDarcL1S;
  fDarcL2A   = event.fDarcL2A;
  fDarcL2R   = event.fDarcL2R;
  fDarcClk   = event.fDarcClk;
  fDarcHold  = event.fDarcHold;
  fDarcSpare = event.fDarcSpare;

  for (Int_t i = 0; i < 4; i++)
    fGlobalInput[i] = event.fGlobalInput[i];

  for (Int_t i = 0; i < 6; i++)
    fGlobalScaler[i] = event.fGlobalScaler[i];

  fRegHeaderArray = new TClonesArray("AliMUONRegHeader", 8);
  for (Int_t index = 0; index < (event.fRegHeaderArray)->GetEntriesFast(); index++) {
    new ((*fRegHeaderArray)[fRegHeaderArray->GetEntriesFast()]) 
        AliMUONRegHeader(*(AliMUONRegHeader*)(event.fRegHeaderArray)->At(index));
  }

  return *this;
}

//___________________________________________
AliMUONDarcHeader::~AliMUONDarcHeader()
{
  /// 
  /// dtor
  ///
  fRegHeaderArray->Delete();
  delete fRegHeaderArray;
}

//___________________________________________
void AliMUONDarcHeader::SetScalersNumbers()
{
  /// set numbers for scaler events for Darc header
  /// since this is provided by the experiment
  /// put dummy numbers to check the monitoring
  
  fGlobalL0    = 1000;
  fGlobalClk   = 10000;
  fGlobalHold  = 100;    
  fGlobalSpare = 1;    

  fDarcL0R   = 1000;
  fDarcL1P   = 900;
  fDarcL1S   = 800;
  fDarcL2A   = 700;
  fDarcL2R   = 700;
  fDarcClk   = 10000;
  fDarcHold  = 100;
  fDarcSpare = 0;

   for (Int_t i = 0; i < 6; i++)
    fGlobalScaler[i] = i;

}

//___________________________________________
void AliMUONDarcHeader::Clear(Option_t* )
{
  /// Clear TClones arrays
  /// instead of deleting
  ///
  fRegHeaderArray->Clear("C");
 
}
