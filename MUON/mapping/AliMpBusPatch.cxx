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

// $Id$
// $MpId: AliMpBusPatch.cxx,v 1.3 2006/03/17 11:51:34 ivana Exp $
// Category: management

// Class AliMpBusPatch
// ---------------
// Class that manages the maps buspatch<>DDL<>DE 
// for the mapping
// Calculates also the maximum DSP and buspatch numbers for a given DE
//
// Author: Ch. Finck; Subatech Nantes

#include "AliMpBusPatch.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"

#include "AliLog.h"

#include "TArrayI.h"
#include "Riostream.h"

ClassImp(AliMpBusPatch)

//////////////////////////////////////////////////////////
//
// This class contains the informations about buspatch vs (DE-DDL)
//
//////////////////////////////////////////////////////////


//_____________________________________________________________________________
AliMpBusPatch::AliMpBusPatch()
  : TObject(),
    fDetElemIdToBusPatch(300),
    fBusPatchToDetElem(300),
    fBusPatchToDDL(300)
{
/// Default constructor

  for (Int_t i = 0; i < 10; i++)
    fMaxBusPerCh[i] = 0;

}


//_____________________________________________________________________________
AliMpBusPatch::AliMpBusPatch(const AliMpBusPatch& rhs)
  : TObject(rhs)
{
/// Copy constructor

 *this = rhs;
}

//_____________________________________________________________________________
AliMpBusPatch::~AliMpBusPatch() 
{
/// Destructor

  fDetElemIdToBusPatch.Delete();
  fBusPatchToDetElem.Delete();
  fBusPatchToDDL.Delete();
}

//_____________________________________________________________________________
AliMpBusPatch& AliMpBusPatch::operator = (const AliMpBusPatch& /*rhs*/) 
{
/// Assignment operator
 
  AliFatal("= operator not implemented");

  return *this;
}

//____________________________________________________________________
Int_t AliMpBusPatch::GetDEfromBus(Int_t busPatchId)
{
 /// getting DE id from bus patch
  Long_t it = fBusPatchToDetElem.GetValue(busPatchId);

 if ( it ) 
   return (Int_t)it;
 else 
   return -1;
}

//____________________________________________________________________
TArrayI*  AliMpBusPatch::GetBusfromDE(Int_t idDE)
{
/// getting bus patch from DE id 

  return (TArrayI*)fDetElemIdToBusPatch.GetValue(idDE);
}
//____________________________________________________________________
Int_t AliMpBusPatch::GetDDLfromBus(Int_t busPatchId)
{
/// getting DE id from bus patch
  Long_t it = fBusPatchToDDL.GetValue(busPatchId);

 if ( it ) 
   return (Int_t)it;
 else 
   return -1;
}

//____________________________________________________________________
void AliMpBusPatch::GetDspInfo(Int_t iCh, Int_t& iDspMax, Int_t* iBusPerDSP) 
const
{
/// calculates the number of DSP & buspatch per block

  Int_t iBusPerBlk = fMaxBusPerCh[iCh]/4; //per half chamber; per block

  iDspMax =  iBusPerBlk/5; //number max of DSP per block
  if (iBusPerBlk % 5 != 0)
    iDspMax += 1;
  
  for (Int_t i = 0; i < iDspMax; i++) {
    if ((iBusPerBlk -= 5) > 0) 
      iBusPerDSP[i] = 5;
    else 
      iBusPerDSP[i] = iBusPerBlk + 5;
  }
  
}
//____________________________________________________________________
void AliMpBusPatch::ReadBusPatchFile()
{
/// idDE <> buspatch <> iDDL map's
  
   TString infile = AliMpFiles::BusPatchFilePath();

   ifstream in(infile, ios::in);
   if (!in) AliError("DetElemIdToBusPatch.dat not found.");
       
   char line[80];

   Int_t iChprev = 1;
   Int_t maxBusPatch = 0;

   while ( in.getline(line,80) ) {

      if ( line[0] == '#' ) continue;

      TString tmp(AliMpHelper::Normalize(line));

      Int_t blankPos  = tmp.First(' ');
      Int_t blankPos1 = tmp.Last(' ');

      TString sDE(tmp(0, blankPos));

      Int_t idDE = atoi(sDE.Data());
      
      if (idDE/100 != iChprev) {
	fMaxBusPerCh[iChprev-1] = maxBusPatch-iChprev*100+1;
	iChprev = idDE/100;
      }

      TString sDDL(tmp(blankPos1 + 1, tmp.Length()-blankPos1));

      Int_t iDDL = atoi(sDDL.Data());

      TString busPatch(tmp(blankPos + 1,blankPos1-blankPos-1));
      AliDebug(3,Form("idDE %d buspatch %s iDDL %d\n", idDE, busPatch.Data(), iDDL));

      TArrayI busPatchList;
      // decoding range of buspatch
      AliMpHelper::DecodeName(busPatch,';',busPatchList);
      
      // filling buspatch -> idDE
      for (Int_t i = 0; i < busPatchList.GetSize(); i++) {
	fBusPatchToDetElem.Add((Long_t)busPatchList[i],(Long_t)idDE);
	fBusPatchToDDL.Add((Long_t)busPatchList[i],(Long_t)iDDL);
	maxBusPatch = busPatchList[i];
      }
   
      // filling idDE -> buspatch list (vector)
      fDetElemIdToBusPatch.Add((Long_t)idDE, (Long_t)(new TArrayI(busPatchList))); 

    }
   
   fMaxBusPerCh[iChprev-1] = maxBusPatch-iChprev*100+1;

  in.close();

}
