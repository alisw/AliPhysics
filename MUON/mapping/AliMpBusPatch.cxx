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
// $MpId: AliMpBusPatch.cxx,v 1.5 2006/05/24 13:58:34 ivana Exp $
// Category: management

// Class AliMpBusPatch
// ---------------
// Class that manages the maps buspatch<>DDL<>DE 
// for the mapping
// Calculates also the maximum DSP and buspatch numbers for a given DDL
// Create a bus Iterator for DDL, needed especially for station 3
// Implementing a Sort method for Iterator, not really needed for the moment
//
// Author: Ch. Finck; Subatech Nantes

#include "AliMpBusPatch.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliLog.h"
#include "TArrayI.h"
#include "Riostream.h"

/// \cond CLASSIMP
ClassImp(AliMpBusPatch)
/// \endcond

const Int_t  AliMpBusPatch::fgkOffset = 100;

//_____________________________________________________________________________
AliMpBusPatch::AliMpBusPatch()
  : TObject(),
    fDetElemIdToBusPatch(300),
    fBusPatchToDetElem(300),
    fBusPatchToDDL(300)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpBusPatch::~AliMpBusPatch() 
{
/// Destructor

  fDetElemIdToBusPatch.Delete();
  fBusPatchToDetElem.Delete();
  fBusPatchToDDL.Delete();

}
//____________________________________________________________________
Int_t AliMpBusPatch::GetGlobalBusID(Int_t localID, Int_t ddlID)
{
  /// return global bus id from local bus and ddl id

  return ddlID*fgkOffset + localID;

}
//____________________________________________________________________
Int_t AliMpBusPatch::GetLocalBusID(Int_t globalID, Int_t ddlID)
{
  /// return local bus id from local bus id

  return globalID - ddlID*fgkOffset;

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
void AliMpBusPatch::GetDspInfo(Int_t iDDL, Int_t& iDspMax, Int_t* iBusPerDSP) 
const
{
/// calculates the number of DSP & buspatch per block

  Int_t iBusPerBlk = fBusInDDL[iDDL].GetSize()/2; //per block

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

   while ( in.getline(line,80) ) {

      if ( line[0] == '#' ) continue;

      TString tmp(AliMpHelper::Normalize(line));

      Int_t blankPos  = tmp.First(' ');
      Int_t blankPos1 = tmp.Last(' ');

      TString sDE(tmp(0, blankPos));

      Int_t idDE = atoi(sDE.Data());

      TString sDDL(tmp(blankPos1 + 1, tmp.Length()-blankPos1));

      Int_t iDDL = atoi(sDDL.Data());

      // always working local DDL number... for the moment.

      // not really needed remove for stand alone purpose (Ch.F)
      //      if (iDDL >= AliDAQ::DdlIDOffset("MUONTRK"))
      // 	iDDL -= AliDAQ::DdlIDOffset("MUONTRK");


      TString busPatch(tmp(blankPos + 1,blankPos1-blankPos-1));
      AliDebug(3,Form("idDE %d buspatch %s iDDL %d\n", idDE, busPatch.Data(), iDDL));
      AddDetElem(iDDL, idDE);

      TArrayI busPatchList;
      // decoding range of buspatch
      AliMpHelper::DecodeName(busPatch,';',busPatchList);
      
      // filling buspatch -> idDE

      for (Int_t i = 0; i < busPatchList.GetSize(); i++) {
	fBusPatchToDetElem.Add((Long_t)busPatchList[i],(Long_t)idDE);
	fBusPatchToDDL.Add((Long_t)busPatchList[i],(Long_t)iDDL);
	AddBus(iDDL, busPatchList[i]);
      }
   
      // filling idDE -> buspatch list (vector)
      fDetElemIdToBusPatch.Add((Long_t)idDE, (Long_t)(new TArrayI(busPatchList))); 

    }
   
  in.close();

}
//____________________________________________________________________
void AliMpBusPatch::AddBus(Int_t iDDL, Int_t busPatch)
{
/// add bus patch number per DDL

  fBusInDDL[iDDL].Set(fBusInDDL[iDDL].GetSize() + 1);
  fBusInDDL[iDDL].AddAt(busPatch, fBusInDDL[iDDL].GetSize() - 1);

}

//____________________________________________________________________
void AliMpBusPatch::AddDetElem(Int_t iDDL, Int_t detElem)
{
/// add DE per DDL

  fDeInDDL[iDDL].Set(fDeInDDL[iDDL].GetSize() + 1);
  fDeInDDL[iDDL].AddAt(detElem, fDeInDDL[iDDL].GetSize() - 1);

}
//____________________________________________________________________
Int_t AliMpBusPatch::NextBusInDDL(Int_t iDDL)
{
/// Next bus patch number in DDL

  if (fBusItr[iDDL] >= fBusInDDL[iDDL].GetSize())
    return -1;

  return fBusInDDL[iDDL].At(fBusItr[iDDL]++);

}


//____________________________________________________________________
void AliMpBusPatch::ResetBusItr(Int_t iDDL)
{
/// reset bus iterator for the given DDL

  fBusItr[iDDL] = 0;
}

//____________________________________________________________________
void AliMpBusPatch::Sort()
{
/// sort bus patch number for all DDL

  // put it hardware wise
  // this method is not used for the moment.
  Int_t numberOfDdls = 20;

  for (Int_t j = 0; j < numberOfDdls; j++) {
    Sort(fBusInDDL[j], 0, fBusInDDL[j].GetSize() - 1);

    if (AliLog::GetGlobalDebugLevel() == 1) {
      printf("DDL %d\n",j);
      for (Int_t i = 0; i <  fBusInDDL[j].GetSize(); i++)
	printf("buspatch %d index %d\n",fBusInDDL[j].At(i), i);
    }
  } 

}

//____________________________________________________________________
void AliMpBusPatch::Sort(TArrayI& arr, Int_t start, Int_t end)
{
/// sort bus patch number per DDL
/// not really needed, but for future developments ?
/// quicksort method, not in Root ?

  Int_t pivot;
  Int_t starth;
  Int_t endh; // store pivot # keep start & end in memory for split

  starth = start;
  endh = end;
  pivot = arr[start];

  while(start < end) {
      while((arr[end] >= pivot) && (start < end))
	end--;

      if (start != end) {
	  arr[start] = arr[end];
	  start++;
      }
      while ((arr[start] <= pivot) && (start < end))
	start++;

      if (start != end) {
	  arr[end] = arr[start];
	  end--;
      }
  }

  arr[start] = pivot;
  pivot = start;
  start = starth;
  end   = endh;

  if(start < pivot)
    Sort(arr, start, pivot-1);

  if(end > pivot)
    Sort(arr, pivot+1, end);

}
