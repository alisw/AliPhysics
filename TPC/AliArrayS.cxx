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

/*
$Log$
Revision 1.1.4.2  2000/04/10 11:32:37  kowal2

"ROOT"-based class with some extra functionality

*/

///////////////////////////////////////////////////////////////////////
//   Added additional functionality  to original TArrayS              //
//   multiple inheritance from TObject to be possible use automatic   //
//   branch mechanism for tree
//   function Expand to be possible expand array without deleting     //
//   array contents                                                  //                  //
//                                                                   //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////
#include "AliArrayS.h"
ClassImp(AliArrayS)
void AliArrayS::Expand(Int_t n)
{
  //
  // Set array size of TArrayS object to n integers and copy old array
  // If n<0 leave array unchanged.
  // user are responsible for apropriate size of array
  // 
  if (n < 0) return;  
  fArray = (Short_t*)  TStorage::ReAlloc(fArray, n * sizeof(UShort_t),fN * sizeof(UShort_t));
  if (fArray!=0) fN= n; 
  else fN =0;
}



void AliArrayS::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTPC.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);     
      //read pad parameters
      R__b >> fN;
      if (fArray!=0){
	delete [] fArray;
	fArray =0;
      }
      if (fN>0){
	fArray = new Short_t[fN];
	R__b.ReadFastArray(fArray,fN); 
      }
   } else {
      R__b.WriteVersion(AliArrayS::IsA());
      TObject::Streamer(R__b);   
      R__b << fN;      
      if (fN>0) R__b.WriteFastArray(fArray,fN); 
   }
}
 
