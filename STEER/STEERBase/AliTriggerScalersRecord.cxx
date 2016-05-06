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

/* $Id: AliTriggerScalersRecord.cxx 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//
// Class to define the ALICE Trigger Scalers Record 
//
// Each record consists of 1 time stamp (4 words)  (AliTimeStamp)
// and an array with the scalers (AliTriggerScalers) for each trigger class 
// in partition  
//
//////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TObjArray.h>
#include "AliLog.h"  
#include "AliTriggerScalers.h"
#include "AliTriggerScalersRecord.h"

using std::endl;
using std::cout;
ClassImp( AliTriggerScalersRecord )
//_____________________________________________________________________________
AliTriggerScalersRecord::AliTriggerScalersRecord():
  fTimestamp(),
  fScalers(),
  fTimeGroup(0)
{
 //Default constructor
}

//_____________________________________________________________________________
void AliTriggerScalersRecord::SetTimeStamp( UInt_t orbit, UInt_t period, 
                                            UInt_t seconds, UInt_t microsecs )
{
   fTimestamp.SetTimeStamp( orbit, period, seconds, microsecs );
}

//_____________________________________________________________________________
void AliTriggerScalersRecord::AddTriggerScalers( AliTriggerScalers* scaler ) 
{ 
  fScalers.AddLast( scaler ); 
  fScalers.Sort(); 
}

//_____________________________________________________________________________
void AliTriggerScalersRecord::AddTriggerScalers( UChar_t classIndex, UInt_t LOCB, UInt_t LOCA,        
                                         UInt_t L1CB, UInt_t L1CA, UInt_t L2CB, UInt_t L2CA )
{
    AddTriggerScalers( new AliTriggerScalers( classIndex, LOCB, LOCA, L1CB, L1CA, L2CB, L2CA ) );
    fScalers.Sort();
} 
//_____________________________________________________________________________
void AliTriggerScalersRecord::AddTriggerScalers( UChar_t classIndex, UInt_t LOCB, UInt_t LOCA,        
                                         UInt_t L1CB, UInt_t L1CA, UInt_t L2CB, UInt_t L2CA,
					 UInt_t LMCB, UInt_t LMCA)
{
    AddTriggerScalers( new AliTriggerScalers( classIndex, LOCB, LOCA, L1CB, L1CA, L2CB, L2CA, LMCB, LMCA ) );
    fScalers.Sort();
} 

//_____________________________________________________________________________
Int_t AliTriggerScalersRecord::Compare( const TObject* obj ) const
{
  // Compare  timestamps
  
  return fTimestamp.Compare( &(((AliTriggerScalersRecord*)obj)->fTimestamp) );
}
//_____________________________________________________________________________
const AliTriggerScalers* AliTriggerScalersRecord::GetTriggerScalersForClass( const Int_t classindex ) const
{
   // Find Trigger scaler with class ID = classindex using a brutal force 

   Int_t   position, last;
   AliTriggerScalers *op2 = 0;
   position = 0;
   last = fScalers.GetEntriesFast();
   while (position < last) {
      op2 = (AliTriggerScalers *)fScalers.At(position);
      if( op2 && (op2->GetClassIndex() == classindex )) break;
      op2=0;
      position++;
   }
   return op2;   
}

//_____________________________________________________________________________
AliTriggerScalers* AliTriggerScalersRecord::GetTriggerScalersForClassBinary( const Int_t classindex )
{
   // Find Trigger scaler with class ID = classindex using a binary search. 

   Int_t   base, position, last, result = 0;
   AliTriggerScalers *op2 = NULL;
   
   fScalers.Sort(); 
   
   base = 0;
   last = fScalers.GetEntriesFast();

   while (last >= base) {
      result = 0;
      position = (base+last) / 2;
      op2 = (AliTriggerScalers *)fScalers.At(position);
      if( op2 && op2->GetClassIndex() > classindex ) result = -1;
      if( op2 && op2->GetClassIndex() < classindex ) result = 1;
  
      if (op2 && result == 0)
         return op2;
      if (!op2 || result < 0)
         last = position-1;
      else
         base = position+1;
      op2 = NULL;   
   }
   return op2;   
}
                                      
//_____________________________________________________________________________
void AliTriggerScalersRecord::Print( const Option_t* ) const
{
   // Print
  cout << "Trigger Scalers Record, time group: "<< fTimeGroup << endl;
  fTimestamp.Print();
  for( Int_t i=0; i<fScalers.GetEntriesFast(); ++i ) 
     ((AliTriggerScalers*)fScalers.At(i))->Print();
}
