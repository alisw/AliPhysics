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

/* $Id: AliTriggerScalers.cxx 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define the ALICE Trigger Scalers  
//
//  For each trigger class there are six scalers:
//
//    LOCB       L0 triggers before any vetos 
//    LOCA       L0 triggers after all vetos 
//    L1CB       L1 triggers before any vetos 
//    L1CA       L1 triggers after all vetos 
//    L2CB       L2 triggers before any vetos 
//    L2CA       L2 triggers after all vetos 
//    LMCB       LM triggers before any vetos 
//    LMCA       LM triggers after all vetos 
//
//////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliLog.h"
#include "AliTriggerScalers.h"

using std::endl;
using std::cout;
ClassImp( AliTriggerScalers )

//_____________________________________________________________________________
AliTriggerScalers::AliTriggerScalers(): 
  TObject(),
  fClassIndex(0),
  fLOCB(0),     
  fLOCA(0),     
  fL1CB(0),     
  fL1CA(0),     
  fL2CB(0),     
  fL2CA(0),      
  fLMCB(0),     
  fLMCA(0)      
{
  // Default constructor
}

//_____________________________________________________________________________
AliTriggerScalers::AliTriggerScalers( UChar_t classIndex, UInt_t LOCB, UInt_t LOCA,        
                                      UInt_t L1CB, UInt_t L1CA, UInt_t L2CB, UInt_t L2CA ):   
  TObject(),
  fClassIndex( classIndex ),
  fLOCB(LOCB),     
  fLOCA(LOCA),     
  fL1CB(L1CB),     
  fL1CA(L1CA),     
  fL2CB(L2CB),     
  fL2CA(L2CA),      
  fLMCB(0),     
  fLMCA(0)      
{
  // Default constructor for L0,L1,L2 levels 
}
//_____________________________________________________________________________
AliTriggerScalers::AliTriggerScalers( UChar_t classIndex, UInt_t LOCB, UInt_t LOCA,        
                                      UInt_t L1CB, UInt_t L1CA, UInt_t L2CB, UInt_t L2CA, UInt_t LMCB,UInt_t LMCA ):   
  TObject(),
  fClassIndex( classIndex ),
  fLOCB(LOCB),     
  fLOCA(LOCA),     
  fL1CB(L1CB),     
  fL1CA(L1CA),     
  fL2CB(L2CB),     
  fL2CA(L2CA),      
  fLMCB(LMCB),     
  fLMCA(LMCA)      
{
  // Constructor for L0,L1,L2,LM levels 
}

//_____________________________________________________________________________
Int_t  AliTriggerScalers::Compare( const TObject* obj ) const
{
  // Compare Scaler by class index (to sort in Scaler Record by class index)
  if( fClassIndex < ((AliTriggerScalers*)obj)->fClassIndex ) return -1;
  if( fClassIndex > ((AliTriggerScalers*)obj)->fClassIndex ) return 1;
  return 0;
}
//____________________________________________________________________________
void AliTriggerScalers::GetAllScalers(UInt_t *scalers) const
{
 scalers[0]=fLOCB;
 scalers[1]=fLOCA;
 scalers[2]=fL1CB;
 scalers[3]=fL1CA;
 scalers[4]=fL2CB;
 scalers[5]=fL2CA;
}
//____________________________________________________________________________
void AliTriggerScalers::GetAllScalersM012(UInt_t *scalers) const
{
 scalers[2]=fLOCB;
 scalers[3]=fLOCA;
 scalers[4]=fL1CB;
 scalers[5]=fL1CA;
 scalers[6]=fL2CB;
 scalers[7]=fL2CA;
 scalers[0]=fLMCB;
 scalers[1]=fLMCA;
}
//_____________________________________________________________________________
void AliTriggerScalers::Print( const Option_t* ) const
{
   // Print
  cout << "Trigger Scalers for Class: " << (Int_t)fClassIndex << endl;
  cout << "  LOCB: " << fLOCB << " LOCA: " << fLOCA; //<< endl;
  cout << "  L1CB: " << fL1CB << " L1CA: " << fL1CA; //<< endl;
  cout << "  L2CB: " << fL2CB << " L2CA: " << fL2CA << endl;
  cout << "  LMCB: " << fLMCB << " LMCA: " << fLMCA << endl;
}
