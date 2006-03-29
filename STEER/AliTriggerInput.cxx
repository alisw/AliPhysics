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

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define a Trigger Input from an specific detector                                                                                           //
//
//
//                        name         description     id mask
//    Ej:
//      AliTriggerInput( "V0_MB_L0", "VO minimum bias", 0x01 );
//      AliTriggerInput( "V0_SC_L0", "VO semi central", 0x02 );
//      AliTriggerInput( "V0_C_L0",  "VO central",      0x04 );
//
//    The name must be globaly unique. Spaces are not allowed.
//    As convention should start with detector name then an id
//    and the trigger level (L0, L1, L2)
//
//    A maximun of 60 inputs trigger are allow.
//    So, the id mask should set only bit from the position 1 to 60.
//
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliTriggerInput.h"

ClassImp( AliTriggerInput )

//_____________________________________________________________________________
void AliTriggerInput::Print( const Option_t* ) const
{
   // Print
   cout << "Trigger Input:" << endl; 
   cout << "  Name:        " << GetName() << endl;
   cout << "  Description: " << GetTitle() << endl;
   cout << "  Value:       " << hex << "Ox" << fValue << dec << endl;
}
