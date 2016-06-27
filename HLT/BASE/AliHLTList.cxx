/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Mikolaj Krzewicki                                              *
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

/* $Id: AliHLTList.cxx 2016-08-25 13:09:59Z mkrzewic $ */

//-----------------------------------------------------------------
//           Implementation of the AliHLTList class
//   TList that becomes owner of its contents after streaming
//   so it is safer to stream
//
// Origin: Mikolaj Krzewicki, mkrzewic@cern.ch
//-----------------------------------------------------------------

#include "AliHLTList.h"
#include "TObject.h"
#include "TBuffer.h"
#include <string.h>

using namespace std;

//_______________________________________________________________________
void AliHLTList::Streamer(TBuffer &b)
{
  // Stream all objects in the array to or from the I/O buffer.
  if (b.IsReading()) {
    b.ReadClassBuffer(AliHLTList::Class(),this);
    this->SetOwner(kTRUE);
  } else {
    b.WriteClassBuffer(AliHLTList::Class(),this);
  }
}

