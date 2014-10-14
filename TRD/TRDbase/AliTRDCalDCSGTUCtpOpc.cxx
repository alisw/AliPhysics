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

/* $Id: AliTRDCalDCSGTUCtpOpc.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD DCS GTU parameters                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalDCSGTUCtpOpc.h"

ClassImp(AliTRDCalDCSGTUCtpOpc)

//_____________________________________________________________________________
AliTRDCalDCSGTUCtpOpc::AliTRDCalDCSGTUCtpOpc()
  :TNamed()  
    ,fId(0)
    ,fOpcode(0)
    ,fDirection(0)
    ,fInverted(0)
    ,fDelay(0)
    ,fConnected(0)
{
  //
  // AliTRDCalDCSGTU default constructor
  //



}

//_____________________________________________________________________________
AliTRDCalDCSGTUCtpOpc::AliTRDCalDCSGTUCtpOpc(const char *name, const char *title)
  :TNamed(name,title)
    ,fId(0)
    ,fOpcode(0)
    ,fDirection(0)
    ,fInverted(0)
    ,fDelay(0)
    ,fConnected(0)
{
  //
  // AliTRDCalDCSGTU constructor
  //



}


