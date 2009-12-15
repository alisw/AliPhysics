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

/* $Id: AliTRDCalFEE.cxx 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD FEE parameters                             //
//  Empty dummy class for backward compability                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalFEE.h"

ClassImp(AliTRDCalFEE)

//_____________________________________________________________________________
AliTRDCalFEE::AliTRDCalFEE()
  :TNamed()
{
  //
  // AliTRDCalFEE default constructor
  //

}

//_____________________________________________________________________________
AliTRDCalFEE::AliTRDCalFEE(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
{
  //
  // AliTRDCalFEE constructor
  //

}
