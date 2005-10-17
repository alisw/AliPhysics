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
//                                                                           //
//  Calibration base class for a single ROC                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalROCVdrift.h"

ClassImp(AliTRDCalROCVdrift)

//_____________________________________________________________________________
AliTRDCalROCVdrift::AliTRDCalROCVdrift():AliTRDCalROC()
{
  //
  // Default constructor
  //

  fNchannels    = 0;
  fVdrift       = 0;

}

//_____________________________________________________________________________
AliTRDCalROCVdrift::AliTRDCalROCVdrift(Int_t p, Int_t c)
                   :AliTRDCalROC(p,c)
{
  //
  // Constructor that initializes a given pad plane type
  //

  fNchannels = fNrows * fNcols;
  fVdrift    = new Float_t[fNchannels];

}

//_____________________________________________________________________________
AliTRDCalROCVdrift::AliTRDCalROCVdrift(const AliTRDCalROCVdrift &c)
                   :AliTRDCalROC(c)
{
  //
  // AliTRDCalROCVdrift copy constructor
  //

  ((AliTRDCalROCVdrift &) c).Copy(*this);

}

//_____________________________________________________________________________
AliTRDCalROCVdrift::~AliTRDCalROCVdrift()
{
  //
  // AliTRDCalROCVdrift destructor
  //

  if (fVdrift) {
    delete [] fVdrift;
    fVdrift = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalROCVdrift &AliTRDCalROCVdrift::operator=(const AliTRDCalROCVdrift &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalROCVdrift &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalROCVdrift::Copy(TObject &c) const
{
  //
  // Copy function
  //

  Int_t iBin = 0;

  ((AliTRDCalROCVdrift &) c).fNchannels = fNchannels;

  if (((AliTRDCalROCVdrift &) c).fVdrift) delete [] ((AliTRDCalROCVdrift &) c).fVdrift;
  ((AliTRDCalROCVdrift &) c).fVdrift = new Float_t[fNchannels];
  for (iBin = 0; iBin < fNchannels; iBin++) {
    ((AliTRDCalROCVdrift &) c).fVdrift[iBin] = fVdrift[iBin];
  }                                                                             

  AliTRDCalROC::Copy(c);

}

