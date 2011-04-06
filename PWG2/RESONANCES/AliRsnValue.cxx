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

////////////////////////////////////////////////////////////////////////////////
//
//  This class contains all code which is used to compute any of the values
//  which can be of interest within a resonance analysis. Besides the obvious
//  invariant mass, it allows to compute other utility values on all possible
//  targets, in order to allow a wide spectrum of binning and checks.
//  When needed, this object can also define a binning in the variable which
//  it is required to compute, which is used for initializing axes of output
//  histograms (see AliRsnFunction).
//  The value computation requires this object to be passed the object whose
//  informations will be used. This object can be of any allowed input type
//  (track, pair, event), then this class must inherit from AliRsnTarget.
//  Then, when value computation is attempted, a check on target type is done
//  and computation is successful only if expected target matches that of the
//  passed object.
//  In some cases, the value computation can require a support external object,
//  which must then be passed to this class. It can be of any type inheriting
//  from TObject.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliCentrality.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"
#include "AliRsnDaughterDef.h"

#include "AliRsnValue.h"

ClassImp(AliRsnValue)

//_____________________________________________________________________________
AliRsnValue::AliRsnValue() :
   AliRsnTarget(),
   fComputedValue(0),
   fBinArray(0)
{
//
// Default constructor without arguments.
// Initialize data members to meaningless values.
// This method is provided for ROOT streaming,
// but should never be used directly by a user.
//
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, Int_t nbins, Double_t min, Double_t max) :
   AliRsnTarget(name),
   fComputedValue(0.0),
   fBinArray(0)
{
//
// Main constructor (version 1).
// This constructor defines in meaningful way all data members,
// and defined a fixed binnings, subdividing the specified interval
// into that many bins as specified in the integer argument.
// ---
// This method is also the entry point for all instances
// of this class which don't need to do binning (e.g.: TNtuple inputs),
// since arguments 3 to 5 have default values which don't create any
// binning array, in order not to allocate memory when this is useless.
//

   SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, Double_t min, Double_t max, Double_t step) :
   AliRsnTarget(name),
   fComputedValue(0.0),
   fBinArray(0)
{
//
// Main constructor (version 2).
// This constructor defines in meaningful way all data members
// and creates enough equal bins of the specified size to cover
// the required interval.
//

   SetBins(min, max, step);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue
(const char *name, Int_t nbins, Double_t *array) :
   AliRsnTarget(name),
   fComputedValue(0.0),
   fBinArray(0)
{
//
// Main constructor (version 3).
// This constructor defines in meaningful way all data members
// and creates a set of variable bins delimited by the passed array.
//

   SetBins(nbins, array);
}

//_____________________________________________________________________________
AliRsnValue::AliRsnValue(const AliRsnValue& copy) :
   AliRsnTarget(copy),
   fComputedValue(copy.fComputedValue),
   fBinArray(copy.fBinArray)
{
//
// Copy constructor.
// Duplicates the binning array and copies all settings.
//
}

//_____________________________________________________________________________
AliRsnValue& AliRsnValue::operator=(const AliRsnValue& copy)
{
//
// Assignment operator.
// Works like copy constructor.
//

   AliRsnTarget::operator=(copy);

   fComputedValue = copy.fComputedValue;
   fBinArray = copy.fBinArray;

   return (*this);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t min, Double_t max)
{
//
// Set binning for the axis in equally spaced bins
// where the number of bins, minimum and maximum are given.
//

   if (!nbins) {
      fBinArray.Set(0);
      return;
   }

   fBinArray.Set(nbins + 1);

   Double_t mymax = TMath::Max(min, max);
   Double_t mymin = TMath::Min(min, max);

   Int_t    k = 0;
   Double_t binSize = (mymax - mymin) / ((Double_t)nbins);

   fBinArray[0] = mymin;
   for (k = 1; k <= nbins; k++) fBinArray[k] = fBinArray[k - 1] + binSize;
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Double_t min, Double_t max, Double_t step)
{
//
// Set binning for the axis in equally spaced bins
// where the bin size, minimum and maximum are given.
//

   Double_t dblNbins = TMath::Abs(max - min) / step;
   Int_t    intNbins = ((Int_t)dblNbins) + 1;

   SetBins(intNbins, min, max);
}

//_____________________________________________________________________________
void AliRsnValue::SetBins(Int_t nbins, Double_t *array)
{
//
// Set binning for the axis in unequally spaced bins
// using the same way it is done in TAxis
//

   if (!nbins) {
      fBinArray.Set(0);
      return;
   }

   Int_t i;
   fBinArray.Set(nbins);
   for (i = 0; i < nbins; i++) fBinArray[i] = array[i];
}

//_____________________________________________________________________________
Bool_t AliRsnValue::Eval(TObject *, Bool_t)
{
//
// Evaluation of the required value.
// Checks that the passed object is of the right type
// and if this check is successful, computes the required value.
// The output of the function tells if computing was successful,
// and the values must be taken with GetValue().
//

   AliWarning("This method must be overridden by derived classes");
   return kTRUE;
}

//_____________________________________________________________________________
void AliRsnValue::Print(Option_t *option) const
{
//
// Print informations about this object
//

   AliInfo("=== VALUE INFO =================================================");
   AliInfo(Form(" Name                  : %s", GetName()));
   AliInfo(Form(" Current computed value: %f", fComputedValue));
   if (!strcmp(option, "BINS")) {
      Int_t i;
      for (i = 0; i < fBinArray.GetSize(); i++) {
         AliInfo(Form(" Bin limit #%03d        = %f", i, fBinArray[i]));
      }
   }
}
