//
// Class AliRsnCutPIDNSigma
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnPIDRange.h"

ClassImp(AliRsnPIDRange)

//_____________________________________________________________________________
AliRsnPIDRange::AliRsnPIDRange(Double_t nsigma, Double_t pmin, Double_t pmax) : TObject() ,
   fPMin(pmin),
   fPMax(pmax),
   fNSigmaCut(nsigma)
{
//
// Default constructor
//

}

//_____________________________________________________________________________
AliRsnPIDRange::AliRsnPIDRange(const AliRsnPIDRange &copy) : TObject(copy),
   fPMin(copy.fPMin),
   fPMax(copy.fPMax),
   fNSigmaCut(copy.fNSigmaCut)
{
//
// Copy constructor
//

}

//_____________________________________________________________________________
AliRsnPIDRange &AliRsnPIDRange::operator=(const AliRsnPIDRange &copy)
{
//
// Assignment operator.
//

   TObject::operator=(copy);
   if (this == &copy)
      return *this;
   fPMin = copy.fPMin;
   fPMax = copy.fPMax;
   fNSigmaCut = copy.fNSigmaCut;

   return (*this);
}
