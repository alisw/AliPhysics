// $Id$
//--------------------------------------------------------------------------
// implementation of the AliTkEtaPhiVector class
//--------------------------------------------------------------------------

#include <Riostream.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "AliTkChargedJet.h"
#include "AliTkEtaPhiVector.h"

ostream& operator<<(ostream& s, const AliTkEtaPhiVector& v) {
  return s << "eta=" << v.Eta() << " phi=" << v.Phi();
}
