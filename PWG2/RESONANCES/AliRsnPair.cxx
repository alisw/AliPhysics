//
// *** Class AliRsnPair ***
//
// "Core" method for defining the work on a pari of particles.
// For one analysis, one must setup one of this for each pair he wants to analyze,
// adding to it all analysis which he desires to do.
// Here he defines the cuts, and the particle types and charges, and can add
// functions which do different operations on the same pair, and some binning
// with respect to some kinematic variables (eta, momentum)
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <Riostream.h>
#include <TList.h>

#include "AliLog.h"

#include "AliRsnMother.h"
#include "AliRsnEvent.h"
#include "AliRsnFunction.h"
#include "AliRsnCutSet.h"
#include "AliRsnCutStd.h"
#include "AliRsnValue.h"
#include "AliRsnCutManager.h"

#include "AliRsnPair.h"

ClassImp(AliRsnPair)

//_____________________________________________________________________________
AliRsnPair::AliRsnPair(const char *name, AliRsnPairDef *def) :
  TNamed(name, ""),
  fOnlyTrue(kFALSE),
  fCheckDecay(kFALSE),
  fIsMixed(kFALSE),
  fPairDef(def),
  fCutManager(),
  fMother(),
  fEvent(0x0)
{
//
// Default constructor
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnPair::AliRsnPair(const AliRsnPair& copy) :
  TNamed(copy),
  fOnlyTrue(copy.fOnlyTrue),
  fCheckDecay(copy.fCheckDecay),
  fIsMixed(copy.fIsMixed),
  fPairDef(copy.fPairDef),
  fCutManager(copy.fCutManager),
  fMother(copy.fMother),
  fEvent(0x0)
{
//
// Default constructor
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnPair& AliRsnPair::operator=(const AliRsnPair& copy)
{
  fOnlyTrue = copy.fOnlyTrue;
  fCheckDecay = copy.fCheckDecay;
  fIsMixed = copy.fIsMixed;
  fPairDef = copy.fPairDef;
  fMother = copy.fMother;
  fCutManager = copy.fCutManager;
  fEvent = 0x0;

  return (*this);
}

//_____________________________________________________________________________
AliRsnPair::~AliRsnPair()
{
//
// Destructor
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnPair::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliInfo(Form("PDG %d %d", AliPID::ParticleCode(fPairDef->GetPID(0)), AliPID::ParticleCode(fPairDef->GetPID(1))));
  AliInfo(Form("Masses %f %f", fPairDef->GetMass(0), fPairDef->GetMass(1)));

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
Bool_t AliRsnPair::Fill
(AliRsnDaughter *daughter0, AliRsnDaughter *daughter1, AliRsnEvent *ev0, AliRsnEvent *ev1)
{
//
// Sets the two passed daughters to the AliRsnMother data member of this object
// which is used to perform all computations to fill the value list.
// This operation is done successfully only when the first passed object matches
// the required object type (track/V0) and the required charge for first element in pair def,
// and the second passed object does the same w.r. to the second element in pair def.
// Moreover, all cuts are checked and the operation fails if a cut check is unsuccessful.
// Finally, if a true pair is required, this is checked at the end.
//

  AliDebug(AliLog::kDebug+2,"<-");
  
  // first of all, compute the 4-momenta of the daughters
  // and that of the mother, according to current pair def
  // this could be needed for some cuts
  fMother.SetDaughters(daughter0, fPairDef->GetMass(0), daughter1, fPairDef->GetMass(1));
  
  // check for correct type-charge match for first element
  if (daughter0->RefType() != fPairDef->GetDaughterType(0)) return kFALSE;
  if (daughter0->ChargeChar() != fPairDef->GetCharge(0)) return kFALSE;
  
  // check for correct type-charge match for second element
  if (daughter1->RefType() != fPairDef->GetDaughterType(1)) return kFALSE;
  if (daughter1->ChargeChar() != fPairDef->GetCharge(1)) return kFALSE;
    
  // cuts on track #1 & common
  fCutManager.SetEvent(ev0);
  if (!fCutManager.PassDaughter1Cuts(daughter0)) 
  {
    AliDebug(AliLog::kDebug+2, "Specific cuts for track #1 not passed");
    return kFALSE;
  }
  if (!fCutManager.PassCommonDaughterCuts(daughter0))
  {
    AliDebug(AliLog::kDebug+2, "Common cuts for track #1 not passed");
    return kFALSE;
  }
  
  // cuts on track #2 & common
  fCutManager.SetEvent(ev1);
  if (!fCutManager.PassDaughter2Cuts(daughter1))
  {
    AliDebug(AliLog::kDebug+2, "Specific cuts for track #2 not passed");
    return kFALSE;
  }
  if (!fCutManager.PassCommonDaughterCuts(daughter1))
  {
    AliDebug(AliLog::kDebug+2, "Common cuts for track #2 not passed");
    return kFALSE;
  }
  
  // point to first event as reference
  fEvent = ev0;
  
  // define pair & check
  fCutManager.SetEvent(fEvent);
  if (!fCutManager.PassMotherCuts(&fMother)) return kFALSE;
  
  // if required a true pair, check this here and eventually return a fail message
  if (fOnlyTrue)
  {
    // are we in a MonteCarlo?
    if (!daughter0->GetParticle() || !daughter1->GetParticle()) return kFALSE;
    
    // are the daughters really secondaries (in MC)?
    Int_t m0 = daughter0->GetParticle()->GetFirstMother();
    Int_t m1 = daughter1->GetParticle()->GetFirstMother();
    if (m0 < 0 || m1 < 0) return kFALSE;
    
    // if they are, do they come from the same mother?
    if (m0 != m1) return kFALSE;
    
    //cout << "Checking a true pair..." << endl;
    
    // if they do, is this mother the correct type?
    Int_t mpdg0 = TMath::Abs(daughter0->GetMotherPDG());
    Int_t mpdg1 = TMath::Abs(daughter1->GetMotherPDG());
    Int_t mpdg  = TMath::Abs(fPairDef->GetMotherPDG());
    if (mpdg0 != mpdg) return kFALSE; //{cout << "Mother of d0 is " << mpdg0 << " instead of " << mpdg << endl; return kFALSE;}
    if (mpdg1 != mpdg) return kFALSE; //{cout << "Mother of d1 is " << mpdg1 << " instead of " << mpdg << endl; return kFALSE;}
    
    // do they match the expected decay channel (that is, are they the expected types)?
    if (fCheckDecay)
    {
      //cout << "Checking decay tree..." << endl;
      Int_t pdg0 = TMath::Abs(daughter0->GetParticle()->GetPdgCode());
      Int_t pdg1 = TMath::Abs(daughter1->GetParticle()->GetPdgCode());
      if (AliPID::ParticleCode(fPairDef->GetPID(0)) != pdg0) return kFALSE; // {cout << "PDG0 is " << pdg0 << " instead of " << fPairDef->GetPID(0) << endl; return kFALSE;};
      if (AliPID::ParticleCode(fPairDef->GetPID(1)) != pdg1) return kFALSE; // {cout << "PDG1 is " << pdg1 << " instead of " << fPairDef->GetPID(1) << endl; return kFALSE;};
      //cout << "Decay tree accepted!" << endl;
    }
    
    // ok... if we arrive here that must really be a true pair! :-)
    //cout << "Pair accepted!" << endl;
  }

  AliDebug(AliLog::kDebug+2,"->");
  
  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnPair::Compute()
{
//
// Virtual method to compute pair quantities of interest
//

  AliWarning("Implement this method in derived classes");
}

//_____________________________________________________________________________
void AliRsnPair::Init(const char* /*prefix*/, TList* /*list*/)
{
//
// Virtual method to compute pair quantities of interest
//

  AliWarning("Implement this method in derived classes");
}
