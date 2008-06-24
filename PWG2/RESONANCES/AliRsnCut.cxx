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

//=========================================================================
// Class AliRsnCut
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()]
// - a value equal to a given reference     [--> MatchesValue()  ]
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//=========================================================================

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnMCInfo.h"
#include "AliRsnPairParticle.h"
#include "AliRsnPairDef.h"
#include "AliRsnCut.h"

const Double_t AliRsnCut::fgkDSmallNumber = 1e-100;
const Double_t AliRsnCut::fgkDBigNumber = 1e10;
const Int_t    AliRsnCut::fgkIBigNumber = 32767;

ClassImp (AliRsnCut)

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut() :
  TNamed(),
  fDMin(-fgkDBigNumber),
  fDMax( fgkDBigNumber),
  fIMin(-fgkIBigNumber),
  fIMax( fgkIBigNumber),
  fUIMin(0),
  fUIMax(2 * (UInt_t)fgkIBigNumber),
  fRsnCutType (kLastCutType),
  fRsnCutVarType (kDouble_t)
{
//
// Constructor
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut (const char *name, const char *title, ERsnCutType type) :
  TNamed (name,title),
  fDMin(-fgkDBigNumber),
  fDMax( fgkDBigNumber),
  fIMin(-fgkIBigNumber),
  fIMax( fgkIBigNumber),
  fUIMin(0),
  fUIMax(2 * (UInt_t)fgkIBigNumber),
  fRsnCutType (type),
  fRsnCutVarType (kDouble_t)
{
//
// Constructor with arguments but not limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut (const char *name, const char *title, ERsnCutType type, Double_t min, Double_t max) :
  TNamed (name,title),
  fDMin(min),
  fDMax(max),
  fIMin(-fgkIBigNumber),
  fIMax( fgkIBigNumber),
  fUIMin(0),
  fUIMax(2 * (UInt_t)fgkIBigNumber),
  fRsnCutType (type),
  fRsnCutVarType (kDouble_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut (const char * name, const char * title, ERsnCutType type, Int_t min, Int_t max) :
  TNamed (name,title),
  fDMin(-fgkDBigNumber),
  fDMax( fgkDBigNumber),
  fIMin(min),
  fIMax(max),
  fUIMin(0),
  fUIMax(2 * (UInt_t)fgkIBigNumber),
  fRsnCutType (type),
  fRsnCutVarType (kInt_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::AliRsnCut (const char * name, const char * title, ERsnCutType type, UInt_t min, UInt_t max) :
  TNamed (name,title),
  fDMin(-fgkDBigNumber),
  fDMax( fgkDBigNumber),
  fIMin(-fgkIBigNumber),
  fIMax( fgkIBigNumber),
  fUIMin(min),
  fUIMax(max),
  fRsnCutType (type),
  fRsnCutVarType (kUInt_t)
{
//
// Constructor with arguments and limits
//
}

//________________________________________________________________________________________________________________
AliRsnCut::~ AliRsnCut()
{
//
// Destructor.
// Does absolutely nothing.
//
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsBetween (const Double_t & theValue)
{
//
// Interval check.
// Question: "Is the argument included between fDMin and fDMax?"
// (not implemented for integer values because usually it is not used with them)
//
    return ((theValue >= fDMin) && (theValue <= fDMax));
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue (const Int_t &theValue)
{
//
// Reference check.
// Question: "Is the argument equal to fIMin?" (fIMax is assumed never used)
//
    return (theValue == fIMin);
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue (const UInt_t &theValue)
{
//
// Reference check.
// Question: "Is the argument equal to fUIMin?" (fUIMax is assumed never used)
//
    return (theValue == fUIMin);
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::MatchesValue (const Double_t &theValue)
{
//
// Reference check.
// Question: "Is the argument reasonably close to fDMin?" (fDMax is assumed never used)
// Here, "reasonably close" means that the difference is smaller than the
// 'fgkSmallNumber' global static data member of this class
//
    return (TMath::Abs (theValue - fDMin) < fgkDSmallNumber);
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues (ERsnCutType type, const Double_t & theValue, const Double_t & theValue2)
{
//
// (Re)assignment of cut values
//
    fRsnCutType = type;
    fDMin = theValue;
    fDMax = theValue2;
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues (ERsnCutType type, const Int_t& theValue, const Int_t& theValue2)
{
//
// (Re)assignment of cut values
//
    fRsnCutType = type;
    fIMin = theValue;
    fIMax = theValue2;
}

//________________________________________________________________________________________________________________
void AliRsnCut::SetCutValues (ERsnCutType type, const UInt_t& theValue, const UInt_t& theValue2)
{
//
// (Re)assignment of cut values
//
    fRsnCutType = type;
    fUIMin = theValue;
    fUIMax = theValue2;
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ECutSetType type, AliRsnDaughter *daughter)
{
//
// Core of the whole class.
// According to the kind of cut selected in the enumeration,
// checks the cut taking the right values from the argument.
// Depending on the second argument type, only some cuts are checked
// (the ones for that type of object), otherwise kTRUE is returned in order
// not to act as a cleaning factor for an AND with other cuts.
//
    AliDebug (AliLog::kDebug, "<-");
    AliRsnMCInfo *mcinfo = daughter->GetMCInfo();

    // check type
    if (type != kParticle) {
        AliWarning(Form("Mismatch: type = %d (expected %d), class type = %s (expected AliRsnDaughter)", type, kParticle, daughter->ClassName()));
        return kTRUE;
    }

    switch (fRsnCutType) {
        case kMomentum:
            return IsBetween (daughter->P());
        case kTransMomentum:
            return IsBetween (daughter->Pt());
        case kEta:
            return IsBetween (daughter->Eta());
        case kRadialImpactParam:
            return IsBetween (daughter->Vt());
        case kMomentumMC:
            if (mcinfo) return IsBetween (mcinfo->P());
            else return kTRUE;
        case kTransMomentumMC:
            if (mcinfo) return IsBetween (mcinfo->P());
            else return kTRUE;
        case kStatus:
            return daughter->CheckFlag(fUIMin);
        case kChargePos:
            return (daughter->Charge() > 0);
        case kChargeNeg:
            return (daughter->Charge() < 0);
        case kPIDType:
            return MatchesValue((Int_t)daughter->PIDType());
        /*
        case kEtaMC:
            if (mcinfo) return IsBetween (mcinfo->Eta());
            else return kTRUE;
        case kMcVt:
            if (mcinfo) return IsBetween (mcinfo->Vt());
            else return kTRUE;
        case kEsdNSigma:
            return IsBetween (daughter->GetNSigma());
        case kEsdNSigmaCalculate:
            return IsBetween (daughter->GetESDInfo()->GetNSigmaCalculate());
        */
        default:
            AliWarning("Requested a cut which cannot be applied to a single track");
            return kTRUE;
    }

    return kTRUE;
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::IsSelected(ECutSetType type, AliRsnPairParticle * pair)
{
    AliDebug (AliLog::kDebug, "<-");

    // check type
    if (type != kPair) {
        AliWarning(Form("Mismatch: type = %d (expected %d), class type = %s (expected AliRsnPairParticle)", type, kPair, pair->ClassName()));
        return kTRUE;
    }

    switch (fRsnCutType) {
        case kMomentum:
            return IsBetween (pair->GetP());
        case kTransMomentum:
            return IsBetween (pair->GetPt());
        /*
        case kEta:
            return IsBetween (daughter->Eta());
        */
        case kMomentumMC:
            return IsBetween (pair->GetPMC());
        case kTransMomentumMC:
            return IsBetween (pair->GetPtMC());
        case kRestMomentum:
            return CheckRestMomentum(pair);
        case kIsPdgEqual:
            return pair->IsPDGEqual();
        case kIsLabelEqual:
            return pair->IsLabelEqual();
        case kIsTruePair:
            return pair->IsTruePair(fIMin);
        default:
            AliWarning("Requested a cut which cannot be applied to a pair");
            return kTRUE;
    }

    return kTRUE;
}

//________________________________________________________________________________________________________________
void AliRsnCut::PrintAllValues()
{
  AliInfo (Form ("fRsnCutType=%d fRsnCutVarType=%d",fRsnCutType,fRsnCutVarType));
  AliInfo (Form ("fDMin=%.2e fDMax=%.2e",fDMin,fDMax));
  AliInfo (Form ("fIMin=%d fIMax=%d",fIMin,fIMax));
  AliInfo (Form ("fUIMin=%d fUIMax=%d",fUIMin,fUIMax));
}

//________________________________________________________________________________________________________________
Bool_t AliRsnCut::CheckRestMomentum(AliRsnPairParticle *pair)
{
//
// Check the cut on daughter momenta in rest reference frame of mother
//

    Double_t beta = pair->GetP() / pair->GetMass();
    Double_t gamma = 1. / TMath::Sqrt(1. - beta*beta);

    Double_t p1labP = 0.0, p1labT, p1restP, p1restTot;
    p1labP += pair->GetDaughter(0)->Px() * pair->GetP(0);
    p1labP += pair->GetDaughter(0)->Py() * pair->GetP(1);
    p1labP += pair->GetDaughter(0)->Pz() * pair->GetP(2);
    p1labP /= pair->GetP();

    p1labT = TMath::Sqrt(pair->GetDaughter(0)->P2() - p1labP*p1labP);

    p1restP = gamma*p1labP - beta*gamma*pair->GetDaughter(0)->E();

    p1restTot = TMath::Sqrt(p1restP*p1restP + p1labT*p1labT);

    return IsBetween(p1restTot);
}



