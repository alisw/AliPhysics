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
//  This class is a simple set of definitions which are used to select among
//  all daughter candidates in term of the object type (track, V0, cascade),
//  charge and eventually PID.
//  The PID assigned to the definition is the one which is 'assigned' to the
//  candidate daughter for whatever implies a mass hypothesis. When MC is 
//  available and one requires this, thie PID is also used to select only the
//  daughters which are of that species in MC.
// 
//  NOTE: charge is a single character. If one wants to select a well-defined 
//        charge, he must use '+', '-' or '0', while any other character is
//        interpreted as 'no charge selection'.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliRsnDaughterDef.h"

ClassImp(AliRsnDaughterDef)


//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef() :
   fPID(AliRsnDaughter::kUnknown),
   fMass(0.0),
   fCharge(0),
   fRefType(AliRsnDaughter::kNoType)
{
//
// This version of constructor leaves everything undefined:
// this will cause all daughters to be accepted.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(AliRsnDaughter::ESpecies type, Char_t sign) :
   fPID(type),
   fMass(AliRsnDaughter::SpeciesMass(type)),
   fCharge(sign),
   fRefType(AliRsnDaughter::RefType(type))
{
//
// This version of constructor initializes the PID type
// and the charge (optional, leave 2nd argument to default to include both),
// and calls 'SetPID()' to assign the object type accordingly.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(EPARTYPE type, Char_t sign) :
   fPID(AliRsnDaughter::FromAliPID(type)),
   fMass(AliRsnDaughter::SpeciesMass(AliRsnDaughter::FromAliPID(type))),
   fCharge(sign),
   fRefType(AliRsnDaughter::RefType(AliRsnDaughter::FromAliPID(type)))
{
//
// This version of constructor initializes the PID type
// and the charge (optional, leave 2nd argument to default to include both),
// and calls 'SetPID()' to assign the object type accordingly.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(AliRsnDaughter::ERefType refType, Char_t sign) :
   fPID(AliRsnDaughter::kUnknown),
   fMass(0.0),
   fCharge(sign),
   fRefType(refType)
{
//
// This version of constructor initializes the object type
// and the charge (optional, leave 2nd argument to default to include both),
// and leaves the PID type undefined.
// This is useful when one is interested in all tracks/V0s/cascades without
// requiring them to be identified as a certain species, but if one then requires
// an object linked to this definition to compute a rapidity or a transverse mass,
// this will not work.
//
}

//_____________________________________________________________________________
AliRsnDaughterDef::AliRsnDaughterDef(const AliRsnDaughterDef &copy) :
   TObject(copy),
   fPID(copy.fPID),
   fMass(copy.fMass),
   fCharge(copy.fCharge),
   fRefType(copy.fRefType)
{
//
// Copy constructor has standard behavior.
//
}

//_____________________________________________________________________________
const AliRsnDaughterDef& AliRsnDaughterDef::operator=(const AliRsnDaughterDef &copy)
{
//
// Assignment operator has standard behavior.
//

   fMass = copy.fMass;
   fCharge = copy.fCharge;
   fPID = copy.fPID;
   fRefType = copy.fRefType;

   return (*this);
}

//_____________________________________________________________________________
Bool_t AliRsnDaughterDef::MatchesDaughter(AliRsnDaughter *checked, Bool_t truePID)
{
//
// Checks if the argument matches the definitions, by combining the other
// inline methods, and using the same philosophy.
// The only exception is for the PID matching, which can be disabled
// by second argument. In this case, a track is considered matched
// if it is matched just in object type and charge.
//

   Bool_t chargeMatch = MatchesCharge(checked);
   Bool_t objMatch    = MatchesRefType(checked);
   Bool_t pidMatch    = (truePID ? MatchesPID(checked) : kTRUE);
      
   // return the AND of all
   return (chargeMatch && objMatch && pidMatch);
}
