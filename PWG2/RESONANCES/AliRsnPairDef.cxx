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
//  This class is a simple set of definitions which are used to define a
//  decay tree to be studied for a resonance, in terms of the PID and charge
//  of its candidate daughters, which in turn determins what kind of objects
//  the analysis must take into account.
//  This object contains two AliRsnDaughterDef which define a model for each
//  of the two expected daughters (see also AliRsnDaughterDef class) plus a
//  mass hypothesis for the resonance, which is used for computin quantities
//  which need it (like rapidity or Mt), and a PDG code, which is used to 
//  check for true pairs, when needed. In all other cases, these two additional
//  values can be left to their default (meaningless) value.
//  Since this object must define a decay channel, the only provided constructor
//  allow to set a PID and a charge.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliRsnMother.h"
#include "AliRsnPairDef.h"

ClassImp(AliRsnPairDef)

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef() :
   fMotherMass(0.0), 
   fMotherPDG(0), 
   fDef1(), 
   fDef2()
{
//
// Constructor.
// If the two pointers are well initialized, they are used to init the members.
//
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef
(EPARTYPE type1, Char_t ch1, EPARTYPE type2, Char_t ch2, Int_t pdg, Double_t mass) :
   fMotherMass(mass), 
   fMotherPDG(pdg), 
   fDef1(type1, ch1), 
   fDef2(type2, ch2)
{
//
// Constructor.
// If the two pointers are well initialized, they are used to init the members.
//
}

//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef
(AliRsnDaughter::ESpecies type1, Char_t ch1, AliRsnDaughter::ESpecies type2, Char_t ch2, Int_t pdg, Double_t mass) :
   fMotherMass(mass), 
   fMotherPDG(pdg), 
   fDef1(type1, ch1), 
   fDef2(type2, ch2)
{
//
// Constructor.
// If the two pointers are well initialized, they are used to init the members.
//
}
   
//_____________________________________________________________________________
AliRsnPairDef::AliRsnPairDef(const AliRsnPairDef &copy) :
   TObject(copy),
   fMotherMass(copy.fMotherMass),
   fMotherPDG(copy.fMotherPDG),
   fDef1(copy.fDef1),
   fDef2(copy.fDef2)
{
//
// Copy constructor with standard behavior
//
}

//_____________________________________________________________________________
const AliRsnPairDef& AliRsnPairDef::operator=(const AliRsnPairDef &copy)
{
//
// Assignment operator with standard behavior.
//

   fMotherMass = copy.fMotherMass;
   fMotherPDG = copy.fMotherPDG;
   fDef1 = copy.fDef1;
   fDef2 = copy.fDef2;

   return (*this);
}
