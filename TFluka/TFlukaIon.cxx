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

/* $Id$*/

//
// Stores user defined ion properties. 
// Fluka allows only one user defined ion which can be used as a beam particle.
// Author:
// A. Morsch 
// andreas.morsch@cern.ch
//

#include "TFlukaIon.h"
#include <TDatabasePDG.h>
   
ClassImp(TFlukaIon)


TFlukaIon::TFlukaIon() : 
    TNamed("", "Ion"),
    fZ(0), 
    fA(0),
    fQ(0),
    fExEnergy(0.),
    fMass(0.)
{
// Default constructor
}


TFlukaIon::TFlukaIon(const char* name, Int_t z, Int_t a, Int_t q, Double_t exE, Double_t mass) :
    TNamed(name, "Ion"),
    fZ(z), 
    fA(a),
    fQ(q),
    fExEnergy(exE),
    fMass(mass)

{
// Constructor
    AddIon(a, z);
}

Int_t TFlukaIon::GetIonPdg(Int_t z, Int_t a, Int_t i)
{
// Acording to
// http://cepa.fnal.gov/psm/stdhep/pdg/montecarlorpp-2006.pdf

  return 1000000000 + 10*1000*z + 10*a + i;
}  

Int_t TFlukaIon::GetZ(Int_t pdg)
{
// Acording to
// http://cepa.fnal.gov/psm/stdhep/pdg/montecarlorpp-2006.pdf

  return (pdg - 1000000000)/10000; 
}  


Int_t TFlukaIon::GetA(Int_t pdg)
{
// Acording to
// http://cepa.fnal.gov/psm/stdhep/pdg/montecarlorpp-2006.pdf

    Int_t a = pdg - 1000000000;
    a %= 10000;
    a /= 10;
    return (a);
}  

Int_t TFlukaIon::GetIsomerNumber(Int_t pdg)
{
// Acording to
// http://cepa.fnal.gov/psm/stdhep/pdg/montecarlorpp-2006.pdf

    Int_t is = pdg - 1000000000;
    is %= 10000;
    is %= 10;
    return (is);
}  

void TFlukaIon::AddIon(Int_t a, Int_t z)
{

    // Add a new ion
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    const Double_t kAu2Gev   = 0.9314943228;
    Int_t pdg =  GetIonPdg(z, a);
    if (pdgDB->GetParticle(pdg)) return;
    
    pdgDB->AddParticle(Form("Iion A  = %5d Z = %5d", a, z),"Ion", Float_t(a) * kAu2Gev + 8.071e-3, kTRUE,
		       0, 3 * z, "Ion", pdg);
}

void  TFlukaIon::AddIon(const char* name, Int_t z, Int_t a, Int_t q,
			Double_t /*exE*/, Double_t mass)
{
// User defined ion
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    const Double_t kAu2Gev   = 0.9314943228;
    Int_t pdg =  GetIonPdg(z, a, q);
    if (pdgDB->GetParticle(pdg)) return;
    if (mass = 0.) mass = Float_t(a) * kAu2Gev + 8.071e-3;

    pdgDB->AddParticle(name, "User Ion", mass, kTRUE, 0, 3 * z, "Ion", pdg);
}

void TFlukaIon::WriteUserInputCard(FILE* pFlukaVmcInp)
{
    // Write the user input card
    // EVENTYPE          0.        0.        2.        0.        0.        0.DPMJET    
    fprintf(pFlukaVmcInp,"EVENTYPE          0.        0.        2.        0.        0.        0.DPMJET\n");
}
