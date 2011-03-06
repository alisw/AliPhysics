/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


//-------------------------------------------------------------------------
//     Realisation of AliVParticle for MC Particles
//     Basically a stripped down AliMCParicle / TParticle
//     with minimum information on MC tracks
//     Author: Christian Klein-Bösing, CERN
//-------------------------------------------------------------------------


#include "AliAODMCParticle.h"
#include "AliAODEvent.h"

#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TClonesArray.h"


ClassImp(AliAODMCParticle)

TString AliAODMCParticle::fgkStdBranchName("mcparticles");

AliAODMCParticle::AliAODMCParticle():
AliVParticle(),
  fPdgCode(0),
  fFlag(0),
  fLabel(0),
  fMother(0),
  fPx(0),
  fPy(0),
  fPz(0),
  fE(0),
  fVx(0),
  fVy(0),
  fVz(0),
  fVt(0)
{
  // Default Constructor
  fDaughter[0] =   fDaughter[1] = 0;
}

    
AliAODMCParticle::AliAODMCParticle(const AliMCParticle* mcpart, Int_t label,Int_t flag):
    AliVParticle(),
    fPdgCode(mcpart->Particle()->GetPdgCode()),
    fFlag(flag),
    fLabel(label),
    fMother(mcpart->GetMother()),
    fPx(mcpart->Particle()->Px()),
    fPy(mcpart->Particle()->Py()),
    fPz(mcpart->Particle()->Pz()),
    fE(mcpart->Particle()->Energy()),
    fVx(mcpart->Particle()->Vx()),
    fVy(mcpart->Particle()->Vy()),
    fVz(mcpart->Particle()->Vz()),
    fVt(mcpart->Particle()->T())
{
    // Constructor
  fDaughter[0] =  mcpart->GetFirstDaughter(); 
  fDaughter[1] =  mcpart->GetLastDaughter();
}
    
    
AliAODMCParticle::AliAODMCParticle(const AliAODMCParticle& mcPart) :
    AliVParticle(mcPart),
    fPdgCode(mcPart.fPdgCode),
    fFlag(mcPart.fFlag),
    fLabel(mcPart.fLabel),
    fMother(mcPart.fMother),
    fPx(mcPart.fPx),
    fPy(mcPart.fPy),
    fPz(mcPart.fPz),
    fE(mcPart.fE),
    fVx(mcPart.fVx),
    fVy(mcPart.fVy),
    fVz(mcPart.fVz),
    fVt(mcPart.fVt)
{
  // Copy constructor
  fDaughter[0] = mcPart.fDaughter[0]; 
  fDaughter[1] = mcPart.fDaughter[1]; 
}

AliAODMCParticle& AliAODMCParticle::operator=(const AliAODMCParticle& mcPart)
{ 

  if (this!=&mcPart) { 
    AliVParticle::operator=(mcPart);
    fPdgCode    = mcPart.fPdgCode;
    fFlag       = mcPart.fFlag;
    fLabel      = mcPart.fLabel;
    fMother     = mcPart.fMother;
    fPx         = mcPart.fPx;
    fPy         = mcPart.fPy;
    fPz         = mcPart.fPz;
    fE          = mcPart.fE;
    fVx         = mcPart.fVx;
    fVy         = mcPart.fVy;
    fVz         = mcPart.fVz;
    fVt         = mcPart.fVt;
    fDaughter[0] = mcPart.fDaughter[0]; 
    fDaughter[1] = mcPart.fDaughter[1]; 
  }  
  
  return *this;

}

Double_t AliAODMCParticle::M()         const
{
    TParticlePDG* pdg =  TDatabasePDG::Instance()->GetParticle(fPdgCode);
    if (pdg) {
	return (pdg->Mass());
    } else {
	return GetCalcMass();
    }
}


Short_t AliAODMCParticle::Charge()     const
{
    TParticlePDG* pdg =  TDatabasePDG::Instance()->GetParticle(fPdgCode);
    if (pdg) {
	return (Short_t (pdg->Charge()));
    } else {
	return -99;
    }
}

void AliAODMCParticle::Print(const Option_t */*opt*/) const {
// Print particle information
  if(TDatabasePDG::Instance()->GetParticle(fPdgCode)){
    Printf(">>> PDG (%d) : %s",fPdgCode,TDatabasePDG::Instance()->GetParticle(fPdgCode)->GetName());
  }
  else{
    Printf(">>> PDG (%d) : %s",fPdgCode,"Unknown");
  }
  Printf(">>  P(%3.3f,%3.3f,%3.3f) V((%3.3f,%3.3f,%3.3f,%3.3f)",fPx,fPy,fPz,fVx,fVy,fVz,fVt);  
  Printf(">   Mother %d, First Daughter %d Last Daughter %d , Status %d, PhysicalPrimary %d",
	 fMother,fDaughter[0],fDaughter[1],GetStatus(),
	 IsPhysicalPrimary());
}
