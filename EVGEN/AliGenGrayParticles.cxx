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

/*
$Log$
Revision 1.3  2002/12/02 10:02:40  morsch
Corrections introduced by F. Silker:
- SetBetaSource
- Particle type according to charge.

Revision 1.2  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.1.2.1  2002/10/10 16:40:08  hristov
Updating VirtualMC to v3-09-02

Revision 1.1  2002/10/08 13:53:17  morsch
Gray particle generator, first commit.

*/

/*
  Generator for gray nucluons in pA interactions. 
  Source is modelled by a relativistic Maxwell distributions.
  Original code by  Ferenc Sikler  <sikler@rmki.kfki.hu>
 */
#include "AliGenGrayParticles.h"
#include "AliGrayParticleModel.h"
#include "AliPDG.h"
#include <TDatabasePDG.h>

 ClassImp(AliGenGrayParticles)
    
 AliGenGrayParticles::AliGenGrayParticles():AliGenerator(-1)
{
// Default constructor
}

AliGenGrayParticles::AliGenGrayParticles(Int_t npart)
    :AliGenerator(npart)
{
// Constructor
    fName  = "GrayParticles";
    fTitle = "Generator for gray particles in pA collisions";
    SetPmax();
    SetTarget();
    SetNominalCmsEnergy();
    SetCharge();
    SetTemperature();
    SetBetaSource();
}

//____________________________________________________________
AliGenGrayParticles::~AliGenGrayParticles()
{
// Destructor
}


void AliGenGrayParticles::Init()
{
  //
  // Initialization
  //
    Float_t kMass  = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
    fMomentum = fCMS/2. * fZTarget / fATarget;
    fBeta     = fMomentum / TMath::Sqrt(kMass * kMass + fMomentum * fMomentum);
}


void AliGenGrayParticles::Generate()
{
  //
  // Generate one event
  //
    Float_t p[3];
    Float_t origin[3] = {0., 0., 0.};
    Float_t polar [3] = {0., 0., 0.};    
    Int_t nt, i;
    for(i = 0;i < fNpart; i++) {
	Int_t kf;
        if(fCharge==1) kf = kProton;
                  else kf = kNeutron;
	GenerateSlow(fCharge, fTemperature, fBetaSource, p);
	
	SetTrack(fTrackIt, -1, kf, p, origin, polar,
		 0., kPNoProcess, nt, 1.);

	KeepTrack(nt);
    }
}




void AliGenGrayParticles::GenerateSlow(Int_t charge, Double_t T, Double_t beta, Float_t* q)
/* 
   Emit a slow nucleon with "temperature" T [GeV], 
   from a source moving with velocity beta         
   Three-momentum [GeV/c] is given back in q[3]    
*/

{
 Double_t m, pmax, p, f, theta, phi;
 
 TDatabasePDG * pdg = TDatabasePDG::Instance();
 const Double_t kMassProton  = pdg->GetParticle(kProton) ->Mass();
 const Double_t kMassNeutron = pdg->GetParticle(kNeutron)->Mass();
 
 /* Select nucleon type */
 if(charge == 0) m = kMassNeutron;
 else m = kMassProton;

 /* Momentum at maximum of Maxwell-distribution */

 pmax = TMath::Sqrt(2*T*(T+sqrt(T*T+m*m)));

 /* Try until proper momentum                                  */
 /* for lack of primitive function of the Maxwell-distribution */
 /* a brute force trial-accept loop, normalized at pmax        */

 do
 {
     p = Rndm() * fPmax;
     f = Maxwell(m, p, T) / Maxwell(m , pmax, T);
 }
 while(f < Rndm());

 /* Spherical symmetric emission */
 theta = TMath::ACos(2. * Rndm() - 1.);
 phi   = 2. * TMath::Pi() * Rndm();

 /* Determine momentum components in system of the moving source */
 q[0] = p * TMath::Sin(theta) * TMath::Cos(phi);
 q[1] = p * TMath::Sin(theta) * TMath::Sin(phi);
 q[2] = p * TMath::Cos(theta);

 /* Transform to system of the target nucleus                             */
 /* beta is passed as negative, because the gray nucleons are slowed down */
 Lorentz(m, -beta, q);

 /* Transform to laboratory system */
 Lorentz(m, fBeta, q);
}

Double_t AliGenGrayParticles::Maxwell(Double_t m, Double_t p, Double_t T)
{
/* Relativistic Maxwell-distribution */
    Double_t ekin;
    ekin = TMath::Sqrt(p*p+m*m)-m;
    return (p*p * exp(-ekin/T));
}


void AliGenGrayParticles::Lorentz(Double_t m, Double_t beta, Float_t* q)
{
/* Lorentz transform in the direction of q[2] */
 
    Double_t gamma  = 1/sqrt(1-beta*beta); 
    Double_t energy = sqrt(m*m + q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    q[2] = gamma * (q[2] + beta*energy);
}







