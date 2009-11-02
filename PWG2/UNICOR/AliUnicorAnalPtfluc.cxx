/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
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

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

//=============================================================================
// pt-fluctuations analyzer
// Loop over pairs and fill pair histograms. The first particle is always 
// taken from ev0 and the second from ev1 so for ev0=ev1 one is getting true 
// pairs, otherwise mixed ones. 
//=============================================================================

#include <TROOT.h>
#include <TRandom2.h>
#include <TVector2.h>
#include <TMath.h>
#include "AliUnicorEvent.h"
#include "AliUnicorHN.h"
#include "AliUnicorAnalPtfluc.h"

ClassImp(AliUnicorAnalPtfluc)

//=============================================================================
AliUnicorAnalPtfluc::AliUnicorAnalPtfluc(Char_t *nam, Int_t pid0, Int_t pid1) : 
  AliUnicorAnal(nam), fPid0(pid0), fPid1(pid1) 
{
  // constructor

  TAxis *ax[5];
  ax[0] = new TAxis(2,-0.5,1.5);   ax[0]->SetTitle("trumix");
  ax[1] = new TAxis(9,0,0.9);      ax[1]->SetTitle("centrality");
  ax[2] = new TAxis(6,-0.5,5.5);   ax[2]->SetTitle("n-pt0-pt1-pt00-pt11-pt01");
  ax[3] = new TAxis(48,-180,180);  ax[3]->SetTitle("dphi (deg)");
  ax[4] = new TAxis(40,-2,2);      ax[4]->SetTitle("deta");
  AliUnicorHN *pair = new AliUnicorHN("pair",5,ax);
  for (int i=0; i<5; i++) delete ax[i];
  fHistos.Add(pair);
  gROOT->cd();
  printf("%s object named %s created\n",ClassName(),GetName());
}
//=============================================================================
void AliUnicorAnalPtfluc::Process(Int_t tmr, AliUnicorEvent * const ev0, AliUnicorEvent * const ev1) 
{
  // process pairs from one or two (if mixing) events

  double ptmin=0.1;  // GeV
  double ptmax=1.5;  // GeV
  double etamin=-9;  
  double etamax=9;  

  // mixing-and-rotating-proof centrality

  double cent = (ev0->Centrality()+ev1->Centrality())/2.0;

  // loop over pairs 

  AliUnicorHN *pair = (AliUnicorHN*) fHistos.At(0);
  static TRandom2 ran;
  for (int i=0; i<ev0->NParticles(); i++) {
    if (!ev0->ParticleGood(i,fPid0)) continue;
    double eta0 = ev0->ParticleEta(i);
    double phi0 = ev0->ParticlePhi(i);
    double pt0 = ev0->ParticlePt(i);
    if (eta0 < etamin) continue;
    if (eta0 > etamax) continue;
    if (pt0 < ptmin) continue;
    if (pt0 > ptmax) continue;
    for (int j=0; j<ev1->NParticles(); j++) {
      if (ev0==ev1 && j==i) continue;
      if (ev0==ev1 && j<i && fPid0==fPid1) continue;
      if (!ev1->ParticleGood(j,fPid1)) continue;
      double eta1 = ev1->ParticleEta(j);
      double phi1 = ev1->ParticlePhi(j);
      double pt1 = ev1->ParticlePt(j);
      if (eta1 < etamin) continue;
      if (eta1 > etamax) continue;
      if (pt1 < ptmin) continue;
      if (pt1 > ptmax) continue;
      double deta = eta1-eta0;
      double dphi = phi1-phi0;
      // randomize order
      if (ran.Rndm()<0.5) { 
	double buf = pt0;
	pt0 = pt1;
	pt1 = buf;
	deta = -deta;
	dphi = -dphi;
      }
      dphi = TVector2::Phi_mpi_pi(dphi);
      dphi*=TMath::RadToDeg();
      pair->Fill((double) tmr, cent, 0.0, dphi, deta, 1.0);   // number of pairs
      pair->Fill((double) tmr, cent, 1.0, dphi, deta, pt0);
      pair->Fill((double) tmr, cent, 2.0, dphi, deta, pt1);
      pair->Fill((double) tmr, cent, 3.0, dphi, deta, pt0*pt0);
      pair->Fill((double) tmr, cent, 4.0, dphi, deta, pt1*pt1);
      pair->Fill((double) tmr, cent, 5.0, dphi, deta, pt0*pt1);
    }
  }
}
//=============================================================================
