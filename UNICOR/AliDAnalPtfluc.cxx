// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

//=============================================================================
// pt-fluctuations analyzer
//=============================================================================

#include <TROOT.h>
#include <TRandom2.h>
#include <TVector2.h>
#include <TMath.h>
#include "AliDEvent.h"
#include "AliDHN.h"
#include "AliDAnalPtfluc.h"

ClassImp(AliDAnalPtfluc)

//=============================================================================
AliDAnalPtfluc::AliDAnalPtfluc(Char_t *nam, Int_t pid0, Int_t pid1) : 
  AliDAnal(nam), fPid0(pid0), fPid1(pid1) 
{
  // constructor

  TAxis *ax[5];
  ax[0] = new TAxis(2,-0.5,1.5);   ax[0]->SetTitle("trumix");
  ax[1] = new TAxis(9,0,0.9);      ax[1]->SetTitle("centrality");
  ax[2] = new TAxis(6,-0.5,5.5);   ax[2]->SetTitle("n-pt0-pt1-pt00-pt11-pt01");
  ax[3] = new TAxis(48,-180,180);  ax[3]->SetTitle("dphi (deg)");
  ax[4] = new TAxis(40,-2,2);      ax[4]->SetTitle("deta");
  AliDHN *pair = new AliDHN("pair",5,ax);
  for (int i=0; i<5; i++) delete ax[i];
  fHistos.Add(pair);
  gROOT->cd();
  printf("%s object named %s created\n",ClassName(),GetName());
}
//=============================================================================
void AliDAnalPtfluc::Process(Int_t tmr, AliDEvent *ev0, AliDEvent *ev1) 
{
  // process pairs from one or two (if mixing) events

  double ptmin=0.1;  // GeV
  double ptmax=1.5;  // GeV
  double etamin=-9;  
  double etamax=9;  

  // mixing-and-rotating-proof centrality

  double cent = (ev0->Centrality()+ev1->Centrality())/2.0;

  // loop over pairs 

  AliDHN *pair = (AliDHN*) fHistos.At(0);
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
      pair->Fill(tmr, cent, 0, dphi, deta, 1);   // number of pairs
      pair->Fill(tmr, cent, 1, dphi, deta, pt0);
      pair->Fill(tmr, cent, 2, dphi, deta, pt1);
      pair->Fill(tmr, cent, 3, dphi, deta, pt0*pt0);
      pair->Fill(tmr, cent, 4, dphi, deta, pt1*pt1);
      pair->Fill(tmr, cent, 5, dphi, deta, pt0*pt1);
    }
  }
}
//=============================================================================
