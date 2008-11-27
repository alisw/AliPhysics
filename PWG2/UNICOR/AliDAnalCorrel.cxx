// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2005

//=============================================================================
// two-particle correlation analyzer
//=============================================================================

#include <TROOT.h>
#include <TMath.h>
#include <TRandom2.h>
#include "AliDEvent.h"
#include "AliDHN.h"
#include "AliDAnalCorrel.h"

ClassImp(AliDAnalCorrel)
 
//=============================================================================
AliDAnalCorrel::AliDAnalCorrel(Char_t *nam, Double_t emi, Double_t ema, 
			 Int_t pid0, Int_t pid1): 
  AliDAnal(nam), fPid0(pid0), fPid1(pid1), fMass0(0), fMass1(0), fPa()
{
  // constructor
  // emi and ema define the rapidity range for histogram

  TParticlePDG *part0 = AliDAnal::fgPDG.GetParticle(fPid0);
  TParticlePDG *part1 = AliDAnal::fgPDG.GetParticle(fPid1);
  fMass0 = part0? part0->Mass() : 0;
  fMass1 = part1? part1->Mass() : 0;

  double pi = TMath::Pi();
  TAxis *ax[8];
  ax[0] = new TAxis(3,-0.5,2.5);ax[0]->SetTitle("trumixrot");
  ax[1] = new TAxis(5,0,0.5);   ax[1]->SetTitle("centrality");
  ax[2] = new TAxis(3,emi,ema); ax[2]->SetTitle("pair y");
  ax[3] = new TAxis(8,-pi,pi);  ax[3]->SetTitle("pair phi"); // wrt event plane
  double a0[]={0,0.1,0.2,0.3,0.4,0.5,0.7,1.0};
  ax[4] = new TAxis(7,a0);      ax[4]->SetTitle("(pair pt)/2 (GeV)");
  ax[5] = new TAxis(8,0,pi);    ax[5]->SetTitle("q-theta");
  ax[6] = new TAxis(16,-pi,pi); ax[6]->SetTitle("q-phi");
  double a1[100];
  for (int i=0;i<20;i++) a1[i]=i*0.005;
  for (int i=0;i<45;i++) a1[20+i]=0.1+i*0.02;
  ax[7] = new TAxis(64,a1);     ax[7]->SetTitle("q (GeV/c)");
  AliDHN *pair = new AliDHN("pair",8,ax);
  for (int i=0; i<8; i++) delete ax[i];
  fHistos.Add(pair);
  gROOT->cd();
  printf("%s object named %s created\n",ClassName(),GetName());
}
//=============================================================================
void AliDAnalCorrel::Process(Int_t tmr, AliDEvent *ev0, AliDEvent *ev1, Double_t phirot) 
{
  // process pairs from one or two (if mixing) events
  // tmr tells which histogram (bins) to fill: tru,mix,rot

  static TRandom2 ran;
  AliDHN *pair = (AliDHN*) fHistos.At(0);

  // mixing-and-rotating-proof centrality and reaction plane angle
  // (but not rotation-proof for rotation angles much different from 0 and 180)
  // true and rotated pairs are within the triangle (j<i), mixed - all
  // thus, proper rotation is either by 180, or by 170 AND 190, etc. 

  double cent = (ev0->Centrality()+ev1->Centrality())/2.0;
  double q0x,q0y,q1x,q1y;
  ev0->RP(q0x,q0y);
  ev1->RP(q1x,q1y); 
  double rpphi = atan2(q0y+q1y,q0x+q1x);

  // loop over pairs 

  for (int i=0; i<ev0->NParticles(); i++) {
    if (!ev0->ParticleGood(i,fPid0)) continue;
    for (int j=0; j<ev1->NParticles(); j++) {
      if (ev0==ev1 && j<i && fPid0==fPid1 ) continue; 
      if (ev0==ev1 && j==i) continue; // beware, not even when rotated or non-identical
      if (!ev1->ParticleGood(j,fPid1)) continue;
      if (!ev0->PairGood(ev0->ParticleP(i),ev0->ParticleTheta(i),ev0->ParticlePhi(i),
			 ev1->ParticleP(j),ev1->ParticleTheta(j),ev1->ParticlePhi(j)+phirot)) continue;
      fPa.Set0(fMass0,ev0->ParticleP(i),ev0->ParticleTheta(i),ev0->ParticlePhi(i));
      fPa.Set1(fMass1,ev1->ParticleP(j),ev1->ParticleTheta(j),ev1->ParticlePhi(j)+phirot);
      if (ev0==ev1 && fPid0==fPid1 && ran.Rndm()>=0.5) fPa.Swap();
      fPa.CalcLAB();
      fPa.CalcPairCM();
      if (fPa.QCM()==0) return; // should not be too frequent
      double phi = TVector2::Phi_mpi_pi(fPa.Phi()-rpphi);
      pair->Fill(tmr,                    // 0 for tru, 1 for mix, 2 for rot
		 cent,                   // centrality
		 fPa.Rapidity(),         // pair rapidity
		 phi,                    // pair phi wrt reaction plane
		 fPa.Pt()/2.0,           // half of pair pt
		 fPa.QCMTheta(),         // polar angle of Q
		 fPa.QCMPhiOut(),        // azimuthal angle of Q w.r.t. out
		 fPa.QCM(),              // |p2-p1| in c.m.s.
		 1.0);                   // weigth
    }
  }
}
//=============================================================================
