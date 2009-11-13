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
// high-pt correlation analyzer
// Loop over pairs and fill pair histograms. The first particle is always 
// taken from ev0 and the second from ev1 so for ev0=ev1 one is getting true 
// pairs, otherwise mixed ones. Overlaps a bit with AliUnicorAnalCorrel. 
//=============================================================================

#include <TROOT.h>
#include <TMath.h>
#include <TVector2.h>
#include "AliUnicorEvent.h"
#include "AliUnicorHN.h"
#include "AliUnicorAnalHighpt.h"

ClassImp(AliUnicorAnalHighpt)
 
//=============================================================================
AliUnicorAnalHighpt::AliUnicorAnalHighpt(Char_t *nam, Double_t emi, Double_t ema, Int_t pid0, 
			 Int_t pid1): AliUnicorAnal(nam), fPid0(pid0), fPid1(pid1)
{
  // constructor
  // emi and ema define the rapidity range for histogram

  Double_t ewi = ema-emi; // width of the pseudorapidity range
  double pi = TMath::Pi();
  TAxis *ax[9];
  ax[0] = new TAxis(2,-0.5,1.5);   ax[0]->SetTitle("trumix");
  ax[1] = new TAxis(2,-0.5,1.5);   ax[1]->SetTitle("weigth"); // 1 or ass pt
  ax[2] = new TAxis(5,0,0.5);      ax[2]->SetTitle("centrality");
  ax[3] = new TAxis(3,emi,ema);    ax[3]->SetTitle("trig eta");
  ax[4] = new TAxis(1,-pi,pi);     ax[4]->SetTitle("trig phi"); // w.r.t. RP
  //  ax[4] = new TAxis(8,-pi,pi);     ax[4]->SetTitle("trig phi"); // w.r.t. RP
  double a5[]={0,2.5,3,4,6,8,10,15,20,30,50,100};
  ax[5] = new TAxis(11,a5);        ax[5]->SetTitle("trig pt (GeV)");
  ax[6] = new TAxis(20,-ewi,ewi);  ax[6]->SetTitle("ass eta - trig eta"); 
  ax[7] = new TAxis(36,-1,2*pi-1); ax[7]->SetTitle("ass phi - trig phi"); 
  ax[8] = new TAxis(10,0,1);       ax[8]->SetTitle("ass pt / trig pt");

  AliUnicorHN *pair = new AliUnicorHN("pair",9,ax);
  for (int i=0; i<9; i++) delete ax[i];
  fHistos.Add(pair);
  gROOT->cd();
  printf("%s object named %s created\n",ClassName(),GetName());
}
//=============================================================================
void AliUnicorAnalHighpt::Process(const AliUnicorEvent * const ev0, const AliUnicorEvent * const ev1) 
{
  // process pairs from one or two (if mixing) events
  // true pairs are within the triangle (j<i), mixed - all

  // mixing-proof centrality and reaction plane angle

  double cent = (ev0->Centrality()+ev1->Centrality())/2.0;
  double q0x,q0y,q1x,q1y;
  ev0->RP(q0x,q0y);
  ev1->RP(q1x,q1y); 
  double rpphi = atan2(q0y+q1y,q0x+q1x);

  // loop over pairs 

  int tmr = (ev0!=ev1); // 0 for true, 1 for mixed
  AliUnicorHN *pair = (AliUnicorHN*) fHistos.At(0);
  for (int i=0; i<ev0->NParticles(); i++) {
    if (!ev0->ParticleGood(i,fPid0)) continue;
    double eta0 = ev0->ParticleEta(i);
    double phi0 = ev0->ParticlePhi(i);
    double pt0 = ev0->ParticlePt(i);
    for (int j=0; j<ev1->NParticles(); j++) {
      if (ev0==ev1 && j==i) continue; // same particle
      if (ev0==ev1 && j<i && fPid0==fPid1) continue;
      if (!ev1->ParticleGood(j,fPid1)) continue;
      double eta1 = ev1->ParticleEta(j);
      double phi1 = ev1->ParticlePhi(j);
      double pt1 = ev1->ParticlePt(j);
      int order = (pt0>pt1); 
      double deta = order*(eta1-eta0); // ass-trig
      double dphi = order*(phi1-phi0);
      dphi = TVector2::Phi_mpi_pi(dphi);
      if (dphi<-1) dphi+=2*TMath::Pi(); 
      double relphi = TVector2::Phi_mpi_pi(phi0-rpphi);  // wrt. reaction plane
      pair->Fill((double)tmr,0.0,cent,eta0,relphi,pt0,deta,dphi,pt1/pt0,1.0); // number of pairs
      pair->Fill((double)tmr,1.0,cent,eta0,relphi,pt0,deta,dphi,pt1/pt0,pt1); // weigthed with ass pt
    }
  }
}
//=============================================================================
