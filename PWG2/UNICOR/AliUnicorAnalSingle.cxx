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

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// single particle analyzer
// Loop over tracks of one event and, for tracks that fulfill the quality 
// and pid cuts, fill single track histograms. 
//=============================================================================

#include <cmath>
#include <TROOT.h>
#include <TMath.h>
#include <TAxis.h>
#include <TParticlePDG.h>
#include "AliUnicorHN.h"
#include "AliUnicorEvent.h"
#include "AliUnicorAnalSingle.h"

ClassImp(AliUnicorAnalSingle)

//=============================================================================
AliUnicorAnalSingle::AliUnicorAnalSingle(Char_t *nam, Double_t emi, Double_t ema, Int_t pid) : 
  AliUnicorAnal(nam), fPid(pid), fMass(0.0) 
{
  // constructor
  // emi and ema define the rapidity range for histograms

  fPid = pid;
  TParticlePDG *part = AliUnicorAnal::fgPDG.GetParticle(fPid);
  fMass = part? part->Mass() : 0;

  double pi = TMath::Pi();
  TAxis *ax[10];
  ax[0] = new TAxis(30,-1,1);    ax[0]->SetTitle("vertex z");
  ax[1] = new TAxis(80,emi,ema); ax[1]->SetTitle("eta");
  ax[2] = new TAxis(90,-pi,pi);  ax[2]->SetTitle("phi");
  AliUnicorHN *zep = new AliUnicorHN("zep",3,ax);
  for (int i=0; i<3; i++) delete ax[i];

  ax[0] = new TAxis(20,0,1);     ax[0]->SetTitle("centrality");
  ax[1] = new TAxis(80,emi,ema); ax[1]->SetTitle("y");
  ax[2] = new TAxis(80,0,2);     ax[2]->SetTitle("pt (GeV)");
  AliUnicorHN *cyp = new AliUnicorHN("cyp",3,ax);
  for (int i=0; i<3; i++) delete ax[i];

  ax[0] = new TAxis(10,emi,ema); ax[0]->SetTitle("eta");
  ax[1] = new TAxis(150,0,3);    ax[1]->SetTitle("p (GeV)");
  ax[2] = new TAxis(150,0.5,3.5);ax[2]->SetTitle("sqrt(dedx (mips))");
  AliUnicorHN *epd = new AliUnicorHN("epd",3,ax);
  for (int i=0; i<3; i++) delete ax[i];

  fHistos.Add(zep);
  fHistos.Add(cyp);
  fHistos.Add(epd);
  gROOT->cd();
}
//=============================================================================
void AliUnicorAnalSingle::Process(AliUnicorEvent *ev) 
{
  // fill single particle histograms

  AliUnicorHN *zep = (AliUnicorHN*) fHistos.At(0);
  AliUnicorHN *cyp = (AliUnicorHN*) fHistos.At(1);
  AliUnicorHN *epd = (AliUnicorHN*) fHistos.At(2);
  for (int i=0; i<ev->NParticles(); i++) {
    if (!ev->ParticleGood(i,fPid)) continue;
    zep->Fill(ev->Zver(),ev->ParticleEta(i),ev->ParticlePhi(i),1.0);
    double y = fMass>0? ev->ParticleY(i,fMass) : ev->ParticleEta(i);
    cyp->Fill(ev->Centrality(),y,ev->ParticlePt(i),1.0);
    epd->Fill(ev->ParticleEta(i),ev->ParticleP(i),sqrt(ev->ParticleDedx(i)),1.0);
  }
}
//=============================================================================
