// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// single particle analyzer
//=============================================================================

#include <cmath>
#include <TROOT.h>
#include <TMath.h>
#include <TAxis.h>
#include <TParticlePDG.h>
#include "AliDHN.h"
#include "AliDEvent.h"
#include "AliDAnalSingle.h"

ClassImp(AliDAnalSingle)

//=============================================================================
AliDAnalSingle::AliDAnalSingle(Char_t *nam, Double_t emi, Double_t ema, Int_t pid) : 
  AliDAnal(nam), fPid(pid), fMass(0.0) 
{
  // constructor
  // emi and ema define the rapidity range for histograms

  fPid = pid;
  TParticlePDG *part = AliDAnal::fgPDG.GetParticle(fPid);
  fMass = part? part->Mass() : 0;

  double pi = TMath::Pi();
  TAxis *ax[10];
  ax[0] = new TAxis(30,-1,1);    ax[0]->SetTitle("vertex z");
  ax[1] = new TAxis(80,emi,ema); ax[1]->SetTitle("eta");
  ax[2] = new TAxis(90,-pi,pi);  ax[2]->SetTitle("phi");
  AliDHN *zep = new AliDHN("zep",3,ax);
  for (int i=0; i<3; i++) delete ax[i];

  ax[0] = new TAxis(20,0,1);     ax[0]->SetTitle("centrality");
  ax[1] = new TAxis(80,emi,ema); ax[1]->SetTitle("y");
  ax[2] = new TAxis(80,0,2);     ax[2]->SetTitle("pt (GeV)");
  AliDHN *cyp = new AliDHN("cyp",3,ax);
  for (int i=0; i<3; i++) delete ax[i];

  ax[0] = new TAxis(10,emi,ema); ax[0]->SetTitle("eta");
  ax[1] = new TAxis(150,0,3);    ax[1]->SetTitle("p (GeV)");
  ax[2] = new TAxis(150,0.5,3.5);ax[2]->SetTitle("sqrt(dedx (mips))");
  AliDHN *epd = new AliDHN("epd",3,ax);
  for (int i=0; i<3; i++) delete ax[i];

  fHistos.Add(zep);
  fHistos.Add(cyp);
  fHistos.Add(epd);
  gROOT->cd();
  printf("%s object named %s created\n",ClassName(),GetName());
}
//=============================================================================
void AliDAnalSingle::Process(AliDEvent *ev) 
{
  // fill single particle histograms

  AliDHN *zep = (AliDHN*) fHistos.At(0);
  AliDHN *cyp = (AliDHN*) fHistos.At(1);
  AliDHN *epd = (AliDHN*) fHistos.At(2);
  for (int i=0; i<ev->NParticles(); i++) {
    if (!ev->ParticleGood(i,fPid)) continue;
    zep->Fill(ev->Zver(),ev->ParticleEta(i),ev->ParticlePhi(i),1.0);
    double y = fMass>0? ev->ParticleY(i,fMass) : ev->ParticleEta(i);
    cyp->Fill(ev->Centrality(),y,ev->ParticlePt(i),1.0);
    epd->Fill(ev->ParticleEta(i),ev->ParticleP(i),sqrt(ev->ParticleDedx(i)),1.0);
  }
}
//=============================================================================
