// $Id$
//--------------------------------------------------------------------------
// AliTkChargedJetFinder.cxx
// implementation of a simple jet finder for charged particles
// based on CDF PRD 65, 092002 (2002)
// T. Kollegger <kollegge@ikf.physik.uni-frankfurt.de> 10/12/2002
//--------------------------------------------------------------------------

#include <Riostream.h>
//-----------------------------------------------------------------------
// STL includes
#include <list>
#include <map>
//-----------------------------------------------------------------------
// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TVector2.h>
#include <TClonesArray.h>
#include <TParticle.h>
//--------------------------------------------------------------------------
#include "AliTkEtaPhiVector.h"
#ifdef DOCHAERGED
#include "TkChargedJet.h"
#else
#include "AliTkConeJet.h"
#include "AliTkConeJetEvent.h"
#include "AliTkTowerV2.h"
#endif
//--------------------------------------------------------------------------
#ifdef ALICEINTERFACE
#include <AliJetParticle.h>
#include <AliJetEvent.h>
#include <AliJetEventParticles.h>
#endif
//-------------------------------------------------------------------------
#include "AliTkChargedJetFinder.h"

//------------------------------------------------------------------------
// implementation of the jet helper class
//------------------------------------------------------------------------

jet::jet(){
  fPt=-999;
  fEta=-999;
  fPhi=-999;
  fNParticles=0;
}

AliTkEtaPhiVector jet::getCentroid() const
{
  return AliTkEtaPhiVector(fEta,fPhi);
}

void jet::addParticle(TParticle *particle) 
{
  // add a particle to the jet
  fParticles.push_back(particle);
  if(fNParticles==0) {
    fEta=particle->Eta();
    fPhi=particle->Phi();
    fPt=0;
  }
  fPt+=particle->Pt();
  fNParticles++;
}

Double_t jet::getPt() const
{
#if 1
  return fPt;
#else
  // returns the pt of the jet
  Double_t pt = 0.0;

  // loop over all particles which are assigned to the jet...
  // and sum the pt of them
  for(list<TParticle *>::const_iterator i = particles.begin();
      i != particles.end(); ++i) {
    TParticle *particle = *i;
    pt += particle->Pt();
  }

  return pt;
#endif
}

Int_t jet::getNParticles() const
{
  return fNParticles;
}

TParticle *jet::getParticle(Int_t i)  const
{
  if ((i < 0) || ((UInt_t)i > fParticles.size()-1)) {
    cerr << "jet:: out of range" << i << endl;
    return NULL;
  }
  list<TParticle *>::const_iterator iter = fParticles.begin();
  for(Int_t pos = 0; pos < i+1; pos++) {
    ++iter;
  }
  return (*iter);
}
    
TClonesArray *jet::getParticles() const
{
  TClonesArray *parts = new TClonesArray("TParticle",this->getNParticles());
  Int_t i = 0;
  for(list<TParticle *>::const_iterator iter = fParticles.begin();
      iter != fParticles.end(); ++iter, i++) {
    TParticle *particle = *iter;
    if (particle == NULL) {
      cerr << "jet: What's wrong? NULL pointer in particles" << endl;
      continue;
    }
    new ((*parts)[i]) TParticle(*particle);
  }
  return parts;
}

Double_t jet::getDiff(TParticle *particle) const
{
  // calculate the difference between jet center and particle in eta-phi space
  return sqrt(getDiffSq(particle));
}

Double_t jet::getDiffSq(TParticle *particle) const
{
  // calculate the square of the difference between jet center
  // and particle in eta-phi space
  AliTkEtaPhiVector jetCenter = getCentroid();
  AliTkEtaPhiVector particlePos(particle->Eta(),particle->Phi());
  return jetCenter.diffSq(particlePos);
}

ostream& operator<< (ostream& s,jet& j) 
{
  return s << "jet: position "  << j.getCentroid() 
	   << " #particles="    << j.getNParticles() 
	   << " pT="            << j.getPt();
}


//------------------------------------------------------------------------
// implementation of the jet finder class
//------------------------------------------------------------------------

ClassImp(AliTkChargedJetFinder)

AliTkChargedJetFinder::AliTkChargedJetFinder() : TObject() 
{  
  fOutput=0;
  fR=0;
  fRSq=0;
  fEtaMin=0;
  fEtaMax=0;
  fPtCut=0;
  fPtSeed=0;
  fMinJetPt=0;
  
#ifdef DOCHARGED
  fMyTJet=0;
#else
  fEvoutevent=0;
#endif

#ifdef DOHISTOS
  // histograms...
  fOutput_name=0;
  fHistos=0;
  fHistoutfile=0;
#endif

  fEvout_name=0;
  fEvoutfile=0;
  fEvouttree=0;

#ifdef ALICEINTERFACE
  fAliParticles=0;
#endif

  defaultSettings();
}

AliTkChargedJetFinder::~AliTkChargedJetFinder()
{
  if(fEvoutevent) delete fEvoutevent;
  if(fEvoutfile) delete fEvoutfile;
  if(fEvout_name) delete[] fEvout_name;
#ifdef ALICEINTERFACE
  if(fAliParticles) delete fAliParticles;
#endif
#ifdef DOHISTOS
  if(fHistos) delete fHistos;
  if(fHistoutFile) delete fHistoutFile;
  if(fOutput_name) delete[] fOutput_name;
#endif
}

void AliTkChargedJetFinder::defaultSettings() 
{
  // default settings for a first study...
  // jet finder radius
  setFinderR(0.7);
  // lower bound of eta range
  setEtaMin(-1);
  // high bound of eta range
  setEtaMax(1);
  // seed value
  setPtSeed(5);
  // minimum jet value
  setMinJetPt(5);

  fEvout_name=0;
  setEvOutFilename("$JF_DATADIR/charged_jets.evout.root");
#ifdef DOHISTOS
  fOutput_name=0;
  setHistFilename("$JF_DATADIR/charged_jets.hist.root");
#endif
}

void AliTkChargedJetFinder::init() 
{
#ifdef DOHISTOS
  // initalization of the jet finder
  fHistos = new TObjArray(10);
  TH1F *hist = new TH1F("ncharged","charged tracks",31,-0.5,30.5);
  hist->GetXaxis()->SetTitle("n_{charged} tracks");
  hist->GetYaxis()->SetTitle("entries");
  fHistos->Add(hist);
  hist = new TH1F("jet_pt","jet p_{T}",100,0,100);
  hist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hist->GetYaxis()->SetTitle("entries");
  fHistos->Add(hist);

  // declare all histograms
  if (fOutput_name) {
    if(fOutput) cout << "Writing summary output to " << getHistFilename() << endl;
    fHistoutfile = new TFile(getHistFilename(),"RECREATE");
  }
#endif

#ifdef DOCHARGED
  if (fEvout_name) {
    if(fOutput) cout << "Writing event output to " << getEvOutFilename() << endl;
    fEvoutfile = new TFile(getEvOutFilename(),"RECREATE");
    fMyTJet = new TkChargedJet();
    fEvouttree = new TTree("jets","");
    fEvouttree->Branch("jet","TkChargedJet",&(this->fMyTJet));
  }
#else
  fEvoutfile = 0;
  fEvouttree = 0;
  fEvoutevent = new AliTkConeJetEvent();

  if (getEvOutFilename()) {
    if(fOutput) cout << "Writing finder output to " << getEvOutFilename() << endl;
    fEvoutfile = new TFile(getEvOutFilename(),"RECREATE");
    fEvouttree = new TTree("jets","TKConeJetFinderV2 jets");
    fEvouttree->Branch("ConeFinder","AliTkConeJetEvent",&(this->fEvoutevent),32000,0);
  }
#endif
}

void AliTkChargedJetFinder::initEvent(TClonesArray *newParticles,Int_t type) 
{
  // initalizes the jet finder for a new event
#ifndef DOCHARGED
  if(fEvoutevent){
    fEvoutevent->Clear();
    fEvoutevent->setRadius(fR);
    fEvoutevent->setEtCut(fPtSeed);
    fEvoutevent->setPtCut(fPtCut);
    fEvoutevent->setDesc(*new TString("AliTkChargedJetFinder"));
  }
#endif

  // delete old particle list and old jets
  fParticles.erase(fParticles.begin(),fParticles.end());
  fJets.erase(fJets.begin(),fJets.end());

  TParticle *particle;
  TIterator *iter = newParticles->MakeIterator();
#ifdef testptr
  list<TParticle *>::const_iterator siter;
#endif

  // create new particle list
  switch (type) {
  case 1:
    // TParticles in TClonesArray
    while ((particle = (TParticle *) iter->Next()) != NULL) {
      if (isTParticleAccepted(particle)) {
	addParticle(particle);
      }
    }
    if(fOutput) cout << "found particles for jets: " << fParticles.size() << endl;
    break;
  case 2:
    cerr << "AliTkChargedJetFinder: don't know how to fill Et hist from TClonesArray with type " << type << endl;
    ;
    break;
  }
#ifdef testptr
    for (siter=particles.begin(); siter != particles.end(); ++siter) {
      cout << "particle eta=" << (*siter)->Eta()
	   << " phi=" << (*siter)->Phi()
	   << " pT=" << (*siter)->Pt()
	   << " pointer=" << (*siter)
	   << endl;
    }
#endif
}

#ifdef ALICEINTERFACE
void AliTkChargedJetFinder::initEvent(const AliJetEventParticles *p,TString desc)
{
  if(fEvoutevent){
    fEvoutevent->Clear();
    fEvoutevent->setRadius(fR);
    fEvoutevent->setEtCut(fPtSeed);
    fEvoutevent->setPtCut(fPtCut);
    fEvoutevent->setDesc(desc);
    //fEvoutevent->setJetParticles(p);
  }

  // delete old particle list and old jets
  fParticles.erase(fParticles.begin(),fParticles.end());
  fJets.erase(fJets.begin(),fJets.end());

  const TClonesArray *particles=p->GetParticles();
  if(fAliParticles) fAliParticles->Clear();
  else fAliParticles=new TClonesArray("TParticle",0);
  fAliParticles->Expand(particles->GetEntriesFast());

  // loop over all particles...
  Int_t i=0;
  AliJetParticle *aliparticle = NULL;
  TIterator *iter = particles->MakeIterator();
  while ((aliparticle = (AliJetParticle *) iter->Next()) != NULL) {
    Float_t pt=aliparticle->Pt();
    if(pt<fPtCut) continue;
    // particle is accepted through reader in JETAN
    TParticle *particle = new((*fAliParticles)[i]) TParticle(0,0,0,0,0,0,aliparticle->Px(),aliparticle->Py(),
			                           aliparticle->Pz(),aliparticle->Energy(),0,0,0,0);
    //mark particle (-123 for Pythia)
    particle->SetWeight(aliparticle->GetType()); 
    //particle is accepted through reader
    addParticle(particle);
    i++;
  }
  delete iter;
}
#endif

void AliTkChargedJetFinder::run() 
{
  // loop over particles as long as there are 
  // some which are not assigned to a jet
  if (fParticles.empty()) {
    return;
  }

  list<TParticle *>::iterator iter;
  while (!fParticles.empty()) {
    // get particle with highest pT...
    iter = fParticles.begin();
    if((*iter)->Pt()<fPtSeed) {
      fParticles.erase(iter);
      continue;
    }
    // create a new jet...
    jet myJet;
    // and add the particle with highest pT...
    myJet.addParticle(*iter);	
    // delete it from the list of unassigned particles...
    fParticles.erase(iter);
    // loop over all particles to see if they are within the new jet...
    for(iter = fParticles.begin(); iter != fParticles.end(); ++iter) {  
      if (myJet.getDiffSq(*iter) < fRSq) {
	// add this particle to the jet...
	myJet.addParticle(*iter);
	// remove this particle from the list of unassigned ones...
	fParticles.erase(iter);
	// and start again from the beginnig of the list
	// -- the pT-weigthed center of the jet might have changed,
	//    so one might add now a particle which was before rejected
	iter = fParticles.begin();
	iter--;
      }
    }
    // we have our final jet - add it to the jet list...
    fJets.push_back(myJet);
  }
}

void AliTkChargedJetFinder::finishEvent() 
{
  // analyse event...
  checkJets();

#ifdef DOCHARGED
  list<jet>::iterator highjet = findHighestJet();
  if (highjet == NULL) {
    cerr << "AliTkChargedJetFinder:: no jet found - no output! " << endl;
    return;
  }
  if(fOutput) {
    cout << "high jet ncharged=" << (*highjet).getNParticles()
	 << " pt="               << (*highjet).getPt() << endl;
    cout << *highjet << endl;
  }

#ifdef DOHISTOS
  // write out event histograms
  if ((*highjet).getPt() > 30) {
    ((TH1F *)fHistos->At(0))->Fill((*highjet).getNParticles());
  }
  ((TH1F *)fHistos->At(1))->Fill((*highjet).getPt());
#endif

  TkChargedJet *jet = new TkChargedJet(*highjet);
  if ((fEvoutfile) && (fEvoutfile->IsOpen())) {
    delete fMyTJet;
    if (jet) {
      fMyTJet = jet;
    } else {
      fMyTJet = new TkChargedJet();
    }
    fEvoutfile->cd();
    fEvouttree->Fill();
  }
#else /******************************************************/
  if (fEvoutevent){ 
    AliTkConeJet *conejet = new AliTkConeJet();
    AliTkTowerV2 *myTower = new AliTkTowerV2();
    list<jet>::iterator iter;
    for (iter = fJets.begin();
	 iter != fJets.end(); ++iter) {

      Float_t jetet=(*iter).getPt();
      if(jetet<fMinJetPt) continue;
      conejet->Clear();
      conejet->setEta((*iter).Eta());
      conejet->setPhi((*iter).Phi());
      conejet->setEt(jetet);

      myTower->Clear();
      myTower->setEta((*iter).Eta());
      myTower->setPhi((*iter).Phi());
      myTower->setEt(jetet);
      myTower->setParticleList((*iter).getParticles());
      conejet->addTower(myTower);
      conejet->calculateValues();
      fEvoutevent->addJet(conejet);
    }
    fEvoutevent->sortJets();
    fEvoutevent->Print("");
    if (fEvouttree)
      fEvouttree->Fill();
    delete conejet;
    delete myTower;
  }
#endif
}

void AliTkChargedJetFinder::finish() 
{
#ifdef DOHISTOS
  // write out histograms
  if ((fHistoutfile) && (fHistoutfile->IsOpen())) {
    fHistoutfile->cd();
    ((TH1F *)fHistos->At(0))->Write();
    ((TH1F *)fHistos->At(1))->Write();
    fHistoutfile->Write();
    fHistoutfile->Close();
  }
#endif
  if ((fEvoutfile) && (fEvoutfile->IsOpen())) {
    fEvoutfile->cd();
    fEvouttree->Write();
#ifdef DOHISTOS
    ((TH1F *)fHistos->At(0))->Write();
#endif
    fEvoutfile->Write();
    fEvoutfile->Close();
  }
}

void AliTkChargedJetFinder::setEvOutFilename(const Char_t *filename) 
{
  if (!fEvout_name) {
    fEvout_name = new Char_t[4096];
  }
  strcpy(fEvout_name,filename);
}

#ifdef DOHISTOS
void AliTkChargedJetFinder::setHistFilename(const Char_t *filename) 
{
  if (!fOutput_name) {
    fOutput_name = new Char_t[4096];
  }
  strcpy(fOutput_name,filename);
}
#endif

bool AliTkChargedJetFinder::isTParticleAccepted(TParticle *particle) 
{
  // particle cuts...

  // check if particle is stable
  if (particle->GetStatusCode() != 1) return false;

  // check for charge
  TParticlePDG *partPDG;
  partPDG = particle->GetPDG();
  if (partPDG->Charge() == 0) {
    return false;
  }

  // check eta range
  if ((particle->Eta() < fEtaMin) || (particle->Eta() > fEtaMax)) {
    return false;
  }

  return true;
}

void AliTkChargedJetFinder::addParticle(TParticle *particle) 
{
  // add particle at right position
  // first point to speed up program, use a sorted tree for example...

  // if first particle, store at beginning
  if (fParticles.empty()) {
    fParticles.push_front(particle);
    return;
  }
  
  // search insert position
  list<TParticle *>::iterator iter;
  for (iter = fParticles.begin(); iter != fParticles.end(); ++iter) {
    if ((*iter)->Pt() < particle->Pt()) {
      fParticles.insert(iter,particle);
      break;
    }
  }
}

list<jet>::iterator AliTkChargedJetFinder::findHighestJet() 
{
  list<jet>::iterator iter;
  list<jet>::iterator high;

  if (fJets.empty()) {
    return fJets.begin();
  }

  high = fJets.begin();
  for (iter = fJets.begin(); iter != fJets.end(); ++iter) {
    //cout << *iter << endl;
    if ((*iter).getPt() > (*high).getPt()) {
      high = iter;
    }
  }
  return high;
}

void AliTkChargedJetFinder::checkJets() 
{
  list<jet>::iterator i;
  list<jet>::iterator j;
  for (i = fJets.begin(); i != fJets.end(); ++i) {
    for (j = i, ++j; j != fJets.end(); ++j) {
      AliTkEtaPhiVector v1 = (*i).getCentroid();
      AliTkEtaPhiVector v2 = (*j).getCentroid();
      Double_t diff = v1.diff(v2);
      if (diff < fR) {
	cerr << "AliTkChargedJetFinder: Something is wrong - jets to close..." << endl;
	cerr << v1 << " pt=" << (*i).getPt() <<  " -- " 
	     << v2 << " pt=" << (*j).getPt() << endl;
      }
    }
  }
}
