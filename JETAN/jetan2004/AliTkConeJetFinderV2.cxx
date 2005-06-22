// $Id$

#include <Riostream.h>
#include <vector>
#include <list>
#include <map>

#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "AliTkConeJet.h"
#include "AliTkTowerV2.h"
#include "AliTkConeJetEvent.h"
#include "AliTkConeJetFinderV2.h"

#ifdef ALICEINTERFACE
#include <AliJetParticle.h>
#include <AliJetEvent.h>
#include <AliJetEventParticles.h>
#endif

//numerical precision
#define __epsilon__ 1e-6

//==== Implementation of tower =====
tower::tower() : fEta_min(-999),fEta_max(-999),fEta_center(-999),
                 fPhi_min(-999),fPhi_max(-999),fPhi_center(-999),
                 fEt(0), fParticles(0)
{
}

tower::tower(const tower& t) 
{
  fEta_min = t.getEtaMin();
  fEta_max = t.getEtaMax();
  fEta_center = t.getEta();
  fPhi_min = t.getPhiMin();
  fPhi_max = t.getPhiMax();
  fPhi_center = t.getPhi();
  fEt = t.getEt();
}

tower::tower(Float_t phimin, Float_t phimax, Float_t etamin, Float_t etamax) 
{
  fEta_min = etamin;
  fEta_max = etamax;
  fEta_center = (etamax+etamin)/2.;
  fPhi_min = phimin;
  fPhi_max = phimax;
  fPhi_center = (phimax+phimin)/2.;
  fEt = 0;
}

tower::~tower()
{
  clearParticles();
}

tower& tower::operator+=(const Float_t E) 
{
  fEt += E;
  return *this;
}

tower& tower::operator+=(const TParticle *part) 
{
  fEt += part->Pt();   // Et for a massless particle...
  addParticle(part);
  return *this;
}

ostream& operator<<(ostream& s,const tower& t) 
{
  return s << "tower info: etamin=" << t.getEtaMin()
	   << " etamax=" << t.getEtaMax() 
	   << " phimin=" << t.getPhiMin()
	   << " phimax=" << t.getPhiMax()
	   << " Et=" << t.getEt();
}

list<const TParticle*> *tower::getParticles() 
{
  list<const TParticle*> *newList = new list<const TParticle*>;
  list<const TParticle*> &myList = *newList;
  copy(fParticles.begin(),fParticles.end(),back_inserter(myList));
  return newList;
}

//==== Implementation of protojet =====
protojet:: protojet() : fCentroid(-999,-999),
			fEtWeightedCentroid(-999,-999),
			fEt(-999),fUpdate(kFALSE),fTowers(0)
{
}
  
protojet::protojet(const protojet& p) : 
                        fCentroid(-999,-999),
			fEtWeightedCentroid(-999,-999),
			fEt(-999),fUpdate(kFALSE),fTowers(0)
{
  setCentroidPosition(p.getCentroidPosition());
}
  
protojet::protojet(const protojet *p) : 
                        fCentroid(-999,-999),
			fEtWeightedCentroid(-999,-999),
			fEt(-999),fUpdate(kFALSE),fTowers(0)
{
  setCentroidPosition(p->getCentroidPosition());
}

bool protojet::operator==(const protojet &p1) 
{
  AliTkEtaPhiVector otherCenter = p1.getCentroidPosition();
  if (fCentroid.diffSq(otherCenter) > __epsilon__) 
    return false;
  
  if ((getEt() - p1.getEt())*(getEt() - p1.getEt()) > __epsilon__)
   return false;

  return true;
}

void protojet::update()
{
  Float_t Et = 0;
  Float_t phi = 0;
  Float_t eta = 0;

  for (list<tower*>::const_iterator iter = fTowers.begin();
       iter != fTowers.end(); ++iter) {
    Et += (*iter)->getEt();
    phi += (*iter)->getEt() * (*iter)->getPhi();
    eta += (*iter)->getEt() * (*iter)->getEta();
  }
  if (Et != 0) {
    phi /= Et;
    eta /= Et;
  }
  fEtWeightedCentroid.setVector(eta,phi);
  fEt=Et;
  fUpdate=kFALSE;
}

list<tower*> protojet::getTowerList() const 
{
  list<tower*> newList;
  copy(fTowers.begin(),fTowers.end(),back_inserter(newList));
  return newList;
}

list<tower*> protojet::getSharedTowerList(const protojet *other) const 
{
  list<tower*> newList;
  for (list<tower*>::const_iterator i = fTowers.begin();
       i != fTowers.end(); ++i) {
    if (other->hasTower((*i))) {
      newList.push_back((*i));
    }
  }
  return newList;
}

bool protojet::hasTower(tower *pTower) const 
{
  bool bHasTower = false;
  for (list<tower*>::const_iterator iter = fTowers.begin();
       iter != fTowers.end(); ++iter) {
    if (pTower == (*iter)) {
      bHasTower = true;
      break;
    }
  }
  return bHasTower;
}

bool protojet::shareTowers(const protojet *other) const 
{
  bool bShareTowers = false;
  for (list<tower*>::const_iterator iter = fTowers.begin();
       iter != fTowers.end(); ++iter) {
    if (other->hasTower(*iter)) {
      bShareTowers = true;
      break;
    }
  }
  return bShareTowers;
}

Float_t protojet::diffToCenter(tower *pTower) const 
{
  AliTkEtaPhiVector v(pTower->getEta(),pTower->getPhi());
  AliTkEtaPhiVector center = getCentroidPosition();
  
  return v.diff(center);
}

ostream& operator<<(ostream& s,const protojet& p) 
{
  return s << "Protojet info: eta=" << p.Eta()
	   << " phi=" << p.Phi() << " Et=" << p.getEt();
}


//================== implementation of AliTkConeJetFinderV2 ==============
ClassImp(AliTkConeJetFinderV2)

AliTkConeJetFinderV2::AliTkConeJetFinderV2() 
                  : TObject(),
		    fOutput(0),fNTowers(0),fEtaBins(0),fEtaMin(0),
		    fEtaMax(0),fEtaWidth(0),fPhiBins(0),fPhiMin(0),fPhiMax(0),
		    fPhiWidth(0),fEtCut(0),fEtMinJet(0),fPtCut(0),fRadius(0),fTowers(0),
		    fSeedPointsNew(),fProtojets(),fJets(0),
		    fEvoutfile(0),fEvout_name(0),fEvoutevent(0),fEvouttree(0)
#ifdef ALICEINTERFACE
		    ,fAliParticles(0)
#endif

#ifdef DOHISTOS
		    ,fHistos(),fEventHistos(),fHistFile(0)
#endif
{  
  defaultSettings();
}


AliTkConeJetFinderV2::~AliTkConeJetFinderV2()
{
  if(fEvoutevent) delete fEvoutevent;
  if(fEvoutfile) delete fEvoutfile;
  if(fEvout_name) delete[] fEvout_name;
#ifdef ALICEINTERFACE
  if(fAliParticles) delete fAliParticles;
#endif

#ifdef DOHISTOS
  if(fHistFile) delete fHistFile;
#endif
}

void  AliTkConeJetFinderV2::defaultSettings() 
{
  fOutput = kFALSE;
  fRadius = 0.7;
  // some more or less usefull default settings
  // only possibilty to control the jet finder without setters/getters
  fPhiBins = (Int_t) (2 * TMath::Pi() / 0.1);
  fPhiMin = 0;
  fPhiMax = 2 * TMath::Pi();
  fEtaBins = 20;
  fEtaMin = -1;
  fEtaMax = 1;
  fEtCut = 0;
  fPtCut = 0;
  fEtMinJet = 0;
  fNTowers = fPhiBins*fEtaBins;
}

void AliTkConeJetFinderV2::setSettings(Int_t phibins,Int_t etabins)
{
  fPhiBins = phibins;
  fPhiMin = 0;
  fPhiMax = 2 * TMath::Pi();
  fEtaBins = etabins;
  fEtaMin = -1;
  fEtaMax = 1;
  fNTowers = fPhiBins*fEtaBins;
}

void AliTkConeJetFinderV2::init() 
{
  createTowers();

  fEvoutevent = new AliTkConeJetEvent();
  fEvoutfile = 0;
  fEvouttree = 0;

#ifdef DOHISTOS
  createHistos();
  createEventHistos();
#endif

  if (getEvOutFilename()) {
    if(fOutput)cout << "Writing finder output to " << getEvOutFilename() << endl;
    fEvoutfile = new TFile(getEvOutFilename(),"RECREATE");
    fEvouttree = new TTree("jets","TKConeJetFinderV2 jets");
    fEvouttree->Branch("ConeFinder","AliTkConeJetEvent",&(this->fEvoutevent),32000,0);
  }
}

#ifdef ALICEINTERFACE
void AliTkConeJetFinderV2::initEvent(const AliJetEventParticles *p,TString desc) 
{
  if(fEvoutevent){
    fEvoutevent->Clear();
    fEvoutevent->setDesc(desc);
    //fEvoutevent->setJetParticles(p);
    fEvoutevent->setJetParticles(0);
  }
  const TClonesArray *parts=p->GetParticles();
  initEvent_(parts,2);
}
#endif

void AliTkConeJetFinderV2::initEvent(const TClonesArray *particles,Int_t type,TString desc) 
{
  if(fEvoutevent){
    fEvoutevent->Clear();
    fEvoutevent->setDesc(desc);
  }
  initEvent_(particles,type);
}
	
void AliTkConeJetFinderV2::initEvent(const TClonesArray *particles,Int_t type) 
{
#ifdef DOHISTOS
  // clear event histograms for new event
  clearEventHistos();
#endif

  if(fEvoutevent){
    fEvoutevent->Clear();
  }

  initEvent_(particles,type);
}


void AliTkConeJetFinderV2::initEvent_(const TClonesArray *particles,Int_t type) 
{
  fEvoutevent->setRadius(fRadius);
  fEvoutevent->setEtCut(fEtCut);
  fEvoutevent->setPtCut(fPtCut);

  // reset all towers
  for (vector<tower>::iterator i = fTowers->begin(); i != fTowers->end();++i) {
     (*i).clear();
  }
  // reset seed points
  fSeedPointsNew.erase(fSeedPointsNew.begin(),fSeedPointsNew.end());
  // reset protojet list
  for (list<protojet*>::iterator i = fProtojets.begin();i != fProtojets.end(); ++i) {
    protojet *p = (*i);
    if (p) {
      delete p;
    }
  }
  fProtojets.erase(fProtojets.begin(),fProtojets.end());

  // reset jet list
  for (list<protojet*>::iterator i = fJets.begin();
       i != fJets.end(); ++i) {
    protojet *p = (*i);
    if (p) {
      delete p;
    }
  }
  fJets.erase(fJets.begin(),fJets.end());

  // fill Et towers from particles
  switch (type) {
  case 1:
    fillTowersFromTParticles(particles);
    break;
#ifdef ALICEINTERFACE
  case 2:
    fillTowersFromAliParticles(particles);
    break;
#endif
  default:
    cerr << "AliTkConeJetFinderV2: don't know how to fill Et hist from TClonesArray with type " << type << endl;
  }
  createSeedPoints();

#ifdef DOHISTOS
  fillTowerHist();
#endif
}

Bool_t AliTkConeJetFinderV2::isJetEnergy(Float_t min, Float_t max) 
{ // Checks if there was an jet with Et between min and max
  Float_t men=maxJetEnergy();

  if((men>min)&&(men<max)) return kTRUE;
  return kFALSE;
}

Float_t AliTkConeJetFinderV2::maxJetEnergy() 
{
  Float_t men=-1;
  for (list<protojet*>::const_iterator iter = fJets.begin(); 
       iter != fJets.end(); ++iter) {   

    Float_t et=(*iter)->getEt();
    if(men<et) men=et;
  }
  return men;
}

void AliTkConeJetFinderV2::run() 
{
  findProtojets();
  findJets();
}

void AliTkConeJetFinderV2::finishEvent() 
{
  if (fEvoutevent){ 
    AliTkConeJet *jet = new AliTkConeJet();
    AliTkTowerV2 *myTower = new AliTkTowerV2();
    for (list<protojet*>::const_iterator iter = fJets.begin();
	 iter != fJets.end(); ++iter) {

      Float_t jetet=(*iter)->getEt();
      if(jetet<fEtMinJet) continue;
      jet->Clear();
      jet->setEta((*iter)->Eta());
      jet->setPhi((*iter)->Phi());
      jet->setEt(jetet);

      list<tower*> jettowers = (*iter)->getTowerList();
      for (list<tower*>::const_iterator twiter = jettowers.begin();
	   twiter != jettowers.end(); ++twiter) {
	myTower->Clear();
	myTower->setEta((*twiter)->getEta());
	myTower->setPhi((*twiter)->getPhi());
	myTower->setEt((*twiter)->getEt());
	list<const TParticle*> &twparts = *((*twiter)->getParticles());
	for (list<const TParticle*>::const_iterator partiter = twparts.begin();
	    partiter != twparts.end(); ++partiter) {
	  myTower->addParticle(*partiter);
	}
	delete &twparts;
	jet->addTower(myTower);
      }
      jet->calculateValues();
      fEvoutevent->addJet(jet);
    }
    fEvoutevent->sortJets();

    if (fEvouttree)
      fEvouttree->Fill();
    delete jet;
    delete myTower;
  }
#ifdef DOHISTOS  
  writeEventHistos();
#endif
}

void AliTkConeJetFinderV2::finish() 
{
  if ((fEvoutfile) && (fEvoutfile->IsOpen())) {
    fEvoutfile->cd();
    fEvouttree->Write();
    fEvoutfile->Close();
  }

#ifdef DOHISTOS
  writeHistos();
  fHistFile->Close();
#endif
}

void AliTkConeJetFinderV2::createTowers()
{
  if(fOutput) cout << "Creating " << fNTowers << " tower" << endl;
  if (fPhiMax > fPhiMin) {
    fPhiWidth = (fPhiMax-fPhiMin)/(Float_t)fPhiBins;
  } else {
    fPhiWidth = (fPhiMin-fPhiMax)/(Float_t)fPhiBins;
  }
  if (fEtaMax > fEtaMin) {
    fEtaWidth = (fEtaMax-fEtaMin)/(Float_t)fEtaBins;
  } else {
    fEtaWidth = (fEtaMin-fEtaMax)/(Float_t)fEtaBins;
  }
  if(fOutput) cout << "Delta(eta)=" << fEtaWidth << " Delta(phi)=" << fPhiWidth << endl;

  // let's create the container
  fTowers = new vector<tower>(fNTowers);
  
  // set tower boundaries...
  Float_t eta = fEtaMin;
  Float_t phi = fPhiMin;
  for (Int_t e = 0; e < fEtaBins; e++) {
    for (Int_t p = 0; p < fPhiBins; p++) {
      Int_t tower = e*fPhiBins + p;
      (*fTowers)[tower].setPhiMin(phi);
      (*fTowers)[tower].setPhiMax(phi + fPhiWidth);
      (*fTowers)[tower].setPhi(phi + fPhiWidth/2.);
      (*fTowers)[tower].setEtaMin(eta);
      (*fTowers)[tower].setEtaMax(eta + fEtaWidth);
      (*fTowers)[tower].setEta(eta + fEtaWidth/2.);
				   
      phi += fPhiWidth;
    }
    eta += fEtaWidth;
    phi = fPhiMin;
  }
}

void AliTkConeJetFinderV2::fillTowersFromTParticles(const TClonesArray *particles) 
{
  // input TClonesArray with TParticles, e.g. from PYTHIA
  // fills the Et grid from these particles
  if (!particles) {
    return;
  }

  // loop over all particles...
  TParticle *particle = NULL;
  TIterator *iter = particles->MakeIterator();
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    // check if particle is accepted
    if (1) { //isTParticleAccepted(particle)) {
      Float_t pt=particle->Pt();
      if(pt<fPtCut) continue;
      Int_t tower = findTower(particle->Phi(),particle->Eta());
      if (tower >= 0) {
	(*fTowers)[tower] += particle;
#ifdef DOHISTOS
      // calculate Et for a massless particle...
      // idea for Et check: Et = E *sin(theta);
      Float_t Et = TMath::Sqrt(pt*pt);
	// just as check - fill particle Et in histogram
	((TH2F *)fEventHistos.At(1))->Fill(particle->Eta(),particle->Phi(),Et);
#endif
      }
    }
  }
  delete iter;
}

#ifdef ALICEINTERFACE
void AliTkConeJetFinderV2::fillTowersFromAliParticles(const TClonesArray *particles)
{
  // input TClonesArray with AliParticles
  // fills the Et grid from these particles
  if (!particles) {
    return;
  }

  if(fAliParticles) fAliParticles->Clear();
  else fAliParticles=new TClonesArray("TParticle",0);
  fAliParticles->Expand(particles->GetEntriesFast());

  // loop over all particles...
  AliJetParticle *aliparticle = NULL;
  TIterator *iter = particles->MakeIterator();

  Int_t i=0;
  //Float_t meanet=0,totet=0;
  while ((aliparticle = (AliJetParticle *) iter->Next()) != NULL) {
    Float_t pt=aliparticle->Pt();
    if(pt<fPtCut) continue;
    // particle is accepted through reader in JETAN
    TParticle *particle = new((*fAliParticles)[i]) TParticle(0,0,0,0,0,0,aliparticle->Px(),aliparticle->Py(),
			                           aliparticle->Pz(),aliparticle->Energy(),0,0,0,0);
    //mark particle (-123 for Pythia)
    particle->SetWeight(aliparticle->GetType()); 
    Int_t tower = findTower(particle->Phi(),particle->Eta());
    if (tower >= 0) {
      (*fTowers)[tower] += particle;
      //meanet+=pt;totet+=pt;
#ifdef DOHISTOS
      // calculate Et for a massless particle...
      // idea for Et check: Et = E *sin(theta);
      Float_t Et = TMath::Sqrt(particle->Pt()*particle->Pt());
      // just as check - fill particle Et in histogram
      ((TH2F *)fEventHistos.At(1))->Fill(particle->Eta(),particle->Phi(),Et);
#endif
    }
    i++;
  }
  delete iter;
  //cout << "Mean particle " << meanet/i << " " << totet << endl;
}

#endif

Int_t AliTkConeJetFinderV2::findTower(Float_t phi,Float_t eta) 
{
  if ((phi < fPhiMin) || (phi > fPhiMax) ||
      (eta < fEtaMin) || (eta > fEtaMax)) {
    return -1;
  }
  Int_t phibins = (Int_t) ((phi - fPhiMin) / fPhiWidth);
  Int_t etabins = (Int_t) ((eta - fEtaMin) / fEtaWidth);
  return (etabins * fPhiBins + phibins);
}

void AliTkConeJetFinderV2::createSeedPoints() 
{
  // function should decide if it makes sense 
  // to use fEtCut for seed towers
  // + midpoints or simply tower list
  // uses only tower list so far...

#if 0
  Float_t met=0.;Float_t set=0.;
  Int_t counter=0;
  for(vector<tower>::iterator iter = fTowers->begin();
      iter != fTowers->end(); ++iter) {
    met+=iter->getEt();
    set+=iter->getEt()*iter->getEt();
    counter++;
  }
  if(counter>1){
    met/=counter;
    set=set/counter-met*met;
    if(set>0) set=TMath::Sqrt(set)/(counter-1);
    else set=0;
    if(fOutput) cout << "Tower Mean Et: " << met << "+-" << set << " " << counter << endl;
  }  
  Float_t meanEt=met-set; //store mean tower et

  for(vector<tower>::iterator pos = fTowers->begin();
      pos != fTowers->end(); ++pos) {
    if(pos->getEt()<meanEt){
      pos->clear();
    }
  }
#endif  

  for(vector<tower>::iterator iter = fTowers->begin();
      iter != fTowers->end(); ++iter) {
    if(iter->getEt()>fEtCut){
      AliTkEtaPhiVector seedPoint((*iter).getEta(),(*iter).getPhi());
      fSeedPointsNew.push_back(seedPoint);
    }
  }
  if(fOutput) cout << "created " << fSeedPointsNew.size() << " seed points" << endl;
}

void AliTkConeJetFinderV2::findProtojets() 
{
  const Float_t radiusSq = fRadius*fRadius;
  const Int_t maxIterations = 100; 
  const Bool_t geoCut = kTRUE;
  const Float_t geoCutRadiusSq = 2*(fPhiWidth*fPhiWidth + fEtaWidth*fEtaWidth);

  // loop over all seedpoints
  for(list<AliTkEtaPhiVector>::iterator iter = fSeedPointsNew.begin();
      iter != fSeedPointsNew.end(); ++iter) {
    // create a new protojet at seedpoint position
    protojet *pj = new protojet();
    pj->setCentroidPosition(*iter);
    // loop over all towers and add all within "radius" to the protojet...
    for (vector<tower>::iterator tower = fTowers->begin();
	 tower != fTowers->end(); ++tower) {
      if(tower->getEt()<=0) continue; 
      AliTkEtaPhiVector TwCenter((*tower).getEta(),(*tower).getPhi());
      AliTkEtaPhiVector protojetCenter(pj->getCentroidPosition());
      if (TwCenter.diffSq(protojetCenter) < radiusSq) {
	pj->addTower(&(*tower));
      }
    }
    //update mean values
    pj->update();
    // have added all towers within "radius" to protojet...
    // iterate seedpoint until stable
    // lets get the pT-weigthed center of the protojet...
    AliTkEtaPhiVector centerGeo = pj->getCentroidPosition();
    AliTkEtaPhiVector centerEt = pj->getEtWeightedPosition();
    Int_t iteration = 0;
    Bool_t wasGeoCut = kFALSE;
    while ((centerGeo.diffSq(centerEt) > __epsilon__) &&
	   (iteration < maxIterations)) {
      iteration++;
      pj->eraseTowers();
      pj->setCentroidPosition(centerEt);
      centerGeo = pj->getCentroidPosition();
      
      // if geoCut == kTRUE, break if it leaves original bin...
      if ((geoCut == kTRUE) &&
	  ((*iter).diffSq(centerGeo) > geoCutRadiusSq)) {
	wasGeoCut = kTRUE;
	break;
      }
      // loop over all towers and add all within "radius" to the protojet...
      for (vector<tower>::iterator tower = fTowers->begin();
	   tower != fTowers->end(); ++tower) {
	if(tower->getEt()<=0) continue;
	AliTkEtaPhiVector TwCenter((*tower).getEta(),(*tower).getPhi());
	if (TwCenter.diffSq(centerEt) < radiusSq) {
	  pj->addTower(&(*tower));
	}
      }
      // have added all towers within "radius" to protojet...
      centerEt = pj->getEtWeightedPosition();
    }

    // we have a stable protojet (or run out of iterations...)
#if 0
    if((fOutput) &&(iteration==maxIterations))
      cout << " FindProtojets warning: max iterations of " 
      << maxIterations << " reached !!!" << endl;
#endif
    // let's add it to the protojet list...
    if (!wasGeoCut) {
      addProtojet(pj);
    } else {
      delete pj;
    }
  }
  //if(fOutput) cout << "found " << fProtojets.size() << "proto jets" << endl;
}

void AliTkConeJetFinderV2::addProtojet(protojet *pj) 
{
  if (!pj) {
    return;
  }
  for(list<protojet*>::const_iterator iter = fProtojets.begin();
      iter != fProtojets.end(); ++iter) {
    if ((*pj) == (*(*iter))) {
      delete pj;
      return;
    }
  }
  fProtojets.push_back(pj);
}

void AliTkConeJetFinderV2::dumpProtojets(Float_t etmin) 
{
  for(list<protojet*>::const_iterator iter = fProtojets.begin();
      iter != fProtojets.end(); ++iter) {
    if ((*iter)->getEt() > etmin) {
      cout << (*(*iter)) << endl;
    }
  }
}

void AliTkConeJetFinderV2::findJets() 
{
  // loop over all protojets until list is empty
  while(!fProtojets.empty()) {
    // find the protojet with maximum Et
    // should be easy to have a sorted list...
    list<protojet*>::iterator maxEtProtojet = fProtojets.begin();
    Float_t maxEt = 0;
    for (list<protojet*>::iterator iter = fProtojets.begin();
	 iter != fProtojets.end(); ++iter) {
      if ((*iter)->getEt() > maxEt) {
	maxEt = (*iter)->getEt();
	maxEtProtojet = iter;
      }
    }
    // we've found the protojet with the highest Et - remove it from the list
    protojet *jet1 = *maxEtProtojet;
    fProtojets.erase(maxEtProtojet);
    // loop again over all protojets to find sharing jet with highest Et
    list<protojet*>::iterator maxEtNeighbor = fProtojets.begin();
    maxEt = 0;
    for (list<protojet*>::iterator iter = fProtojets.begin();
	 iter != fProtojets.end(); ++iter) {
      if (((*iter)->getEt() > maxEt) &&
	  (jet1->shareTowers(*iter))) {
	maxEt = (*iter)->getEt();
	maxEtNeighbor = iter;
      }
    }
    if (maxEt>0) {
      // jet's share towers
      // merging splitting step...
      protojet *jet2 = (*maxEtNeighbor);
      fProtojets.erase(maxEtNeighbor);
      splitMergeJets(jet1,jet2);
    } else {
      // protojet 1 doesn't share towers with other protojets, make it a jet...
      addJet(jet1);
    }
  }
  if(fOutput){
    cout << "found " << fJets.size() << " jets" << endl;
    dumpJets();
  }
}

void AliTkConeJetFinderV2::splitMergeJets(protojet *jet1,protojet *jet2) 
{
  const Float_t EtRatioCut = 0.5;

  // let's calcualte the shared energy...
  Float_t fEtShared = 0.0;
  list<tower*> sharedTowers = jet1->getSharedTowerList(jet2);
  for (list<tower*>::const_iterator iter = sharedTowers.begin();
       iter != sharedTowers.end(); ++iter) {
    fEtShared += (*iter)->getEt();
  }
  Float_t fEtJet2 = jet2->getEt();
  // calculate the ratio of the shared energy and decided if split or merge
  Float_t fEtRatio = EtRatioCut + 1.;
  if(fEtJet2) fEtRatio = fEtShared / fEtJet2;
  if (fEtRatio > EtRatioCut) {
    mergeJets(jet1,jet2);
  } else {
    splitJets(jet1,jet2);
  }
}

void AliTkConeJetFinderV2::mergeJets(protojet *jet1,protojet *jet2) {
  // merge protojets...
  list<tower*> jet2Towers = jet2->getTowerList();
  for (list<tower*>::const_iterator iter = jet2Towers.begin();
       iter != jet2Towers.end(); ++iter) {
    // add towers from protojet2 to protojet1...
    // don't add shared towers...
    if (!jet1->hasTower((*iter))) {
      jet1->addTower((*iter));
    }
  }
  addProtojet(jet1);
  delete jet2;
}

void AliTkConeJetFinderV2::splitJets(protojet *jet1,protojet *jet2) 
{
  // split protojets...
  list<tower*> sharedTowers = jet1->getSharedTowerList(jet2);
  for (list<tower*>::const_iterator iter = sharedTowers.begin();
       iter != sharedTowers.end(); ++iter) {
    // search nearest jet
    if (jet1->diffToCenter(*iter) < jet2->diffToCenter(*iter)) {
      // tower closer to jet1
      jet2->eraseTower(*iter);
    } else {
      // tower closer to jet2
      jet1->eraseTower(*iter);
    }
  }
  if (jet1->shareTowers(jet2)) {
    cerr << "!!! SplitJets: Something is wrong !!!" << endl;
    addProtojet(jet1);
    delete jet2;
    return;
  }
  addProtojet(jet1);
  addProtojet(jet2);
}

void AliTkConeJetFinderV2::dumpJets() 
{
  cout << "----- found jets > " << fEtMinJet << " GeV -----" << endl;
  for (list<protojet*>::const_iterator iter = fJets.begin();
       iter != fJets.end(); ++iter) {
    if((*iter)->getEt()>fEtMinJet) cout << (*(*iter)) << endl;
  }
}

#ifdef DOHISTOS
void AliTkConeJetFinderV2::createHistos() 
{
}
#endif

#ifdef DOHISTOS
void AliTkConeJetFinderV2::createEventHistos() 
{
  fHistFile = new TFile("$JF_DATADIR/ConeFinderV2.root","RECREATE");
  TH2F *h = new TH2F("Etbins","Etbins",
                     fEtaBins,fEtaMin,fEtaMax,
		     fPhiBins,fPhiMin,fPhiMax);
  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("#phi");
  h->GetZaxis()->SetTitle("E_{t} (GeV)");
  fEventHistos.Add(h);

  h = new TH2F("Etbin_check","Etbin_check",
               fEtaBins,fEtaMin,fEtaMax,
	       fPhiBins,fPhiMin,fPhiMax);
  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("#phi");
  h->GetZaxis()->SetTitle("E_{t} (GeV)");
  fEventHistos.Add(h);

  h = new TH2F("seedpoints","seedpoints",
	       fEtaBins,fEtaMin,fEtaMax,
	       fPhiBins,fPhiMin,fPhiMax);
  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("#phi");
  fEventHistos.Add(h);
}
#endif
#ifdef DOHISTOS
void AliTkConeJetFinderV2::clearEventHistos() 
{
  ((TH2F *)fEventHistos.At(0))->Reset();
  ((TH2F *)fEventHistos.At(1))->Reset();
  ((TH2F *)fEventHistos.At(2))->Reset();
}
#endif

#ifdef DOHISTOS
void AliTkConeJetFinderV2::writeEventHistos() 
{
  fHistFile->cd();
  ((TH2F *)fEventHistos.At(0))->Write();
  ((TH2F *)fEventHistos.At(1))->Write();
  ((TH2F *)fEventHistos.At(2))->Write();
}
#endif

#ifdef DOHISTOS
void AliTkConeJetFinderV2::writeHistos() 
{
}
#endif

#ifdef DOHISTOS
void AliTkConeJetFinderV2::fillTowerHist() 
{
  for (vector<tower>::iterator i=(*fTowers).begin(); i!=(*fTowers).end(); ++i) {
    tower t = *i;
    Float_t eta = (t.getEtaMax() + t.getEtaMin()) / 2.;
    Float_t phi = (t.getPhiMax() + t.getPhiMin()) / 2.;
    ((TH2F *)fEventHistos.At(0))->Fill(eta,phi,t.getEt());
  }
}
#endif

void AliTkConeJetFinderV2::setEvOutFilename(const Char_t *filename) 
{
  if (!fEvout_name) {
    fEvout_name = new Char_t[4096];
  }

  strcpy(fEvout_name,filename);
}

