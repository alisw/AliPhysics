// $Id$

// includes 
#include <Riostream.h>
#include <list>
#include <queue>
#include <iterator>

#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TStopwatch.h>

#include "AliTkKtJetFinder.h"

void AliTkKtJetFinder::init() {
  DebugOutput("TkJetFinder::init() called");

  cout << "TkJetFinder::init() - finder initialized" << endl;
  cout << "precluster parameter:" << endl;
  cout << " - phi binning: " << phiBins << " bins from "
       << phiMin << " to " << phiMax << ", bin width "
       << (phiMax-phiMin)/(Float_t)phiBins << endl;
  cout << " - theta binning: " << thetaBins << " bins from "
       << thetaMin << " to " << thetaMax << ", bin width "
       << (thetaMax-thetaMin)/(Float_t)thetaBins << endl;
  cout << "jetfinder parameter:" << endl;
  cout << " - D:" << finder_D << endl
       << " - D_Cut:" << finder_DCut << endl;

  DebugOutput("TkJetFinder::init() finished");
}

void AliTkKtJetFinder::makeTParticles(TClonesArray *particles) {
  DebugOutput("TkJetFinder::makeTParticle() called");
  initEvent();
  preclusterTParticles(particles);
  DebugOutput("TkJetFinder::makeTParticle() preclustering done");
  findJets();
  DebugOutput("TkJetFinder::makeTParticle() jetfinding done");
  DebugOutput("TkJetFinder::makeTParticle() finished");
}

void AliTkKtJetFinder::clear() {
  DebugOutput("TkJetFinder::clear() called");

  DebugOutput("TkJetFinder::clear() finished");
}

void AliTkKtJetFinder::finish() {
  DebugOutput("TkJetFinder::finish() called");
  
  DebugOutput("TkJetFinder::finish() finished");
}

void AliTkKtJetFinder::initEvent() {
  DebugOutput("TkJetFinder::initEvent() called");

  // clean precluster list
  if (firstPreCluster) {
    TPreCluster *act = firstPreCluster;
    TPreCluster *next;
    while (act != NULL) {
      next = act->next;
      delete act;
      act = next;
    }
  }

  // delete precluster array
  delete[] preClusterArray;

  // reset list
  firstPreCluster = NULL;
  lastPreCluster = NULL;
  preClusterUID = 0;
  preClusterArray = NULL;

  // delete dlist;
  for (list<Tdlist *>::iterator i = myDList.begin(); i != myDList.end(); ++i) {
    Tdlist *dlist = *i;
    delete dlist;
    dlist = NULL;
  }
  myDList.clear();
  
  
  DebugOutput("TkJetFinder::initEvent() finished");
}

void AliTkKtJetFinder::preclusterTParticles(TClonesArray *particles) {
  Char_t buffer[1024];
  Int_t acceptedParticles = 0;
  Int_t acceptedPreCluster = 0;

  DetailedOutput("TkJetFinder::preclusterTParticle() called");

  TStopwatch timer;
  timer.Reset();

  sprintf(buffer,"TkJetFinder::preclusterTParticle() - found %d particles",
	  particles->GetEntries());
  DetailedOutput(buffer);

  cout << "WARNING - TkJetFinder::preclusterTParticle()" << endl
       << "          no preclustering implemented - each particle accepted!"
       << endl;

  timer.Start();
  // start the particle loop
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  TPreCluster *precluster;
  while((particle = (TParticle *) iter->Next()) != NULL) {
    if (isTParticleAccepted(particle)) {
      acceptedParticles++;
      precluster = new TPreCluster;
      // to be kind of realistic use maseless particles 
      // - TPC doesn't measure energy and there is no hadron calorimeter...
      precluster->id = preClusterUID;
      precluster->E = TMath::Sqrt(particle->P()*particle->P());
      precluster->px = particle->Px();
      precluster->py = particle->Py();
      precluster->pz = particle->Pz();
      precluster->next = NULL;
      if (true) {
	// add precluster to list of preclusters
	addPreCluster(precluster);
	acceptedPreCluster++;
	preClusterUID++;
      } else {
	delete precluster;
      }
    }
  }

  // build an array of pointers for direct access on the precluster...
  preClusterArray = new TPreCluster*[preClusterUID];
  precluster = firstPreCluster;
  while (precluster) {
    preClusterArray[precluster->id] = precluster;
    precluster = precluster->next;
  }

  timer.Stop();

  // if you want to be really sure what the precluster algorithm does...
  //dumpPreClusters();
  //dumpPreClusterArray();

  // status information at the end...
  sprintf(buffer,"TkJetFinder::preclusterTParticle() - used %d particles",
	  acceptedParticles);
  DebugOutput(buffer);
  sprintf(buffer,
	  "TkJetFinder::preclusterTParticle() - created %d preclusters",
	  acceptedPreCluster);
  DebugOutput(buffer);
  
  sprintf(buffer,
	  "TkJetFinder::preclusterTParticle() - timing CPU: %f, real %f",
	  timer.CpuTime(),timer.RealTime());
  TimingOutput(buffer);
  
  cout << "some stupid tests..." << endl;
  Tdlist *p1 = new Tdlist;
  Tdlist *p2 = new Tdlist;
  Tdlist *p3 = new Tdlist;
  p1->d = 3;
  p2->d = 2;
  p3->d = 1;
  if (p1==p2) {
    cout << "p1==p2 (" << p1->d << "==" << p2->d << ")" << endl;
  } else {
    cout << "p1!=p2 (" << p1->d << "!=" << p2->d << ")" << endl;
  }
  if (p1 < p2) {
    cout << "p1 < p2 (" << p1->d << ">" << p2->d << ")" << endl;
  } else {
    cout << "p1 > p2 (" << p1->d << "<" << p2->d << ")" << endl;
  }
  priority_queue<Tdlist *> myHeap;
  myHeap.push(p2);
  myHeap.push(p1);
  myHeap.push(p3);

  cout << myHeap.size() << endl;
  while (!myHeap.empty()) {
    Tdlist *myEntry = myHeap.top();
    myHeap.pop();
    cout << "d = " << myEntry->d << endl;
    delete myEntry;
  }

  DetailedOutput("TkJetFinder::preclusterTParticle() finished");
}

void AliTkKtJetFinder::findJets() {
  Char_t buffer[1024];

  TStopwatch timer;
  timer.Reset();

  DetailedOutput("TkJetFinder::findJets() called");

  timer.Start();
  // let's start...
  buildNewDList();
  // find the precluster/pair with lowest dmin...
  //list<Tdlist *>::iterator i;
  list<Tdlist *>::iterator min;

  myDList.sort();
  while(!myDList.empty()) {
    min = myDList.begin();

    // we have the entry with the lowest min...
    //cout << (*min)->d << endl;
    delete *min;
    myDList.erase(min);
  }
    

  timer.Stop();

  sprintf(buffer,
	  "TkJetFinder::findJets() - timing CPU: %f, real %f",
	  timer.CpuTime(),timer.RealTime());
  TimingOutput(buffer);

  DetailedOutput("TkJetFinder::findJets() finish");
}

void AliTkKtJetFinder::setPhiBins(Int_t nPhiBins, 
			       Float_t fPhiMin, Float_t fPhiMax) {
  setNPhiBins(nPhiBins);
  setPhiMin(fPhiMin);
  setPhiMax(fPhiMax);
}

void AliTkKtJetFinder::setThetaBins(Int_t nThetaBins,
				 Float_t fThetaMin, Float_t fThetaMax) {
  setNThetaBins(nThetaBins);
  setThetaMin(fThetaMin);
  setThetaMax(fThetaMax);
}


void AliTkKtJetFinder::setDefaultOptions() {
  // precluster options
  setPhiBins((Int_t)(TMath::Pi()*2./0.1), 0. ,2.*TMath::Pi());
  setThetaBins(10,(TMath::Pi()/2. - 0.5),(TMath::Pi()/2. + 0.5));

  // finder options
  setD(1.);
  setDCut(-1.);

  // debug options
  setDebugLevel(0);
  setDebugFilename("KtJetFinderDebug.root");

  cout << "TkJetFinder::setDefaultOptions() - DEFAULT options set" << endl;
}

Bool_t AliTkKtJetFinder::isTParticleAccepted(TParticle *particle) {
  // check if particle is accepted
  // makes sense to write this into a own class, but now I'm lazy

  // check if particle is stable -> a detectable particle
  if (particle->GetStatusCode() != 1) {
    return kFALSE;
  }
  // and now just for fun
  TParticlePDG *partPDG;
  partPDG = particle->GetPDG();
//   if (partPDG->Charge() == 0) {
//     cout << partPDG->GetName() << endl;
//     return kFALSE;
//   }

  // default case: accept
  return kTRUE;
}

void AliTkKtJetFinder::addPreCluster(TPreCluster *precluster) {
  if (firstPreCluster == NULL) {
    firstPreCluster = precluster;
    lastPreCluster = precluster;
    return;
  }

  lastPreCluster->next = precluster;
  precluster->next = NULL;
  lastPreCluster = precluster;
}

void AliTkKtJetFinder::deletePreCluster(Int_t UID) {
  TPreCluster *p1 = NULL;
  TPreCluster *p2 = firstPreCluster;
  while (p2!=NULL) {
    if (p2->id == UID) {
      p1->next = p2->next;
      delete p2;
      break;
    }
    p1 = p2;
    p2 = p2->next;
  }
}

void AliTkKtJetFinder::dumpPreClusters() {
  TPreCluster *pre = firstPreCluster;
  while (pre) {
    cout << "UID: " << pre->id
	 << " E: " << pre->E
	 << " px: " << pre->px
	 << " py: " << pre->py
	 << " pz: " << pre->pz
	 << " next: " << pre->next
	 << endl;
    pre = pre->next;
  }
}

void AliTkKtJetFinder::dumpPreClusterArray() {
  for (Int_t i = 0; i < preClusterUID; i++) {
    cout << "UID: " << preClusterArray[i]->id
	 << " E: " << preClusterArray[i]->E
	 << " px: " << preClusterArray[i]->px
	 << " py: " << preClusterArray[i]->py
	 << " pz: " << preClusterArray[i]->pz
	 << " next: " << preClusterArray[i]->next
	 << endl;
  }
}

void AliTkKtJetFinder::addD(Tdlist *newD) {
    myDList.push_back(newD);
}

void AliTkKtJetFinder::buildNewDList() {
  Char_t buffer[1024];

  // step 1 - add all preclusters to dlist
  struct Tdlist *dlist;
  for (Int_t i = 0; i < preClusterUID; i++) {
    dlist = new Tdlist;
    dlist->d = calcD((TPreCluster *)preClusterArray[i]);
    dlist->id1 = preClusterArray[i]->id;
    dlist->id2 = -1;
    dlist->prev = NULL;
    dlist->next = NULL;
    addD(dlist);
  }

  sprintf(buffer,
	  "AliTkKtJetFinder::buildNewDList - added %d preclusters to new list",
	  (Int_t)myDList.size());
  DetailedOutput(buffer);

  // and now for all combinations...
  for (Int_t i = 0; i < preClusterUID; i++) {
    for (Int_t j = i + 1; j < preClusterUID; j++) {
      if (i == j) {
	// makes no sense...
	continue;
      }
      dlist = new Tdlist;
      dlist->d = calcD((TPreCluster *)preClusterArray[i],
		       (TPreCluster *)preClusterArray[j]);
      dlist->id1 = preClusterArray[i]->id;
      dlist->id2 = preClusterArray[j]->id;
      dlist->prev = NULL;
      dlist->next = NULL;
      addD(dlist);
    }
  }

  sprintf(buffer,
	  "AliTkKtJetFinder::buildNewDList - created list with %d entries",
	  (Int_t)myDList.size());
  DetailedOutput(buffer);
}

Float_t AliTkKtJetFinder::calcD(TPreCluster *p1, TPreCluster *p2) {
  if (!p1) {
    // something is wrong - we need a precluster
    return -9999;
  }
  if (!p2) {
    // called with one particle only
    // return pT^2
    return (p1->px*p1->px + p1->py*p1->py);
  } 
  // we have two particles... 
  // calculate the minimum pt...
  Float_t pt1Sq = p1->px*p1->px + p1->py*p1->py;
  Float_t pt2Sq = p2->px*p2->px + p2->py*p2->py;
  Float_t ptMinSq;
  if (pt1Sq < pt2Sq) {
    ptMinSq = pt1Sq;
  } else {
    ptMinSq = pt2Sq;
  }
  // calculate rapidities - !move this into precluster structure!!!
  Float_t y1 = 0.5 * TMath::Log((p1->E+p1->pz)/(p1->E-p1->pz));
  Float_t y2 = 0.5 * TMath::Log((p2->E+p2->pz)/(p2->E-p2->pz));
  // calculate phi - !move this into precluster structure!!!
  Float_t phi1 = 0;
  Float_t phi2 = 0;
  // calculate ((y1-y2)^2 + (phi1-phi2)^2)/(D^2)
  // !save D^2 instead of D to save time!!!
  Float_t R = ((y1-y2)*(y1-y2) + (phi1-phi2)*(phi1-phi2))/(finder_D*finder_D);
  
  return ptMinSq * R;
}

void AliTkKtJetFinder::setDebugLevel(Int_t nDebugLevel) {
  debugLevel = nDebugLevel;

  // description of debug levels
  // levels are additive!!!
  //  0 -> no debug output - only normal output (default)
  //  1 -> normal debug output
  //       major function calls, summary of major results
  //  2 -> timing output
  //       timing information for major steps
  //  4 -> detailed results
  //       output of much more results
  //  8 -> result histos
  //       histograms of event results for internal consitency checks
  // 16 -> timing histos
  //       histograms with timing info - not cleared per event!
  // 32 -> write debug histos to file
}

void AliTkKtJetFinder::DebugOutput(Char_t *output) {
  if ((debugLevel & 1) == 1) {
    cout << output << endl;
  }
}

void AliTkKtJetFinder::TimingOutput(Char_t *output) {
  if ((debugLevel & 2) == 2) {
    cout << output << endl;
  }
}

void AliTkKtJetFinder::DetailedOutput(Char_t *output) {
  if ((debugLevel & 4) == 4) {
    cout << output << endl;
  }
}


ClassImp(AliTkKtJetFinder)
