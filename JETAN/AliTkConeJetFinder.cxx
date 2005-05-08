//$Id$
/**************************************************************************
 * AliTkConeJetFinder.cxx                                                    *
 * Thorsten Kollegger <kollegge@ikf.physik.uni-frankfurt.de               *
 * Jet finder based on a cone algorithm with seeds and addition of        *
 * midpoints, for a description of the algorithm see hep-ex/0005012       *
 *************************************************************************/

// all includes in AliTkConeJetFinder.h
#include "AliTkConeJetFinder.h"

ClassImp(AliTkConeJetFinder)

AliTkConeJetFinder::AliTkConeJetFinder() : TObject() {
  setEtaGrid(-999,999,-999);
  setPhiGrid(-999,999,-999);
  hEtHist = NULL;
}

AliTkConeJetFinder::~AliTkConeJetFinder() {
  if (hEtHist) {
    delete hEtHist;
  }
}

void AliTkConeJetFinder::setDefaultSettings() {
  setEtaGrid(100,-5.,5.);
  setPhiGrid((Int_t) (2*TMath::Pi()/0.1),
	     0,
	     2*TMath::Pi());
  setJetConeRadius(0.7);
}

void AliTkConeJetFinder::Init() {
  // create an array for the tower info
  // and at this stage we know how many towers we have...
  towers = NULL;
  nTower = nEtaBins*nPhiBins;
  towers = new AliTkTower[nTower];
  if (!towers) {
    cerr << "Couldn't allocate space for tower info!!!" << endl;
  }
  // let's create the towers;
  Float_t etaBinSize = (nEtaMax-nEtaMin)/(Float_t)nEtaBins;
  cout << "Init: etaBinSize=" << etaBinSize << endl;
  Float_t phiBinSize = (nPhiMax-nPhiMin)/(Float_t)nPhiBins;
  cout << "Init: phiBinSize=" << phiBinSize << endl;
  //AliTkTower *tower;
  Float_t eta = nEtaMin;
  Float_t phi = nPhiMin;
  for (Int_t uid = 0; uid < nTower; uid++) {
    towers[uid].uid = uid;
    towers[uid].etaMin = eta;
    towers[uid].phiMin = phi;
    towers[uid].etaMax = eta + etaBinSize;
    towers[uid].phiMax = phi + phiBinSize;
    towers[uid].Et = 0;
    phi += phiBinSize;
    if (phi > nPhiMax) {
      phi = nPhiMin;
      eta += etaBinSize;
    }
  }

  // create an array for the seed list...
  SeedTower = new Int_t[nTower];

  // create the array to store the protojets
  const Int_t expectedNumberOfProtoJets = 1000;
  protojets = new TObjArray(expectedNumberOfProtoJets);
  protojets->SetOwner(kTRUE);

  // called once to initalize the finder...
  // create the histogram which stores the Et information
  hEtHist = new TH2D("Tower_EtHist","Tower_EtHist",
		     nEtaBins,nEtaMin,nEtaMax,
		     nPhiBins,nPhiMin,nPhiMax);
  hEtHist->GetXaxis()->SetTitle("#eta");
  hEtHist->GetYaxis()->SetTitle("#phi (rad)");
  hEtHist->GetZaxis()->SetTitle("E_{t} in bins (GeV)");

  // create histogram to hold Et cone centeres around tower center
  hEtConeHist = new TH2D("Cone_EtHist","Cone_EtHist",
			 nEtaBins,nEtaMin,nEtaMax,
			 nPhiBins,nPhiMin,nPhiMax);
  hEtConeHist->GetXaxis()->SetTitle("#eta");
  hEtConeHist->GetYaxis()->SetTitle("#phi (rad)");
  hEtConeHist->GetZaxis()->SetTitle("E_{t} in cones (GeV)");

  // create histogram to show the stable protojets before merging
  hStableProtoJetHist = new TH2D("StableProtoJetHist","StableProtoJetHist",
				 nEtaBins,nEtaMin,nEtaMax,
				 nPhiBins,nPhiMin,nPhiMax);
  hStableProtoJetHist->GetXaxis()->SetTitle("#eta");
  hStableProtoJetHist->GetYaxis()->SetTitle("#phi (rad)");
  hStableProtoJetHist->GetZaxis()->SetTitle("E_{t} in stable cones (GeV)");

  hTowerEt = new TH2D("tower_EtCheck","tower EtCheck",
			nEtaBins,nEtaMin,nEtaMax,
			nPhiBins,nPhiMin,nPhiMax);
  hTowerEt->GetXaxis()->SetTitle("#eta");
  hTowerEt->GetYaxis()->SetTitle("#phi (rad)");
  hTowerEt->GetZaxis()->SetTitle("E_{t} in bins (GeV)");

}

void AliTkConeJetFinder::InitEvent() {
  // called at each events
  // clears the jet finder for the next event

  // clear tower array - only Et
  for (Int_t i = 0; i < nTower; i++) {
    towers[i].Et = 0.;
  }

  // clear protojet list
  protojets->Clear();

  // clear control histograms
  hEtHist->Reset();
  hEtConeHist->Reset();
  hTowerEt->Reset();
  hStableProtoJetHist->Reset();
  
}

void AliTkConeJetFinder::FillEtHistFromTParticles(TClonesArray *particles) {
  // input TClonesArray with TParticles, e.g. from PYTHIA
  // fills the Et grid from these particles
  if (!particles) {
    return;
  }
  TParticle *particle = NULL;
  Float_t Et = 0;

  // loop over all particles...
  TIterator *iter = particles->MakeIterator();
  while ((particle = (TParticle *) iter->Next()) != NULL) {
    // check if particle is accepted
    if (isTParticleAccepted(particle)) {
      // calculate Et for a massless particle...
      Et = TMath::Sqrt(particle->Pt()*particle->Pt());
      // idea for Et check: Et = E *sin(theta);
      Int_t tower = findTower(particle->Eta(),particle->Phi());
      if (tower >= 0) {
	towers[tower].Et += Et;
      } else {
// 	cerr << "Couldn't find tower!!! eta="
// 	     << particle->Eta() << " phi=" << particle->Phi()
// 	     << endl;
      }
      // and to the histogram for control reasons
      AddEtHist(particle->Eta(),particle->Phi(),Et);
    }
  }
  
  // at this point we have filled all towers...
  // let's make sure tower info is right
  for (Int_t i = 0; i < nTower; i++) {
    Float_t etacenter = (towers[i].etaMin + towers[i].etaMax)/2.;
    Float_t phicenter = (towers[i].phiMin + towers[i].phiMax)/2.;
    hTowerEt->Fill(etacenter,phicenter,towers[i].Et);
  }

}

Bool_t AliTkConeJetFinder::isTParticleAccepted(TParticle *particle) {
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

void AliTkConeJetFinder::AddEtHist(Float_t eta, Float_t phi, Double_t Et) {
  // add particle/tower withEt at position eta,phi to the current Et grid...
  if (hEtHist) {
    Double_t oldEt;
    // get the Et already stored in the bin eta,phi
    oldEt = hEtHist->GetBinContent(hEtHist->GetXaxis()->FindFixBin(eta),
				   hEtHist->GetYaxis()->FindFixBin(phi));
    // add the new one
    Et += oldEt;
    // and set the bin to the sum of old+new Et.
    hEtHist->SetBinContent(hEtHist->GetXaxis()->FindFixBin(eta),
			   hEtHist->GetYaxis()->FindFixBin(phi),
			   Et);
  }
}

Int_t AliTkConeJetFinder::run() {
  Int_t status;
  status = FindProtoJets();
  cout << "Found " << protojets->GetEntries() << " stable protojets"
       << endl;
  return status;
}

Int_t AliTkConeJetFinder::FindProtoJets() {
  Float_t etaTower;
  Float_t phiTower;

  Float_t eta;
  Float_t phi;
  Float_t Et;

  Float_t etaDiff;
  Float_t phiDiff;

  AliTkProtoJet protojet;

  // let's take each tower as seed
  // this means seedless search...
  // introduce later step with seed+midpoints...
  for (Int_t i=0; i < nTower; i++) {
    // loop over all towers
    etaTower = (towers[i].etaMax+towers[i].etaMin)/2.;
    phiTower = (towers[i].phiMax+towers[i].phiMin)/2.;

    eta = etaTower;
    phi = phiTower;

    // let's calculate a protojet with centered around tower center
    if ((CalculateConeCentroid(&eta,&phi,&Et) > 0)) {
      // ok, we found a valid centroid around tower pos
      // let's fill the initial seed into the EtConeHist
      hEtConeHist->Fill(etaTower,phiTower,Et);
      // let's calculate the vector to the weighted cone center
      etaDiff = (eta - etaTower)*(eta - etaTower);
      phiDiff = calcPhiDiff(phi,phiTower);

      // let's iterate the proto-jet while it's in the tower boundaries
      // or it is stable
      // or the maximum number of iterations is reached
      Bool_t isStable = kFALSE;
      Bool_t isInTower = isProtoJetInTower(i,eta,phi);
      //Bool_t isInTower = isProtoJetInTower(etaBin,phiBin,eta,phi);
      Int_t nIterations = 0;
      Float_t etaOld;
      Float_t phiOld;
      
      // make const variables private members and set in default...
      const Int_t nMaxIterations = 100;
      while ((nIterations < nMaxIterations) &&
	     (!isStable) &&
	     (isInTower)) {
	etaOld = eta;
	phiOld = phi;
	CalculateConeCentroid(&eta,&phi,&Et);
	etaDiff = (eta - etaOld)*(eta - etaOld);
	phiDiff = calcPhiDiff(phi,phiOld)*calcPhiDiff(phi,phiOld);
	// and check protojet again...
	isStable = isProtoJetStable(etaDiff,phiDiff);
	// we allow the jet to flow away...
	//isInTower = kTRUE;
	// we force protojets to stay in tower
	isInTower = isProtoJetInTower(i,eta,phi);
	nIterations++;
      }

      if ((isStable) || (nIterations == nMaxIterations)) {
// 	cout << "Found a stable protojet at eta=" << eta
// 	     << " phi=" << phi << " with Et=" << Et << endl;
	if (nIterations == nMaxIterations) {
	  cout << "To many iterations on protojet" << endl;
	}
	// create a new protojet object
	AliTkProtoJet *protojet = new AliTkProtoJet();
	protojet->eta = eta;
	protojet->phi = phi;
	protojet->Et = Et;

	// let's check if we allready found this protojet
	TIter *iter = new TIter(protojets);
	AliTkProtoJet *pj;
	Bool_t wasFound = kFALSE;
	while ((pj = (AliTkProtoJet *) iter->Next()) && !wasFound) {
	  wasFound = protojet->IsEqual(pj);
	}

	// and add it to the list of protojets...
	if (!wasFound) {
	  protojets->Add(protojet);
	  //cout << protojet->Et << endl;
	  hStableProtoJetHist->Fill(eta,phi,Et);
	} else {
	  delete protojet;
	}
      }
    }
  }
  return 0;
}

Int_t AliTkConeJetFinder::FindJets() {
  //let's look for the highest Et-Jet
  //AliTkProtoJet *pj;

  return 0;
}


Int_t AliTkConeJetFinder::CalculateConeCentroid(Float_t *eta, Float_t *phi,
					     Float_t *Et) {

  //====================================================================
  // WORKS STILL ON HISTOGRAM
  // MUST WORK ON TOWERS
  //====================================================================


  // let's copy some stuff
  Float_t etaCenter = *eta;
  Float_t phiCenter = *phi;
  
  /*
  cout << "CalculateConeCentroid: Initial Cone Center eta="
       << etaCenter 
       << " phi="
       << phiCenter
       << endl;
  */

  Float_t etaTower = 0;
  Float_t phiTower = 0;
  Float_t EtTower = 0;

  Float_t etaDiff = 0;
  Float_t phiDiff = 0;

  Float_t etaJet = 0;
  Float_t phiJet = 0;
  Float_t EtJet = 0;

  // lets calculate boundaries from jet cone
  Float_t etaMin = etaCenter - jetRadius;
  Float_t etaMax = etaCenter + jetRadius;
  Float_t phiMin = phiCenter - jetRadius;

  if (phiMin < 0) {
    phiMin = 2*TMath::Pi() + phiMin;
  }
  Float_t phiMax = phiCenter + jetRadius;
  if (phiMax > 2*TMath::Pi()) {
    phiMax = phiMax - 2*TMath::Pi();
  }

  Int_t etaBinMin = hEtHist->GetXaxis()->FindFixBin(etaMin);
  if (etaBinMin == 0) {
    return -1;
  }
  Int_t etaBinMax = hEtHist->GetXaxis()->FindFixBin(etaMax);
  if (etaBinMax == hEtHist->GetXaxis()->GetNbins()+1) {
    return -1;
  }
  Int_t phiBinMin = hEtHist->GetYaxis()->FindFixBin(phiMin);
  Int_t phiBinMax = hEtHist->GetYaxis()->FindFixBin(phiMax);
 
  // let's loop over all bins/towers...
  for (Int_t etaBin = etaBinMin; etaBin != etaBinMax+1; etaBin++) {
    for (Int_t phiBin = phiBinMin; phiBin != phiBinMax+1; phiBin++) {
      // roll over in phi...
      if (phiBin == hEtHist->GetYaxis()->GetNbins()+1) {
	phiBin = 0;
	continue;
      }
      // let's get the tower
      etaTower = hEtHist->GetXaxis()->GetBinCenter(etaBin);
      phiTower = hEtHist->GetYaxis()->GetBinCenter(phiBin);
      EtTower = hEtHist->GetCellContent(etaBin,phiBin);

      // let's check if the bin is in the cone
      etaDiff = (etaTower - etaCenter)*(etaTower - etaCenter);
      phiDiff = calcPhiDiff(phiTower,phiCenter)*
	calcPhiDiff(phiTower,phiCenter);
      if (TMath::Sqrt(etaDiff+phiDiff) < jetRadius) {
	// tower is in cone - add it to the cone
	EtJet += EtTower;
	if (EtTower != 0) {
	  etaJet += EtTower * etaTower;
	  phiJet += EtTower * phiTower;
	}
      }
    }
  }
  
  // normalize eta and phi
  if (EtJet != 0) {
    etaJet /= EtJet;
    phiJet /= EtJet;
  }

  /*
  cout << "CalculateConeCentroid: Final Cone Center eta="
       << etaJet
       << " phi="
       << phiJet
       << " Et="
       << EtJet
       << endl;
  */

  // return the jet cone parameters
  *eta = etaJet;
  *phi = phiJet;
  *Et = EtJet;
  

  return 1;
}

Float_t AliTkConeJetFinder::calcPhiDiff(Float_t phi1, Float_t phi2) {
  Float_t phidiff1 = TMath::Sqrt((phi1-phi2)*(phi1-phi2));
  Float_t phidiff2 = TMath::Sqrt((phi1-phi2-2*TMath::Pi())*
				 (phi1-phi2-2*TMath::Pi()));
  Float_t phidiff3 = TMath::Sqrt((phi1-phi2+2*TMath::Pi())*
				 (phi1-phi2+2*TMath::Pi()));
  // return the minmum;
  Float_t phiret = phidiff1;
  if (phiret > phidiff2) {
    phiret = phidiff2;
  }
  if (phiret > phidiff3) {
    phiret = phidiff3;
  }
  return phiret;
}

Bool_t AliTkConeJetFinder::isProtoJetStable(Float_t etaDiff, Float_t phiDiff) {
  const Float_t epsilon = 0.00001;
  if (TMath::Sqrt(etaDiff*phiDiff) < epsilon) {
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliTkConeJetFinder::isProtoJetInTower(Int_t etaBin, Int_t  phiBin,
					  Float_t eta, Float_t phi) {
  Float_t etaBinMax = hEtHist->GetXaxis()->GetBinUpEdge(etaBin);
  Float_t etaBinMin = hEtHist->GetXaxis()->GetBinLowEdge(etaBin);
  Float_t phiBinMax = hEtHist->GetYaxis()->GetBinUpEdge(phiBin);
  Float_t phiBinMin = hEtHist->GetYaxis()->GetBinLowEdge(phiBin);
  if ((eta < etaBinMax) &&
      (eta > etaBinMin) &&
      (phi < phiBinMax) &&
      (phi > phiBinMin)) {
    return kTRUE;
  }

  return kFALSE;
}

Bool_t AliTkConeJetFinder::isProtoJetInTower(Int_t tower, Float_t eta, 
					  Float_t phi) {
  if ((eta > towers[tower].etaMin) &&
      (eta < towers[tower].etaMax) &&
      (phi > towers[tower].phiMin) &&
      (phi < towers[tower].phiMax)) {
    return kTRUE;
  }
  return kFALSE;
}

Int_t AliTkConeJetFinder::findTower(Float_t eta, Float_t phi) {
  Float_t etaBinSize = (nEtaMax - nEtaMin)/(Float_t)nEtaBins;
  Float_t phiBinSize = (nPhiMax - nPhiMin)/(Float_t)nPhiBins;

  if ((eta < nEtaMin) || (eta > nEtaMax)) {
    return -1;
  }
  if ((phi < nPhiMin) || (phi > nPhiMax)) {
    return -1;
  }

  Int_t etaBin = (Int_t) ((eta - nEtaMin)/etaBinSize);
  Int_t phiBin = (Int_t) ((phi - nPhiMin)/phiBinSize);

  Int_t bin = etaBin * nPhiBins + phiBin;

  if (bin >= nTower) {
    return -2;
  }
 
  return bin;
}

// lot of setters/getters... nothing special, no comments...
void AliTkConeJetFinder::setEtaNBins(Int_t nbins) {
  nEtaBins = nbins;
}

Int_t AliTkConeJetFinder::getEtaNBins() {
  return nEtaBins;
}

void AliTkConeJetFinder::setEtaRange(Float_t min, Float_t max) {
  nEtaMin = min;
  nEtaMax = max;
}

Float_t AliTkConeJetFinder::getEtaRangeMin() {
  return nEtaMin;
}

Float_t AliTkConeJetFinder::getEtaRangeMax() {
  return nEtaMax;
}

void AliTkConeJetFinder::setEtaGrid(Int_t nbins, Float_t min, Float_t max) {
  setEtaNBins(nbins);
  setEtaRange(min,max);
}

void AliTkConeJetFinder::setPhiNBins(Int_t nbins) {
  nPhiBins = nbins;
}

Int_t AliTkConeJetFinder::getPhiNBins() {
  return nPhiBins;
}

void AliTkConeJetFinder::setPhiRange(Float_t min, Float_t max) {
  nPhiMin = min;
  nPhiMax = max;
}

Float_t AliTkConeJetFinder::getPhiRangeMin() {
  return nPhiMin;
}

Float_t AliTkConeJetFinder::getPhiRangeMax() {
  return nPhiMax;
}

void AliTkConeJetFinder::setPhiGrid(Int_t nbins, Float_t min, Float_t max) {
  setPhiNBins(nbins);
  setPhiRange(min,max);
}

void AliTkConeJetFinder::setJetConeRadius(Float_t r) {
  jetRadius = r;
}

Float_t AliTkConeJetFinder::getJetConeRadius() {
  return jetRadius;
}

TObjArray *AliTkConeJetFinder::getProtoJetList() {
  return protojets;
}

//==========================================================================
ClassImp(AliTkTower)

//==========================================================================
ClassImp(AliTkProtoJet)

Bool_t AliTkProtoJet::IsEqual(AliTkProtoJet *other) {
  const Float_t epsilon = 0.0000001;
  
  Float_t etaDiff = (other->eta - eta)*(other->eta - eta);
  Float_t phiDiff = (other->phi - phi)*(other->phi - phi);
  Float_t EtDiff = (other->Et - Et)*(other->Et - Et);

  Bool_t isEqual = kFALSE;
  if ((etaDiff < epsilon) &&
      (phiDiff < epsilon) &&
      (EtDiff < epsilon)) {
    isEqual = kTRUE;
  }
  return isEqual;
}

bool SProtoJet::operator==(const SProtoJet &s1) {
  Float_t epsilon = 0.0000001;
  bool b;
  if (eta > s1.eta) {
    b = ((eta - s1.eta) < epsilon);
  } else {
    b = ((s1.eta - eta) < epsilon);
  }
  if (phi > s1.phi) {
    b = b && ((phi - s1.phi) < epsilon);
  } else {
    b = b && ((s1.phi - phi) < epsilon);
  }
  if (Et > s1.Et) {
    b = b && ((Et - s1.Et) < epsilon);
  } else {
    b = b && ((s1.Et - Et) < epsilon);
  }
  return b;
}
