/*

       July 2008 BW mass is limited by "PYTHIA method", by I. Lokhtin and L. Malinina
                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#include <functional>
#include <algorithm>
#include <vector>
#include <iostream>
#include <TRandom.h>
#include <TError.h>
#include <TMath.h>

#ifndef DATABASE_PDG
#include "DatabasePDG.h"
#endif
#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif
#ifndef DECAY_CHANNEL
#include "DecayChannel.h"
#endif
#ifndef HADRONDECAYER_INCLUDED
#include "HadronDecayer.h"
#endif
#ifndef UKUTILITY_INCLUDED
#include "UKUtility.h"
#endif
#ifndef PARTICLE_INCLUDED
#include "Particle.h"
#endif
#include "HYJET_COMMONS.h"

//calculates decay time in fm/c
//calculates 1,2 and 3 body decays

using std::cout;
using std::endl;

Double_t GetDecayTime(const Particle &parent, Double_t weakDecayLimit) {
  ParticlePDG *pDef = parent.Def(); 
  Double_t width = pDef->GetWidth(); //GeV
  if(width > weakDecayLimit) {
    const Double_t slope =  parent.E() * 0.1973 / (pDef->GetMass() * width);
//    cout<<"get decay time"<<pDef->GetPDG()<<" "<<slope<<endl;
    return -slope * TMath::Log(gRandom->Rndm());//in fm/c
  }

  return 0.;
}


extern "C" void mydelta_();
extern SERVICEEVCommon SERVICEEV;

void Decay(List_t &output, Particle &parent, ParticleAllocator &allocator, DatabasePDG* database) {
  // check if the parent particle has been decayed already
  
  Int_t daughters = parent.GetNDaughters();

//cout<<"in Decay pdg"<< parent.Def()->GetPDG()<<"daughters "<<daughters<<endl;

  if(daughters>0)  // particle decayed already
    return;

  // Get the PDG properties of the particle
  ParticlePDG *pDef = parent.Def();

  // Get the number of posible decay channels
  Int_t nDecayChannel = pDef->GetNDecayChannels();

  // check the full branching of this specie
  Double_t fullBranching = pDef->GetFullBranching(); // Only 3 or less body decays

  // return if particle has no branching
  if(fullBranching < 0.00001)
    return;

  // get the PDG mass of the specie  
  Double_t PDGmass = pDef->GetMass();
  int ComprCodePyth=0;
  Float_t Delta =0;

  Bool_t success = kFALSE;
  Int_t iterations = 0;
  // Try to decay the particle
  while(!success) {
    // get a random mass using the Breit Wigner distribution 
    Double_t BWmass = gRandom->BreitWigner(PDGmass, pDef->GetWidth()); 
    //!!!!    
    //      BWmass = PDGmass;
    // Try to cut the Breit Wigner tail of the particle using the cuts from pythia
    // The Delta variable is obtained from pythia based on the specie
    int encoding =pDef->GetPDG();
    SERVICEEV.ipdg = encoding;
    mydelta_();      
    ComprCodePyth=SERVICEEV.KC;
    Delta = SERVICEEV.delta;// PYDAT2.PMAS[KC][3];

    //if there are no such particle in PYTHIA particle table, we take Delta=0.4
    if(ComprCodePyth==0){
      BWmass=PDGmass; 
      Delta=0.0;
    } 

    //bad delta - an exception
    if(ComprCodePyth==254){
      BWmass=PDGmass; 
      Delta=0.0;
    } 
      
    //for particles from PYTHIA table only, if the BW mass is outside the cut range then quit this iteration and generate another BW mass
    if(ComprCodePyth!=0 && Delta>0 && (BWmass<PDGmass-Delta || BWmass>PDGmass+Delta)){
      //      std::cout<<"encoding"<<encoding<<"delta"<<Delta<<"width "<<pDef->GetWidth()<<"mass"<<BWmass<<std::endl;
      continue;
    }    
    //----    
    
    if(BWmass>5)
      std::cout<<" > 5 encoding"<<encoding<<" pdgmass "<<PDGmass<<" delta "<<Delta<<"width "<<pDef->GetWidth()<<" mass "<<BWmass<<"CC"<<ComprCodePyth<<std::endl;
    
    // check how many decay channels are allowed with the generated mass
    Int_t nAllowedChannels = database->GetNAllowedChannels(pDef, BWmass);
    // if no decay channels are posible with this mass, then generate another BW mass
    if(nAllowedChannels==0) {    
      iterations++;
      continue;
    }

    std::vector<Particle> apDaughter;
    std::vector<Double_t> dMass; //daughters'mass
    std::vector<Double_t> dMom;
    std::vector<Double_t> sm;
    std::vector<Double_t> rd;

    // we need to choose an allowed decay channel
    Double_t randValue = gRandom->Rndm() * fullBranching;
    Int_t chosenChannel = 1000;
    Bool_t found = kFALSE;
    Int_t channelIterations = 0;
    while(!found) {
      for(Int_t nChannel = 0; nChannel < nDecayChannel; ++nChannel) {
	randValue -= pDef->GetDecayChannel(nChannel)->GetBranching();
	if(randValue <= 0. && database->IsChannelAllowed(pDef->GetDecayChannel(nChannel), BWmass)) {
	  chosenChannel = nChannel;
	  found = kTRUE;
	  break;
	}
      }
      channelIterations++;
    }

    // get the PDG information for the chosen decay channel
    DecayChannel *dc = pDef->GetDecayChannel(chosenChannel);
    Int_t nSec = dc->GetNDaughters();

    // Adjust the parent momentum four-vector for the MC generated Breit-Wigner mass
    Particle parentBW(database->GetPDGParticle(parent.Encoding()));
    parentBW.Pos(parent.Pos());
    Double_t BWenergy = TMath::Sqrt(parent.Mom().X()*parent.Mom().X() + 
				    parent.Mom().Y()*parent.Mom().Y() +
				    parent.Mom().Z()*parent.Mom().Z() +
				    BWmass*BWmass);

    Int_t NB = (Int_t)parent.GetType(); //particle from jets

    TLorentzVector MomparentBW(parent.Mom().X(), parent.Mom().Y(), parent.Mom().Z(), BWenergy); 
    parentBW.Mom(MomparentBW);
    // take into account BW when calculating boost velocity (for wide resonances it matters)
    TVector3 velocityBW(parentBW.Mom().BoostVector());

    // now we have an allowed decay
    // first case: one daughter particle
    if(nSec == 1) {
      // initialize the daughter particle
      Particle p1(database->GetPDGParticle(dc->GetDaughterPDG(0)));
      p1.SetLastMotherPdg(parentBW.Encoding());
      p1.SetLastMotherDecayCoor(parentBW.Pos());
      p1.SetLastMotherDecayMom(parentBW.Mom());

      // link the parent and daughters trough their indexes in the list
      Int_t parentIndex = -1;
      if(parent.GetMother()==-1) parentIndex = parent.SetIndex();   // parents which are primaries don't have yet an index
      else parentIndex = parent.GetIndex();                         // parents which are secondaries have an index
      Int_t p1Index = p1.SetIndex();                           // set the daughter index
      p1.SetMother(parentIndex);                               // set the mother index for this daughter 
      parent.SetDaughter(p1Index);                             // set p1 as daughter to the parent 
      if(parent.GetMother()==-1) allocator.AddParticle(parent, output);   // add it only if its a primary particle
      allocator.AddParticle(p1, output);
      success = kTRUE;  
    }
    // second case: two daughter particles
    else if(nSec == 2) {
      // initialize the daughter particles
      Particle p1(database->GetPDGParticle(dc->GetDaughterPDG(0)));
      p1.Pos(parentBW.Pos());
      Particle p2(database->GetPDGParticle(dc->GetDaughterPDG(1)));
      p2.Pos(parentBW.Pos());
      
      // calculate the momenta in rest frame of mother for the two particles (theta and phi are isotropic)
      MomAntiMom(p1.Mom(), p1.TableMass(), p2.Mom(), p2.TableMass(), BWmass);
    
      // boost to the laboratory system (to the mother velocity)
      p1.Mom().Boost(velocityBW);
      p2.Mom().Boost(velocityBW);

      //store information about mother
      p1.SetLastMotherPdg(parentBW.Encoding());
      p1.SetLastMotherDecayCoor(parentBW.Pos());
      p1.SetLastMotherDecayMom(parentBW.Mom());
      p2.SetLastMotherPdg(parentBW.Encoding());
      p2.SetLastMotherDecayCoor(parentBW.Pos());
      p2.SetLastMotherDecayMom(parentBW.Mom());
      // std::cout<<"2d NB="<<NB<<std::endl;
      //set to daughters the same type as has mother
      p1.SetType(NB);
      p2.SetType(NB);


      // check the kinematics in the lab system
      Double_t deltaS = TMath::Sqrt((parentBW.Mom().X()-p1.Mom().X()-p2.Mom().X())*(parentBW.Mom().X()-p1.Mom().X()-p2.Mom().X())+
				    (parentBW.Mom().Y()-p1.Mom().Y()-p2.Mom().Y())*(parentBW.Mom().Y()-p1.Mom().Y()-p2.Mom().Y())+
				    (parentBW.Mom().Z()-p1.Mom().Z()-p2.Mom().Z())*(parentBW.Mom().Z()-p1.Mom().Z()-p2.Mom().Z())+
				    (parentBW.Mom().E()-p1.Mom().E()-p2.Mom().E())*(parentBW.Mom().E()-p1.Mom().E()-p2.Mom().E()));
      // if deltaS is too big then repeat the kinematic procedure
 
 
      if(deltaS>0.001) {

	cout << "2-body decay kinematic check in lab system: " << pDef->GetPDG() << " --> " << p1.Encoding() << " + " << p2.Encoding() << endl;
	cout << "Mother    (e,px,py,pz): " << parentBW.Mom().E() << "\t" << parentBW.Mom().X() << "\t" << parentBW.Mom().Y() << "\t" << parentBW.Mom().Z() << endl;
	cout << "Mother    (x,y,z,t): " << parentBW.Pos().X() << "\t" << parentBW.Pos().Y() << "\t" << parentBW.Pos().Z() << "\t" << parentBW.Pos().T() << endl;

	cout << "Daughter1 (e,px,py,pz): " << p1.Mom().E() << "\t" << p1.Mom().X() << "\t" << p1.Mom().Y() << "\t" << p1.Mom().Z() << endl;
	cout << "Daughter2 (e,px,py,pz): " << p2.Mom().E() << "\t" << p2.Mom().X() << "\t" << p2.Mom().Y() << "\t" << p2.Mom().Z() << endl;	
	cout << "2-body decay delta(sqrtS) = " << deltaS << endl;
	cout << "Repeating the decay algorithm ..." << endl;

	iterations++;
	continue;
      }
      // push particles to the list of secondaries
      Int_t parentIndex = -1;
      if(parent.GetMother()==-1) parentIndex = parent.SetIndex();   // parents which are primaries don't have yet an index
      else parentIndex = parent.GetIndex();                         // parents which are secondaries have an index
      p1.SetIndex(); 
      p2.SetIndex();
      p1.SetMother(parentIndex); 
      p2.SetMother(parentIndex);
      parent.SetDaughter(p1.GetIndex());
      parent.SetDaughter(p2.GetIndex());
      if(parent.GetMother()==-1) allocator.AddParticle(parent, output);      // add it only if its a primary
      allocator.AddParticle(p1, output);
      allocator.AddParticle(p2, output);
      success = kTRUE;
    }

    // third case: three daughter particle
    else if(nSec == 3) {
      // initialize the daughter particle
      Particle p1(database->GetPDGParticle(dc->GetDaughterPDG(0)));
      p1.Pos(parentBW.Pos());
      Particle p2(database->GetPDGParticle(dc->GetDaughterPDG(1)));
      p2.Pos(parentBW.Pos());
      Particle p3(database->GetPDGParticle(dc->GetDaughterPDG(2)));
      p3.Pos(parentBW.Pos());
      // calculate the momenta in the rest frame of the mother particle
      Double_t pAbs1 = 0., pAbs2 = 0., pAbs3 = 0., sumPabs = 0., maxPabs = 0.;
      Double_t mass1 = p1.TableMass(), mass2 = p2.TableMass(), mass3 = p3.TableMass();
      TLorentzVector &mom1 = p1.Mom(), &mom2 = p2.Mom(), &mom3 = p3.Mom(); 
      Double_t deltaMass = BWmass - mass1 - mass2 - mass3;

      do {
	Double_t rd1 = gRandom->Rndm();
	Double_t rd2 = gRandom->Rndm();
	if (rd2 > rd1)
	  std::swap(rd1, rd2);
	// 1
	Double_t e = rd2*deltaMass;
	pAbs1 = TMath::Sqrt(e*e + 2*e*mass1);
	sumPabs = pAbs1;
	maxPabs = sumPabs;
	// 2
	e = (1-rd1)*deltaMass;
	pAbs2 = TMath::Sqrt(e*e + 2*e*mass2);
	
	if(pAbs2 > maxPabs)
	  maxPabs = pAbs2;
	
	sumPabs += pAbs2;
	// 3
	e = (rd1-rd2)*deltaMass;
	pAbs3 = TMath::Sqrt(e*e + 2*e*mass3);
	
	if (pAbs3 > maxPabs)
	  maxPabs =  pAbs3;
	sumPabs  +=  pAbs3;
      } while(maxPabs > sumPabs - maxPabs);
      
      // isotropic sample first particle 3-momentum
      Double_t cosTheta = 2*(gRandom->Rndm()) - 1;
      Double_t sinTheta = TMath::Sqrt(1 - cosTheta*cosTheta);
      Double_t phi      = TMath::TwoPi()*(gRandom->Rndm());
      Double_t sinPhi   = TMath::Sin(phi);
      Double_t cosPhi   = TMath::Cos(phi);
      
      mom1.SetPxPyPzE(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta, 0);
      mom1 *= pAbs1;
      // sample rest particle 3-momentum
      Double_t cosThetaN = (pAbs2*pAbs2 - pAbs3*pAbs3 - pAbs1*pAbs1)/(2*pAbs1*pAbs3);
      Double_t sinThetaN = TMath::Sqrt(1 - cosThetaN*cosThetaN);
      Double_t phiN      = TMath::TwoPi()*(gRandom->Rndm());
      Double_t sinPhiN   = TMath::Sin(phiN);
      Double_t cosPhiN   = TMath::Cos(phiN);
      
      mom3.SetPxPyPzE(sinThetaN*cosPhiN*cosTheta*cosPhi - sinThetaN*sinPhiN*sinPhi + cosThetaN*sinTheta*cosPhi,
		      sinThetaN*cosPhiN*cosTheta*sinPhi + sinThetaN*sinPhiN*cosPhi + cosThetaN*sinTheta*sinPhi,
		      -sinThetaN*cosPhiN*sinTheta + cosThetaN*cosTheta,
		      0.);
      
      mom3 *= pAbs3*mom3.P();
      mom2 = mom1;
      mom2 += mom3;
      mom2 *= -1.;
      // calculate energy
      mom1.SetE(TMath::Sqrt(mom1.P()*mom1.P() + mass1*mass1));
      mom2.SetE(TMath::Sqrt(mom2.P()*mom2.P() + mass2*mass2));
      mom3.SetE(TMath::Sqrt(mom3.P()*mom3.P() + mass3*mass3));
      
      // boost to Lab system
      mom1.Boost(velocityBW);
      mom2.Boost(velocityBW);
      mom3.Boost(velocityBW);
      
      p1.SetLastMotherPdg(parentBW.Encoding());
      p1.SetLastMotherDecayCoor(parentBW.Pos());
      p1.SetLastMotherDecayMom(parentBW.Mom());
      p2.SetLastMotherPdg(parentBW.Encoding());
      p2.SetLastMotherDecayCoor(parentBW.Pos());
      p2.SetLastMotherDecayMom(parentBW.Mom());
      p3.SetLastMotherPdg(parentBW.Encoding());
      p3.SetLastMotherDecayCoor(parentBW.Pos());
      p3.SetLastMotherDecayMom(parentBW.Mom());

      //set to daughters the same type as has mother  
      p1.SetType(NB);
      p2.SetType(NB);
      p3.SetType(NB);
 //      std::cout<<"3d NB="<<NB<<std::endl;


            
      // energy conservation check in the lab system
      Double_t deltaS = TMath::Sqrt((parentBW.Mom().X()-p1.Mom().X()-p2.Mom().X()-p3.Mom().X())*(parentBW.Mom().X()-p1.Mom().X()-p2.Mom().X()-p3.Mom().X()) +
				    (parentBW.Mom().Y()-p1.Mom().Y()-p2.Mom().Y()-p3.Mom().Y())*(parentBW.Mom().Y()-p1.Mom().Y()-p2.Mom().Y()-p3.Mom().Y()) +
				    (parentBW.Mom().Z()-p1.Mom().Z()-p2.Mom().Z()-p3.Mom().Z())*(parentBW.Mom().Z()-p1.Mom().Z()-p2.Mom().Z()-p3.Mom().Z())	+
				    (parentBW.Mom().E()-p1.Mom().E()-p2.Mom().E()-p3.Mom().E())*(parentBW.Mom().E()-p1.Mom().E()-p2.Mom().E()-p3.Mom().E()));
      // if deltaS is too big then repeat the kinematic procedure
      if(deltaS>0.001) {

	cout << "3-body decay kinematic check in lab system: " << pDef->GetPDG() << " --> " << p1.Encoding() << " + " << p2.Encoding() << " + " << p3.Encoding() << endl;
	cout << "Mother    (e,px,py,pz): " << parentBW.Mom().E() << "\t" << parentBW.Mom().X() << "\t" << parentBW.Mom().Y() << "\t" << parentBW.Mom().Z() << endl;
	cout << "Daughter1 (e,px,py,pz): " << p1.Mom().E() << "\t" << p1.Mom().X() << "\t" << p1.Mom().Y() << "\t" << p1.Mom().Z() << endl;
	cout << "Daughter2 (e,px,py,pz): " << p2.Mom().E() << "\t" << p2.Mom().X() << "\t" << p2.Mom().Y() << "\t" << p2.Mom().Z() << endl;
	cout << "Daughter3 (e,px,py,pz): " << p3.Mom().E() << "\t" << p3.Mom().X() << "\t" << p3.Mom().Y() << "\t" << p3.Mom().Z() << endl;
	cout << "3-body decay delta(sqrtS) = " << deltaS << endl;
	cout << "Repeating the decay algorithm..." << endl;

	iterations++;
	continue;
      }

      Int_t parentIndex = -1;
      if(parent.GetMother()==-1) parentIndex = parent.SetIndex();   // parents which are primaries don't have yet an index
      else parentIndex = parent.GetIndex();                         // parents which are secondaries have an index
      p1.SetIndex();
      p2.SetIndex();
      p3.SetIndex();
      p1.SetMother(parentIndex); 
      p2.SetMother(parentIndex);
      p3.SetMother(parentIndex);
      parent.SetDaughter(p1.GetIndex());
      parent.SetDaughter(p2.GetIndex());
      parent.SetDaughter(p3.GetIndex());
      if(parent.GetMother()==-1) allocator.AddParticle(parent, output);    // add it only if its a primary
      allocator.AddParticle(p1, output);
      allocator.AddParticle(p2, output);
      allocator.AddParticle(p3, output);
      success = kTRUE;
    }
  }
  return;
}
