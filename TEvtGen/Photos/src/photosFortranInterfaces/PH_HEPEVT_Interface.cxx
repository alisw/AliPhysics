#include <vector>
#include <cmath>
#include "PhotosBranch.h"
#include "PhotosParticle.h"
#include "PH_HEPEVT_Interface.h"
#include "Log.h"
using namespace std;

namespace Photospp
{

vector<PhotosParticle*> PH_HEPEVT_Interface::m_particle_list;
int PH_HEPEVT_Interface::ME_channel=0;
int PH_HEPEVT_Interface::decay_idx=0;

void PH_HEPEVT_Interface::clear(){

  m_particle_list.clear();

  hep.nevhep=0; 
  hep.nhep=0;
  

  /**  for(int i=0; i < NMXHEP; i++){

    hep.isthep[i]=0;
    hep.idhep[i]=0;
    
    for(int j=0; j<2; j++){
      hep.jmohep[i][j]=0;
      hep.jdahep[i][j]=0;
    }
    
    for(int j=0; j<5; j++)
      hep.phep[i][j]=0;
    
    for(int j=0; j<4; j++)
      hep.vhep[i][j]=0;
  
      ph_phoqed_.qedrad[i]=0;
  
      }**/
}

void PH_HEPEVT_Interface::add_particle(int i,PhotosParticle * particle,
				       int first_mother, int last_mother,
				       int first_daughter, int last_daughter){

  if(i>0)
    i--; //account for fortran indicies begining at 1
  else
    Log::Warning()<<"Index given to PH_HEPEVT_Interface::add_particle "
	          <<"is too low (it must be > 0)."<<endl;

  //add to our internal list of pointer/index pairs
  m_particle_list.push_back(particle);

  //now set the element of PH_HEPEVT
  hep.nevhep=0; //dummy
  hep.nhep = hep.nhep + 1; // 1.II.2014: fix for gcc 4.8.1. Previously it was:
                           // hep.nhep = hep.nhep++; which technically is an undefined operation
  hep.isthep[i]=particle->getStatus();
  hep.idhep[i]=particle->getPdgID();

  hep.jmohep[i][0]=first_mother;
  hep.jmohep[i][1]=last_mother;

  hep.jdahep[i][0]=first_daughter;
  hep.jdahep[i][1]=last_daughter;

  hep.phep[i][0]=particle->getPx();
  hep.phep[i][1]=particle->getPy();
  hep.phep[i][2]=particle->getPz();
  hep.phep[i][3]=particle->getE();
  
  // if massFrom4Vector=true (default) - get sqrt(e^2-p^2)
  // otherwise - get mass from event record
  if(!Photos::massFrom4Vector) hep.phep[i][4]=particle->getMass();
  else                         hep.phep[i][4]=particle->getVirtuality();

  int pdgid = abs(particle->getPdgID());

  // if 'forceMass' for this PDGID was used - overwrite mass
  if(Photos::forceMassList)
  {
    for(unsigned int j=0;j<Photos::forceMassList->size();j++)
    {
      if(pdgid == abs(Photos::forceMassList->at(j)->first))
      {
        double mass = Photos::forceMassList->at(j)->second;
        
        // when 'forceMass' is used the mass provided is larger than 0.0
        // when 'forceMassFromEventRecord' is used mass is -1.0
        // in this case - get mass from event record
        if(mass<0.0) mass = particle->getMass();
        hep.phep[i][4] = mass;
      }
    }
  }

  hep.vhep[i][0]=0;
  hep.vhep[i][1]=0;
  hep.vhep[i][2]=0;
  hep.vhep[i][3]=0;

  hep.qedrad[i]=1;

}

int PH_HEPEVT_Interface::set(PhotosBranch *branch)
{
	PH_HEPEVT_Interface::clear();
	int idx=1;

	//get mothers
	vector<PhotosParticle *> mothers = branch->getMothers();
	int nmothers=mothers.size();

	//check if mid-particle exist
	decay_idx=0;
	PhotosParticle *decay_particle = branch->getDecayingParticle();
	if(decay_particle) decay_idx=nmothers+1;

	//get daughters
	vector<PhotosParticle *> daughters = branch->getDaughters();
	int ndaughters=daughters.size();

	for(int i=0;i<nmothers;i++)
	{
		if(decay_idx)
			add_particle(idx++,mothers.at(i),
			             0,0, //mothers
			             decay_idx,decay_idx); //daughters
		else
			add_particle(idx++,mothers.at(i),
			             0,0, //mothers
			             nmothers+1,nmothers+ndaughters); //daughters
	}

	if(decay_particle)
		add_particle(idx++,decay_particle,
		             1,nmothers, //mothers
		             nmothers+2,nmothers+1+ndaughters); //daughters

	for(int i=0;i<ndaughters;i++)
	{
		if(decay_idx)
			add_particle(idx++,daughters.at(i),
			             decay_idx,decay_idx, //mothers
			             0,0); //daughters
		else
			add_particle(idx++,daughters.at(i),
			             1,nmothers, //mothers
			             0,0); //daughters
	}
	//Log::RedirectOutput( phodmp_ , Log::Debug(1000) );
	Log::Debug(1000,false)<<"PH_HEPEVT returning: "<<( (decay_idx) ? decay_idx : 1 )<<" from "<<idx-1<<" particles."<<endl;
	return (decay_idx) ? decay_idx : 1;
}

void PH_HEPEVT_Interface::get(){

  int index = 0;

  //if no photons have been added to the event record, do nothing.
  if(hep.nhep == (int) m_particle_list.size())
    return;

  //phodmp_();

  int  particle_count  = m_particle_list.size();
  int  daughters_start = hep.jmohep[hep.nhep-1][0];
  int  photons         = hep.nhep - m_particle_list.size();
  bool isPhotonCreated = (photons>0);
  
  std::vector<PhotosParticle*> photon_list; // list of added photons
                                            // which need kinematical treatment
                                            // in special case

  // we decipher daughters_start from  last entry 
  // that is last daughter in  ph_hepevt_
  // another option of this functionality may be 
  // hep.jdahep[ hep.jmohep[hep.nhep-1][0]-1][0];
  // Update daughters_start if there are two mothers
  // NOTE: daughters_start is index for C++ arrays, while hep.jmohep
  //       contains indices for Fortran arrays.
  if(hep.jmohep[hep.nhep-1][1]>0)
    daughters_start = hep.jmohep[hep.nhep-1][1];
  
  index = particle_count;

  // Add extra photons
  for(;photons>0; photons--, index++){
    
    if(hep.idhep[index]!=PhotosParticle::GAMMA)
      Log::Fatal("PH_HEPEVT_Interface::get(): Extra particle added to the PH_HEPEVT common block in not a photon!",6);
    
    //create a new particle
    PhotosParticle * new_photon;
    new_photon = m_particle_list.at(0)->createNewParticle(hep.idhep[index],
							  hep.isthep[index],
							  hep.phep[index][4],
							  hep.phep[index][0],
							  hep.phep[index][1],
							  hep.phep[index][2],
							  hep.phep[index][3]);
    
    //add into the event record
    //get mother particle of photon
    PhotosParticle * mother =  m_particle_list.at(hep.jmohep[index][0]-1);
    mother->addDaughter(new_photon);
    
    //add to list of photons
    photon_list.push_back(new_photon);
  }

  // Before we update particles, we check for special cases
  // At this step, particles are yet unmodified
  // but photons are already in the event record
  bool special=false;
  PhotosParticle *p1 = NULL;
  PhotosParticle *p2 = NULL;

  if( isPhotonCreated )
  {
    std::vector<PhotosParticle*> daughters;

    // in the following we create list of   daughters,
    // later  we calculate bool special which is true only if all
    // daughters self-decay
    // at peresent warning for  mixed self-decay and not self decay 
    // daughters is not printed.

    for(int i=daughters_start;i<particle_count;i++)
    {
      PhotosParticle *p = m_particle_list.at(i);

      daughters.push_back(p);
    }

    // Check if this is a special case
    special = true;
    
    if(daughters.size()==0) special = false;
    
    // special = false if there is a stable particle on the list
    //                 or there is a particle without self-decay
    for(unsigned int i=0;i<daughters.size();i++)
    {
      if(daughters[i]->getStatus()==1)
      {
        special = false;
        break;
      }
      
      // NOTE: We can use 'getDaughters' here, because vertices
      //       of daughters are not being modified by Photos right now
      //       (so there will be no caching)
      std::vector<PhotosParticle*> daughters2 = daughters[i]->getDaughters();
      
      if(daughters2.size()!=1 || 
         daughters2[0]->getPdgID() != daughters[i]->getPdgID() )
      {
        special = false;
        break;
      }
    }

    if( special )
    {
      double px1=0.0, py1=0.0, pz1=0.0, e1=0.0;
      double px2=0.0, py2=0.0, pz2=0.0, e2=0.0;

      // get sum of 4-momenta of unmodified particles
      for(unsigned int i=0;i<daughters.size();i++)
      {
        // ignore photons
        if(daughters[i]->getPdgID()==22) continue;

        px1+=daughters[i]->getPx();
        py1+=daughters[i]->getPy();
        pz1+=daughters[i]->getPz();
        e1 +=daughters[i]->getE();
      }

      // get sum of 4-momenta of particles in self-decay vertices
      for(unsigned int i=0;i<daughters.size();i++)
      {
        // ignore photons
        if(daughters[i]->getPdgID()==22) continue;

        // since 'allDaughtersSelfDecay()' is true
        // each of these particles has exactly one daughter
        px2 += daughters[i]->getDaughters().at(0)->getPx();
        py2 += daughters[i]->getDaughters().at(0)->getPy();
        pz2 += daughters[i]->getDaughters().at(0)->getPz();
        e2  += daughters[i]->getDaughters().at(0)->getE();
      }

      //cout<<"ORIG: "<<px1<<" "<<py1<<" "<<pz1<<" "<<e1<<endl;
      //cout<<"SELF: "<<px2<<" "<<py2<<" "<<pz2<<" "<<e2<<endl;

      p1 = m_particle_list.at(0)->createNewParticle(0,-1,0.0,px1,py1,pz1,e1);
      p2 = m_particle_list.at(0)->createNewParticle(0,-2,0.0,px2,py2,pz2,e2);

      // Finaly, boost photons to appropriate frame
      for(unsigned int i=0;i<photon_list.size();i++)
      {
        PhotosParticle *boosted = photon_list[i]->createNewParticle( 22, 1,
                                    0.0,
                                    photon_list[i]->getPx(),
                                    photon_list[i]->getPy(),
                                    photon_list[i]->getPz(),
                                    photon_list[i]->getE()   );
                
        boosted->boostToRestFrame(p1);
        boosted->boostFromRestFrame(p2);
        
        photon_list[i]->createSelfDecayVertex(boosted);
        
        delete boosted;
      }

      Log::Warning()<<"Hidden interaction, all daughters self decay."
          <<"Potentially over simplified solution applied."<<endl;
    }
  }

  //otherwise loop over particles which are already in the
  //event record and modify their 4 momentum
  //4.03.2012: Fix to prevent kinematical trap in vertex of simultaneous:
  //           z-collinear and non-conservation pf E,p for dauthters of grandmothers
  for(index=daughters_start; index < particle_count && index < (int) m_particle_list.size(); index++){

    PhotosParticle * particle = m_particle_list.at(index);

    if(hep.idhep[index]!=particle->getPdgID())
      Log::Fatal("PH_HEPEVT_Interface::get(): Something is wrong with the PH_HEPEVT common block",5);

    // If photons were added - for each daughter create a history entry
    if(isPhotonCreated && Photos::isCreateHistoryEntries)
    {
      particle->createHistoryEntry();
    }
    
    //check to see if this particle's 4-momentum has been modified
    bool   update=false;

    // don't update particle if difference lower than THRESHOLD * particle energy (default threshold = 10e-8)
    double threshold = NO_BOOST_THRESHOLD*hep.phep[index][3];
    if( fabs(hep.phep[index][0]-particle->getPx()) > threshold ||
        fabs(hep.phep[index][1]-particle->getPy()) > threshold ||
        fabs(hep.phep[index][2]-particle->getPz()) > threshold ||
        fabs(hep.phep[index][3]-particle->getE())  > threshold    ) update=true;

    if(update)
    {

      //modify this particle's momentum and it's daughters momentum
      //Steps 1., 2. and 3. must be executed in order.

      //1. boost the particles daughters into it's (old) rest frame
      particle->boostDaughtersToRestFrame(particle);

      //2. change this particles 4 momentum
      particle->setPx(hep.phep[index][0]);
      particle->setPy(hep.phep[index][1]);
      particle->setPz(hep.phep[index][2]);
      particle->setE(hep.phep[index][3]);

      //3. boost the particles daughters back into the lab frame
      particle->boostDaughtersFromRestFrame(particle);

      if(special && particle->getDaughters().size()>0){

      // Algorithm for special case:
      // a. get self-daughter of 'particle'
      PhotosParticle *particled = particle->getDaughters().at(0);

      // b. boost 'particled' daughters to rest frame
      particled->boostDaughtersToRestFrame(particled);

      // c. copy four momentum of 'particle' into four momentum of
      //    its self-daughter 'particled'

      particled->setPx( particle->getPx() );
      particled->setPy( particle->getPy() );
      particled->setPz( particle->getPz() );
      particled->setE ( particle->getE()  );

      // d. boost self daughter to rest-frame of <e1>
      //    boost self daughter from rest-frame of <e2>

      particled->boostToRestFrame(p1);
      particled->boostFromRestFrame(p2);

      // e. boost the 'particled' daughters back into the lab frame
      particled->boostDaughtersFromRestFrame(particled);
      }

    }
  }
  
  // cleanup
  if(p1) delete p1;
  if(p2) delete p2;
}

void PH_HEPEVT_Interface::prepare()
{
	check_ME_channel();
}

void PH_HEPEVT_Interface::complete()
{

}

void PH_HEPEVT_Interface::check_ME_channel()
{
	ME_channel=0;

// Check mothers:

	if(decay_idx==2)                              return; // Only one mother present
	if(hep.idhep[0]*hep.idhep[1]>0) return; // Mothers have same sign

	Log::Debug(900)<<"ME_channel: Mothers PDG:  "<<hep.idhep[0]<<" "<<hep.idhep[1]<<endl;
	if(decay_idx)
		Log::Debug(900,false)<<"            Intermediate: "<<hep.idhep[decay_idx-1]<<endl;

	int              firstDaughter=3;
	if(decay_idx==0) firstDaughter=2; // if no intermediate particle - daughters start at idx 2

	// Are mothers in range +/- 1-6; +/- 11-16?
	int mother1 = abs(hep.idhep[0]);
	int mother2 = abs(hep.idhep[1]);
	if( mother1<1 || (mother1>6 && mother1<11) || mother1>16 ) return;
	if( mother2<1 || (mother2>6 && mother2<11) || mother2>16 ) return;

//Check daughters

	// Z: check for pairs 11 -11 ; 13 -13 ; 15 -15
	// -------------------------------------------
	int firstPDG =0;
	int secondPDG=0;
	for(int i=firstDaughter; i<hep.nhep;i++)
	{
		int pdg = abs(hep.idhep[i]);
		if(pdg==11 || pdg==13 || pdg==15)
		{
			if(firstPDG==0) firstPDG=hep.idhep[i];
			else
			{
				secondPDG=hep.idhep[i];
				// Just in case two pairs are genereted - verify that we have a pair with oposite signs
				if(firstPDG*secondPDG>0) secondPDG=0;
				break;
			}
		}
	}

	if( ME_channel==0 && firstPDG!=0 && secondPDG!=0 &&
	    firstPDG==-secondPDG ) ME_channel=1;

	// W: check for pairs 11 -12; -11 12; 13 -14; -13 14; 15 -16; -15 16
	// -----------------------------------------------------------------
	firstPDG =0;
	secondPDG=0;
	for(int i=firstDaughter; i<hep.nhep;i++)
	{
		int pdg = abs(hep.idhep[i]);
		if(pdg>=11 && pdg<=16)
		{
			if(firstPDG==0) firstPDG=hep.idhep[i];
			else
			{
				secondPDG=hep.idhep[i];
				// Just in case two pairs are genereted - verify that we have a pair with oposite signs
				if(firstPDG*secondPDG>0) secondPDG=0;
				break;
			}
		}
	}

	firstPDG =abs(firstPDG);
	secondPDG=abs(secondPDG);

	if(  ME_channel==0 && firstPDG!=0 && secondPDG!=0 &&
	   ( ( firstPDG==11 && secondPDG==12 ) || (firstPDG == 12 && secondPDG == 11) ||
	     ( firstPDG==13 && secondPDG==14 ) || (firstPDG == 14 && secondPDG == 13) ||
	     ( firstPDG==15 && secondPDG==16 ) || (firstPDG == 16 && secondPDG == 15)
	   )
	  ) ME_channel=2;

	Log::Debug(901)<<"ME_channel: Found ME_channel: "<<ME_channel<<endl;

// Check intermediate particle (if exists):

	// Verify that intermediate particle PDG matches ME_channel found
	if(ME_channel>0 && decay_idx)
	{
		int pdg=hep.idhep[decay_idx-1];

		if(ME_channel==1 && !(pdg==22 || pdg==23) ) ME_channel=0; //gamma/Z
		if(ME_channel==2 && !(pdg==24 || pdg==-24)) ME_channel=0; //W+/W-

		if(ME_channel==0)
			Log::Debug(901,false)<<"            but set to 0: wrong intermediate particle: "<<pdg<<endl;
	}

// Check flags

	switch(ME_channel)
	{
		case  0: break; // Ok - no channel found
		case  1: if(!Photos::meCorrectionWtForZ) ME_channel=0; break;
		case  2: if(!Photos::meCorrectionWtForW) ME_channel=0; break;
		default: Log::Error()<<"Unimplemented ME channel: "<<ME_channel<<endl; break;
	}
	Log::Debug(902)<<"ME_channel: Finally, after flag check, ME_channel is: "<<ME_channel<<endl;
}


} // namespace Photospp
