//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtParticle.hh
//
// Description:Class to describe all particles
//
// Modification history:
//
//    DJL/RYD     Sept. 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTPARTICLE_HH
#define EVTPARTICLE_HH

//#include <iostream.h>
#include <assert.h>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include <string>
#include <vector>

class EvtDiracSpinor;
class EvtVector4C;
class EvtTensor4C;
class EvtStdHep;
class EvtSecondary;
class EvtRaritaSchwinger;

const int MAX_DAUG =100;
const int MAX_LEVEL=10;
const int MAX_TRIES=10000;

class EvtParticle{

public:

  /**
  * Default constructor.
  */  
  EvtParticle();

  /**
  * Destructor.
  */
  virtual ~EvtParticle();

  /**
  * Returns polarization vector in the parents restframe.
  */
  virtual EvtVector4C epsParent(int i) const; 

  /**
  * Returns polarization vector in the particles own restframe.
  */
  virtual EvtVector4C eps(int i) const; 


  /**
  * Returns polarization vector in the parents restframe for a photon.
  */
  virtual EvtVector4C epsParentPhoton(int i); 

  /**
  * Returns polarization vector in the particles own restframe for a photon.
  */
  virtual EvtVector4C epsPhoton(int i); 


  /**
  * Returns Dirac spinor in the parents restframe for a Dirac particle.
  */
  virtual EvtDiracSpinor spParent(int) const;

  /**
  * Returns Dirac spinor in the particles own restframe for a Dirac particle.
  */
  virtual EvtDiracSpinor sp(int) const;



  /**
  * Returns Dirac spinor in the parents restframe for a Neutrino particle.
  */
  virtual EvtDiracSpinor spParentNeutrino() const;

  /**
  * Returns Dirac spinor in the particles own restframe for a
  * Neutrino particle.
  */
  virtual EvtDiracSpinor spNeutrino() const;
  

  /**
  * Returns tensor in the parents restframe for a spin 2 particle.
  */
  virtual EvtTensor4C epsTensorParent(int i) const; 

  /**
  * Returns tensor in the particles own restframe for a spin 2 particle.
  */
  virtual EvtTensor4C epsTensor(int i) const; 
  
  /**
   * Returns Rarita-Schwinger spinor in the parents restframe for a 
   * Rarita-Schwinger particle.
   */
  virtual EvtRaritaSchwinger spRSParent(int) const;
  
  /**
   * Returns Rarita-Schwinger spinor in the particles own restframe for a 
   * Rarita-Schwinger particle.
   */
  virtual EvtRaritaSchwinger spRS(int) const;
  


  /**
  * Initialiaze particle with id and 4momentum.
  */ 
  virtual void init(EvtId part_n,const EvtVector4R& p4)=0;
  
  /**
  * Add another daughter to the particle
  */
  void addDaug(EvtParticle *node); 
  
  /**
  * Decay particle
  */
  void decay();

  /** 
  * Delete a decay chain
  */
  void deleteTree();
  void deleteDaughters(bool keepChannel=false);

  /**
  * Should only be used internally.
  */ 
  void setChannel( int i );
  
  /**
  * Creates the daughters in the list of ids and 
  * adds them to the parent. Note that momentum
  * is left uninitialized, this is _only_ creation.
  */ 
  void makeDaughters(unsigned int ndaug,EvtId *id);

  /**
  * Creates the daughters in the list of ids and 
  * adds them to the parent. Note that momentum
  * is left uninitialized, this is _only_ creation.
  */ 
  void makeDaughters(unsigned int ndaug, std::vector<EvtId> idVector);

  /**
  * Similar to the routine above except that here 
  * momentum is generated according to phase space 
  * daughters are filled with this momentum.
  */ 
  double initializePhaseSpace(unsigned int numdaughter,EvtId *daughters,
			      bool forceResetMasses=false,
			      double poleSize=-1., int whichTwo1=0, 
			      int whichTwo2=1); 
  
  /**
  * Get pointer the the i:th daugther.
  */  
  EvtParticle *getDaug(int i);

  /**
  * Iterates over the particles in a decay chain.
  */
  EvtParticle *nextIter(EvtParticle *rootOfTree=0);

  /**
  * Makes stdhep list
  */
  void makeStdHep(EvtStdHep& stdhep,
		 EvtSecondary& secondary,EvtId *stable_parent_ihep);
  void makeStdHep(EvtStdHep& stdhep);

  /**
  * Gets 4vector in the labframe, i.e., the frame in which the root
  * particles momentum is measured.
  */
  EvtVector4R getP4Lab() const;


  /**
  * Gets 4vector in the labframe for the 4-momentum befor FSR was 
  * generated in the parents decay. The lab frame is where the root
  * particles momentum is measured.
  */
  EvtVector4R getP4LabBeforeFSR();
 
  /**
  * Gets 4vector in the particles restframe, i.e. this functiont will
  * return (m,0,0,0)
  */
  EvtVector4R getP4Restframe() const;

  /**
  * Returns the 4position of the particle in the lab frame.
  */
  EvtVector4R get4Pos() const;

  /**
  * Returns pointer to parent particle.
  */
  EvtParticle *getParent() const;
  
  /**
  * Makes partptr the idaug:th daugther.
  */ 
  void insertDaugPtr(int idaug,EvtParticle* partptr){ _daug[idaug]=partptr;
                                                 partptr->_parent=this; }
  /**
  * Returns mass of particle.
  */
  double mass() const;

  /**
  * Used internally to decide if first time particle is decayed.
  */
  int firstornot() const;
  void setFirstOrNot();
  void resetFirstOrNot();

  /**
  * Returns Id of particle.
  */
  EvtId getId() const;

  /**
  * Returns the PDG id of the particle
  */
  
  int getPDGId() const;

  /** 
  * Returns particle type.
  */

  EvtSpinType::spintype getSpinType() const;

  /**
  * Returns number of spin states of the particle.
  */
  int getSpinStates() const;

  /**
  * Returns 4momentum in parents restframe.
  */
  const EvtVector4R& getP4() const;

  /**
  * Sets the 4momentum in the parents restframe.
  */
  void setP4(const EvtVector4R& p4){
    _p=p4;
    _pBeforeFSR=p4;
  }

  void setP4WithFSR(const EvtVector4R& p4){
    _p=p4;
  }

  void setFSRP4toZero(){
    _pBeforeFSR.set(0.0,0.0,0.0,0.0);
  }

  /**
  * Retunrs the decay channel.
  */
  int getChannel() const;
  
  /**
  * Returns number of daugthers.
  */ 
  size_t getNDaug() const;
  void resetNDaug() {_ndaug=0; return;}

  /**
  * Prints out the particle "tree" of a given particle.  The
  * tree consists of all daughters (and their daughters, etc)
  * and their properties.
  */
  void printTree() const;

  void printTreeRec(unsigned int level) const;

  std::string treeStr() const;
  std::string treeStrRec(unsigned int level) const;

  /**
  * Prints information for the particle.
  */
  void printParticle() const;

  /**
  * Set lifetime of the particle in parents restframe.
  */  
  void setLifetime(double tau);
  
  /**
  * Generate lifetime according to pure exponential.
  */
  void setLifetime();

  /**
  * Returns the lifetime.
  */
  double getLifetime();

  /** 
  * Set diagonal spindensity matrix.
  */
  void setDiagonalSpinDensity();
  
  /** 
  * Set spindensity matrix for e+e- -> V
  */
  void setVectorSpinDensity(); 
  
  /**
  * Set forward spin density matrix.
  */
  void setSpinDensityForward(const EvtSpinDensity& rho){_rhoForward=rho;}

  /**
  * Set forward spin density matrix according to the density matrix
  * rho in the helicity amplitude basis.
  */
  void setSpinDensityForwardHelicityBasis(const EvtSpinDensity& rho);
  void setSpinDensityForwardHelicityBasis(const EvtSpinDensity& rho,
					  double alpha,
					  double beta,
					  double gamma);
  
  /**
  * Returns a rotation matrix need to rotate the basis state
  * to the helicity basis. The EvtSpinDensity matrix is just use
  * as a matrix here. This function is to be implemented in each
  * derived class.
  */
  virtual EvtSpinDensity rotateToHelicityBasis() const=0;
  virtual EvtSpinDensity rotateToHelicityBasis(double alpha,
					       double beta,
					       double gamma) const=0;

  /**
  * Get forward spin density matrix.
  */
  EvtSpinDensity getSpinDensityForward(){return _rhoForward;}

  /**
  * Set backward spin density matrix.
  */
  void setSpinDensityBackward(const EvtSpinDensity& rho){_rhoBackward=rho;}

  /**
  * Get backward spin density matrix.
  */
  EvtSpinDensity getSpinDensityBackward(){return _rhoBackward;}

  //Hacks will be removed when better solutions are thought of!
  //This is used to suppress use of random numbers when doing initialization
  //of some models.
  void noLifeTime() { _genlifetime=0; }

  //lange - April 29, 2002
  void setId(EvtId id) { _id=id;}
  void initDecay(bool useMinMass=false);
  bool generateMassTree();

  double compMassProb();

  //setMass will blow away any existing 4vector
  void setMass(double m) { _p=EvtVector4R(m,0.0,0.0,0.0);}

  //void setMixed() {_mix=true;}
  //void setUnMixed() {_mix=false;}
  //bool getMixed() {return _mix;}

  //void takeCConj() {report(Severity::Info,"EvtGen") << "should take conj\n";}

  //this means that the particle has gone through initDecay
  // and thus has a mass
  bool isInitialized() {return _isInit;}
  bool hasValidP4() {return _validP4;}
  bool isDecayed() {return _isDecayed;}


  // decay prob - only relevent if already decayed
  // and is a scalar particle
  // returned is a double* that should be prob/probMax
  double* decayProb() {return _decayProb;}
  void setDecayProb( double p);

  // Return the name of the particle (from the EvtId number)
  std::string getName();
   
protected:

  void setp( double e, double px, double py, double pz) { 
    _p.set(e,px,py,pz); 
    _pBeforeFSR=_p;
  }

  void setp( const EvtVector4R& p4 ) { 
    _p =p4; 
    _pBeforeFSR=_p;
  }

  void setpart_num(EvtId particle_number ) { 
    assert(_channel==-10||
	   _id.getId()==particle_number.getId()||
	   _id.getId()==-1);
    _id = particle_number; 
  }
  bool           _validP4;

private:

  EvtParticle*   _daug[MAX_DAUG];
  size_t         _ndaug; 
  EvtParticle*   _parent;
  int            _channel; 
  int            _first;
  EvtId          _id;
  EvtVector4R    _p;
  EvtVector4R    _pBeforeFSR;
  double         _t;
  bool            _isInit;
  bool           _isDecayed;

  //bool _mix;

  EvtSpinDensity _rhoForward;
  EvtSpinDensity _rhoBackward;

  void makeStdHepRec(int firstparent,int lastparent,EvtStdHep& stdhep,
		 EvtSecondary& secondary,EvtId *stable_parent_ihep);
  void makeStdHepRec(int firstparent,int lastparent,EvtStdHep& stdhep);


  //This is a hack until things gets straightened out. (Ryd)
  int         _genlifetime;
  
  //should never be used, therefor is private.
  //these does _not_ have an implementation
  EvtParticle& operator=(const EvtParticle& p);
  EvtParticle(const EvtParticle& p);

  double *_decayProb;

};

#endif

