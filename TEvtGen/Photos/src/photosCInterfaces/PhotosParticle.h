#ifndef _PhotosParticle_h_included_
#define _PhotosParticle_h_included_

/**
 * @class PhotosParticle
 *
 * @brief Abstract base class for particle in the event. This class also
 * handles boosting.
 *
 * PhotosParticle is a Photos representation of a particle. It has virtual
 * getter and setter methods that need to be implemented by a derived class.
 * An example of this is PhotosHepMCParticle. In this way it provides an
 * interface to the information in the Event Record.
 * 
 * @author Nadia Davidson
 * @date 16 June 2008
 */

#include <vector>
#include "Photos.h"

namespace Photospp
{

class PhotosParticle
{
public:
	/** Stable particle status */
	static const int STABLE=1;

	/** Decayed particle status */
	static const int DECAYED=2;

	/** History particle status */
	static const int HISTORY=3;

	/** X Axis */
	static const int X_AXIS=1;

	/** Y Axis */
	static const int Y_AXIS=2;

	/** Z Axis */
	static const int Z_AXIS=3;

	/** Z0 particle */
	static const int Z0 = 23;

	/** H particle */
	static const int HIGGS = 25;

	/** H0 particle */
	static const int HIGGS_H = 35;

	/** A0 particle */
	static const int HIGGS_A = 36;

	/** H+ particle */
	static const int HIGGS_PLUS = 37;

	/** H- particle */
	static const int HIGGS_MINUS = -37;

	/** W+ particle */
	static const int W_PLUS = 24;

	/** W- particle */
	static const int W_MINUS = -24;

	/** photon */
	static const int GAMMA = 22;

	/** tau+ particle */
	static const int TAU_PLUS = -15;

	/** tau- particle */
	static const int TAU_MINUS = 15;

	/** tau neutrino particle */
	static const int TAU_NEUTRINO = 16;

	/** tau antineutrino particle */
	static const int TAU_ANTINEUTRINO = -16;

	/** muon+ particle */
	static const int MUON_PLUS = -13;

	/** muon- particle */
	static const int MUON_MINUS = 13;

	/** muon neutrino particle */
	static const int MUON_NEUTRINO = 14;

	/** muon antineutrino particle */
	static const int MUON_ANTINEUTRINO = -14;

	/** e+ particle */
	static const int POSITRON = -11;

	/** e- particle */
	static const int ELECTRON = 11;

	/** e neutrino particle */
	static const int ELECTRON_NEUTRINO = 12;

	/** e antineutrino particle */
	static const int ELECTRON_ANTINEUTRINO = -12;

	/** up quark */
	static const int UP = 2;

	/** anti-up quark */
	static const int ANTIUP = -2;

	/** down quark */
	static const int DOWN = 1;

	/** anti-down quark */
	static const int ANTIDOWN = -1;

	/** All other particle types*/
	static const int OTHER = 0;

public:
	virtual ~PhotosParticle(){};

	/** Return whether the particle has any chidren */
	bool hasDaughters();

	/** Traverse the event structure and find the final version
	    of this particle which does not have a particle of it's own type
	    as it's daughter. eg. Generally the final stable copy */
	PhotosParticle * findLastSelf();

	/** Traverse the event structure and find the first set of mothers
	    which are not of the same type as this particle. */
	std::vector<PhotosParticle *> findProductionMothers();

	/** Return whole decay tree starting from this particle */
	std::vector<PhotosParticle *> getDecayTree();

	/** Transform this particles four momentum from the lab frome
	    into the rest frame of the paramter PhotosParticle. */
	void boostToRestFrame(PhotosParticle * boost);

	/** Transform the four momentum of all the daughters recursively
	    into the frame of the "particle" PhotosParticle. */
	void boostDaughtersToRestFrame(PhotosParticle * boost);

	/** Transform this particles four momentum from the rest frame of
	    the paramter PhotosParticle, back into the lab frame. */
	void boostFromRestFrame(PhotosParticle * boost);

	/** Transform this particles four momentum from the lab frame to
	    the rest frame of the parameter PhotosParticle. */
	void boostDaughtersFromRestFrame(PhotosParticle * boost);

	/** Do a Lorenz transformation along the Z axis. */
	void boostAlongZ(double pz, double e);

	/** rotate this particles 4-momentum by an angle phi from
	    the axisis "axis" towards the axis "second_axis". */
	void rotate(int axis, double phi, int second_axis=Z_AXIS);

	/** rotate 4-momentum of daughters of this particle by an angle phi from
	    the axisis "axis" towards the axis "second_axis". */
	void rotateDaughters(int axis, double phi, int second_axis=Z_AXIS);

	/** Returns the angle around the axis "axis" needed to rotate
	    the four momenum is such a way that the non-Z component
	    disappears and Z>0. This is used to in rotating the coordinate
	    system into a frame with only a Z component before calling
	    boostAlongZ(). */
	double getRotationAngle(int axis, int second_axis=Z_AXIS);

	/** Get scalar momentum */
	double getP();

	/** Get momentum component in the direction of "axis" (x,y,z) */
	double getP(int axis);

	/** Set momentum component in the direction of "axis" (x,y,z) */
	void  setP(int axis, double p_component);

	/** Get sqrt(e^2-p^2) */
	virtual double getVirtuality();

public:
	/** check that the 4 momentum in conserved at the vertices producing
	    and ending this particle */
	virtual bool checkMomentumConservation()=0;

	/** Returns the px component of the four vector */
	virtual double getPx()=0;

	/** Returns the py component of the four vector */
	virtual double getPy()=0;

	/** Returns the pz component of the four vector */
	virtual double getPz()=0;

	/** Returns the energy component of the four vector */
	virtual double getE()=0;
	
	/** Get the invariant mass from the event record*/
	virtual double getMass() = 0;

	/** Set the px component of the four vector */
	virtual void setPx( double px )=0;

	/** Set the px component of the four vector */
	virtual void setPy( double py )=0;

	/** Set the pz component of the four vector */
	virtual void setPz( double pz )=0;

	/** Set the energy component of the four vector */
	virtual void setE( double e )=0;

	/** Set the mothers of this particle via a vector of PhotosParticle */
	virtual void setMothers(std::vector<PhotosParticle*> mothers)=0;

	/** Set the daughters of this particle via a vector of PhotosParticle */
	virtual void setDaughters(std::vector<PhotosParticle*> daughters)=0;

	/** Add a new daughter to this particle */
	virtual void addDaughter(PhotosParticle* daughter)=0;

	/** Returns the mothers of this particle via a vector of PhotosParticle */
	virtual std::vector<PhotosParticle*> getMothers()=0;

	/** Returns the daughters of this particle via a vector of PhotosParticle */
	virtual std::vector<PhotosParticle*> getDaughters()=0;

	/** Returns all particles in the decay tree of this particle
	    via a vector of PhotosParticle */
	virtual std::vector<PhotosParticle*> getAllDecayProducts()=0;

	/** Set the PDG ID code of this particle */
	virtual void setPdgID(int pdg_id)=0;

	/** Set the mass of this particle */
	virtual void setMass(double mass)=0;

	/** Set the status of this particle */
	virtual void setStatus(int status)=0;

	/** Get the PDG ID code of this particle */
	virtual int getPdgID()=0;

	/** Get the status of this particle */
	virtual int getStatus()=0;

	/** Get the barcode of this particle */
	virtual int getBarcode()=0;

	/** Create a new particle of the same type, with the given
	    properties. The new particle bares no relations to this
	    particle, but it provides a way of creating a intance of
	    the derived class. eg. createNewParticle() is used inside
	    filhep_() so that an eg. PhotosHepMCParticle is created without
	    the method having explicit knowledge of the PhotosHepMCParticle
	    class */
	virtual PhotosParticle * createNewParticle(int pdg_id, int status,
	                                           double mass, double px,
	                                           double py, double pz,
	                                           double e)=0;

  /** Create history entry of this particle before modifications
      of PHOTOS. Implementation of this method depends strongly
      on the event record. */
  virtual void createHistoryEntry()=0;
  
	/** Create a self-decay vertex for this particle
	    with 'out' being the outgoing particle in new vertex */
	virtual void createSelfDecayVertex(PhotosParticle *out)=0;

	/** Print some information about this particle to standard output */
	virtual void print()=0;
};

} // namespace Photospp
#endif
