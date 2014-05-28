#ifndef _PhotosBranch_h_included_
#define _PhotosBranch_h_included_

/**
 * @class PhotosBranch
 *
 * @brief Single branching point
 *
 * Contains information about daughters and mothers of a single branch.
 * Each branch will be converted to HEPEVT and processed by photos.
 *
 * @author Tomasz Przedzinski
 * @date 8 July 2010
 */

#include <vector>
#include "PhotosParticle.h"
using std::vector;

namespace Photospp
{

class PhotosBranch
{
public:
	/** Create branch out of decaying particle */
	PhotosBranch(PhotosParticle* p);

	/** Return decaying particle. NULL if branching does not have mid-particle */
	PhotosParticle*          getDecayingParticle() { return particle;  }

	/** Get list of mothers */
	vector<PhotosParticle *> getMothers()          { return mothers;   }

	/** Get list of daughters */
	vector<PhotosParticle *> getDaughters()        { return daughters; }

	/** Get list of all particles used by branch */
	vector<PhotosParticle *> getParticles();

	/** Check if branch is suppressed */
	int getSuppressionStatus() { return suppression; }

	/** Check if branch is forced */
	int getForcingStatus()     { return forcing; }

	/** Checks momentum conservation of decaying particle.
	    If it does not exist, checks momentum of first mother passed to photos */
	bool checkMomentumConservation();

	/** Process single branch */
	void process();

	/** Create branches from particles list */
	static vector<PhotosBranch *> createBranches(vector<PhotosParticle *> particles);
private:
	/** Checks if branching is suppressed by PHOTOS. */
	int checkSuppressionLevel() { return checkList(false); }

	/** Checks if branching is forced by PHOTOS. */
	int checkForcingLevel()     { return checkList(true);  }

	/** Algorithm used for suppression/forcing check */
	int checkList(bool forceOrSuppress);
private:
	/** State of branching suppression*/
	int suppression;
	/** State of branching forcing*/
	int forcing;
	/** Decaying particle */
	PhotosParticle          *particle;
	/** List of mothers   */
	vector<PhotosParticle *> mothers;
	/** List of daughters */
	vector<PhotosParticle *> daughters;
};

} // namespace Photospp
#endif
