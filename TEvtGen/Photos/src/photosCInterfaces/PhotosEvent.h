#ifndef _PhotosEvent_h_included_
#define _PhotosEvent_h_included_

/**
 * @class PhotosEvent
 *
 * @brief Abstract base class for containing the event information.
 *
 * PhotosEvent contains virtual methods, which need to be implemented
 * by the appropriate interface class to the event record. An object of
 * PhotosEvent type should be created by the user and processed
 * via the process() method.
 *
 * @author Nadia Davidson
 * @date 16 June 2008
 */

#include <vector>
#include "PhotosBranch.h"
#include "PhotosParticle.h"
using std::vector;

namespace Photospp
{

class PhotosEvent
{
public:
	virtual ~PhotosEvent();

	/** Get an unfiltered list of particles from the event record */
	virtual vector<PhotosParticle*> getParticleList() = 0;

	/** Print informations about the event */
	virtual void print() = 0;

	/** Process event */
	void process();
private:
	/** Filter suppressed and invalid particles. */
	vector<PhotosParticle *> filterParticles(vector<PhotosParticle *> particles);

	/** branch points which should be given to PHOTOS */
	vector<PhotosBranch *> m_branch_points;
};

} // namespace Photospp
#endif
