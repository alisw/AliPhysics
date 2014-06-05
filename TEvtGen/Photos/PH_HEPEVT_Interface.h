#ifndef _PH_HEPEVT_Interface_included_
#define _PH_HEPEVT_Interface_included_

#include <vector>
#include "PhotosBranch.h"
#include "PhotosParticle.h"
#include "f_Init.h"

namespace Photospp
{
  // const static int NMXHEP = 10000; at present NMXHEP is defined in f_Init.h
const static double NO_BOOST_THRESHOLD=1.0e-8;

class PH_HEPEVT_Interface
{
public:
	/** Convert PhotosBranch to HEPEVT */
	static int  set(PhotosBranch* branch);

	/** Update event record with data from HEPEVT */
	static void get();

	/** Prepare particles for processing */
	static void prepare();

	/** Check channel for complete matrix element calculation */
	static void check_ME_channel();

	/** Finalize processing */
	static void complete();

	/** Clear HEPEVT */
	static void clear();
public:
	/** Index of decaying particle*/
	static int decay_idx;
	/** Number of channel to be used - flag for fortran routines */
	static int ME_channel;
private:
	/** Add single particle to HEPEVT */
	static void add_particle(int i, PhotosParticle * particle,
	                         int first_mother, int last_mother,
	                         int first_daughter, int last_daughter);

	/** List of particles added to HEPEVT */
	static std::vector<PhotosParticle*> m_particle_list;
};

} // namespace Photospp
#endif
