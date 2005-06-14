////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRIGGER_RECORD_HPP
#define dHLT_TRIGGER_RECORD_HPP

#include "BasicTypes.hpp"
#include "Point.hpp"

namespace dHLT
{


/* The sign of the particle as given by L0.
 */
enum ParticleSign
{
    Minus       = -1,
    UnknownSign = 0,
    Plus        = 1
};


typedef UInt TriggerRecordID;


/* Data structure containing information about L0 validated trigger hit.
 */
struct TriggerRecord
{
	ParticleSign sign;     // The sign of the particle.
	Float pt;              // Transverse momentum of the particle.
	Point station1impact;  // Impact point of particle on trigger station 1.
	Point station2impact;  // Impact point of particle on trigger station 2.


	/* Default constructor.
	   Sets the fSign to UnknownSign, fPt to -1 and the impact points are set to zero.
	 */
	TriggerRecord()
		: station1impact(), station2impact()
	{
		sign = UnknownSign;
		pt = -1.0;
	};

	/* Creates a trigger record with the specifed particle sign, pt and impact points.
	   The impactpoint1 corresponds to trigger station 1 and simmilarly impactpoint2
	   corresponds to station 2.
	 */
	TriggerRecord(
			const ParticleSign particle_sign, const Float particle_pt,
			const Point impactpoint1, const Point impactpoint2
		)
	{
		sign = particle_sign;
		pt = particle_pt;
		station1impact = impactpoint1;
		station2impact = impactpoint2;
	};
};


} // dHLT

#endif // dHLT_TRIGGER_RECORD_HPP
