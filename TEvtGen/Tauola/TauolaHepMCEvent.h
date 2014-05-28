#ifndef _TauolaHepMCEvent_h_included_
#define _TauolaHepMCEvent_h_included_

/**
 * @class TauolaHepMCEvent
 *
 * @brief Interface to HepMC::GenEvent objects
 *
 * This class implements the virtual methods of
 * TauolaEvent. In this way it provides an
 * interface between the generic TauolaEvent class
 * and a HepMC::GenEvent object.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 *
 * This code is licensed under GNU General Public Licence.
 * For more informations, see: http://www.gnu.org/licenses/
 */

#include <iostream>
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "TauolaEvent.h"
#include "TauolaParticle.h"
#include "TauolaHepMCParticle.h"

namespace Tauolapp
{

class TauolaHepMCEvent : public TauolaEvent{

 public:

  /** Constructor which keeps a pointer to the HepMC::GenEvent*/
  TauolaHepMCEvent(HepMC::GenEvent * event);

  ~TauolaHepMCEvent();

  /** Returns the HepMC::GenEvent */
  HepMC::GenEvent * getEvent();

  /** Implementation of TauolaEvent virtual method.
      This returns a list of particles in the event with 
      pdg id = "pdgID". */
  std::vector<TauolaParticle*> findParticles(int pdgID);

  /** Implementation of TauolaEven virtual method.
      This returns a list of particles in the event with
      pdg id = "pdgID" and stable status code. */
  std::vector<TauolaParticle*> findStableParticles(int pdgID);

  /** Overriding of TauolaEvent decayEndgame method.
      Converts the momentum and length units */
  void eventEndgame();

 private:

  /** The event */
  HepMC::GenEvent * m_event;
  /** List of particles to be decayed */
  std::vector<TauolaParticle*> m_tau_list;
  /** Momentum unit name */
  string m_momentum_unit;
  /** Length unit name */
  string m_length_unit;

};

} // namespace Tauolapp
#endif  

