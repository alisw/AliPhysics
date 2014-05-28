#ifndef _TauolaEvent_h_included_
#define _TauolaEvent_h_included_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TauolaParticlePair.h"

/**
 * @class TauolaEvent
 *
 * @brief Abstract base class for containing the event information.
 *
 * TauolaEvent contains virtual methods, which need to be implemented
 * by the appropriate interface class to the event record. Currently only
 * TauolaHepMCEvent does this. An object of TauolaEvent type should be
 * created by the user and can be decayed via the decayTaus() method.
 *
 * This class is responsible for finding taus, (or tau and 
 * it's neutrino) and creating TauolaParticlePairs out of them.
 *
 * @author Nadia Davidson
 * @date 16 June 2008
 */

namespace Tauolapp
{

class TauolaEvent{

 public:
  virtual ~TauolaEvent(){};

   /** create TauolaParticlePairs */
   std::vector<TauolaParticle*> findPairs();

  /** Decay taus in this event.*/
   void decayTaus();

   /** Undecay taus in this event but removing their daughters and
    returning the status cods to 1.*/
   void undecayTaus();

   /** Final touches to event record after all decays are finished.
       Some event records (e.g. HepMC) need it. */
   virtual void eventEndgame() {}

   /** return a list of all particle with pdg_id = absolute value of pdg_id.
       This method must be implemented by a derived class. eg. 
       TauolaHepMCEvent */
   virtual std::vector<TauolaParticle*> findParticles(int pdg_id)=0;

   /** return a list of all particle with pdg_id = absolute value of pdg_id
       and stable status code. This method must be implemented by a derived class. 
       eg. TauolaHepMCEvent */
   virtual std::vector<TauolaParticle*> findStableParticles(int pdg_id)=0;


 private:    

};

} // namespace Tauolapp
#endif  

