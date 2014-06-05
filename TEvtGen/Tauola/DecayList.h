#ifndef _DecayList_h_included_
#define _DecayList_h_included_

/** 
 * @class DecayList
 *
 * @brief Storage class for keeping track of TauolaParticles and
 * their index (as passed to Tauola).
 *
 * This class contains a list of TauolaParticles. The index of the
 * TauolaParticle in the list is passed to FORTRAN TAUOLA (as we can 
 * not pass TauolaParticles directly). A static copy of the class is 
 * used for any run of Tauola and is cleared after each decay generation.
 * To be compatible with how indicies are used in Tauola, those that are:
 * - < 0 are relative to some other given position in the list
 * - = 0 is the position one greater than some given position in the list.
 * - > 0 are absolute. 1 is the first item (which should generally be 
 * the tau's position) 
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TauolaParticle.h"

using namespace std;

namespace Tauolapp
{

class TauolaParticle;

class DecayList {
  
 public:
  /** Return the TauolaParticle corresponding to the index (absolute)
      in the list of particle */
  static TauolaParticle * getParticle(int index);

  /** Adds the new particle into the list and delete the previous
     particle at the same position if it exists */
  static void updateList(TauolaParticle * new_particle,
                         int index);
  
  /** Adds the new particle to the end of list */ 
  static void addToEnd(TauolaParticle * new_particle);

  /** clear all entries from the list */
  static void clear();

  /** Translates index (absolute and relative) to
      absolute index. If a relative index is given (negative integer)
      it is taken relative from the end of the list */
  static int getAbsoluteIndex(int index);

  /** Translates index (absolute and relative) to
      absolute index. If a relative index is given (negative integer)
      it is taken relative to the parameter "neg_index_relative_to" */
  static int getAbsoluteIndex(int index, 
                              int neg_index_relative_to);

  /** Return index (absolute) of "particle" */
  static int getAbsoluteIndex(TauolaParticle * particle);

  /** Print the contents of the list */
  static void print();

 private:
  /** vector used for TauolaParticle mapping */
  static vector<TauolaParticle*> m_particle_list;

};

} // namespace Tauolapp
#endif
