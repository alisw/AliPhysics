#ifndef __PLOTS_CLASS_HEADER__
#define __PLOTS_CLASS_HEADER__

/**
 * This class contains routines for plot generation.
 *
 * The output files are adopted for the draw.C root script. 
 *
 * @author Tomasz Przedzinski
 * @date 20 December 2009
 */

#include "TauolaParticlePair.h"

namespace Tauolapp
{

class Plots
{
public:
  /**   SANC tables plots
  Writes the data for plots for the 1-1 table or 2-2 table
  using 11-11 table as the born-level table, so before
  running, the 11-11 table must be substituted with born-level
  table for either 1-1 or 2-2. */
  Plots();

  /** Sets cosTheta (for plots 1 and 2)
            and incoming particle pdgid (for all plots). */
  void setSancVariables(int inc, double cos);

  /** SANC test - three functions - table, born level and plzap0 for selected cosTheta */
  void SANCtest1();

  /** Weights - three functions - w, w0 and w/w0 for selected cosTheta */
  void SANCtest2();

  /** Error check - one function - table vs born for all cosTheta */
  void SANCtest3();

  /** cross-section - three functions - w, w0 and w/w0 for all cosTheta*/
  void SANCtest4();
private:
  /* Incoming particle PDG ID */
  int    m_incoming_pdg_id;
  /* cos(theta) used for plots */
  double m_cosTheta;
  /* Number of points in plot */
  int    m_n_plot_points;
  /* TauolaParticlePair class object */
  TauolaParticlePair t_pair;
};

} // namespace Tauolapp
#endif
