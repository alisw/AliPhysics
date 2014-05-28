#ifndef _PhotosDebugRandom_included_
#define _PhotosDebugRandom_included_
/**
 * @class PhotosDebugRandom
 *
 * @brief Random generator for debugging purposes
 *
 * Random generator with ability to save and restore state, set custom
 * state and print current state. This class extends PhotosRandom class.
 * PhotosRandom is a static class, therefore this extension works
 * automatically for FORTRAN code without any changes.
 * Our random generator is the C++ follow-up of FORTRAN generator by 
 * F. James DD-Report, November 1988. and  G. Marsaglia and A. Zaman,  
 * FSU-SCR-87-50,
 *
 * @author Tomasz Przedzinski
 * @date 28th June 2012
 */
 
#include "PhotosRandom.h"

namespace Photospp
{

class PhotosDebugRandom : PhotosRandom
{
public:
  /* Save current state */
  static void saveState();
  
  /* Restore state from save */
  static void restoreState();

  /* Set current state provided by user */
  static void setState(int i, int j, double c, double list[97]);

  /* Save state provided by user */
  static void setSaveState(int i, int j, double c, double list[97]);
  
  /* Print state in a form that can be easily copied into the code */
  static void print();

private:
  static int i97_saved;
  static int j97_saved;
  static double uran_saved[97];
  static double cran_saved;
};

} // namespace Photospp
#endif
