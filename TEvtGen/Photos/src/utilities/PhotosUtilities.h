#ifndef _PhotosUtilities_h_included_
#define _PhotosUtilities_h_included_
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;
/**
 * @class PhotosUtilities
 *
 * @brief Support functions
 *
 * Functions for boosting, rotation, ...
 *
 * @author Tomasz Przedzinski, Zbigniew Was
 * @date 29 June 2013
 */
namespace Photospp
{

namespace PhotosUtilities
{
  /** PHOton radiation in decays calculation of TRIangle fie */
  double PHOTRI(double A,double B,double C);

  /** PHOton radiation in decays Calculate ANgle from X and Y */
  double PHOAN1(double X,double Y);

  /** PHOton radiation in decays Calculate ANgle from X and Y equiv to PHOAN1 */
  double PHOAN2(double X,double Y);

  /** PHOton radiation in decays ROtation routine around 2-nd axis */
  void PHORO2(double ANGLE,double PVEC[4]);


  /** PHOton radiation in decays ROtation routine around 3-rd axis */
  void PHORO3(double ANGLE,double PVEC[4]);

  /** Boot to-from restr frame of PBOOS1 */
  void PHOB(int MODE,double PBOOS1[4],double vec[4]);

  /** PHOton radiation in decays BOost routine along arbitrary axis axis */
  void bostdq(int mode,double qq[4],double pp[4],double r[4]);


  /** PHOton radiation in decays BOost routine along 3-rd axis */
  void PHOBO3(double ANGLE,double PVEC[4]);

  /** trivial method to fill value into array on positions beg to end */
  void fill_val(int beg, int end, double* array, double value); 

 /** PHOeps vector product (normalized to unity) */
  void PHOEPS(double vec1[4], double vec2[4], double eps[4]);

 /** PHOton radiation  in decays function for SPIn determination */
  double PHOSPI(int idhep);

 /**  PHOton radiation in decays CHArge determination */
  double PHOCHA(int idhep);

} // namespace PhotosUtilities

} // namespace Photospp
#endif

