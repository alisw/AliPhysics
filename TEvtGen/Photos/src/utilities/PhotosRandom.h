#ifndef _PhotosRandom_included_
#define _PhotosRandom_included_

/**
 * @class PhotosRandom
 *
 * @brief Photos random number generator rewritten from PHOTOS FORTRAN
 *
 * Generates uniformly distributed random numbers between 0 and 1.
 * Must be initialized by call to PhotosRandom::initialize().
 * Original authors: B. van Eijk, F. James, G. Marsaglia and A. Zaman
 *
 * @author Tomasz Przedzinski
 * @date 18th October 2010
 */

namespace Photospp
{

class PhotosRandom
{
public:
	/* Change the seed. Default is s1=1802 and s2=9373
	   These values must be in range [0,31327] and [0,30080] respectively */
	static void setSeed(int s1,int s2);

	/* Initialization routine. Must be called at least once before
	   the generator can be used. */
	static void initialize();

	/* Uniform distribution between 0 and 1 */
	static double randomReal();

protected:
	static bool init;
	static int iseed[2];
	static int i97;
	static int j97;
	static double uran[97];
	static double cran;
	static const double cdran;
	static const double cmran;
};

} // namespace Photospp
#endif
