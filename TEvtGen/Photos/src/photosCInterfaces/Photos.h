#ifndef _Photos_h_included_
#define _Photos_h_included_

/**
 * @class Photos
 *
 * @brief Controls the configuration and initialization of Photos.
 *
 * This is the main configuration class for Photos C++ Interface.
 * It is also used for invoking methods for processing single particle or branch.
 *
 * @author Nadia Davidson
 * @date 16th June 2008
 */
#include <stdarg.h>
#include <vector>
#include "PhotosParticle.h"
#include "PhotosRandom.h"
#include "f_Init.h"
// WARNING: VARIANT B of phase space generation was not tested in C++ for all
// options of the program initialization.. 
//#define VARIANTB true
using std::vector;
using std::pair;

namespace Photospp
{

class PhotosParticle;

class Photos
{
public:
	static const int VER_MAJOR=3, VER_MINOR=56;
	static const int DAT_DAY  =3, DAT_MONTH=5, DAT_YEAR=14;
public:

	/** Initalize Photos with the parameters previously set via the
	   setter methods */
	static void initialize();

	/** Prints info on  Photos initialization (reinitialization)
	   status */
	static void iniInfo();

	/** Process decay of single particle */
	static void processParticle(PhotosParticle *p);
	/** Process decay of whole decay branch starting from given particle */
	static void processBranch(PhotosParticle *p);

	/** Suppress processing of a single decay */
	static void suppressBremForDecay (int count, int motherID, ... );
	/** Suppress processing of whole decay branch */
	static void suppressBremForBranch(int count, int motherID, ... );

	/** Suppress all processing. Only forced decays will be processed. */
	static void suppressAll()                      { isSuppressed=true; }

	/** Force processing of a single decay */
	static void forceBremForDecay (int count, int motherID, ... );

	/** Force processing of a whole decay branch */
	static void forceBremForBranch(int count, int motherID, ... );

  /** Block emissions id decays pi0 and K_L -> gamma e+ e- 
      1           = no suppression
      2 (default) = suppressed emissions in K_L -> gamma e+ e- ... and all pi0 decays */
  static void IPHEKL_setPi0KLnoEmission(int m);

  static bool IPHQRK_setQarknoEmission(int MODCOR, int PDGID);
  
  /** If event record allows it, create history entries of particles
      before Photos processing */
  static void createHistoryEntries(bool flag, int status);

  /** Ignore particles with given status code */
  static void ignoreParticlesOfStatus(int status);

  /** Remove 'status' from the list of ignored status codes */
  static void deIgnoreParticlesOfStatus(int status);
  
  /** Returns 'true' if status code is ignored */
  static bool isStatusCodeIgnored(int status);
public:
  /** Substitute build-in generator with external one */
  static void setRandomGenerator( double (*gen)() );

	/** Seed for RANMAR used by fortran part of the Photos */
	static void setSeed(int iseed1, int iseed2)    { PhotosRandom::setSeed(iseed1,iseed2); }

	/** Maximum interference weight */
	static void maxWtInterference(double interference) { phokey_.fint=interference; }

	/** Minimal energy (in units of decaying particle mass) for photons to be explicitly generated */
	static void setInfraredCutOff(double cut_off)  { phocop_.xphcut=cut_off; }

	/** Coupling constant alpha QED */
	static void setAlphaQED(double alpha)          { phocop_.alpha=alpha; }

	/** Key for interference, matrix element weight */
	static void setInterference(bool interference) { phokey_.interf=(int)interference; }

	/** Set double bremsstrahlung generation */
	static void setDoubleBrem(bool doub)           { phokey_.isec=(int)doub; }

	/** Set bremsstrahlung generation up to multiplicity of 4 */
	static void setQuatroBrem(bool quatroBrem)     { phokey_.itre=(int)quatroBrem; }

	/* Key for partial effects of  matrix element (in leptonic W decays) */
	static void setCorrectionWtForW(bool corr) { phokey_.ifw=(int)corr; }

	/** Set exponentiation mode */
	static void setExponentiation(bool expo);

	/** Switch for complete effects of matrix element (in  scalar  to 2 scalars decays) */
	static void setMeCorrectionWtForScalar(bool corr);

	/** Switch for complete effects of matrix element (in leptonic W decays) */
	static void setMeCorrectionWtForW(bool corr);

	/** Switch for complete effects of matrix element (in leptonic Z decays) */
	static void setMeCorrectionWtForZ(bool corr);

	/** Set photon emission in top pair production in quark (gluon) pair annihilation */
	static void setTopProcessRadiation(bool top)         { phokey_.iftop=(int)top; }

	/* Set if PHOTOS should stop at critical error. True by default.
	   WARNING: These stops are an essential source of debugging information flow
                    from event record to PHOTOS algorithm. Never switch it off! The only exception:
	            you have checked your set-up including particular physics initialization
                    with the substantially large sample and you submit large production. */
	static void setStopAtCriticalError(bool stop);

	/** Initialize kinematic corrections */
	static void initializeKinematicCorrections(int flag) { PHCORK(flag); }

	/** Force mass value to be sqrt(e^2-p^2) for all particle momenta
	    taken from event record. May be important for numerical stability.
		May lead to faulty results due to rounding errors for
		hiper-relativistic electron, for example. */
	static void forceMassFrom4Vector(bool flag) { massFrom4Vector=flag; }
	
  /** When particles with PDGID and -PDGID will be processed by Photos,
      their mass value will be taken from event record instead of being
      calculated from 4-vector.

      This works only if 'forceMassFrom4Vector' is set to 'true' (default)      
      This routine may be executed several times with different PDGID values. */
  static void forceMassFromEventRecord(int pdgid);

  /** When particles with PDGID and -PDGID will be processed by Photos,
      their mass value will be given by user instead of being calculated
      from 4-vector.

      This works only if 'forceMassFrom4Vector' is set to 'true' (default)
      This routine may be executed several times with different PDGID values. */  
  static void forceMass(int pdgid, double mass);

	/** set energy momentum conservation threshold */
	static void setMomentumConservationThreshold(double threshold){momentum_conservation_threshold=threshold; }

public:
	/** Is in suppressed mode */
	static bool isSuppressed;

	/** Is mass from 4-vector or from event record */
	static bool massFrom4Vector;
	
	/** List of suppressed decays */
	static vector<vector<int>* >    *supBremList;

	/** List of forced decays */
	static vector<vector<int>* >    *forceBremList;

	/** List of forced mass values */
	static vector<pair<int,double>* > *forceMassList;
  
  /** List of ignored status codes */
	static vector<int >             *ignoreStatusCodeList;

 	/** Threshold for momentum conservation check */
	static double momentum_conservation_threshold;

	/** Flag for complete effects of matrix element (in scalars decays) */
	static bool meCorrectionWtForScalar;

	/** Flag for complete effects of matrix element (in leptonic Z decays) */
	static bool meCorrectionWtForZ;
	
	/** Flag for complete effects of matrix element (in leptonic W decays) */
	static bool meCorrectionWtForW;
  
  /** Flag for creating historic entries */
  static bool isCreateHistoryEntries;

  /** Status of history entries */
  static int  historyEntriesStatus;

  /** Pointer to random generator function */
  static double (*randomDouble)();
public:
	/** Get instance of Photos */
	Photos& getInstance() { return _instance; }
private:
	/* Singleton: only one instance allowed.
	   Constructor sets default values of PHOTOS parameters */
	 Photos();
	~Photos() {}
	Photos(const Photos&);
	Photos& operator=(const Photos&);
	static Photos _instance;
};

} // namespace Photospp
#endif

