////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliFemtoNSigmaFilter:                                                                                  //
//                                                                                                        //
//  Authors:  Jesse Buxton (jesse.thomas.buxton@cern.ch)                                                  //
//            Andrew Kubera (andrew.michael.kubera@cern.ch)                                               //
//            Jai Salzwedel (Jai.Salzwedel@cern.ch)                                                       //
//                                                                                                        //
//  The purpose of this class is to aid in the selection of particles. In particular, this class is       //
//  intended to work with functions such as AliFemtoV0TrackCut::IsKaonNSigma etc.  This class is designed //
//  to work specifically with the AliFemtoV0TrackCut class to help find the correct daughters for V0      //
//  reconstruction.  If desired, these methods could easily be added to classes such as                   //
//  AliFemtoESDTrackCut in AliFemtoUser.  An updated AliFemtoV0TrackCut including these methods will      //
//  likely later be merged, but for now one should use the AliFemtoV0TrackCutNSigmaFilter class in the    //
//  AliFemtoUser directory.                                                                               //
////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // Note:  When I speak of NSigma values being available/valid, this simply means NSigma > -999.
  // Note:  Momentum variables are named with P (or Mom), but AliFemtoV0TrackCutNSigmaFilter will allow the user to cut on NSigma in p or pT. (14/12/2015 not yet implemented)
  // Note:  Struct variables have a lower case first letter (excluding "a").
  // Note:  Function parameters begin with "a".

  // ------------------------------------------------------------ General overview: ----------------------------------------------------------------------------------------------

  // The user creates a number of (struct) NSigmaCut objects through the Add(TPCAndTOF/TPC/TOF)Cut methods.  Each NSigmaCut object represents an interval in momentum, [pMin,pMax),
  // and contains NSigma cut values for the TPC and/or TOF.  The cuts contained in fTPCAndTOFCutCollection (and created by the AddTPCAndTOFCut method) are implemented when both
  // the experimental TPC and TOF NSigma values for a given test particle are available.  The cuts contained in fTPC(TOF)CutCollection (and created by AddTPC(TOF)Cut) are
  // implemented when the experimental TPC(TOF) NSigma value is available.  The TPCAndTOF cuts are preferred over TPC(TOF) cuts (i.e. if the user sets up a TPCAndTOF cut
  // and a TPC cut in the same momentum region, if both experimental NSigmaTPC and NSigmaTOF values are available, the TPCAndTOF cut will be utilized).

  // Note: When setting the NSigma cuts, both pMin and pMax must be supplied.  If, for instance, you would like to extend your last cut object out to infinity, simply choose
  //       an appropriately large number.  Also, pMin and pMax must be positive numbers.


  // -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  // --*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*

  // Note:  A TPCAndTOFCut can only pass if both the TPC and TOF NSigma values are available.  More specifically: If, for a particular region in momentum, the user sets
  //        ONLY a TPCAndTOFCut, the test particle will fail if one of the experimental NSigma values (TPC or TOF) is invalid (and, obviously, when both values are valid 
  //        but one or both fail the NSigma cuts).  See Example 3 (below).
  //        Therefore, if you would like to, for instance, use the TPC(TOF) when the TOF(TPC) NSigma value is not available, you must also set a TPC(TOF)Cut for the desired
  //        momentum regime.  For instance, see Example 1 and Example 3.

  // Note:  Constructing a TPCAndTOFCut, TPCCut, and TOFCut for a particular region in momemtum will have the following behavior:
  //        1.) If both experimental NSigmaTPC and NSigmaTOF are available, evaluate the TPCAndTOFCut.  If BOTH NSigmaTPC and NSigmaTOF pass their respective cuts, the test
  //            particle passes.  Otherwise, it fails.
  //        2.) If only NSigmaTPC is valid, evaluate the TPCCut.  The test particle either passes or fails.
  //        3.) If only NSigmaTOF is valid, evaluate the TOFCut.  The test particle either passes or fails.
  //        4.) If neither NSigmaTPC nor NSigmaTOF are available, the test particle fails.

  // --*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*


  // --@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@

  // Note:  fMomMaxPreferTPC information:
  //          For a particular region in momentum, if the user defines an individual TPCCut and an individual TOFCut (instead of using a TPCAndTOFCut), there exists an 
  //          ambiguity in which cut to implement when both experimental NSigmaTPC and NSigmaTOF values are valid.  In this case, the member variable fMomMaxPreferTPC 
  //          dictates which detector to use in the cut: for momenta <= fMomMaxPreferTPC, use the TPCCut, otherwise use the TOFCut.  See Example 4.

  //          The default value will be set to 0.5 GeV/c, but this can be adjusted by the user.  For instance, when dealing with protons, this value can
  //          probably be extended up to ~1 GeV/c. 
  //          In most cases, this variable will go unused.  This is mainly used as a type of safeguard to resolve ambiguity, but I suppose some users could 
  //          desire this additional functionality.  

  // Note:  Constructing a TPCCut and TOFCut for a particular region in momentum will have the following behavior:
  //        1.) If both experimental NSigmaTPC and NSigmaTOF values are valid
  //             i.) If the test particle's momentum is <= fMomMaxPreferTPC, evaluate the TPCCut.  The test particle either passes or fails.
  //             ii.)If the test particle's momentum is > fMomMaxPreferTPC, evaluate the TOFCut.  The test particle either passes of fails.
  //        2.) If NSigmaTPC value is valid (but NSigmaTOF is not), evaluate the TPCCut.
  //        3.) If NSigmaTOF value is valid (but NSigmaTPC is not), evaluate the TOFCut.
  //        4.) If neither the NSigmaTPC nor the NSigmaTOF value is valid, the particle fails the cut.

  //          In this case,
  //          if both the NSigmaTPC and NSigmaTOF values are valid, the TPC cut will be evaulated if the particle has a momentum <= fMomMaxPreferTPC,
  //          and the TOF cut will be evaluated if the momentum is > fMomMaxPreferTPC.
  //          If the NSigmaTPC value is valid by NSigmaTOF is not, the TPC cut will be evaluated
  //          If the NSigmaTOF value is valid by NSigmaTPC is not, the TOF cut will be evaluated
  //          If neither the NSigmaTPC nor the NSigmaTOF value is valid, the particle will fail the cut
  // --@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@

  //!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!
  // Improper configuration:
  // When building the various cuts, the method CheckForOverlap checks the user's configuration for a few simple errors.  If an error is found, the member fImproperConfig is set the true,
  // and all candidates will fail.  For now, each cut collection is checked for:
  //   1. A negative pMin value
  //   2. A negative pMax value
  //   3. pMin > pMax
  //   4. Overlapping bins in a cut collection
  //      NOTE:  This checks for overlapping bins within a single cut collection (for instance, overlapping bins in the fTPCCutCollection).  Overlapping bins amongst the (3) different
  //             cut collections is allowed, and in many cases, is desired. 

  // The configurations are obviously defined improper by the authors.  However, the user may notice some novel way to use the filter with a configuration we deem "bad".  In this case, the user
  // may set fOverrideImproperConfig=true (via SetOverrideImproperConfig(true)), which will negate the fImproperConfig guard.  We have implemented this guard for a reason.
  // Therefore, the user should have a deep understanding of the functionality of this class, and must understand completely how his ("improper") configuration will function, before overriding 
  // the improper configuration guard.

  // Note:  SetOverrideImproperConfig(bool aOverride) should be called before the cut collections are made.  However, if it is called after, the override will still function properly, but the 
  //        error messages will read as if fOverrideImproperConfig=false
  //!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!--!!


  //--------------------------------------------------------- EXAMPLES ----------------------------------------------------------------------------------
  //Methods to add TPC and TOF cuts
  //Example 1:  For the momentum range 0 < p < 0.4 GeV/c
  //            If you would like nsigmaTPC < 3 && nsigmaTOF < 3 when both values are available, 
  //            but nsigmaTPC < 1 if no TOF available (and implicitly fail if experimental nsigmaTPC is unavailable), you would set the following:
  //              AddTPCAndTOFCut(0.,0.4,3.,3.);
  //              AddTPCCut(0.,0.4,1.);

  //Example 2:  For the momentum range 0.6 < p < 5.0 GeV/c
  //            If you would like nsigmaTPC < 3 && nsigmaTOF < 3 when both values are available, 
  //            but fail if no TPC and/or TOF available, you would set the following:
  //              AddTPCAndTOFCut(0.6,5.0,3.,3.);

  //Example 3:  For the momentum range 0.4 < p < 1.0 GeV/c
  //		If you would like nsigmaTPC < 3 && nsigmaTOF < 4 when both experimental values are valid,
  //		nsigmaTPC < 1 if no TOF
  //		nsigmaTOF < 2 if no TPC
  //		  AddTPCAndTOFCut(0.4,1.0,3.,4.);
  //		  AddTPCCut(0.4,1.0,1.);
  //		  AddTOFCut(0.4,1.0,2.);

  //Example 4:  For the momentum range 0.4 < p < 1.0 GeV/c
  //		nsigmaTPC < 1
  //		nsigmaTOF < 2
  //		For p <= 0.6 GeV/c, use the TPC if available, otherwise use TOF if available (otherwise fail)
  //		For p > 0.6 GeV/c, use the TOF if available, otherwise use the TPC if available (otherwise fail)
  //		  AddTPCCut(0.4,1.0,1.);
  //		  AddTOFCut(0.4,1.0,2.);
  //		  SetMomMaxPreferTPC(0.6);
  //-----------------------------------------------------------------------------------------------------------------------------------------------------



//_______________________________________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________________________________

#ifndef ALIFEMTONSIGMAFILTER_H
#define ALIFEMTONSIGMAFILTER_H

//includes and any constant variable declarations
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>  //std::sort
#include <assert.h>

#include "TString.h"

using std::cerr;
using std::endl;
using std::vector;
using std::sort;

class AliFemtoNSigmaFilter {

public:

  enum DetectorType {kTPC=0, kTOF=1, kTPCAndTOF=2};

  struct NSigmaCut
  {
    DetectorType detectorType;
    double pMin;  //later, user can choose to feed p or pT
    double pMax;  //later, user can choose to feed p or pT
    double nSigmaValueTPC;
    double nSigmaValueTOF;
  };

  typedef vector<NSigmaCut> VecNSigmaCuts;

  //Constructor, destructor, copy constructor, assignment operator
  AliFemtoNSigmaFilter();
  virtual ~AliFemtoNSigmaFilter();

  void CheckForOverlap(VecNSigmaCuts &aCollection);  //if overlapping values, there will exist a region with multiple NSigma values.  This is bad
  static bool ComparePMin(const NSigmaCut &aFirst, const NSigmaCut &aSecond);  //used by SortCollection to sort the cut collections by pMin
  void SortCollection(VecNSigmaCuts &aCollection);  //sort the collections to make iteration easier

  void AddTPCAndTOFCut(double aPMin, double aPMax, double aNSigmaValueTPC, double aNSigmaValueTOF);
  void AddTPCCut(double aPMin, double aPMax, double aNSigmaValueTPC);
  void AddTOFCut(double aPMin, double aPMax, double aNSigmaValueTOF);

  void SetMomMaxPreferTPC(double aMaxMom);
  void SetOverrideImproperConfig(bool aOverride);  //DO NOT set to true without a complete understanding of your "improper" configuration

  int FindBinOfInterest(double aMom, VecNSigmaCuts &aCollection);
  bool Pass(double aMom, double aNSigmaTPC, double aNSigmaTOF);


private:

  double fMomMaxPreferTPC;  // The default value will be set to 0.5 GeV/c.  See documentation at top of current file for more information
  VecNSigmaCuts fTPCAndTOFCutCollection;
  VecNSigmaCuts fTPCCutCollection;
  VecNSigmaCuts fTOFCutCollection;

  double fAbsoluteMomMax; // the maximum momentum defined in all cut collections (the absolute maximum pMax)
                          // this variable will help speed things up by immediately failing any particle with 
                          // a momentum outside of the range of all of the cut collections.

  bool fImproperConfig;  //internal member used to detect improper configurations
  bool fOverrideImproperConfig;  //allows the user (via SetOverrideImproperConfig) to override the fImproperConfig member


#ifdef __ROOT__
  ClassDef(AliFemtoNSigmaFilter, 1)
#endif
};


//inline stuff
inline void AliFemtoNSigmaFilter::SetMomMaxPreferTPC(double aMaxMom) {fMomMaxPreferTPC = aMaxMom;}
inline void AliFemtoNSigmaFilter::SetOverrideImproperConfig(bool aOverride) {fOverrideImproperConfig = aOverride;}

#endif


