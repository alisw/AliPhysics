///
/// \file AliFemtoNSigmaFilter.cxx
///


#include "AliFemtoNSigmaFilter.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoNSigmaFilter);
  /// \endcond
#endif

#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <algorithm>  // std::sort

using std::cout;
using std::cerr;
using std::endl;

//________________________________________________________________________________________________________________
AliFemtoNSigmaFilter::AliFemtoNSigmaFilter():

  fMomMaxPreferTPC(0.5)

  , fTPCAndTOFCutCollection()
  , fTPCCutCollection()
  , fTOFCutCollection()

  , fAbsoluteMomMax(0.0)
  , fImproperConfig(false)
  , fOverrideImproperConfig(false)
{
  // no-op
}


//________________________________________________________________________________________________________________
AliFemtoNSigmaFilter::AliFemtoNSigmaFilter(const AliFemtoNSigmaFilter& aFilter) :
  fMomMaxPreferTPC(aFilter.fMomMaxPreferTPC)

  , fTPCAndTOFCutCollection(aFilter.fTPCAndTOFCutCollection)
  , fTPCCutCollection(aFilter.fTPCCutCollection)
  , fTOFCutCollection(aFilter.fTOFCutCollection)

  , fAbsoluteMomMax(aFilter.fAbsoluteMomMax)
  , fImproperConfig(aFilter.fImproperConfig)
  , fOverrideImproperConfig(aFilter.fOverrideImproperConfig)
{
  //no-op
}

//________________________________________________________________________________________________________________
AliFemtoNSigmaFilter& AliFemtoNSigmaFilter::operator=(const AliFemtoNSigmaFilter& aFilter)
{
  if(this == &aFilter) {return *this;}

  fMomMaxPreferTPC = aFilter.fMomMaxPreferTPC;

  fTPCAndTOFCutCollection = aFilter.fTPCAndTOFCutCollection;
  fTPCCutCollection = aFilter.fTPCCutCollection;
  fTOFCutCollection = aFilter.fTOFCutCollection;

  fAbsoluteMomMax = aFilter.fAbsoluteMomMax;
  fImproperConfig = aFilter.fImproperConfig;
  fOverrideImproperConfig = aFilter.fOverrideImproperConfig;

  return *this;
}

//________________________________________________________________________________________________________________
AliFemtoNSigmaFilter::~AliFemtoNSigmaFilter()
{
  //no-op
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::CheckForOverlap(VecNSigmaCuts &aCollection)
{
  bool tImproperConfig = false;

  //called within SortCollection

  //First, make sure all defined momentum values are positive
  //Also, make sure that pMin < pMax
  for (VecNSigmaCuts::iterator cut = aCollection.begin(); cut != aCollection.end(); ++cut) {
    if (cut->pMin < 0.0) {
      tImproperConfig = true;
      cerr << "W-AliFemtoNSigmaFilter: WARNING - pMin is set to a negative value!!!!!\n";
    }
    if (cut->pMax < 0.0) {
      tImproperConfig = true;
      cerr << "W-AliFemtoNSigmaFilter: WARNING - pMax is set to a negative value!!!!!\n";
    }
    if (cut->pMin >= cut->pMax) {
      tImproperConfig = true;
      cerr << "W-AliFemtoNSigmaFilter: WARNING - pMin >= pMax, making this bin useless!!!!!" << endl;
    }
  }

  // Now, check for overlapping bins (remembering that sorting has already occurred)
  for (UInt_t i = 1; i < aCollection.size(); i++) { //note, starting at i=1, not i=0
    if (aCollection[i].pMin < aCollection[i-1].pMax) {
      tImproperConfig = true;
      cerr << "WARNING:  Overlapping bins created!!!!!" << endl;
      if (fOverrideImproperConfig) {
        cerr << "User set fOverrideImproperConfig = true.  Therefore, in the overlap region, the bin"
                "with smallest pMin will be used" << endl;
      }
    }
  }

  if (tImproperConfig && !fOverrideImproperConfig) {
    cerr << "Improper configuration detected, and user did not set fOverrideImproperConfig=true."
            "Therefore, all candidates will fail the cut" << endl;
  }

  fImproperConfig = tImproperConfig;
}



//________________________________________________________________________________________________________________
bool AliFemtoNSigmaFilter::ComparePMin(const AliFemtoNSigmaFilter::NSigmaCut &aFirst, const AliFemtoNSigmaFilter::NSigmaCut &aSecond)
{
  return (aFirst.pMin < aSecond.pMin);
}
//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::SortCollection(VecNSigmaCuts &aCollection)
{
  sort(aCollection.begin(), aCollection.end(), ComparePMin);
  CheckForOverlap(aCollection);
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::AddTPCAndTOFCut(double aPMin, double aPMax, double aNSigmaValueTPC, double aNSigmaValueTOF)
{
  NSigmaCut tNSigmaCut = {
    kTPCAndTOF,
    aPMin,
    aPMax,
    aNSigmaValueTPC,
    aNSigmaValueTOF
  };

  fTPCAndTOFCutCollection.push_back(tNSigmaCut);
  SortCollection(fTPCAndTOFCutCollection);

  //if the maximum pMax in fTPCAndTOFCutCollection is greater than fAbosoluteMomMax, reset fAbosoluteMomMax equal to it
  if (fTPCAndTOFCutCollection.back().pMax > fAbsoluteMomMax) {
    fAbsoluteMomMax = fTPCAndTOFCutCollection.back().pMax;
  }
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::AddTPCCut(double aPMin, double aPMax, double aNSigmaValueTPC)
{
  NSigmaCut tNSigmaCut = {kTPC, aPMin, aPMax, aNSigmaValueTPC, 0.0};
  fTPCCutCollection.push_back(tNSigmaCut);
  SortCollection(fTPCCutCollection);

  //if the maximum pMax in fTPCCutCollection is greater than fAbosoluteMomMax, reset fAbosoluteMomMax equal to it
  if (fTPCCutCollection.back().pMax > fAbsoluteMomMax) {
    fAbsoluteMomMax = fTPCCutCollection.back().pMax;
  }
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::AddTOFCut(double aPMin, double aPMax, double aNSigmaValueTOF)
{
  NSigmaCut tNSigmaCut = {kTOF, aPMin, aPMax, 0.0, aNSigmaValueTOF};
  fTOFCutCollection.push_back(tNSigmaCut);
  SortCollection(fTOFCutCollection);

  //if the maximum pMax in fTOFCutCollection is greater than fAbosoluteMomMax, reset fAbosoluteMomMax equal to it
  if (fTOFCutCollection.back().pMax > fAbsoluteMomMax) {
    fAbsoluteMomMax = fTOFCutCollection.back().pMax;
  }
}



//_______________________________________________________________________________________________________
int AliFemtoNSigmaFilter::FindBinOfInterest(double aMom, const VecNSigmaCuts &aCollection) const
{
  for (unsigned int i = 0; i < aCollection.size(); i++) {
    if ( (aCollection[i].pMin <= aMom) && (aMom < aCollection[i].pMax) ) {
      return i;
    }
  }

  return -1; // implies no cut exists in this collection for this particular momentum value
}


//_______________________________________________________________________________________________________
bool AliFemtoNSigmaFilter::Pass(double aMom, double aNSigmaTPC, double aNSigmaTOF)
{
  // If the user configures the AliFemtoNSigmaFilter incorrectly (in the opinion of the authors), all
  // candidates fail.
  // However, if the user finds some novel way of using the filter with a configuration we deem bad, the
  // fOverrideImproperConfig allows one to negate this check
  if (fImproperConfig && !fOverrideImproperConfig) {
    return false;
  }

  // Cut out any candidates with momentum greater than the range of all of the cut collection
  if (aMom > fAbsoluteMomMax) {
    return false;
  }

  // Check if N-Sigma values are available (i.e. > -999)
  bool tExistTPC = (aNSigmaTPC > -999),
       tExistTOF = (aNSigmaTOF > -999);

  // If both TPC and TOF nsigma values do not exist, the cut should fail
  if (!tExistTPC && !tExistTOF) {
    return false;
  }

  // Next, find which collections need to be considered for this momentum, and which particular momentum
  // bins are needed. If -1 is returned, there is no cut for this momentum & detector pair
  const Int_t iTPCAndTOF = FindBinOfInterest(aMom, fTPCAndTOFCutCollection),
                    iTPC = FindBinOfInterest(aMom, fTPCCutCollection),
                    iTOF = FindBinOfInterest(aMom, fTOFCutCollection);

  // Check if TPCAndTOF cut exists for this momentum (iTPCAndTOF > -1)
  // and check if nsigma information is available for both (i.e. tExistTPC=true && tExistTOF=true)
  // If so, apply the cuts
  if (iTPCAndTOF > -1 && tExistTPC && tExistTOF) {
    return (TMath::Abs(aNSigmaTPC) < fTPCAndTOFCutCollection[iTPCAndTOF].nSigmaValueTPC
         && TMath::Abs(aNSigmaTOF) < fTPCAndTOFCutCollection[iTPCAndTOF].nSigmaValueTOF);
  }

  // If user has implemented separate TPC and TOF nsigma cuts for this momentum value, and both
  // experimental values are valid, there exists an ambiguity in which detector to use.  The ambiguity is
  // resolved by the fMomMaxPreferTPC member.
  // See the documentation in the header file for a more complete description
  else if (iTPC > -1 && tExistTPC && iTOF > -1 && tExistTOF) {
    if (aMom <= fMomMaxPreferTPC) { //prefer TPC
      return TMath::Abs(aNSigmaTPC) < fTPCCutCollection[iTPC].nSigmaValueTPC;
    }
    else { // prefer TOF
      return TMath::Abs(aNSigmaTOF) < fTOFCutCollection[iTOF].nSigmaValueTOF;
    }

  }

  // Check if TPC cut exists for this momentum (iTPC > -1) and if nsigmaTPC information is available
  // If so, apply the cut
  if (iTPC > -1 && tExistTPC) {
    return TMath::Abs(aNSigmaTPC) < fTPCCutCollection[iTPC].nSigmaValueTPC;
  }

  // Check if TOF cut exists for this momentum (iTOF > -1) and if nsigmaTOF information is available
  // If so, apply the cut
  else if (iTOF > -1 && tExistTOF) {
    return TMath::Abs(aNSigmaTOF) < fTOFCutCollection[iTOF].nSigmaValueTOF;
  }

  // If user does not define any nsigma cuts for the particular momentum of interest, return false
  return false;

}
