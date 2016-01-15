///////////////////////////////////////////////////////////////////////////
// AliFemtoNSigmaFilter:                                                 //
///////////////////////////////////////////////////////////////////////////


#include "AliFemtoNSigmaFilter.h"

#ifdef __ROOT__
ClassImp(AliFemtoNSigmaFilter)
#endif




//________________________________________________________________________________________________________________
//****************************************************************************************************************
//________________________________________________________________________________________________________________


//________________________________________________________________________________________________________________
AliFemtoNSigmaFilter::AliFemtoNSigmaFilter():

  fMomMaxPreferTPC(0.5),

  fTPCAndTOFCutCollection(0),
  fTPCCutCollection(0),
  fTOFCutCollection(0),

  fAbsoluteMomMax(0.),

  fImproperConfig(false),
  fOverrideImproperConfig(false)


{


}



//________________________________________________________________________________________________________________
AliFemtoNSigmaFilter::~AliFemtoNSigmaFilter()
{
  /* noop */
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::CheckForOverlap(VecNSigmaCuts &aCollection)
{
  bool tImproperConfig = false;

  //called within SortCollection

  //First, make sure all defined momentum values are positive
  //Also, make sure that pMin < pMax
  for(unsigned int i=0; i<aCollection.size(); i++)
  {
    if(aCollection[i].pMin < 0.) 
    {
      tImproperConfig = true;
      cerr << "WARNING:  pMin is set to a negative value!!!!!" << endl;
    }
    if(aCollection[i].pMax < 0.) 
    {
      tImproperConfig = true;
      cerr << "WARNING:  pMax is set to a negative value!!!!!" << endl;
    }
    if(aCollection[i].pMin > aCollection[i].pMax) 
    {
      tImproperConfig = true;
      cerr << "WARNING:  pMin > pMax, making this bin useless!!!!!" << endl;
    }
  }

  //Now, check for overlapping bins (remembering that sorting has already occurred)
  for(unsigned int i=1; i<aCollection.size(); i++) //note, starting at i=1, not i=0
  {
    if(aCollection[i].pMin < aCollection[i-1].pMax) 
    {
      tImproperConfig = true;
      cerr << "WARNING:  Overlapping bins created!!!!!" << endl;
      if(fOverrideImproperConfig)
      {
        cerr << "User set fOverrideImproperConfig = true.  Therefore, in the overlap region, the bin with smallest pMin will be used" << endl;
      }
    }
  }

  if(tImproperConfig && !fOverrideImproperConfig) 
  {
    cerr << "Improper configuration detected, and user did not set fOverrideImproperConfig=true.  Therefore, all candidates will fail the cut" << endl;
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
  NSigmaCut tNSigmaCut = {kTPCAndTOF,aPMin,aPMax,aNSigmaValueTPC,aNSigmaValueTOF};
  fTPCAndTOFCutCollection.push_back(tNSigmaCut);
  SortCollection(fTPCAndTOFCutCollection);

  //if the maximum pMax in fTPCAndTOFCutCollection is greater than fAbosoluteMomMax, reset fAbosoluteMomMax equal to it
  if(fTPCAndTOFCutCollection.back().pMax > fAbsoluteMomMax) fAbsoluteMomMax = fTPCAndTOFCutCollection.back().pMax;
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::AddTPCCut(double aPMin, double aPMax, double aNSigmaValueTPC)
{
  NSigmaCut tNSigmaCut = {kTPC,aPMin,aPMax,aNSigmaValueTPC,0.};
  fTPCCutCollection.push_back(tNSigmaCut);
  SortCollection(fTPCCutCollection);

  //if the maximum pMax in fTPCCutCollection is greater than fAbosoluteMomMax, reset fAbosoluteMomMax equal to it
  if(fTPCCutCollection.back().pMax > fAbsoluteMomMax) fAbsoluteMomMax = fTPCCutCollection.back().pMax;
}


//________________________________________________________________________________________________________________
void AliFemtoNSigmaFilter::AddTOFCut(double aPMin, double aPMax, double aNSigmaValueTOF)
{
  NSigmaCut tNSigmaCut = {kTOF,aPMin,aPMax,0.,aNSigmaValueTOF};
  fTOFCutCollection.push_back(tNSigmaCut);
  SortCollection(fTOFCutCollection);

  //if the maximum pMax in fTOFCutCollection is greater than fAbosoluteMomMax, reset fAbosoluteMomMax equal to it
  if(fTOFCutCollection.back().pMax > fAbsoluteMomMax) fAbsoluteMomMax = fTOFCutCollection.back().pMax;
}



//________________________________________________________________________________________________________________
int AliFemtoNSigmaFilter::FindBinOfInterest(double aMom, VecNSigmaCuts &aCollection)
{
  for(unsigned int i=0; i<aCollection.size(); i++)
  {
    if( (aCollection[i].pMin <= aMom) && (aCollection[i].pMax > aMom) ) return i;
  }

  return -1; //implies no cut exists in this collection for this particular momentum value
}


//________________________________________________________________________________________________________________
bool AliFemtoNSigmaFilter::Pass(double aMom, double aNSigmaTPC, double aNSigmaTOF)
{
  //If the user configures the AliFemtoNSigmaFilter incorrectly (in the opinion of the authors), all candidates fail
  //However, if the user finds some novel way of using the filter with a configuration we deem bad, the fOverrideImproperConfig
  //allows one to negate this check
  if(fImproperConfig && !fOverrideImproperConfig) return false;

  //Cut out any candidates with momentum greater than the range of all of the cut collection
  if(aMom > fAbsoluteMomMax) return false;

  //First, check if nsigma values are available (i.e. > -999)
  bool tExistTPC = true;
  bool tExistTOF = true;

  if(aNSigmaTPC < -999) tExistTPC = false;
  if(aNSigmaTOF < -999) tExistTOF = false;

  //If both TPC and TOF nsigma values do not exist, the cut should fail
  if(!tExistTPC && !tExistTOF) return false;

  //Next, find which collections need to be considered for this momentum, and which particular momentum bins are needed
  //If the returned value = -1, this implies no cut exists for this particular momentum value and detector
  int iTPCAndTOF = FindBinOfInterest(aMom,fTPCAndTOFCutCollection);
  int iTPC = FindBinOfInterest(aMom,fTPCCutCollection);
  int iTOF = FindBinOfInterest(aMom,fTOFCutCollection);

  // Check if TPCAndTOF cut exists for this momentum (iTPCAndTOF > -1)
  // and check if nsigma information is available for both (i.e. tExistTPC=true && tExistTOF=true)
  // If so, apply the cuts
  if(iTPCAndTOF > -1 && tExistTPC && tExistTOF )
  {
    if(TMath::Abs(aNSigmaTPC) < fTPCAndTOFCutCollection[iTPCAndTOF].nSigmaValueTPC 
       && TMath::Abs(aNSigmaTOF) < fTPCAndTOFCutCollection[iTPCAndTOF].nSigmaValueTOF) {return true;}
    else {return false;}
  }

  // If user has implemented separate TPC and TOF nsigma cuts for this momentum value, and both experimental values are valid, there 
  // exists an ambiguity in which detector to use.  The ambiguity is resolved by the fMomMaxPreferTPC member.
  // See the documentation in the header file for a more complete description
  else if(iTPC > -1 && tExistTPC && iTOF > -1 && tExistTOF)
  {
    if(aMom <= fMomMaxPreferTPC) //prefer TPC
    {
      if(TMath::Abs(aNSigmaTPC) < fTPCCutCollection[iTPC].nSigmaValueTPC) {return true;}
      else {return false;}
    }
    else //prefer TOF
    {
      if(TMath::Abs(aNSigmaTOF) < fTOFCutCollection[iTOF].nSigmaValueTOF) {return true;}
      else {return false;}
    }

  }

  // Check if TPC cut exists for this momentum (iTPC > -1) and if nsigmaTPC information is available
  // If so, apply the cut
  if(iTPC > -1 && tExistTPC)
  {
    if(TMath::Abs(aNSigmaTPC) < fTPCCutCollection[iTPC].nSigmaValueTPC) {return true;}
    else {return false;}
  }

  // Check if TOF cut exists for this momentum (iTOF > -1) and if nsigmaTOF information is available
  // If so, apply the cut
  else if(iTOF > -1 && tExistTOF)
  {
    if(TMath::Abs(aNSigmaTOF) < fTOFCutCollection[iTOF].nSigmaValueTOF) {return true;}
    else {return false;}
  }


  // If user does not define any nsigma cuts for the particular momentum of interest, return false
  return false;

}

