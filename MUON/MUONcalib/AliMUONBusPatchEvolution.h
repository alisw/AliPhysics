#ifndef ALIMUONBUSPATCHEVOLUTION_H
#define ALIMUONBUSPATCHEVOLUTION_H

/// \ingroup calib
/// \class AliMUONBusPatchEvolution
/// \brief Utility class to massage the output of the MCHBPEVO DA
///
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#include <map>
#include <vector>

class AliMergeableCollection;
class TH1;

class AliMUONBusPatchEvolution : public TObject
{
public:

	AliMUONBusPatchEvolution(AliMergeableCollection& hc);
	AliMUONBusPatchEvolution(AliMergeableCollection& hc, const std::map<int,int>& nofPadsPerBusPatch);

	void Augment();

	Bool_t GetFaultyBusPatches(int timeResolution,
			int requiredEvents,
			float occupancyThreshold,
			std::map<int,double>& faultyBusPatchOccupancies);

	void GetTimeResolutions(std::vector<int>& timeResolutions);

    void ShrinkTimeAxis();

private:

   void ComputeNumberOfPads();

   Bool_t FillNumberOfPads();

	int GetTimeResolution(TH1& h);

   void GroupByStation(int timeResolution);
   void GroupByChamber(int timeResolution);
   void GroupByDDL(int timeResolution);
   void GroupByDE(int timeResolution);

   void Normalize();

  void NumberOfPadsFromHistosToMap(std::map<int,int>& nofPadsPerBusPatch);

private:
  AliMergeableCollection& fBPEVO;
  std::map<int,int> fNofPads;

  ClassDef(AliMUONBusPatchEvolution,1)

};

#endif
