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
#include "TTimeStamp.h"

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

    static TH1* ExpandTimeAxis(const TH1& h, Int_t expansionTime, Int_t timeResolution=-1);

    static Bool_t GetTimeOffset(const TH1& h, TTimeStamp& origin);

private:

   void ComputeNumberOfPads();

   Bool_t FillNumberOfPads();

	static int GetTimeResolution(const TH1& h);

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
