/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>

#include "AliLog.h"

#include "AliTimeRangeMasking.h"

templateClassImp(AliTimeRangeMask)

template<typename time_type, typename bitmap_type>
std::vector<std::string> AliTimeRangeMask<time_type, bitmap_type>::fMaskReasonNames = {"BadTPCPID", "T1", "T2"};

template<typename time_type, typename bitmap_type>
AliTimeRangeMask<time_type, bitmap_type>::AliTimeRangeMask()
  : TObject(),
    fStart(),
    fEnd(),
    fMaskReasons()
{
}

template<typename time_type, typename bitmap_type>
AliTimeRangeMask<time_type, bitmap_type>::AliTimeRangeMask(AliTimeRangeMask const& other)
  : TObject(other),
    fStart(other.fStart),
    fEnd(other.fEnd),
    fMaskReasons(other.fMaskReasons)
{ }

template<typename time_type, typename bitmap_type>
AliTimeRangeMask<time_type, bitmap_type>::AliTimeRangeMask(time_type start, time_type end, bitmap_type reasons)
  : TObject(),
    fStart(start),
    fEnd(end),
    fMaskReasons(reasons)
{ }

template<typename time_type, typename bitmap_type>
AliTimeRangeMask<time_type, bitmap_type>::AliTimeRangeMask(time_type time)
  : TObject(),
    fStart(time),
    fEnd(time),
    fMaskReasons()
{
}

template<typename time_type, typename bitmap_type>
std::string AliTimeRangeMask<time_type, bitmap_type>::CollectMaskReasonNames() const
{
  std::string reasons;
  for (int iReason = 0; iReason < kNReasons; ++iReason) {
    if (HasReason(MaskReason(iReason))) {
      if (reasons.length()) {
        reasons += " | ";
      }
      reasons += fMaskReasonNames[iReason];
    }

  }

  return reasons;
}

template<typename time_type, typename bitmap_type>
void AliTimeRangeMask<time_type, bitmap_type>::Print(Option_t* option) const
{
  std::cout << fStart << " - " << fEnd << ": " << CollectMaskReasonNames() << "\n";
}

template class AliTimeRangeMask<ULong64_t, UShort_t>;

templateClassImp(AliTimeRangeMasking)

template<typename time_type, typename bitmap_type>
AliTimeRangeMasking<time_type, bitmap_type>::AliTimeRangeMasking()
  : TObject(),
    fArrTimeRanges("AliTimeRangeMask<ULong64_t, UShort_t>", 10)
{
}


template<typename time_type, typename bitmap_type>
AliTimeRangeMask<time_type, bitmap_type>* AliTimeRangeMasking<time_type, bitmap_type>::AddTimeRangeMask(time_type start, time_type end, bitmap_type reasons)
{
  if ( const auto* range = FindTimeRangeMask(start))  {
    const std::string reasonsString = range->CollectMaskReasonNames();
    AliErrorF("Start time %llu already in range [%llu, %llu]: %s", 
        start, range->GetStart(), range->GetEnd(), reasonsString.data());
    return nullptr;
  }

  if ( const auto* range = FindTimeRangeMask(end))  {
    const std::string reasonsString= range->CollectMaskReasonNames();
    AliErrorF("End time %llu already in range [%llu, %llu]: %s", 
        end, range->GetStart(), range->GetEnd(), reasonsString.data());
    return nullptr;
  }

  return new(fArrTimeRanges[fArrTimeRanges.GetEntriesFast()]) AliTimeRangeMask<time_type, bitmap_type>(start, end, reasons);
}

template<typename time_type, typename bitmap_type>
AliTimeRangeMask<time_type, bitmap_type>* AliTimeRangeMasking<time_type, bitmap_type>::FindTimeRangeMask(time_type time) const
{
  for (auto o : fArrTimeRanges) {
    auto const val = (AliTimeRangeMask<time_type, bitmap_type>*)o;
    if ( *val == time ) return val;
  }
  return nullptr;
}

template<typename time_type, typename bitmap_type>
void AliTimeRangeMasking<time_type, bitmap_type>::Print(Option_t* option) const
{
  for (const auto& range : fArrTimeRanges) {
    range->Print();
  }
}

template class AliTimeRangeMasking<ULong64_t, UShort_t>;
