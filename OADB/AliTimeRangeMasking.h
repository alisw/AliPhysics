/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#ifndef ALITIMERANGEMASKING_H
#define ALITIMERANGEMASKING_H

/// \file AliTimeRangeMasking.h
/// \brief classes for defining time ranges with a certain mask to be able to cut on
/// \author Jens Wiechula, jens.wiechula@ikf.uni-frankfurt.de

#include <vector>
#include <string>

#include "Rtypes.h"
#include "TObject.h"
#include "TClonesArray.h"

/// \class AliTimeRangeMask
/// A Class for assiciating a bit mask with a time range
template<typename time_type, typename bitmap_type>
class AliTimeRangeMask : public TObject {
  public:
    enum MaskReason { kBadTPCPID = 0, kNReasons };

    AliTimeRangeMask();
    AliTimeRangeMask(AliTimeRangeMask const&);
    AliTimeRangeMask(time_type start, time_type end, bitmap_type reasons = 0);
    AliTimeRangeMask(time_type time);

    void SetReason(MaskReason reason) { SETBIT(fMaskReasons, reason); }
    void ResetReason(MaskReason reason) { CLRBIT(fMaskReasons, reason); }
    Bool_t HasReason(MaskReason reason) const { return TESTBIT(fMaskReasons, reason); }

    time_type GetStart() const { return fStart; }
    time_type GetEnd()   const { return fEnd; }

    bitmap_type GetMaskReasons() const { return fMaskReasons; }

    /// equal operator
    ///
    /// defined as range is contained in
    bool operator== (const AliTimeRangeMask& other) const
    {
      return other.fStart >= fStart && other.fEnd<=fEnd;
    }

    std::string CollectMaskReasonNames() const;

    /// equal operator
    ///
    /// defined as time is contained in range
    bool operator== (const time_type time) const
    {
      return time >= fStart && time<=fEnd;
    }

    bool operator< (const AliTimeRangeMask& other) const
    {
      return fEnd<other.fStart;
    }

    virtual void Print(Option_t* option = "") const;

  private:
    time_type fStart;      ///< start time of masked range
    time_type fEnd;        ///< end time of masked range
    bitmap_type fMaskReasons;  ///< bitmap of masking reasons

    static std::vector<std::string> fMaskReasonNames; //!< mask reason names

    ClassDef(AliTimeRangeMask, 1);
};

/// \class AliTimeRangeMasking
/// A Class for keeping several time ranges with mask of type AliTimeRangeMask
template<typename time_type, typename bitmap_type>
class AliTimeRangeMasking : public TObject {
  public:
    AliTimeRangeMasking();

    AliTimeRangeMask<time_type, bitmap_type>* AddTimeRangeMask(time_type start, time_type end, bitmap_type reasons = {});


    AliTimeRangeMask<time_type, bitmap_type>* FindTimeRangeMask(time_type time) const;

    virtual void Print(Option_t* option = "") const;

  private:
    TClonesArray fArrTimeRanges;

    ClassDef(AliTimeRangeMasking, 1);
};

#endif
