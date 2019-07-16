/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
#ifndef ALITIMERANGECUT_H
#define ALITIMERANGECUT_H

/// \file AliTimeRangeCut.h
/// \brief A class for cutting on AliTimeRangeMasking definitions
/// \author Jens Wiechula, jens.wiechula@ikf.uni-frankfurt.de

#include "TString.h"

#include "AliTimeRangeMasking.h"

class AliVEvent;

/// \class AliTimeRangeCut
/// \brief A class for cutting on AliTimeRangeMasking definitions
///
/// This class can be used to cut events which are masked by their global id
/// (AliVHeader::GetEventIdAsLong()) following a certain mask defined in AliTimeRangeMasking.
/// The masking regions are defined in the OADB object COMMON/PHYSICSSELECTION/data/TimeRangeMasking.root
/// They are set by the macro $ALICE_PHYSICS_SRC/OADB/macros/SetupTimeRangeMasking.C
///
/// The class should be used as follows:
/// * Use AliEventCuts with calling UseTimeRangeCut
/// * Use Standalone
///   - Add `AliTimeRangeCut fTimeRangeCut;` as data member in the task
///   - In UserExec call
///     `fTimeRangeCut.InitFromEvent(InputEvent());`
///     `const Bool_t cutThisEvent = fTimeRangeCut.CutEvent(InputEvent());`
///     or in case the time range has been masked for different reasons, but one is only interested in a specific reason
///     `const Bool_t cutThisEvent = fTimeRangeCut.CutEvent(InputEvent(), bitmask);`
///     for the bit definitions see [AliTimeRangeMask](@ref AliTimeRangeMask)
class AliTimeRangeCut : public TObject {
  public:
    AliTimeRangeCut() : fOADBPath(), fTimeRangeMasking(0x0), fLastRun(-1) {}
    ~AliTimeRangeCut() { delete fTimeRangeMasking; }

    void InitFromEvent(const AliVEvent* event); 
    void InitFromRunNumber(const Int_t run);

    UShort_t GetMask(const AliVEvent* event) const;
    UShort_t GetMask(const ULong64_t gid) const;

    Bool_t CutEvent(const AliVEvent* event, const UShort_t mask = 0) const;
    Bool_t CutEvent(const ULong64_t gid, const UShort_t mask = 0) const;

    void SetOADBPath(const TString& path) { fOADBPath = path; }
    
    const TString& GetOADPath() const { return fOADBPath; }

  private:
    AliTimeRangeCut(const AliTimeRangeCut&);
    AliTimeRangeCut& operator= (const AliTimeRangeCut&);

    TString fOADBPath; ///< OADB path
    AliTimeRangeMasking<ULong64_t, UShort_t>* fTimeRangeMasking; //!< Time Range masksking object
    Int_t fLastRun; //!< last set run number

    ClassDef(AliTimeRangeCut, 1)
};

#endif
