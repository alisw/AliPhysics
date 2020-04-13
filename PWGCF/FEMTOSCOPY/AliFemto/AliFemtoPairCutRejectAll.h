///
/// \file AliFemtoPairCutRejectAll.h
///


#pragma once

#ifndef ALIFEMTOPAIRCUTREJECTALL_H
#define ALIFEMTOPAIRCUTREJECTALL_H

#include "AliFemtoPairCut.h"


/// \class AliFemtoPairCutRejectAll
/// \brief Reject all pairs - contains no state
///
/// This is useful for "analyses" only using track monitors, and
/// can ignore the results of pairs.
///
/// \author Andrew Kubera, Ohio State Univeristy
///
class AliFemtoPairCutRejectAll : public AliFemtoPairCut {
public:

  AliFemtoPairCutRejectAll()
    : AliFemtoPairCut()
    { }

  AliFemtoPairCutRejectAll(const AliFemtoPairCutRejectAll &orig)
    : AliFemtoPairCut(orig)
    { }

  virtual ~AliFemtoPairCutRejectAll()
    { }

  AliFemtoPairCutRejectAll& operator=(const AliFemtoPairCutRejectAll &rhs)
    {
      AliFemtoPairCut::operator=(rhs);
      return *this;
    }

  bool Pass(const AliFemtoPair *) override
    {
      return false;
    }

  TList* ListSettings() override
    {
      return new TList();
    }

  AliFemtoString Report() override
    {
      AliFemtoString report("AliFemtoPairCutRejectAll report:");
      return report;
    }

  AliFemtoPairCutRejectAll* Clone() const override
    {
      return new AliFemtoPairCutRejectAll(*this);
    }
};


#endif
