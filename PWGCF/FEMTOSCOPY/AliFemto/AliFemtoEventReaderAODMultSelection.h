///
/// \file AliFemtoEventReaderAODMultSelection.h
///

#pragma once

#ifndef ALIFEMTOEVENTREADER_MULTSELECTION_H_
#define ALIFEMTOEVENTREADER_MULTSELECTION_H_

#include "AliFemtoEventReaderAODChain.h"

class AliFemtoEvent;

/// \class AliFemtoEventReaderAODMultSelection
/// \brief UNSTABLE! Event reader designed to read AOD files using the MultiplicitySelection
///        task, required for AOD files taken during run2 (2015)
///
/// \author Andrew Kubera, The Ohio State University <andrew.kubera@cern.ch>
///
class AliFemtoEventReaderAODMultSelection : public AliFemtoEventReaderAODChain {
public:
  /// Default Constructor
  AliFemtoEventReaderAODMultSelection()
    : AliFemtoEventReaderAODChain()
  {}

  AliAODEvent* GetAodSource()
    { return fEvent; }

protected:
  /// Returns a new AliFemtoEvent from the AliAODEvent pointed to by fEvent data member.
  virtual AliFemtoEvent* CopyAODtoFemtoEvent();

};


#endif
