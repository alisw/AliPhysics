///
/// \file AliFemtoEventReaderAODKinematicsMultSelection.h
///

#pragma once

#ifndef ALIFEMTOEVENTREADER_KINEMATICSMULTSELECTION_H_
#define ALIFEMTOEVENTREADER_KINEMATICSMULTSELECTION_H_

#include "AliFemtoEventReaderAODKinematicsChain.h"

class AliFemtoEvent;

/// \class AliFemtoEventReaderAODMultSelection
/// \brief UNSTABLE! Event reader designed to read AOD files using the MultiplicitySelection
///        task, required for AOD files taken during run2 (2015)
///
/// \author Andrew Kubera, The Ohio State University <andrew.kubera@cern.ch>
///
class AliFemtoEventReaderAODKinematicsMultSelection : public AliFemtoEventReaderAODKinematicsChain {
public:
  /// Default Constructor
  AliFemtoEventReaderAODKinematicsMultSelection()
    : AliFemtoEventReaderAODKinematicsChain()
  {}

  AliAODEvent* GetAodSource()
    { return fEvent; }

protected:
  /// Returns a new AliFemtoEvent from the AliAODEvent pointed to by fEvent data member.
  virtual AliFemtoEvent* ReturnHbtEvent();

};


#endif
