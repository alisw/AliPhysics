//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      Caltech, UCSB
//
// Module: EvtGen/EvtAmpIndex.hh
//
// Description:This class keeps track of indices on amplitude objects.
//
// Modification history:
//
//    Ryd     Nov 22, 2002         Module created
//
//------------------------------------------------------------------------

#ifndef EVTAMPINDEX_HH
#define EVTAMPINDEX_HH

#include <vector>

class EvtAmpIndex {

  friend class EvtAmpSubIndex;

public:

  EvtAmpIndex(std::vector<int> ind);
  virtual ~EvtAmpIndex() {}

  void reset();
  bool next();

  int index();

private:

  std::vector<int> _ind;
  int _size;
  std::vector<int> _state;
  std::vector<int> _nstate;

};


#endif

