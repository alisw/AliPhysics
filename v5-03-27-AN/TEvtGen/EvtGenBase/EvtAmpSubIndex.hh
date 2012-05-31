//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      Caltech
//
// Module: EvtAmpSubIndex.hh
//
// Description:This class keeps track of indices on amplitude objects.
//             Used for a subset of indices in an EvtAmpIndex object.
//
// Modification history:
//
//    Ryd     Nov 22, 2002         Module created
//
//------------------------------------------------------------------------

#ifndef EVTAMPSUBINDEX_HH
#define EVTAMPSUBINDEX_HH

#include <vector>
class EvtAmpIndex;

class EvtAmpSubIndex {

public:

  EvtAmpSubIndex(EvtAmpIndex* ind,std::vector<int> sub);
  virtual ~EvtAmpSubIndex() {}

  int index();

private:

  EvtAmpIndex* _ind;
  std::vector<int> _sub;
  int _size;
  std::vector<int> _nstate;

};


#endif

