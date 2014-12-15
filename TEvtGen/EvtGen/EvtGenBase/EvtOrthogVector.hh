//--------------------------------------------------------------------------
// 
// Environment: 
// This software is part of the EvtGen package developed jointly 
// for the BaBar and CLEO collaborations.  If you use all or part 
// of it, please give an appropriate acknowledgement.
// 
// Copyright Information: See EvtGen/COPYRIGHT 
// Copyright (C) 2000 Caltech, LLNL
// 
// Module: EvtGen/EvtOrthogVector.hh
// 
// Description:
// 
// Modification history: 
//
// Lange August 11, 2000 Created
//
//------------------------------------------------------------------------


#ifndef EVTORTHOGVECTOR_HH
#define EVTORTHOGVECTOR_HH

#include <vector>

class EvtOrthogVector {

public:

  EvtOrthogVector(int n,
		  std::vector<double> *vectors);

  ~EvtOrthogVector();

  std::vector<double> getOrthogVector() {return _orthogVector;}


private:

  int _dimen;
  std::vector<int> _holder;
  void findOrthog(int dim, std::vector<int> invect,
		  std::vector<double> *vectors);

  std::vector<double> _orthogVector;

  int findEvenOddSwaps();


};

#endif

