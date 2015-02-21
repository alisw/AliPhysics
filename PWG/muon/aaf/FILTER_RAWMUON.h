#ifndef FILTER_RAWMUON_H
#define FILTER_RAWMUON_H

class TTree;

namespace AAF {
  
  int FILTER_RAWMUON(const char* from, const char* to);

  namespace RAWMUON {

  	  int CheckFile(const char* from, int maxevents=-1);
  	  void DisableBranches(TTree* tree);

  }

}

#endif
