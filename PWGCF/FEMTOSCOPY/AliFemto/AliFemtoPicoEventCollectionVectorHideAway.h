///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoPicoEventCollectionVectorHideAway: a helper class for         //
// managing many mixing buffers with up to three variables used for      //
// binning.                                                              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOPICOEVENTCOLLECTIONVECTORHIDEAWAY_H
#define ALIFEMTOPICOEVENTCOLLECTIONVECTORHIDEAWAY_H
#include "AliFemtoPicoEvent.h"
#include "AliFemtoPicoEventCollection.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include <vector>
#include <list>
#include <float.h>
#include <limits.h>

#if !defined(ST_NO_NAMESPACES)
using std::vector;
using std::list;
#endif

class AliFemtoPicoEventCollectionVectorHideAway {
public:
  AliFemtoPicoEventCollectionVectorHideAway(int bx=1, double lx=-FLT_MAX, double ux=FLT_MAX,
					    int by=1, double ly=-FLT_MAX, double uy=FLT_MAX,
					    int bz=1, double lz=-FLT_MAX, double uz=FLT_MAX);
  AliFemtoPicoEventCollectionVectorHideAway(const AliFemtoPicoEventCollectionVectorHideAway& aColl);
  ~AliFemtoPicoEventCollectionVectorHideAway();
  AliFemtoPicoEventCollectionVectorHideAway& operator=(const AliFemtoPicoEventCollectionVectorHideAway& aColl);

  AliFemtoPicoEventCollection* PicoEventCollection(int bx, int by, int bz);
  AliFemtoPicoEventCollection* PicoEventCollection(double x, double y=0, double z=0);
  unsigned int GetBinXNumber(double x);
  unsigned int GetBinYNumber(double y);
  unsigned int GetBinZNumber(double z);
private:
  int fBinsTot;                                        // Total number of bins 
  int fBinsx,fBinsy,fBinsz;                            // Number of bins on x, y, z axis
  double fMinx,fMiny,fMinz;                            // Minima on x, y, z axis
  double fMaxx,fMaxy,fMaxz;                            // Maxima on x, y, z axis
  double fStepx,fStepy,fStepz;                         // Steps on x, y, z axis
  AliFemtoPicoEventCollection* fCollection;            // Pico event collection
  AliFemtoPicoEventCollectionVector fCollectionVector; // Collection vector
};

#endif
