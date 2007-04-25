/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.3  2002/02/04 18:58:33  laue
 * *** empty log message ***
 *
 * Revision 1.2  2001/11/11 18:34:13  laue
 * AliFemtoPicoEventCollectionVectorHideAway: updated for 3d grid
 * AliFemtoVertexMultAnalysis: new
 *
 * Revision 1.1  2000/07/16 21:44:11  laue
 * Collection and analysis for vertex dependent event mixing
 *
 *
 **************************************************************************/
#ifndef AliFemtoPicoEventCollectionVectorHideAway_hh
#define AliFemtoPicoEventCollectionVectorHideAway_hh
#include "Infrastructure/AliFemtoPicoEvent.h"
#include "Infrastructure/AliFemtoPicoEventCollection.h"
#include "Infrastructure/AliFemtoPicoEventCollectionVector.h"
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
  AliFemtoPicoEventCollection* PicoEventCollection(int, int, int);
  AliFemtoPicoEventCollection* PicoEventCollection(double x, double y=0, double z=0);
private:
  int fBinsTot;
  int fBinsx,fBinsy,fBinsz;
  double fMinx,fMiny,fMinz;
  double fMaxx,fMaxy,fMaxz;
  double fStepx,fStepy,fStepz;
  AliFemtoPicoEventCollection* fCollection;
  AliFemtoPicoEventCollectionVector fCollectionVector;
};

#endif
