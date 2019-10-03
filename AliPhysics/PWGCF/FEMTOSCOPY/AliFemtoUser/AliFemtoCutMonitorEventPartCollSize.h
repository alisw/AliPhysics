////////////////////////////////////////////////////////////////////////////////
//                                                                            //  
// AliFemtoCutMonitorEventPartCollSize - the cut monitor for events to study  //
// the particle multiplicities per event                                      //
//                                                                            //
// author: Jesse Buxton jesse.thomas.buxton@cern.ch                           //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCUTMONITOREVENTPARTCOLLSIZE_H
#define ALIFEMTOCUTMONITOREVENTPARTCOLLSIZE_H

class AliFemtoEvent;
class AliFemtoTrack;
class AliFemtoV0;
class AliFemtoKink;
class AliFemtoPair;
class TH1D;
class TH2D;
class TList;

#include "AliFemtoString.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoCutMonitor.h"

class AliFemtoCutMonitorEventPartCollSize: public AliFemtoCutMonitor
{
public:

  AliFemtoCutMonitorEventPartCollSize();
  AliFemtoCutMonitorEventPartCollSize(const char *aName1, const int& aNBins1, const int& aNLow1, const int& aNHigh1, const char *aName2, const int& aNBins2, const int& aNLow2, const int& aNHigh2);
  virtual ~AliFemtoCutMonitorEventPartCollSize();

  AliFemtoCutMonitorEventPartCollSize(const AliFemtoCutMonitorEventPartCollSize &aCut);
  AliFemtoCutMonitorEventPartCollSize& operator=(const AliFemtoCutMonitorEventPartCollSize& aCut);

  virtual AliFemtoString Report();

  virtual void Fill(const AliFemtoParticleCollection *aCollection1, const AliFemtoParticleCollection *aCollection2);

  void Write();

  virtual TList *GetOutputList();

protected:
  TH1D *fCollSizePerEvent1;
  TH1D* fCollSizePerEvent2;

  TH2D* fCollSizePerEvent1vsEvent2;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCutMonitorEventPartCollSize, 1);
  /// \endcond
#endif

};

#endif
