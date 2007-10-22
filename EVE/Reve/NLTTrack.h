// $Header$

#ifndef REVE_NLTTrack_H
#define REVE_NLTTrack_H

#include <Reve/Track.h>
#include <Reve/NLTBases.h>

namespace Reve {

class NLTProjection;

class NLTTrack : public Track,
		 public NLTProjected
{
  friend class NLTTrackGL;

private:
  NLTTrack(const NLTTrack&);            // Not implemented
  NLTTrack& operator=(const NLTTrack&); // Not implemented

  Vector*            fOrigPnts;
  Int_t              GetBreakPointIdx(Int_t start);
  void               GetBreakPoint(Int_t N, Bool_t back, Float_t& x, Float_t& y, Float_t& z);

protected:
  std::vector<Int_t> fBreakPoints;
  NLTProjection*     fProjection;

public:
  NLTTrack();
  virtual ~NLTTrack(){}

  virtual void SetProjection(NLTProjector* proj, NLTProjectable* model);

  virtual void UpdateProjection();
  virtual void MakeTrack(Bool_t recurse=kTRUE);

  void         PrintLineSegments();

  virtual void CtrlClicked(Reve::Track*); // marked as signal in Track

  ClassDef(NLTTrack, 1);
}; // endclass NLTTrack


/**************************************************************************/
// TrackRnrStyle
/**************************************************************************/

class NLTTrackList : public TrackList,
                     public NLTProjected
{
private:
  NLTTrackList(const NLTTrackList&);            // Not implemented
  NLTTrackList& operator=(const NLTTrackList&); // Not implemented

public:
  NLTTrackList();
  virtual ~NLTTrackList() {}

  virtual void SetProjection(NLTProjector* proj, NLTProjectable* model);
  virtual void UpdateProjection(){};

  ClassDef(NLTTrackList, 1); 
};// endclass NLTTrackList

}

#endif
