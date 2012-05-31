// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveAliEVEHOMERSource_H
#define AliEveAliEVEHOMERSource_H

#include <TEveElement.h>

#include <TNamed.h>

class AliHLTHOMERSourceDesc;

class AliEveHOMERSource : public TEveElement,
			  public TNamed
{
public:
  struct SourceId
  {
    TString fDet, fSDet, fSSDet, fType;

    SourceId(): fDet(), fSDet(), fSSDet(), fType() {}

    struct CmpByDet
    {
      bool operator()(const SourceId& s1, const SourceId& s2) const
      {
	Int_t r;
	if ((r = s1.fDet  .CompareTo(s2.fDet)  )) return r < 0;
	if ((r = s1.fSDet .CompareTo(s2.fSDet) )) return r < 0;
	if ((r = s1.fSSDet.CompareTo(s2.fSSDet))) return r < 0;
	if ((r = s1.fType .CompareTo(s2.fType) )) return r < 0;
	return false;
      }
    };

    struct CmpByType
    {
      bool operator()(const SourceId& s1, const SourceId& s2) const
      {
	Int_t r;
	if ((r = s1.fType .CompareTo(s2.fType) )) return r < 0;
	if ((r = s1.fDet  .CompareTo(s2.fDet)  )) return r < 0;
	if ((r = s1.fSDet .CompareTo(s2.fSDet) )) return r < 0;
	if ((r = s1.fSSDet.CompareTo(s2.fSSDet))) return r < 0;
	return false;
      }
    };

  };

  struct SourceState
  {
    Bool_t                  fState;
    AliHLTHOMERSourceDesc  *fHandle;

    SourceState() : fState(kFALSE), fHandle(0) {}
    SourceState(Bool_t state) : fState(state), fHandle(0) {}
  };


  AliEveHOMERSource(const Text_t* n="HOMER Source", const Text_t* t="");
  virtual ~AliEveHOMERSource() {}

  const SourceId* GetSourceId() const  { return fSrcId; }
  void SetSourceId(const SourceId* id) { fSrcId = id; }

  SourceState* GetSourceState() const  { return fSrcState; }
  void SetSourceState(SourceState* st) { fSrcState = st; TEveElement::SetRnrState(st->fState); }

  void SetSource(const SourceId* id, SourceState* st) { fSrcId = id; fSrcState = st; TEveElement::SetRnrState(st->fState); }

  virtual Bool_t SingleRnrState() const { return kTRUE; }
  virtual Bool_t SetRnrState(Bool_t rnr);

protected:
  const SourceId    *fSrcId;
        SourceState *fSrcState;

private:
  AliEveHOMERSource(const AliEveHOMERSource&);            // Not implemented
  AliEveHOMERSource& operator=(const AliEveHOMERSource&); // Not implemented

  ClassDef(AliEveHOMERSource, 0); // Description of an HOMER source.
};

#endif
