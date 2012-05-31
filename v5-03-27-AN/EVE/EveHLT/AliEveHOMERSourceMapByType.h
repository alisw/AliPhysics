// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveAliEveHOMERSourceMapByType_H
#define AliEveAliEveHOMERSourceMapByType_H

#include "AliEveHOMERSourceMap.h"

//______________________________________________________________________________
//
// AliEveHOMERSourceMap is an abstract container for HLT HOMER sources.
// The concrete implementations AliEveHOMERSourceMapByDet and
// AliEveHOMERSourceMapByType allow retrieval of HOMER sources in proper
// order as required for their display in EVE object browser.

class AliEveHOMERSourceMapByType : public AliEveHOMERSourceMap
{
public:
  AliEveHOMERSourceMapByType(ESourceGrouping_e grouping);
  virtual ~AliEveHOMERSourceMapByType() {}

  virtual void FillMap(const TList* handles, Bool_t def_state);

protected:
  typedef std::map<AliEveHOMERSource::SourceId,
		   AliEveHOMERSource::SourceState,
		   AliEveHOMERSource::SourceId::CmpByType> Map_t;
  typedef Map_t::iterator                                  Map_i;

  Map_t fMap;

  struct iterator_imp : public iterator_imp_base,
			public Map_i
  {
    iterator_imp()         : Map_i()   {}
    iterator_imp(Map_i mi) : Map_i(mi) {}
    virtual ~iterator_imp() {}

    virtual const AliEveHOMERSource::SourceId&    id()    const { return Map_i::operator*().first; }
    virtual const AliEveHOMERSource::SourceState& state() const { return Map_i::operator*().second; }
    virtual       AliEveHOMERSource::SourceState& state() { return Map_i::operator*().second; }

    virtual iterator_imp* clone() const
    { return new iterator_imp(*this); }

    virtual iterator_imp& operator++()
    { Map_i::operator++(); return *this; }

    virtual bool operator!=(const iterator_imp_base& rhs) const
    { const Map_i &lhs = *this; return lhs != dynamic_cast<const Map_i&>(rhs); }

    virtual TString description() const;
  };

  void insert(AliEveHOMERSource::SourceId& sid, AliEveHOMERSource::SourceState& sst, Bool_t def_state);

  virtual iterator_imp_base* iterator_imp_new()   { return new iterator_imp; }
  virtual iterator_imp_base* iterator_imp_begin() { return new iterator_imp(fMap.begin()); }
  virtual iterator_imp_base* iterator_imp_end()   { return new iterator_imp(fMap.end());   }

  // ClassDef(AliEveHOMERSourceMapByType, 0);
};

#endif
