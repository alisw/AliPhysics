// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveAliEveHOMERSourceMap_H
#define AliEveAliEveHOMERSourceMap_H

#include "AliEveHOMERSource.h"

#include <TNamed.h>

#include <map>

class AliEveHOMERSourceMap : public TNamed
{
protected:
  struct iterator_imp_base
  {
    virtual ~iterator_imp_base() {}

    virtual const AliEveHOMERSource::SourceId&    id()    = 0;
    virtual       AliEveHOMERSource::SourceState& state() = 0;

    virtual iterator_imp_base* clone() = 0;

    virtual iterator_imp_base& operator++() = 0;
    virtual bool operator!=(const iterator_imp_base& o) = 0;

    virtual TString description() = 0;
  };

public:
  struct iterator
  {
    iterator_imp_base* m_imp;

    iterator()                       :  m_imp(0)   {}
    iterator(iterator_imp_base* imp) :  m_imp(imp) {}
    ~iterator() { delete m_imp; }

    iterator& operator= (const iterator& o) { delete m_imp; m_imp = o.m_imp->clone(); return *this; }
    bool      operator!=(const iterator& o) { return m_imp->operator!=(*o.m_imp); }

    const AliEveHOMERSource::SourceId&    id()    { return m_imp->id(); }
          AliEveHOMERSource::SourceState& state() { return m_imp->state(); }

    Int_t   level();
    TString description() { return m_imp->description(); }

    iterator& operator++() { m_imp->operator++(); return *this; }
    iterator& skip_children()
    {
      Int_t lvl = level();
      do operator++(); while (level() > lvl);
      return *this;
    }
  };

  iterator begin() { return iterator(iterator_imp_begin()); }
  iterator end()   { return iterator(iterator_imp_end()); }
  
  enum ESourceGrouping_e { kSG_ByDet, kSG_ByType };

protected:
  ESourceGrouping_e fGrouping; // Not used so far ...

  virtual iterator_imp_base* iterator_imp_new()   = 0; // Not used so far ...
  virtual iterator_imp_base* iterator_imp_begin() = 0;
  virtual iterator_imp_base* iterator_imp_end()   = 0;

public:
  AliEveHOMERSourceMap(ESourceGrouping_e grouping);
  virtual ~AliEveHOMERSourceMap() {}

  static AliEveHOMERSourceMap* Create(ESourceGrouping_e grouping);

  virtual void FillMap(TList* handles, Bool_t def_state) = 0;

  void PrintXXX();

  ClassDef(AliEveHOMERSourceMap, 0);
};

/******************************************************************************/

class AliEveHOMERSourceMapByDet : public AliEveHOMERSourceMap
{
protected:
  typedef std::map<AliEveHOMERSource::SourceId,
		   AliEveHOMERSource::SourceState,
		   AliEveHOMERSource::SourceId::CmpByDet>  Map_t;
  typedef Map_t::iterator                                  Map_i;

  Map_t fMap;

  struct iterator_imp : public iterator_imp_base,
			public Map_i
  {
    iterator_imp()         : Map_i()   {}
    iterator_imp(Map_i mi) : Map_i(mi) {}
    virtual ~iterator_imp() {}

    virtual const AliEveHOMERSource::SourceId&    id()    { return Map_i::operator*().first; }
    virtual       AliEveHOMERSource::SourceState& state() { return Map_i::operator*().second; }

    virtual iterator_imp* clone()
    { return new iterator_imp(*this); }

    virtual iterator_imp& operator++()
    { Map_i::operator++(); return *this; }

    virtual bool operator!=(const iterator_imp_base& o)
    { return Map_i::operator!=(dynamic_cast<const Map_i&>(o)); }

    virtual TString description();
  };

  void insert(AliEveHOMERSource::SourceId& sid, AliEveHOMERSource::SourceState& sst, Bool_t def_state);

  virtual iterator_imp_base* iterator_imp_new()   { return new iterator_imp; }
  virtual iterator_imp_base* iterator_imp_begin() { return new iterator_imp(fMap.begin()); }
  virtual iterator_imp_base* iterator_imp_end()   { return new iterator_imp(fMap.end());   }

public:
  AliEveHOMERSourceMapByDet(ESourceGrouping_e grouping);
  virtual ~AliEveHOMERSourceMapByDet() {}

  virtual void FillMap(TList* handles, Bool_t def_state);

  // ClassDef(AliEveHOMERSourceMapByDet, 1);
};

/******************************************************************************/

class AliEveHOMERSourceMapByType : public AliEveHOMERSourceMap
{
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

    virtual const AliEveHOMERSource::SourceId&    id()    { return Map_i::operator*().first; }
    virtual       AliEveHOMERSource::SourceState& state() { return Map_i::operator*().second; }

    virtual iterator_imp* clone()
    { return new iterator_imp(*this); }

    virtual iterator_imp& operator++()
    { Map_i::operator++(); return *this; }

    virtual bool operator!=(const iterator_imp_base& o)
    { return Map_i::operator!=(dynamic_cast<const Map_i&>(o)); }

    virtual TString description();
  };

  void insert(AliEveHOMERSource::SourceId& sid, AliEveHOMERSource::SourceState& sst, Bool_t def_state);

  virtual iterator_imp_base* iterator_imp_new()   { return new iterator_imp; }
  virtual iterator_imp_base* iterator_imp_begin() { return new iterator_imp(fMap.begin()); }
  virtual iterator_imp_base* iterator_imp_end()   { return new iterator_imp(fMap.end());   }

public:
  AliEveHOMERSourceMapByType(ESourceGrouping_e grouping);
  virtual ~AliEveHOMERSourceMapByType() {}

  virtual void FillMap(TList* handles, Bool_t def_state);

  // ClassDef(AliEveHOMERSourceMapByType, 1);
};

#endif
