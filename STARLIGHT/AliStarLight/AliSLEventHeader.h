// -*- C++ -*-
#ifndef ALI_SL_EVENT_HEADER_H
#define ALI_SL_EVENT_HEADER_H

#include <TMap.h>
#include <TVector2.h>

#include "AliLog.h"
#include "AliGenEventHeader.h"

template<typename T>
struct type2val {
  typedef T type;
} ;

class AliSLEventHeader : public AliGenEventHeader {
public:
  AliSLEventHeader(const char* name="STARLIGHT Event Header")
    : AliGenEventHeader(name)
    , fEventInfo(new TMap) {
    fEventInfo->SetOwnerKeyValue(kTRUE);
  }
  virtual ~AliSLEventHeader() { delete fEventInfo; fEventInfo = NULL; }

  template<typename T>
  T GetEventInfo(const char* name) const {
    return GetEventInfo(name, type2val<T>());
  }
  TMap* GetEventInfo() { return fEventInfo; }

  virtual void Print(Option_t* ) const;

protected:
  Double_t GetEventInfo(const char* name, const type2val<Double_t>&) const;
  TVector2 GetEventInfo(const char* name, const type2val<TVector2>&) const;
private:
  AliSLEventHeader(const AliSLEventHeader& );
  AliSLEventHeader& operator=(const AliSLEventHeader& );

  TMap *fEventInfo;

  ClassDef(AliSLEventHeader, 2);
} ;

#endif // ALI_SL_EVENT_HEADER_H
