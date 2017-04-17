// -*- C++ -*-
#ifndef ALI_SL_EVENT_HEADER_H
#define ALI_SL_EVENT_HEADER_H

#include <THashList.h>
#include <TParameter.h>

#include "AliLog.h"
#include "AliGenEventHeader.h"

class AliSLEventHeader : public AliGenEventHeader {
public:
  AliSLEventHeader(const char* name="STARLIGHT Event Header")
    : AliGenEventHeader(name)
    , fEventInfo(new THashList) {
    fEventInfo->SetOwner(kTRUE);
  }
  virtual ~AliSLEventHeader() { delete fEventInfo; fEventInfo = NULL; }

  Double_t GetEventInfo(const char* name) const;
  TList* GetEventInfo() { return fEventInfo; }
  void ClearEventInfo() { if (fEventInfo) fEventInfo->Clear(); }

  virtual void Print(Option_t* ) const;

protected:
private:
  AliSLEventHeader(const AliSLEventHeader& );
  AliSLEventHeader& operator=(const AliSLEventHeader& );

  THashList *fEventInfo;

  ClassDef(AliSLEventHeader, 1);
} ;

#endif // ALI_SL_EVENT_HEADER_H
