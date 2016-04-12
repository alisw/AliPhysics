#ifndef ALIHLTCDHWRAPPER_H
#define ALIHLTCDHWRAPPER_H

#include <assert.h>
#include "AliRawDataHeader.h"
#include "AliRawDataHeaderV3.h"

#define CHECK_AND_CALL(func, args...)					\
  ( GetVersion() == 2 ?							\
    reinterpret_cast<const AliRawDataHeader*>(fCDH)->func(args) :	\
    reinterpret_cast<const AliRawDataHeaderV3*>(fCDH)->func(args) )


class AliHLTCDHWrapper {
 public:
  AliHLTCDHWrapper() : fCDH(NULL) {}
  AliHLTCDHWrapper(const AliHLTCDHWrapper& other) : fCDH(other.fCDH) { CheckVersion(); }
  AliHLTCDHWrapper(const void* cdh) : fCDH(cdh) { CheckVersion(); }

  ~AliHLTCDHWrapper() {}

  inline AliHLTCDHWrapper& operator=(const AliHLTCDHWrapper& other) { 
    fCDH = other.fCDH;
    CheckVersion(); 
    return *this;
  }

  inline AliHLTCDHWrapper& operator=(const void*& cdh) { 
    fCDH = cdh;
    CheckVersion(); 
    return *this;
  }

  inline void CheckVersion() {
#ifdef DEBUG
    if(fCDH)
      assert(GetVersion() == 2 || GetVersion() == 3);
#endif
  }
  
  inline UChar_t GetVersion() const { 
    return (reinterpret_cast<const AliRawDataHeader*>(fCDH))->GetVersion();
  }

  inline UInt_t GetHeaderSize() {
    return (GetVersion() == 2 ? 
	    sizeof(AliRawDataHeader) : sizeof(AliRawDataHeaderV3) );
  }

  inline const void* GetHeader() const {
    return fCDH;
  }

  inline UInt_t GetDataSize() const {
    //first word, independent of Version
    return *((UInt_t*)fCDH);
  }

  inline UShort_t GetEventID1() const {
    return CHECK_AND_CALL(GetEventID1);
  }

  inline UInt_t GetEventID2() const {
    return CHECK_AND_CALL(GetEventID2);
  }

  inline UChar_t GetL1TriggerMessage() const {
    return CHECK_AND_CALL(GetL1TriggerMessage);
  }

  inline UChar_t GetAttributes() const {
    return CHECK_AND_CALL(GetAttributes);
  }

  inline Bool_t TestAttribute(Int_t index) const {
    return CHECK_AND_CALL(TestAttribute, index);
  }

  /*
  inline void SetAttribute(Int_t index) {
    CHECK_AND_CALL(SetAttribute, index);
  }
  */

  /*
  inline void ResetAttribute(Int_t index) {
    CHECK_AND_CALL(ResetAttribute, index);
  }
  */

  inline UInt_t GetSubDetectors() const {
    return CHECK_AND_CALL(GetSubDetectors);
  }

  inline UInt_t GetStatus() const {
    return CHECK_AND_CALL(GetStatus);
  }

  inline UInt_t GetMiniEventID() const {
    return CHECK_AND_CALL(GetMiniEventID);
  }

  inline ULong64_t GetTriggerClasses() const {
    return CHECK_AND_CALL(GetTriggerClasses);
  }

  inline ULong64_t GetTriggerClassesNext50() const {
    return CHECK_AND_CALL(GetTriggerClassesNext50);
  }

  inline ULong64_t GetROI() const {
    return CHECK_AND_CALL(GetROI);
  }

  /*
  inline void SetTriggerClass(ULong64_t mask) {
    CHECK_AND_CALL(SetTriggerClass, mask);
  }
  */

 private:
  const void* fCDH;

};


#endif
