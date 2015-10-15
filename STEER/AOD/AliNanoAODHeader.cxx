#include "AliNanoAODHeader.h"
#include "AliLog.h"

ClassImp(AliNanoAODHeader)

AliNanoAODHeader& AliNanoAODHeader::operator=(const AliNanoAODHeader& evt) {

    AliVHeader::operator=(evt); // FIXME: ok?
    AliNanoAODStorage::operator=(evt);

    return *this;
}


void  AliNanoAODHeader::Clear(Option_t * /*opt*/) {
  // empty storage
  fVars.clear();
  fNVars = 0;
}

void AliNanoAODHeader::NotImplemented(void) const {
  AliError("Not implemented");
}
