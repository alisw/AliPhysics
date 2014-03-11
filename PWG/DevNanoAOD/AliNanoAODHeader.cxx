#include "AliNanoAODHeader.h"


ClassImp(AliNanoAODHeader)

AliNanoAODHeader& AliNanoAODHeader::operator=(const AliNanoAODHeader& evt) {

    AliVHeader::operator=(evt); // FIXME: ok?
    AliNanoAODStorage::operator=(evt);

    return *this;
}


