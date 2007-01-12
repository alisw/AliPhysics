#ifndef FLIPPCB_H
#define FLIPPCB_H

class TString;
class AliMpPCB;
class AliMpMotifPosition;
class AliMpMotifType;
class AliMpSlatMotifMap;

const char* NameIt(const TString& baseName);

AliMpPCB* Duplicate(const AliMpPCB& src, AliMpSlatMotifMap& motifMap);

void flipPCB(const char* srcName);

#endif