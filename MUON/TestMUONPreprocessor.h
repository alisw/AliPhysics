#ifndef TESTMUONPREPROCESSOR_H
#define TESTMUONPREPROCESSOR_H

class TMap;

TMap* CreateDCSAliasMap(const char* inputCDB, Int_t runNumber);

#endif
