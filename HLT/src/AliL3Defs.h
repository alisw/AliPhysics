#ifndef _ALIL3DEFS_H_
#define _ALIL3DEFS_H_

#include "AliL3RootTypes.h"

const Int_t NPatches = 6;
const Int_t NRowsSlice = 176;
const Int_t NRows[6][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
const Int_t NumRows[6] = {32,32,28,28,24,32};
const Double_t Pi = 3.14159265358979323846;
const Double_t ToRad = Pi/180.;
const Int_t MaxNPads = 256;
const Int_t MaxNTimeBins = 512;
const Double_t BField = 0.2;

#endif /* _ALIL3DEFS_H_ */
