#ifndef _ALIL3DEFS_H_
#define _ALIL3DEFS_H_

#include <Rtypes.h>

const Int_t NRowsSlice = 175;
const Int_t NRows[6][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
//const Int_t NRows[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
const Double_t Pi = 3.14159265358979323846;
const Double_t ToRad = Pi/180.;
const Int_t MaxNPads = 256;
const Int_t MaxNTimeBins = 512;

/*fRow[0][0] = 0;     // first row
  fRow[0][1] = 31;
  fRow[1][0] = 32;
  fRow[1][1] = 63;
  fRow[2][0] = 64;
  fRow[2][1] = 91;
  fRow[3][0] = 92;
  fRow[3][1] = 119;
  fRow[4][0] = 120;
  fRow[4][1] = 143;   
  fRow[5][0] = 144;
  fRow[5][1] = 175;   // last row*/ 
#endif /* _ALIL3DEFS_H_ */
