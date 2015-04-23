#include "ACORDE/MakeACORDEZeroMisAlignment.C"
#include "ACORDE/macros/MakeACORDEOCDBCalib.C"
void MakeACORDECDBObjects()
{
  MakeACORDEZeroMisAlignment();  //ACORDE/Align/Data               
  MakeACORDEOCDBCalib();  //ACORDE/Calib/Data               
}
