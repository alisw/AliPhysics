// --
// --
// Implementation for TTree output in PHOS DA
// for calibrating energy by pi0 and MIP.
// --
// -- Author: Hisayuki Torii (Hiroshima Univ.)
// --
#include <iostream>
#include <math.h>
#include <Rtypes.h>
#include "AliPHOSDATreeDigit.h"

ClassImp(AliPHOSDATreeDigit)
//------------------------------------------------------------------------
bool AliPHOSDATreeDigit::IsNeighbor(const AliPHOSDATreeDigit& digit) const{
  // Check wether the given digit is neighboring to this digit.
  // Return true if yes.
  if( fabs(this->GetRow() - digit.GetRow()) < 2 && fabs(this->GetCol() - digit.GetCol()) < 2 ){
    return true;
  }
  return false;
}
//------------------------------------------------------------------------
void AliPHOSDATreeDigit::Print(char* opt){
  // Print out
  std::cout<<" AliPHOSDATreeDigit:: "<<opt<<" E="<<fEnergy<<" (row,col)=("
	   <<GetRow()<<","<<GetCol()<<") absid="<<fAbsId<<std::endl;
}
//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const AliPHOSDATreeDigit& digit){
  // Print out
  //std::cout<<" AliPHOSDATreeDigit:: E="<<digit.GetEnergy()<<" (row,col)=("
  //<<digit.GetRow()<<","<<digit.GetCol()<<") absid="<<digit.GetAbsId()<<std::endl;
  out<<" AliPHOSDATreeDigit:: E="<<digit.fEnergy<<" (row,col)=("
  <<(int)(digit.fAbsId/56)<<","<<digit.fAbsId%56<<") absid="<<digit.fAbsId;
  return out;
}
//------------------------------------------------------------------------
