#include "f_Decay.h"

namespace Tauolapp
{

void TauolaDecay(int sign_type, double *polx, double* poly, 
                 double * polz, double * poln){

  taupos_.npa=1; //tau position in particle list
  taupos_.npb=1;

  double pol[4]={0};
  dekay_(&sign_type, pol);
  *polx=pol[0];
  *poly=pol[1];
  *polz=pol[2];
  *poln=pol[3];

  //  std::cout << "Polarimetric: "<<pol[0]<<","<<pol[1]
  //            <<","<<pol[2]<<","<<pol[3]<<std::endl;
  

}

void TauolaWriteDecayToEventRecord(int sign_type){
  taupos_.npa=1; //tau position in particle list
  taupos_.npb=1;

  double pol[4]={0};

  sign_type+=10;
  dekay_(&sign_type, pol); //write to event record

}

} // namespace Tauolapp
