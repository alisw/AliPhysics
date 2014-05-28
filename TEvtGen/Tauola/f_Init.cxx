#include "f_Init.h"
#include "Log.h"
#include <iostream>
using std::endl;

namespace Tauolapp
{

void f_interface_tauolaInitialize(int pdg_id, int firstDecayMode, 
                                  int secondDecayMode, bool rad,
                                  double rad_cut_off, double iniphy){

  double pol[4]={0,0,0,0}; //dummy

  jaki_.jak1=firstDecayMode;
  jaki_.jak2=secondDecayMode;

  taurad_.itdkrc=rad;
  taurad_.xk0dec=rad_cut_off; //<-this one was causing nan's
  idfc_.idff=pdg_id;
  
  inimas_();
  initdk_();
  float iniphy_param = iniphy;
  iniphy_(&iniphy_param);

  int init_state=-1;
  //  dexay_(&init_state,pol);
  dekay_(&init_state,pol);
}

double f_getTauMass(){
  return (double) parmas_.amtau;
}

void f_interface_tauolaInitialise(int pdg_id, int firstDecayMode, 
                                  int secondDecayMode, bool rad,
                                  double rad_cut_off, double iniphy)
{
  Log::Warning() <<"Deprecated routine 'f_interface_tauolaInitialise'"<<endl;
  Log::Warning(0)<<"Use 'f_interface_tauolaInitialize' instead."<<endl;

  f_interface_tauolaInitialize(pdg_id, firstDecayMode, 
                               secondDecayMode, rad,
                               rad_cut_off, iniphy);

  // Deprecated routines:  initialise, setInitialisePhy,
  //                       f_interface_tauolaInitialise
}

} // namespace Tauolapp
