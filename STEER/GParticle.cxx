///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Particle class                                                           //
//  This class contains the information on the particles being transported   //
//  by the MonteCarlo. For the moment the numbering convention is the one    //
//  of GEANT 3.21                                                            //
//                                                                           //
//Begin_Html
/*
<img src="picts/GParticleClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GParticle.h"
 
ClassImp(GParticle)
 
//_____________________________________________________________________________
const Text_t *GParticle::GetName() const
{
  //
  // Return name of this particle (Geant convention)
  //
  const char *gnames[50] = {
    "Gamma","Positron","Electron","Neutrino","Muon +"
    ,"Muon -","Pion 0","Pion +","Pion -","Kaon 0 long"
    ,"Kaon +","Kaon -","Neutron","Proton","AntiProton"
    ,"Kaon 0 short","Eta","Lambda","Sigma +","Sigma 0"
    ,"Sigma -","Xi 0","Xi -","Omega -","Antineutron"
    ,"AntiLambda","AntiSigma -","AntiSigma 0","AntiSigma +","AntiXi 0"
    ,"AntiXi +","AntiOmega +","Tau +","Tau -","D +"
    ,"D -","D 0","Anti D 0","DS+","DS-"
    ,"Lambda C +","W +","W -","Z 0","Deuteron"
    ,"Triton","Alpha","Geantino","He3","Cerenkov"};
  
  if (fKF < 1 || fKF >50) return "GParticle";
  return gnames[fKF-1];
}

//_____________________________________________________________________________
const Text_t *GParticle::GetTitle() const
{
  //
  // Return title of this particle (Geant convention)
  //
  static char name[20];
  sprintf(name,"GParticle %d",fKF);
  return name;
}
