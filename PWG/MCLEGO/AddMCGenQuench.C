// Author: Leticia Cunqueiro

//Some explanations: 
//if fpythia==1, quenched pythia (qpythia) is configured.
//if fpythia==2, pyquen afterburner applyed
  
// The list of avaliable tunes can be checked in  $ALICE_ROOT/PYTHIA6/pythiaX.f , inside subroutine PYTUNE 
// Note that if you are using QPYTHIA , it only works with Q2-ordered showers, that is tune<300.
//More modern tunes are pt-ordered and use a different FSR shower. 

//In the implementation of qpythia in aliroot I set as default the following computation of the geometry:
//   -Glauber to determine the overlapping region of the nuclei and the impact parameter of the collision
//  -Random sampling of the initial coordinates of the hard scattering in the overlapping region
//  -Given the coordinates and direction of each parton in the shower, calculate the path length to the "end"
//   of the medium and the integrated qhat along this path length, alla PQM.

//In the PQM approach, you integrate the qhat(dx,dy) along the path length and this is purely geometrical. See formula (13) in  http://arxiv.org/pdf/hep-ph/0406201.pdf
//There is a free parameter k (in fm), that sets the scale of the transport coefficient.  In 0-10% PbPb collisions,
//the average <qhat>  and k are related via a number:
//  <qhat>/k =5.87 e^-5 fm^-4
//If you want a <qhat> of 10 GeV2/fm,
//         10 /5.87e^-5 (GeV2 fm^3)=k
//If you want k in fm then you have to divide by the squared of hbarc (hbar c=0.197 GeV fm)
// This gives k=4.4e^6 fm, which is the quench value we set in SetQhat

R__LOAD_LIBRARY(liblhapdf)
R__LOAD_LIBRARY(libqpythia)
R__LOAD_LIBRARY(libEGPythia6)
R__LOAD_LIBRARY(libAliPythia6)

#include "AliGenerator.h"
#include "AliGenPythia.h"

AliGenerator* CreatePythia6Gen(Float_t e_cms, Int_t ptHardMin, Int_t ptHardMax, Int_t fpythia, Double_t quench = 4.4e6, Int_t ianglepyquen=2,Float_t ptWeight=0) {
    
  AliGenPythia* genP = new AliGenPythia(1);

  //   vertex position and smearing 
  genP->SetVertexSmear(kPerEvent);
  
  //   charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb
  if (ptHardMin>0.) {
    genP->SetProcess(kPyJets);
    genP->SetPtHard((float)ptHardMin,(float)ptHardMax);
    genP->SetWeightPower(ptWeight);
   } else
    genP->SetProcess(kPyMb); // Minimum Bias  

  //   Centre of mass energy 
  genP->SetEnergyCMS(e_cms); // in GeV
    
  //for jet quenching with QPYTHIA
  if (fpythia == 1){
    genP->SetTune(103); //tune DW, standard choice for Q2 showers
    genP->SetQuench(4);
    genP->SetQhat(quench); 
  }

  //for pyquen afterburner
  if (fpythia == 2){
    genP->SetTune(103); //tune DW, standard choice for Q2 showers
    genP->SetQuench(2);
    genP->SetPyquenPar(1,0.1,0,0,ianglepyquen);
  }

  return genP;
}

AliGenerator* AddMCGenQuench(Float_t e_cms = 2760., Double_t ptHardMin = 0., Double_t ptHardMax = 0., Int_t fpythia = 1, Double_t quench=4.4e6, Int_t ianglepyquen = 2,Float_t ptWeight=0) 
{
  //Add Pythia generator: pt-hard bin or min bias
 
  return CreatePythia6Gen(e_cms, ptHardMin, ptHardMax, fpythia, quench, ianglepyquen,ptWeight);
}
