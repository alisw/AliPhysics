#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

#include "TString.h"

#include "EMCAL/AliEMCALReconstructioner.h"
#endif

void Go(TString deb = ""){
  AliEMCALReconstructioner * a ;   

  cout << "AliEMCAL:> Single File default reconstruction started" << endl ;
  a = new AliEMCALReconstructioner("galice.root") ;  //first -single file default reconstruction
  a->ExecuteTask(deb.Data()) ;
  cout << "AliEMCAL:> Single File default reconstruction finished" << endl ;
  // delete a ; 

  cout << "AliEMCAL:> Single File branch TEST reconstruction started" << endl ;
  a = new AliEMCALReconstructioner("galice.root","test") ;  //another branch single file recontruction
  a->ExecuteTask(deb.Data()) ;
  cout << "AliEMCAL:> Single File branch TEST reconstruction ended" << endl ;
  //delete a ; 
  
  cout << "AliEMCAL:> Split File default reconstruction started" << endl ;
  a = new AliEMCALReconstructioner("galice.root","Default",kTRUE) ; //Split file default reconstruction
  a->ExecuteTask(deb.Data()) ;
  cout << "AliEMCAL:> Split File default reconstruction ended" << endl ;
  //delete a ; 

  cout << "--------AliEMCAL:> Reconstruction OK------------------"<< endl ;
}
