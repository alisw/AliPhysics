#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

#include "TString.h"

#include "PHOS/AliPHOSReconstructioner.h"
#endif

void Go(TString deb = ""){
  AliPHOSReconstructioner * a ;   

  cout << "AliPHOS:> Single File default reconstruction started" << endl ;
  a = new AliPHOSReconstructioner("galice.root") ;  //first -single file default reconstruction
  a->ExecuteTask(deb.Data()) ;
  cout << "AliPHOS:> Single File default reconstruction finished" << endl ;
  // delete a ; 

  cout << "AliPHOS:> Single File branch TEST reconstruction started" << endl ;
  a = new AliPHOSReconstructioner("galice.root","test") ;  //another branch single file recontruction
  a->ExecuteTask(deb.Data()) ;
  cout << "AliPHOS:> Single File branch TEST reconstruction ended" << endl ;
  //delete a ; 
  
  cout << "AliPHOS:> Split File default reconstruction started" << endl ;
  a = new AliPHOSReconstructioner("galice.root","Default",kTRUE) ; //Split file default reconstruction
  a->ExecuteTask(deb.Data()) ;
  cout << "AliPHOS:> Split File default reconstruction ended" << endl ;
  //delete a ; 

  cout << "--------AliPHOS:> Reconstruction OK------------------"<< endl ;
}
