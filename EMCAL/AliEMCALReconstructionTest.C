void AliEMCALReconstructionTest(){

  cout << "AliEMCAL:> Single File default reconstruction started" << endl ;

  AliEMCALReconstructioner a("galice.root") ;  //first -single file default reconstruction
  a.ExecuteTask() ;
  cout << "AliEMCAL:> Single File default reconstruction finished" << endl ;
 
  cout << "AliEMCAL:> Single File branch TEST reconstruction started" << endl ;
  AliEMCALReconstructioner a("galice.root","test") ;  //another branch single file recontruction
  a.ExecuteTask() ;
  cout << "AliEMCAL:> Single File branch TEST reconstruction ended" << endl ;
  
  cout << "AliEMCAL:> Split File default reconstruction started" << endl ;
  AliEMCALReconstructioner a("galice.root","Default",kTRUE) ; //Split file default reconstruction
  a.ExecuteTask() ;
  cout << "AliEMCAL:> Split File default reconstruction ended" << endl ;
  
  cout << "--------AliEMCAL:> Reconstruction OK------------------"<< endl ;
  
}
