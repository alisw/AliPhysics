void AliPHOSReconstructionTest(){
cout << "AliPHOS:> Single File default reconstruction started" << endl ;
AliPHOSReconstructioner a("galice.root") ;  //first -single file default reconstruction
a.ExecuteTask() ;
cout << "AliPHOS:> Single File default reconstruction finished" << endl ;
cout << "AliPHOS:> Single File branch TEST reconstruction started" << endl ;
AliPHOSReconstructioner a("galice.root","test") ;  //another branch single file recontruction
a.ExecuteTask() ;
cout << "AliPHOS:> Single File branch TEST reconstruction ended" << endl ;
cout << "AliPHOS:> Split File default reconstruction started" << endl ;
AliPHOSReconstructioner a("galice.root","Default",kTRUE) ; //Split file default reconstruction
a.ExecuteTask() ;
cout << "AliPHOS:> Split File default reconstruction ended" << endl ;
cout << "--------AliPHOS:> Reconstruction OK------------------"<< endl ;

}
