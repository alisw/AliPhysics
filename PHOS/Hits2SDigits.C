void Hits2SDigits(){
 gSystem->Setenv("CONFIG_SPLIT_FILE","1") ;
 if (gSystem->Getenv("CONFIG_SPLIT_FILE")) 
  cout << "SPLIT" << endl; 
 else
  cout << "NO SPLIT" << endl ; 
 TFile * f = TFile::Open("galice.root","update");
 gAlice = (AliRun*) f->Get("gAlice") ;  
 AliPHOSv1*    fPHOS  = (AliPHOSv1 *)gAlice->GetDetector("PHOS") ;  
 gAlice->SetDebug(10) ; 
 gAlice->Hits2SDigits("PHOS") ;



}
