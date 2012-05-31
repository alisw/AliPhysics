void SDigits2Digits(){
 gSystem->Setenv("CONFIG_SPLIT_FILE","1") ;
 
 TFile * f = TFile::Open("galice.root","update") ;
 gAlice = (AliRun*) f->Get("gAlice") ;  
 AliPHOSv1*    fPHOS  = (AliPHOSv1 *)gAlice->GetDetector("PHOS") ;  
 gAlice->SetDebug(10) ; 
 gAlice->SDigits2Digits("PHOS") ;



}
