void cphisto(){

  TFile * f = new TFile("histos.root") ;
  TList * l=(TList*)f->Get("histESD") ;
  cout << l->GetSize() << endl ;

  TFile fout("sum.root","recreate") ;
  for(Int_t i=0;i<l->GetSize();i++){
    l->At(i)->Write() ;
  }
  fout.Close() ;

}
