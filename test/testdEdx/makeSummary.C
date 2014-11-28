void DrawSysWatchTime(){
  TString prefix="/hera/alice/marsland/MAF/MAFbenchmark/mcmaker/workdir/test_6_375000_25_20140713_20/";
  TString listSyswatch = gSystem->GetFromPipe(Form("ls %s/*/mult_10000/event_1/simwatch.log",prefix.Data()));  
  TObjArray * arraySyswatch= listSyswatch.Tokenize("\n");
  Int_t nfiles= arraySyswatch->GetEntries();
  TTree * treeSyswatch[] = new TTree*[nfiles];

  
  for (Int_t ifile=0; ifile<nfiles; ifile++){
    treeSyswatch[ifile] = AliSysInfo::MakeTree(arraySyswatch->At(ifile)->GetName());    
    treeSyswatch[0]->AddFriend( treeSyswatch[ifile],Form("T%d",ifile));
  }

  treeSyswatch[0]->Draw("deltaT:sname","deltaT>1","");
  for (Int_t ifile=1; ifile<nfiles; ifile++){
    treeSyswatch[0]->SetMarkerStyle(25);
    treeSyswatch[0]->SetMarkerColor(1+ifile);
    treeSyswatch[0]->Draw(Form("T%d.stampSec-T%d.stampOldSec:sname",ifile,ifile),"deltaT>1","same");    
  }  
}

 
