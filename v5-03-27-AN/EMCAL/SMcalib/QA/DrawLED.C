void DrawLED(int runno=1890,int mod=0,int col=0,int row=0,int gain = 0){
  
  Char_t fname[256];
  sprintf(fname,"LED_%09d.root",runno);
  
  TFile *f = new TFile(fname,"read");
  TTree *tree = f->Get("fTreeAmpVsTime");
  f.ls();
  
  AliCaloCalibSignal *calib = new AliCaloCalibSignal();
       
  int channo = calib->GetChannelNum(mod,col,row,gain);
  char arg[64];
  sprintf(arg,"fChannelNum==%d",channo);
  cout<<arg<<endl;
  tree->Draw("fAmp:fHour",arg);

}



