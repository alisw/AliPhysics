void make_init(char *file){

  TFile * rootf = new TFile(file,"READ");

  if(!rootf->IsOpen()){
    cerr<<"no file: "<<file<<endl;
    return;
  }

  AliTPCParam* par = (AliTPCParam*)rootf->Get("75x40_100x60");

  if(!par){
    cerr<<"no AliTPCParam 75x40_100x60 in file: "<<file<<endl;
    return;
  }

  int fNRowLow = par->GetNRowLow();
  int fNRowUp  = par->GetNRowUp();
  int fNRow= fNRowLow+ fNRowUp;
  int fNSectorLow = par->GetNInnerSector();
  int fNSectorUp = par->GetNOuterSector();
  int fNSector = fNSectorLow + fNSectorUp;
  int fNSlice = fNSectorLow;

  FILE *f = fopen("Init.cxx","w");
  fprintf(f,"void AliL3Transform::Init(){\n");
  fprintf(f,"  //sector:\n");
  fprintf(f,"  fNRowLow = %d;\n",fNRowLow);
  fprintf(f,"  fNRowUp = %d;\n",fNRowUp);
  fprintf(f,"  fNSectorLow = %d;\n",fNSectorLow);
  fprintf(f,"  fNSectorUp = %d;\n",fNSectorUp);
  fprintf(f,"  fNSector = %d;\n",fNSector);
  fprintf(f,"  fPadPitchWidthLow = %f;\n",par->GetPadPitchWidth(0));
  fprintf(f,"  fPadPitchWidthUp = %f;\n",par->GetPadPitchWidth(fNSectorLow));
  fprintf(f,"  fZWidth = %.20f;\n",par->GetZWidth());
  fprintf(f,"  fZSigma = %.20f;\n",par->GetZSigma());
  fprintf(f,"\n  //slices:\n");
  fprintf(f,"  fNSlice = %d;\n",fNSectorLow);
  fprintf(f,"  fNRow = %d;\n",fNRow);
//  fprintf(f,"  fPi = 3.14159265358979323846;\n");
  fprintf(f,"  fPi = %.15f;\n",TMath::Pi());
  fprintf(f,"  for(Int_t i=0;i<36;i++){\n");
  fprintf(f,"    fCos[i] = cos(fPi/9*i);\n");
  fprintf(f,"    fSin[i] = sin(fPi/9*i);\n");
  fprintf(f,"  }\n\n");
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
    fprintf(f,"  fX[%d] = %3.15f;\n",i,par->GetPadRowRadii(sec,row));
  }
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
    fprintf(f,"  fNPads[%d] = %d;\n",i,par->GetNPads(sec,row));
  }

  fprintf(f,"}\n");
  fclose(f);
}
