/** $Id$ 

Important macro to get certain Aliroot parameters. They are stored
in a file "Init.cxx". Compare the contents of the class AliL3Transform
with the result of this macro to check that there are no differences.
*/

void Make_Init(char *file, char *tofile="Init.cxx"){

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

  int fNTimeBins = par->GetMaxTBin()+1;
  int fNRowLow = par->GetNRowLow();
  int fNRowUp  = par->GetNRowUp();
  int fNRow= fNRowLow + fNRowUp;
  int fNSectorLow = par->GetNInnerSector();
  int fNSectorUp = par->GetNOuterSector();
  int fNSector = fNSectorLow + fNSectorUp;
  int fNSlice = fNSectorLow;

  FILE *f = fopen(tofile,"w");
  fprintf(f,"void AliL3Transform::Init(){\n");

  fprintf(f,"  //sector:\n");
  fprintf(f,"  fNTimeBins = %d ;\n",fNTimeBins);
  fprintf(f,"  fNRowLow = %d ;\n",fNRowLow);
  fprintf(f,"  fNRowUp = %d ;\n",fNRowUp);
  fprintf(f,"  fNSectorLow = %d ;\n",fNSectorLow);
  fprintf(f,"  fNSectorUp = %d ;\n",fNSectorUp);
  fprintf(f,"  fNSector = %d ;\n",fNSector);
  fprintf(f,"  fPadPitchWidthLow = %f ;\n",par->GetPadPitchWidth(0));
  fprintf(f,"  fPadPitchWidthUp = %f ;\n",par->GetPadPitchWidth(fNSectorLow));
  fprintf(f,"  fZWidth = %.20f ;\n",par->GetZWidth());
  fprintf(f,"  fZSigma = %.20f ;\n",par->GetZSigma());

  fprintf(f,"\n  //slices:\n");
  fprintf(f,"  fNSlice = %d ;\n",fNSectorLow);
  fprintf(f,"  fNRow = %d ;\n",fNRow);

  //rotation shift put in by hand -> Constantin 
  fprintf(f,"  fNRotShift = 0.5 ;\n");

  fprintf(f,"  fPi = %.15f ;\n",TMath::Pi());
  fprintf(f,"  for(Int_t i=0;i<36;i++){\n");
  fprintf(f,"    fCos[i] = cos(2*fPi/9*(i+0.5));\n");
  fprintf(f,"    fSin[i] = sin(2*fPi/9*(i+0.5));\n");
  fprintf(f,"  }\n\n");

  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
    fprintf(f,"  fX[%d] = %3.15f ;\n",i,par->GetPadRowRadii(sec,row));
  }
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
    fprintf(f,"  fNPads[%d] = %d ;\n",i,par->GetNPads(sec,row));
  }

  fprintf(f,"}\n");
  fclose(f);
}

