// $Id$ 

/**
   Important macro to get certain Aliroot parameters. They are stored
   in a file "Init.cxx". New init of AliL3Transform uses output to read certain
   TPC parameters.
*/

void Make_Init(char *file, char *tofile="Init.cxx"){

  TFile * rootf = new TFile(file,"READ");

  if(!rootf->IsOpen()){
    cerr<<"no file: "<<file<<endl;
    return;
  }

  AliRun *gAlice = (AliRun*)rootf->Get("gAlice");
  if(!gAlice){
    cerr<<"no gAlice in file: "<<file<<endl;
    return;
  }  

  AliTPCParam* par = (AliTPCParam*)rootf->Get("75x40_100x60");
  if(!par){
    cerr<<"no AliTPCParam 75x40_100x60 in file: "<<file<<endl;
    return;
  }

  AliTPCParamSR *param=(AliTPCParamSR*)par;
  AliTPCPRF2D    * prfinner   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter   = new AliTPCPRF2D;
  AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE);
  rf->SetGauss(param->GetZSigma(),param->GetZWidth(),1.);
  rf->SetOffset(3*param->GetZSigma());
  rf->Update();
  
  TDirectory *savedir=gDirectory;
  TFile *if=TFile::Open("$ALICE_ROOT/TPC/AliTPCprf2d.root");
  if (!if->IsOpen()) { 
    cerr<<"Can't open $ALICE_ROOT/TPC/AliTPCprf2d.root !\n" ;
    exit(3);
  }
  prfinner->Read("prf_07504_Gati_056068_d02");
  prfouter->Read("prf_10006_Gati_047051_d03");
  if->Close();
  savedir->cd();
  
  param->SetInnerPRF(prfinner);
  param->SetOuterPRF(prfouter); 
  param->SetTimeRF(rf);
  
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

  fprintf(f,"  fBFieldFactor = %d ;\n",gAlice->Field()->Factor());
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
  fprintf(f,"  fZLength = %.20f ;\n",par->GetZLength());
  fprintf(f,"  fZOffset = %.20f ;\n",par->GetZOffset());
  fprintf(f,"  fDiffT = %.20f ;\n",par->GetDiffT());
  fprintf(f,"  fDiffL = %.20f ;\n",par->GetDiffL());
  fprintf(f,"  fInnerPadLength = %f ;\n",par->GetInnerPadLength());
  fprintf(f,"  fOuterPadLength = %f ;\n",par->GetOuterPadLength());
  fprintf(f,"  fInnerPRFSigma = %.20f ;\n",param->GetInnerPRF()->GetSigmaX());
  fprintf(f,"  fOuterPRFSigma = %.20f ;\n",param->GetOuterPRF()->GetSigmaX());
  fprintf(f,"  fTimeSigma = %.20f ;\n",param->GetTimeRF()->GetSigma());
  
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

void Make_Default(char *file,char *tofile)
{
  /*
    Macro to write out default values, which should be used to initialize
    the static data members of the AliL3Transform class. Macro does more
    or less the same as the above, only the output syntax is changed in order
    to use it for static data member initialization.
  */
  
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
  fprintf(f,"Int_t AliL3Transform::fNTimeBins = %d ;\n",fNTimeBins);
  fprintf(f,"Int_t AliL3Transform::fNRowLow = %d ;\n",fNRowLow);
  fprintf(f,"Int_t AliL3Transform::fNRowUp = %d ;\n",fNRowUp);
  fprintf(f,"Int_t AliL3Transform::fNSectorLow = %d ;\n",fNSectorLow);
  fprintf(f,"Int_t AliL3Transform::fNSectorUp = %d ;\n",fNSectorUp);
  fprintf(f,"Int_t AliL3Transform::fNSector = %d ;\n",fNSector);
  fprintf(f,"Double_t AliL3Transform::fPadPitchWidthLow = %f ;\n",par->GetPadPitchWidth(0));
  fprintf(f,"Double_t AliL3Transform::fPadPitchWidthUp = %f ;\n",par->GetPadPitchWidth(fNSectorLow));
  fprintf(f,"Double_t AliL3Transform::fZWidth = %.20f ;\n",par->GetZWidth());
  fprintf(f,"Double_t AliL3Transform::fZSigma = %.20f ;\n",par->GetZSigma());
  fprintf(f,"Int_t AliL3Transform::fNSlice = %d ;\n",fNSectorLow);
  fprintf(f,"Int_t AliL3Transform::fNRow = %d ;\n",fNRow);
  fprintf(f,"Double_t AliL3Transform::fNRotShift = 0.5 ;\n");
  fprintf(f,"Double_t AliL3Transform::fPi = %.15f ;\n",TMath::Pi());
  fprintf(f,"Double_t AliL3Transform::fX[176] = {\n");
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
    
  fprintf(f,"                                    %3.15f,\n",par->GetPadRowRadii(sec,row));
  }
  fprintf(f,"};\n\n");
  
  fprintf(f,"Int_t AliL3Transform::fNPads[176] = {\n");
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
  fprintf(f,"                                     %d,\n",par->GetNPads(sec,row));
  }
  fprintf(f,"};\n");
  fclose(f);
}
