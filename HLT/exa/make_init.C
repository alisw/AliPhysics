// $Id$ 

/**
   Macro to get the default parameters from AliTPCParam. 
   Output is written to a file, which can be inserted directly into 
   AliL3Transform.cxx as the default parameters to be set at top.
*/

void make_init(char *tofile="Init.cxx")
{
  AliTPCParamSR *param = new AliTPCParamSR();
  param->SetDefault();
  if(!param)
    {
      cerr<<"AliL3Transform::MakeInitFile : No TPC parameters found"<<endl;
      return;
    }
  
  AliTPCPRF2D    * prfinner   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter1   = new AliTPCPRF2D;
  AliTPCPRF2D    * prfouter2   = new AliTPCPRF2D;  
  AliTPCRF1D     * rf    = new AliTPCRF1D(kTRUE);
  rf->SetGauss(param->GetZSigma(),param->GetZWidth(),1.);
  rf->SetOffset(3*param->GetZSigma());
  rf->Update();
  
  TDirectory *savedir=gDirectory;
  TFile *f1=TFile::Open("$ALICE_ROOT/TPC/AliTPCprf2d.root");
  if (!f1->IsOpen()) 
    { 
      cerr<<"Can't open $ALICE_ROOT/TPC/AliTPCprf2d.root !\n" ;
      exit(3);
    }
  prfinner->Read("prf_07504_Gati_056068_d02");
  prfouter1->Read("prf_10006_Gati_047051_d03");
  prfouter2->Read("prf_15006_Gati_047051_d03");  
  f1->Close();
  savedir->cd();
  
  param->SetInnerPRF(prfinner);
  param->SetOuter1PRF(prfouter1); 
  param->SetOuter2PRF(prfouter2);
  param->SetTimeRF(rf);
  
  Int_t fNTimeBins = param->GetMaxTBin()+1;
  Int_t fNRowLow = param->GetNRowLow();
  Int_t fNRowUp  = param->GetNRowUp();
  Int_t fNRowUp1 = param->GetNRowUp1();
  Int_t fNRowUp2 = param->GetNRowUp2();
  Int_t fNRow= fNRowLow + fNRowUp;
  Int_t fNSectorLow = param->GetNInnerSector();
  Int_t fNSectorUp = param->GetNOuterSector();
  Int_t fNSector = fNSectorLow + fNSectorUp;
  Int_t fNSlice = fNSectorLow;
  
  FILE *f = fopen(tofile,"w");
  if(!f)
    {
      cerr<<"Error opening file "<<tofile<<endl;
      return;
    }
  fprintf(f,"const Double_t AliL3Transform::fBFACT = 0.0029980;\n");
  fprintf(f,"Double_t AliL3Transform::fBField = 0.2;\n");
  fprintf(f,"Int_t AliL3Transform::fVersion = 0;\n");
  fprintf(f,"Int_t AliL3Transform::fBFieldFactor = %d ;\n",gAlice->Field()->Factor());
  fprintf(f,"Int_t AliL3Transform::fNTimeBins = %d ;\n",fNTimeBins);
  fprintf(f,"Int_t AliL3Transform::fNRowLow = %d ;\n",fNRowLow);
  fprintf(f,"Int_t AliL3Transform::fNRowUp = %d ;\n",fNRowUp);
  fprintf(f,"Int_t AliL3Transform::fNRowUp1 = %d ;\n",fNRowUp1);
  fprintf(f,"Int_t AliL3Transform::fNRowUp2 = %d ;\n",fNRowUp2);
  fprintf(f,"Int_t AliL3Transform::fNSectorLow = %d ;\n",fNSectorLow);
  fprintf(f,"Int_t AliL3Transform::fNSectorUp = %d ;\n",fNSectorUp);
  fprintf(f,"Int_t AliL3Transform::fNSector = %d ;\n",fNSector);
  fprintf(f,"Double_t AliL3Transform::fPadPitchWidthLow = %f ;\n",param->GetInnerPadPitchWidth());
  fprintf(f,"Double_t AliL3Transform::fPadPitchWidthUp = %f ;\n",param->GetOuterPadPitchWidth());
  fprintf(f,"Double_t AliL3Transform::fZWidth = %.20f ;\n",param->GetZWidth());
  fprintf(f,"Double_t AliL3Transform::fZSigma = %.20f ;\n",param->GetZSigma());
  fprintf(f,"Double_t AliL3Transform::fZLength = %.20f ;\n",param->GetZLength());
  fprintf(f,"Double_t AliL3Transform::fZOffset = %.20f ;\n",param->GetZOffset());
  fprintf(f,"Double_t AliL3Transform::fDiffT = %.20f ;\n",param->GetDiffT());
  fprintf(f,"Double_t AliL3Transform::fDiffL = %.20f ;\n",param->GetDiffL());
  fprintf(f,"Double_t AliL3Transform::fInnerPadLength = %f ;\n",param->GetInnerPadLength());
  fprintf(f,"Double_t AliL3Transform::fOuter1PadLength = %f ;\n",param->GetOuter1PadLength());
  fprintf(f,"Double_t AliL3Transform::fOuter2PadLength = %f ;\n",param->GetOuter2PadLength());
  fprintf(f,"Double_t AliL3Transform::fInnerPRFSigma = %.20f ;\n",param->GetInnerPRF()->GetSigmaX());
  fprintf(f,"Double_t AliL3Transform::fOuter1PRFSigma = %.20f ;\n",param->GetOuter1PRF()->GetSigmaX());
  fprintf(f,"Double_t AliL3Transform::fOuter2PRFSigma = %.20f ;\n",param->GetOuter2PRF()->GetSigmaX());
  fprintf(f,"Double_t AliL3Transform::fTimeSigma = %.20f ;\n",param->GetTimeRF()->GetSigma());
  fprintf(f,"Int_t AliL3Transform::fNSlice = %d ;\n",fNSectorLow);
  fprintf(f,"Int_t AliL3Transform::fNRow = %d ;\n",fNRow);
  fprintf(f,"Int_t AliL3Transform::fNRotShift = 0.5 ;\n");
  fprintf(f,"Double_t AliL3Transform::fPi = %.15f ;\n",TMath::Pi());
    
  fprintf(f,"Double_t AliL3Transform::fX[159] = {");
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
    
    fprintf(f,"                                    %3.15f,\n",param->GetPadRowRadii(sec,row));
  }
  fprintf(f,"};\n\n");
  
  fprintf(f,"Int_t AliL3Transform::fNPads[159] = {");
  for(Int_t i=0;i<fNRow;i++){
    int sec,row;
    if( i < fNRowLow){sec =0;row =i;}
    else{sec = fNSectorLow;row =i-fNRowLow;}
  fprintf(f,"                                     %d,\n",param->GetNPads(sec,row));
  }
    
  fprintf(f,"};\n");
  fclose(f);
}
