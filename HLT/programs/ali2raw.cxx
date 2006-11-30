// $Id$ 

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#ifndef use_aliroot
BOMB THE COMPILE: USEPACKAGE=ALIROOT
#endif

#include "AliHLTStandardIncludes.h"

#include "AliHLTRootTypes.h"
#include "AliHLTMemHandler.h"
#include "AliHLTFileHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"

#if __GNUC__ == 3
using namespace std;
#else
#include <stream.h>
#endif

/**
   Converts from AliRoot Digits to L3-RawDigits
 */
Int_t Make_Init(Char_t *file, Char_t *tofile);

int main (Int_t argc, Char_t** argv)
{
  Int_t first=0;
  Int_t last=0;

  if(argc<3)
  {
    cout<<"Usage: ali2raw datafile path [firstsector] [lastsector]"<<endl;
    return -1;
  }
  if(argc<4) first=atoi(argv[3]);
  if(argc<5) first=atoi(argv[4]);

  Char_t transname[1000];
  const Int_t npatch = 6;
  Int_t row[npatch][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
  
  strcpy(transname,argv[1]);
  strcat(transname,"/.transform.config");
  if(Make_Init(argv[1],transname)<0) exit(1);


  /*
  AliHLTFileHandler *fFileHandler = new AliHLTFileHandler(); 
  fFileHandler->SetAliInput(in);
  AliHLTTransform *fTransformer = new AliHLTTransform(path);


  fFileHandler->Init(fTransformer);
  for(int slice=first; slice<=last; slice++){
    for(int patch=0;patch<npatch;patch++){
      cerr<<"reading slice: "<<slice<<" patch: "<<patch;
      fFileHandler->Free();
      fFileHandler->Init(slice,patch,row[patch]);      
      sprintf(name,"%sdigits_%d_%d.raw",path,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary();
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }      
  }
  fFileHandler->CloseAliInput();
  */
}

int Make_Init(char *file, char *tofile=NULL){

  TFile *rootf = new TFile(file,"READ");

  if(!rootf->IsOpen()){
    cerr<<"could not open file: "<<file<<endl;
    return -1;
  }

  AliTPCParam* par = (AliTPCParam*)rootf->Get("75x40_100x60");

  if(!par){
    cerr<<"no AliTPCParam 75x40_100x60 in file: "<<file<<endl;
    return -2;
  }

  Int_t fNTimeBins = par->GetMaxTBin()+1;
  Int_t fNRowLow = par->GetNRowLow();
  Int_t fNRowUp  = par->GetNRowUp();
  Int_t fNRow= fNRowLow + fNRowUp;
  Int_t fNSectorLow = par->GetNInnerSector();
  Int_t fNSectorUp = par->GetNOuterSector();
  Int_t fNSector = fNSectorLow + fNSectorUp;

  Char_t tofilename[1000];
  if(tofile){
    strcpy(tofilename,tofile);
    strcat(tofilename,".transform-config");
  } else strcpy(tofilename,file);

  FILE *f = fopen(tofilename,"w");
  fprintf(f,"void AliHLTTransform::Init(){\n");

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
  return 0;
}



#if 0
void singlepatch(Char_t* in,Int_t first, Int_t last,Char_t *path="",Int_t event=0)
{
  AliHLTLogger l;
  //l.UnSet(AliHLTLogger::kDebug);
  //l.UnSet(AliHLTLogger::kAll);
  //l.Set(AliHLTLogger::kInformational);
  //l.UseStdout();
  l.UseStream();
  
  Char_t fname[100];
  sprintf(fname,"%sevent_%d/",path,event);
  Char_t name[256];
  AliHLTFileHandler *fFileHandler = new AliHLTFileHandler(); 
  fFileHandler->SetAliInput(in);
  Int_t srow[2] = {0,175};
  Int_t patch=0;
  for(Int_t slice=first; slice<=last; slice++)
    {
      cerr<<"reading slice: "<<slice;
      fFileHandler->Free();
      fFileHandler->Init(slice,patch,srow);
      sprintf(name,"%sdigits_%d_%d.raw",fname,slice,patch);
      fFileHandler->SetBinaryOutput(name);
      fFileHandler->AliDigits2CompBinary(event);
      fFileHandler->CloseBinaryOutput();      
      cerr<<" done"<<endl;
    }
  fFileHandler->CloseAliInput();
  
}
#endif







