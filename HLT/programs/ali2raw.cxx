/* $Id$ 
   Author: Constantin Loizides <loizides@ikf.physik.uni-frankfurt.de>
*/

#ifndef use_aliroot
BOMB THE COMPILE: USEPACKAGE=ALIROOT
#endif

#include <stream.h>
#include "AliL3RootTypes.h"
#include "AliL3MemHandler.h"
#include "AliL3FileHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Transform.h"

/**
   Converts from AliRoot Digits to L3-RawDigits
 */
int Make_Init(char *file, char *tofile=NULL);

int main (int argc, char** argv)
{
  int first=0;
  int last=0;

  if(argc<3)
  {
    cout<<"Usage: ali2raw datafile path [firstsector] [lastsector]"<<endl;
    return -1;
  }
  if(argc<4) first=atoi(argv[3]);
  if(argc<5) first=atoi(argv[4]);

  char transname[1000];
  const Int_t npatch = 6;
  Int_t row[npatch][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
  
  strcpy(transname,argv[1]);
  strcat(transname,"/.transform.config");
  if(Make_Init(argv[1],transname)<0) exit(1);


  /*
  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);
  AliL3Transform *fTransformer = new AliL3Transform(path);


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

  int fNTimeBins = par->GetMaxTBin()+1;
  int fNRowLow = par->GetNRowLow();
  int fNRowUp  = par->GetNRowUp();
  int fNRow= fNRowLow + fNRowUp;
  int fNSectorLow = par->GetNInnerSector();
  int fNSectorUp = par->GetNOuterSector();
  int fNSector = fNSectorLow + fNSectorUp;
  //  int fNSlice = fNSectorLow;

  char tofilename[1000];
  if(tofile){
    strcpy(tofilename,tofile);
    strcat(tofilename,".transform-config");
  } else strcpy(tofilename,file);

  FILE *f = fopen(tofilename,"w");
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
  return 0;
}



#if 0
void singlepatch(char* in,int first, int last,char *path="",int event=0)
{
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  //l.UseStdout();
  l.UseStream();
  
  char fname[100];
  sprintf(fname,"%sevent_%d/",path,event);
  char name[256];
  AliL3FileHandler *fFileHandler = new AliL3FileHandler(); 
  fFileHandler->SetAliInput(in);
  Int_t srow[2] = {0,175};
  int patch=0;
  for(int slice=first; slice<=last; slice++)
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







