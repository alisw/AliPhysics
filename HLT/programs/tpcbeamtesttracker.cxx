// $Id$

#include "AliL3StandardIncludes.h"
#include "AliL3RootTypes.h"
#include "AliLevel3.h"
#include "AliL3Transform.h"
#include "AliL3RawDataFileHandler.h"
#include "AliL3SpacePointData.h"
#include "AliL3ClustFinderNew.h"
#include "AliL3ConfMapper.h"
#include "AliL3Vertex.h"

#if __GNUC__== 3
using namespace std;
#else
#include <stream.h>
#include <string.h>
#include <stdlib.h>
#endif


#include <sys/stat.h>
struct stat stat_results;


//This program does the clusterfinding and the tracking.
int main(Int_t argc,Char_t **argv)
{

  Char_t cfile[1024];
  Char_t path[1024]; 
  Int_t slice=0;
  Int_t patch=-1;
  
  if(argc<3){
    cout<<"Usage: tpcbeamtesttracker filename path_to_cosmics"<<endl;
    return -1;
  }
  if (argc>2) {
    sprintf(cfile,"%s",argv[1]);
    sprintf(path,"%s",argv[2]);
  }

  AliL3Transform::Init("./l3-cosmics-transform.config");
  AliL3Transform::SetZeroSup(5);
  AliL3RawDataFileHandler *f=new AliL3RawDataFileHandler();

  f->Init(slice,patch);
  
  f->SetMappingFile("remap.txt");
  f->ReadMappingFile();

  Char_t fname[1024];
  Char_t pname[1024];

  f->SetRawPedestalsInput("./pedout.out");
  f->ReadRawPedestalsInput();
  //f->SetPedVal(85);

  sprintf(pname,"%s/%s",path,cfile);
#if 0
  f->SetRawInput(pname);
  f->ReadRawInput();
#else
  if (stat(pname, &stat_results) != 0){
    exit(1);
  }
  Char_t *t=new Char_t[stat_results.st_size];
  ifstream *fin=new ifstream(pname,fstream::binary);
  fin->read(t,stat_results.st_size);
  delete fin;
  f->ReadRawInputPointer(t);
#endif

  UInt_t nrows;
  AliL3DigitRowData *data=(AliL3DigitRowData*)f->RawData2Memory(nrows);

  AliL3MemHandler *out = new AliL3MemHandler();
  sprintf(fname,"./%s-digits_%d_%d.raw",cfile,slice,patch);
  out->SetBinaryOutput(fname);
  out->Memory2Binary(nrows,data);
  out->CloseBinaryOutput();
  out->Free();

  AliL3MemHandler *mem=new AliL3MemHandler();

  UInt_t maxclusters=100000;
  UInt_t pointsize = maxclusters*sizeof(AliL3SpacePointData);
  AliL3SpacePointData *points = (AliL3SpacePointData*)mem->Allocate(pointsize);

  Bool_t rawsp=kTRUE;
  AliL3ClustFinderNew *cf = new AliL3ClustFinderNew();
  cf->InitSlice(slice,patch,maxclusters);
  cf->SetMatchWidth(10);
  cf->SetSTDOutput(kTRUE);
  cf->SetRawSP(rawsp);
  cf->SetThreshold(10);
  cf->SetDeconv(kFALSE);
  cf->SetCalcErr(kTRUE);
  cf->SetOutputArray(points);
  cf->Read(nrows,data);
  cf->ProcessDigits();
  Int_t npoints = cf->GetNumberOfClusters();

  sprintf(fname,"./%s-points_%d_%d.raw",cfile,slice,patch);
  out->Transform(npoints,points,slice);
  out->SetBinaryOutput(fname);
  out->Memory2Binary(npoints,points);
  out->CloseBinaryOutput();
  out->Free();

  delete out;
  delete mem;
  delete cf;
  delete f;


}
