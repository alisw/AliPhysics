// $Id$

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group


/**
 Example program how to run the "standalone" clusterfinder.

 Important: give patch=-1 for one-patch slices.
*/


#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3Logger.h"
#include "AliL3RootTypes.h"
#include "AliL3Transform.h"
#include "AliL3ClustFinderNew.h"
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"

#if __GNUC__ == 3
using namespace std;
#else
#include <stream.h>
#endif

#include <libgen.h>

#define MAXCLUSTER 25000

int main(Int_t argc,Char_t **argv)
{
  Int_t slice=0;
  Int_t patch=0;
  Int_t fm=4;
  Int_t th=10;
  Bool_t de=kFALSE;

  AliL3Logger l;
  l.Set(AliL3Logger::kAll);
  l.UseStderr();
  //l.UseStdout();
  //l.UseStream();

  if(argc<2){
    cout<<"Usage: runit datafile [slice] [patch] [match] [threshold] [deconv]"<<endl;
    return -1;
  }
  if (argc>2) {
    slice=atoi(argv[2]);
  }
  if (argc>3) {
    patch=atoi(argv[3]);
  }
  if (argc>4) {
    fm=atoi(argv[4]);
  }
  if (argc>5) {
    th=atoi(argv[5]);
  }
  if (argc>6) {
    de=kTRUE;
  }

  AliL3DigitRowData *digits = 0;
  unsigned int nrows=0;
  
  //reading transformer init file
  Char_t fname[1024];
  strcpy(fname,argv[1]);
  AliL3Transform::Init(dirname(fname)); 
  strcpy(fname,argv[1]);

  //Does all the file/data handling  
  AliL3MemHandler file; 

  //Give slice and patch information (see filename convention)
  if((patch>=0)&&(patch<6)) file.Init(slice,patch);
  else {
    Int_t srows[2]={0,AliL3Transform::GetLastRow(5)};
    patch=0;
    file.Init(slice,patch,srows);
  }

  //Open the data file:
  if(!file.SetBinaryInput(fname))
    {
      cerr<<"Error opening file "<<fname<<endl;
      return -1;
    }

  //Store the data in memory, and get the pointer to it:
  digits = file.CompBinary2Memory(nrows);
  file.CloseBinaryInput();

  //The cluster finder itself.
  AliL3ClustFinderNew cf; 

  //Init cluster finder
  cf.InitSlice(slice,patch,0,nrows-1,MAXCLUSTER);
  cf.SetMatchWidth(fm);
  cf.SetThreshold(th);
  //cf.SetXYError(0.2);
  //cf.SetZError(0.3);
  cf.SetSTDOutput(kTRUE);
  cf.SetCalcErr(kTRUE);
  cf.SetDeconv(de); //standard is false

  //Allocate memory to store found spacepoints 
  AliL3MemHandler fpoints;
  AliL3SpacePointData *points=(AliL3SpacePointData*)fpoints.Allocate(MAXCLUSTER*sizeof(AliL3SpacePointData));
  cf.SetOutputArray(points);

  //Give the data pointer to the cluster finder
  cf.Read(nrows,digits);

  //Start processing:
  cf.ProcessDigits();

  exit(0);
}
