/* $Id$
   Author: Constantin Loizides <loizides@ikf.physik.uni-frankfurt.de>
*/

#include <stream.h>
#include <libgen.h>
#include "AliL3RootTypes.h"
#include "AliL3Transform.h"
#include "AliL3ClustFinderNew.h"
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"

/**
 Example program how to run the "standalone" clusterfinder.
*/

int main(int argc,char **argv)
{
  Int_t slice=0;
  Int_t patch=0;

  /*
    AliL3Logger l;
    l.Set(AliL3Logger::kAll);
    l.UseStdout();
    //l.UseStream();
  */

  if(argc<2){
    cout<<"Usage: runit datafile [slice] [patch]"<<endl;
    return -1;
  }
  if (argc>2) {
    slice=atoi(argv[2]);
  }
  if (argc>3) {
    patch=atoi(argv[3]);
  }

  AliL3DigitRowData *digits = 0;
  unsigned int ndigits=0;
  
  //Storing all specific quantities, needed by the Cluster Finder.
  Char_t fname[1024];
  strcpy(fname,argv[1]);
  AliL3Transform::Init(dirname(fname)); 
  strcpy(fname,argv[1]);

  //Does all the file/data handling  
  AliL3MemHandler file; 

  //Give slice and patch information (see filename convention)
  if((patch>=0)&&(patch<6)) file.Init(slice,patch);
  else {
    Int_t srows[2]={0,175};
    file.Init(slice,0,srows);
  }

  //Open the data file:
  if(!file.SetBinaryInput(fname))
    {
      cerr<<"Error opening file "<<fname<<endl;
      return -1;
    }

  //Store the data in memory, and get the pointer to it:
  digits = file.CompBinary2Memory(ndigits);
  file.CloseBinaryInput();

  //The cluster finder itself.
  AliL3ClustFinderNew cf; 

  //Init cluster finder
  cf.InitSlice(0,0,0,ndigits-1,10000);
  //cf.SetXYError(0.2);
  //cf.SetZError(0.3);
  cf.SetSTDOutput(kTRUE);
  cf.SetCalcErr(kTRUE);

  //Switch off deconvolution:
  cf.SetDeconv(kFALSE);
  
  //Allocate memory to store found spacepoints 
  AliL3MemHandler fpoints;
  AliL3SpacePointData *points=(AliL3SpacePointData*)fpoints.Allocate(10000*sizeof(AliL3SpacePointData));
  cf.SetOutputArray(points);

  //Give the data pointer to the cluster finder
  cf.Read(ndigits,digits);

  //Start processing:
  cf.ProcessDigits();
  
  return 0;
}
