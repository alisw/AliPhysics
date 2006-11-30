// $Id$

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group


/**
   This program extracts parameters and lookup tables needed for the
   vhdl implementation of the Hough transform. 
*/

#include <AliHLTStandardIncludes.h>
#include <AliHLTRootTypes.h>
#include <AliHLTTransform.h>
#include <AliHLTLogging.h>
#include <AliHLTLogger.h>
#include <AliHLTMemHandler.h>
#include <AliHLTHoughTransformerVhdl.h>

#if __GNUC__ == 3
using namespace std;
#endif

int main(Int_t argc,Char_t **argv)
{
  Int_t patch=0;
  Int_t slice=0;
  Char_t path[1000];

  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStderr();
  //l.UseStdout();
  //l.UseStream();
    
  if (argc>1) {
    slice=atoi(argv[1]);
  }
  if (argc>2) {
    patch=atoi(argv[2]);
  }  
  if (argc>3) {
    strcpy(path,argv[3]);
  } else strcpy(path,"/tmp/data/RawData/slice0");
  if(argc>4){
    cout<<"Usage: transform [slice] [patch] [path]"<<endl;
    exit(1);
  }

  AliHLTTransform::Init(path);
  cerr << "Transform version: " << AliHLTTransform::GetVersion() << endl;

  AliHLTHoughTransformerVhdl vtest(slice,patch,100,10);
  vtest.CreateHistograms(64,0.1,64,-30,30);
  vtest.PrintVhdl();
  exit(0);
}
  
