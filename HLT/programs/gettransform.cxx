// $Id$

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group


/**
   This program extracts parameters and lookup tables needed for the
   vhdl implementation of the Hough transform. 
*/

#include <AliL3StandardIncludes.h>
#include <AliL3RootTypes.h>
#include <AliL3Transform.h>
#include <AliL3Logging.h>
#include <AliL3Logger.h>
#include <AliL3MemHandler.h>
#include <AliL3HoughTransformerVhdl.h>

#if __GNUC__ == 3
using namespace std;
#endif

int main(Int_t argc,Char_t **argv)
{
  Int_t patch=0;
  Int_t slice=0;
  Char_t path[1000];

  AliL3Logger l;
  l.Set(AliL3Logger::kAll);
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

  AliL3Transform::Init(path);
  cerr << "Transform version: " << AliL3Transform::GetVersion() << endl;

  AliL3HoughTransformerVhdl vtest(slice,patch,100,10);
  vtest.CreateHistograms(64,0.1,64,-30,30);
  vtest.PrintVhdl();
  exit(0);
}
  
