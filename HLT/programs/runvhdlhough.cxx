// $Id$
   
// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group


#include "AliHLTStandardIncludes.h"

#include "AliHLTRootTypes.h"
#include "AliHLTTransform.h"
#include "AliHLTVHDLClusterFinder.h"
#include "AliHLTAltroMemHandler.h"
#include "AliHLTLogging.h"
#include "AliHLTLogger.h"

#if __GNUC__ == 3
using namespace std;
#else
#include <stream.h>
#endif


/**
Example program how to run the vhdl hough code.
*/

int main(int argc,char **argv)
{
#if 0
  Int_t slice=0;
  Int_t patch=0;
  Int_t fm=4;
  Int_t th=10;

  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStderr();
  //l.UseStdout();
  //l.UseStream();

  if(argc<2){
    cout<<"Usage: runvhdlcf altrodatafile [slice] [patch] [matchwidth] [threshold]"<<endl;
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

  //Storing all specific quantities, needed by the Cluster Finder.
  //Char_t fname[1024];
  //strcpy(fname,argv[1]);
  //AliHLTTransform::Init(dirname(fname)); 
  //strcpy(fname,argv[1]);

  FILE *afile=fopen(argv[1],"r");
  if(!afile) {
    cerr << "Can't open file " << argv[1] << endl;
    exit(1);
  }

  AliHLTVHDLClusterFinder cf;
  cf.SetASCIIInput(afile);

  //set cluster finder parameters
  cf.SetMatchWidth(fm);
  cf.SetThreshold(th);
  //cf.SetXYError(0.2);
  //cf.SetZError(0.3);
  cf.SetSTDOutput(kTRUE);
  cf.SetCalcErr(kTRUE);
  //cf.SetDeconv(kFALSE);
  
  //start processing data
  cf.ProcessDigits();

  if(afile) fclose(afile);
  exit(1);
#endif
}





