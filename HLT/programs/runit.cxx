/* $Id$
   Author: Constantin Loizides <loizides@ikf.physik.uni-frankfurt.de>
*/

#include <stream.h>
#include <libgen.h>
#include "AliL3Transform.h"
#include "AliL3ClustFinderNew.h"
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"

/**
 Example program how to run the "standalone" clusterfinder.
*/

int main(int argc,char **argv)
{
  if(argc!=2)
    {
      cout<<"Usage: runit datafile"<<endl;
      return -1;
    }

  AliL3DigitRowData *digits = 0;
  unsigned int ndigits=0;
  
  //Does all the file/data handling  
  AliL3MemHandler file; 

  //Open the data file:
  if(!file.SetBinaryInput(argv[1]))
    {
      cerr<<"Error opening file "<<argv[1]<<endl;
      return -1;
    }

  //Store the data in memory, and get the pointer to it:
  digits = file.CompBinary2Memory(ndigits);
  file.CloseBinaryInput();

  //Storing all detector-spesific quantities, needed by the clusterfinder.
  AliL3Transform::Init(dirname(argv[1])); 

  //The cluster finder itself.
  AliL3ClustFinderNew cf; 

  //Init cluster finder
  cf.InitSlice(0,0,0,ndigits-1,10000);
  cf.SetXYError(0.2);
  cf.SetZError(0.3);
  cf.SetSTDOutput(kTRUE);

  //Switch off deconvolution:
  cf.SetDeconv(false);
  
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

