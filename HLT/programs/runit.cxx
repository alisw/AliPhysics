#include "stream.h"
#include "AliL3ClustFinderNew.h"
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"

//Example program how to run the "standalone" clusterfinder.
int main(int argc,char **argv)
{
  if(argc!=2)
    {
      cout<<"Usage: runit datafile"<<endl;
      return -1;
    }
  
  AliL3DigitRowData *digits = 0;
  unsigned int ndigits=0;
  
  AliL3MemHandler file; //Does all the file/data handling
  //AliL3Transform::Init(path); //Storing all detector-spesific quantities, needed by the clusterfinder.
  AliL3ClustFinderNew cf; //The cluster finder itself.

  //Open the data file:
  if(!file.SetBinaryInput(argv[1]))
    {
      cerr<<"Error opening file "<<argv[1]<<endl;
      return -1;
    }

  //Allocate memory to store found spacepoints 
  AliL3MemHandler fpoints;
  AliL3SpacePointData *points=(AliL3SpacePointData*)fpoints.Allocate(10000*sizeof(AliL3SpacePointData));
  cf.SetOutputArray(points);

  //Store the data in memory, and get the pointer to it:
  digits = file.CompBinary2Memory(ndigits);
  file.CloseBinaryInput();
  
  //Switch off deconvolution:
  cf.SetDeconv(false);
  
  //Init cluster finder
  cf.InitSlice(0,0,0,20,10000);
  cf.SetXYError(0.2);
  cf.SetZError(0.3);

  //Give the data pointer to the cluster finder
  cf.Read(ndigits,digits);
  
  //Start processing:
  cf.ProcessDigits();
  
  return 0;
}





