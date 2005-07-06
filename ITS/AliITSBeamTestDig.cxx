////////////////////////////////////////////////////
//  Class to define                               //
//  beam test raw 2 dig conv.                     //
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                //
////////////////////////////////////////////////////

#include "AliITSBeamTestDig.h"

ClassImp(AliITSBeamTestDig)



//_____________________________________________________________
  AliITSBeamTestDig::AliITSBeamTestDig(): TTask()
{
  //
  // Default constructor
  //
  fReader=0;
  fTreeD=0;
  fITSHeader=0;
  fITSgeom=0;
}

//_____________________________________________________________
  AliITSBeamTestDig::AliITSBeamTestDig(const Text_t* name, const Text_t* title): TTask(name,title)
{
  //
  // Standard constructor
  //

  fReader=0;
  fTreeD=0;
  fITSHeader=0;
  fITSgeom = 0;
}

//______________________________________________________________________
AliITSBeamTestDig::AliITSBeamTestDig(const AliITSBeamTestDig &bt):TTask(bt){
    // Copy constructor. 

  fReader=bt.fReader;
  fTreeD=bt.fTreeD;
  fITSHeader=bt.fITSHeader;
  fITSgeom = bt.fITSgeom;
}
//______________________________________________________________________
AliITSBeamTestDig& AliITSBeamTestDig::operator=(const AliITSBeamTestDig &source){
    // Assignment operator. This is a function which is not allowed to be
    // done to the ITS beam test dig. It exits with an error.
    // Inputs:
    if(this==&source) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITSBeamTestDig");
    exit(1);
    return *this; //fake return
}


  
 
