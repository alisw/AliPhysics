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
  fReaderDate=0;
  fTreeD=0;
  fBt=0;
  fITSHeader=0;
}

//_____________________________________________________________
  AliITSBeamTestDig::AliITSBeamTestDig(const Text_t* name, const Text_t* title): TTask(name,title)
{
  //
  // Standard constructor
  //

  fReaderDate=0;
  fTreeD=0;
  fBt=0;
  fITSHeader=0;
 
}

//______________________________________________________________________
AliITSBeamTestDig::AliITSBeamTestDig(const AliITSBeamTestDig &bt):TTask(bt){
    // Copy constructor. 

  fReaderDate=bt.fReaderDate;
  fTreeD=bt.fTreeD;
  fBt=bt.fBt;
  fITSHeader=bt.fITSHeader;
}
//______________________________________________________________________
AliITSBeamTestDig& AliITSBeamTestDig::operator=(AliITSBeamTestDig &bt){
    // Assignment operator. This is a function which is not allowed to be
    // done to the ITS beam test dig. It exits with an error.
    // Inputs:
    if(this==&bt) return *this;
    Error("operator=","You are not allowed to make a copy of the AliITSBeamTestDig");
    exit(1);
    return *this; //fake return
}


  
 
