////////////////////////////////////////////////////
//  Class to define                               //
//  beam test raw 2 dig conv.                     //
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                //
////////////////////////////////////////////////////

#include "AliITSBeamTestDig.h"

ClassImp(AliITSBeamTestDig)



//_____________________________________________________________
  AliITSBeamTestDig::AliITSBeamTestDig(): TTask(),
fITSHeader(0),
fReader(0),
fTreeD(0),
fITSgeom(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________
AliITSBeamTestDig::AliITSBeamTestDig(const Text_t* name, const Text_t* title): TTask(name,title),
fITSHeader(0),
fReader(0),
fTreeD(0),
fITSgeom(0)
{
  //
  // Standard constructor
  //

}

//______________________________________________________________________
AliITSBeamTestDig::AliITSBeamTestDig(const AliITSBeamTestDig &bt):TTask(bt),
fITSHeader(bt.fITSHeader),
fReader(bt.fReader),
fTreeD(bt.fTreeD),
fITSgeom(bt.fITSgeom){
    // Copy constructor. 

}
//______________________________________________________________________
AliITSBeamTestDig& AliITSBeamTestDig::operator=(const AliITSBeamTestDig &source){
    // Assignment operator
  this->~AliITSBeamTestDig();
  new(this) AliITSBeamTestDig(source);
  return *this;
}


  
 
