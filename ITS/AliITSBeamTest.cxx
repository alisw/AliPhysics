////////////////////////////////////////////////////
//  Base class to define                          //
//  ITS beam test                                 //
//                                                //   
//  Origin: E. Crescio crescio@to.infn.it         //
////////////////////////////////////////////////////

#include "AliITSBeamTest.h"

const Int_t AliITSBeamTest::fgkNumberOfSPD=6;
const Int_t AliITSBeamTest::fgkNumberOfSDD=2;
const Int_t AliITSBeamTest::fgkNumberOfSSD=4;


ClassImp(AliITSBeamTest)


//_____________________________________________________________
AliITSBeamTest::AliITSBeamTest() : AliITS()
{
  //
  // Constructor
  //

 
  SetNumberOfSPD(fgkNumberOfSPD);
  SetNumberOfSDD(fgkNumberOfSDD);
  SetNumberOfSSD(fgkNumberOfSSD);
}

//_____________________________________________________________
AliITSBeamTest::AliITSBeamTest(const char* name,const char *title) : AliITS(name,title)
{
  //
  // Constructor
  //

 
  SetNumberOfSPD(fgkNumberOfSPD);
  SetNumberOfSDD(fgkNumberOfSDD);
  SetNumberOfSSD(fgkNumberOfSSD);
}


//__________________________________________________________________
AliITSBeamTest::~AliITSBeamTest()
{
  //
  // Destructor
  //

}

//____________________________________________________________________
Int_t AliITSBeamTest::GetNumberOfSubDet(const TString& det) const{

  if(det.Contains("SPD")) return fNspd;
  if(det.Contains("SDD")) return fNsdd;
  if(det.Contains("SSD")) return fNssd;
  return 0;
}
