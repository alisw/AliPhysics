#include "AliHBTCorrectOSLCorrelFctn.h"
/******************************************************************/
/******************************************************************/
/******************************************************************/

//____________________
///////////////////////////////////////////////////////
//                                                   //
// AliHBTCorrectOSLCorrelFctn                        //
//                                                   //
// Class for calculating Q Invariant correlation     //
// taking to the account resolution of the           //
// detector and coulomb effects.                     //
//                                                   //
///////////////////////////////////////////////////////


AliHBTCorrectOSLCorrelFctn::AliHBTCorrectOSLCorrelFctn(const char* name, const char* title):
 AliHBTOnePairFctn3D(name,title),
 fMeasCorrelFctn(0x0),
 fSmearedNumer(0x0),
 fSmearedDenom(0x0),
 fMeasNumer(0x0),
 fMeasDenom(0x0)
{
//ctor
}
/******************************************************************/

AliHBTCorrectOSLCorrelFctn::AliHBTCorrectOSLCorrelFctn(const AliHBTCorrectOSLCorrelFctn& in):
 AliHBTOnePairFctn3D(in),
 fMeasCorrelFctn(0x0),
 fSmearedNumer(0x0),
 fSmearedDenom(0x0),
 fMeasNumer(0x0),
 fMeasDenom(0x0)
{
//cpy constructor
 in.Copy(*this);
}
/******************************************************************/

AliHBTCorrectOSLCorrelFctn::~AliHBTCorrectOSLCorrelFctn()
{
 //dtor
 delete fMeasCorrelFctn;
 delete fSmearedNumer;
 delete fSmearedDenom;
 delete fMeasNumer;
 delete fMeasDenom;
}

/******************************************************************/
