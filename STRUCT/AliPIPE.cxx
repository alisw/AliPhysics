///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Beam pipe class                                                          //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliPIPEClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliPIPE.h"
#include "AliRun.h"
 
ClassImp(AliPIPE)
 
//_____________________________________________________________________________
AliPIPE::AliPIPE()
{
  //
  // Default constructor for beam pipe
  //
}
 
//_____________________________________________________________________________
AliPIPE::AliPIPE(const char *name, const char *title)
       : AliModule(name,title)
{
  //
  // Standard constructor for beam pipe
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
