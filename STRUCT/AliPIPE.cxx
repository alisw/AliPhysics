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
#include "TGeant3.h"
 
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
       : AliDetector(name,title)
{
  //
  // Standard constructor for beam pipe
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliPIPE::BuildGeometry()
{
  //
  // Create ROOT TNode geometry. Only for sensitive detectors
  //
}
 


//_____________________________________________________________________________
void AliPIPE::StepManager()
{
  //
  // Called at every step in the beam pipe
  //
}
