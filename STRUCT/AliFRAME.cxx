///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Space Frame                                                              //
//  This class contains the description of the Space Frame geometry          //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliFRAMEClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliFRAME.h"
#include "AliRun.h"
#include "AliMC.h"
 
ClassImp(AliFRAME)
 
//_____________________________________________________________________________
AliFRAME::AliFRAME()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliFRAME::AliFRAME(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor
  //
    
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliFRAME::BuildGeometry()
{
  //
  // ROOT TNode geometry is built only for sensitive detectors
  // and not for structural elements
  //
}
 
//_____________________________________________________________________________

void AliFRAME::StepManager()
{
  //
  // Called at every step in the space frame
  //
}
