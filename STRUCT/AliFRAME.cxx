///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Space Frame                                                              //
//  This class contains the description of the Space Frame geometry          //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliFRAMEClass.gif">
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
  : AliModule(name,title)
{
  //
  // Standard constructor
  //
    
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
