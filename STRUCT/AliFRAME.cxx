/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

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
AliFRAME::AliFRAME():
    fRefVolumeId(0)
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliFRAME::AliFRAME(const char *name, const char *title)
    : AliModule(name,title),
      fRefVolumeId(0)
{
  //
  // Standard constructor
  //
    
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
