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
//  Magnetic Dipole                                                          //
//  This class contains the description of the base class for the            //
//  magnetic dipole                                                          //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliDIPOClass.gif">
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
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDIPO.h"
#include "AliRun.h"
#include "AliLog.h"
 
ClassImp(AliDIPO)
 
//_____________________________________________________________________________
AliDIPO::AliDIPO()
{
  //
  // Default constructor for the Dipole Magnet
  //
}
 
//_____________________________________________________________________________
AliDIPO::AliDIPO(const char *name, const char *title)
  : AliModule(name,title)
{
  //
  // Standard constructor for the Magnetic Dipole
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliDIPO::Init()
{
  //
  // Initialise Magnetic Dipole after it has been built
  //
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" DIPO_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the ABSO initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}
 
