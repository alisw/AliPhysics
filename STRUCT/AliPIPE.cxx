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
//  Beam pipe class                                                          //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPIPEClass.gif">
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
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
