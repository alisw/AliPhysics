/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////
//
//  Base Class for Detector specific Merging/Digitization   
//                  
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

/*
$Log$
Revision 1.5  2002/10/14 14:57:32  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.3.6.1  2002/07/24 10:08:13  alibrary
Updating VirtualMC

Revision 1.4  2002/07/17 07:29:53  jchudoba
Add private method GetNInputStreams(). Do not use it, it's just a temporary fix the PHOS and EMCAL code.

Revision 1.3  2001/11/14 14:50:33  jchudoba
Pass custom name and title to the TTask base class

Revision 1.2  2001/10/04 15:56:34  jchudoba
TTask inheritance

Revision 1.1  2001/07/27 13:02:06  jchudoba
ABC for detector digits merging/digitization

*/

// system includes
#include <Riostream.h>

// ROOT includes

// AliROOT includes
#include "AliDigitizer.h"
#include "AliRunDigitizer.h"

ClassImp(AliDigitizer)

AliDigitizer::AliDigitizer(const Text_t* name, const Text_t* title)
  :TTask(name,title) 
{
//
// dummy default ctor with name and title
//
  fManager = 0;
}

AliDigitizer::AliDigitizer(AliRunDigitizer *manager, 
			   const Text_t* name, const Text_t* title)
  :TTask(name,title)
{
//
// ctor with name and title
//
  fManager = manager;
  manager->AddDigitizer(this);
}



AliDigitizer::~AliDigitizer() {;}
////////////////////////////////////////////////////////////////////////
Int_t AliDigitizer::GetNInputStreams() const
{
//
// return number of input streams
//
  Int_t nInputStreams = 0 ;
  if (fManager)
    nInputStreams = fManager->GetNinputs() ;
  return nInputStreams ; 
}
