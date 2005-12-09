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

//____________________________________________________________________
// 
// Base class for FMD reconstruction algorithms. 
//
// Derived classes will implement various ways of reconstructing the
// charge particle multiplicity in the FMD.  
// 
//      +---------------------+       +---------------------+
//      | AliFMDReconstructor |<>-----| AliFMDMultAlgorithm |
//      +---------------------+       +---------------------+
//                                               ^
//                                               |
//                                   +-----------+---------+
//                                   |                     |
//                         +-------------------+   +------------------+
//                         | AliFMDMultPoisson |   | AliFMDMultNaiive |
//                         +-------------------+   +------------------+
//
// AliFMDReconstructor acts as a manager class.  It contains a list of
// AliFMDMultAlgorithm objects.  The call graph looks something like 
//
//
//       +----------------------+            +----------------------+
//       | :AliFMDReconstructor |            | :AliFMDMultAlgorithm |
//       +----------------------+            +----------------------+
//                  |                                  |
//    Reconstruct  +-+                                 |
//    ------------>| |                         PreRun +-+
//                 | |------------------------------->| |   
//                 | |                                +-+
//                 | |-----+ (for each event)          |
//                 | |     | *ProcessEvent             |
//                 |+-+    |                           |
//                 || |<---+                 PreEvent +-+
//                 || |------------------------------>| |      
//                 || |                               +-+
//                 || |-----+                          |
//                 || |     | ProcessDigits            |
//                 ||+-+    |                          |
//                 ||| |<---+                          |
//                 ||| |         *ProcessDigit(digit) +-+
//                 ||| |----------------------------->| |
//                 ||| |                              +-+
//                 ||+-+                               |
//                 || |                     PostEvent +-+
//                 || |------------------------------>| |
//                 || |                               +-+
//                 |+-+                                |
//                 | |                        PostRun +-+
//                 | |------------------------------->| |
//                 | |                                +-+
//                 +-+                                 |
//                  |                                  |
//
//
#include "AliFMDMultAlgorithm.h"	// ALIFMDMULTALGORITHM_H
#include "AliFMDDigit.h"		// ALIFMDDIGIT_H
#include <TClonesArray.h>               // ROOT_TClonesArray

//____________________________________________________________________
ClassImp(AliFMDMultAlgorithm)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDMultAlgorithm::AliFMDMultAlgorithm(const char* name, const char* title)
  : TNamed(name, title), 
    fTreeR(0), 
    fMult(0), 
    fFMD(0)
{
  // Default CTOR
}

//____________________________________________________________________
AliFMDMultAlgorithm::~AliFMDMultAlgorithm()
{
  // DTOR
  if (fMult) {
    fMult->Delete();
    delete fMult;
  }
}


//____________________________________________________________________
void
AliFMDMultAlgorithm::PreEvent(TTree* treeR, Float_t /* ipZ */) 
{
  // Executed before each event.
  if (fMult) fMult->Clear();
  fNMult = 0;
  fTreeR = treeR;
}



//____________________________________________________________________
// 
// EOF
//
