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

/*
$Log$
Revision 1.2  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/
#include <stdlib.h>

#include "AliMC.h"

ClassImp(AliMC)

AliMC* AliMC::fgMC=0;

AliMC* gMC;

AliMC::AliMC(const char *name, const char *title) : TNamed(name,title)
{
  if(fgMC) {
    printf("Cannot initialise twice MonteCarlo class\n");
    exit(1);
  } else {
    fgMC=this;
    gMC=this;
  }
}

