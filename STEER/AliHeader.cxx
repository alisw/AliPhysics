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
*/

#include "AliHeader.h"
#include "stdio.h"
 
ClassImp(AliHeader)

AliHeader::AliHeader()
{
  fRun=0;	
  fNvertex=0;
  fNprimary=0;
  fNtrack=0;
  fEvent=0;
}

AliHeader::AliHeader(Int_t run, Int_t event)
{
  fRun=run;	
  fNvertex=0;
  fNprimary=0;
  fNtrack=0;
  fEvent=event;
}

void AliHeader::Reset(Int_t run, Int_t event)
{
  fRun=run;	
  fNvertex=0;
  fNprimary=0;
  fNtrack=0;
  fEvent=event;
}

void AliHeader::Dump()
{
  printf(
"\n=========== Header for run %d Event %d = beginning ======================================\n",
  fRun,fEvent);
  printf("              Number of Vertex %d\n",fNvertex);
  printf("              Number of Primary %d\n",fNprimary);
  printf("              Number of Tracks %d\n",fNtrack);
  printf(
  "=========== Header for run %d Event %d = end ============================================\n\n",
  fRun,fEvent);
}
