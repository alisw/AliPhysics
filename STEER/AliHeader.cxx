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
Revision 1.4  2000/10/02 21:28:14  fca
Removal of useless dependecies via forward declarations

Revision 1.3  2000/07/12 08:56:25  fca
Coding convention correction and warning removal

Revision 1.2  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

#include "AliHeader.h"
#include <stdio.h>
 
ClassImp(AliHeader)

AliHeader::AliHeader()
{
  //
  // Default constructor
  //
  fRun=0;	
  fNvertex=0;
  fNprimary=0;
  fNtrack=0;
  fEvent=0;
}

AliHeader::AliHeader(Int_t run, Int_t event)
{
  //
  // Standard constructor
  //
  fRun=run;	
  fNvertex=0;
  fNprimary=0;
  fNtrack=0;
  fEvent=event;
}

void AliHeader::Reset(Int_t run, Int_t event)
{
  //
  // Resets the header with new run and event number
  //
  fRun=run;	
  fNvertex=0;
  fNprimary=0;
  fNtrack=0;
  fEvent=event;
}

void AliHeader::Print(const char* option)
{
  //
  // Dumps header content
  //
  printf(
"\n=========== Header for run %d Event %d = beginning ======================================\n",
  fRun,fEvent);
  printf("              Number of Vertex %d\n",fNvertex);
  printf("              Number of Primary %d\n",fNprimary);
  printf("              Number of Tracks %d\n",fNtrack);
  printf(
  "=========== Header for run %d Event %d = end ============================================\n\n",
  fRun,fEvent);
  
  // print  particle file map
  char* oMap = strstr(option,"Map");
  if (oMap) {
    printf("\nParticle file map: \n");
    for (Int_t i=0; i<fNtrack; i++) 
      printf("   %d th entry: %d \n",i,fParticleFileMap[i]);
  }    
}
