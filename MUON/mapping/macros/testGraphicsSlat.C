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

// $Id$
// $MpId: $

//
// Test macro for drawing slat motifs with real contours
// Christian Finck, Subatech
//
#if !defined(__CINT__) || defined(__MAKECINT__)

// MUON includes
#include "AliMpSt345Reader.h"
#include "AliMpSlat.h"
#include "AliMpVPainter.h"
#include "AliMpMotifReader.h"
#include "AliMpMotifType.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotif.h"
#include "TVector2.h"
#include "TCanvas.h"
#endif

void testGraphicsMotif(Option_t* motifType = "R43", const TVector2& padSizes = TVector2(2.5,0.5))
{
  // Warning : this function leaks memory. But should be fine as only used 
  // interactively to check a few motifs at once...
  //
  AliMpMotifReader reader(kStation345,kBendingPlane);
  AliMpMotifType* type = reader.BuildMotifType(motifType);
  if (!type)
  {
    cerr << "Motif not found" << endl;
    return;
  }
  AliMpMotif* motif = new AliMpMotif(motifType,type,padSizes);
  AliMpMotifPosition* pos = new AliMpMotifPosition(0,motif,TVector2(0,0));
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(pos);
  if (!painter)
  {
    cerr << "Could not get a painter !" << endl;
    return;
  }
  TCanvas* c = new TCanvas();
  painter->Draw("MP");
}

//112230N
//112233NR3
//220000N
//122000NR1
//112200NR2
void testGraphicsSlat(AliMpPlaneType planeType = kBendingPlane, 
	              Option_t* option = "PMCI",
		      Bool_t saveJPG = false)
{
  // P plane
  // M motif
  // P pad
  // I indices

  Char_t *slatName[19] = {"122000SR1", "112200SR2", "122200S", "222000N", "220000N",
			  "122000NR1", "112200NR2", "122200N",
			  "122330N", "112233NR3", "112230N", "222330N", "223300N", "333000N", "330000N",
			  "112233N", "222333N", "223330N", "333300N"};
			  
  TCanvas *c1[19];
  Char_t c1Name[255];
  Char_t c1NameJpg[255];

  for (Int_t i = 0; i < 19; i++) {
    sprintf(c1Name, "%s%d", "c1", i);
    c1[i]= new TCanvas(c1Name,slatName[i],10,10,1200,800);     

    Char_t* slatType = slatName[i];
    AliMpSt345Reader* reader = new AliMpSt345Reader();
    AliMpSlat* slat = reader->ReadSlat(slatType, planeType);
    AliMpVPainter* painter = AliMpVPainter::CreatePainter(slat);
    painter->Draw(option);
  
    if (planeType == kNonBendingPlane)
      sprintf(c1NameJpg, "%s%s", slatName[i], "_NonBending.jpg");
    else 
      sprintf(c1NameJpg, "%s%s", slatName[i], "_Bending.jpg");
    
    if (saveJPG) c1[i]->Print(c1NameJpg);
  }
 
}
