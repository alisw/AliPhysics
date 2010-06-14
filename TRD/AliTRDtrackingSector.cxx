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

/* $Id: AliTRDtrackingSector.cxx 23810 2008-02-08 09:00:27Z hristov $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Tracking data container for one sector                                   //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDtrackingSector.h"
#include "AliTRDtrackingChamber.h"


ClassImp(AliTRDtrackingSector)

//_____________________________________________________________________________
AliTRDtrackingSector::AliTRDtrackingSector()
  :fSector(-1)
  ,fN(0)
  ,fGeom(NULL)
{
  // Default constructor
  
  memset(fChamber, 0, AliTRDgeometry::kNdets*sizeof(AliTRDtrackingChamber*));
  memset(fIndex, -1, AliTRDgeometry::kNdets*sizeof(Char_t));
  memset(fX0, 0, AliTRDgeometry::kNlayer*sizeof(Float_t));
}

//_____________________________________________________________________________
AliTRDtrackingSector::AliTRDtrackingSector(AliTRDgeometry *geo, Int_t gs)
  :fSector(gs)
  ,fN(0)
  ,fGeom(geo)
{
  //
  // AliTRDtrackingSector Constructor
  //

  memset(fChamber, 0, AliTRDgeometry::kNdets*sizeof(AliTRDtrackingChamber*));
  memset(fIndex, -1, AliTRDgeometry::kNdets*sizeof(Char_t));
  memset(fX0, 0, AliTRDgeometry::kNlayer*sizeof(Float_t));
}

    
//_____________________________________________________________________________
void AliTRDtrackingSector::Init(const AliTRDReconstructor *rec)
{		
// 	Steer building of tracking chambers and build tracking sector.
// 	Propagate radial position information (calibration/alignment aware) from chambers to sector level
//
  
  AliTRDchamberTimeBin *tb = NULL;
  AliTRDtrackingChamber **tc = &fChamber[0];
  for(Int_t ic = 0; (ic<AliTRDgeometry::kNdets) && (*tc); ic++, tc++){
    for(Int_t itb=0; itb<AliTRDseedV1::kNtb; itb++){
      if(!(tb = (*tc)->GetTB(itb))) continue;
      tb->SetReconstructor(rec);
    }
    (*tc)->Build(fGeom, rec->IsHLT());
  }
    
  Int_t nl;
  for(int il=0; il<AliTRDgeometry::kNlayer; il++){
    fX0[il] = 0.; nl = 0;
    for(int is=0; is<AliTRDgeometry::kNstack; is++){
      Int_t idx = is*AliTRDgeometry::kNlayer + il;
      if(fIndex[idx]<0) continue;
      fX0[il] += GetChamber(fIndex[idx])->GetX(); 
      nl++; 
    }
    if(!nl){
      //printf("Could not estimate radial position  of plane %d in sector %d.\n", ip, fSector);
      continue;
    }
    fX0[il] /= Float_t(nl);
  }
}



//_____________________________________________________________________________
void AliTRDtrackingSector::Clear(const Option_t *opt)
{
// Reset counters and steer chamber clear

  AliTRDtrackingChamber **tc = &fChamber[0];
  for(Int_t ich=0; ich<fN; ich++, tc++){ 
    (*tc)->Clear(opt);
    delete (*tc); (*tc) = NULL;   // I would avoid
  }	
  memset(fIndex, -1, AliTRDgeometry::kNdets*sizeof(Char_t));
  fN = 0;
}

//_____________________________________________________________________________
AliTRDtrackingChamber* AliTRDtrackingSector::GetChamber(Int_t stack, Int_t layer, Bool_t build)
{
// Return chamber at position (stack, plane) in current 
// sector or build a new one if it is not already created
  
  Int_t ch = stack*AliTRDgeometry::kNlayer + layer;
  if(fIndex[ch] >= 0) return fChamber[Int_t(fIndex[ch])];
  else if(!build) return NULL;
  
  // CHAMBER HAS TO BE BUILD
  Int_t rch = ch;do rch--; while(rch>=0 && fIndex[rch]<0);
  fIndex[ch] = rch >=0 ? fIndex[rch]+1 : 0; 
  fN++;
  
  memmove(&fChamber[Int_t(fIndex[ch])+1], &fChamber[Int_t(fIndex[ch])], (AliTRDgeometry::kNdets-fIndex[ch]-1)*sizeof(void*));
  for(Int_t ic = ch+1; ic<AliTRDgeometry::kNdets; ic++) fIndex[ic] += fIndex[ic] >= 0 ? 1 : 0;
  
  AliTRDtrackingChamber *chmb = fChamber[Int_t(fIndex[ch])] = new AliTRDtrackingChamber();
  chmb->SetDetector(AliTRDgeometry::GetDetector(layer, stack, fSector));
  return chmb;
}

//_____________________________________________________________________________
AliTRDtrackingChamber** AliTRDtrackingSector::GetStack(Int_t stack)
{
// Return chamber at position (stack, plane) in current 
// sector or build a new one if it is not already created
  
  if(stack<0 || stack>=AliTRDgeometry::kNstack) return NULL;
  
  Int_t ich, n = 0;
  for(int il=0; il<AliTRDgeometry::kNlayer; il++){
    ich = stack*AliTRDgeometry::kNlayer + il;
    if(fIndex[ich] < 0) fStack[il] = NULL; 
    else{
      fStack[il] = fChamber[Int_t(fIndex[ich])];
      n++;
    }
  }
  
  return n ? &fStack[0] : NULL;
}

//_____________________________________________________________________________
void AliTRDtrackingSector::Print(Option_t *opt) const
{
// Dump info about this tracking sector and the tracking chamber within
// 

  printf("\n\tSector[%2d]\n", fSector);
  for(int il=0; il<AliTRDgeometry::kNlayer; il++){
    for(int is =0; is<AliTRDgeometry::kNstack; is++){
      Int_t ch = is*AliTRDgeometry::kNlayer + il;
      if(opt) fChamber[Int_t(fIndex[ch])]->Print(opt);
      else printf("%2d[%2d] ", fIndex[ch], fIndex[ch]>=0 ? fChamber[Int_t(fIndex[ch])]->GetNClusters() : 0);
    }
    if(!opt) printf("\n");
  }

}
