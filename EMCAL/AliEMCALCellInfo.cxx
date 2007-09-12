/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Log$ */

// Aug 1, 2007 - cells information in one place
#include "AliEMCALCellInfo.h"
#include "AliEMCALGeometry.h"

#include <TObjArray.h>
 
ClassImp(cellInfo) 

cellInfo::cellInfo() : absId(-1), nSupMod(-1), nModule(-1), nIphi(-1), 
nIeta(-1), iPhi(-1), iEta(-1), iPhim(-1), iEtam(-1)   
{
}

ClassImp(AliEMCALCellInfo)
// ------------------------------------------------------------------------------

AliEMCALCellInfo::AliEMCALCellInfo() : TNamed("",""), fTable(0), fCurrentInd(0)
{
}

AliEMCALCellInfo::AliEMCALCellInfo(const char* name, const Int_t nrow) : TNamed(name,"table of cell information") , fTable(0), fCurrentInd(0)
{
  fTable = new TObjArray(nrow);
}

void AliEMCALCellInfo::AddAt(cellInfo* r)
{
  (*fTable)[fCurrentInd] = new cellInfo(*r);
  fCurrentInd++;
}

AliEMCALCellInfo::~AliEMCALCellInfo()
{
  if(fTable) {
    fTable->Delete();
    delete fTable;
  }
}

cellInfo* AliEMCALCellInfo::GetTable(Int_t i) const
{
  return (cellInfo*)fTable->At(i);
}

void AliEMCALCellInfo::PrintTable(int ind1, int ind2) const
{
  printf(" %s : %s : #rows %i \n", GetName(), GetTitle(), fTable->GetSize());
  if(ind1==-1 && ind2==-1) return;
  printf(" absId nSupMod nModule nIphi nIeta iPhi iEta iPhim iEtam\n");
  if(ind1 < 0) ind1 = 0;
  if(ind2 >= fTable->GetSize()) ind2 = fTable->GetSize();
  for(int i=ind1; i<ind2; i++) {
    cellInfo* r = GetTable(i);
    if(r==0) break;
    printf(" %5.5i    %2.2i      %3.3i    %1.1i    %1.1i    %2.2i    %2.2i    %2.2i    %2.2i\n", 
    r->absId, r->nSupMod, r->nModule, r->nIphi, r->nIeta, r->iPhi, r->iEta, r->iPhim,r->iEtam);
  }
}


AliEMCALCellInfo *AliEMCALCellInfo::GetTableForGeometry(const char* geoName)
{
  if(geoName==0) return 0;
  AliEMCALGeometry *g=AliEMCALGeometry::GetInstance(geoName);  
  return GetTableForGeometry(g);
}

AliEMCALCellInfo *AliEMCALCellInfo::GetTableForGeometry(AliEMCALGeometry *g)
{
  if(g==0) return 0;
  AliEMCALCellInfo *t = new AliEMCALCellInfo("CellInfo",g->GetNCells());

  for(Int_t absid=0; absid<g->GetNCells(); absid++){
    cellInfo r;
    r.absId   = absid;

    g->GetCellIndex(r.absId,  r.nSupMod, r.nModule, r.nIphi, r.nIeta);
    g->GetCellPhiEtaIndexInSModule(r.nSupMod,r.nModule,r.nIphi,r.nIeta, r.iPhi,r.iEta);
    g->GetModulePhiEtaIndexInSModule(r.nSupMod,r.nModule, r.iPhim, r.iEtam); 

    t->AddAt(&r);
  }
  return t;
}
