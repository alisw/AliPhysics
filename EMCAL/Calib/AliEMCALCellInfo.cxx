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

/* History of cvs commits:
 *
* $Log$
* Revision 1.2  2007/10/16 14:36:39  pavlinov
* fixed code violation (almost)
*
* Revision 1.1  2007/09/12 11:19:24  pavlinov
* added pi0 calibration, linearity, shower profile
*
*/

// Aug 1, 2007 - cells information in one place
#include "AliEMCALCellInfo.h"
#include "AliEMCALGeometry.h"

#include <TObjArray.h>
 
ClassImp(AliEMCALCellIndexes) 

// ------------------------------------------------------------------------------
AliEMCALCellIndexes::AliEMCALCellIndexes() : fAbsId(-1), fNSupMod(-1), fNModule(-1), fNIphi(-1), 
fNIeta(-1), fIPhi(-1), fIEta(-1), fIPhim(-1), fIEtam(-1)   
{
}


ClassImp(AliEMCALCellInfo)

// ------------------------------------------------------------------------------
AliEMCALCellInfo::AliEMCALCellInfo() : TNamed("",""), fTable(0), fCurrentInd(0)
{
  //default constructor
}

// ------------------------------------------------------------------------------
AliEMCALCellInfo::AliEMCALCellInfo(const AliEMCALCellInfo& info) 
  : TNamed(info.GetName(),info.GetTitle()), 
    fTable(info.fTable), fCurrentInd(info.fCurrentInd)
{
  //copy constructor
}

// ------------------------------------------------------------------------------
AliEMCALCellInfo::AliEMCALCellInfo(const char* name, const Int_t nrow) : TNamed(name,"table of cell information") , fTable(0), fCurrentInd(0)
{
  fTable = new TObjArray(nrow);
}

// ------------------------------------------------------------------------------
void AliEMCALCellInfo::AddAt(AliEMCALCellIndexes* r)
{
  (*fTable)[fCurrentInd] = new AliEMCALCellIndexes(*r);
  fCurrentInd++;
}

// ------------------------------------------------------------------------------
AliEMCALCellInfo::~AliEMCALCellInfo()
{
  if(fTable) {
    fTable->Delete();
    delete fTable;
  }
}

// ------------------------------------------------------------------------------
AliEMCALCellIndexes* AliEMCALCellInfo::GetTable(Int_t i) const
{
  // Oct 16, 2007
  return (AliEMCALCellIndexes*)fTable->At(i);
}

// ------------------------------------------------------------------------------
void AliEMCALCellInfo::PrintTable(int ind1, int ind2) const
{
  // Oct 16, 2007
  printf(" %s : %s : #rows %i \n", GetName(), GetTitle(), fTable->GetSize());
  if(ind1==-1 && ind2==-1) return;
  printf(" fAbsId fNSupMod fNModule fNIphi fNIeta fIPhi fIEta fIPhim fIEtam\n");
  if(ind1 < 0) ind1 = 0;
  if(ind2 >= fTable->GetSize()) ind2 = fTable->GetSize();
  for(int i=ind1; i<ind2; i++) {
    AliEMCALCellIndexes* r = GetTable(i);
    if(r==0) break;
    printf(" %5.5i    %2.2i      %3.3i    %1.1i    %1.1i    %2.2i    %2.2i    %2.2i    %2.2i\n", 
    r->fAbsId, r->fNSupMod, r->fNModule, r->fNIphi, r->fNIeta, r->fIPhi, r->fIEta, r->fIPhim,r->fIEtam);
  }
}

// ------------------------------------------------------------------------------
AliEMCALCellInfo *AliEMCALCellInfo::GetTableForGeometry(const char* geoName)
{
  // Oct 16, 2007
  if(geoName==0) return 0;
  AliEMCALGeometry *g=AliEMCALGeometry::GetInstance(geoName);  
  return GetTableForGeometry(g);
}

// ------------------------------------------------------------------------------
AliEMCALCellInfo *AliEMCALCellInfo::GetTableForGeometry(AliEMCALGeometry *g)
{
  // Oct 16, 2007
  if(g==0) return 0;
  AliEMCALCellInfo *t = new AliEMCALCellInfo("CellInfo",g->GetNCells());

  for(Int_t absid=0; absid<g->GetNCells(); absid++){
    AliEMCALCellIndexes r;
    r.fAbsId   = absid;

    g->GetCellIndex(r.fAbsId,  r.fNSupMod, r.fNModule, r.fNIphi, r.fNIeta);
    g->GetCellPhiEtaIndexInSModule(r.fNSupMod,r.fNModule,r.fNIphi,r.fNIeta, r.fIPhi,r.fIEta);
    g->GetModulePhiEtaIndexInSModule(r.fNSupMod,r.fNModule, r.fIPhim, r.fIEtam); 

    t->AddAt(&r);
  }
  return t;
}
