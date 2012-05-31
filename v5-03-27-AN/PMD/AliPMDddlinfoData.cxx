/***************************************************************************
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
//
// Author : 
//
#include "TNamed.h"
#include "AliCDBEntry.h"
#include "AliPMDddlinfoData.h"


ClassImp(AliPMDddlinfoData)

AliPMDddlinfoData::AliPMDddlinfoData()
{
  // Default constructor
  Reset();
}
// ----------------------------------------------------------------- //
AliPMDddlinfoData::AliPMDddlinfoData(const char* name)
{
  //constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  
}
// ----------------------------------------------------------------- //
AliPMDddlinfoData::AliPMDddlinfoData(const AliPMDddlinfoData& ddlinfoda) :
  TNamed(ddlinfoda)
{
  // copy constructor
  SetName(ddlinfoda.GetName());
  SetTitle(ddlinfoda.GetName());
  Reset();

  for(Int_t iddl = 0; iddl < 6; iddl++)
    {
      fModules[iddl] = ddlinfoda.GetNoOfModulePerDdl(iddl);
      for(Int_t imod = 0; imod < 12; imod++)
	{
	  fModuleNo[iddl][imod] = ddlinfoda.GetModulesPerDdl(iddl,imod);
	}
    }
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartRowA[idet][ismn] = ddlinfoda.GetStartRowA(idet,ismn);
	  fStartRowB[idet][ismn] = ddlinfoda.GetStartRowB(idet,ismn);
	  fEndRowA[idet][ismn]   = ddlinfoda.GetEndRowA(idet,ismn);
	  fEndRowB[idet][ismn]   = ddlinfoda.GetEndRowB(idet,ismn);
	  fStartColA[idet][ismn] = ddlinfoda.GetStartColA(idet,ismn);
	  fStartColB[idet][ismn] = ddlinfoda.GetStartColB(idet,ismn);
	  fEndColA[idet][ismn]   = ddlinfoda.GetEndColA(idet,ismn);
	  fEndColB[idet][ismn]   = ddlinfoda.GetEndColB(idet,ismn);
	}
    }

}
// ----------------------------------------------------------------- //
AliPMDddlinfoData &AliPMDddlinfoData::operator =(const AliPMDddlinfoData& ddlinfoda)
{
  //asignment operator
  SetName(ddlinfoda.GetName());
  SetTitle(ddlinfoda.GetName());
  Reset();

  for(Int_t iddl = 0; iddl < 6; iddl++)
    {
      fModules[iddl] = ddlinfoda.GetNoOfModulePerDdl(iddl);
      for(Int_t imod = 0; imod < 12; imod++)
	{
	  fModuleNo[iddl][imod] = ddlinfoda.GetModulesPerDdl(iddl,imod);
	}
    }
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartRowA[idet][ismn] = ddlinfoda.GetStartRowA(idet,ismn);
	  fStartRowB[idet][ismn] = ddlinfoda.GetStartRowB(idet,ismn);
	  fEndRowA[idet][ismn]   = ddlinfoda.GetEndRowA(idet,ismn);
	  fEndRowB[idet][ismn]   = ddlinfoda.GetEndRowB(idet,ismn);
	  fStartColA[idet][ismn] = ddlinfoda.GetStartColA(idet,ismn);
	  fStartColB[idet][ismn] = ddlinfoda.GetStartColB(idet,ismn);
	  fEndColA[idet][ismn]   = ddlinfoda.GetEndColA(idet,ismn);
	  fEndColB[idet][ismn]   = ddlinfoda.GetEndColB(idet,ismn);
	}
    }

  return *this;
}
// ----------------------------------------------------------------- //
AliPMDddlinfoData::~AliPMDddlinfoData()
{
  //destructor
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData::Reset()
{

  for(Int_t iddl = 0; iddl < 6; iddl++)
    {
      fModules[iddl] = -1;
      for(Int_t imod = 0; imod < 12; imod++)
	{
	  fModuleNo[iddl][imod] = -1;
	}
    }
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartRowA[idet][ismn] = -1;
	  fStartRowB[idet][ismn] = -1;
	  fEndRowA[idet][ismn]   = -1;
	  fEndRowB[idet][ismn]   = -1;
	  fStartColA[idet][ismn] = -1;
	  fStartColB[idet][ismn] = -1;
	  fEndColA[idet][ismn]   = -1;
	  fEndColB[idet][ismn]   = -1;
	}
    }
  
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetNoOfModulePerDdl(Int_t iddl) const
{
  return fModules[iddl];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetModulesPerDdl(Int_t iddl, Int_t imod) const
{
  return fModuleNo[iddl][imod];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetStartRowA(Int_t idet, Int_t ismn) const
{
  return fStartRowA[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetStartRowB(Int_t idet, Int_t ismn) const
{
  return fStartRowB[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetEndRowA(Int_t idet, Int_t ismn) const
{
  return fEndRowA[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetEndRowB(Int_t idet, Int_t ismn) const
{
  return fEndRowB[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetStartColA(Int_t idet, Int_t ismn) const
{
  return fStartColA[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetStartColB(Int_t idet, Int_t ismn) const
{
  return fStartColB[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetEndColA(Int_t idet, Int_t ismn) const
{
  return fEndColA[idet][ismn];
}
// ----------------------------------------------------------------- //
Int_t AliPMDddlinfoData:: GetEndColB(Int_t idet, Int_t ismn) const
{
  return fEndColB[idet][ismn];
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetNoOfModulePerDdl(Int_t iddl, Int_t nmod)
{
  fModules[iddl] = nmod;
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetModuleNoPerDdl(Int_t iddl, Int_t mod[])
{
  for(Int_t i = 0; i < 12; i++) fModuleNo[iddl][i] = mod[i];
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetStartRowA(Int_t srowa[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartRowA[idet][ismn] = srowa[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetStartRowB(Int_t srowb[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartRowB[idet][ismn] = srowb[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetEndRowA(Int_t erowa[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fEndRowA[idet][ismn] = erowa[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetEndRowB(Int_t erowb[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fEndRowB[idet][ismn] = erowb[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetStartColA(Int_t scola[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartColA[idet][ismn] = scola[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetStartColB(Int_t scolb[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fStartColB[idet][ismn] = scolb[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetEndColA(Int_t ecola[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fEndColA[idet][ismn] = ecola[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //
void AliPMDddlinfoData:: SetEndColB(Int_t ecolb[][24])
{
  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  fEndColB[idet][ismn] = ecolb[idet][ismn];
	}
    }
}
// ----------------------------------------------------------------- //

void AliPMDddlinfoData::Print(Option_t *) const
{
  printf("\n ######ddlinfo File for each ddl and patchbus ####\n");

  for(Int_t iddl=0; iddl<6; iddl++)
    {
      printf("%d %d \n",iddl, fModules[iddl]);
      for(Int_t imod = 0; imod < 12; imod++)
	{
	  printf("%d \n",fModuleNo[iddl][imod]);
	}
    }

  for(Int_t idet = 0; idet < 2; idet++)
    {
      for(Int_t ismn = 0; ismn < 24; ismn++)
	{
	  printf("%d %d %d %d %d %d %d %d %d %d \n",idet, ismn,
		 fStartRowA[idet][ismn], fEndRowA[idet][ismn],
		 fStartColA[idet][ismn], fEndColA[idet][ismn],
		 fStartRowB[idet][ismn], fEndRowB[idet][ismn],
		 fStartColB[idet][ismn], fEndColB[idet][ismn]);
	}
      printf("\n");
    }

}
