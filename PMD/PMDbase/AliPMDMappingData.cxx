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
#include "AliPMDMappingData.h"


ClassImp(AliPMDMappingData)

AliPMDMappingData::AliPMDMappingData()
{
  // Default constructor
  Reset();
}
// ----------------------------------------------------------------- //
AliPMDMappingData::AliPMDMappingData(const char* name)
{
  //constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
  
}
// ----------------------------------------------------------------- //
AliPMDMappingData::AliPMDMappingData(const AliPMDMappingData& mapda) :
  TNamed(mapda)
{
  // copy constructor
  SetName(mapda.GetName());
  SetTitle(mapda.GetName());
  Reset();
  for(Int_t iddl = 0; iddl < kDdl; iddl++)
    {
      for(Int_t imod = 0; imod < 48; imod++)
	{
	  fBeginPatchBus[iddl][imod] = mapda.GetBeginPatchBus(iddl,imod);
	  fEndPatchBus[iddl][imod]   = mapda.GetEndPatchBus(iddl,imod);
	}
      for(Int_t ibus = 0; ibus < kBus; ibus++)
	{
	  fModuleNo[iddl][ibus]    = mapda.GetModuleNo(iddl,ibus);
	  fMcmperBus[iddl][ibus]   = mapda.GetMcmperBus(iddl,ibus);
	  fStartRowBus[iddl][ibus] = mapda.GetStartRowBus(iddl,ibus);
	  fEndRowBus[iddl][ibus]   = mapda.GetEndRowBus(iddl,ibus);
	  fStartColBus[iddl][ibus] = mapda.GetStartColBus(iddl,ibus);
	  fEndColBus[iddl][ibus]   = mapda.GetEndColBus(iddl,ibus);
	}
    }

}
// ----------------------------------------------------------------- //
AliPMDMappingData &AliPMDMappingData::operator =(const AliPMDMappingData& mapda)
{
  //asignment operator
  SetName(mapda.GetName());
  SetTitle(mapda.GetName());
  Reset();

  for(Int_t iddl = 0; iddl < kDdl; iddl++)
    {
      for(Int_t imod = 0; imod < 48; imod++)
	{
	  fBeginPatchBus[iddl][imod] = mapda.GetBeginPatchBus(iddl,imod);
	  fEndPatchBus[iddl][imod]   = mapda.GetEndPatchBus(iddl,imod);
	}
      for(Int_t ibus = 0; ibus < kBus; ibus++)
	{
	  fModuleNo[iddl][ibus]    = mapda.GetModuleNo(iddl,ibus);
	  fMcmperBus[iddl][ibus]   = mapda.GetMcmperBus(iddl,ibus);
	  fStartRowBus[iddl][ibus] = mapda.GetStartRowBus(iddl,ibus);
	  fEndRowBus[iddl][ibus]   = mapda.GetEndRowBus(iddl,ibus);
	  fStartColBus[iddl][ibus] = mapda.GetStartColBus(iddl,ibus);
	  fEndColBus[iddl][ibus]   = mapda.GetEndColBus(iddl,ibus);
	}
    }

  return *this;
}
// ----------------------------------------------------------------- //
AliPMDMappingData::~AliPMDMappingData()
{
  //destructor
}
// ----------------------------------------------------------------- //
void AliPMDMappingData::Reset()
{
  for(Int_t i = 0; i < 6; i++)
    {
      for(Int_t j = 0; j < 48; j++)
	{
	  fBeginPatchBus[i][j] = -1;
	  fEndPatchBus[i][j] = -1;
	}
    }

  for(Int_t iddl = 0; iddl < kDdl; iddl++)
    {
      for(Int_t ibus = 0; ibus < kBus; ibus++)
	{
	  fModuleNo[iddl][ibus]    = -1;
	  fMcmperBus[iddl][ibus]   = -1;
	  fStartRowBus[iddl][ibus] = -1;
	  fEndRowBus[iddl][ibus]   = -1;
	  fStartColBus[iddl][ibus] = -1;
	  fEndColBus[iddl][ibus]   = -1;
	}
    }

}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetBeginPatchBus(Int_t iddl, Int_t imod) const
{
  return fBeginPatchBus[iddl][imod];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetEndPatchBus(Int_t iddl, Int_t imod) const
{
  return fEndPatchBus[iddl][imod];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetModuleNo(Int_t iddl, Int_t ibus) const
{
  return fModuleNo[iddl][ibus];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetMcmperBus(Int_t iddl, Int_t ibus) const
{
  return fMcmperBus[iddl][ibus];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetStartRowBus(Int_t iddl, Int_t ibus) const
{
  return fStartRowBus[iddl][ibus];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetEndRowBus(Int_t iddl, Int_t ibus) const
{
  return fEndRowBus[iddl][ibus];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetStartColBus(Int_t iddl, Int_t ibus) const
{
  return fStartColBus[iddl][ibus];
}
// ----------------------------------------------------------------- //
Int_t AliPMDMappingData:: GetEndColBus(Int_t iddl, Int_t ibus) const
{
  return fEndColBus[iddl][ibus];
}
// ----------------------------------------------------------------- //
void AliPMDMappingData::SetPatchBus(Int_t iddl, Int_t imod,
				    Int_t bpatchbus, Int_t epatchbus)
{
  fBeginPatchBus[iddl][imod] = bpatchbus;
  fEndPatchBus[iddl][imod]   = epatchbus;
}
// ----------------------------------------------------------------- //
void AliPMDMappingData::SetModuleNo(Int_t iddl, Int_t ibus, Int_t modno)
{
  fModuleNo[iddl][ibus] = modno;
}
// ----------------------------------------------------------------- //
void AliPMDMappingData::SetMcmperBus(Int_t iddl, Int_t ibus, Int_t totmcm)
{
  fMcmperBus[iddl][ibus] = totmcm;
}

// ----------------------------------------------------------------- //
void AliPMDMappingData::SetRowBus(Int_t iddl,Int_t ibus,Int_t rows,Int_t rowe)
{
  fStartRowBus[iddl][ibus] = rows;
  fEndRowBus[iddl][ibus]   = rowe;
}

// ----------------------------------------------------------------- //
void AliPMDMappingData::SetColBus(Int_t iddl,Int_t ibus,Int_t cols,Int_t cole)
{
  fStartColBus[iddl][ibus] = cols;
  fEndColBus[iddl][ibus]   = cole;
}
// ----------------------------------------------------------------- //

void AliPMDMappingData::Print(Option_t *) const
{
  printf("\n ######Mapping File for each ddl and patchbus ####\n");


  for(Int_t iddl = 0; iddl < kDdl; iddl++)
    {
      for(Int_t imod = 0; imod < 48; imod++)
	{
	  printf("%d %d %d %d \n",iddl, imod, fBeginPatchBus[iddl][imod],
		 fEndPatchBus[iddl][imod]);
	}
      for(Int_t ibus = 0; ibus < kBus; ibus++)
	{
	  printf("%d %d %d %d %d %d %d %d\n",iddl, ibus,
		 fModuleNo[iddl][ibus], fMcmperBus[iddl][ibus],
		 fStartRowBus[iddl][ibus], fEndRowBus[iddl][ibus],
		 fStartColBus[iddl][ibus], fEndColBus[iddl][ibus]);
	}
      printf("\n");
    }

}
