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


//Root includes

//AliRoot includes
#include "AliMUONLoader.h"
#include "AliMUONConstants.h"
#include "AliLog.h"

ClassImp(AliMUONLoader)
//___________________________________________________________________
AliMUONLoader::AliMUONLoader()
  : AliLoader(),
    fMUONData(0)
{
//default constructor
}
//_______________________________________________________________________________
AliMUONLoader::AliMUONLoader(const Char_t* detname,const Char_t* eventfoldername)
  : AliLoader(detname,eventfoldername),
    fMUONData(0)
{
}
//_______________________________________________________________________________
AliMUONLoader::AliMUONLoader(const Char_t * detname,TFolder* eventfolder)
  : AliLoader(detname,eventfolder),
  fMUONData(0)
{
//constructor
}
//___________________________________________
AliMUONLoader::AliMUONLoader(const AliMUONLoader& rhs)
  : AliLoader(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//_______________________________________________________________________________
AliMUONLoader::~AliMUONLoader()
{
//detructor 
}
//-------------------------------------------------------------------
AliMUONLoader&  
AliMUONLoader::operator=(const AliMUONLoader& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//_______________________________________________________________________________
void AliMUONLoader::SetMUONData(AliMUONData * MUONData)
{
  fMUONData = MUONData;
}
//_______________________________________________________________________________
AliMUONData * AliMUONLoader::GetMUONData()
{
  return fMUONData;
}

