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

 
#include "AliPMDddldata.h"

ClassImp(AliPMDddldata)

//------------------------------------------------------------
AliPMDddldata::AliPMDddldata():
  fDetector(-1),
  fSMN(-1),
  fModule(-1),
  fPatchBus(-1),
  fMCM(-1),
  fChannel(-1),
  fRow(-1),
  fCol(-1),
  fSignal(-1),
  fBit(-1)
{
  //
  // ctor
  //

}

//___________________________________________
AliPMDddldata::~AliPMDddldata()
{
  // 
  // dtor
  //
}

//___________________________________________
AliPMDddldata::AliPMDddldata(const AliPMDddldata & ddl) :
  TObject(),
  fDetector(ddl.fDetector),
  fSMN(ddl.fSMN),
  fModule(ddl.fModule),
  fPatchBus(ddl.fPatchBus),
  fMCM(ddl.fMCM),
  fChannel(ddl.fChannel),
  fRow(ddl.fRow),
  fCol(ddl.fCol),
  fSignal(ddl.fSignal),
  fBit(ddl.fBit)
{
  //
  // copy ctor
  //
}

//___________________________________________
AliPMDddldata& AliPMDddldata::operator=(const AliPMDddldata &ddl)
{
  // 
  // assignment operator
  //
  if (this != &ddl)
    {
      fDetector = ddl.fDetector;
      fSMN   = ddl.fSMN;
      fModule   = ddl.fModule;
      fPatchBus = ddl.fPatchBus;
      fMCM      = ddl.fMCM;
      fChannel  = ddl.fChannel;
      fRow      = ddl.fRow;
      fCol      = ddl.fCol;
      fSignal   = ddl.fSignal;
      fBit      = ddl.fBit;
    }
  return *this;
}
