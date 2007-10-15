/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSModuleXCoordinate.h"
#include "AliHLTPHOSConstants.h"
#include "Rtypes.h"
#include <iostream>

using namespace PhosHLTConst;
using namespace std;

AliHLTPHOSModuleXCoordinate::AliHLTPHOSModuleXCoordinate(int x) : fX(x)
{
  //  SetX(x);
}


AliHLTPHOSModuleXCoordinate::AliHLTPHOSModuleXCoordinate() : fX(0)
{
  //never to be called
}



AliHLTPHOSModuleXCoordinate::~AliHLTPHOSModuleXCoordinate()
{

}


void
AliHLTPHOSModuleXCoordinate::SetX(int x)
{
  if( (x >= 0)  && (x < N_XCOLUMNS_MOD) )
    {
      fX = x;
    }
  else
    {
      cout << "ERROR, x value out of range" << endl;
      cout << "Attemt to set x to " << x <<endl;
      cout << "Allowd arne is " << 0 << " to " << N_XCOLUMNS_MOD <<endl;
    }

}



int
AliHLTPHOSModuleXCoordinate::GetX()
{
  return fX;
}
