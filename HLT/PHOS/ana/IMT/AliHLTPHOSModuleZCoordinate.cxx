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
#include "AliHLTPHOSModuleZCoordinate.h"
#include "AliHLTPHOSConstants.h"
#include "Rtypes.h"
#include <iostream>

//  const unsigned int N_ZROWS_MOD      =  56;                    /**<Number of rows per module*/       
//  const unsigned int N_XCOLUMNS_MOD   =  64;  

using namespace PhosHLTConst;
using namespace std;


AliHLTPHOSModuleZCoordinate::AliHLTPHOSModuleZCoordinate(int z) : fZ(0)
{
  SetZ(z);
}


AliHLTPHOSModuleZCoordinate::AliHLTPHOSModuleZCoordinate() : fZ(0)
{
  //never to be called
}

AliHLTPHOSModuleZCoordinate::~AliHLTPHOSModuleZCoordinate()
{

}

void
AliHLTPHOSModuleZCoordinate::SetZ(int z)
{
  if( (z >= 0)  && (z < N_ZROWS_MOD) )
    {
      fZ = z;
    }
  else
    {
      cout << "ERROR, x value out of range" << endl;
      cout << "Attemt to set z to " << z <<endl;
      cout << "Allowd arne is " << 0 << " to " << N_ZROWS_MOD <<endl;
    }

}

int
AliHLTPHOSModuleZCoordinate::GetZ()
{
  return fZ;
}
