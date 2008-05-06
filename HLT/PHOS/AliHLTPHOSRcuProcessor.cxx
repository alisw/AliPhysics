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
#include "AliHLTPHOSRcuProcessor.h"

AliHLTPHOSRcuProcessor::AliHLTPHOSRcuProcessor():AliHLTPHOSProcessor(), 
						 AliHLTPHOSRcuProperties()
			// 			 fkEquippmentID(0),
// 						 fRcuX(0), 
// 						 fRcuZ(0),
// 						 fRcuZOffset(0),
// 						 fRcuXOffset(0),
// 						 fIsSetEquippmentID(kFALSE)    
{
  //  cout << "AliHLTPHOSRcuProcessor::AliHLTPHOSRcuProcessor() "<< endl;
}

AliHLTPHOSRcuProcessor::~AliHLTPHOSRcuProcessor()
{

}


/*
const AliHLTUInt16_t
AliHLTPHOSRcuProcessor::GetEquippmentID() const
{
  return fkEquippmentID;
}
*/


 /*
void 
AliHLTPHOSRcuProcessor::SetEquippmentID(AliHLTUInt16_t id)
{
  AliHLTUInt16_t  &ref = const_cast<AliHLTUInt16_t&>(fkEquippmentID); 
  ref = id;
}
*/



 
  //void 
  //AliHLTPHOSRcuProcessor::SetCoordinates(AliHLTUInt16_t /*equippmentID*/)
/*
{
  int rcuIndex =  (fkEquippmentID - 1792)%N_RCUS_PER_MODULE;
  fModID = (fkEquippmentID  -1792 -rcuIndex)/N_RCUS_PER_MODULE;
  
  if(rcuIndex == 0)
    {
      fRcuX = 0; 
      fRcuZ = 0;
    }

  if(rcuIndex == 1)
    {
      fRcuX = 0; 
      fRcuZ = 1;
    }
 
  if(rcuIndex == 2)
    {
      fRcuX = 1; 
      fRcuZ = 0;
    }

  if(rcuIndex == 3)
    {
      fRcuX = 1; 
      fRcuZ = 1;
    }

  fRcuZOffset =  N_ZROWS_RCU*fRcuZ;
  fRcuXOffset =  N_XCOLUMNS_RCU*fRcuX;

*/

//   cout <<"********InitInfo************"<< endl;
//   cout <<"AliHLTPHOSRawAnalyzerComponent::SetCoordinate casted"<< endl;
//   cout <<"Equpippment ID =\t"<< fkEquippmentID <<endl;
//   cout <<"Mod ID =\t"<<  (int)fModID<<endl;
//   cout <<"RCUX =\t\t" << (int)fRcuX << endl;
//   cout <<"RCUZ =\t\t" << (int)fRcuZ << endl;
//   cout <<"RcuZOffset = \t" <<  (int)fRcuZOffset << endl;
//   cout <<"RcuXOffset = \t" <<  (int)fRcuXOffset << endl << endl;

//}


