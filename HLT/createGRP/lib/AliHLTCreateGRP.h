#ifndef ALI_HLT_CREATE_GRP_H
#define ALI_HLT_CREATE_GRP_H

// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: David Rohr <drohr@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include <TObject.h>

class AliDCSSensor;

class AliHLTCreateGRP : public TObject
{
public:
	int CreateGRP(Int_t runNumber, TString detectorList, TString beamType, TString runType, UInt_t startShift, UInt_t endShift, Bool_t defaults);
	
private:
	UInt_t createDetectorMask(TObjArray* listOfDetectors);
	template <typename T> AliDCSSensor* CreateSensor(const char* id, T value, UInt_t starttime, UInt_t endtime);
	
	ClassDef(AliHLTCreateGRP, 0);
};
#endif 
