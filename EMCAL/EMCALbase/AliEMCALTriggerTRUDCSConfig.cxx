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

#include "AliEMCALTriggerTRUDCSConfig.h"
#include <bitset>
#include <iomanip>
#include <iostream>
#include <sstream>

ClassImp(AliEMCALTriggerTRUDCSConfig) ;

AliEMCALTriggerTRUDCSConfig::AliEMCALTriggerTRUDCSConfig() : TObject()
,fSELPF(0x1e1f)
,fL0SEL(0x1)
,fL0COSM(0)
,fGTHRL0(0)
,fRLBKSTU(0)
,fFw(0x21)
{
	for (Int_t i=0;i<6;i++) fMaskReg[i] = 0;
}

Int_t AliEMCALTriggerTRUDCSConfig::GetSegmentation()
{
	if (fL0SEL & 0x0001)
		return 2;
	else
		return 1;
}

bool AliEMCALTriggerTRUDCSConfig::operator==(const AliEMCALTriggerTRUDCSConfig &other) const {
	return (fSELPF == other.fSELPF) && (fL0SEL == other.fL0SEL) && (fL0COSM == other.fL0COSM)
					&& (fGTHRL0 == other.fGTHRL0) && (fRLBKSTU == other.fRLBKSTU) && (fFw == other.fFw)
					&& !memcmp(fMaskReg, other.fMaskReg, sizeof(UInt_t) * 6);
}

std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerTRUDCSConfig &conf){
	stream << "SELPF: " << std::hex << conf.fSELPF << ", L0SEL: " << conf.fL0SEL << ", L0COSM: " << std::dec 
				 << conf.fL0COSM << ", GTHRL0: " << conf.fGTHRL0 << ", RLBKSTU: " << conf.fRLBKSTU << ", FW: " << std::hex
				 << conf.fFw << std::dec << std::endl;
	for(int ireg = 0; ireg < 6; ireg++){
		stream << "Reg" << ireg << ": " << std::bitset<sizeof(UInt_t) *8>(conf.fMaskReg[ireg]) << " (" << conf.fMaskReg[ireg] << ")" << std::endl;
	}
	return stream;
}

std::string AliEMCALTriggerTRUDCSConfig::ToJSON() const {
	std::stringstream jsonstring;
	jsonstring << "{"
						 << "\"fSELPF\":" << fSELPF << ","
						 << "\"fL0SEL\":" << fL0SEL << ","
						 << "\"fL0COSM\":" << fL0COSM << ","
						 << "\"fGTHRL0\":" << fGTHRL0 << ","
						 << "\"fRLBKSTU\":" << fRLBKSTU << ","
						 << "\"fFw\":" << fFw << ","
						 << "\"fMaskReg\":[" << fMaskReg[0] << "," << fMaskReg[1] << "," << fMaskReg[2] << "," << fMaskReg[3] << "," << fMaskReg[4] << "," << fMaskReg[5] << "]"
						 << "}";
	return jsonstring.str();
}