#ifndef ALIHLTEMCALRAWANALYZERLMSCOMPONENT_H
#define ALIHLTEMCALRAWANALYZERLMSCOMPONENT_H

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

//Evaluation of Amplitude and Peak position using Least Mean Square (LMS) fit
//-----------
//-----------
//-----------
//-----------

#include "AliHLTEMCALRawAnalyzerComponent.h"

class  AliCaloRawAnalyzerLMS;

//      AliHLTEMCALRawAnalyzerLMSComponent.h:
class  AliHLTEMCALRawAnalyzerLMSComponent : public AliHLTEMCALRawAnalyzerComponent
{

 public:
  AliHLTEMCALRawAnalyzerLMSComponent();
  virtual ~AliHLTEMCALRawAnalyzerLMSComponent();
  virtual const char* GetComponentID() ;
  virtual AliHLTComponent* Spawn();
  virtual int Deinit();

 private:
  AliHLTEMCALRawAnalyzerLMSComponent              (   const AliHLTEMCALRawAnalyzerLMSComponent &,  TString det ); 
  AliHLTEMCALRawAnalyzerLMSComponent & operator = (   const AliHLTEMCALRawAnalyzerLMSComponent );
  //  AliEMCALRawAnalyzerLMS *fRawAnalyzer;

};

#endif
