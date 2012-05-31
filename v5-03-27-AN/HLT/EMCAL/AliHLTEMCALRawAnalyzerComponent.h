#ifndef ALIHLTEMCALRAWANALYZERCOMPONENT_H
#define ALIHLTEMCALRAWANALYZERCOMPONENT_H

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

// Base class fro anlyzing EMCAL raww data
// Further documentation found in base class
// --------------
// --------------


#include "AliHLTCaloRawAnalyzerComponentv3.h"
 

class  AliHLTEMCALRawAnalyzerComponent : public AliHLTCaloRawAnalyzerComponentv3
{
 public:
  AliHLTEMCALRawAnalyzerComponent( fitAlgorithm algo );
  virtual ~AliHLTEMCALRawAnalyzerComponent();
  void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  virtual const char* GetComponentID() = 0;
  virtual AliHLTComponent* Spawn() = 0; 
  
 private:
  AliHLTEMCALRawAnalyzerComponent();
  AliHLTEMCALRawAnalyzerComponent(const AliHLTEMCALRawAnalyzerComponent & );
  AliHLTEMCALRawAnalyzerComponent & operator = (const AliHLTEMCALRawAnalyzerComponent  &);
  virtual void InitMapping( const int specification );
};

#endif
