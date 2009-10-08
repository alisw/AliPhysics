// $Id: AliHLTPHOSRawAnalyzerChiSquareFit.cxx 30036 2008-11-24 16:43:44Z odjuvsla $


/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
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


#include "AliHLTCaloRawAnalyzerChiSquareFit.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliHLTCaloRawAnalyzerChiSquareFit) 

AliHLTCaloRawAnalyzerChiSquareFit::AliHLTCaloRawAnalyzerChiSquareFit(const AliHLTCaloRawAnalyzerChiSquareFit&):AliHLTCaloRawAnalyzer()
{
  //comment
}

AliHLTCaloRawAnalyzerChiSquareFit::AliHLTCaloRawAnalyzerChiSquareFit():AliHLTCaloRawAnalyzer()
{
  //comment
}


AliHLTCaloRawAnalyzerChiSquareFit::~AliHLTCaloRawAnalyzerChiSquareFit()
{

} //end AliHLTCaloRawAnalyzerChiSquareFit


void 
AliHLTCaloRawAnalyzerChiSquareFit::Evaluate(int /*start*/, int /*length*/)
{
  /*

  */

  //thats all for the moment 
} //end FitChiSquareFit





