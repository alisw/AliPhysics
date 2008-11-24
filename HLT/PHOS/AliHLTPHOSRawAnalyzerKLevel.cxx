// $Id$


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

#include "AliHLTPHOSRawAnalyzerKLevel.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliHLTPHOSRawAnalyzerKLevel) 

//________________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerKLevel::AliHLTPHOSRawAnalyzerKLevel(const AliHLTPHOSRawAnalyzerKLevel&):AliHLTPHOSRawAnalyzer(), fTKLevel(0)
{
  //comment  
}


//________________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerKLevel::AliHLTPHOSRawAnalyzerKLevel():AliHLTPHOSRawAnalyzer(), fTKLevel(0) 
{
  //comment
  cout <<"You cannot invoke the Fitter without arguments"<<endl;;
}


//________________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerKLevel::~AliHLTPHOSRawAnalyzerKLevel()
{
  //comment
} //end AliHLTPHOSRawAnalyzerKLevel



//________________________________________________________________________________________________________
void 
AliHLTPHOSRawAnalyzerKLevel::Evaluate(int /*start*/, int /*length*/)
{
  //thats all 
} //end FitKLevel




