// -*- mode: c++ -*-

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

#include "AliCaloRawAnalyzerFitter.h"
#include "TF1.h"
#include <iostream>

using std::cout;
using std::endl;


AliCaloRawAnalyzerFitter::AliCaloRawAnalyzerFitter(const char *name, const char *nameshort ) :AliCaloRawAnalyzer( name, nameshort), 
											      fkEulerSquared(7.389056098930650227),
											      fTf1(0),
											      fFixTau(true)
{
  for(int i=0; i < ALTROMAXSAMPLES; i++)
    {
      fXaxis[i] = i;
    }

  fTf1 = new TF1( "myformula", "[0]*((x - [1])/[2])^2*exp(-2*(x -[1])/[2])",  0, 30 ); 
 
  if (fFixTau) 
    {
      fTf1->FixParameter(2, fTau);
    }
  else 
    {
      fTf1->ReleaseParameter(2); // allow par. to vary
      fTf1->SetParameter(2, fTau);
    }
}


AliCaloRawAnalyzerFitter::~AliCaloRawAnalyzerFitter()
{
  delete fTf1;
}


void 
AliCaloRawAnalyzerFitter::PrintFitResult(const TF1 *f) const
{
  //shutting up the rule checker
  cout << endl;
  cout << __FILE__ << __LINE__ << "Using this samplerange we get" << endl;
  cout << __FILE__ << __LINE__ << "AMPLITUDE = " << f->GetParameter(0)/fkEulerSquared << ",.. !!!!" << endl;
  cout << __FILE__ << __LINE__ << "TOF = " << f->GetParameter(1) << ",.. !!!!" << endl;
  cout << __FILE__ << __LINE__ << "NDF = " << f->GetNDF() << ",.. !!!!" << endl;
  cout << endl << endl;
}


// Bool_t 
// AliCaloRawAnalyzerFitter::GetFixTau() const 
// { 
//   return fFixTau; 
// }; 
