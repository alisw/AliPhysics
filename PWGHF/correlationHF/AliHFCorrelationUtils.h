#ifndef ALIHFCORRELATIONUTILS_H
#define ALIHFCORRELATIONUTILS_H

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

/////////////////////////////////////////////////////////////
//
// Collections of methods used in several points of correlation analyses
//
// Author: A. Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>

#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>

#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>


class AliHFCorrelationUtils : public TObject {

 public:
  
  enum ETypeOfMeson{kDzero=0,kDstar=1,kDplus=2,kDaverage=4};
  enum ETypeOfSystem{kpp=0,kpPb=1};

  static TH1D* ReflectHisto(TH1D *h,Double_t scale);
  static TH1D* DuplicateHistoTo2piRange(TH1D *h,Double_t scale);
  static void GetMCClosureModulation(Double_t ptD, Double_t ptTrmin, Double_t ptTrmax, Double_t mod[]);
  
 private:
  
  ClassDef(AliHFCorrelationUtils,1);

};

#endif

