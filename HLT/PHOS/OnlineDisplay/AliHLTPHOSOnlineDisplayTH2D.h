//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSONLINEDISPLAYTH2D_H
#define ALIHLTPHOSONLINEDISPLAYTH2D_H

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

#define ZROWS 3
#define XCOLS 3
#define NHISTOGRAMS  ZROWS*XCOLS +4

#include <TH2D.h>
#include <TCanvas.h>

#include <iostream>

//#include "AliHLTPHOSBase.h"

#include "AliHLTPHOSConstants.h" 

using namespace PhosHLTConst;
using namespace std;

class AliHLTPHOSOnlineDisplay;

//class  AliHLTPHOSOnlineDisplayTH2D : public TH2D, public AliHLTPHOSBase
class  AliHLTPHOSOnlineDisplayTH2D : public TH2D
{
public:
  AliHLTPHOSOnlineDisplayTH2D(AliHLTPHOSOnlineDisplay  *onlineDisplayPtr, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, 
			      Int_t nbinsy, Double_t ylow, Double_t yup);

  virtual ~AliHLTPHOSOnlineDisplayTH2D();
  virtual void ExecuteEvent(Int_t event, Int_t pz, Int_t px);
  
  void SetRunNumber(const int runnumber) 
  { 
    fRunNumber = runnumber;  
    cout << __FILE__ <<":"<< __LINE__ << "RunNumber was set to "<< fRunNumber  <<endl; ;
  };

  int GetZBin(Int_t pz);
  int GetXBin(Int_t px);

  void EvaluateBinPosition(const char *info, int *z, int *x);


private:
  AliHLTPHOSOnlineDisplayTH2D(); 
  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr;
  TCanvas  *fgRawDataCanvasPtr[NGAINS];
  TH1D     *fgRawDataPlotsPtr[NHISTOGRAMS][NGAINS]; 
  TCanvas  *fgRawDataCanvasSinglePtr[NGAINS];
  TH1D     *fgRawDataPlotsSinglePtr[NGAINS]; 
  

  int fRunNumber;
  //  bool fIsSetRunNumber;

};

#endif
