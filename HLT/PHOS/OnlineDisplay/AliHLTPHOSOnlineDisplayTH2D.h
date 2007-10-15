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

#define Z_ROWS 3
#define X_COLS 3
#define N_HISTOGRAMS  Z_ROWS*X_COLS +4

#include <TH2D.h>
#include <TCanvas.h>
#include "AliHLTPHOSBase.h"

class AliHLTPHOSOnlineDisplay;

class  AliHLTPHOSOnlineDisplayTH2D : public TH2D, public AliHLTPHOSBase
{
public:
  AliHLTPHOSOnlineDisplayTH2D(AliHLTPHOSOnlineDisplay  *onlineDisplayPtr, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);
  virtual ~AliHLTPHOSOnlineDisplayTH2D();
 
  virtual void ExecuteEvent(Int_t event, Int_t pz, Int_t px);


  int GetZBin(Int_t pz);
  int GetXBin(Int_t px);

  void SetGain(int gain);

private:
  AliHLTPHOSOnlineDisplayTH2D(); 
  AliHLTPHOSOnlineDisplay *fOnlineDisplayPtr;

  //  TH1D *fHistPtr;
  TCanvas  *fgRawDataCanvas;
  TH1D     *fgRawDataPlotsPtr[N_HISTOGRAMS]; 

  int fGain;

};

#endif
