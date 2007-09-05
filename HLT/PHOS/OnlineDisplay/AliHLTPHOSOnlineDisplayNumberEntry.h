#ifndef ALIHLTPHOSONLINEDISPLAYNUMBERENTRY_H
#define ALIHLTPHOSONLINEDISPLAYNUMBERENTRY_H 

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE DCS Project.  *
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

#include <TGNumberEntry.h>


class AliHLTPHOSOnlineDisplayNumberEntry : public TGNumberEntry
{
 public:
  AliHLTPHOSOnlineDisplayNumberEntry();
  ~AliHLTPHOSOnlineDisplayNumberEntry();
  AliHLTPHOSOnlineDisplayNumberEntry(const TGWindow* parent = 0, Double_t val = 0, Int_t digitwidth = 5, 
		  Int_t id = -1, TGNumberFormat::EStyle style = kNESReal, 
		  TGNumberFormat::EAttribute attr = kNEAAnyNumber, 
		  TGNumberFormat::ELimit limits = kNELNoLimits, Double_t min = 0, Double_t max = 1);

  virtual void ValueChanged(Long_t t);

  virtual void ValueSet(Long_t t);

  void SetLimits(int low, int high);
  void SetButtonType(char c);

 private:
  char buttonType;
  int lowLimit;
  int highLimit;
};


#endif  
