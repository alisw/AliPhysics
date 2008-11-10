//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSALTROCONFIG_H
#define ALIHLTPHOSALTROCONFIG_H

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



class  AliHLTPHOSAltroConfig
{
public:
  AliHLTPHOSAltroConfig();
  virtual ~AliHLTPHOSAltroConfig();
  void SetNPresSamples(int presamples);
  void SetNSamples(int samples);
  void SetAltroZeroSupression(bool isZerosupressed);
  void SetAltroBaselineSubtraction(bool isAltroBaselineSubtraction);
  //  void SetSoftwareBaselineSubtraction(bool isSoftwareBaselineSubtraction);
  inline int  GetNPresSamples(){return  fNPresamples;}; 
  inline int  GetNSamples(){return  fNSamples;}; 
  inline bool GetIsAltroZroSupresses(){return   fIsAltroZeroSupressed;}; 
  inline bool GetIsAltroBaselineSubtraction(){return fIsAltroBaselineSubtraction;};
  void PrintAltroDefaultValues();

protected:
  //Altro Config
  int fNPresamples;
  int fNSamples;
  int fNTotalSamples;
  bool fIsAltroZeroSupressed;
  bool fIsAltroBaselineSubtraction;
};

#endif
