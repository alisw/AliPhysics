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
  void SetNPresSamples(const int presamples);
  void SetNSamples(const int samples);
  void SetAltroZeroSupression(const bool isZerosupressed);
  void SetAltroBaselineSubtraction(const bool isAltroBaselineSubtraction);
  //  void SetSoftwareBaselineSubtraction(bool isSoftwareBaselineSubtraction);
  int  GetNPresSamples() const {return  fNPresamples;}; 
  int  GetNSamples() const {return  fNSamples;}; 
  bool GetIsAltroZroSupresses() const {return   fIsAltroZeroSupressed;}; 
  bool GetIsAltroBaselineSubtraction() const {return fIsAltroBaselineSubtraction;};
  void PrintAltroDefaultValues() const;

protected:
  //Altro Config
  int fNPresamples; //comment
  int fNSamples; //comment
  int fNTotalSamples; //comment
  bool fIsAltroZeroSupressed; //comment
  bool fIsAltroBaselineSubtraction; //comment
};

#endif
