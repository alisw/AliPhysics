//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTest.h
    @author Matthias Richter
    @date   
    @brief  Abstract interface for the AliHLTTestProcessor test methods
 */

class AliHLTTest
{
public:
  AliHLTTest() {};
  virtual ~AliHLTTest() {};

  virtual bool CheckRunNo(unsigned runNo) const = 0;
  virtual bool CheckChainId(const char* chainId) const = 0;
  virtual bool CheckDataType(const char* id, const char* origin) const = 0;
  virtual bool CheckMagneticField(float bz) = 0;

};
