//-*- Mode: C++ -*-
// $Id$

 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSHISTOGRAMPRODUCER_H
#define ALIHLTPHOSHISTOGRAMPRODUCER_H


/**
 * Class does 
 *
 * @file   AliHLTPHOSHistogramProducer.h
 * @author Oystein Djuvsland
 * @date
 * @brief  
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSBase.h"

class TH1D;
class TNtuple;
class AliHLTPHOSCaloClusterContainerStruct;

/** 
 * @class AliHLTPHOSHistogramProducer
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSHistogramProducer : public AliHLTPHOSBase
{
  
public:
  
  AliHLTPHOSHistogramProducer();
  ~AliHLTPHOSHistogramProducer();
  

  AliHLTPHOSHistogramProducer(const AliHLTPHOSHistogramProducer &) :
    AliHLTPHOSBase(),
    fClusterEnergiesHistPtr(0),
    fMultiplicitiesHistPtr(0),
    fClusterNtuplePtr(0),
    fFillClusterEnergies(false),
    fFillMultiplicities(false),
    fFillNtuple(false),
    fMaxNtupleEntries(1000000000)
  {
    //comment
  }
  
  AliHLTPHOSHistogramProducer & operator = (const AliHLTPHOSHistogramProducer)
  {
    //Assignment
    return *this;
  }

  Int_t Fill(AliHLTPHOSCaloClusterContainerStruct* clusterContainerPtr);

  Int_t InitializeObjects();

  TH1D* GetClusterEnergiesHistogram() { return fClusterEnergiesHistPtr; }
  TH1D* GetMultiplicitiesHistogram() { return fMultiplicitiesHistPtr; }
  TNtuple* GetClusterNtuple() { return fClusterNtuplePtr; }

  void SetFillClusterEnergies(bool val) { fFillClusterEnergies = val; }
  void SetFillMultiplicities(bool val) { fFillMultiplicities = val; }
  void SetFillClusterNtuple(bool val) { fFillNtuple = val; }
  void SetMaxNtupleEntries(Int_t n) { fMaxNtupleEntries = n; }
  
private:

  TH1D* fClusterEnergiesHistPtr;
  TH1D* fMultiplicitiesHistPtr;
  TNtuple* fClusterNtuplePtr;

  bool fFillClusterEnergies;
  bool fFillMultiplicities;
  bool fFillNtuple;
  
  Int_t fMaxNtupleEntries;

};
#endif

