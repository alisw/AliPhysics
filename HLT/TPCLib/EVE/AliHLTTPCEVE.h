//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCEVE_H
#define ALIHLTTPCEVE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCEVE.h
    @author Matthias Richter
    @date   2008-11-22
    @brief  AliEVE bindings for the HLT TPC.
*/

#include "AliHLTLogging.h"

class TEveManager;
class TEvePointSet;
class TEveElement;
class AliRawReader;
class AliHLTOUT;
class AliHLTTPCClusterData;
class AliHLTTPCSpacePointContainer;

/**
 * @class AliHLTTPCEVE
 * AliEVE bindings for the HLT TPC.
 * The class provides methods to convert HLT TPC proprietary data structures
 * to EVE objects.
 * 
 * @ingroup alihlt_tpc
 */
class AliHLTTPCEVE : public AliHLTLogging {
 public:
  /** standard constructor */
  AliHLTTPCEVE();
  /** standard destructor */
  virtual ~AliHLTTPCEVE();

  /**
   * Create EVE point collection from HLTOUT of HLT.Digits.root
   * Calls ::MakePointSetFromHLTOUT(AliHLTOUT*, TEveElement*, Float_t)
   * @param path        path of the digit file
   * @param eventNo     event number to be displayed
   * @param cont        EVE element collection
   * @param maxR        geometrical cut, maximum radius of clusters
   */
  TEvePointSet* MakePointSetFromHLTDigits(const char* path, int eventNo, TEveElement* cont=0, Float_t maxR=270) const;

  /**
   * Create EVE point collection from HLTOUT from RawReader.
   * Calls ::MakePointSetFromHLTOUT(AliHLTOUT*, TEveElement*, Float_t)
   * @param pRawReader  the RawReader
   * @param cont        EVE element collection
   * @param maxR        geometrical cut, maximum radius of clusters
   */
  TEvePointSet* MakePointSetFromHLTOUT(AliRawReader* pRawReader, TEveElement* cont=0, Float_t maxR=270) const;

  /**
   * Create EVE point collection from HLTOUT instance.
   * Base method to create EVE point list from HLTOUT.
   * @param pHLTOUT     the HLTOUT instance
   * @param cont        EVE element collection
   * @param maxR        geometrical cut, maximum radius of clusters
   */
  TEvePointSet* MakePointSetFromHLTOUT(AliHLTOUT* pHLTOUT, TEveElement* cont=0, Float_t maxR=270) const;

  /**
   * Add clusters from AliHLTTPCClusterData set to EVE point collection.
   * @param clusters    EVE point collection
   * @param data        AliHLTTPCClusterData set
   * @param sizeInByte  size of the buffer of the data set in byte
   * @param slice       TPC slice number, determines ratotion of space points
   * @param maxR        geometrical cut, maximum radius of clusters
   */
  int AddClusters(TEvePointSet* clusters, const AliHLTTPCClusterData* data, unsigned int sizeInByte, int slice, Float_t maxR) const;

  int AddClusters(TEvePointSet* clusters, const AliHLTTPCSpacePointContainer* points, Float_t maxR=5000) const;

  int AddPointSet(TEveManager* pEve, const AliHLTTPCSpacePointContainer* points, Float_t maxR=5000, const char* title=NULL) const;

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTTPCEVE(const AliHLTTPCEVE&);
  /** assignment operator prohibited */
  AliHLTTPCEVE& operator=(const AliHLTTPCEVE&);

  ClassDef(AliHLTTPCEVE, 0)
};
#endif
