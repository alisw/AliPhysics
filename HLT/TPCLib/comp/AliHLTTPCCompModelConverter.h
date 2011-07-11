// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELCONVERTER_H
#define ALIHLTTPCCOMPMODELCONVERTER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCompModelConverter.h
    @author Timm Steinbeck
    @date   
    @brief  Declaration of a copy component. */

#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCCompModelAnalysis.h"
#include "AliHLTLogging.h"
#include "AliHLTGlobalBarrelTrack.h"
#include <vector>

class AliHLTSpacePointContainer;

/**
 * @class AliHLTTPCCompModelConverter
 * @brief
 *
 */
class AliHLTTPCCompModelConverter: public AliHLTLogging
    {
    public:
      /** standard constructor */
      AliHLTTPCCompModelConverter();

      /** constructor including model/track analysis */
      AliHLTTPCCompModelConverter(AliHLTTPCCompModelAnalysis* modelanalysis);

      /** standard destructor */
      virtual ~AliHLTTPCCompModelConverter();
      
      /** initialisation function 
       * @return zero upon success
       */
      int Init();
      
      /** function to set input tracks
       * @param tracklets pointer to AliHLTTracksData
       * @param size of the buffer in byte
       * @return zero upon success
       */
      int SetInputTracks( const AliHLTTracksData* pTracks, unsigned sizeInByte );

      /** function to set input tracks
       * @param clusters pointer to AliHLTTPCClusterData
       * @param slice UInt_t slice number
       * @param patch UInt_t patch number
       * @return zero upon success
       */
      int SetInputClusters( AliHLTTPCClusterData* clusters, UInt_t slice, UInt_t patch );
      
      /** function to convert input to Vestbo-model (-> compression can then follow) */
      void Convert();
      
      /** function to get output model data size
       * @return unsigned long value of output model data size
       */
      unsigned long GetOutputModelDataSize() const;
      
      /** function to output model data
       * @param [in] data       pointer to target buffer
       * @param [in,out] dataSize   [in] size of buffer, [out] size of data
       * @return zero upon success
       */
      int OutputModelData( AliHLTUInt8_t* data, unsigned long& dataSize ) const;
      
      /** function to select remaining clusters */
      void SelectRemainingClusters();
      
      /** function to get remaining clusters output data size
       * @return unsigned long value = size
       */
      unsigned long GetRemainingClustersOutputDataSize() const;

      /** function to get remaining clusters
       * @param [in] data       pointer to target buffer
       * @param [in,out] dataSize   [in] size of buffer, [out] size of data
       */
      int GetRemainingClusters( AliHLTUInt8_t* const data, unsigned long& dataSize ) const;
      
      /** function to define minimal hits for one track
       * @param minHits unsigned
       */
      void SetMinHits( unsigned minHits )
      {
	fMinHits = minHits;
      }

    protected:

      /** function to expand track data */
      void ExpandTrackData();

      /** member variable input track array */
      vector<AliHLTGlobalBarrelTrack> fInputTrackArray;
      /** member variable output track array */
      AliHLTTPCTrackArray fOutputTrackArray;
      
      /** pointer to instance of analysis model class */
      AliHLTTPCCompModelAnalysis* fModelAnalysisInstance; // pointer to instance of analysis model class

      /** array of cluster data */
      AliHLTTPCClusterData* fClusters[36][6];

      /** temporary array how of used cluster sizes */
      unsigned long fClusterUsedSizes[36][6]; //! temporary, how many entries do we have in the following clusters used arrays
      /** array to determine the used clusters in tracking */
      bool* fClusterUsed[36][6]; //! keep track of which clusters have been used, can be transient.
      /** number of minimal hit for one track */
      unsigned fMinHits;
      
    private:
      /** copy constructor prohibited */
      AliHLTTPCCompModelConverter(const AliHLTTPCCompModelConverter&);
      /** assignment operator prohibited */
      AliHLTTPCCompModelConverter& operator=(const AliHLTTPCCompModelConverter&);
     
      ClassDef(AliHLTTPCCompModelConverter, 0)

    };
#endif
