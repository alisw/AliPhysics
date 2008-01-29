// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELDECONVERTER_H
#define ALIHLTTPCCOMPMODELDECONVERTER_H
/* TPCCompModelDeconverterright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full TPCCompModelDeconverterright notice                               */

/** @file   AliHLTTPCCompModelDeconverter.h
    @author Timm Steinbeck
    @date   
    @brief  Declaration of a copy component. */

#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTLogging.h"

/**
 * @class AliHLTTPCCompModelDeconverter
 * @brief A dummy HLT processing component. 
 *
 * An implementiation of a copy component that just copies its input data
 * to debug a components input data
 * @ingroup alihlt_tutorial
 */
class AliHLTTPCCompModelDeconverter: public AliHLTLogging
    {
    public:

      /** standard constructor */
      AliHLTTPCCompModelDeconverter();

      /** standard destructor */
      virtual ~AliHLTTPCCompModelDeconverter();
      
      /** initialisation function */
      int Init();
      
      /** function to set track cluster model input data
       * @param data AliHLTUInt8_t* pointer to input data
       * @param size UInt_t input size
       * @return zero upon success
       */      
      int SetTrackClusterModelInputData( AliHLTUInt8_t* data, UInt_t size );

      /** function to set remaining cluster model input data
       * @param data AliHLTUInt8_t* pointer to input data
       * @param size UInt_t input size
       * @return zero upon success
       */    
      int SetRemainingClustersModelInputData( AliHLTUInt8_t* data, UInt_t size );

      /** function to deconvert tracks
       * @param data AliHLTUInt8_t* pointer to input data
       * @param size UInt_t& input size
       * @return zero upon success
       */    
      int DeconvertTracks( AliHLTUInt8_t* data, UInt_t& size );

      /** function to deconvert clusters
       * @param slice UInt_t slice number
       * @param patch UInt_t patch number
       * @param data AliHLTUInt8_t* pointer to input data
       * @param size UInt_t& input size
       * @return zero upon success
       */    
      int DeconvertClusters( UInt_t slice, UInt_t patch, AliHLTUInt8_t* data, UInt_t& size );

    protected:

      /** member variable of input track array */
      AliHLTTPCTrackArray fInputTrackArray;
	
      /** member variable pointer to track cluster model data */
      AliHLTUInt8_t* fTrackClusterModelData;
      /** member variable for track cluster model data size */
      UInt_t fTrackClusterModelDataSize;

      /** member variable pointer to remaining clusters model data */
      AliHLTUInt8_t* fRemainingClustersModelData;
      /** member variable for remaining clusters model data size */
      UInt_t fRemainingClustersModelDataSize;

    private:
      /** copy constructor prohibited */
      AliHLTTPCCompModelDeconverter(const AliHLTTPCCompModelDeconverter&);
      /** assignment operator prohibited */
      AliHLTTPCCompModelDeconverter& operator=(const AliHLTTPCCompModelDeconverter&);

	ClassDef(AliHLTTPCCompModelDeconverter, 1)

    };
#endif
