// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCCALIBSEEDMAKERCOMPONENT_H
#define ALIHLTTPCCALIBSEEDMAKERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCalibSeedMakerComponent.h
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief A component to create TPC seeds from HLT clusters/tracks
*/

#include "AliHLTProcessor.h"

//forward declarations
class AliHLTTPCSpacePointData;
class AliTPCParam;
//class TClonesArray;
class TObjArray;
class AliTPCclusterMI;
class TH2F;

/**
 * @class AliHLTTPCCalibSeedMakerComponent
 * 
 *
 * @ingroup alihlt_tpc_components
 */ 

class AliHLTTPCCalibSeedMakerComponent : public AliHLTProcessor {
    
   public:
   
   /** standard constructor */    
   AliHLTTPCCalibSeedMakerComponent();           
   /** destructor */
   virtual ~AliHLTTPCCalibSeedMakerComponent();

      // Public functions to implement AliHLTComponent's interface.
      // These functions are required for the registration process
      
      /** interface function, see AliHLTComponent for description */
      const char* GetComponentID();							     
      /** interface function, see AliHLTComponent for description */
      void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			     
      /** interface function, see AliHLTComponent for description */
      AliHLTComponentDataType GetOutputDataType();					     
      /** interface function, see AliHLTComponent for description */
      int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   
      /** interface function, see AliHLTComponent for description */
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); 
      /** interface function, see AliHLTComponent for description */
      AliHLTComponent* Spawn(); 							   
      /** function for acting on the saving and cleaning histograms, after they are filled */
        
   protected:
	
      // Protected functions to implement AliHLTComponent's interface.
      // These functions provide initialization as well as the actual processing capabilities of the component. 

      int DoInit( int argc, const char** argv );
      int DoDeinit();
      int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );
      
      using AliHLTProcessor::DoEvent;
      
   private:
             
      /** copy constructor prohibited */
      AliHLTTPCCalibSeedMakerComponent(const AliHLTTPCCalibSeedMakerComponent&);

      /** assignment operator prohibited */
      AliHLTTPCCalibSeedMakerComponent& operator=(const AliHLTTPCCalibSeedMakerComponent&);
      AliTPCParam    *fTPCGeomParam;  //!transient
      
      //AliHLTTPCSpacePointData *fClustersArray[36][6]; //! transient
      //UInt_t                   fNSpacePoints[36][6];  //! transient
      //TClonesArray            *fSeedArray; // array of AliTPCseed objects
      //TObjArray               *seedArray; // array of AliTPCseed arrays
      TObjArray *fSeedArray;
      
      static const Int_t fkNPartition = 36*6;           // number of patches in TPC
      AliTPCclusterMI   *fPartitionClusters[fkNPartition];  //! arrays of cluster data for each TPC partition
      Int_t              fNPartitionClusters[fkNPartition]; //! number of clusters for each TPC partition
      
      TH2F *fdEdx; // energy loss vs. momentum
                    
      ClassDef(AliHLTTPCCalibSeedMakerComponent, 6)
    };

#endif
