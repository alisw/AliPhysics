////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_CLUSTERING_IO_INTERFACE_HPP
#define dHLT_CLUSTERING_IO_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "EventID.hpp"
#include "TriggerRecord.hpp"
#include "RegionOfInterest.hpp"
#include "Cluster.hpp"
#include "Track.hpp"
#include "ADCStream.hpp"

namespace dHLT
{
  namespace Clustering
  {
    
    
    class IOInterface
    {
    public:
      
      // Do not free memory of ADC streams until the corresponding
      // ReleaseADCStream method is called.
      
      /* When a new ADC stream is available to the framework it is passed on
	 to a cluster finder with this method.
      */
      virtual void AddADCStream(const EventID event, const ADCStream* stream) = 0;
      
      /* When no more ADC streams are expected then this methods is called.
       */
      virtual void EndOfADCStreams(const EventID event) = 0;
      
    };
    
    
    class IOCallback
    {
    public:
      
      // Do not release the returned memory block until ReturnClusters is
      // called with the same memory pointer.
      virtual ClusterPoint* AllocateClusterBlock(const UInt size) = 0;
      
      /* Once clusters are found they are returned with this method to the framework.
       */
      virtual void ReturnClusters(const EventID event, const ROI region, ClusterPoint* clusters, const UInt count) = 0;
      
      /* When no more clusters can be found by the cluster finders then this
	 method is called.
      */
      virtual void EndOfClusters(const EventID event) = 0;
      
      // Do not release the specified ADC stream memory block until this method is called.
      virtual void ReleaseADCStream(const ADCStream* stream) = 0;
      
    };
    
    
    
  } // Clustering
} // dHLT

#endif // dHLT_CLUSTERING_IO_INTERFACE_HPP
