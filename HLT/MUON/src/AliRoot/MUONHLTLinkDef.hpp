////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedefs;

#pragma link C++ nestedclass;    // Makes the namespaces work properly.
#pragma link C++ nestedtypedef;  // Makes the namespaces work properly.
#pragma link C++ nestedfunction; // Makes the namespaces work properly.

#pragma link C++ namespace AliMUONHLT;

#pragma link C++ function AliMUONHLT::Version;
#pragma link C++ function AliMUONHLT::MajorVersion;
#pragma link C++ function AliMUONHLT::MinorVersion;
#pragma link C++ function AliMUONHLT::BuildNumber;

#pragma link C++ class AliMUONHLT::Region+;
#pragma link C++ class AliMUONHLT::Point+;
#pragma link C++ class AliMUONHLT::TriggerRecord+;
#pragma link C++ class AliMUONHLT::ADCStream+;
#pragma link C++ class AliMUONHLT::Track+;

#pragma link C++ class AliMUONHLT::MicrodHLT+;

#pragma link C++ class AliMUONHLT::ADCStreamSource+;
#pragma link C++ struct AliMUONHLT::ADCStreamSource::DataBlock+;

#pragma link C++ class AliMUONHLT::TriggerSource+;
#pragma link C++ enum AliMUONHLT::TriggerSource::AreaType+;
#pragma link C++ enum AliMUONHLT::TriggerSource::SourceType+;
#pragma link C++ class AliMUONHLT::TriggerSource::EventData+;

#pragma link C++ class AliMUONHLT::ClusterSource+;
#pragma link C++ enum AliMUONHLT::ClusterSource::AreaType+;
#pragma link C++ enum AliMUONHLT::ClusterSource::SourceType+;
#pragma link C++ class AliMUONHLT::ClusterSource::BlockData+;
#pragma link C++ class AliMUONHLT::ClusterSource::EventData+;

#pragma link C++ class AliMUONHLT::TrackSink+;
#pragma link C++ class AliMUONHLT::TrackSink::EventData+;

#pragma link C++ class AliMUONHLT::TrackerInterface+;
#pragma link C++ class AliMUONHLT::TrackerCallback+;

#pragma link C++ class AliMUONHLT::ClusterFinderInterface+;
#pragma link C++ class AliMUONHLT::ClusterFinderCallback+;

#endif // __CINT__

