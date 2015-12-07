// rm -f galice.root ;  aliroot -l -q -b $ALICE_SOURCE/HLT/global/physics/macros/testConfigOnlineCalib.C $ALICE_ROOT/HLT/exa/recraw-local.C'("raw.root","local:///opt/HLT-DEV/HCDB", 23, 23, "HLT", "chains=RootWriter ignore-hltout")'

void testConfigOnlineCalib()
{
	AliHLTPluginBase::InitInstance();

	AliHLTConfiguration calib1("myCalibration", "HLTAnalysisManagerComponent", "GLOBAL-flat-esd-converter", "-fPushEventModulo=3 -MinTracks=5 -QueueDepth=0 AddTaskMacro=$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/AddTaskTPCCalib.C(\"TPCCalib:CalibTimeDrift\")");
/*	AliHLTConfiguration calibmerge("myCalibrationMerger" , "RootObjectMerger" , "myCalibration" , "-QueueDepth 0 -cumulative");
	AliHLTConfiguration preproc("myTPCOfflinePreprocessor" , "TPCOfflinePreprocessorWrapper" , "myCalibrationMerger" , "-QueueDepth 0");

	AliHLTConfiguration eventTrigger("myCustomTrigger", "Zero", "TPC-DP", "");

	AliHLTConfiguration zmqsink("myZMQsink", "ZMQsink", "myTPCOfflinePreprocessor", "out=PUB@tcp://*:60203");
	AliHLTConfiguration zmqsource("myZMQsource", "ZMQsource", "myCustomTrigger", "in=SUB+tcp://localhost:60203");

	AliHLTConfiguration mapPrepare1("myMapPrepare1", "TPCClusterTransformationPrepare", "myZMQsource", "-QueueDepth 0 -MinSector 0 -MaxSector 35 -NoInitialObject");
	AliHLTConfiguration mapPrepare2("myMapPrepare2", "TPCClusterTransformationPrepare", "myZMQsource", "-QueueDepth 0 -MinSector 36 -MaxSector 71 -NoInitialObject");
	AliHLTConfiguration mapPreparemerge("myMapPrepare", "RootObjectMerger", "myMapPrepare1 myMapPrepare2", "-QueueDepth 0");*/
	
	TString clusterTransformation = "TPC-ClusterTransformation";
	AliHLTConfiguration overrideClusterTransformation(clusterTransformation.Data(), "TPCClusterTransformation", "TPC-HWCFDecoder", "-update-object-on-the-fly");
	/*AliHLTConfiguration overrideClusterTransformation(clusterTransformation.Data(), "TPCClusterTransformation", "TPC-HWCFDecoder myMapPrepare", "-update-object-on-the-fly");*/

	AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "myCalibration GLOBAL-esd-converter", "-directory testDir -datafile test.root");
	/*AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "myCalibrationMerger myTPCOfflinePreprocessor myZMQsink GLOBAL-esd-converter", "-directory testDir -datafile test.root");*/
}
