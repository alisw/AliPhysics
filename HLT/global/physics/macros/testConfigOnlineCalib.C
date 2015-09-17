// rm -f galice.root ;  aliroot -l -q -b $ALICE_SOURCE/HLT/global/physics/macros/testConfigOnlineCalib.C $ALICE_ROOT/HLT/exa/recraw-local.C'("raw.root","local:///opt/HLT-DEV/HCDB", 23, 23, "HLT", "chains=RootWriter ignore-hltout")'

void testConfigOnlineCalib()
{
	AliHLTSystem* pHLT = AliHLTPluginBase::GetInstance();

	if (1)
	{
		AliHLTConfiguration calib1("myCalibration1", "HLTAnalysisManagerComponent", "GLOBAL-flat-esd-converter", "-fPushEventModulo=3  -MinTracks=5 -QueueDepth=0 AddTaskMacro=$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/AddTaskTPCCalib.C('TPCCalib:CalibTimeDrift')");
//		AliHLTConfiguration calib2("myCalibration2", "TPCCalibManagerComponent", "GLOBAL-flat-esd-converter", "-fPushEventModulo=3");
		AliHLTConfiguration calibmerge("myCalibrationMerger" , "RootObjectMerger" , "myCalibration1" , "-QueueDepth 0 -cumulative");
		AliHLTConfiguration preproc("myTPCOfflinePreprocessor" , "TPCOfflinePreprocessorWrapper" , "myCalibrationMerger" , "-QueueDepth 0");

		AliHLTConfiguration eventTrigger("myCustomTrigger", "Zero", "TPC-DP", "");

		AliHLTConfiguration zmqsink("myZMQsink", "ZMQsink", "myTPCOfflinePreprocessor", "-ZMQsocketMode PUB -ZMQendpoint @tcp://*:60201");
		AliHLTConfiguration zmqsource("myZMQsource", "ZMQsource", "myCustomTrigger", "-ZMQsocketMode SUB -ZMQendpoint >tcp://localhost:60201");

		AliHLTConfiguration mapPrepare1("myMapPrepare1", "TPCClusterTransformationPrepare", "myZMQsource", "-QueueDepth 0 -MinSector 0 -MaxSector 35");
		AliHLTConfiguration mapPrepare2("myMapPrepare2", "TPCClusterTransformationPrepare", "myZMQsource", "-QueueDepth 0 -MinSector 36 -MaxSector 71");
		AliHLTConfiguration mapPreparemerge("myMapPrepare", "RootObjectMerger", "myMapPrepare1 myMapPrepare2", "-QueueDepth 0");

		TString clusterTransformation = "TPC-ClusterTransformation";
		AliHLTConfiguration overrideClusterTransformation(clusterTransformation.Data(), "TPCClusterTransformation", "TPC-HWCFDecoder myMapPrepare", "-initialize-on-the-fly");

		AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "myCalibrationMerger myTPCOfflinePreprocessor myZMQsink GLOBAL-esd-converter", "-directory testDir -datafile test.root");
	}
	else if (0)
	{
		AliHLTConfiguration calib("myCalib", "TPCCalibManagerComponent", "GLOBAL-esd-converter", "");
		AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "myCalib", "-directory testDir -datafile test-calib.root");
	}
	else
	{
		AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "GLOBAL-esd-converter", "-directory testDir -datafile test-esdevent.root");
	}
}
