void testConfigOnlineCalib()
{
	AliHLTSystem* pHLT = AliHLTPluginBase::GetInstance();

	if (1)
	{
//		AliHLTConfiguration calib1("myCalibration1" , "AsyncCalibration" , "GLOBAL-esd-converter" , "-QueueDepth 0");
//		AliHLTConfiguration calib2("myCalibration2" , "AsyncCalibration" , "GLOBAL-esd-converter" , "-QueueDepth 0");
		AliHLTConfiguration calib1("myCalibration1", "TPCCalibManagerComponent", "GLOBAL-flat-esd-converter", "-fPushEventModulo=3");
//		AliHLTConfiguration calib2("myCalibration2", "TPCCalibManagerComponent", "GLOBAL-flat-esd-converter", "-fPushEventModulo=3");
		AliHLTConfiguration calibmerge("myCalibrationMerger" , "TPCClusterTransformationMerger" , "myCalibration1" , "-cumulative");

		AliHLTConfiguration eventTrigger("myCustomTrigger", "Zero", "TPC-DP", "");

		AliHLTConfiguration zmqsink("myZMQsink", "ZMQsink", "myCalibrationMerger", "-ZMQsocketMode PUB -ZMQendpoint @tcp://*:60201");
		AliHLTConfiguration zmqsource("myZMQsource", "ZMQsource", "myCustomTrigger", "-ZMQsocketMode SUB -ZMQendpoint >tcp://localhost:60201");

		AliHLTConfiguration mapPrepare1("myMapPrepare1", "TPCClusterTransformationPrepare", "myZMQsource", "-QueueDepth 0 -MinSector 0 -MaxSector 35");
		AliHLTConfiguration mapPrepare2("myMapPrepare2", "TPCClusterTransformationPrepare", "myZMQsource", "-QueueDepth 0 -MinSector 36 -MaxSector 71");
		AliHLTConfiguration mapPreparemerge("myMapPrepare", "TPCClusterTransformationMerger", "myMapPrepare1 myMapPrepare2", "");

		TString clusterTransformation = "TPC-ClusterTransformation";
		AliHLTConfiguration overrideClusterTransformation(clusterTransformation.Data(), "TPCClusterTransformation", "TPC-HWCFDecoder myMapPrepare", "-initialize-on-the-fly");

		AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "myCalibrationMerger myZMQsink GLOBAL-esd-converter", "-directory testDir -datafile test.root");
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
