void ZMQHLTchain(int configId=2, const char* outSocket="PUB@tcp://*:60324")
{
  enum
  {
    kGlobalHistoChainConfig = 1,
    kCompressionHistoChainConfig = 2,
    kCorrelationMinBiasHistoChainConfig = 3,
    kCorrelationSemiCentralHistoChainConfig = 4,
    kCorrelationCentralHistoChainConfig = 5,
    kCorrelationAllTriggerHistoChainConfig = 6,
    kEmcalHistoChainConfig = 7,
    kMuonHistoChainConfig = 8,
    kCorrelationsHistoChainConfig = 9,
    kClusterAttachmentConfig = 10
  };

	AliHLTSystem* system = AliHLTPluginBase::GetInstance();
	TString listOfChains = "";
	
	// Add additional libraries that need to be loaded here:
	system->LoadComponentLibraries("libAliHLTUtil.so");
	system->LoadComponentLibraries("libAliHLTGlobal.so");
	system->LoadComponentLibraries("libAliHLTTPC.so");
	system->LoadComponentLibraries("libAliHLTEMCAL.so");
	system->LoadComponentLibraries("libAliHLTMUON.so");
	
	AliHLTLogging log;
	
	switch (configId)
  {
    case kCompressionHistoChainConfig:
      {
        AliHLTConfiguration blockfilter("blockfilter", "BlockFilter", "source", "-datatype 'HWCLUST1' 'TPC ' -datatype 'REMCLSCM' 'TPC ' -datatype 'COMPDESC' 'TPC '");
        AliHLTConfiguration compressionHisto("compressionHisto", "TPCDataCompressorMonitor", "blockfilter", "-pushback-period=20");
        listOfChains += " compressionHisto";
      }
      break;
    case kClusterAttachmentConfig:
      {
        AliHLTConfiguration clusterAttachment("promptQA" , "PromptRecoQA" , "source", "-axis=tpcTrackPt,0,0,0 -axis=tpcClusterCharge,0,0,0 -reset -pushback-period=20 -ResetAfterPush=1");
        listOfChains += " promptQA";
      }
      break;
    case kGlobalHistoChainConfig:
      {
        const char* stdCuts="Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9";
        TString config;
        config += "-interval 20 -maxentries 250000 -maxmemory 4000000";
        config += " -histogram TrackPt(100,0,10) -size 1000 -expression Track_pt -cut Track_TPCclus>0";
        config += " -histogram TrackPhi(180,0,360) -size 1000 -expression Track_phi -cut Track_TPCclus>0";
        config += " -histogram TrackMultiplicity(60,1,600) -size 1000 -expression trackcount";
        config += " -histogram TrackMultiplicityTrend -size 1000 -expression trackcount:event";
        config += Form(" -histogram TrackMultiplicityPrimary(60,1,600) -size 1000 -expression Sum$(%s)",stdCuts);
        config += Form(" -histogram TrackMultiplicityPrimaryTrend -size 1000 -expression Sum$(%s):event",stdCuts);
        config += Form(" -histogram TrackMultiplicityBackground(60,1,600) -size 1000 -expression Sum$(!(%s))",stdCuts);
        config += Form(" -histogram TrackMultiplicityBackgroundTrend -size 1000 -expression Sum$(!(%s)):event",stdCuts);
        config += Form(" -histogram TrackEta(100,-2,2) -size 1000 -expression Track_eta -cut %s",stdCuts);
        config += Form(" -histogram TrackTPCclus(200,0,200) -size 1000 -expression Track_TPCclus -cut %s",stdCuts);
        config += " -histogram TrackITSclus(7,-0.5,6.5) -size 1000 -expression Track_ITSclus";
        config += Form(" -histogram TrackTheta(90,0,180) -size 1000 -expression Track_theta -cut %s",stdCuts);
        config += " -histogram TrackDCAr(200,-50,50) -size 1000 -expression Track_DCAr -cut Track_TPCclus>0";
        config += " -histogram TrackCharge(7,-1.75,1.75) -size 1000 -expression Track_charge -cut Track_TPCclus>0";
        config += " -histogram VertexXY(100,-2.5,2.5,100,-2.5,2.5) -size 1000 -expression vertexY:vertexX -cut nContributors>3 -opt colz";
        config += " -histogram VertexX(100,-2.5,2.5) -size 1000 -expression vertexX -cut nContributors>3";
        config += " -histogram VertexY(100,-2.5,2.5) -size 1000 -expression vertexY -cut nContributors>3";
        config += " -histogram VertexZ(600,-30,30) -size 1000 -expression vertexZ -cut nContributors>3";
        config += " -histogram VertexTrendX -size 1000 -expression vertexX:event -cut nContributors>3";
        config += " -histogram VertexTrendY -size 1000 -expression vertexY:event -cut nContributors>3";
        AliHLTConfiguration tpcGlobalHisto("tpcGlobalHisto", "GlobalHisto", "source", config.Data());
        listOfChains += " tpcGlobalHisto";
      }
      break;
    case kCorrelationMinBiasHistoChainConfig:
      {
        AliHLTConfiguration correlationHisto1("correlationHisto1", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2");
        AliHLTConfiguration renameHisto1("renameHisto1", "ObjectRenamer", "correlationHisto1", "-suffix _minbias");
        listOfChains += " renameHisto1";
      }
      break;
    case kCorrelationSemiCentralHistoChainConfig:
      {
        AliHLTConfiguration correlationHisto2("correlationHisto2", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVLN");
        AliHLTConfiguration renameHisto2("renameHisto2", "ObjectRenamer", "correlationHisto2", "-suffix _semicentral");
        listOfChains += " renameHisto2";
      }
      break;
    case kCorrelationCentralHistoChainConfig:
      {
        AliHLTConfiguration correlationHisto3("correlationHisto3", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVHN");
        AliHLTConfiguration renameHisto3("renameHisto3", "ObjectRenamer", "correlationHisto3", "-suffix _central");
        listOfChains += " renameHisto3";
      }
      break;
    case kCorrelationAllTriggerHistoChainConfig:
      {
        AliHLTConfiguration correlationHisto4("correlationHisto4", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2 -addTrigger CVLN -addTrigger CVHN");
        listOfChains += " correlationHisto4";
      }
      break;
    case kEmcalHistoChainConfig:
      {
        AliHLTConfiguration emcalHisto("emcalHisto", "EmcalClusterMonitor", "source", "-pushback-period=20");
        listOfChains += " emcalHisto";
      }
      break;
    case kMuonHistoChainConfig:
      {
        AliHLTConfiguration muonHisto("muonHisto", "MUONClusterHistogrammer", "source", "-pushback-period=20");
        listOfChains += " muonHisto";
      }
      break;
    case kCorrelationsHistoChainConfig:
      {
        AliHLTConfiguration correlationHisto1("correlationHisto1", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2");
        AliHLTConfiguration renameHisto1("renameHisto1", "ObjectRenamer", "correlationHisto1", "-suffix _minbias");
        listOfChains += " renameHisto1";
        AliHLTConfiguration correlationHisto2("correlationHisto2", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVLN");
        AliHLTConfiguration renameHisto2("renameHisto2", "ObjectRenamer", "correlationHisto2", "-suffix _semicentral");
        listOfChains += " renameHisto2";
        AliHLTConfiguration correlationHisto3("correlationHisto3", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CVHN");
        AliHLTConfiguration renameHisto3("renameHisto3", "ObjectRenamer", "correlationHisto3", "-suffix _central");
        listOfChains += " renameHisto3";
        AliHLTConfiguration correlationHisto4("correlationHisto4", "MultiplicityCorrelations", "source", "-pushback-period=20 -addTrigger CPBI1 -addTrigger CPBI2 -addTrigger CVLN -addTrigger CVHN");
        listOfChains += " correlationHisto4";
      }
      break;

    default:
      log.LoggingVarargs(kHLTLogError, "", FUNCTIONNAME(), __FILE__, __LINE__,
          "Could not configure the monitoring chain, unknown configId = %d.",
          configId
          );
      return;
  }
	
  //send blocks via ZMQ
  AliHLTConfiguration zmqSink("sink","ZMQsink", listOfChains.Data(), Form("out=%s",outSocket));
  printf("done configuring chain\n");
}

