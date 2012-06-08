// $Id$
// The class categories definitions for Doxygen

/** @defgroup STEER STEER
 *  Category of AliRoot steering classes
 *  @{
 */

/** @defgroup STEERbase STEERbase
 *  Category of AliRoot steering classes
 *  @ingroup STEER
 *  @{
 */
// From STEERBaseLinkDef.h 56494 2012-05-15 20:58:42Z morsch

    enum   AliLog::EType_t {};
 
    class AliVParticle {};
    class AliVTrack {};
    class AliVCluster {};
    class AliVCaloCells {};
    class AliVVertex {};
    class AliVEvent {};
    class AliVHeader {};
    class AliVEventHandler {};
    class AliVEventPool {};
    class AliVCuts {};
    class AliVVZERO {};
    class AliVZDC {};
    class AliCentrality {};
    class AliEventplane {};

    class AliMixedEvent {};

    class AliPID {};
    class AliLog {};

    class AliRunTag {};
    class AliLHCTag {};
    class AliDetectorTag {};
    class AliEventTag {};
    class AliFileTag {};

    class AliRunTagCuts {};
    class AliLHCTagCuts {};
    class AliDetectorTagCuts {};
    class AliEventTagCuts {};

    class AliTagCreator {};

    class AliHeader {};
    class AliGenEventHeader {};
    class AliDetectorEventHeader {};
    class AliGenCocktailEventHeader {};
    class AliGenPythiaEventHeader {};
    class AliGenHijingEventHeader {};
    class AliCollisionGeometry {};
    class AliGenDPMjetEventHeader {};
    class AliGenHerwigEventHeader {};
    class AliGenGeVSimEventHeader {};
    class AliGenEposEventHeader {};
    class AliStack {};
    class AliMCEventHandler {};
    class AliInputEventHandler {};

    class AliTrackReference {};
    class AliSysInfo {};

    class AliMCEvent {};
    class AliMCParticle {};
    class AliMCVertex {};

    class  AliMagF {};
    class  AliMagWrapCheb {};
    class  AliCheb3DCalc {};
    class  AliCheb3D {};

    class  AliNeutralTrackParam {};

    class AliCodeTimer {};
    class AliCodeTimer::AliPair {};

    class  AliPDG {};

    class AliTimeStamp {};
    class AliTriggerScalers {};
    class AliTriggerScalersRecord {};

    class  AliExternalTrackParam {};
    class AliQA {};

    class AliTRDPIDReference {};
    class AliTRDPIDParams {};
    class AliTRDPIDParams::AliTRDPIDThresholds {};
    class AliITSPidParams {};
    class AliPIDResponse {};
    class AliITSPIDResponse {};
    class AliTPCPIDResponse {};
    class AliTPCdEdxInfo {};
    class AliTOFPIDResponse {};
    class AliTRDPIDResponse {};
    class AliEMCALPIDResponse {};
    class AliPIDCombined {};
    class AliTOFHeader {};

    class AliDAQ {};
    class AliRefArray {};

    class AliOADBContainer {};

    class AliMathBase {};
    class  TTreeDataElement {};
    class  TTreeStream {};
    class  TTreeSRedirector {};

    class AliVMFT {};
    class AliCounterCollection {};
    
    class AliVCaloTrigger {};

    class AliTOFPIDParams {};

/** @} */

/** @defgroup STEER0 STEER0
 *  Category of AliRoot steering classes
 *  @ingroup STEER
 *  @{
 */
// From STEERLinkDef.h 54207 2012-01-27 19:17:40Z hristov

    enum VertexSmear_t {};
    enum VertexSource_t {};

    class  AliGenerator {};
    class  AliVertexGenerator {};
    class  AliRun {};
    class  AliModule {};
    class  AliDetector {};
    class  AliDigit {};
    class  AliHit {};
    class  AliLego {};
    class  AliLegoGenerator {};
    class  AliLegoGeneratorXYZ {};
    class  AliLegoGeneratorPhiZ {};
    class  AliLegoGeneratorEta {};
    class  AliLegoGeneratorEtaR {};
    class  AliDigitNew {};
    class  AliGeometry {};
    class  AliRecPoint {};
    class  AliHitMap {};
    class  AliRndm {};
    class  AliDebugVolume {};
    class  AliConfig {};
    class  AliDigitizer {};
    class  AliDigitizationInput {};
    class  AliStream {};
    class  AliMergeCombi {};
    class  AliGausCorr {};
    class  AliLoader {};
    class  AliDataLoader {};
    class  AliBaseLoader {};
    class  AliObjectLoader {};
    class  AliTreeLoader {};
    class  AliRunLoader {};
    class  AliReconstructor {};
    class  AliMC {};
    class  AliSimulation {};
    class  AliReconstruction {};
    class  AliRecoInputHandler {};
    class  AliVertexGenFile {};
    class  AliVertexer {};

    class AliTriggerDetector {};
    class AliCentralTrigger {};
    class AliTriggerUtils {};

    class AliGeomManager {};
    class AliAlignObj {};
    class AliAlignObjParams {};
    class AliAlignObjMatrix {};
    class AliMisAligner {};

    class AliTrackFitter {};
    class AliTrackFitterRieman {};
    class AliTrackFitterKalman {};
    class AliTrackFitterStraight {};
    class AliTrackResiduals {};
    class AliTrackResidualsChi2 {};
    class AliTrackResidualsFast {};
    class AliTrackResidualsLinear {};
    class AliAlignmentTracks {};

    class  AliRieman {};

    class AliTriggerDetector {};
    class AliCentralTrigger {};
    class AliCTPRawStream {};
    class AliSignalProcesor {};
    class  AliHelix {};
    class  AliCluster {};
    class  AliCluster3D {};
    class  AliTracker {};
    class  AliTrackleter {};
    class  AliV0 {};
    class  AliKink {};

    class  AliSelectorRL {};

    class AliSurveyObj {};
    class AliSurveyPoint {};
    class AliSurveyToAlignObjs {};

    class AliFstream {};
    class AliCTPRawData {};

    class AliQADataMaker {};
    class AliQADataMakerSim {};
    class AliQADataMakerRec {};
    class AliCorrQADataMakerRec {};
    class AliGlobalQADataMaker {};
    class AliQAManager {};
    class AliQAChecker {};
    class AliCorrQAChecker {};
    class AliGlobalQAChecker {};
    class AliQACheckerBase {};
    class AliQAThresholds {};
    class AliMillepede {};

    class AliPlaneEff {};

    class AliTriggerRunScalers {};
    class AliGRPPreprocessor {};
    class AliGRPRecoParam {};

    class AliRelAlignerKalman {};

    class AliESDTagCreator {};

    class AliGRPObject {};

    class AliQAv1 {};

    class AliRunInfo {};
    class AliEventInfo {};
    class AliDetectorRecoParam {};
    class AliRecoParam {};

    class AliMillePede2 {};
    class AliMillePedeRecord {};
    class AliMinResSolve {};
    class AliMatrixSparse {};
    class AliVectorSparse {};
    class AliMatrixSq {};
    class AliSymMatrix {};
    class AliSymBDMatrix {};
    class AliRectMatrix {};
    class AliParamSolver {};

    class AliGRPManager {};
    class AliDCSArray {}; 	 
    class AliLHCReader {};
    class AliCTPTimeParams {};
    class AliCTPInputTimeParams {};

    class AliLHCDipValT<Double_t> {}; 	 
    class AliLHCDipValT<Int_t> {}; 	 
    class AliLHCDipValT<Float_t> {}; 	 
    class AliLHCDipValT<Char_t> {}; 	 
    class AliLHCData {};
    class AliLHCClockPhase {};

    class AliLTUConfig {};

    typedef AliLHCDipValD {}; 	 
    typedef AliLHCDipValI {}; 	 
    typedef AliLHCDipValF {}; 	 
    typedef AliLHCDipValC {};

/** @} */

/** @defgroup ESD ESD
 *  Category of AliRoot event sumary data classes
 *  @ingroup STEER
 *  @{
 */
// From ESDLinkDef.h 54829 2012-02-25 20:47:28Z morsch

    enum   AliESDEvent::ESDListIndex {};

    class  AliESD {};
    class  AliESDEvent {};
    class  AliESDInputHandler {};
    class  AliESDInputHandlerRP {};
    class  AliESDRun {};
    class  AliESDHeader {};
    class  AliESDHLTDecision {};
    class  AliESDZDC {};
    class  AliESDCaloTrigger {};
    class  AliESDfriend {};                                                                                                           
    class  AliESDtrack {};
    class  AliESDfriendTrack {};
    class  AliESDMuonTrack {};
    class  AliESDPmdTrack {};
    class  AliESDTrdTrigger {};
    class  AliESDTrdTrack {};
    class  AliESDTrdTracklet {};
    class  AliESDHLTtrack {};
    class  AliESDv0 {};
    class  AliESDcascade {};
    class  AliVertex {};
    class  AliESDVertex {};
    class  AliESDpid {};
    class  AliESDkink {};
    class  AliESDV0Params {};
    class  AliESDCaloCluster {};
    class  AliESDMuonCluster {};
    class  AliESDMuonPad {};

    class  AliKFParticleBase {};
    class  AliKFParticle {};
    class  AliKFVertex {};

    class  AliKalmanTrack {};
    class  AliVertexerTracks {};
    class  AliStrLine {};
    class  AliTrackPointArray {};
    class  AliTrackPoint {};

    class AliTrackPointArray {};
    class AliTrackPoint {};

    class  AliESDFMD {};
    class  AliFMDMap {};
    class  AliFMDFloatMap {};

    class  AliESDVZERO {};
    class  AliESDTZERO {};
    class  AliESDACORDE {};

    class  AliESDMultITS {};
    class  AliMultiplicity {};

    class  AliSelector {};

    class  AliRawDataErrorLog {};

    class  AliMeanVertex {};
    class  AliESDCaloCells {};

    class  AliESDVZEROfriend {};
    class  AliESDTZEROfriend {};

    class  AliESDHandler {};
    class  AliTrackerBase {};

    namespace AliESDUtils {};

    class  AliTriggerIR {};
    class  AliTriggerScalersESD {};
    class  AliTriggerScalersRecordESD {};
    class AliTriggerCluster {};
    class AliTriggerDescriptor {};
    class AliTriggerInput {};
    class AliTriggerInteraction {};
    class AliTriggerPFProtection {};
    class AliTriggerBCMask {};
    class AliTriggerClass {};
    class AliTriggerConfiguration {};
    class AliExpression {};
    class AliVariableExpression {};
    class AliESDCosmicTrack {};

    class  AliV0vertexer {};
    class  AliCascadeVertexer+ {};
    
/** @} */

/** @defgroup CDB CDB
 *  Category of AliRoot Conditions database classes
 *  @ingroup STEER
 *  @{
 */
// From CDBLinkDef.h 50616 2011-07-17 09:35:46Z hristov

    class AliCDBPath {};
    class AliCDBRunRange {};
    class AliCDBId {};
    class AliCDBMetaData {};
    class AliCDBEntry {};
    class AliCDBStorage {};
    class AliCDBStorageFactory {};
    class AliCDBManager {};
    class AliCDBParam {};
    class AliCDBLocal {};
    class AliCDBLocalFactory {};
    class AliCDBLocalParam {};
    class AliCDBDump {};
    class AliCDBDumpFactory {};
    class AliCDBDumpParam {}; 
    class AliCDBGrid {};
    class AliCDBGridFactory {};
    class AliCDBGridParam {};

    class AliDCSValue {};
    class AliDCSSensor {};
    class AliDCSSensorArray {};
    class AliDCSGenDB {};
    class  AliSplineFit {};

    class AliPreprocessor {};

    class AliShuttleInterface {};

    class AliGRPDCS {};
    class AliCDBHandler {};

    class  AliBaseCalibViewer {};
    class  AliBaseCalibViewerGUI {};
    class  AliCalibViewerGUItime {};

/** @} */

/** @defgroup AOD AOD
 *  Category of AliRoot AOD classes
 *  @ingroup STEER
 *  @{
 */
// From AODLinkDef.h 56945 2012-06-07 14:19:25Z fca

    enum   AliAODVertex::AODVtx_t {};
    enum   AliAODTrack::AODTrk_t {};
    enum   AliAODTrack::AODTrkPID_t {};

    class AliAODEvent {};
    class AliAODHeader {};
    class AliAODTrack {};
    class AliAODPid {};
    class AliAODVertex {};
    class AliAODCluster {};
    class AliAODCaloCluster {};
    class AliAODPmdCluster {};
    class AliAODFmdCluster {};
    class AliAODJet {};
    class AliAODJetEventBackground {};
    class AliAODPhoton {};
    class AliAODRedCov<3> {};
    class AliAODRedCov<4> {};
    class AliAODRedCov<6> {};
    class AliAODRecoDecay {};
    class AliAODv0 {};
    class AliAODcascade {};
    class AliAODHandler {};
    class AliAODExtension {};
    class AliAODBranchReplicator {};
    class AliAODInputHandler {};
    class AliAODTracklets {};
    class AliAODTagCreator {};
    class AliAODCaloCells {};
    class AliAODCaloTrigger {};
    class AliAODDiJet {};
    class AliAODMCParticle {};
    class AliAODMCHeader {};
    class AliAODPWG4Particle {};
    class AliAODPWG4ParticleCorrelation {};
    class AliAODDimuon {};
    class AliAODpidUtil {};
    class AliAODTZERO {};
    class AliAODVZERO {};
    class AliAODZDC {};

/** @} */

/** @} */
