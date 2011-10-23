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
// From STEERBaseLinkDef.h 40405 2010-04-14 14:41:08Z cvetan

    enum  AliLog::EType_t {};
 
    class AliVParticle {};
    class AliVTrack {};
    class AliVVertex {};
    class AliVEvent {};
    class AliVHeader {};
    class AliVEventHandler {};
    class AliVEventPool {};
    class AliVCuts {};

    class AliMixedEvent {};

    class AliPID {};
    class AliLog {};

    class AliRunTag {};
    class AliLHCTag {};
    class AliDetectorTag {};
    class AliEventTag {};

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

    class AliMagF {};
    class AliMagWrapCheb {};
    class AliCheb3DCalc {};
    class AliCheb3D {};


    class AliCodeTimer {};
    class AliCodeTimer::AliPair {};

    class AliPDG {};

    class AliTimeStamp {};
    class AliTriggerScalers {};
    class AliTriggerScalersRecord {};
    
    class AliExternalTrackParam {};
    class AliQA {};

    class AliITSPidParams {};
    class AliITSPIDResponse {};
    class AliTPCPIDResponse {};
    class AliTOFPIDResponse {};
    class AliTRDPIDResponse {};

    class AliDAQ {};

/** @} */

/** @defgroup STEER0 STEER0
 *  Category of AliRoot steering classes
 *  @ingroup STEER
 *  @{
 */
// From STEERLinkDef.h 38305 2010-01-15 14:30:37Z shahoian

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
    class  AliTaskLoader {};
    class  AliRunLoader {};
    class  AliReconstructor {};
    class  AliTrackMap {};
    class  AliMemoryWatcher {};
    class  AliMC {};
    class  AliSimulation {};
    class  AliReconstruction {};
    class  AliVertexGenFile {};
    class  AliVertexer {};
    class  AliV0vertexer {};
    class  AliCascadeVertexer {};

    class AliExpression {};
    class AliVariableExpression {};
    class AliTriggerInput {};
    class AliTriggerDetector {};
    class AliTriggerConfiguration {};
    class AliTriggerBCMask {};
    class AliTriggerInteraction {};
    class AliTriggerDescriptor {};
    class AliTriggerClass {};
    class AliCentralTrigger {};
    class AliTriggerCluster {};
    class AliTriggerPFProtection {};

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

    class  AliRieman;

    class AliExpression {};
    class AliVariableExpression {};
    class AliTriggerInput {};
    class AliTriggerDetector {};
    class AliTriggerConfiguration {};
    class AliTriggerBCMask {};
    class AliTriggerInteraction {};
    class AliTriggerDescriptor {};
    class AliTriggerClass {};
    class AliCentralTrigger {};
    class AliCTPRawStream {};
    class AliMathBase {};
    class AliSignalProcesor {};
    class  AliHelix {};
    class  AliCluster {};
    class  AliCluster3D {};
    class  AliTracker {};
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

/** @} */

/** @defgroup ESD ESD
 *  Category of AliRoot event sumary data classes
 *  @ingroup STEER
 *  @{
 */
// From ESDLinkDef.h 40103 2010-03-31 09:03:39Z belikov

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
    class  AliESDTrdTrack {};
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
    class  AliNeutralTrackParam {};
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

    class  AliTriggerIR {};

    class  AliESDVZEROfriend {};

    class  AliTriggerScalersESD {};
    class  AliTriggerScalersRecordESD {};
    class  AliESDHandler {};
    class  AliTrackerBase {};
    
/** @} */

/** @defgroup CDB CDB
 *  Category of AliRoot Conditions database classes
 *  @ingroup STEER
 *  @{
 */
// From CDBLinkDef.h 31411 2009-03-11 13:13:53Z zampolli

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

    class  TTreeDataElement {};
    class  TTreeStream {};
    class  TTreeSRedirector {};

/** @} */

/** @defgroup AOD AOD
 *  Category of AliRoot AOD classes
 *  @ingroup STEER
 *  @{
 */
// From AODLinkDef.h 37138 2009-11-23 12:19:02Z morsch

    enum   AliAODVertex::AODVtx_t {};
    enum   AliAODTrack::AODTrk_t {};
    enum   AliAODTrack::AODTrkPID_t {};
    enum   AliAODCluster::AODClu_t {};
    enum   AliAODCluster::AODCluPID_t {};

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
    class AliAODInputHandler {};
    class AliAODTracklets {};
    class AliAODTagCreator {};
    class AliAODCaloCells {};
    class AliAODDiJet {};
    class AliAODMCParticle {};
    class AliAODMCHeader {};
    class AliAODPWG4Particle {};
    class AliAODPWG4ParticleCorrelation {};
    class AliAODDimuon {};

/** @} */

/** @} */
