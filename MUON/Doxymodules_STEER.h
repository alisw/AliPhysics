// $Id$
// The class categories definitions for Doxygen

/** @defgroup STEER STEER
 *  Category of AliRoot steering classes
 *  @{
 */
// From STEERLinkDef.h revision 1.95

    class AliPDG {};

    class AliGenerator {};
    class AliVertexGenerator {};
    class AliRun {};
    class AliModule {};
    class AliDetector {};
    class AliDigit {};
    class AliHit {};
    class AliHeader {};
    class AliDisplay {};
    class AliPoints {};
    class AliMagF {};
    class AliMagFC {};
    class AliMagFCM {};
    class AliMagFMaps {};
    class AliMagFDM {};
    class AliLego {};
    class AliLegoGenerator {};
    class AliLegoGeneratorXYZ {};
    class AliLegoGeneratorPhiZ {};
    class AliLegoGeneratorEta {};
    class AliDigitNew {};
    class AliGeometry {};
    class AliRecPoint {};
    class AliSegmentation {};
    class AliHitMap {};
    class AliRndm {};
    class AliMCQA {};
    class AliDebugVolume {};
    class AliStack {};
    class AliConfig {};
    class AliGenEventHeader {};
    class AliDigitizer {};
    class AliRunDigitizer {};
    class AliStream {};
    class AliMergeCombi {};
    class AliFieldMap {};
    class AliGausCorr {};
    class AliLoader {};
    class AliDataLoader {};
    class AliBaseLoader {};
    class AliObjectLoader {};
    class AliTreeLoader {};
    class AliTaskLoader {};
    class AliRunLoader {};
    class AliTrackReference {};
    class AliReconstructor {};
    class AliTrackMap {};
    class AliTrackMapper {};
    class AliCollisionGeometry {};
    class AliMemoryWatcher {};
    class AliMC {};
    class AliSimulation {};
    class AliReconstruction {};
    class AliVertexGenFile {};
    class AliVertexer {};
    class AliVertexerTracks {};
    class AliStrLine {};
    class AliV0vertexer {};
    class AliCascadeVertexer {};

    class AliExpression {};
    class AliVariableExpression {};
    class AliTriggerInput {};
    class AliTriggerDetector {};
    class AliTriggerCondition {};
    class AliTriggerDescriptor {};
    class AliCentralTrigger {};

    class AliDetectorEventHeader {};

    class AliAlignObj {};
    class AliAlignObjAngles {};
    class AliAlignObjMatrix {};

    class AliTrackFitter {};
    class AliTrackFitterRieman {};
    class AliTrackFitterStraight {};
    class AliTrackResiduals {};
    class AliTrackResidualsChi2 {};
    class AliTrackResidualsFast {};
    class AliTrackResidualsLinear {};
    class AliAlignmentTracks {};

    class TTreeDataElement {};
    class  TTreeStream {};
    class  TTreeSRedirector {};

    class  AliRieman;

    class AliExpression {};
    class AliVariableExpression {};
    class AliTriggerInput {};
    class AliTriggerDetector {};
    class AliTriggerCondition {};
    class AliTriggerDescriptor {};
    class AliCentralTrigger {};
    class AliCTPRawData {};
    class AliCTPRawStream {};
    class AliMathBase {};
    class AliSignalProcesor {};
    class  AliHelix {};
    class  AliCluster {};
    class  AliClusterTGeo {};
    class  AliTracker {};
    class  AliV0 {};
    class  AliKink {};

    class  AliSelectorRL {};

    class  AliSplineFit {};


/** @} */

/** @defgroup STEER_ESD STEER_ESD
 *  Category of AliRoot event sumary data classes
 *  @{
 */
// From ESDLinkDef.h revision 1.29

    enum  AliLog::EType_t {};

    class  AliESD {};
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

    class  AliKalmanTrack {};
    class  AliExternalTrackParam {};
    class  AliLog {};
    class  AliPID {};
    class  AliTrackPointArray {};
    class  AliTrackPoint {};

    class AliRunTag {};
    class AliLHCTag {};
    class AliDetectorTag {};
    class AliEventTag {};

    class AliTagCreator {};
    class AliRunTagCuts {};
    class AliEventTagCuts {};

    class AliTrackPointArray {};
    class AliTrackPoint {};

    class  AliESDFMD {};
    class  AliFMDMap {};
    class  AliFMDFloatMap {};

    class  AliESDVZERO {};

    class  AliESDMultITS {};
    class  AliMultiplicity {};

    class  AliSelector {};

    class  AliRawDataErrorLog {};

/** @} */

/** @defgroup STEER_CDB STEER_CDB
 *  Category of AliRoot Conditions database classes
 *  @{
 */
// From CDBLinkDef.h revision 1.5

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

    class AliPreprocessor {};

    class AliShuttleInterface {};

    class AliGRPPreprocessor {};
    class AliGRPDCS {};

/** @} */

/** @defgroup STEER_AOD STEER_AOD
 *  Category of AliRoot AOD classes
 *  @{
 */
// From AODLinkDef.h revision 1.5

    enum   AliAODVertex::AODVtx_t {};
    enum   AliAODTrack::AODTrk_t {};
    enum   AliAODTrack::AODTrkPID_t {};
    enum   AliAODCluster::AODClu_t {};
    enum   AliAODCluster::AODCluPID_t {};

    class  AliAODEvent {};
    class  AliVirtualParticle {};
    class  AliAODHeader {};
    class  AliAODTrack {};
    class  AliAODVertex {};
    class  AliAODCluster {};
    class  AliAODJet {};
    class  AliAODRedCov<Int_t> {};

/** @} */
