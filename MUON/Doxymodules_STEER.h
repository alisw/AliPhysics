// $Id$
// The class categories definitions for Doxygen

/** @defgroup STEER STEER
 *  Category of AliRoot steering classes
 *  @{
 */

/** @defgroup STEER0 STEER0
 *  Category of AliRoot steering classes
 *  @ingroup STEER
 *  @{
 */
// From STEERLinkDef.h revision 1.107

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
    class AliMagFMapsV1 {};
    class AliMagFDM {};
    class AliMagFCheb {};
    class AliCheb3DCalc {};
    class AliCheb3D {};
    class AliLego {};
    class AliLegoGenerator {};
    class AliLegoGeneratorXYZ {};
    class AliLegoGeneratorPhiZ {};
    class AliLegoGeneratorEta {};
    class AliLegoGeneratorEtaR {};
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

    class AliGeomManager {};
    class AliAlignObj {};
    class AliAlignObjParams {};
    class AliAlignObjMatrix {};

    class AliTrackFitter {};
    class AliTrackFitterRieman {};
    class AliTrackFitterKalman {};
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

    class  AliSplineFit {};

    class  AliDCSValue {};
    class  AliDCSSensor {};
    class  AliDCSSensorArray {};

    class  AliSurveyObj {};
    class  AliSurveyPoint {};

    class  AliCodeTimer {};
    class  AliCodeTimer::AliPair {};

    class  AliFstream {};
    class  AliCTPRawData {};

/** @} */

/** @defgroup ESD ESD
 *  Category of AliRoot event sumary data classes
 *  @ingroup STEER
 *  @{
 */
// From ESDLinkDef.h revision 1.35

    enum   AliLog::EType_t {};
    enum   AliESD::ESDListIndex_t {};

    class  AliESD {};
    class  AliESDRun {};
    class  AliESDHeader {};
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

    class  AliKFParticleBase {};
    class  AliKFParticle {};
    class  AliKFVertex {};

    class  AliKalmanTrack {};
    class  AliExternalTrackParam {};
    class  AliVertexerTracks {};
    class  AliStrLine {};
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
    class AliLHCTagCuts {};
    class AliDetectorTagCuts {};
    class AliEventTagCuts {};

    class AliTrackPointArray {};
    class AliTrackPoint {};

    class  AliESDFMD {};
    class  AliFMDMap {};
    class  AliFMDFloatMap {};

    class  AliESDVZERO {};
    class  AliESDTZERO {};

    class  AliESDMultITS {};
    class  AliMultiplicity {};

    class  AliSelector {};

    class  AliRawDataErrorLog {};

/** @} */

/** @defgroup CDB CDB
 *  Category of AliRoot Conditions database classes
 *  @ingroup STEER
 *  @{
 */
// From CDBLinkDef.h revision 1.9

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

    class AliPreprocessor {};

    class AliShuttleInterface {};

    class AliGRPPreprocessor {};
    class AliGRPDCS {};

/** @} */

/** @defgroup AOD AOD
 *  Category of AliRoot AOD classes
 *  @ingroup STEER
 *  @{
 */
// From AODLinkDef.h revision 1.8

    enum   AliAODVertex::AODVtx_t {};
    enum   AliAODTrack::AODTrk_t {};
    enum   AliAODTrack::AODTrkPID_t {};
    enum   AliAODCluster::AODClu_t {};
    enum   AliAODCluster::AODCluPID_t {};

    class  AliVParticle {};
    class  AliVEvent {};
    class  AliVHeader {};
    class  AliVEventHandler {};

    class  AliAODEvent {};
    class  AliAODHeader {};
    class  AliAODTrack {};
    class  AliAODVertex {};
    class  AliAODCluster {};
    class  AliAODJet {};
    class  AliAODPhoton {};
    class  AliAODRedCov<Int_t> {};
    class  AliAODRedCov<Int_t> {};
    class  AliAODRedCov<Int_t> {};
    class  AliAODRecoDecay;
    class  AliAODHandler {};

/** @} */
