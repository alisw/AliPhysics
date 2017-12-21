// #include <Riostream.h>
// #include <TChain.h>
// #include <TVectorT.h> 
// #include <TVector3.h> 
// #include <TFile.h>
// #include <TH1.h> 
// #include <TClonesArray.h>
// #include <TObjArray.h>
// #include <TObjString.h>
// #include <TFormula.h>
// #include <TString.h>
// #include <TRefArray.h>
// #include <TNtuple.h>
// #include <TArrayF.h>


#include "AliJFilter.h" 
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h" 
#include "AliMCEvent.h" 
#include "AliStack.h" 
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliESDCaloCluster.h" 
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliAODVertex.h" 
#include "AliAODTracklets.h" 
#include "AliAODPid.h" 
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliESDUtils.h"
//#include "AliESDVZERO.h" 
#include "AliCentrality.h" 
#include "AliAODTracklets.h"
#include "AliMultiplicity.h"
#include "AliJConst.h"
#include "AliESDRun.h"
#include "AliDAQ.h"
#include "AliESDVZERO.h"
#include "AliExternalTrackParam.h"
#include "AliHeader.h" 
//== EMCAL
#include "AliESDCaloCluster.h"
#include "AliEMCALGeometry.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALPIDUtils.h"

#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJPhoton.h"
//#include "AliJCaloCell.h"
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"

#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliPHOSGeoUtils.h"
#include "AliAnalysisUtils.h"


ClassImp(AliJFilter);

//______________________________________________________________________________
AliJFilter::AliJFilter() :   
	TNamed(),
	fEsdTrackCuts(0x0), 
	fESDFilter(0x0), 
	fIsRealOrMC(0),
	fStoreEventPlaneSource(0),
	fOADBPath(),
	fCaloClustersArr(0),
	fClusterThreshold(0),
	fTrackThreshold(0),
	fEventSuccess(0),
	fMcMap(0),
	fTrackList(0),
	fMCTrackList(0x0),
	fPhotonList(0x0),
	fCaloCellList(0x0),
	fHeaderList(0x0),
	fRunInfoList(0x0),
	fPIDResponse(0x0),
	fPIDCombined(0x0),
	fVZEROData(0x0), 
	fTZEROData(0x0), 
	//fFMDData(0x0), 
	fZDCData(0x0), 
	fEMCLabels(0),
	fEMCTreeLabels(0),
	fAliJRunHeader(0x0),
	fEMCALGeometry(0x0),    
	fEMCALRecoUtils(0x0),    
	fPHOSGeom(0x0),
	fAnaUtils(0x0),
	fMyTask(0x0)
{
	//Default constructor
}

//______________________________________________________________________________
AliJFilter::AliJFilter(const char *name,AliAnalysisTaskSE *task):
	TNamed(name,name), 
	fEsdTrackCuts(0x0), 
	fESDFilter(0x0), 
	fIsRealOrMC(0),
	fStoreEventPlaneSource(0),
	fOADBPath(),
	fCaloClustersArr(0),
	fClusterThreshold(0),
	fTrackThreshold(0),
	fEventSuccess(0),
	fMcMap(0),
	fTrackList(0),
	fMCTrackList(0x0),
	fPhotonList(0x0),
	fCaloCellList(0x0),
	fHeaderList(0x0),
	fRunInfoList(0x0),
	fPIDResponse(0x0),
	fPIDCombined(0x0),
	fVZEROData(0x0), 
	fTZEROData(0x0), 
	//fFMDData(0x0), 
	fZDCData(0x0), 
	fEMCLabels(0),
	fEMCTreeLabels(0),
	fAliJRunHeader(0x0),
	fEMCALGeometry(0x0),    
	fEMCALRecoUtils(0x0),    
	fPHOSGeom(0x0),
	fAnaUtils(0x0),
	fMyTask(0x0)
{
	// Constructor
	if(task->DebugLevel() > 5) cout << "---- AliJFilter Constructor ----"<<endl;

}

//____________________________________________________________________________
AliJFilter::AliJFilter(const AliJFilter& ap) :
	TNamed(ap.GetName(), ap.GetTitle()),
	fEsdTrackCuts(ap.fEsdTrackCuts), 
	fESDFilter(ap.fESDFilter), 
	fIsRealOrMC(ap.fIsRealOrMC),
	fStoreEventPlaneSource(ap.fStoreEventPlaneSource),
	fOADBPath(ap.fOADBPath),
	fCaloClustersArr(ap.fCaloClustersArr),
	fClusterThreshold(ap.fClusterThreshold),
	fTrackThreshold(ap.fTrackThreshold),
	fEventSuccess(ap.fEventSuccess),
	fMcMap(ap.fMcMap),
	fTrackList(ap.fTrackList),
	fMCTrackList(ap.fMCTrackList),
	fPhotonList(ap.fPhotonList),
	fCaloCellList(ap.fCaloCellList),
	fHeaderList(ap.fHeaderList),
	fRunInfoList(ap.fRunInfoList),
	fPIDResponse(ap.fPIDResponse),
	fPIDCombined(ap.fPIDCombined),
	fVZEROData(ap.fVZEROData), 
	fTZEROData(ap.fTZEROData), 
	//fFMDData(ap.fFMDData), 
	fZDCData(ap.fZDCData), 
	fEMCLabels(ap.fEMCLabels),
	fEMCTreeLabels(ap.fEMCTreeLabels),
	fAliJRunHeader(ap.fAliJRunHeader),
	fEMCALGeometry(ap.fEMCALGeometry),    
	fEMCALRecoUtils(ap.fEMCALRecoUtils),    
	fPHOSGeom(ap.fPHOSGeom),
	fAnaUtils(ap.fAnaUtils),
	fMyTask(ap.fMyTask)
{ 
	// cpy ctor
}

//_____________________________________________________________________________
AliJFilter& AliJFilter::operator = (const AliJFilter& ap)
{
	// assignment operator

	this->~AliJFilter();
	new(this) AliJFilter(ap);
	return *this;
}

//______________________________________________________________________________
AliJFilter::~AliJFilter()
{
	// destructor 
	delete fMcMap;
	delete fTrackList;
	delete fMCTrackList;
	delete fPhotonList;
	delete fCaloCellList;
	delete fHeaderList;
	delete fAliJRunHeader;
	delete fRunInfoList;
	delete fPIDResponse;
	delete fPIDCombined;
	delete fEMCALRecoUtils;
	delete fEMCALGeometry;
	delete fPHOSGeom;
	delete fAnaUtils;
	delete fVZEROData;
	delete fTZEROData;
	delete fZDCData;
	//  delete fFMDData;


}

//________________________________________________________________________

void AliJFilter::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fMyTask->DebugLevel() > 1) printf("AliJFilter::UserCreateOutPutData() \n");

	//== RUN HEADER
	cout<<"TEST2 "<<fAliJRunHeader<<endl;
	if(!fAliJRunHeader) fAliJRunHeader = new AliJRunHeader();
	fRunInfoList  = new TList();
	fRunInfoList->SetName("RunInfoList");
	fRunInfoList->SetOwner();
	fRunInfoList->Clear();
	fRunInfoList->Add(fAliJRunHeader);

	//=== Other Objects
	fCaloClustersArr = new TRefArray();
	fEMCALGeometry = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
	fEMCALRecoUtils = new AliEMCALRecoUtils();
	fPHOSGeom = new AliPHOSGeoUtils();
	fAnaUtils = new AliAnalysisUtils();
	fAnaUtils->SetUseOutOfBunchPileUp( kTRUE );
	fMcMap = new TArrayI();

	//=== Set Tree and TClonesArray
	//== TRACKS
	AddList("AliJTrackList", "AliJTrack", &fTrackList, 1000);
	if( fAliJRunHeader->GetStoreEMCalInfo() ){
		AddList("AliJPhotonList", "AliJPhoton", &fPhotonList, 1000);
		//BS AddList("AliJCaloCell", "AliJCaloCell", &fCaloCellList, 1000);
	}
	if( IsMC() ) 
		AddList("AliJMCTrackList", "AliJMCTrack", &fMCTrackList, 1000);
	//== Event Header
	AddList("AliJEventHeaderList", "AliJEventHeader", &fHeaderList, 1000);

	//== EventPlane SRC
	if( fAliJRunHeader->GetStoreEventPlaneSource() ){
		fVZEROData = new AliESDVZERO;
		fTZEROData = new AliESDTZERO;
		fZDCData   = new AliESDZDC;
	}
	//== PID
	//   fPIDCombined = new AliPIDCombined;
	//   fPIDCombined->SetDefaultTPCPriors();
	//   fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
	//   fPIDResponse = ((AliInputEventHandler*) (man->GetInputEventHandler()))->GetPIDResponse();
	//   fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
	//   if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());

	cout << "Add(fAliJRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJFilter::UserExec(Option_t* /*option*/) 
{
	// user loop
	AliJRunHeader *runh = fAliJRunHeader;
	Bool_t hasGoodTrack, hasGoodCluster;

	fEventSuccess = kFALSE;

	// Processing of one event
	DEBUG( 5, 1, "------- AliJFilter Exec-------" );
	if(!((fMyTask->Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  fMyTask->Entry())); 

	//=== Init Variables
	fTrackList->Clear();
	if( IsMC() ){
		fMCTrackList->Clear();
		fEMCLabels.clear();
		fEMCTreeLabels.clear();
	}

	if( fAliJRunHeader->GetStoreEMCalInfo() ){
		fPhotonList->Clear("C");
		fCaloCellList->Clear();
	}
	fHeaderList->Clear();

	hasGoodCluster = kTRUE;
	hasGoodTrack = kTRUE;

	//=== CHECK ESD, AOD, MC event
	if( !Event() ) return;

	if( FromESD() ) {   //Reading ESD  
		DEBUG( 5, 1, "\t------- Start ESD " );
		if( !ESDEvent() ) return;
		if( runh->GetWithoutSDD() && !(ESDEvent()->GetTriggerMask()  & (1<<13)) ) return;

		if( IsMC() ){
			if( ! MCEvent() ) return;
		}
	}

	if( FromAOD() ) {
		DEBUG( 5, 1, "\t------- Start AOD " );
		if( !AODEvent() ) return; 
	}


	// pileup rejection
	if( runh->IsPP() && fAnaUtils->IsPileUpEvent( Event() )) // Only APPLY for P+P
		return;

	//--------------------------------------------------------------- 
	// RUN Header
	//--------------------------------------------------------------- 
	if(!runh->GetRunNumber()){ //new run has started : I suppose no change of run in process
		runh->SetRunNumber( Event()->GetRunNumber() );
		if( FromESD() ){
			//==== General ====//
			runh->SetBeamEnergy( ESDEvent()->GetBeamEnergy() );
			runh->SetBeamType( ESDEvent()->GetBeamType() );
			//==== Detector status ==//
			if( ESDEvent()->GetCurrentL3() > 0 ) runh->SetL3MagnetFieldPolarity(1);
			if( ESDEvent()->GetCurrentL3() < 0 ) runh->SetL3MagnetFieldPolarity(-1);
			runh->SetL3MagnetFieldIntensity( ESDEvent()->GetMagneticField() );
			runh->SetCurrentL3( ESDEvent()->GetCurrentL3() );
			runh->SetCurrentDip( ESDEvent()->GetCurrentDip() );
			runh->SetUniformBMap( ESDEvent()->IsUniformBMap() );
			//==== Triggers ====//
			const AliESDRun* esdRun = ESDEvent()->GetESDRun();
			for(Int_t triggerBit=0; triggerBit<kRangeTriggerTableAlice; triggerBit++){
				runh->SetActiveTriggersAlice( triggerBit, esdRun->GetTriggerClass(triggerBit) );
			}
		}
		else if( FromAOD() ){
			//==== General ====//
			cout << "Run # = "<< AODEvent()->GetRunNumber() << endl;
			runh->SetRunNumber( AODEvent()->GetRunNumber() );
			//TODO runh->SetBeamEnergy( ESDEvent()->GetBeamEnergy() );
			//TODO runh->SetBeamType( ESDEvent()->GetBeamType() );
			//==== Detector status ==//
			//TODO runh->Setl3MgFieldPolarity(1);
			runh->SetL3MagnetFieldIntensity( AODEvent()->GetMagneticField() );
			runh->SetCurrentL3( AODEvent()->GetMagneticField()*30000.0/5.00668 );
			runh->SetCurrentDip( AODEvent()->GetMuonMagFieldScale()*6000.0 );
			runh->SetUniformBMap( kFALSE ); // TODO is this?
		}
		cout << "Add(fAliJRunHeader) is done =============" << endl;
	}

	//--------------------------------------------------------------- 
	// EventHeader and read Others
	//--------------------------------------------------------------- 
	if( FromESD() ){   //Reading ESD  
		DEBUG( 5, 1, "\t------- Start READ ESD " );

		ReadESDHeader( ESDEvent() );
		ReadESDTracks( ESDEvent() );

		if( fAliJRunHeader->GetStoreEMCalInfo() ){
			ReadESDCaloClusters( ESDEvent() );
			ReadESDCaloCells( ESDEvent() );
		}
		if( IsMC() ){
			ReadMCTracksFromESD();
			//RemapMCLabels();
		}
	}else if( FromAOD() ){ 
		DEBUG( 5, 1, "\t------- Start READ AOD " );
		ReadAODHeader( AODEvent() );
		ReadAODTracks( AODEvent() );
		if( fAliJRunHeader->GetStoreEMCalInfo() ){
			ReadAODCaloClusters( AODEvent() );
			ReadAODCaloCells( AODEvent() );
		}
		if( IsMC() ){
			ReadMCTracksFromAOD();
			//RemapMCLabels();
		}

	}else{
		cout << "Error: Not correct InputDataFormat especified " << endl;
		return;
	}

	if( hasGoodCluster || hasGoodTrack ){
		//=== TODO : need this?
		AliAODHandler* outputHandler = 
			(AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
		outputHandler->SetFillAOD(kTRUE);
		outputHandler->SetFillExtension(kTRUE);
		fEventSuccess = kTRUE;
	}
	else{
		fTrackList->Clear();
		if( IsMC() ){
			fMCTrackList->Clear();
			fEMCLabels.clear();
			fEMCTreeLabels.clear();
		}

		if( fAliJRunHeader->GetStoreEMCalInfo() ){
			fPhotonList->Clear("C");
			fCaloCellList->Clear();
		}
		fHeaderList->Clear();
	}

	DEBUG( 5, 1, "\t------- End UserExec " );
}

//______________________________________________________________________________
void AliJFilter::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

	//   TString formula(fEsdTrackCuts->GetMaxDCAToVertexXYPtDep());
	//   if(formula.Length()>0){ // momentum dep DCA cut for AOD
	//     formula.ReplaceAll("pt","x");
	//   }
}

//______________________________________________________________________________
void AliJFilter::Terminate(Option_t *)
{
	// termination
	fTrackList->Clear();
	if( IsMC() ) fMCTrackList->Clear();
	if( fAliJRunHeader->GetStoreEMCalInfo() ){
		fPhotonList->Clear();
		fCaloCellList->Clear();
	}
	fHeaderList->Clear();

	// Processing when the event loop is ended
	cout<<"PWG4JCORRAN Analysis DONE !!"<<endl; 

}

//______________________________________________________________________________
void AliJFilter::ReadESDTracks(AliESDEvent * esd)
	//void AliJFilter::ReadESDTracks(const AliESDEvent * esd)
{
	// Read the AliESDtrack and fill the list of AliJTrack containers
	Int_t nt = esd->GetNumberOfTracks();
	DEBUG( 5, 1 , Form("ESD::NumberOfTracks = %d",nt), "AliJFilter::ReadESDTracks" ); 

	//==== Prepare TPC, GCG track ====//
	Float_t ptMaxTPC = 0;
	Float_t ptMinTPC = 1E10;
	Float_t ptMaxGCG = 0;
	Float_t ptMinGCG = 1E10;
	for(int i = 0;i<32;i++){
		AliESDtrackCuts* cuts = (AliESDtrackCuts*)fESDFilter->GetCuts()->At(i);
		if(!cuts) continue;
		Float_t tmp1= 0,tmp2 = 0;
		cuts->GetPtRange(tmp1,tmp2);
		if( TESTBIT ( fAliJRunHeader->GetStoreTPCTrackBitMask(), i ) ){
			if(tmp1<ptMinTPC)ptMinTPC=tmp1;
			if(tmp2>ptMaxTPC)ptMaxTPC=tmp2;
		}
		if( TESTBIT(fAliJRunHeader->GetStoreGCGTrackBitMask() , i ) ){
			if(tmp1<ptMinGCG)ptMinGCG=tmp1;
			if(tmp2>ptMaxGCG)ptMaxGCG=tmp2;
		}
	} 

	//==== loop over tracks ====//
	for(Int_t it = 0; it < nt; it++) { 

		AliESDtrack *track = esd->GetTrack(it);
		if( !track ) continue;
		UInt_t filterMap = fESDFilter->IsSelected( track );
		if(! filterMap ) continue; // apply track selection criteria

		//====create a new AliJTrack and fill the track info
		AliJTrack * ctrack = new( (*fTrackList)[fTrackList->GetEntriesFast()] ) AliJTrack;
		ctrack->SetPxPyPzE(track->Px(), track->Py(), track->Pz(), 0 );
		Double32_t pos[3];
		track->GetXYZ(pos);
		ctrack->SetTrackPos( pos );
		ctrack->SetTPCdEdx( track->GetTPCsignal()  );
		ctrack->SetParticleType(kJNone);
		ctrack->SetCharge(track->Charge());
		ctrack->SetFilterMap( filterMap );
		ctrack->SetLabel( track->GetLabel() );

		ReadESDPID( track, ctrack );
		//==== TPC Tracks ====//
		if( filterMap & fAliJRunHeader->GetStoreTPCTrackBitMask() ) {
			ConvertESDTPCOnlyTracks( esd, it, ctrack, ptMinTPC, ptMaxTPC );
		}
		//==== GCG Tracks ====//
		if( filterMap & fAliJRunHeader->GetStoreGCGTrackBitMask() ) {
			ConvertESDGCGTracks( esd, it, ctrack, ptMinGCG, ptMaxGCG );
		}

		Float_t b[2];
		Float_t bCov[3];
		track->GetImpactParameters(b,bCov);
		//         ctrack->SetDCAtoVertexXY( b[0] );
		//         ctrack->SetDCAtoVertexZ( b[1] );

		//if( track->P()>1 ) DEBUG( 5, 1, Form("P = %f", track->P() ) ) ;

	} // end tracks loop
}

//______________________________________________________________________________
void AliJFilter::ConvertESDTPCOnlyTracks(AliESDEvent* esd, int iTrack, AliJTrack * ctrack, double ptmin, double ptmax)
{

	const AliESDVertex *vtxSPD = esd->GetPrimaryVertexSPD();

	Double_t pos[3] = { 0. };      
	Double_t covTr[21]={0.};
	//Double_t pid[10]={0.};  

	Double_t p[3] = { 0. };

	Double_t pDCA[3] = { 0. }; // momentum at DCA
	Double_t rDCA[3] = { 0. }; // position at DCA
	Float_t  dDCA[2] = {0.};    // DCA to the vertex d and z
	Float_t  cDCA[3] = {0.};    // covariance of impact parameters


	AliESDtrack* esdTrack = esd->GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 

	// Track selection

	AliESDtrack *track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esd),esdTrack->GetID());
	if(!track) return;

	if(track->Pt()>0.)
	{
		// only constrain tracks above threshold
		AliExternalTrackParam exParam;
		// take the B-field from the ESD, no 3D fieldMap available at this point
		Bool_t relate = false;
		relate = track->RelateToVertexTPC(vtxSPD,esd->GetMagneticField(),kVeryBig,&exParam);
		if(!relate){
			delete track;
			return;
		}
		// fetch the track parameters at the DCA (unconstraint)
		if(track->GetTPCInnerParam()){
			track->GetTPCInnerParam()->GetPxPyPz(pDCA);
			track->GetTPCInnerParam()->GetXYZ(rDCA);
		}
		// get the DCA to the vertex:
		track->GetImpactParametersTPC(dDCA,cDCA);
		// set the constrained parameters to the track
		track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
	}

	track->GetPxPyPz(p);

	double p2[3];
	esdTrack->GetInnerPxPyPz(p2);

	Float_t pT = track->Pt();
	if(pT<ptmin||pT>ptmax){
		delete track;
		return;
	}

	track->GetXYZ(pos);
	track->GetCovarianceXYZPxPyPz(covTr);

	ctrack->SetTPCTrack(p[0], p[1], p[2]);

	delete track;
}


void AliJFilter::ConvertESDGCGTracks(AliESDEvent *esd, int iTrack, AliJTrack *ctrack, double ptMin, double ptMax)
{

	Double_t pos[3] = { 0. };      
	Double_t covTr[21]={0.};
	Double_t p[3] = { 0. };

	Double_t pDCA[3] = { 0. }; // momentum at DCA
	Double_t rDCA[3] = { 0. }; // position at DCA
	Float_t  dDCA[2] = {0.};    // DCA to the vertex d and z
	Float_t  cDCA[3] = {0.};    // covariance of impact parameters


	AliESDtrack* esdTrack = esd->GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
	const AliExternalTrackParam * exParamGC = esdTrack->GetConstrainedParam();
	if(!exParamGC) return;

	// fetch the track parameters at the DCA (unconstrained)
	esdTrack->GetPxPyPz(pDCA);
	esdTrack->GetXYZ(rDCA);
	// get the DCA to the vertex:
	esdTrack->GetImpactParameters(dDCA,cDCA);

	if (!esdTrack->GetConstrainedPxPyPz(p)) return;


	Float_t pT = exParamGC->Pt();
	if(pT<ptMin||pT>ptMax){
		return;
	}

	esdTrack->GetConstrainedXYZ(pos);
	exParamGC->GetCovarianceXYZPxPyPz(covTr);

	ctrack->SetGCGTrack(p[0], p[1], p[2]);
}



//_________________________________________________________________________________-
void AliJFilter::ReadESDPID(AliESDtrack *track, AliJTrack *ctrack)
{
	// To reduce the size of output, the variables which cannot be calculated later are only kept
	// expected TOF signal, TPC momentum for expected TPC signal. Measured values are stored in ReadESDTrack()
	// 1. expected TOF signal
	Double_t times[AliPID::kSPECIES];
	track->GetIntegratedTimes(times);
	for(int ip=0; ip < (AliJTrack::kNAliJTrkPID); ip++) {
		ctrack->SetExpectedTOFsignal(AliJTrack::AliJTrkPID(ip), times[ip]);

	}
	// 2. TPC momentum
	Double_t momTPC = track->GetTPCmomentum();
	ctrack->SetTPCmomentum(momTPC);
}

//______________________________________________________________________________
Bool_t AliJFilter::ReadAODTracks(const AliAODEvent * aod)
{
	// AOD track reader
	Bool_t hasGoodTrack;
	hasGoodTrack = kFALSE;

	// Read the AliAODtrack and fill the list of AliJTrack containers
	Int_t nt = aod->GetNumberOfTracks();
	Int_t listnt = 0;

	DEBUG(5, 1, Form("AOD::NumberOfTracks = %d",nt) );

	//==== loop over tracks ====//
	for(Int_t it = 0; it < nt; it++) { 

        AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
        if(!track) {
            AliFatal("Not a standard AOD track");
            continue;
        }
		//if(track->GetFilterMap() & (1 << 7) ) continue;
		//if(!AcceptAODTrack(track)) continue; 
		//if(! fEsdTrackCuts->IsSelected(track)) continue; //apply loose selection criteria
		//FK//if(track->GetType() != AliAODTrack::kPrimary) continue; // only primaries 
		//

		AliJTrack * ctrack = new( (*fTrackList)[listnt++] ) AliJTrack;
		ctrack->SetID( track->GetID() );
		ctrack->SetPxPyPzE(track->Px(), track->Py(), track->Pz(), 0 );
		Double32_t pos[3];
		track->GetXYZ(pos);
		ctrack->SetTrackPos( pos );
		//TODO if( fStoreTPCTrack )
		ctrack->SetParticleType(kJNone);
		ctrack->SetCharge(track->Charge());
		ctrack->SetStatus(track->GetStatus());//
		ctrack->SetFlags( track->GetFlags() );
		ctrack->SetLabel( track->GetLabel() );
		//     //FilterMap
		//     UInt_t filterMap=0;
		//     for( unsigned int i=0;i<sizeof(filterMap)*8;i++ ){
		//       if( track->TestFilterBit( BIT(i) )){
		//         SETBIT( filterMap ,  i);
		//       }
		//     }
		ctrack->SetFilterMap( track->GetFilterMap() );

		//PID TODO
		double const * pid = track->PID();
		//ctrack->SetPID(AliJTrack::kElectronAliJ,pid[AliAODTrack::kElectron],AliJTrack::kTOF);
		//ctrack->SetPID(AliJTrack::kMuonAliJ,    pid[AliAODTrack::kMuon],    AliJTrack::kTOF);
		//ctrack->SetPID(AliJTrack::kPionAliJ,    pid[AliAODTrack::kPion],    AliJTrack::kTOF);
		//ctrack->SetPID(AliJTrack::kKaonAliJ,    pid[AliAODTrack::kKaon],    AliJTrack::kTOF);
		//ctrack->SetPID(AliJTrack::kProtonAliJ,  pid[AliAODTrack::kProton],  AliJTrack::kTOF);
		//TPC
		ctrack->SetTPCnClust(track->GetTPCNcls());
		ctrack->SetTPCdEdx( track->GetTPCsignal()  );
		ctrack->SetTOFsignal( track->GetTOFsignal() );
		ctrack->SetLabel( track->GetLabel() );
		for( int i=0;i<int(sizeof(UInt_t)*8);i++ ){
			ctrack->SetBit( i, track->TestBit( i ));
		}

		// check track threshold
		if( track->Pt() > fTrackThreshold )
			hasGoodTrack = kTRUE;

		//if(fMyTask->DebugLevel() > 5 && track->P()>1 ) cout << "P = " << track->P() << endl;
	} // end tracks loop

	return hasGoodTrack;
}


//______________________________________________________________________________
AliJEventHeader* AliJFilter::ReadCommonHeader(AliVEvent *event){
	//Read the AliVEvent and fill the list of AliJEventHeader containers
	//create a header and fill it
	AliJEventHeader *hdr = new( (*fHeaderList)[fHeaderList->GetEntriesFast()] ) AliJEventHeader;


	// Get Centrality as a percent from 0% to 100%
	AliCentrality *cent = event->GetCentrality();
	if( cent ){
		hdr->SetCentrality( cent->GetCentralityPercentile("V0M"));
		hdr->SetCentralityArray(AliJEventHeader::kcV0M, cent->GetCentralityPercentile("V0M"));
		hdr->SetCentralityArray(AliJEventHeader::kcFMD, cent->GetCentralityPercentile("FMD"));
		hdr->SetCentralityArray(AliJEventHeader::kcTRK, cent->GetCentralityPercentile("TRK"));
		hdr->SetCentralityArray(AliJEventHeader::kcTKL, cent->GetCentralityPercentile("TKL"));
		hdr->SetCentralityArray(AliJEventHeader::kcCL0, cent->GetCentralityPercentile("CL0"));
		hdr->SetCentralityArray(AliJEventHeader::kcCL1, cent->GetCentralityPercentile("CL1"));
		hdr->SetCentralityArray(AliJEventHeader::kcV0MvsFMD, cent->GetCentralityPercentile("V0MvsFMD"));
		hdr->SetCentralityArray(AliJEventHeader::kcTKLvsV0, cent->GetCentralityPercentile("TKLvsV0"));
		hdr->SetCentralityArray(AliJEventHeader::kcZEMvsZDC, cent->GetCentralityPercentile("ZEMvsZDC"));
		hdr->SetCentralityArray(AliJEventHeader::kcV0A, cent->GetCentralityPercentile("V0A"));
		hdr->SetCentralityArray(AliJEventHeader::kcV0C, cent->GetCentralityPercentile("V0C"));
	}
	hdr->SetTriggerMaskAlice(event->GetTriggerMask()); //ULong64_t
	hdr->SetTriggerMaskJCorran(ConvertTriggerMask()); //UInt_t
	hdr->SetEventType(event->GetEventType());
	hdr->SetBunchCrossNumber(event->GetBunchCrossNumber());

	int ncontributors = 0;
	const AliVVertex * vtxESD = event->GetPrimaryVertex();
	if(vtxESD){
		hdr->SetXVertex(vtxESD->GetX()); //FK// EFF
		hdr->SetYVertex(vtxESD->GetY()); //FK// EFF
		hdr->SetZVertex(vtxESD->GetZ());
		//hdr->SetZVertexErr(vtxESD->GetZRes());
		double covMat[6];
		vtxESD->GetCovarianceMatrix(covMat);
		hdr->SetZVertexErr(TMath::Sqrt(covMat[5])); // GetZRes := TMath::Sqrt(fCovZZ)
		ncontributors = vtxESD->GetNContributors(); // get number of contributors to vertex 
		hdr->SetVtxMult( vtxESD->GetNContributors() );
	}else{
		hdr->SetZVertex(9999);
		hdr->SetZVertexErr(9999);
	}
	hdr->SetVtxMult(ncontributors); //FK// EFF contrib to vertex
	return hdr;
}
//______________________________________________________________________________
void AliJFilter::ReadESDHeader(AliESDEvent *esd)
{
	// Read the AliESDEvent and fill the list of AliJEventHeader containers
	if(!esd) return;
	if( fAliJRunHeader->GetRefitESDVertexTracks() )
		AliESDUtils::RefitESDVertexTracks( esd ); // TODO only for LHC11a right?
	AliJEventHeader *hdr = ReadCommonHeader( esd );
	//   AliMultiplicity *fSPDMult =(AliMultiplicity *) esd->GetMultiplicity();
	//   if(fSPDMult) hdr->SetSPDTrackletMult(fSPDMult->GetNumberOfTracklets());
	// This is moved from ReadCommonHeader. AOD should have same.TODO!!
	AliESDVZERO *v0 = esd->GetVZEROData();

	if( v0 ) hdr->SetV0Mult(v0->GetMTotV0A() + v0->GetMTotV0C());
	if( v0 ) hdr->SetV0AMult(v0->GetMTotV0A());
	if( v0 ) hdr->SetV0CMult(v0->GetMTotV0C());

	const AliESDRun* esdRun = esd->GetESDRun();
	//cout <<"========================"<<endl;
	//cout << (esdRun->GetDetectorsInReco() & AliDAQ::kSPD) << endl;
	//cout << (esdRun->GetDetectorsInReco() & AliDAQ::kSSD) << endl;
	//cout << (esdRun->GetDetectorsInReco() & AliDAQ::kSDD) << endl;
	//cout << (esdRun->GetDetectorsInReco() & AliDAQ::kTPC) << endl;
	if(esdRun->GetDetectorsInReco() & AliDAQ::kSPD) hdr->SetSPDTrackletMult(AliESDtrackCuts::GetReferenceMultiplicity( esd, AliESDtrackCuts::kTracklets, 1.0 ));
	if((esdRun->GetDetectorsInReco() & AliDAQ::kSSD) || (esdRun->GetDetectorsInReco() & AliDAQ::kSDD)) hdr->SetITSSATrackletMult(AliESDtrackCuts::GetReferenceMultiplicity( esd, AliESDtrackCuts::kTrackletsITSSA, 1.0 ));
	if(esdRun->GetDetectorsInReco() & AliDAQ::kTPC) hdr->SetITSTPCTrackletMult(AliESDtrackCuts::GetReferenceMultiplicity( esd, AliESDtrackCuts::kTrackletsITSTPC, 1.0 ));


	//TODO  Store Detector data
	if( fAliJRunHeader->GetStoreEventPlaneSource() ){
		*fVZEROData = *esd->GetVZEROData();
		*fTZEROData = AliESDTZERO(*esd->GetESDTZERO());
		*fZDCData  = *esd->GetESDZDC();
	}
	hdr->SetEventID( esd->GetEventNumberInFile());
	//const AliESDVertex * vtxESD = esd->GetPrimaryVertex();
	//if( vtxESD->GetStatus() == 0 ) hdr->SetVtxMult( 0 );
	// if fNcontributes > 0 then status is always true. do we need this?

	//==== MC ====/
	if( IsMC() ){
		const AliVVertex * primaryMCVertex = MCEvent()->GetPrimaryVertex();
		//cout<<"AliMCEvent = "<<MCEvent()<<endl;
		//cout<<"AliVVertex = "<<primaryMCVertex<<endl;
		if( primaryMCVertex ){
			hdr->SetXVertexMC( primaryMCVertex->GetX() );
			hdr->SetYVertexMC( primaryMCVertex->GetY() );
			hdr->SetZVertexMC( primaryMCVertex->GetZ() );
		}
		AliESDHeader * esdHeader = esd->GetHeader();
		hdr->SetL0TriggerInputs( esdHeader->GetL0TriggerInputs() );
	}
}

//______________________________________________________________________________
void AliJFilter::ReadAODHeader(AliAODEvent *aod)
{  
	//Read the AliAODEvent and fill the list of AliJEventHeader containers
	AliJEventHeader *hdr = ReadCommonHeader( aod );

	const AliAODTracklets *trackletsSPD = aod->GetTracklets();
	if(trackletsSPD){
		hdr->SetSPDTrackletMult(trackletsSPD->GetNumberOfTracklets());
	}

	hdr->SetFiredTriggers( aod->GetFiredTriggerClasses() );
	//TODO hdr->SetEventID( esd->GetEventNumberInFile());
	//==== MC ====//
	if( IsMC() ){
		AliAODMCHeader *aodMCheader = (AliAODMCHeader *) aod->FindListObject(AliAODMCHeader::StdBranchName());
		hdr->SetXVertexMC( aodMCheader->GetVtxX() );
		hdr->SetYVertexMC( aodMCheader->GetVtxY() );
		hdr->SetZVertexMC( aodMCheader->GetVtxZ() );
	}

	AliAODHeader * ah = dynamic_cast<AliAODHeader*>(aod->GetHeader());
    
    if(!ah) {
            AliFatal("Not a standard AOD");
            return;
    }
	hdr->SetESDFileName( ah->GetESDFileName() );
	hdr->SetEventNumberESDFile( ah->GetEventNumberESDFile() );
}

//______________________________________________________________________________
Int_t AliJFilter::GetSuperModuleNumber(bool isemcal, AliVCluster *cluster, AliVCaloCells *cells, Int_t absId)
{
	//get super module number 
	if(isemcal){
		Int_t absIdMax  = -1, iSM =-1, ieta = -1, iphi = -1;
		Bool_t shared = kFALSE;
		fEMCALRecoUtils->GetMaxEnergyCell(fEMCALGeometry, cells, cluster, absIdMax,  iSM, ieta, iphi, shared);

		if(iSM < 0 || iphi < 0 || ieta < 0 ) 
		{
			AliFatal(Form("Negative value for super module: %d, or cell ieta: %d, or cell iphi: %d, check EMCAL geometry name\n",
						iSM,ieta,iphi));
		}

		return iSM ;

	} else {
		Int_t    relId[4];
		if ( absId >= 0) {
			fPHOSGeom->AbsToRelNumbering(absId,relId);
			fPHOSGeom->AbsToRelNumbering(absId,relId);
			return relId[0]-1; 
		} else return -1;
	}//PHOS

	return -1;
}

//______________________________________________________________________________
Double_t * AliJFilter::GetCellsAmplitude( bool isemcal, AliVCluster *cluster, AliVCaloCells *emCells, AliVCaloCells *phoCells )
{
	// cell amplitude reader
	Int_t iCell, nCell;
	UShort_t *cellAddrs;
	Double_t *amps;

	// get cluster cells
	nCell = cluster->GetNCells();

	amps = new Double_t[nCell];

	// get the cell addresses
	cellAddrs = cluster->GetCellsAbsId();

	// get the cell amplitudes
	for( iCell = 0; iCell < nCell; iCell++ ){
		if( isemcal )
			amps[iCell] = emCells->GetCellAmplitude( cellAddrs[iCell] );
		else
			amps[iCell] = phoCells->GetCellAmplitude( cellAddrs[iCell] );

	}

	return amps;
}

//_____________________________________________________________________________

UInt_t AliJFilter::ConvertTriggerMask(){

	//convert alice trigger mask to jcorran trigger mask
	UInt_t triggerMaskJC=0;
	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kMB){
		// minimum bias TBit 0 
		triggerMaskJC |= (1<<kMinBiasTriggerBitJCorran); 
		fAliJRunHeader->SetActiveTriggersJCorran( kMinBiasTriggerBitJCorran, "MinBiasTriggerBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kHighMult){
		//high multiplicity trigger TBit 1 
		triggerMaskJC |= (1<<kHighMultTriggerBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kHighMultTriggerBitJCorran,"HighMultTriggerBitJCorran");
	}

	if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
				->IsEventSelected() & AliVEvent::kEMC1) ||
			(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			 ->IsEventSelected() & AliVEvent::kEMC7 )){
		//EMCAL L0   TBit2
		triggerMaskJC |= (1<<kEmc0TriggerBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kEmc0TriggerBitJCorran,"Emc0TriggerBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kEMCEGA){
		//EMCAL Gamma TBit3
		triggerMaskJC |= (1<<kEmc1GammaTriggerBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kEmc1GammaTriggerBitJCorran,"Emc1GammaTriggerBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kEMCEJE){
		//EMCAL JET TBit4
		triggerMaskJC |= (1<<kEmc1JetTriggerBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kEmc1JetTriggerBitJCorran,"Emc1JetTriggerBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kCentral){
		//central trigger TBit 5 
		triggerMaskJC |= (1<<kCentralTriggerBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kCentralTriggerBitJCorran,"CentralTriggerBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kSemiCentral){
		//semi-central trigger TBit 6 
		triggerMaskJC |= (1<<kSemiCentralTriggerBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kSemiCentralTriggerBitJCorran,"SemiCentralTriggerBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kFastOnly){
		//semi-central trigger TBit 6 
		triggerMaskJC |= (1<<kFastOnlyBitJCorran);
		fAliJRunHeader->SetActiveTriggersJCorran( kFastOnlyBitJCorran ,"FastOnlyBitJCorran");
	}

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & AliVEvent::kINT7){
		// minimum bias TBit 0 
		triggerMaskJC |= (1<<kINT7TriggerBitJCorran); 
		fAliJRunHeader->SetActiveTriggersJCorran( kINT7TriggerBitJCorran, "INT7TriggerBitJCorran");
	}

	return triggerMaskJC;
}


//______________________________________________________________________________
void AliJFilter::ReadMCTracksFromESD(){ 
	//store MC information from AliStack
	if(!MCEvent()) return;
	AliStack *stack = MCEvent()->Stack();
	if(!stack) return;
	Int_t np    = MCEvent()->GetNumberOfTracks();

	//  AliGenEventHeader* genHeader = fMC->GenEventHeader();
	//  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
	//  Double_t ptHard = 0;
	//  Double_t nTrials = 1; // Trials for MC trigger weigth for real data
	//  nTrials = pythiaGenHeader->Trials();
	//  ptHard  = pythiaGenHeader->GetPtHard();
	//  Int_t nprim = stack->GetNtrack();

	Long64_t ntrack = 0;

	for(Long64_t iTrack = 0; iTrack < np; iTrack++){
		AliMCParticle *track = (AliMCParticle*) MCEvent()->GetTrack(iTrack);
		if(!track){
			Printf("ERROR: Could not receive track %d",(int) iTrack);
			continue;
		}
		Bool_t isPrimary = stack->IsPhysicalPrimary(iTrack);
		if(isPrimary){
			//create a new JMCTrack and fill the track info
			AliJMCTrack *ctrack = new( (*fMCTrackList)[ntrack++] ) AliJMCTrack;

			TParticle *partStack = stack->Particle(iTrack);
			Int_t   pdg  = partStack->GetPdgCode();

			Char_t ch     = (Char_t) partStack->GetPDG()->Charge();
			Int_t label    = track->GetLabel();

			ctrack->SetLabel(label);
			ctrack->SetPdgCode(pdg);
			ctrack->SetPxPyPzE( partStack->Px(), partStack->Py(), partStack->Pz(), partStack->Energy());
			ctrack->SetCharge(ch); 
			ctrack->SetFlag(AliJMCTrack::kPrimary, isPrimary);

			ctrack->SetProductionVertex(partStack->Vx(),partStack->Vy(),partStack->Vz());
		}// loop for al primary tracks
	} 
}

//--------------------------------------------------------------------
void AliJFilter::ReadMCTracksFromAOD(){
	//retreive MC particles from event //FKEFF// 
	if(!AODEvent()) return;  TClonesArray *mcArray = (TClonesArray*) AODEvent()->
		FindListObject(AliAODMCParticle::StdBranchName());
	if(!mcArray){
		Printf("No MC particle branch found");
		return;
	}

	Long64_t ntrack = 0;
	Long64_t np = mcArray->GetEntriesFast();

	for(Long64_t it = 0; it < np; it++) {
		AliAODMCParticle *track = (AliAODMCParticle*) mcArray->At(it);
		if(!track){
			Error("ReadEventAODMC", "Could not receive particle %d",(int) it);
			continue;
		}
		bool isPrimary = track->IsPhysicalPrimary();
		if(isPrimary){
			//create a new JMCTrack and fill the track info
			AliJMCTrack *ctrack = new ((*fMCTrackList)[ntrack++]) AliJMCTrack;;

			Int_t   pdg  = track->GetPdgCode();

			Char_t ch     = (Char_t) track->Charge();
			Int_t label    = track->GetLabel();

			ctrack->SetLabel(label);
			ctrack->SetPdgCode(pdg);
			ctrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E());
			ctrack->SetCharge(ch);
			ctrack->SetFlag(AliJMCTrack::kPrimary, isPrimary);

			ctrack->SetProductionVertex(track->Xv(),track->Yv(),track->Zv());
      
      // If the particle has no daughters, it must be still alive when hitting the detector
      Int_t nDaughters = track->GetNDaughters();
      Bool_t isFinal = kFALSE;
      if(nDaughters == 0) isFinal = kTRUE;
      ctrack->SetIsFinal(isFinal);
		}
	}

}


//--------------------------------------------------------------------
void AliJFilter::RemapMCLabels(){
	// remaps all MC labels to the new arrays

	Int_t i, j, label, mother0, mother1;
	AliJTrack *track;
	AliJPhoton *cluster;
	// BS AliJCaloCell *cell;
	AliJMCTrack *mctrack;

	// tracks
	for( i = 0; i < fTrackList->GetEntries(); i++ ){
		track = (AliJTrack*)fTrackList->At( i );

		track->SetLabel( fMcMap->At( track->GetLabel() ));
	}

	// clusters
	if( fAliJRunHeader->GetStoreEMCalInfo() ){
		for( i = 0; i < fPhotonList->GetEntries(); i++ ){
			cluster = (AliJPhoton*)fPhotonList->At( i );
			for( j = 0; j < cluster->GetNEMCLabel(); j++ ){
				label = cluster->GetEMCLabel( j );
				// no label clusters protection
				if( label >= 0 )
					cluster->SetEMCLabel( j, fMcMap->At( label ));
			}
		}

		/*  BS
		// cells
		for( i = 0; i < fCaloCellList->GetEntries(); i++ ){
		cell = (AliJCaloCell*)fCaloCellList->At( i );
		label = cell->GetMcLabel();
// no label cells protection
if( label >= 0 )
cell->SetMcLabel( fMcMap->At( cell->GetMcLabel() ));
}
*/
}

// MC particles
for( i = 0; i < fMCTrackList->GetEntries(); i++ ){
	mctrack = (AliJMCTrack*)fMCTrackList->At( i );

	mother0 = mctrack->GetMother( 0 );
	mother1 = mctrack->GetMother( 1 );

	if( mother0 >= 0 )
		mother0 = fMcMap->At( mother0 );
	if( mother1 >= 0 )
		mother1 = fMcMap->At( mother1 );

	mctrack->SetMother( mother0, mother1 );
}
}

//--------------------------------------------------------------------


void AliJFilter::PrintOut() const {
	//AliJRunHeader * RunInfo = fAliJRunHeader;
}

//********************************************
//    UTILS
//********************************************
void AliJFilter::AddList(const char* aname, const char* cname, TClonesArray **obj, int nlist){
	*obj = new TClonesArray(cname, nlist);
	(*obj)->SetName(aname);
	(*obj)->SetOwner();
}

