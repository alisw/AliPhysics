#include "AliJFilter.h" 
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMCEvent.h" 
#include "AliStack.h" 
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODTrack.h"
#include "AliAnalysisFilter.h"
#include "AliAODVertex.h" 
#include "AliAODTracklets.h" 
#include "AliAODPid.h" 
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h" 
#include "AliMultSelection.h"
#include "AliAODTracklets.h"
#include "AliMultiplicity.h"
#include "AliJConst.h"
#include "AliDAQ.h"
#include "AliExternalTrackParam.h"
#include "AliHeader.h" 
#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"

#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAnalysisUtils.h"
#include "AliJRunTable.h"


ClassImp(AliJFilter);

//______________________________________________________________________________
AliJFilter::AliJFilter() :   
	TNamed(),
	fIsRealOrMC(0),
	fStoreEventPlaneSource(0),
	fOADBPath(),
	fTrackThreshold(0),
	fEventSuccess(0),
	fMcMap(0),
	fTrackList(0),
	fMCTrackList(0x0),
	fHeaderList(0x0),
	fRunInfoList(0x0),
	fPIDResponse(0x0),
	fPIDCombined(0x0),
	fAliJRunHeader(0x0),
	fAnaUtils(0x0),
	fMyTask(0x0),
	pfOutlierLowCut(0x0),
	pfOutlierHighCut(0x0),
	fFirstEvent(kTRUE),
	fRunTable(0)
{
	//Default constructor
}

//______________________________________________________________________________
AliJFilter::AliJFilter(const char *name,AliAnalysisTaskSE *task):
	TNamed(name,name), 
	fIsRealOrMC(0),
	fStoreEventPlaneSource(0),
	fOADBPath(),
	fTrackThreshold(0),
	fEventSuccess(0),
	fMcMap(0),
	fTrackList(0),
	fMCTrackList(0x0),
	fHeaderList(0x0),
	fRunInfoList(0x0),
	fPIDResponse(0x0),
	fPIDCombined(0x0),
	fAliJRunHeader(0x0),
	fAnaUtils(0x0),
	fMyTask(0x0),
	fFirstEvent(kTRUE),
	fRunTable(0)
{
	// Constructor
	if(task->DebugLevel() > 5) cout << "---- AliJFilter Constructor ----"<<endl;
	// Cut for 15o period
	pfOutlierLowCut = new TF1("fLowCut","[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);
	pfOutlierHighCut = new TF1("fHighCut","[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);

	pfOutlierLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
	pfOutlierHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);

}

//____________________________________________________________________________
AliJFilter::AliJFilter(const AliJFilter& ap) :
	TNamed(ap.GetName(), ap.GetTitle()),
	fIsRealOrMC(ap.fIsRealOrMC),
	fStoreEventPlaneSource(ap.fStoreEventPlaneSource),
	fOADBPath(ap.fOADBPath),
	fTrackThreshold(ap.fTrackThreshold),
	fEventSuccess(ap.fEventSuccess),
	fMcMap(ap.fMcMap),
	fTrackList(ap.fTrackList),
	fMCTrackList(ap.fMCTrackList),
	fHeaderList(ap.fHeaderList),
	fRunInfoList(ap.fRunInfoList),
	fPIDResponse(ap.fPIDResponse),
	fPIDCombined(ap.fPIDCombined),
	fAliJRunHeader(ap.fAliJRunHeader),
	fAnaUtils(ap.fAnaUtils),
	fMyTask(ap.fMyTask),
	fFirstEvent(kTRUE),
	fRunTable(ap.fRunTable)

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
	delete fHeaderList;
	delete fAliJRunHeader;
	delete fRunInfoList;
	delete fPIDResponse;
	delete fPIDCombined;
	delete fAnaUtils;


}

//________________________________________________________________________

void AliJFilter::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fMyTask->DebugLevel() > 1) printf("AliJFilter::UserCreateOutPutData() \n");

        pfOutlierLowCut = new TF1("fLowCut","[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);
        pfOutlierHighCut = new TF1("fHighCut","[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);

        pfOutlierLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
        pfOutlierHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);

	//== RUN HEADER
	cout<<"TEST2 "<<fAliJRunHeader<<endl;
	if(!fAliJRunHeader) fAliJRunHeader = new AliJRunHeader();
	fRunInfoList  = new TList();
	fRunInfoList->SetName("RunInfoList");
	//fRunInfoList->SetOwner();
	fRunInfoList->Clear();
	fRunInfoList->Add(fAliJRunHeader);

	fAnaUtils = new AliAnalysisUtils();
	fAnaUtils->SetUseOutOfBunchPileUp( kTRUE );
	fMcMap = new TArrayI();

	DEBUG( 5, 0, Form("IsMC() = %d",IsMC())  );
	//=== Set Tree and TClonesArray
	//== TRACKS
	AddList("AliJTrackList", "AliJTrack", &fTrackList, 1000);
	if( IsMC() ) 
		AddList("AliJMCTrackList", "AliJMCTrack", &fMCTrackList, 1000);
	//== Event Header
	AddList("AliJEventHeaderList", "AliJEventHeader", &fHeaderList, 1000);

	//== EventPlane SRC
	cout << "Add(fAliJRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJFilter::UserExec(Option_t* /*option*/) 
{
	// user loop
	AliJRunHeader *runh = fAliJRunHeader;
	Bool_t hasGoodTrack;

	fEventSuccess = kFALSE;

	// Processing of one event
	DEBUG( 5, 1, "------- AliJFilter Exec-------" );
	if(!((fMyTask->Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  fMyTask->Entry())); 

	//=== Init Variables
	fTrackList->Clear();
	if( IsMC() ){
		fMCTrackList->Clear();
	}

	fHeaderList->Clear();

	hasGoodTrack = kTRUE;

	//=== CHECK ESD, AOD, MC event
	if( !Event() ) return;
	if( !AODEvent() ) return; 


	// pileup rejection
	if( runh->IsPP() && fAnaUtils->IsPileUpEvent( Event() )) // Only APPLY for P+P
		return;

	if( fFirstEvent ) {
		fRunTable = & AliJRunTable::GetSpecialInstance();
		fRunTable->SetRunNumber( Event()->GetRunNumber() );
		fFirstEvent = kFALSE;
	}
	if(!IsGoodEvent( AODEvent() )) // check now for 15o data event selection
		return; 

	//--------------------------------------------------------------- 
	// RUN Header
	//--------------------------------------------------------------- 
	if(!runh->GetRunNumber()){ //new run has started : I suppose no change of run in process
		runh->SetRunNumber( Event()->GetRunNumber() );
		if( fFirstEvent ) {
			fRunTable = & AliJRunTable::GetSpecialInstance();
			fRunTable->SetRunNumber( Event()->GetRunNumber() );
			fFirstEvent = kFALSE;
		}
		//==== General ====//
		cout << "Run # = "<< AODEvent()->GetRunNumber() << endl;
		runh->SetRunNumber( AODEvent()->GetRunNumber() );
		runh->SetL3MagnetFieldIntensity( AODEvent()->GetMagneticField() );
		runh->SetCurrentL3( AODEvent()->GetMagneticField()*30000.0/5.00668 );
		runh->SetCurrentDip( AODEvent()->GetMuonMagFieldScale()*6000.0 );
		runh->SetUniformBMap( kFALSE ); // TODO is this?
		cout << "Add(fAliJRunHeader) is done =============" << endl;
	}

	//--------------------------------------------------------------- 
	// EventHeader and read Others
	//--------------------------------------------------------------- 
	DEBUG( 5, 1, "\t------- Start READ AOD " );
	ReadAODHeader( AODEvent() );
	ReadAODTracks( AODEvent() );
	if( IsMC() ){
		ReadMCTracksFromAOD();
		//RemapMCLabels();
	}

	if( hasGoodTrack ){
		//=== TODO : need this?
		//AliAODHandler* outputHandler = 
		//		(AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
		//	outputHandler->SetFillAOD(kTRUE);
		//	outputHandler->SetFillExtension(kTRUE);
		fEventSuccess = kTRUE;
	} else {
		fTrackList->Clear();
		if( IsMC() ){
			fMCTrackList->Clear();
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
	fHeaderList->Clear();

	// Processing when the event loop is ended
	cout<<"PWG4JCORRAN Analysis DONE !!"<<endl; 

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
		//double const * pid = track->PID();
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
AliJEventHeader* AliJFilter::ReadCommonHeader(AliAODEvent *event){
	//Read the AliVEvent and fill the list of AliJEventHeader containers
	//create a header and fill it
	AliJEventHeader *hdr = new( (*fHeaderList)[fHeaderList->GetEntriesFast()] ) AliJEventHeader;

	float fcent = -999;
	// centrality
	if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
		AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
		fcent = pms->GetMultiplicityPercentile("V0M");
	} else {
		fcent = -1;
		//cout<<"warning: centrality unavailable";
	}

	hdr->SetCentrality( fcent );
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


//--------------------------------------------------------------------
void AliJFilter::ReadMCTracksFromAOD(){
	//retreive MC particles from event //FKEFF// 
	if(!AODEvent()) { 
		DEBUG(5, 0, "No AODEvent"); 
		return;  
	}
	TClonesArray *mcArray = (TClonesArray*) AODEvent()->FindListObject(AliAODMCParticle::StdBranchName());
	if(!mcArray){
		Printf("No MC particle branch found");
		return;
	}

	Long64_t ntrack = 0;
	Long64_t np = mcArray->GetEntriesFast();
	DEBUG(5, 1, Form("MC all::NumberOfTracks = %d",(int) np) );

	for(Long64_t it = 0; it < np; it++) {
		AliAODMCParticle *track = (AliAODMCParticle*) mcArray->At(it);
		if(!track){
			DEBUG(5, 0, Form("eadEventAODMC Could not receive particle %d",(int) it) );
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
	DEBUG(5, 1, Form("MC primary::NumberOfTracks = %d",fMCTrackList->GetEntriesFast()) );

}


//--------------------------------------------------------------------
void AliJFilter::RemapMCLabels(){
	// remaps all MC labels to the new arrays

	Int_t i, mother0, mother1;
	AliJTrack *track;
	// BS AliJCaloCell *cell;
	AliJMCTrack *mctrack;

	// tracks
	for( i = 0; i < fTrackList->GetEntries(); i++ ){
		track = (AliJTrack*)fTrackList->At( i );

		track->SetLabel( fMcMap->At( track->GetLabel() ));
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
// Perodic specific event selection before fill the events
//________________________________________________________________________
Bool_t AliJFilter::IsGoodEvent(AliAODEvent *event) {
	// RunTable to sepecify the run conditions
	// Taken from AliJFFlucTask for basic event selection cuts need to add few more
	if(fRunTable->GetRunNumberToPeriod(fRunTable->GetRunNumber()) == AliJRunTable::kLHC15o){
		const AliVVertex* vtTrc = event->GetPrimaryVertex();
		const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();
		double covTrc[6],covSPD[6];
		vtTrc->GetCovarianceMatrix(covTrc);
		vtSPD->GetCovarianceMatrix(covSPD);
		double dz = vtTrc->GetZ()-vtSPD->GetZ();
		double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
		double errTrc = TMath::Sqrt(covTrc[5]);
		double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
		if(TMath::Abs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20)
			return kFALSE;
		AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
		if(!pms){
			AliError("MultSelection unavailable.");
			return kFALSE;
		}

		Float_t v0mcent = pms->GetMultiplicityPercentile("V0M");
		Float_t cl0cent = pms->GetMultiplicityPercentile("CL0");
		return kTRUE;
		if(cl0cent < pfOutlierLowCut->Eval(v0mcent) || cl0cent > pfOutlierHighCut->Eval(v0mcent))
			return kFALSE;
	}
	return kTRUE;
}

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

