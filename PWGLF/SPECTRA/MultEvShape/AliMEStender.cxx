#include <TROOT.h>
#include <TString.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TTreeStream.h>
#include <TGeoGlobalMagField.h>

#include <AliAnalysisManager.h>
#include <AliGeomManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDtrackCuts.h>
#include <AliPIDCombined.h>
#include <AliESDInputHandler.h>
#include <AliESDInputHandlerRP.h>
#include <AliMCEventHandler.h>

#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliESDRun.h>
#include <AliESDtrack.h>
#include <AliStack.h>
#include <AliMCEvent.h>
#include <AliMCParticle.h>

#include <AliMultiplicity.h>
#include <AliCentrality.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliAnalysisUtils.h>
#include <AliPPVsMultUtils.h>

#include "AliMESbaseTask.h"
#include "AliMEStender.h"
#include "AliMESeventInfo.h"
#include "AliMEStrackInfo.h"

ClassImp(AliMEStender::AliMESconfigTender)
ClassImp(AliMEStender)

//________________________________________________________________________
AliMEStender::AliMESconfigTender::AliMESconfigTender() : TObject(), fTrackCuts(0), fEventCuts(0), fPIDpriors(0) {}

//________________________________________________________________________
void AliMEStender::AliMESconfigTender::Print(Option_t *) const
{
  // Dump config info to stdout
  printf("MES TENDER CONFIGURATION\n   Event cuts : ");
  switch(fEventCuts){
    case kNoEC: printf("No\n"); break;
    case kStandard: printf("Trigger[MB, HM], Vertex[Yes]\n"); break;
    default: printf("Not defined [%d]\n", fEventCuts); break;
  }
  printf("   Track cuts : ");
  switch(fTrackCuts){
    case kNoTC: printf("No\n"); break;
    case kStandardITSTPCTrackCuts2010: printf("StandardITSTPCTrackCuts2010\n"); break;
    default: printf("Not defined [%d]\n", fTrackCuts); break;
  }
  printf("   PID priors : ");
  switch(fPIDpriors){
    case kNoPP: printf("Flat Priors\n"); break;
    case kTPC: printf("DefaultTPCPriors\n"); break;
	case kIterative: printf("LHC10d Iterative Priors\n"); break;
    default: printf("Not defined [%d]\n", fPIDpriors); break;
  }

  printf("\n");
}

//________________________________________________________________________
Bool_t AliMEStender::ConfigTask(AliMESconfigTender::EMESconfigEventCuts ec,
                                AliMESconfigTender::EMESconfigTrackCuts tc,
                                AliMESconfigTender::EMESconfigPIDpriors pp)
{
  // configure tender for current run
  fConfig.fEventCuts = ec;
  fConfig.fTrackCuts = tc;
  fConfig.fPIDpriors = pp;
  fConfig.Print();
  return kTRUE;
}

//________________________________________________________________________
AliMEStender::AliMEStender()
  : AliAnalysisTaskSE()
  ,fConfig()
  ,fTrackFilter(NULL)
  ,fPIDcomb(NULL)
  ,fTracks(NULL)
  ,fEvInfo(NULL)
  ,fMCtracks(NULL)
  ,fMCevInfo(NULL)
  ,fUtils(NULL)
{
  //
  // Constructor
  //
}

//________________________________________________________________________
AliMEStender::AliMEStender(const char *name)
  : AliAnalysisTaskSE(name)
  ,fConfig()
  ,fTrackFilter(NULL)
  ,fPIDcomb(NULL)
  ,fTracks(NULL)
  ,fEvInfo(NULL)
  ,fMCtracks(NULL)
  ,fMCevInfo(NULL)
  ,fUtils(NULL)
{
  //
  // Constructor
  //
  DefineOutput(AliMESbaseTask::kQA,          TList::Class());
  DefineOutput(AliMESbaseTask::kEventInfo+1, AliMESeventInfo::Class());
  DefineOutput(AliMESbaseTask::kTracks+1,    TObjArray::Class());
}

//________________________________________________________________________
void AliMEStender::SetMCdata(Bool_t mc)
{
  // prepare task for MC processing
  SetBit(kMCdata, mc);
  if(mc){
    DefineOutput(AliMESbaseTask::kMCeventInfo+1, AliMESeventInfo::Class());
    DefineOutput(AliMESbaseTask::kMCtracks+1,    TObjArray::Class());
  }
}

//________________________________________________________________________
AliMEStender::~AliMEStender()
{
  //
  // Destructor
  //
  if(fTrackFilter) delete fTrackFilter;
  if(fPIDcomb) delete fPIDcomb;

  if(fEvInfo) delete fEvInfo;
  if(fTracks->GetEntries()){ fTracks->Delete(); delete fTracks; }
  if(fMCevInfo) delete fMCevInfo;
  if(fMCtracks){ fMCtracks->Delete(); delete fMCtracks; }

  if(fUtils) delete fUtils;

  if(DebugLevel()) AliMESbaseTask::CloseDebugStream();

}

//________________________________________________________________________
void AliMEStender::UserCreateOutputObjects()
{
  // Build user objects
  BuildQAHistos(); PostData(AliMESbaseTask::kQA, fHistosQA);

// ------- track cuts
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts *lTrackCuts(NULL);
  switch(fConfig.fTrackCuts){
  case AliMESconfigTender::kStandardITSTPCTrackCuts2010:
    lTrackCuts = new AliESDtrackCuts("std10TC", "Standard 2010");
    lTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,0);
    fTrackFilter->AddCuts(lTrackCuts);
    break;
  case AliMESconfigTender::kNoTC:
  default:
    AliDebug(2, "No track cuts selected");
    break;
  }

  // PID priors
//   fPIDcomb = new AliPIDCombined();
//   fPIDcomb->SetSelectedSpecies(AliPID::kSPECIES);
/*
  switch(fConfig.fPIDpriors){
	case AliMESconfigTender::kTPC:
		// default aliroot priors
		fPIDcomb->SetDefaultTPCPriors();
		break;
	case AliMESconfigTender::kIterative:
	{  // data priors identified @ 15.04.2015 by Cristi for LHC10d
		AliInfo("Loading iterative data priors ...");
		fPIDcomb->SetPriorDistribution(AliPID::kMuon, fPriorsDist[0]);
		fPIDcomb->SetPriorDistribution(AliPID::kElectron, fPriorsDist[0]);
		fPIDcomb->SetPriorDistribution(AliPID::kPion, fPriorsDist[1]);
		fPIDcomb->SetPriorDistribution(AliPID::kKaon, fPriorsDist[2]);
		fPIDcomb->SetPriorDistribution(AliPID::kProton, fPriorsDist[3]);
		AliInfo("Done loading iterative data priors.");
		break;
	}
	case AliMESconfigTender::kNoPP:
	{ // flat priors
		fPIDcomb->SetEnablePriors(kFALSE);  // FLAT priors
		break;
	}
	default:
		AliDebug(2, "No PID priors selected");
		break;
  }
*/

  fTracks = new TObjArray(200);
  fTracks->SetOwner(kTRUE);
  fEvInfo = new AliMESeventInfo;

  fUtils = new AliPPVsMultUtils();

  if(!HasMCdata()) return;
  fMCtracks = new TObjArray(200);
  fMCtracks->SetOwner(kTRUE);
  fMCevInfo = new AliMESeventInfo;

}

#include "AliGRPManager.h"
//________________________________________________________________________
void AliMEStender::UserExec(Option_t */*opt*/)
{
  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    AliError("ESD event not available");
    return;
  }

//   AliMESeventInfo *fEvInfo   = dynamic_cast<AliMESeventInfo*>(GetOutputData(AliMESbaseTask::kEventInfo+1));
  if(!fEvInfo){
    AliError("REC event info missing. Processing skipped");
    return;
  }
//   TObjArray *fTracks   = dynamic_cast<TObjArray*>(GetOutputData(AliMESbaseTask::kTracks+1));
  if(!fTracks){
    AliError("REC track array missing. Processing skipped");
    return;
  }
  fEvInfo->Clear("");
  fTracks->Delete();
  // check QA container
  if(!fHistosQA){
    AliError("No QA container defined.");
    return;
  }

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  AliPIDResponse *pidResponse = inputHandler->GetPIDResponse();
  if (!pidResponse) AliFatal("This Task needs the PID response attached to the inputHandler");

  //init magnetic field
  if(!TGeoGlobalMagField::Instance()->GetField() && !TGeoGlobalMagField::Instance()->IsLocked()){
    AliGRPManager grpManager;
    if(!grpManager.ReadGRPEntry()) AliError("Cannot get GRP entry for magnetic field");
    if(!grpManager.SetMagField()) AliError("Problem with magnetic field setup");
  }


  ((TH1*)fHistosQA->At(kEfficiency))->Fill(0);  // all events
  if( AliPPVsMultUtils::IsMinimumBias(fESD) ) {
	  ((TH1*)fHistosQA->At(kEfficiency))->Fill(1);  // events after Physics Selection (for MB normalisation to INEL)
  }
  if( !AliPPVsMultUtils::IsEventSelected(fESD) ) {
	  return;
  }
  ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);  // analyzed events


  // TRIGGER SELECTION
  // MB & HM triggers
  Bool_t triggerMB = (inputHandler->IsEventSelected()& AliVEvent::kMB),
         triggerHM = (inputHandler->IsEventSelected()& AliVEvent::kHighMult);
  if(!triggerHM && !triggerMB){
    AliDebug(2, "Miss trigger");
//     ((TH1*)fHistosQA->At(kEfficiency))->Fill(1);
//     return;
  }
  if(triggerMB) fEvInfo->SetTriggerMB();
  if(triggerHM) fEvInfo->SetTriggerHM();

  // vertex selection
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()>=1){ fEvInfo->SetVertex(); fEvInfo->SetVertexGlob(); }
  else {    // try SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()>=1) fEvInfo->SetVertex();
    else {
      AliDebug(2, "Miss vertex");
//       ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);
//       return;
    }
  }
/*
  if(!AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(fESD)){
	  ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);
// 	  return;
  }
  if(!AliPPVsMultUtils::IsINELgtZERO(fESD)){
	  ((TH1*)fHistosQA->At(kEfficiency))->Fill(2);
// 	  return;
  }

// 	((TH1*)fHistosQA->At(kEfficiency))->Fill(0);
*/

  fEvInfo->SetVertexZ(vertex->GetZ());

  // pile-up selection
  if(fESD->IsPileupFromSPDInMultBins()) fEvInfo->SetPileUp();

  // multiplicity
  AliESDtrackCuts *tc(NULL);
  if((tc = dynamic_cast<AliESDtrackCuts*>(fTrackFilter->GetCuts()->At(0)/*FindObject("std10TC")*/))){
    //MakeMultiplicityESD(fESD, "Combined")
	fEvInfo->SetMultiplicity(AliMESeventInfo::kComb, tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8));
	// fEvInfo->SetMultiplicity(AliMESeventInfo::kComb, AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD));
    // MakeMultiplicityESD(fESD, "Global")
    fEvInfo->SetMultiplicity(AliMESeventInfo::kGlob08, tc->CountAcceptedTracks(fESD));
	// V0M
	fEvInfo->SetMultiplicity(AliMESeventInfo::kV0M, fUtils->GetMultiplicityPercentile(fESD, "V0M"));
	// Combined multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
	fEvInfo->SetMultiplicity(AliMESeventInfo::kComb0408, (tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8) - tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.4)));
  } else {
    AliWarning("No track cuts defined. No multiplicity computed for REC data.");
    fTrackFilter->GetCuts()->ls();
  }


  Double_t val[7] = {0.};
  THnSparse *H(NULL);
  H = (THnSparse*)fHistosQA->At(kTrkInfo);
  AliMEStrackInfo *tmes(NULL);

  // Int_t counter = 0;
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {

    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    if(!(fTrackFilter->IsSelected(track)) && !DebugLevel() ) continue;
    tmes = new AliMEStrackInfo(track, pidResponse, fPIDcomb);

    // TOF matching (old version - Alex)
// 	if(pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk) tmes->SetTOFmisProb(fPIDcomb->GetTOFmismatchProb());
// 	printf("mismatch prob from AliMEStender = %g\n", fPIDcomb->GetTOFmismatchProb());


    // TOF matching (Cristi)
	Int_t fTOFout;
	if ((!(track->GetStatus()&AliESDtrack::kTOFout))==0) {fTOFout=1;}else {fTOFout=0;}
	Int_t ftime;
	if ((!(track->GetStatus()&AliESDtrack::kTIME))==0) {ftime=1;}else {ftime=0;}
	Double_t flength;
	flength=track->GetIntegratedLength();
	Double_t ftimetof;
	ftimetof=track->GetTOFsignal();
	Double_t inttime[5];
	track->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
	Double_t fexptimepi;
	Double_t fexptimeka;
	Double_t fexptimepr;
	fexptimepi=inttime[2];
	fexptimeka=inttime[3];
	fexptimepr=inttime[4];

	// Barbara
	if((fTOFout==1)&&(ftime==1)&&(flength>350)&&(ftimetof>10000)&&(fexptimepi>10000)&&(fexptimeka>10000)&&(fexptimepr>10000)&&(ftimetof<80000)){
		tmes->SetTOFmisProb(1);
	}
	else{
		tmes->SetTOFmisProb(0);
	}

    // counter ++;

    // fill tracks QA
    val[0] = tmes->Pt();
    val[1] = tmes->Eta();
    val[2] = tmes->Phi();
    val[3] = tmes->Rv();
    val[4] = tmes->Zv();
    const AliMEStrackInfo::AliMESpid* pid=tmes->GetPID();
    if(pid){
      val[5] = pid->GetRaw(AliMEStrackInfo::kTPC);
      val[6] = pid->GetRaw(AliMEStrackInfo::kTOF);
    }
    if(H) H->Fill(val);
    fTracks->Add(tmes);

    // printf("tender: counter = %i eta = %f \t pT = %g\n", counter, val[1], val[0]);
  }

  // printf("tender: combined08 = %g \t old combined08 = %i \t global08 = %g\n", fEvInfo->GetMultiplicity(AliMESeventInfo::kComb), (tc->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8)), fEvInfo->GetMultiplicity(AliMESeventInfo::kGlob08));
  // printf("counter = %i \t fTracks = %i\n\n", counter, fTracks->GetEntries());
  // printf("event index = %i\n", fESD->GetEventNumberInFile());

  // leading particle
  fEvInfo->FindLeadingParticle(fTracks);
  // shape
  fEvInfo->MakeDirectivity(fTracks);
  fEvInfo->MakeSphericity(fTracks);
  // fill event QA
  H = (THnSparse*)fHistosQA->At(kEvInfo);
  val[0] = fEvInfo->GetVertexZ();
  val[1] = fEvInfo->GetMultiplicity(AliMESeventInfo::kComb);
  val[2] = fEvInfo->GetEventShape()->GetDirectivity(kTRUE);
  val[3] = fEvInfo->GetEventShape()->GetDirectivity(kFALSE);
  if(H) H->Fill(val);



  // fill debug
  if(DebugLevel()>0){
    Int_t run = fESD->GetRunNumber();
    (*AliMESbaseTask::DebugStream()) << "evInfo"
      <<"run=" << run
      <<"ev.=" << fEvInfo
      << "\n";
  }

  PostData(AliMESbaseTask::kEventInfo+1, fEvInfo);
  PostData(AliMESbaseTask::kTracks+1, fTracks);

  //____ _________________________________
  if(!HasMCdata()) return;
//   AliMESeventInfo *fMCevInfo   = dynamic_cast<AliMESeventInfo*>(GetOutputData(AliMESbaseTask::kMCeventInfo+1));
  if(!fMCevInfo){
    AliError("MC event info missing. MC processing skipped");
    return;
  }
//   TObjArray *fMCtracks         = dynamic_cast<TObjArray*>(GetOutputData(AliMESbaseTask::kMCtracks+1));
  if(!fMCtracks){
    AliError("MC track array missing. MC processing skipped");
    return;
  }
  fMCevInfo->Clear("");
  fMCtracks->Delete();

  AliMCEvent *fMC = dynamic_cast<AliMCEvent*>(MCEvent());
  if (!fMC) {
    AliError("MC event not available.");
    return;
  }

  AliStack *fMCStack = fMC->Stack();
  if(!fMCStack){
    AliError("MC stack not available.");
    return;
  }

  // multiplicity
  // multiplicity for eta (-0.8, 0.8)
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kGlob08, MakeMultiplicityMC(fMC));
  // multiplicity for eta (-0.8,-0.4) & (0.4, 0.8)
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kComb0408, MakeMultiplicity0408MC(fMC));
  // multiplicity for eta (-3.7,-1.7) & (2.8, 5.1)  -> V0M
  fMCevInfo->SetMultiplicity(AliMESeventInfo::kV0M, MakeMultiplicityV0MMC(fMC));

  memset(val, 0, 7*sizeof(Double_t));
  H = (THnSparse*)fHistosQA->At(kMCtrkInfo);
  AliMCParticle *particle(NULL); AliMEStrackInfo *tmesRec(NULL);
  for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++) {

    if(!(particle  = dynamic_cast<AliMCParticle*>(fMC->GetTrack(ipart)))) {
      AliWarning("MC particle pointer is null !!!");
      continue;
    }
    if(particle->E()-TMath::Abs(particle->Pz()) < 0.){
      AliWarning("pz > E !!");
      continue;
    }
    if(TMath::Abs(particle->Charge()) < 3) continue;
    tmes = new AliMEStrackInfo(particle, fMCStack);
    fMCtracks->AddLast(tmes);

    // fill tracks QA
    val[0] = tmes->Pt();
    val[1] = tmes->Eta();
    val[2] = tmes->Phi();
    val[3] = tmes->Rv();
    val[4] = tmes->Zv();
    const AliMEStrackInfo::AliMESpid* pid=tmes->GetPID();
    if(pid){
      const Double_t *prob = pid->GetProb(AliMEStrackInfo::kITS);
      for(Int_t is(0); is<AliPID::kSPECIES; is++) if(prob[is]>0.){ val[5] = is; break;}
    }
    val[6] = 0.;
    if(H) H->Fill(val);

    // define matching with ESD track array
    for(Int_t iesd(0); iesd<fTracks->GetEntries(); iesd++){
      if(!( tmesRec = (AliMEStrackInfo*)fTracks->At(iesd))) continue;
      if(tmesRec->GetLabel()!=ipart) continue;
      tmesRec->SetLabel(fMCtracks->GetEntries()-1);
      tmes->SetLabel(iesd);
    }
  }

  // leading particle
  fMCevInfo->FindLeadingParticle(fMCtracks);
  // shape
  fMCevInfo->MakeDirectivity(fMCtracks);
  fMCevInfo->MakeSphericity(fMCtracks);
  // fill event QA
  H = (THnSparse*)fHistosQA->At(kMCevInfo);
  val[0] = 0.;
  val[1] = fMCevInfo->GetMultiplicity(AliMESeventInfo::kGlob08);
  val[2] = fMCevInfo->GetEventShape()->GetDirectivity(kTRUE);
  val[3] = fMCevInfo->GetEventShape()->GetDirectivity(kFALSE);
  if(H) H->Fill(val);


  // fill debug
  if(DebugLevel()>0){
	  (*AliMESbaseTask::DebugStream()) << "evInfoMC"
      <<"ev.=" << fMCevInfo
      << "\n";
  }
  AliDebug(2, Form("Tracks REC[%d] MC[%d]", fTracks->GetEntries(), fMCtracks?fMCtracks->GetEntries():0));

  PostData(AliMESbaseTask::kQA, fHistosQA);
  PostData(AliMESbaseTask::kMCeventInfo+1, fMCevInfo);
  PostData(AliMESbaseTask::kMCtracks+1, fMCtracks);
}

//_____________________________________________________________________
void AliMEStender::SetDebugLevel(Int_t level)
{
  //
  // Init debug stream and register user task
  //
  AliAnalysisTaskSE::SetDebugLevel(level);
  if(level>=1 && !AliMESbaseTask::DebugStream()) AliMESbaseTask::OpenDebugStream();
  if(level>=1) AliMESbaseTask::AddDebugUser(GetName());
}

//_____________________________________________________________________
void AliMEStender::SetPriors(){

	fPIDcomb = new AliPIDCombined();
	fPIDcomb->SetSelectedSpecies(AliPID::kSPECIES);

	switch(fConfig.fPIDpriors){
		case AliMESconfigTender::kTPC:
			// default aliroot priors
			AliInfo("Setting default priors ...");
			fPIDcomb->SetDefaultTPCPriors();
			AliInfo("Done setting default priors.");
			break;
		case AliMESconfigTender::kIterative:
		{  // data priors identified @ 15.04.2015 by Cristi for LHC10d
			AliInfo("Getting iterative data priors from file...");
			TDirectory *cwd(gDirectory);
			TFile *lPriors=TFile::Open("$ALICE_PHYSICS/PWGLF/SPECTRA/MultEvShape/priorsDist_data_LHC10d_newAliroot.root");
			if (lPriors->IsZombie()) {
				AliError("Could not open the priors file");
				return;
			}
			const Char_t *pname[] = {"e", "e", "pi", "K", "p"};
			AliInfo("Setting iterative data priors ...");
			for(Int_t i=0; i<5; i++){
				fPIDcomb->SetPriorDistribution(AliPID::EParticleType(i), (TH1F*)lPriors->Get(Form("priors_%s_final", pname[i])));
				fPIDcomb->GetPriorDistribution(AliPID::EParticleType(i))->SetDirectory(cwd);
			}
			lPriors->Close(); delete lPriors;
			for(Int_t i=0; i<5; i++){
				TH1 *h(fPIDcomb->GetPriorDistribution(AliPID::EParticleType(i)));
				printf("%s -> %s\n", h->GetName(), h->GetDirectory()->GetName());
			}
			fPIDcomb->SetEnablePriors(kTRUE);
			AliInfo("Done setting iterative data priors.");
			break;
		}
		case AliMESconfigTender::kNoPP:
		{ // flat priors
			AliInfo("Setting flat priors ...");
			fPIDcomb->SetEnablePriors(kFALSE);  // FLAT priors
			AliInfo("Done setting flat priors.");
			break;
		}
		default:
			AliWarning("No PID priors selected");
			break;
	}


}

//________________________________________________________________________
Bool_t AliMEStender::PostProcess()
{
  return kTRUE;
}

//________________________________________________________
Bool_t AliMEStender::BuildQAHistos()
{
  // Make QA sparse histos for
  // - task configuration
  // - efficiency
  // - event info
  // - track info

  if(fHistosQA){  // QA histos already created. Check consistency
    if(fHistosQA->GetEntries()==(HasMCdata()?kNqa:kMCevInfo) && dynamic_cast<AliMEStender::AliMESconfigTender*>(fHistosQA->At(kConfig))){
      AliInfo("QA histos already created");
      return kTRUE;
    } else {
      AliError("QA histos wrongly defined.");
      fHistosQA->ls();
      return kFALSE;
    }
  }

  // build QA histos
  fHistosQA = new TList(); fHistosQA->SetOwner(kTRUE);
  fHistosQA->AddAt(&fConfig, kConfig);

  TH1 *hEff = new TH1I("hEff", "Cuts efficiency;cut;entries", 3, -0.5, 2.5);
  TAxis *ax = hEff->GetXaxis();
  ax->SetBinLabel(1, "OK");
  ax->SetBinLabel(2, "Trigger");
  ax->SetBinLabel(3, "Vertex");
//   ax->SetBinLabel(4, "PileUp");

  fHistosQA->AddAt(hEff, kEfficiency);

  TString st;
  THnSparseI *H(NULL);
  if(!(H = (THnSparseI*)gROOT->FindObject("EventInfo"))){
    const Int_t ndim(4);
    const Char_t *cldTitle[ndim] = {"Z_{vertex}", "Multiplicity Combined", "D_{+}", "D_{-}"};
    const Int_t cldNbins[ndim]   = {100, 100, 20, 20};
    const Double_t cldMin[ndim]  = {-15., 0., 0., 0.},
                   cldMax[ndim]  = {15., 100., 1., 1.};
    st = "Event Info;";
    for(Int_t idim(0); idim<ndim; idim++){ st += cldTitle[idim]; st+=";";}
    H = new THnSparseI("EventInfo", st.Data(), ndim, cldNbins, cldMin, cldMax);
  } else H->Reset();
  fHistosQA->AddAt(H, kEvInfo);

  if(!(H = (THnSparseI*)gROOT->FindObject("TrackInfo"))){
    const Int_t ndim(7);
    const Char_t *cldTitle[ndim] = {"p_{t}", "eta", "phi", "DCA_{r}", "DCA_{z}", "dEdx", "p_Beta"};
    const Int_t cldNbins[ndim]   = {  100,     20,    20,     61,        41,      100,      100};
    const Double_t cldMin[ndim]  = {    0.,  -1.0,     0.,    -3.,     -2.0,       0.,       0.},
                   cldMax[ndim]  = {   10.,     1., (2*TMath::Pi()), 3.,  2.,      200.,     1.1};
    st = "Track Info;";
    for(Int_t idim(0); idim<ndim; idim++){ st += cldTitle[idim]; st+=";";}
    H = new THnSparseI("TrackInfo", st.Data(), ndim, cldNbins, cldMin, cldMax);
  } else H->Reset();
  fHistosQA->AddAt(H, kTrkInfo);

  if(HasMCdata()){
    if((H = (THnSparseI*)fHistosQA->At(kEvInfo))){
      H = (THnSparseI*)H->Clone("MCeventInfo");
      H->Reset();
      fHistosQA->AddAt(H, kMCevInfo);
    } else AliError("Missing EventInfo sparse for MC");
    if((H = (THnSparseI*)fHistosQA->At(kTrkInfo))){
      H = (THnSparseI*)H->Clone("MCtrackInfo");
      H->Reset();
      fHistosQA->AddAt(H, kMCtrkInfo);
    } else AliError("Missing TrackInfo sparse for MC");
  }
  AliInfo("Succesfully build QA");
  fHistosQA->ls();
  return kTRUE;
}

// //________________________________________________________________________
// Int_t AliMEStender::MakeMultiplicityESD(AliESDEvent* const esd, const char *opt)
// {
//   if(strcmp(opt, "Combined")==0) return fTrackCuts->GetReferenceMultiplicity(esd, AliESDtrackCuts::kTrackletsITSTPC,0.8);
//   else if(strcmp(opt, "Global")==0) return fTrackCuts->CountAcceptedTracks(esd);
//   else AliError(Form("Wrong option \"%s\"", opt));
//
//   return -1;
// }

//________________________________________________________________________
Int_t AliMEStender::MakeMultiplicityMC(AliMCEvent * const mc)
{
  AliStack *stack(NULL);
  if(!(stack=mc->Stack())) return -1;

//     Int_t nPrim = stack->GetNprimary();
  Int_t charged(0);
  AliMCParticle *particle(NULL);
  for (Int_t ipart=0; ipart<mc->GetNumberOfTracks(); ipart++) {
    if(!( particle = dynamic_cast<AliMCParticle*>(mc->GetTrack(ipart)))) continue;

    if(particle->E()-TMath::Abs(particle->Pz()) < 0.){
      printf(" - E - AliMEStender::MakeMultiplicityMC : pz > E !!\n");
      continue;
    }

    if(!(stack->IsPhysicalPrimary(particle->GetLabel()))) continue;

    //  ---------  Charged  ----------
    if(TMath::Abs(particle->Charge()) < 3) continue;

    if(TMath::Abs(particle->Eta()) > 0.8) continue;

    charged++;

  }//end track loop

  return charged;
}

Int_t AliMEStender::MakeMultiplicity0408MC(AliMCEvent * const mc)
{
	AliStack *stack(NULL);
	if(!(stack=mc->Stack())) return -1;

	//     Int_t nPrim = stack->GetNprimary();
	Int_t charged(0);
	AliMCParticle *particle(NULL);
	for (Int_t ipart=0; ipart<mc->GetNumberOfTracks(); ipart++) {
		if(!( particle = dynamic_cast<AliMCParticle*>(mc->GetTrack(ipart)))) continue;

		if(particle->E()-TMath::Abs(particle->Pz()) < 0.){
			printf(" - E - AliMEStender::MakeMultiplicityMC : pz > E !!\n");
			continue;
		}

		if(!(stack->IsPhysicalPrimary(particle->GetLabel()))) continue;

		//  ---------  Charged  ----------
		if(TMath::Abs(particle->Charge()) < 3) continue;

		if(TMath::Abs(particle->Eta()) > 0.8) continue;
		if(TMath::Abs(particle->Eta()) < 0.4) continue;

		charged++;

	}//end track loop

	return charged;
}

Int_t AliMEStender::MakeMultiplicityV0MMC(AliMCEvent * const mc)
{
	AliStack *stack(NULL);
	if(!(stack=mc->Stack())) return -1;

	//     Int_t nPrim = stack->GetNprimary();
	Int_t charged(0);
	AliMCParticle *particle(NULL);
	for (Int_t ipart=0; ipart<mc->GetNumberOfTracks(); ipart++) {
		if(!( particle = dynamic_cast<AliMCParticle*>(mc->GetTrack(ipart)))) continue;

		if(particle->E()-TMath::Abs(particle->Pz()) < 0.){
			printf(" - E - AliMEStender::MakeMultiplicityMC : pz > E !!\n");
			continue;
		}

		if(!(stack->IsPhysicalPrimary(particle->GetLabel()))) continue;

		//  ---------  Charged  ----------
		if(TMath::Abs(particle->Charge()) < 3) continue;

		if( particle->Eta() < 5.1 &&  particle->Eta() > 2.8 ) charged++;   // V0A
		if( particle->Eta() < -1.7 &&  particle->Eta() > -3.7 ) charged++;   // V0C

	}//end track loop

	return charged;
}
