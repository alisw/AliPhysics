#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObjString.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDTZERO.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliESDTZEROfriend.h"
#include "AliT0CalibAnalysisTask.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliVEvent.h"

ClassImp(AliT0CalibAnalysisTask)

//________________________________________________________________________
AliT0CalibAnalysisTask::AliT0CalibAnalysisTask(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0), 
  fT0OutTree(0x0), fEvent(-99999),   fOrbit(-99999), fBC(-99999),
  fTrackletSPD(-99999), fNcont(-99999),
  fVertex(-99999), fVertexSPD(-99999), 
  fMeanAC(-99999), fMeanA(-99999),fMeanC(-99999), 
  fMultV0A(-99999), fMultV0C(-99999),fTimeV0A(-99999),fTimeV0C(-99999), fSumampA(-99999), fSumampC(-99999),
  ftimestamp(0), fSep2(0),
  fZDCcut(kFALSE), fT0Trigger(0), fpileup(kFALSE), fTrigger(0x0),
  fT0_amplitude(0x0), fT0_time(0x0),
  fcentralityV0M(0), fcentralityZDC(0), fcentralityTRK(0),
  fESDtracks(0),
  fMultiplicity(-99999),
  fTriggerinput(0x0), fZDCbg(kFALSE),
  fTOFtracks(0), fT0tofTrack(0),
  fESDpid(new AliESDpid()), fPFPbit(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliT0CalibAnalysisTask::~AliT0CalibAnalysisTask() 
{
// Destructor
   if (fT0OutTree) delete fT0OutTree ;
}

//_____________________________________________________________________________
Bool_t AliT0CalibAnalysisTask::UserNotify()
{
//
// Calls the mother class Notify()
//

//  AliDebug(AliLog::kDebug,"<-");
  // AliDebug(AliLog::kDebug,"->");

  return AliAnalysisTaskSE::UserNotify();
}

//________________________________________________________________________
void AliT0CalibAnalysisTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
   fOutputList = new TList();
  fOutputList->SetOwner(); // Will delete the histos on cleanup

  fT0OutTree=new TTree("t0tree","None here");
  fT0OutTree->Branch("fNevent", &fEvent);
  fT0OutTree->Branch("fOrbit", &fOrbit);
  fT0OutTree->Branch("fBC", &fBC);
  fT0OutTree->Branch("fNcont", &fNcont, "fNcont/I");
  fT0OutTree->Branch("vertexT0", &fVertex, "fVertex/F");
  fT0OutTree->Branch("vertexSPD", &fVertexSPD, "vertexSPD/F");
  fT0OutTree->Branch("meanAC", &fMeanAC,"meanAC/F");
  fT0OutTree->Branch("meanA", &fMeanA,"meanA/F");
  fT0OutTree->Branch("meanC", &fMeanC,"meanC/F");
  fT0OutTree->Branch("trackletSPD", &fTrackletSPD,"trackletSPD/I");
  fT0OutTree->Branch("TOFtracks", &fTOFtracks);
  fT0OutTree->Branch("ESDtracks", &fESDtracks);
  fT0OutTree->Branch("multV0A", &fMultV0A,"multV0A/F");
  fT0OutTree->Branch("multV0C", &fMultV0C,"multV0C/F");
  fT0OutTree->Branch("timeV0A", &fTimeV0A,"timeV0A/F");
  fT0OutTree->Branch("timeV0C", &fTimeV0C,"timeV0C/F");
  fT0OutTree->Branch("sumampA", &fSumampA);
  fT0OutTree->Branch("sumampC", &fSumampC);
  fT0OutTree->Branch("pileup", &fpileup);
  fT0OutTree->Branch("trigger", &fTrigger);
  fT0OutTree->Branch("triggerinput", &fTriggerinput);
  fT0OutTree->Branch("T0trigger", &fT0Trigger);
  fT0OutTree->Branch("t0tofTrack",  &fT0tofTrack);
  fT0OutTree->Branch("ZDCcut", &fZDCcut);
  fT0OutTree->Branch("ftimestamp", &ftimestamp);
  fT0OutTree->Branch("centralityV0M", &fcentralityV0M);
  fT0OutTree->Branch("centralityZDC", &fcentralityZDC);
  fT0OutTree->Branch("centralityTRK", &fcentralityTRK);
  for ( Int_t i=0; i<24; i++) {
    fT0OutTree->Branch(Form("amp%i", i+1), &famp[i]);
    fT0OutTree->Branch(Form("amp%i_new", i+1), &famp_new[i]);
    fT0OutTree->Branch(Form("time%i", i+1), &ftime[i]);
    fT0OutTree->Branch(Form("fRawTime%i", i+1), &fRawTime[i]);
  } 
  for ( Int_t i=0; i<5; i++) {
    fT0OutTree->Branch(Form("fOrA%i", i+1), &fOrA[i]);
    fT0OutTree->Branch(Form("fOrC%i", i+1), &fOrC[i]);
    fT0OutTree->Branch(Form("fTVDC%i", i+1), &fTVDC[i]);  
  }
  TString pilename[3] = {"T0pileup", "T0background", "T0satellite"};
  for ( Int_t i=0; i<3; i++) 
    fT0OutTree->Branch(pilename[i].Data(), &fT0pileup[i]);
  fT0OutTree->Branch("multEstimator", &fMultiplicity);
  //ZDC background
 fT0OutTree->Branch("zbg", &fZDCbg);
 fESDpid = new AliESDpid();
 //PFP
 fT0OutTree->Branch("bitPFP", &fPFPbit);

  fOutputList->Add(fT0OutTree);
 
   PostData(1, fOutputList);
  
}

//________________________________________________________________________
void AliT0CalibAnalysisTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  fVertex= fVertexSPD =  fMeanAC = fMeanA =  fMeanC =  fTrackletSPD = fNcont = -99999;
  fMultV0A=fMultV0C=fTimeV0A=fTimeV0C=fMultiplicity = -99999;
  
  for (Int_t i=0; i<24; i++)	{
    famp[i] = famp_new[i] = ftime[i] = fRawTime[i] = -9999;
  }
  for ( Int_t i0=0; i0<5; i0++) {
    fTVDC[i0]=-9999;
    fOrC[i0]=-9999;
    fOrA[i0]=-9999;
  }
  
  Float_t orA , orC, tvdc;
  orA = orC = tvdc = -9999;
  fBC = -9999;
  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  if(fESD)
    {
      //    if(fIsSelected) 
      {
	fZDCbg = kFALSE;
	/*
	  AliTriggerAnalysis *fTriggerAnalysis;
	  Bool_t	fZDCbgA = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kZNABG);
	  Bool_t	fZDCbgC =  fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kZNCBG);
	  fZDCbg = (fZDCbgA);
	  cout<<" ZDC back "<<fZDCbgA<<endl;
	*/
	AliESDVZERO* esdV0 = fESD->GetVZEROData();
	const AliMultiplicity *mult = fESD->GetMultiplicity();
	// printf(" VZERO MULT \n");
	AliESDTZERO* tz= (AliESDTZERO*) fESD->GetESDTZERO();
	
	const Double32_t *amp, *time, *mean,*amp_new ;
	fVertex= fVertexSPD =  fMeanAC = fMeanA =  fMeanC = fTrackletSPD = fNcont =-99999;
	fMultV0A=fMultV0C=fTimeV0A=fTimeV0C=fMultiplicity = -99999;
	
	fSumampA=0;
	fSumampC=0;
	fNcont=0;
	TString triggers = fESD->GetFiredTriggerClasses();
	fTrigger.SetString(triggers.Data());
	TString inputtriggers =fESD-> GetHeader()->GetFiredTriggerInputs();
	fTriggerinput.SetString(inputtriggers.Data());
	fEvent=fESD->GetEventNumberInFile();
	fT0Trigger=fESD->GetT0Trig();
	
	fESDtracks=fESD->GetNumberOfTracks();
	ftimestamp=fESD->GetTimeStamp(); // - 1301217882; 
	//	  printf (" timestamp %d \n", ftimestamp);
	
	fMultiplicity =  AliESDtrackCuts::GetReferenceMultiplicity(fESD);
	//	  cout<<" fMultiplicity "<<fMultiplicity<<endl;
	fOrbit=fESD->GetOrbitNumber();
	fBC=fESD->GetBunchCrossNumber();
	fNcont = fESD->GetPrimaryVertexSPD()->GetNContributors();
	if(fNcont>=0 ) {
	  fVertexSPD = fESD->GetPrimaryVertex()->GetZ();
	}
	fTrackletSPD = mult->GetNumberOfTracklets();
	fMultV0A =  esdV0->GetMTotV0A();
	fMultV0C =  esdV0->GetMTotV0C();
	//TOF hit  
	Int_t ntracksMatchedToTOF = 0; 
	Int_t ntracks = fESD->GetNumberOfTracks();
	for(Int_t itrk=0;itrk<ntracks;itrk++){
	  AliESDtrack* track = fESD->GetTrack(itrk);
	  if (!track) {
	    Printf("ERROR: Could not receive track %d", itrk);
	    continue;
	  }
	  //no track selection just TOF hit
	  if (track->IsOn(AliESDtrack::kTOFout)) ntracksMatchedToTOF++;
	}
	
	fTOFtracks = ntracksMatchedToTOF;
	if(fESDpid){ //get T0_TOF 
	  fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0);
	  fT0tofTrack =(Float_t) (fESDpid->GetTOFResponse().GetStartTime(10.0)); //Get start time "from all tracks 
	  //   printf("@@TOF T0 %f\n",fT0tofTrack);
	}	
	fpileup = fESD->IsPileupFromSPD();
	fTimeV0C = esdV0->GetV0CTime();
	fTimeV0A = esdV0->GetV0ATime();
	//T0 info
	mean = fESD->GetT0TOF();
	fVertex=fESD->GetT0zVertex();
	fMeanA = mean[1];
	fMeanC = mean[2];
	if (TMath::Abs(fMeanA)<9999 && TMath::Abs(fMeanC)<9999) 
	  fMeanAC = mean[0] ;
	fSumampC = tz->GetMultC();
	fSumampA = tz->GetMultA();
	amp=fESD->GetT0amplitude();
	amp_new = tz->GetT0NewAmplitude();
	time=fESD->GetT0time();
	for (Int_t i=0; i<24; i++){ 
	  if( time[i] !=0 ) {
	    ftime[i] = time[i];
	    if( amp[i]>0.1)  famp[i] = amp[i];
	    if(amp_new[i]>0) famp_new[i] = amp_new[i];	
	  }
	}
	//new raw OrA OrC TVDC all CFD
	for (Int_t ii=0; ii<5; ii++){ 
	  orA   = tz->GetOrA(ii);
	  orC   = tz->GetOrC(ii);
	  tvdc  = tz->GetTVDC(ii);
	  if(ii==0) {
	    fOrA[ii] = orA;
	    fOrC[ii] = orC;
	    fTVDC[ii] = tvdc;
	  } else {
	    if ( ii>0 && fOrA[ii-1] !=orA)	fOrA[ii] = orA;
	    else fOrA[ii]=-9999;
	    if ( ii>0 && fOrC[ii-1] !=orC)	fOrC[ii] = orC;
	    else fOrC[ii]=-9999;
	  if ( ii>0 && fTVDC[ii-1] != tvdc)	fTVDC[ii] = tvdc;
	  else fTVDC[ii]=-9999;
	  }
	}
	for (Int_t iii=0; iii<24; iii++)
	  fRawTime[iii] =  tz->GetTimeFull(iii,0);
	
	fT0pileup[0] = tz->GetPileupFlag();
	fT0pileup[1] = tz->GetBackgroundFlag();
	fT0pileup[2] = tz->GetSatellite();
	
	AliCentrality *centrality = fESD->GetCentrality();
	//The centrality function you can use
	if(centrality) {
	  fcentralityV0M = centrality->GetCentralityPercentile("V0M"); // returns the centrality 
	  fcentralityTRK = centrality->GetCentralityPercentile("TRK"); // returns the centrality 
	  fcentralityZDC = centrality->GetCentralityPercentile("ZEMvsZDC"); // returns the 
	}	
	// PFP bit
	fPFPbit = tz->GetT0PileupBits();
       
      } //selected
    } //ESD

  fT0OutTree->Fill();
    
  PostData(1, fOutputList);
  
}      

//________________________________________________________________________
void AliT0CalibAnalysisTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  /* fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
    */  

}
