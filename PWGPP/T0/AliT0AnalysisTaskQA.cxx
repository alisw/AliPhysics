#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliT0AnalysisTaskQA.h"
#include "AliESDpid.h"

//#include "AliCDBMetaData.h"
//#include "AliCDBId.h"
//#include "AliCDBEntry.h"
//#include "AliCDBManager.h"
//#include "AliCDBStorage.h"

// Task should calculate channels offset 
// Authors: Alla
//last change 23 Feb 2012 FK

ClassImp(AliT0AnalysisTaskQA)
//________________________________________________________________________
AliT0AnalysisTaskQA::AliT0AnalysisTaskQA() 
  : AliAnalysisTaskSE(),  fESD(0x0), fTzeroObject(0x0),
  fTzeroORA(0x0), fTzeroORC(0x0), fResolution(0x0), fTzeroORAplusORC(0x0), fTzeroTof(0x0),
    fRunNumber(0),fTimeVSAmplitude(0x0),fCFDVSPmtId(0x0),fSPDVertexVST0Vertex(0x0),
    fOrAvsNtracks(0x0), fOrCvsNtracks(0x0), fT0vsNtracks(0x0),fT0TimevsT0Tof(0x0),
  fESDpid(new AliESDpid()), f0TVX(0), fTriggerCounter(0x0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain

// ########### NEVER define slots in the IO constructor
//  DefineInput(0,  TChain::Class());
//  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliT0AnalysisTaskQA::AliT0AnalysisTaskQA(const char *name) 
  : AliAnalysisTaskSE(name),  fESD(0x0), fTzeroObject(0x0),
  fTzeroORA(0x0), fTzeroORC(0x0), fResolution(0x0), fTzeroORAplusORC(0x0), fTzeroTof(0x0),
    fRunNumber(0),fTimeVSAmplitude(0x0),fCFDVSPmtId(0x0),fSPDVertexVST0Vertex(0x0),
    fOrAvsNtracks(0x0), fOrCvsNtracks(0x0), fT0vsNtracks(0x0),fT0TimevsT0Tof(0x0),
    fESDpid(new AliESDpid()), f0TVX(0), fTriggerCounter(0x0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TObjArray::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
}

//________________________________________________________________________
AliT0AnalysisTaskQA::~AliT0AnalysisTaskQA() 
{
  // Destructor
  // printf("AliT0CalibOffsetChannels~AliT0CalibOffsetChannels() ");
  delete fTzeroORA;
  delete fTzeroORC;
  delete fResolution;
  delete fTzeroORAplusORC;
  delete fTzeroTof;
  delete [] fTimeVSAmplitude;
  delete fCFDVSPmtId;
  delete fSPDVertexVST0Vertex;
  delete fOrAvsNtracks;
  delete fOrCvsNtracks;
  delete fT0vsNtracks;
  delete fT0TimevsT0Tof;
  delete fTriggerCounter;

  delete fESDpid;
  delete fTzeroObject;
}

//------------------------------------------------------------------
void AliT0AnalysisTaskQA::UserCreateOutputObjects()
{
  // Create histograms
 fTimeVSAmplitude = new TH2F*[kNPMT0];

 for (Int_t i=0; i<kNPMT0; i++) {
    fTimeVSAmplitude[i]= new TH2F (Form("fTimeVSAmplitude%d",i+1),"fTimeVsAmplitude",600, -10, 50,500,9000,11000);
  }

  fTzeroORAplusORC = new TH1F("fTzeroORAplusORC","ORA+ORC /2",100,-2000,2000);   //or A plus or C 
  fTzeroTof        = new TH1F("fTzeroTof","t0 from TOF",100,-2000,2000);   //t0 start time from TOF
  fResolution      = new TH1F("fResolution","fResolution",100,-500,500);// or A minus or C spectrum
  fTzeroORA        = new TH1F("fTzeroORA","fTzeroORA",100,-2000,2000);// or A spectrum
  fTzeroORC        = new TH1F("fTzeroORC","fTzeroORC",100,-2000,2000);// or C spectrum
  fCFDVSPmtId      = new TH2F("fCFDVSPmtId","fCFDVSPmtId",24,0,24,1000,2000,11000);  // 
  fSPDVertexVST0Vertex = new TH2F("fSPDVertexVST0Vertex","fSPDVertexVST0Vertex",30,-30,30,30,-30,30);
  fOrAvsNtracks = new TH2F("fAvstracks", "A vs tracks",100, 0, 100, 200, -1000, 1000);
  fOrCvsNtracks = new TH2F("fCvstracks", "C vs tracks",100, 0, 100, 200, -1000, 1000);
  fT0vsNtracks  = new TH2F("fT0ACvstrackes", "T0AC vs tracks",100, 0, 100, 200, -1000, 1000); 
  fT0TimevsT0Tof = new TH2F("fT0TimevsT0Tof", "fT0TimevsT0Tof",50, -1000,1000, 50, -1000,1000); 
  f0TVX        = new TH1F("f0TVX","0TVX position [channels]",200,-500,500);// or C spectrum
  fTriggerCounter = new TH1I("fTriggerCounter", "Trigger counter for calculating efficiencies; Trigger; Nevents", 7,0,7);// counter of events satisfying different triggers
  fTriggerCounter->GetXaxis()->SetBinLabel(1,"All");
  fTriggerCounter->GetXaxis()->SetBinLabel(2,"NoPileup");
  fTriggerCounter->GetXaxis()->SetBinLabel(3,"SPD_Vertex");
  fTriggerCounter->GetXaxis()->SetBinLabel(4,"V0_Time");
  fTriggerCounter->GetXaxis()->SetBinLabel(5,"CINT7-B"); 
  fTriggerCounter->GetXaxis()->SetBinLabel(6,"C0TVX-B");
  fTriggerCounter->GetXaxis()->SetBinLabel(7,"CADAND-B");
  
  fTzeroObject     = new TObjArray(0);
  fTzeroObject->SetOwner(kTRUE);
  
  for (Int_t i=0; i<kNPMT0; i++)
    fTzeroObject->AddAtAndExpand(fTimeVSAmplitude[i],i);

  fTzeroObject->AddAtAndExpand(fCFDVSPmtId,24);
  fTzeroObject->AddAtAndExpand(fSPDVertexVST0Vertex,25);
  fTzeroObject->AddAtAndExpand(fTzeroORAplusORC, 26);
  fTzeroObject->AddAtAndExpand(fResolution, 27);
  fTzeroObject->AddAtAndExpand(fTzeroORA, 28);
  fTzeroObject->AddAtAndExpand(fTzeroORC, 29);
  fTzeroObject->AddAtAndExpand(fT0vsNtracks, 30);
  fTzeroObject->AddAtAndExpand(fOrAvsNtracks,31);
  fTzeroObject->AddAtAndExpand(fOrCvsNtracks, 32);
  fTzeroObject->AddAtAndExpand(fTzeroTof, 33);
  fTzeroObject->AddAtAndExpand(fT0TimevsT0Tof, 34);
  fTzeroObject->AddAtAndExpand(f0TVX, 35);
  fTzeroObject->AddAtAndExpand(fTriggerCounter, 36);

  PostData(1, fTzeroObject);
  // Called once
}



//________________________________________________________________________
void AliT0AnalysisTaskQA::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  fRunNumber =  fESD->GetRunNumber() ; 

  const Double32_t* time = fESD->GetT0time();
  const Double32_t* amplitude = fESD->GetT0amplitude();

 
  for (Int_t i=0; i<kNPMT0; i++) {
    if(time[i]<99999  )
      {
	if(amplitude[i]<50){
	  fTimeVSAmplitude[i]->Fill(amplitude[i],time[i]);
	}
	else {
	  fTimeVSAmplitude[i]->Fill(amplitude[i]/1000.,time[i]); //in RUN2 we don't convert to MIPs
	}	
	fCFDVSPmtId->Fill(i,time[i]);
      }
  }

  const Double32_t* mean = fESD->GetT0TOF();
  Double32_t orA = mean[1];
  Double32_t orC = mean[2];
  Int_t ntracks = fESD->GetNumberOfTracks();

  Int_t ntracksMatchedToTOF = 0; 
  for(Int_t itrk=0;itrk<ntracks;itrk++){
    AliESDtrack* track = fESD->GetTrack(itrk);
    if (!track) {
      Printf("ERROR: Could not receive track %d", itrk);
      continue;
    }
    //no track selection just TOF hit
    if (track->IsOn(AliESDtrack::kTOFout)) ntracksMatchedToTOF++;
  }

  if(orA<9999){
    fTzeroORA->Fill(orA);
    fOrAvsNtracks->Fill(ntracksMatchedToTOF, orA);
  }
  if(orC<9999) {
    fTzeroORC->Fill(orC);
    fOrCvsNtracks->Fill(ntracksMatchedToTOF, orC);
  }
  if(orA<9999 && orC<9999) {
    fResolution->Fill((orA-orC)/2.);
    fTzeroORAplusORC->Fill(mean[0]);
    fT0vsNtracks->Fill(ntracksMatchedToTOF, mean[0]);
  }

  if(fESDpid){ //get T0_TOF 
    fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0);
    Float_t t0tofTrack =(Float_t) (fESDpid->GetTOFResponse().GetStartTime(10.0)); //Get start time from all tracks 
    if (t0tofTrack !=0) fTzeroTof->Fill(t0tofTrack);

    if(orA<9999 && orC<9999 && t0tofTrack !=0){ // T0 time  and  TOF time simultaneously
      fT0TimevsT0Tof->Fill(t0tofTrack, mean[0]);
    }
  }
  
  Double32_t t0vertex = fESD->GetT0zVertex();
  Double32_t esdzvertex;
  const AliESDVertex * esdvertex = fESD->GetPrimaryVertex();
  Int_t nofcontrib=-1;
  if(esdvertex && t0vertex<999)
    {
      nofcontrib=esdvertex->GetNContributors();
      if(nofcontrib>1)
	{
	  esdzvertex=esdvertex->GetZ();
	  fSPDVertexVST0Vertex->Fill(t0vertex,esdzvertex);
	}
    }
  //0TVX position
  AliESDTZERO* tz= (AliESDTZERO*) fESD->GetESDTZERO();
  Float_t tvdc = tz->GetTVDC(0);
  if (tvdc!=0) f0TVX->Fill(tvdc);
  
  Bool_t fIs0TVX_B,fIsCINT7_B,fIsV0Cut,fIsCADAND_B;
  ULong64_t fMask0TVX_B,fMaskCINT7_B,fMaskCADAND_B;
  //Trigger efficiencies
  Bool_t isPileup = fESD->IsPileupFromSPD();
  Float_t vertexSPD=-9999.;
  Int_t ncontSPD = fESD->GetPrimaryVertexSPD()->GetNContributors();
  if(ncontSPD>=0 ) {
      vertexSPD = fESD->GetPrimaryVertex()->GetZ();
  }
  Float_t timeV0A = fESD->GetVZEROData()->GetV0ATime();
  Float_t timeV0C = fESD->GetVZEROData()->GetV0CTime();
  TString triggerName = fESD->GetFiredTriggerClasses();			//! trigger name
  
  fTriggerCounter->Fill("All",1);
  if(!isPileup)fTriggerCounter->Fill("NoPileup",1);
  {
      if(vertexSPD>-10. && vertexSPD<10.)
      {
          printf("VertexSPD: %f\n",vertexSPD);
          fTriggerCounter->Fill("SPD_Vertex",1);
          if((timeV0A+timeV0C)>11.5 && (timeV0A+timeV0C)<17.5 && (timeV0A-timeV0C)>5.5 && (timeV0A-timeV0C)<11.5)
          {
	    fTriggerCounter->Fill("V0_Time",1);
	    if(triggerName.Contains("CINT7-B-NOPF-CENT"))
		fTriggerCounter->Fill("CINT7-B",1);
	    if(triggerName.Contains("CADAND-B-NOPF-CENT"))
	      fTriggerCounter->Fill("CADAND-B",1);
	    if(triggerName.Contains("C0TVX-B-NOPF-CENT"))
	      fTriggerCounter->Fill("C0TVX-B",1);
	    
          }
      }
  }
  PostData(1, fTzeroObject);
  
}      
 //________________________________________________________________________
void AliT0AnalysisTaskQA::Terminate(Option_t *) 
{
  
   // Called once at the end of the query
}
