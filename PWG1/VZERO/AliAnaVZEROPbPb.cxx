#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnaVZEROPbPb.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliCentrality.h"

// VZERO includes
#include "AliESDVZERO.h"
#include "AliESDVZEROfriend.h"

ClassImp(AliAnaVZEROPbPb)


AliAnaVZEROPbPb::AliAnaVZEROPbPb() 
  : AliAnalysisTaskSE("AliAnaVZEROPbPb"), fESD(0), fEsdV0(0), fOutputList(0), fClassesNames(0),
  fNFlags(0),
  fhAdcPMTNoTime(0),
  fhAdcPMTWithTime(0),
  fhTimePMT(0),
  fhWidthPMT(0),
  fhTimeCorr(0),
  fhPmtMult(0),
  fhV0ampl(0),
  fhEvents(0),
  fhVtxXYBB(0),
  fhVtxZBB(0),
  fhVtxXYBGA(0),
  fhVtxZBGA(0),
  fhVtxXYBGC(0),
  fhVtxZBGC(0),
  fhL2Triggers(0),
fhOnlineCharge(0),
fhRecoMult(0),
fhV0vsSPDCentrality(0),
fhTriggerBits(0),
fhTotRecoMult(0),
fhCentrality(0)
{
  // Constructor
  for(Int_t i = 0; i < 2; ++i) {
	fhAdcNoTime[i] = fhAdcWithTime[i] =  fhWidth[i] =  fhTime[i] = NULL;
	fhAdcTime[i] =  fhAdcWidth[i] = NULL;
  }

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnaVZEROPbPb::AliAnaVZEROPbPb(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fEsdV0(0), fOutputList(0),fNClasses(0),fClassesNames(0),
  fNFlags(0),
  fhAdcPMTNoTime(0),
  fhAdcPMTWithTime(0),
  fhTimePMT(0),
  fhWidthPMT(0),
  fhTimeCorr(0),
  fhPmtMult(0),
  fhV0ampl(0),
  fhEvents(0),
  fhVtxXYBB(0),
  fhVtxZBB(0),
  fhVtxXYBGA(0),
  fhVtxZBGA(0),
  fhVtxXYBGC(0),
  fhVtxZBGC(0),
  fhL2Triggers(0),
fhOnlineCharge(0),
fhRecoMult(0),
fhV0vsSPDCentrality(0),
fhTriggerBits(0),
fhTotRecoMult(0),
fhCentrality(0)
{
  // Constructor
  for(Int_t i = 0; i < 2; ++i) {
	fhAdcNoTime[i] = fhAdcWithTime[i] =  fhWidth[i] =  fhTime[i] = NULL;
	fhAdcTime[i] =  fhAdcWidth[i] = NULL;
  }
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
void AliAnaVZEROPbPb::SetClassesNames(const Char_t * nameList){
	TString  names("AllClasses,");
	names += nameList;
	fClassesNames = names.Tokenize(",");
	fNClasses = fClassesNames->GetEntriesFast();
}
//________________________________________________________________________
TH1F * AliAnaVZEROPbPb::CreateHisto1D(const char* name, const char* title,Int_t nBins, 
				    Double_t xMin, Double_t xMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  result->SetOption("hist");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  result->SetMarkerStyle(kFullCircle);
  return result;
}

//________________________________________________________________________
TH2F * AliAnaVZEROPbPb::CreateHisto2D(const char* name, const char* title,Int_t nBinsX, 
				    Double_t xMin, Double_t xMax,
				    Int_t nBinsY,
				    Double_t yMin, Double_t yMax,
				    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH2F* result = new TH2F(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  return result;
}


//________________________________________________________________________
void AliAnaVZEROPbPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
	if(fNClasses==0) {
		AliFatal("No Classes Defined");
		return;
	}

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);


    CreateQAHistos();
	CreateHistosPerL2Trigger();

  PostData(1, fOutputList);
 }
//________________________________________________________________________
void AliAnaVZEROPbPb::Init()
{

}

//________________________________________________________________________
void AliAnaVZEROPbPb::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  fEsdV0 = fESD->GetVZEROData();
  if (!fEsdV0) {
    Printf("ERROR: esd V0  not available");
    return;
  }

  FillQAHistos();

  FillPerL2TriggerHistos();

        

  PostData(1, fOutputList);
}      
//________________________________________________________________________
void AliAnaVZEROPbPb::CreateHistosPerL2Trigger(){
	
	fhOnlineCharge= new TH2F*[fNClasses];
	fhCentrality= new TH1F*[fNClasses];
	fhV0vsSPDCentrality= new TH2F*[fNClasses];
	fhRecoMult= new TH2F*[fNClasses];
	fhTotRecoMult= new TH1F*[fNClasses];
	fhTriggerBits= new TH1F*[fNClasses];

	fhL2Triggers = CreateHisto1D("hL2Triggers","L2 Triggers",fNClasses,0,fNClasses);
  	fOutputList->Add(fhL2Triggers);	  


	TIter iter(fClassesNames);
	TObjString* name;
	Int_t iClass =0;
	while((name = (TObjString*) iter.Next())){
		fhL2Triggers->GetXaxis()->SetBinLabel(iClass+1,name->String().Data());
		
		fhOnlineCharge[iClass] = CreateHisto2D(Form("hOnlineCharge_%s",name->String().Data()),Form("Online Charge for %s",name->String().Data()),100,0.,32000.,100,0.,32000,"V0A","V0C");
		fOutputList->Add(fhOnlineCharge[iClass]);	  

		fhCentrality[iClass] = CreateHisto1D(Form("hV0Centrality_%s",name->String().Data()),Form("V0 centrality for %s",name->String().Data()),100,0.,100.);
		fOutputList->Add(fhCentrality[iClass]);	  

		fhV0vsSPDCentrality[iClass] = CreateHisto2D(Form("hV0vsSPDCentrality_%s",name->String().Data()),Form("Centrality for %s",name->String().Data()),100,0.,100.,100,0.,100,"SPD Centrality (CL1)","V0 Centrality (V0M)");
		fOutputList->Add(fhV0vsSPDCentrality[iClass]);	  

		fhRecoMult[iClass] = CreateHisto2D(Form("hRecoMult_%s",name->String().Data()),Form("Reco Multiplicity for %s",name->String().Data()),100,0.,10000.,100,0.,20000,"V0A Mult","V0C Mult");
	  	fOutputList->Add(fhRecoMult[iClass]);	  

		fhTotRecoMult[iClass] = CreateHisto1D(Form("hTotRecoMult_%s",name->String().Data()),Form("Total Reco Multiplicity for %s",name->String().Data()),200,0.,20000.,"V0A + V0C Mult");
	  	fOutputList->Add(fhTotRecoMult[iClass]);	  

		fhTriggerBits[iClass] = CreateHisto1D(Form("hTriggerBits_%s",name->String().Data()),Form("Trigger Bits for %s",name->String().Data()),16,-0.5,15.5);
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(1,"BBA_AND_BBC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(2,"BBA_OR_BBC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(3,"BGA_AND_BBC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(4,"BGA");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(5,"BBA_AND_BGC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(6,"BGC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(7,"CTA1_AND_CTC1");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(8,"CTA1_OR_CTC1");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(9,"CTA2_AND_CTC2");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(10,"CTA1_OR_CTC1");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(11,"MTA_AND_MTC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(12,"MTA_OR_MTC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(13,"BBA");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(14,"BBC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(15,"BGA_OR_BGC");
		fhTriggerBits[iClass]->GetXaxis()->SetBinLabel(16,"All True BG");
	  	fOutputList->Add(fhTriggerBits[iClass]);	  

		iClass++;
	}
	
}
//________________________________________________________________________
void AliAnaVZEROPbPb::CreateQAHistos(){

  fNFlags    = CreateHisto2D("hNFlags","BB Flags",33,-0.5,32.5,33,-0.5,32.5,"V0A","V0C");
  fOutputList->Add(fNFlags);

	for(int iSide = 0; iSide < 2; ++iSide){
		TString side;
		if(iSide) side = "V0A";
		else side = "V0C";
		
		fhAdcNoTime[iSide] = CreateHisto1D(Form("hAdcNoTime%s",side.Data()),Form("ADC (no Leading Time) %s",side.Data()),200,0,200,"ADC charge","Entries");
	  	fOutputList->Add(fhAdcNoTime[iSide]);	  
		fhAdcWithTime[iSide] = CreateHisto1D(Form("hAdcWithTime%s",side.Data()),Form("ADC (with Leading Time) %s",side.Data()),200,0,200,"ADC charge","Entries");
	  	fOutputList->Add(fhAdcWithTime[iSide]);
	  	fhTime[iSide] = CreateHisto1D(Form("htimepmt%s",side.Data()),Form("Time measured by TDC %s",side.Data()),400,-100,100,"Leading time (ns)","Entries");
	  	fOutputList->Add(fhTime[iSide]);
	  	fhWidth[iSide] = CreateHisto1D(Form("hwidth%s",side.Data()),Form("Signal width measured by TDC %s",side.Data()),128,0,200,"Signal width (ns)","Entries");
	  	fOutputList->Add(fhWidth[iSide]);
	  	fhAdcWidth[iSide] = CreateHisto2D(Form("hadcwidth%s",side.Data()),Form("Time width vs ADC %s",side.Data()),200,0,200,128,0,200,"ADC charge","Width (ns)");
	  	fOutputList->Add(fhAdcWidth[iSide]);
	  	fhAdcTime[iSide] = CreateHisto2D(Form("hAdcTime%s",side.Data()),Form("ADC vs Time %s",side.Data()),1000,-100,100,200,0,200,"Time (ns)","ADC charge");
	  	fOutputList->Add(fhAdcTime[iSide]);
	}

  fhAdcPMTNoTime = CreateHisto2D("hadcpmtnotime","ADC vs PMT index (no leading time)",64,-0.5,63.5,200,0,200,"PMT index","ADC charge");
  fhAdcPMTWithTime = CreateHisto2D("hadcpmtwithtime","ADC vs PMT index (with leading time)",64,-0.5,63.5,200,0,200,"PMT index","ADC charge");

  fhTimePMT = CreateHisto2D("htimepmt","Time measured by TDC vs PMT index",64,-0.5,63.5,200,0,100,"PMT Index","Leading time (ns)");
  fhWidthPMT = CreateHisto2D("hwidthpmt","Time width vs PMT index",64,-0.5,63.5,128,0,200,"PMT Index","Signal width (ns)");

  fhTimeCorr = CreateHisto2D("htimecorr","Average time C side vs. A side",200,0,100,200,0,100,"Time V0A (ns)","Time V0C (ns");

  fhV0ampl  = CreateHisto1D("hV0ampl","V0 multiplicity in single channel (all V0 channels)",500,-0.5,499.5);

  fhEvents = CreateHisto2D("hTriggerDecision","Trigger Decision",3,-0.5,2.5,3,-0.5,2.5,"V0A Decision","V0C Decision");
  fhEvents->GetXaxis()->SetBinLabel(1,"Fake");
  fhEvents->GetXaxis()->SetBinLabel(2,"BB");
  fhEvents->GetXaxis()->SetBinLabel(3,"BG");
  fhEvents->GetYaxis()->SetBinLabel(1,"Fake");
  fhEvents->GetYaxis()->SetBinLabel(2,"BB");
  fhEvents->GetYaxis()->SetBinLabel(3,"BG");

  fhVtxXYBB = CreateHisto2D("fhVtxXYBB","XY SPD vertex (bb)",200,-2,2,200,-2,2);
  fhVtxZBB = CreateHisto1D("fhVtxZBB","Z SPD vertex (bb)",400,-50,50);
  fhVtxXYBGA = CreateHisto2D("fhVtxXYBGA","XY SPD vertex (bga)",200,-2,2,200,-2,2);
  fhVtxZBGA = CreateHisto1D("fhVtxZBGA","Z SPD vertex (bga)",400,-50,50);
  fhVtxXYBGC = CreateHisto2D("fhVtxXYBGC","XY SPD vertex (bgc)",200,-2,2,200,-2,2);
  fhVtxZBGC = CreateHisto1D("fhVtxZBGC","Z SPD vertex (bgc)",400,-50,50);

  fhPmtMult = CreateHisto2D("hV0CellMult","Number of fired PMTs (V0C vs V0A)",33,-0.5,32.5,33,-0.5,32.5,"# Cell (V0A)","# Cell (V0C)");
	fOutputList->Add(fhPmtMult);


  fOutputList->Add(fhAdcPMTNoTime);
  fOutputList->Add(fhAdcPMTWithTime);

  fOutputList->Add(fhTimePMT);
  fOutputList->Add(fhWidthPMT);

  fOutputList->Add(fhTimeCorr);
  fOutputList->Add(fhV0ampl);

  fOutputList->Add(fhEvents);

  fOutputList->Add(fhVtxXYBB);
  fOutputList->Add(fhVtxZBB);
  fOutputList->Add(fhVtxXYBGA);
  fOutputList->Add(fhVtxZBGA);
  fOutputList->Add(fhVtxXYBGC);
  fOutputList->Add(fhVtxZBGC);
  
}
//________________________________________________________________________
Bool_t AliAnaVZEROPbPb::IsGoodEvent(){
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if (!isSelected) return kFALSE;
  const AliESDVertex *primaryVtx = fESD->GetPrimaryVertex();
  if (!primaryVtx) return kFALSE;
  if (!primaryVtx->GetStatus()) return kFALSE;
  Double_t tPrimaryVtxPosition[3];
  primaryVtx->GetXYZ(tPrimaryVtxPosition);
  if (TMath::Abs(tPrimaryVtxPosition[2]) > 10.0) return kFALSE;

  AliCentrality *centrality = fESD->GetCentrality();
  if (centrality->GetQuality()) return kFALSE;
  
  return kTRUE;
}
//________________________________________________________________________
void AliAnaVZEROPbPb::FillPerL2TriggerHistos(){
   TString trigStr(fESD->GetFiredTriggerClasses());
  
//  Bool_t isGoodEvt = IsGoodEvent();
	TIter iter(fClassesNames);
	TObjString* name;
	Int_t iClass =0;
  while((name = (TObjString*) iter.Next())){
	
	if (!trigStr.Contains(name->String().Data()) && iClass>0) continue;
	
	fhOnlineCharge[iClass]->Fill(fEsdV0->GetTriggerChargeA(),fEsdV0->GetTriggerChargeC());
 	
	fhL2Triggers->Fill(iClass);
		
  	fhRecoMult[iClass]->Fill(fEsdV0->GetMTotV0A(),fEsdV0->GetMTotV0C());
  	fhTotRecoMult[iClass]->Fill(fEsdV0->GetMTotV0A()+fEsdV0->GetMTotV0C());

	for(int iTrig = 0; iTrig < 16; ++iTrig){
		if(fEsdV0->GetTriggerBits() & (1<<iTrig)) fhTriggerBits[iClass]->Fill(iTrig);	
	}
	
    AliCentrality *centrality = fESD->GetCentrality();
  	Float_t percentile = centrality->GetCentralityPercentile("V0M");
  	Float_t spdPercentile = centrality->GetCentralityPercentile("CL1");
  	if (spdPercentile < 0) spdPercentile = 0;
  	if (percentile < 0) percentile = 0;
	fhCentrality[iClass]->Fill(percentile);
	fhV0vsSPDCentrality[iClass]->Fill(spdPercentile,percentile);
	iClass++;
  }
}
//________________________________________________________________________
void AliAnaVZEROPbPb::FillQAHistos(){
  // Count V0 flags
  Int_t nV0A = 0, nV0C = 0;
  for(Int_t i = 0; i < 32; ++i) {
    if (fEsdV0->GetBBFlag(i)) nV0C++;
    if (fEsdV0->GetBBFlag(i+32)) nV0A++;
  }
  fNFlags->Fill((Float_t)nV0A,(Float_t)nV0C);

  for (Int_t i=0; i<64; ++i) {
	Int_t side = i/32;
	if (fEsdV0->GetTime(i) < 1e-6) {
		fhAdcNoTime[side]->Fill(fEsdV0->GetAdc(i));
		fhAdcPMTNoTime->Fill(i,fEsdV0->GetAdc(i));
    } else {
	  	fhAdcWithTime[side]->Fill(fEsdV0->GetAdc(i));
		fhAdcPMTWithTime->Fill(i,fEsdV0->GetAdc(i));
    }

	fhTime[side]->Fill(fEsdV0->GetTime(i));
	fhWidth[side]->Fill(fEsdV0->GetWidth(i));
	fhAdcWidth[side]->Fill(fEsdV0->GetAdc(i),fEsdV0->GetWidth(i));
	fhAdcTime[side]->Fill(fEsdV0->GetTime(i),fEsdV0->GetAdc(i));

    fhTimePMT->Fill(i,fEsdV0->GetTime(i));
    fhWidthPMT->Fill(i,fEsdV0->GetWidth(i));

  }

  fhTimeCorr->Fill(fEsdV0->GetV0ATime(),fEsdV0->GetV0CTime());

  AliESDVZERO::Decision flaga = fEsdV0->GetV0ADecision();
  AliESDVZERO::Decision flagc = fEsdV0->GetV0CDecision();

  fhEvents->Fill(flaga,flagc);

  const AliESDVertex *vtx = fESD->GetPrimaryVertexSPD();

  if (flaga <= 1 && flagc <=1) {
    fhVtxXYBB->Fill(vtx->GetXv(),vtx->GetYv());
    fhVtxZBB->Fill(vtx->GetZv());
  }
  else {
    if (flaga == 2) {
      fhVtxXYBGA->Fill(vtx->GetXv(),vtx->GetYv());
      fhVtxZBGA->Fill(vtx->GetZv());
    }
    if (flagc == 2) {
      fhVtxXYBGC->Fill(vtx->GetXv(),vtx->GetYv());
      fhVtxZBGC->Fill(vtx->GetZv());
    }
  }

  fhPmtMult->Fill(fEsdV0->GetNbPMV0A(),fEsdV0->GetNbPMV0C());
  for(Int_t i = 0; i < 64; i++) {
    fhV0ampl->Fill(fEsdV0->GetMultiplicity(i));
  }
  
}
//________________________________________________________________________
void AliAnaVZEROPbPb::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
