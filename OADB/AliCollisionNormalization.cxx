//-------------------------------------------------------------------------
//                      Implementation of   Class AliCollisionNormalization
//
//  This class is used to store the vertex ditributions in the data
//  and in Monte Carlo, needed to compute the real number of
//  collisions a given sample is corresponding to.
//  The strategy matches what described in CERN-THESIS-2009-033 p 119.
//
//    Author:     Michele Floris, CERN
//-------------------------------------------------------------------------

#include "AliCollisionNormalization.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "TFile.h"
#include "TCanvas.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"

ClassImp(AliCollisionNormalization)

const char * AliCollisionNormalization::fgkProcLabel[] = {"SD","DD","ND", "Unknown"};

AliCollisionNormalization::AliCollisionNormalization() :
  TObject(),
  fNbinsVz(0),
  fMinVz(0),
  fMaxVz(0),
  fZRange(9.99),
  fIsMC(0),
  fReferenceXS(0),
  fVerbose(0),
  fEnergy(900),
  fHistVzData       (0),
  fHistProcTypes    (0),
  fHistStatBin0     (0),
  fHistStat     (0),
  fInputEvents(0),   
  fPhysSelEvents(0),
  fBgEvents(0),      
  fAllEvents(0),     
  fAllEventsZRange(0),
  fAllEventsZRangeMult1(0), 
  fAllEventsInBin0ZRange(0),
  fTrigEffBin0(0)  
{

  // ctor

  fNbinsVz =  30;
  fMinVz   = -15;
  fMaxVz   = +15;

  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
      fHistVzMCGen[iproc] = 0;
      fHistVzMCRec[iproc] = 0;
      fHistVzMCTrg[iproc] = 0;
  }
  
  
  BookAllHistos();
}

AliCollisionNormalization::AliCollisionNormalization(Int_t nbinz, Float_t minz, Float_t maxz):
  TObject(),
  fNbinsVz(0),
  fMinVz(0),
  fMaxVz(0),
  fZRange(9.99),
  fIsMC(0),
  fReferenceXS(0),
  fVerbose(0),
  fEnergy(900),
  fHistVzData       (0),
  fHistProcTypes    (0),
  fHistStatBin0     (0),
  fHistStat     (0),
  fInputEvents(0),   
  fPhysSelEvents(0),
  fBgEvents(0),      
  fAllEvents(0),     
  fAllEventsZRange(0),
  fAllEventsZRangeMult1(0), 
  fAllEventsInBin0ZRange(0),
  fTrigEffBin0(0)  
{

  // ctor, allows setting binning
  fNbinsVz = nbinz;
  fMinVz   = minz;
  fMaxVz   = maxz;

  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
      fHistVzMCGen[iproc] = 0;
      fHistVzMCRec[iproc] = 0;
      fHistVzMCTrg[iproc] = 0;
  }

  BookAllHistos();
}

AliCollisionNormalization::AliCollisionNormalization(const char * dataFile, const char * dataListName, 
						     const char * mcFile,   const char * mcListName,
						     const char * eventStatFile) :
  TObject(),
  fNbinsVz(0),
  fMinVz(0),
  fMaxVz(0),
  fZRange(9.99),
  fIsMC(0),
  fReferenceXS(0),
  fVerbose(0),
  fEnergy(900),
  fHistVzData       (0),
  fHistProcTypes    (0),
  fHistStatBin0     (0),
  fHistStat     (0),
  fInputEvents(0),   
  fPhysSelEvents(0),
  fBgEvents(0),      
  fAllEvents(0),     
  fAllEventsZRange(0),
  fAllEventsZRangeMult1(0), 
  fAllEventsInBin0ZRange(0),
  fTrigEffBin0(0)  

{

  // ctor, loads histograms from file
  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
      fHistVzMCGen[iproc] = 0;
      fHistVzMCRec[iproc] = 0;
      fHistVzMCTrg[iproc] = 0;
  }


  TFile * fdata = new TFile (dataFile);
  TFile * fmc   = new TFile (mcFile  );
  TFile * fstat = new TFile(eventStatFile);
  
  if (!fdata->IsOpen() || !fmc->IsOpen() || !fstat->IsOpen()) {
    AliFatal("Cannot open input file(s)");
  }
  
  TList * ldata = (TList*) fdata->Get(dataListName);
  TList * lmc   = (TList*) fmc  ->Get(mcListName  );

  AliCollisionNormalization * cndata = (AliCollisionNormalization*) ldata->FindObject("AliCollisionNormalization");
  AliCollisionNormalization * cnmc   = (AliCollisionNormalization*) lmc  ->FindObject("AliCollisionNormalization");


  // Assign or book all histos
  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
    fHistVzMCGen[iproc]= cnmc->GetVzMCGen(iproc);
    fHistVzMCRec[iproc]= cnmc->GetVzMCRec(iproc);
    fHistVzMCTrg[iproc]= cnmc->GetVzMCTrg(iproc);
  }
  fHistVzData        = cndata->GetVzData();
  fHistProcTypes     = cnmc->GetHistProcTypes();
    
  fHistStatBin0      =  (TH1F*) fstat->Get("fHistStatistics_Bin0");
  fHistStat          =  (TH1F*) fstat->Get("fHistStatistics");
  
}


AliCollisionNormalization::~AliCollisionNormalization(){

  // dtor
  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
    if(fHistVzMCGen[iproc]) { delete fHistVzMCGen[iproc]      ; fHistVzMCGen[iproc]      =0;}
    if(fHistVzMCRec[iproc]) { delete fHistVzMCRec[iproc]      ; fHistVzMCRec[iproc]      =0;}
    if(fHistVzMCTrg[iproc]) { delete fHistVzMCTrg[iproc]      ; fHistVzMCTrg[iproc]      =0;}    
  }
  
  if(fHistVzData       ) { delete fHistVzData       ; fHistVzData       =0;}
  if(fHistStatBin0     ) { delete fHistStatBin0     ; fHistStatBin0     =0;}
  if(fHistStat         ) { delete fHistStat         ; fHistStat         =0;}
  if(fHistProcTypes    ) { delete fHistProcTypes    ; fHistProcTypes    =0;}

}
  
void AliCollisionNormalization::BookAllHistos(){
  // books all histos
  // Book histos of vz distributions vs multiplicity
  // if vzOnly == kTRUE, it returns a 1d histo with vz dist only

  // Do not attach histos to the current directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
    fHistVzMCGen [iproc]      = (TH2F*) BookVzHisto(TString("fHistVzMCGen")+ fgkProcLabel[iproc]      ,"Vz distribution of generated events vs rec multiplicity    ");
    fHistVzMCRec [iproc]      = (TH2F*) BookVzHisto(TString("fHistVzMCRec")+ fgkProcLabel[iproc]      ,"Vz distribution of reconstructed events vs rec multiplicity");
    fHistVzMCTrg [iproc]      = (TH2F*) BookVzHisto(TString("fHistVzMCTrg")+ fgkProcLabel[iproc]      ,"Vz distribution of triggered events vs rec multiplicity    ");    
  }
  fHistVzData        = (TH2F*) BookVzHisto("fHistVzData"       ,"Vz distribution of triggered events vs rec multiplicity    ");
  fHistProcTypes     = new TH1F   ("fHistProcTypes", "Number of events in the different process classes", kNProcs, -0.5 , kNProcs-0.5);

  fHistProcTypes->GetXaxis()->SetBinLabel(kProcSD+1,"SD");
  fHistProcTypes->GetXaxis()->SetBinLabel(kProcND+1,"ND");
  fHistProcTypes->GetXaxis()->SetBinLabel(kProcDD+1,"DD");
  fHistProcTypes->GetXaxis()->SetBinLabel(kProcUnknown+1,"Unknown");

  TH1::AddDirectory(oldStatus);

}

TH1 * AliCollisionNormalization::BookVzHisto(const char * name , const char * title, Bool_t vzOnly) {

  // Book histos of vz distributions vs multiplicity
  // if vzOnly == kTRUE, it returns a 1d histo with vz dist only

 
  // Do not attach histos to the current directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TH1 * h;
  Double_t binLimitsVtx[] = {-30,-25,-20,-15,-10,-7,-5.5,-4,-3,-2,-1,0,1,2,3,4,5.5,7,10,15,20,25,30};
  if (vzOnly) {
    h = new TH1F(name,title,22,binLimitsVtx);
  } else {
    h = new TH2F(name,title,22,binLimitsVtx,100,Double_t(-0.5),Double_t(99.5));
  }
 
  h->SetXTitle("V_{z} (cm)");
  h->SetYTitle("n_{trk}");
  h->Sumw2();

  TH1::AddDirectory(oldStatus);

  return h;

}

TH2F * AliCollisionNormalization::GetVzMCGen (Int_t procType) { 

  //  returns MC gen histo. If proc type < 0 sums over all proc types, reweighting the XS

  if(procType>=0)  return fHistVzMCGen[procType] ;

  TH2F * sum = (TH2F*) fHistVzMCGen[kProcSD]->Clone();
  sum->Reset();

  for(Int_t iproc = 0; iproc < kProcUnknown; iproc++){
    sum->Add(fHistVzMCGen[iproc],GetProcessWeight(iproc));    
  }
  
  return sum;
}
TH2F *   AliCollisionNormalization::GetVzMCRec       (Int_t procType) { 

  //  returns MC rec histo. If proc type < 0 sums over all proc types, reweighting the XS

  if(procType>=0)  return fHistVzMCRec[procType] ;

  TH2F * sum = (TH2F*) fHistVzMCRec[kProcSD]->Clone();
  sum->Reset();

  for(Int_t iproc = 0; iproc < kProcUnknown; iproc++){
    sum->Add(fHistVzMCRec[iproc],GetProcessWeight(iproc));    
  }
  
  return sum;

}


TH2F *   AliCollisionNormalization::GetVzMCTrg       (Int_t procType) { 

  //  returns MC trg histo. If proc type < 0 sums over all proc types, reweighting the XS

  if(procType>=0)  return fHistVzMCTrg[procType] ;

  TH2F * sum = (TH2F*) fHistVzMCTrg[kProcSD]->Clone();
  sum->Reset();

  for(Int_t iproc = 0; iproc < kProcUnknown; iproc++){
    sum->Add(fHistVzMCTrg[iproc],GetProcessWeight(iproc));    
  }
  
  return sum;

}

Double_t AliCollisionNormalization::ComputeNint() {

  // Compute the number of collisions applying all corrections
  // TODO: check error propagation  

  TH1 * hVzData = fHistVzData->ProjectionX("_px",2,-1,"E"); // Skip zero bin 
  Int_t allEventsWithVertex = (Int_t) fHistVzData->Integral(0, fHistVzData->GetNbinsX()+1,
							    2, fHistVzData->GetNbinsY()+1); // include under/overflow!, skip 0 bin

  // Assign histos reweighted with measured XS
  TH2F * histVzMCGen = GetVzMCGen(-1);
  TH2F * histVzMCTrg = GetVzMCTrg(-1);
  TH2F * histVzMCRec = GetVzMCRec(-1);

  // Before we start: print number of input events to the physics
  // selection: this allows to crosscheck that all runs were
  // successfully processed (useful if running on grid: you may have a
  // crash without noticing it).

  fInputEvents = fHistStat->GetBinContent(1,1);
  fPhysSelEvents = fHistStat->GetBinContent( fHistStat->GetNbinsX(),1);

  AliInfo(Form("Input Events (No cuts: %d, After Phys. Sel.:%d)",
	       Int_t(fInputEvents),
	       Int_t(fPhysSelEvents))); // Fixme: will this change with a different trigger scheme?

  // Get or compute BG. This assumes the CINT1B suite
  Double_t triggeredEventsWith0MultWithBG = fHistVzData->Integral(0, fHistVzData->GetNbinsX()+1,1, 1);
  //  Double_t bg = 0; // This will include beam gas + accidentals
  if (fHistStatBin0->GetNbinsY() > 4) { // FIXME: we need a better criterion to decide...
    AliInfo("Using BG computed by Physics Selection");
    // WARNING
    // CHECK IF THE COMPUTATION OF BG OFFSET IS STILL OK, IN CASE OF CHANGES TO THE PHYSICS SELECTION
    Int_t bgOffset = 0; 
    Int_t nbiny = fHistStatBin0->GetNbinsY();
    for(Int_t ibiny =1; ibiny <= nbiny; ibiny++){
      bgOffset++;
      printf("%d, %s\n", bgOffset, fHistStatBin0->GetYaxis()->GetBinLabel(ibiny) );
      
      if(!strncmp("All B", fHistStatBin0->GetYaxis()->GetBinLabel(ibiny),5)) break;
      if((ibiny+1)==nbiny) AliFatal("Cannot compute bg offset");
    }
    
    if(fVerbose > 2) AliInfo(Form("BG Offset: %d",bgOffset));
    

    fBgEvents = fHistStatBin0->GetBinContent(fHistStatBin0->GetNbinsX(), bgOffset+AliPhysicsSelection::kStatRowBG);
    fBgEvents += fHistStatBin0->GetBinContent(fHistStatBin0->GetNbinsX(),bgOffset+AliPhysicsSelection::kStatRowAcc);
    Int_t cint1B = (Int_t) fHistStatBin0->GetBinContent(fHistStatBin0->GetNbinsX(),1);	
    if (cint1B != Int_t(triggeredEventsWith0MultWithBG)) {
      AliWarning(Form("Events in bin0 from physics selection and local counter not consistent: %d - %d", cint1B, Int_t(triggeredEventsWith0MultWithBG)));
    }
  } else {
    AliInfo("Computing BG using CINT1A/B/C/E, ignoring intensities");
    AliError("This will only work for early runs!!! USE BG FROM PHYSICS SELECTION INSTEAD");
    Int_t icol = fHistStatBin0->GetNbinsX();
    Int_t cint1B = (Int_t) fHistStatBin0->GetBinContent(icol,1);	
    Int_t cint1A = (Int_t) fHistStatBin0->GetBinContent(icol,2);	
    Int_t cint1C = (Int_t) fHistStatBin0->GetBinContent(icol,3);	
    Int_t cint1E = (Int_t) fHistStatBin0->GetBinContent(icol,4);      
    fBgEvents   = cint1A + cint1C-2*cint1E ; // FIXME: to be changed to take into account ratios of events
    if (cint1B != triggeredEventsWith0MultWithBG) {
      AliWarning(Form("Events in bin0 from physics selection and local counter not consistent: %d - %d", cint1B, (Int_t)triggeredEventsWith0MultWithBG));
    }
  }

  Double_t triggeredEventsWith0Mult = triggeredEventsWith0MultWithBG - fBgEvents; 
  if(fVerbose > 0) AliInfo(Form("Measured events with vertex: %d",allEventsWithVertex));
  Double_t bin0 = fHistVzData->Integral(0, fHistVzData->GetNbinsX()+1,
					1, 1);
  if(fVerbose > 0) AliInfo(Form("Zero Bin, Meas: %2.2f, BG: %2.2f, Meas - BG: %2.2f", bin0, fBgEvents, triggeredEventsWith0Mult)); 

  // This pointers are here so that I use the same names used in Jan Fiete's code (allows for faster/easier comparison)

  TH2* eTrig =    histVzMCTrg; 
  TH2* eTrigVtx = histVzMCRec; 
  TH1* eTrigVtxProjx = eTrigVtx->ProjectionX("eTrigVtxProjx", 2, eTrigVtx->GetNbinsX()+1);

  // compute trigger and vertex efficiency
  TH2 * effVtxTrig = (TH2*) histVzMCRec->Clone("effVtxTrig");
  effVtxTrig->Reset();
  effVtxTrig->Divide(histVzMCRec,histVzMCGen,1,1,"B");
  // Apply correction to data to get corrected events
  TH2 * correctedEvents = (TH2*) fHistVzData->Clone("correctedEvents");
  correctedEvents->Divide(effVtxTrig);

  //  TH1 * correctedEvents = fHistVzData->ProjectionX("eTrigVtxProjx", 2, eTrigVtx->GetNbinsX()+1); 

  // loop over vertex bins
  for (Int_t i = 1; i <= fHistVzData->GetNbinsX(); i++)
    {
      Double_t alpha = (Double_t) hVzData->GetBinContent(i) / allEventsWithVertex;
      Double_t events = alpha * triggeredEventsWith0Mult;
      
      // if (histVzMCRec->GetBinContent(i,1) == 0)
      //   continue;
      if (!eTrigVtxProjx->GetBinContent(i) || !eTrig->Integral(0, eTrig->GetNbinsX()+1, 1, 1))
        continue;

      Double_t fZ = eTrigVtxProjx->Integral(0, eTrigVtxProjx->GetNbinsX()+1) / eTrigVtxProjx->GetBinContent(i) *
        eTrig->GetBinContent(i, 1) / eTrig->Integral(0, eTrig->GetNbinsX()+1, 1, 1);

      events *= fZ;

      // multiply with trigger correction
      Double_t correction = histVzMCGen->GetBinContent(i,1)/histVzMCTrg->GetBinContent(i,1);
      events *= correction;

      if (fVerbose > 1) Printf("  Bin %d, alpha is %.2f%%, fZ is %.3f, number of events with 0 mult.: %.2f (MC comparison: %.2f)", i, alpha * 100., fZ, events, 
			       histVzMCTrg->GetBinContent(i,1));
      correctedEvents->SetBinContent(i, 1, events);
    }

  
  
  // Integrate correctedEvents over full z range
  fAllEvents       = correctedEvents->Integral(0, correctedEvents->GetNbinsX()+1,0, correctedEvents->GetNbinsY()+1);
  // Integrate correctedEvents over needed z range
  fAllEventsZRange = correctedEvents->Integral(correctedEvents->GetXaxis()->FindBin(-fZRange), correctedEvents->GetXaxis()->FindBin(fZRange-0.0001),
						       0, correctedEvents->GetNbinsY()+1);

  fAllEventsZRangeMult1 = correctedEvents->Integral(correctedEvents->GetXaxis()->FindBin(-fZRange), correctedEvents->GetXaxis()->FindBin(fZRange-0.0001), 
						    2,correctedEvents->GetNbinsY()+1);

  fAllEventsInBin0ZRange = correctedEvents->Integral(correctedEvents->GetXaxis()->FindBin(-fZRange), correctedEvents->GetXaxis()->FindBin(fZRange-0.0001),1,1);

  if(fVerbose > 1) AliInfo(Form("Results in |Vz| < %3.3f",fZRange));
  if(fVerbose > 1) AliInfo(Form(" Events in Bin0: %2.2f, With > 1 track: %2.2f, All corrected: %2.2f",
				fAllEventsInBin0ZRange,
				fAllEventsZRangeMult1,
				fAllEventsZRange));
  Int_t nbin = histVzMCTrg->GetNbinsX();
  fTrigEffBin0 =  histVzMCTrg->Integral(0,nbin+1,1,1)/histVzMCGen->Integral(0,nbin+1,1,1);
  if(fVerbose > 1) {    
    AliInfo(Form("Trigger Efficiency in the zero bin: %3.3f", fTrigEffBin0 ));
  }
  

  AliInfo(Form("Number of collisions in full phase space: %2.2f", fAllEvents));
//   effVtxTrig->Draw();
//   new TCanvas();
//   correctedEvents->Draw(); // FIXME: debug

  return fAllEvents;
}

Long64_t AliCollisionNormalization::Merge(TCollection* list)
{
  // Merge a list of AliCollisionNormalization objects with this
  // (needed for PROOF).  
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  // collections of all histograms
  const Int_t nHists = kNProcs*3+5;
  TList collections[nHists];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliCollisionNormalization* entry = dynamic_cast<AliCollisionNormalization*> (obj);
    if (entry == 0) 
      continue;
    Int_t ihist = -1;

    for(Int_t iproc = 0; iproc < kNProcs; iproc++){
      if (entry->fHistVzMCGen[iproc]      ) collections[++ihist].Add(entry->fHistVzMCGen[iproc]      );
      if (entry->fHistVzMCRec[iproc]      ) collections[++ihist].Add(entry->fHistVzMCRec[iproc]      );
      if (entry->fHistVzMCTrg[iproc]      ) collections[++ihist].Add(entry->fHistVzMCTrg[iproc]      );
    }
    if (entry->fHistVzData       ) collections[++ihist].Add(entry->fHistVzData       );
    if (entry->fHistProcTypes    ) collections[++ihist].Add(entry->fHistProcTypes    );
    if (entry->fHistStatBin0     ) collections[++ihist].Add(entry->fHistStatBin0     );
    if (entry->fHistStat         ) collections[++ihist].Add(entry->fHistStat         );

    count++;
  }

  Int_t ihist = -1;
  for(Int_t iproc = 0; iproc < kNProcs; iproc++){
    if (fHistVzMCGen[iproc]      ) fHistVzMCGen[iproc]      ->Merge(&collections[++ihist]);
    if (fHistVzMCRec[iproc]      ) fHistVzMCRec[iproc]      ->Merge(&collections[++ihist]);
    if (fHistVzMCTrg[iproc]      ) fHistVzMCTrg[iproc]      ->Merge(&collections[++ihist]);

  }  
  if (fHistVzData       ) fHistVzData       ->Merge(&collections[++ihist]);
  if (fHistProcTypes    ) fHistProcTypes    ->Merge(&collections[++ihist]);
  if (fHistStatBin0     ) fHistStatBin0     ->Merge(&collections[++ihist]);
  if (fHistStat         ) fHistStat         ->Merge(&collections[++ihist]);
    
  
  delete iter;

  return count+1;
}

void AliCollisionNormalization::FillVzMCGen(Float_t vz, Int_t ntrk, AliMCEvent * mcEvt) {

  // Fill MC gen histo and the process types statistics

  Int_t evtype = GetProcessType(mcEvt);
  // When I fill gen histos, I also fill statistics of process types (used to detemine ratios afterwards).
  fHistProcTypes->Fill(evtype);
  fHistVzMCGen[evtype]->Fill(vz,ntrk);
}

void AliCollisionNormalization::FillVzMCRec(Float_t vz, Int_t ntrk, AliMCEvent * mcEvt){

  // Fill MC rec histo
  Int_t evtype = GetProcessType(mcEvt);
  fHistVzMCRec[evtype]->Fill(vz,ntrk);

}
void AliCollisionNormalization::FillVzMCTrg(Float_t vz, Int_t ntrk, AliMCEvent * mcEvt) {      

  // Fill MC trg histo
  Int_t evtype = GetProcessType(mcEvt);
  fHistVzMCTrg[evtype]->Fill(vz,ntrk);
}


Int_t AliCollisionNormalization::GetProcessType(const AliMCEvent * mcEvt) {

  // Determine if the event was generated with pythia or phojet and return the process type

  // Check if mcEvt is fine
  if (!mcEvt) {
    AliFatal("NULL mc event");
  } 

  // Determine if it was a pythia or phojet header, and return the correct process type
  AliGenPythiaEventHeader * headPy  = 0;
  AliGenDPMjetEventHeader * headPho = 0;
  AliGenEventHeader * htmp = mcEvt->GenEventHeader();
  if(!htmp) {
    AliFatal("Cannot Get MC Header!!");
  }
  if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
    headPy =  (AliGenPythiaEventHeader*) htmp;
  } else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
    headPho = (AliGenDPMjetEventHeader*) htmp;
  } else {
    AliError("Unknown header");
  }

  // Determine process type
  if(headPy)   {
    if(headPy->ProcessType() == 92 || headPy->ProcessType() == 93) {
      // single difractive
      return kProcSD;
    } else if (headPy->ProcessType() == 94) {
      // double diffractive
      return kProcDD;
    }
    else if(headPy->ProcessType() != 92 && headPy->ProcessType() != 93 && headPy->ProcessType() != 94) {    
      // non difractive
      return kProcND; 
    }
  } else if (headPho) {
    if(headPho->ProcessType() == 5 || headPho->ProcessType() == 6 ) {
      // single difractive
      return kProcSD;
    } else if (headPho->ProcessType() == 7) { 
      // double diffractive
      return kProcDD;      
    } else if(headPho->ProcessType() != 5 && headPho->ProcessType() != 6  && headPho->ProcessType() != 7 ) {
      // non difractive
      return kProcND; 
    }       
  }
  

  // no process type found?
  AliError(Form("Unknown header: %s", htmp->IsA()->GetName()));
  return kProcUnknown;
}

Double_t AliCollisionNormalization::GetProcessWeight(Int_t proc) {

  // Return a weight to be used when combining the MC histos to
  // compute efficiency, defined as the ratio XS in generator / XS
  // measured

  Float_t refSD,  refDD,  refND,  errorSD,  errorDD,  errorND;
  GetRelativeFractions(fReferenceXS,refSD,  refDD,  refND,  errorSD,  errorDD,  errorND);

  static Double_t total = fHistProcTypes->Integral();
  if (fHistProcTypes->GetBinContent(fHistProcTypes->FindBin(kProcUnknown)) > 0) {
    AliError("There were unknown processes!!!");
  }
  static Double_t lSD = fHistProcTypes->GetBinContent(fHistProcTypes->FindBin(kProcSD))/total;
  static Double_t lDD = fHistProcTypes->GetBinContent(fHistProcTypes->FindBin(kProcDD))/total;
  static Double_t lND = fHistProcTypes->GetBinContent(fHistProcTypes->FindBin(kProcND))/total;

  if (fVerbose > 2) {
    AliInfo(Form("Total MC evts: %f",total));
    AliInfo(Form(" Frac SD %4.4f", lSD));
    AliInfo(Form(" Frac DD %4.4f", lDD));
    AliInfo(Form(" Frac ND %4.4f", lND));
  }
  
  switch(proc) {
  case kProcSD:
    return refSD/lSD;
    break;
  case kProcDD:
    return refDD/lDD;
    break;
  case kProcND:
    return refND/lND;
    break;
  default:
    AliError("Unknown process");
  }
  
  return 0;

} 

void AliCollisionNormalization::GetRelativeFractions(Int_t origin, Float_t& refSD, Float_t& refDD, Float_t& refND, Float_t& errorSD, Float_t& errorDD, Float_t& errorND)
{
  // Returns fraction of XS (SD, ND, DD) and corresponding error
  // Stolen from Jan Fiete's drawSystematics macro

  // origin: 
  //   -1 = Pythia (test)
  //   0 = UA5
  //   1 = Data 1.8 TeV
  //   2 = Tel-Aviv
  //   3 = Durham
  //

  Double_t epsilon = 0.0001;
  if(TMath::Abs(fEnergy-900)<epsilon) {

    switch (origin)
      {
      case -10: // Pythia default at 7 GeV, 50% error
	AliInfo("PYTHIA x-sections");
	refSD = 0.192637; errorSD = refSD * 0.5;
	refDD = 0.129877; errorDD = refDD * 0.5;
	refND = 0.677486; errorND = 0;
	break;
	
      case -1: // Pythia default at 900 GeV, as test
	AliInfo("PYTHIA x-sections");
      refSD = 0.223788;
      refDD = 0.123315;
      refND = 0.652897;
      break;
      
      case 0: // UA5
	AliInfo("UA5 x-sections a la first paper");
	refSD = 0.153; errorSD = 0.05;
	refDD = 0.080; errorDD = 0.05;
	refND = 0.767; errorND = 0;
	break;
	
      case 10: // UA5
	AliInfo("UA5 x-sections hadron level definition for Pythia"); 
	// Fractions in Pythia with UA5 cuts selection for SD
	// ND: 0.688662
	// SD: 0.188588 --> this should be 15.3
	// DD: 0.122750
	refSD = 0.224 * 0.153 / 0.189; errorSD = 0.023 * 0.224 / 0.189;
	refDD = 0.095;                 errorDD = 0.06; 
	refND = 1.0 - refSD - refDD; errorND = 0;
	break;
	
      case 11: // UA5
	AliInfo("UA5 x-sections hadron level definition for Phojet"); 
	// Fractions in Phojet with UA5 cuts selection for SD
	// ND: 0.783573
	// SD: 0.151601 --> this should be 15.3
	// DD: 0.064827
	refSD = 0.191 * 0.153 / 0.152; errorSD = 0.023 * 0.191 / 0.152;
	refDD = 0.095;                 errorDD = 0.06; 
	refND = 1.0 - refSD - refDD; errorND = 0;
	break;
      case 2: // tel-aviv model
	AliInfo("Tel-aviv model x-sections");
	refSD = 0.171;
	refDD = 0.094;
	refND = 1 - refSD - refDD;
	break;
	
      case 3: // durham model
	AliInfo("Durham model x-sections");
	refSD = 0.190;
	refDD = 0.125;
	refND = 1 - refSD - refDD;
      break;
      default:
	AliFatal(Form("Unknown origin %d, Energy %f", origin, fEnergy));
      }
  }
  else if(TMath::Abs(fEnergy-1800)<epsilon) {     
    switch (origin)
      {
      case 20: // E710, 1.8 TeV
      AliInfo("E710 x-sections hadron level definition for Pythia");
      // ND: 0.705709
      // SD: 0.166590 --> this should be 15.9
      // DD: 0.127701
      refSD = 0.217 * 0.159 / 0.167; errorSD = 0.024 * 0.217 / 0.167;
      refDD = 0.075 * 1.43;          errorDD = 0.02 * 1.43; 
      refND = 1.0 - refSD - refDD; errorND = 0;
      break;
    
      case 21: // E710, 1.8 TeV
	AliInfo("E710 x-sections hadron level definition for Phojet"); 
	// ND: 0.817462
	// SD: 0.125506 --> this should be 15.9
	// DD: 0.057032
	refSD = 0.161 * 0.159 / 0.126; errorSD = 0.024 * 0.161 / 0.126;
	refDD = 0.075 * 1.43;         errorDD = 0.02 * 1.43;
	refND = 1.0 - refSD - refDD; errorND = 0;
	break;
	
      case 1: // data 1.8 TeV
	AliInfo("??? x-sections");
	refSD = 0.152;
	refDD = 0.092;
	refND = 1 - refSD - refDD;
	break;
      default:
	AliFatal(Form("Unknown origin %d, Energy %f", origin, fEnergy));
      }    
  } 
  else {
    AliFatal(Form("Unknown energy %f", fEnergy));
  }
    
}


// this comment serves no purpose
