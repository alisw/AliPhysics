#include <TCanvas.h>
#include <TProfile.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFileMerger.h>
#include <TKey.h>

#include <AliRawReader.h>
#include <AliRawReaderRoot.h>
#include <AliRawReaderDateOnline.h>
#include <AliCaloRawStreamV3.h>
#include <AliDAQ.h>
#include <AliLog.h>

/// Macro for EMCAL raw data analyisis
/// Instructions:
/// Use offline:
/// run the script copyRawFilesFromAlien.sh with the needed modification of local and alien path
/// per each run: root
///.L countRun2.C+
/// countRun2Chunks(nchunks, runnumb, "15000", "011", "./",  1, 9999999 /*number of events*/, 10 /*cut on signal amplitude*/, typeSel /*  less than 0 = take everything; physicsEvent=7, calibrationEvent=8 ($ALICE_ROOT/RAW/event.h)*/, 1, 1,"EMCAL")
///
/// PlotMergedOutput(nchunks, runnumb, "Amp10", "", "./")
/// "Amp10": vary the number to have different cuts on lower bound of the amplitude recorded
/// 
/// MergeTTree(nchunks, runnumb, "15000", "011", "")
/// The tree can be plotted with ReadTreeCountRun2.C
///
/// Use online:
/// .L countRun2.C+
/// countRun2(2, "@gdc-EMCal-00:")

// config for LHC Run 2 EMCal data (2015-): 
// 20 SuperModules total
const int kNSMEMCal = 12;
const int kNSMDCal = 8;
const int kNSM = 20;

// 40 optical fibers/DDLs total; 2 per SuperModule
// 2 of the 1/3 SMs (A5 and C12) actually don't have 2 DDLs, but space is reserved for them anyhow (the missing DDLs are # 21 and 39)
const int kNDDLEMCal = 24; // 0..23
const int kNDDLDCal = 16; // 24..39
const int kNDDL = 40; 
// There are also 2 DDLs (and 2 LDCs) with trigger=STU data; this macro ignores the STU DDLs

// 12 LDCs with EMCal/DCal data
const int kNLDCEMCal = 8; // aldaqpc137..144 a.k.a. ldc-EMCAL-C-0..3 ldc-EMCAL-A-0..3 ; 3 DDLs per LDC except for the last one   
const int kNLDCDCal = 5; // aldaqpc146..150 a.k.a. ldc-DCAL-0..4; 3 DDLs per LDCs
const int kNLDC = 13; 


// how many channels per DDL and SuperModule
const int kNCHAN = 1152; // FEE channels, per DDL
const int kNCHANTRU = 96; // channels per TRU
const int kNCHANLEDREF = 48; // channels for LED Mon Ref. 

const int kNCHANDDL0 = kNCHAN + kNCHANTRU + kNCHANLEDREF;
const int kNCHANDDL1 = kNCHAN + 2*kNCHANTRU;
 
const int kNCHANTot = kNCHANDDL0 + kNCHANDDL1; // per SuperModule

const int kMARGIN = 50; // +/- bin margin for TH2F plotting
 
void countRun2(Int_t onoffline = 1 /*1 = offline, needs a root file; 2 = online, need a GDC name*/, 
   	const char * gdcNameStr = "15000217514011.12.root" /*"@gdc-GLOBAL-01:"*/,
	const int kMaxEvents=100,
	const int ampCut = 5, // only save channel info when they have at least this much signal
	const int typeSel = 7, // >0=do match; <0=take everything; physicsEvent=7, calibrationEvent=8 ($ALICE_ROOT/RAW/event.h)
	// flag for selecting verbosity
	const int debug = 1, // -1=no output, 0=quiet. 1=some printouts, 2=lot of printouts..
	const int saveTree = 1, // flag if we save the TTree or not
	const TString calo="EMCAL",
	Int_t optchunk=-1)
{
  TH2F *hDDLNchan  = new TH2F("hDDLNchan","DDL id vs Nchan;DDL;Nchan",
			      kNDDL, -0.5, kNDDL - 0.5, 
			      100, -kMARGIN, kNCHANDDL1 + kMARGIN);

  TH2F *hSMNchan  = new TH2F("hSMNchan","SM id vs Nchan;SM;Nchan",
			      kNSM, -0.5, kNSM - 0.5,
			      100, -kMARGIN, kNCHANTot + kMARGIN);

  TH2F *hLDCNchan  = new TH2F("hLDCNchan","LDC id vs Nchan;LDC;Nchan",
			      kNLDC, -0.5, kNLDC - 0.5,
			      100, -kMARGIN, kNCHANTot + kMARGIN);

  TProfile *hDDLNsamp  = new TProfile("hDDLNsamp","DDL id vs Nsamp;DDL;Nsamp",
				      kNDDL, -0.5, kNDDL - 0.5); 

  TProfile *hSMNsamp  = new TProfile("hSMNsamp","SM id vs Nsamp;SM;Nsamp",
				      kNSM, -0.5, kNSM - 0.5);

  TProfile *hLDCNsamp  = new TProfile("hLDCNsamp","LDC id vs Nsamp;LDC;Nsamp",
				      kNLDC, -0.5, kNLDC - 0.5);


  TFile destFile(Form("output_tree_%s",gdcNameStr), "recreate");  

  Int_t evno = 0;
  UInt_t type = 0;
  UInt_t period = 0;
  UInt_t orbit = 0;
  UInt_t bc = 0;

  TTree *treeLDC = new TTree("treeLDC","");
  treeLDC->SetAutoSave(1000000);
  // variables for treeLDC
  Int_t nLDC = kNLDC;
  Int_t nLDCChan[kNLDC] = {0};
  Int_t nLDCSamp[kNLDC] = {0};
  // declare the branches
  treeLDC->Branch("evno", &evno, "evno/I");
  treeLDC->Branch("type", &type, "type/i");
  treeLDC->Branch("period", &period, "period/i");
  treeLDC->Branch("orbit", &orbit, "orbit/i");
  treeLDC->Branch("bc", &bc, "bc/i");
  treeLDC->Branch("nLDC", &nLDC, "nLDC/I");
  treeLDC->Branch( "nLDCChan", &nLDCChan, Form("nLDCChan[nLDC]/I") );
  treeLDC->Branch( "nLDCSamp", &nLDCSamp, Form("nLDCChan[nLDC]/I") );

  TTree *treeSM = new TTree("treeSM","");
  treeSM->SetAutoSave(1000000);
  // variables for treeSM
  Int_t nSM = kNSM;
  Int_t nSMChan[kNSM] = {0};
  Int_t nSMSamp[kNSM] = {0};
  // declare the branches
  treeSM->Branch("evno", &evno, "evno/I");
  treeSM->Branch("type", &type, "type/i");
  treeSM->Branch("period", &period, "period/i");
  treeSM->Branch("orbit", &orbit, "orbit/i");
  treeSM->Branch("bc", &bc, "bc/i");
  treeSM->Branch("nSM", &nSM, "nSM/I");
  treeSM->Branch( "nSMChan", &nSMChan, Form("nSMChan[nSM]/I") );
  treeSM->Branch( "nSMSamp", &nSMSamp, Form("nSMChan[nSM]/I") );

  TTree *treeDDL = new TTree("treeDDL","");
  treeDDL->SetAutoSave(1000000);
  // variables for treeDDL
  Int_t iDDL = 0;
  Int_t nChan = 0;
  Int_t nSamp = 0;
  // maximum number of channels in RCU1
  Int_t hwaddress[kNCHANDDL1] = {0};
  Int_t column[kNCHANDDL1] = {0};
  Int_t row[kNCHANDDL1] = {0};
  Int_t caloflag[kNCHANDDL1] = {0};
  Int_t nbunches[kNCHANDDL1] = {0};
  Int_t nsamples[kNCHANDDL1] = {0};
  Int_t min[kNCHANDDL1] = {0};
  Int_t max[kNCHANDDL1] = {0};
  Int_t timeAtMax[kNCHANDDL1] = {0};
  // declare the branches
  treeDDL->Branch("evno", &evno, "evno/I");
  treeDDL->Branch("type", &type, "type/i");
  treeDDL->Branch("period", &period, "period/i");
  treeDDL->Branch("orbit", &orbit, "orbit/i");
  treeDDL->Branch("bc", &bc, "bc/i");
  treeDDL->Branch("iDDL", &iDDL, "iDDL/I");
  treeDDL->Branch("nChan", &nChan, "nChan/I");
  treeDDL->Branch("nSamp", &nSamp, "nSamp/I");
  treeDDL->Branch( "hwaddress", &hwaddress, Form("hwaddress[nChan]/I") );
  treeDDL->Branch( "column", &column, Form("column[nChan]/I") );
  treeDDL->Branch( "row", &row, Form("row[nChan]/I") );
  treeDDL->Branch( "caloflag", &caloflag, Form("caloflag[nChan]/I") );
  treeDDL->Branch( "nbunches", &nbunches, Form("nbunches[nChan]/I") );
  treeDDL->Branch( "nsamples", &nsamples, Form("nsamples[nChan]/I") );
  treeDDL->Branch( "min", &min, Form("min[nChan]/I") );
  treeDDL->Branch( "max", &max, Form("max[nChan]/I") );
  treeDDL->Branch( "timeAtMax", &timeAtMax, Form("timeAtMax[nChan]/I") );

  cout << " tree's created " << endl; 

  // input
  TString rootreadername = Form("%s",gdcNameStr);
  AliRawReader *reader = NULL;
  if(onoffline==1) {
     if(!rootreadername.Contains(".root")) {
     	Printf("%s is not a root file, did you want online mode? -> onoffline = 2", gdcNameStr);
     	return;
     }
     reader = new AliRawReaderRoot( rootreadername );
  }
  if(onoffline==2) {
     if(!rootreadername.Contains("gdc")) {
     	Printf("%s is not a GDC, did you want offline mode? -> onoffline = 1", gdcNameStr);
     	return;
     }
     reader = new AliRawReaderDateOnline( rootreadername );
  }
  if(!reader) {
     Printf("Select online/offline mode 1 or 2");
     return;
  }
  reader->Reset();
  AliCaloRawStreamV3 *stream = new AliCaloRawStreamV3(reader,calo);
  reader->Select("EMCAL", 0, kMaxEvents);

  cout << " reader and stream created " << endl; 
  
  Int_t iev = 0;
  int ircu = -1;

  // variables used in loop over samples below
  int i = 0;
  int sample = 0;
  int time = 0;
  int startBin = 0;
  int runno = 0;

  Int_t emcID = AliDAQ::DetectorID("EMCAL"); // bit 18..
  printf("Max events %d\n",kMaxEvents);
  while (reader->NextEvent() && iev<kMaxEvents) {

    type = reader->GetType();
    if (debug >= 0) {
      printf("event type %d\n",type);
    }

    // make sure EMCal was readout during the event
    const UInt_t *detPattern = reader->GetDetectorPattern(); 
    UInt_t emcInReadout = ( ((1 << emcID) & detPattern[0]) >> emcID);

    if (debug >= 0) {
      printf("emcInReadout %d typeSel %d\n",emcInReadout, typeSel);
    }

    if ( emcInReadout && (typeSel<0 || typeSel==type) ) { // event selection

    evno = iev; // counter among the events read
    if (iev>=kMaxEvents) 
    {
      printf("Stop\n");
      break;
    }
    iev++;

    runno = reader->GetRunNumber();
    period = reader->GetPeriod(); // seems to always be 0, so we may not need to store this?
    orbit = reader->GetOrbitID();
    bc = reader->GetBCID();
      
    if (debug >= 0) {
      printf("Reading event %02d :run %d  type %d period %d orbit %d bc %d\n",
	     runno, evno, type, period, orbit, bc);
    }

    // reset LDC counters
    int ildc = 0;
    int ism = 0;
    for (ildc=0; ildc<kNLDC; ildc++) {
      nLDCChan[ildc] = 0;
      nLDCSamp[ildc] = 0;
    }

    // reset SM counters
    for (ism=0; ism<kNSM; ism++) {
      nSMChan[ism] = 0;
      nSMSamp[ism] = 0;
    }

    while (stream->NextDDL()) {
      // reset DDL counters
      nChan = 0;
      nSamp = 0;
      for (int ichan=0; ichan<kNCHANDDL1; ichan++) {
	hwaddress[ichan] = 0;
	column[ichan] = 0;
	row[ichan] = 0;
	caloflag[ichan] = 0;
	nbunches[ichan] = 0;
	nsamples[ichan] = 0;
	min[ichan] = 1023;
	max[ichan] = 0;
	timeAtMax[ichan] = 0;
      }

      iDDL = stream->GetDDLNumber();
      printf("iDDL = %d\n",iDDL);
      //     continue;
      ism = iDDL / 2; 
      ircu = iDDL % 2; // crate, within SM, DDL0 or DDL1
      int iside = ism % 2; // A=0, C=1
      int isector = ism / 2; // 0..5 + 6..9 (later sparsified for DCal 9..12 SMs)
      int iddlSide = isector*2 + ircu; // index for either A or C side

      // LDC business: lots of magic numbers below for the current (Mar 2015) configuration
      ildc = 0;

      if (iDDL < kNDDLEMCal) { // EMCal, LDCs are split - separate for A and C sides...
	if (iside==1) { // C side
	  ildc = iddlSide / 3;  
	}
	else { // A side
	  ildc = iddlSide / 3 + kNLDCEMCal/2; 
	}
      }
      else { // DCal    
	if (iside==1) { // C side
	  ildc = iddlSide / 3 + kNLDCEMCal/2;    
	}
	else { // A side, shared first DDL w C side LDC
	  int iddlCorr = iddlSide - 1;
	  ildc = iddlCorr / 3 + kNLDCEMCal - 1;      
	}
      }
      // end of LDC magic number section
      if (debug>0) {
	printf("iDDL %d ism %d ildc %d\n", 
	       iDDL, ism, ildc);
      }

      while (stream->NextChannel()) {
	if (debug>1) {
	  printf("RCU %d DDLmod2 %d branch %d\n", 
	       stream->GetRCUId(), ircu, stream->GetBranch());
	}

	//	printf( "stream->GetHWAddress() =%d \n", stream->GetHWAddress());

	if (stream->GetHWAddress() > 3279) {
	   printf("too large hwaddr: %d RCU %d DDLmod2 %d branch %d\n", 
	      stream->GetHWAddress(),
	      iDDL, ircu, stream->GetBranch());
	   
	   continue;
	}

	hwaddress[nChan] = stream->GetHWAddress();
	column[nChan] = stream->GetColumn();
	row[nChan] = stream->GetRow();
	caloflag[nChan] = stream->GetCaloFlag();
	// counters for this channel
	nsamples[nChan] = 0;
	nbunches[nChan] = 0;
	min[nChan] = 1023;
	max[nChan] = 0;

	while (stream->NextBunch()) {
	  nsamples[nChan] += stream->GetBunchLength();
	  nbunches[nChan]++;

	  // loop over samples; check for min/max
	  const UShort_t *sig = stream->GetSignals();
	  startBin = stream->GetStartTimeBin();
	  for (i = 0; i < stream->GetBunchLength(); i++) {
	    sample = sig[i];
	    time = startBin--;

	    // check if it's a min or max value
	    if (sample < min[nChan]) {
	      min[nChan] = sample;
	    }
	    if (sample > max[nChan]) {
	      max[nChan] = sample;
	      timeAtMax[nChan] = time;
	    }

	  } // loop over samples in bunch
	} // loop over bunches

	if ( (max[nChan]-min[nChan]) > ampCut ) {
	  nSamp += nsamples[nChan];
	  nChan++;
	}

      } // channel loop
      treeDDL->Fill();

      if (debug>0) {
	printf("DDL %d nchan %d nsamp %d\n", iDDL, nChan, nSamp);
      }

      nLDCChan[ildc] += nChan;
      nLDCSamp[ildc] += nSamp;
      nSMChan[ism] += nChan;
      nSMSamp[ism] += nSamp;
      hDDLNchan->Fill(iDDL, nChan);
      hDDLNsamp->Fill(iDDL, nSamp);

    } // DDL loop

    for (ildc=0; ildc<kNLDC; ildc++) {
      if (debug>0) {
	printf("LDC %d nchan %d nsamp %d\n", ildc, nLDCChan[ildc], nLDCSamp[ildc]);
      }
      hLDCNchan->Fill(ildc, nLDCChan[ildc]);
      hLDCNsamp->Fill(ildc, nLDCSamp[ildc]);
    }

    for (ism=0; ism<kNSM; ism++) {
      if (debug>0) {
	printf("SM %d nchan %d nsamp %d\n", ism, nSMChan[ism], nSMSamp[ism]);
      }
      hSMNchan->Fill(ism, nSMChan[ism]);
      hSMNsamp->Fill(ism, nSMSamp[ism]);
    }

    treeLDC->Fill();
    treeSM->Fill();

    } // typeSel
    stream->Reset();
  } // event loop

  if (saveTree) {  // save output tree
    destFile.cd();
    treeLDC->Write();
    treeSM->Write();
    treeDDL->Write();
    destFile.Close();
    //    char cmd[200];
    //    sprintf(cmd, ".! mv root/run.root root/run%09d.root", runno);
    //    gROOT->ProcessLine(cmd);
  }

  // some plotting stuff
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  gROOT->ProcessLine(".! mkdir -p png");
  Char_t outname[256];
  TCanvas *c1 = new TCanvas("c1","",0,0,900,600);
  hDDLNchan->Draw("colz");
  sprintf(outname, "png/DDLNchan_%09d.png",runno);
  c1->SaveAs(outname);

  TCanvas *c2 = new TCanvas("c2","",100,100,900,600);
  hLDCNchan->Draw("colz");
  sprintf(outname, "png/LDCNchan_%09d.png",runno);
  c2->SaveAs(outname);

  TCanvas *c3 = new TCanvas("c3","",200,200,900,600);
  hDDLNsamp->SetMarkerStyle(20);
  hDDLNsamp->Draw("colz");
  sprintf(outname, "png/DDLNsamp_%09d.png",runno);
  c3->SaveAs(outname);

  TCanvas *c4 = new TCanvas("c4","",300,300,900,600);
  
  hLDCNsamp->SetMarkerStyle(20);
  hLDCNsamp->Draw("colz");
  sprintf(outname, "png/LDCNsamp_%09d.png",runno);
  c4->SaveAs(outname);

  TCanvas *c5 = new TCanvas("c5","",100,100,900,600);
  hSMNchan->Draw("colz");
  sprintf(outname, "png/SMNchan_%09d.png",runno);
  c5->SaveAs(outname);

  TCanvas *c6 = new TCanvas("c6","",300,300,900,600);  
  hSMNsamp->SetMarkerStyle(20);
  hSMNsamp->Draw("colz");
  sprintf(outname, "png/SMNsamp_%09d.png",runno);
  c6->SaveAs(outname);
  
  TString sufopt = optchunk>0 ? Form("%d", optchunk) : "";
  TFile *fout = new TFile(Form("outputRun%dAmp%d.%s.root", runno, ampCut, sufopt.Data()), "recreate");
  fout->cd();
  hLDCNchan->Write();
  hDDLNsamp->Write();
  hLDCNsamp->Write();
  hSMNchan ->Write();
  hSMNsamp ->Write();

}

void countRun2Chunks(Int_t nchunks, Int_t runnumb, TString filenamebase = "15000", TString endofname = "011",TString pathtofile="", Int_t onoffline = 1, const int kMaxEvents=999999,
const int ampCut = 5, // only save channel info when they have at least this much signal
const int typeSel = 7, // >0=do match; <0=take everything; physicsEvent=7, calibrationEvent=8 ($ALICE_ROOT/RAW/event.h)	// flag for selecting verbosity
	const int debug = 1, // -1=no output, 0=quiet. 1=some printouts, 2=lot of printouts..
	const int saveTree = 1, // flag if we save the TTree or not
	const TString calo="EMCAL"){

   for(Int_t ifiles = 0 ; ifiles < nchunks ; ifiles++){
      TString filename = Form("%s%s%d%s.%d.root", pathtofile.Data(),filenamebase.Data(), runnumb, endofname.Data(), ifiles+10);
      
   
      countRun2(onoffline, filename, kMaxEvents, ampCut, typeSel, debug, saveTree,calo, ifiles+10);
      
   }

}

void MergeFiles(Int_t nfiles, Int_t runnumb, TString filenamebase = "output_tree_"/*"outputRun"*/, TString startofname = "15000"/*""*/, TString endofname = "011"/*"Amp10"*/,TString pathtofile=""){
   //download the files from alien using the script
   TFileMerger merge;
   merge.OutputFile(Form("%s%s%s%d%s.root", pathtofile.Data(),filenamebase.Data(), startofname.Data(), runnumb, endofname.Data()));   
   for(Int_t i=0; i<nfiles; i++){
      merge.AddFile(Form("%s%s%s%d%s.%d.root", pathtofile.Data(),filenamebase.Data(), startofname.Data(), runnumb, endofname.Data(), i+10));
         
   }
   merge.Merge();
   return;
}

void PlotMergedOutput(Int_t nfiles, Int_t runnumb, TString endofname = "Amp5", TString startofname = "", TString pathtofile=""){
   TString filenamebase = "outputRun";
   TFile *file = new TFile(Form("%s%s%s%d%s.root", pathtofile.Data(),filenamebase.Data(), startofname.Data(), runnumb, endofname.Data()));
   if(!file->IsOpen()){
      MergeFiles(nfiles, runnumb, filenamebase, startofname, endofname, pathtofile);
      
      file=TFile::Open(Form("%s%s%s%d%s.root", pathtofile.Data(),filenamebase.Data(), startofname.Data(), runnumb, endofname.Data()));
      if(!file->IsOpen()) {
      	 Printf("Problem in merging");
      	 return;
      }
   } else {
      Printf("File exist, already merged");
   }
   
   TList *list = file->GetListOfKeys();
   list->ls();
   Int_t n = list->GetEntries();
   Printf("Draw %d items", n);
   Char_t outname[256];
   for(Int_t i = 0; i< n; i++){
      TCanvas *c = new TCanvas(Form("c%d",i), Form("c%d",i), 800, 800);
      
      TObject *obj = ((TKey*)list->At(i))->ReadObj();
      TString classname = obj->IsA()->GetName();
      printf("%s, ", classname.Data());
      
      if( classname.Contains("TH1") || classname.Contains("TProfile")){
      	 TH1F *h = (TH1F*)obj;
      	 c->cd();
      	 h->Draw();
      
      }
      if(classname.Contains("TH2")){
      	 TH2F *h = (TH2F*)obj;
      	 c->cd();
      	 h->Draw("colz");
      
      }
      gROOT->ProcessLine(".! mkdir -p png");
      sprintf(outname, "png/%s%d.png",obj->GetName(),runnumb);
      c->SaveAs(outname);
      Printf("%s",((TKey*)list->At(i))->ReadObj()->GetName());
     
   }

}

void MergeTTree(Int_t nfiles, Int_t runnumb, TString startofname = "15000", TString endofname = "011",  TString pathtofile=""){
   
   TString filenamebase = "output_tree_";
   TFile *file = new TFile(Form("%s%s%s%d%s.root", pathtofile.Data(),filenamebase.Data(), startofname.Data(), runnumb, endofname.Data()));
   if(!file->IsOpen()){
      MergeFiles(nfiles, runnumb, filenamebase, startofname, endofname, pathtofile);
      file=TFile::Open(Form("%s%s%s%d%s.root", pathtofile.Data(),filenamebase.Data(), startofname.Data(), runnumb, endofname.Data()));

      if(!file->IsOpen()) {
      	 Printf("Problem in merging");
      	 return;
      }

   } else {
      Printf("File exist, already merged");
   }
   return;

}
