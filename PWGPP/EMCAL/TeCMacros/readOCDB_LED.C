#if !defined( __CINT__) || defined(__MAKECINT__)
#include <AliCDBEntry.h>
#include <AliCDBEntry.h>
#include <AliCDBManager.h>
#include <AliCDBMetaData.h>
#include <AliCaloCalibSignal.h>
#include <AliDCSValue.h>
#include <AliEMCALPreprocessor.h>
#include <AliEMCALSensorTempArray.h>
#include <AliEMCALTriggerDCSConfig.h>
#include <AliEMCALTriggerSTUDCSConfig.h>
#include <AliEMCALTriggerTRUDCSConfig.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TF1.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMap.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include "LInfo.h"
#endif

LInfo *readOCDB_LED(Int_t runNb  = 286350, Bool_t debug=1)
{
  TGrid::Connect("alien://");

  AliCDBManager*  cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);

  AliCDBEntry *en = 0;
  try {
    en=cdb->Get("EMCAL/Calib/LED/");
  } catch (...) { ; }
  if (!en) {
    printf("no entry found for run %d!\n", runNb);
    return 0;
  }

  AliCaloCalibSignal *emcCalibSignal=dynamic_cast<AliCaloCalibSignal*>(en->GetObject());
  if (!emcCalibSignal) {
    printf("can not cast to calib signal!\n");
    return 0;
  }
  //emcCalibSignal->Print();

  LInfo *ret = new LInfo(runNb);

  TTree *treeAmp = emcCalibSignal->GetTreeAvgAmpVsTime();
  TTree *treeLEDAmp = emcCalibSignal->GetTreeLEDAvgAmpVsTime();
  Int_t tea = treeAmp->GetEntries();
  Int_t tem = treeLEDAmp->GetEntries();
  if (debug) {
    cout << " Tree entries " << tea << endl;
    cout << " TreeLEDMon entries " << tem << endl;
  }

  Int_t    fNum     = 0;
  Double_t fAvgAmp  = 0;
  Double_t fHour    = 0;
  Double_t fRMS     = 0;
  treeAmp->SetBranchAddress("fChannelNum",&fNum);
  treeAmp->SetBranchAddress("fHour",&fHour);
  treeAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
  treeAmp->SetBranchAddress("fRMS",&fRMS);

  // associate variables  also for LED
  treeLEDAmp->SetBranchAddress("fRefNum",&fNum);
  treeLEDAmp->SetBranchAddress("fHour",&fHour);
  treeLEDAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
  treeLEDAmp->SetBranchAddress("fRMS",&fRMS);

  Int_t mod   = 0;
  Int_t col   = 0;
  Int_t row   = 0;
  Int_t gain  = 0;
  Int_t strip = 0;

  for (Int_t i=0; i<tem; ++i) {
    treeLEDAmp->GetEntry(i);
    emcCalibSignal->DecodeRefNum(fNum, &mod, &strip, &gain);
    if (debug) 
      Printf("fRefNum = %d, mod = %d, strip = %d, amp = %f, gain = %d", fNum, mod, strip, fAvgAmp, gain);
    ret->FillStrip(mod,gain,strip, fAvgAmp, fRMS);
  }

  for (Int_t i=0; i<tea; ++i) {
    treeAmp->GetEntry(i);
    emcCalibSignal->DecodeChannelNum(fNum, &mod, &col, &row, &gain);
    if (debug)
      Printf("fChannelNum = %d, mod = %d, col = %d, amp = %f, gain = %d", fNum, mod, col, fAvgAmp, gain);
    ret->FillLed(mod,gain,col, row, fAvgAmp, fRMS);
  }

  return ret;
}

void testOCDB_LED(Int_t runNb  = 286350)
{
  LInfo *i = readOCDB_LED(runNb,1);
  i->Print();
  return;
}

void read_LHC18d() 
{
  Int_t runs[] = {285978,285979,285980,286014,286018,286025,286026,286027,286030,286064,286124,286127,286129,286130,286154,286157,286159,286198,286201,286202,286203,286229,286230,286231,286254,286255,286256,286257,286258,286261,286263,286282,286284,286287,286288,286289,286308,286309,286310,286311,286312,286313,286314,286336,286337,286340,286341,286345,286348,286349,286350};
  Int_t nruns = sizeof(runs)/sizeof(Int_t);

  TObjArray arr;
  arr.SetOwner(1);

  for (Int_t i=nruns-3;i<nruns;++i) {
    Int_t rn = runs[i];
    cout << i << " " << rn << endl;
    LInfo *ti = readOCDB_LED(rn,1);
    if (!ti) 
      continue;
    ti->Print();
    arr.Add(ti);
  }

  TFile *outf = TFile::Open("ledinfo.root","update");
  arr.Write("led_lhc18d",TObject::kSingleKey);
  outf->ls();
  outf->Close();
}

#if 0
// global EMCal numbers
const Int_t kNSM = 20; 
const Int_t kNSectors = 5;
const Int_t kNSides = 2;

// these other numbers are per SuperModule
const Int_t kNCol = 48;
const Int_t kNRow = 24;
const Int_t kNStrip = kNCol / 2;

// amplitude cut, in ADCs
const Int_t kAmpCut = -1; // -1 => include everything for now..
 
// which gains to plot?
const Int_t towerGain = 1;
const Int_t stripGain = 0; 

const Bool_t reverseCSstrips = kTRUE;

void createHistosFromLEDfiles_orig(const char *listfile = "filelistb.txt")

{
  char name1[13];
  char name3[20];
  Float_t rmsValues[kNSM]={0.};
  Float_t previousRmsValues[kNSM]={0.};
  Int_t referenceRun=219495;
  Int_t desiredFlag=0;
  Int_t flagOfChangeRef;
  Int_t numberOfSMsWithBadStrips;
  Int_t prevflagOfChangeRef=0;
  Int_t prevnumberOfSMsWithBadStrips=10;
  Double_t ledMonAmp;
  Double_t errorLedMonAmp;

  int nruns = 0;
  ifstream fin;
  fin.open(listfile);

  Int_t runno = 0;

  
  char filename[100];
  char filenameout[100];
  char * sideStr[] = {"A", "C"};
  Int_t nBadStripsInSM[kNSM];

  while ( fin.good() ) {
    fin >> runno >> filename;
    if ( fin.good() ) {
      cout << " filename fn " << filename << " runno " << runno << endl;
      nruns++;

      flagOfChangeRef=0;
      numberOfSMsWithBadStrips=0;

      TFile *f = TFile::Open(filename);

      if (!f || !f->IsOpen()) {
	Printf("Cannot Open file");
	
	continue;
      }
      
      Printf("filename %s  is now opened", filename);
    
      AliCaloCalibSignal *emcCalibSignal=(AliCaloCalibSignal*) AliCDBEntry->GetObject();
      //AliCaloCalibSignal1 *emcCalibSignal=(AliCaloCalibSignal1*) f->Get("AliCDBEntry");
      // Printf("calo calib ok");
      TTree *treeAmp = emcCalibSignal->GetTreeAvgAmpVsTime();
      TTree *treeLEDAmp = emcCalibSignal->GetTreeLEDAvgAmpVsTime();
      //Printf("calo calib ok1"<<  "tree pointer  "<<treeAmp<<  "   "<<treeLEDAmp<<endl;);
      if ((treeAmp->GetEntries()<10)||(treeLEDAmp->GetEntries()<10)) continue;
      //Printf("calo calib ok2");
      cout << " Tree entries " << treeAmp->GetEntries() << endl;
      cout << " TreeLEDMon entries " << treeLEDAmp->GetEntries() << endl;

      // book some histograms
      TH1F *hStrip[kNSM][2];
      TH1F *hRevStrip[kNSM][2];
      TH2F *h[kNSM][2];
      TH2F *hAmpOverMon[kNSM][2];
      char id[100];
      char title[100]; 
      TH2F *hCount[kNSM][2];   
      TH1F *hStripCount[kNSM][2];
      TH2F *hsmAmpRmsOverMean[kNSM]; 
      TH1F *hsmLEDAmpRmsOverMean[kNSM];

      for (int iSM=0; iSM<kNSM; iSM++) {

	int isector = iSM/2;
	int iside = iSM%2; 

	sprintf(id, "hsmAmpRmsOverMean%d", iSM);
	sprintf(title, "rms  over mean: SM %d (%1s%d)", 
		iSM, sideStr[iside], isector);
	hsmAmpRmsOverMean[iSM]= new TH2F(id, title, kNCol, -0.5, kNCol-0.5, kNStrip, -0.5, kNStrip-0.5);

	sprintf(id, "hsmLEDAmpRmsOverMean%d", iSM);
	sprintf(title, "rms  over mean: SM %d (%1s%d)", 
		iSM, sideStr[iside], isector);
	hsmLEDAmpRmsOverMean[iSM]= new TH1F(id, title, kNCol, -0.5, kNCol-0.5);

	for (int igain=0; igain<2; igain++) {

	  sprintf(id, "hAmpOverMon%02d%d", iSM, igain);
	  sprintf(title, "Tower LEDAmplitude over mon: SM %d (%1s%d) gain %d", 
		  iSM, sideStr[iside], isector, igain);
	  hAmpOverMon[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5,
					     kNRow, -0.5, kNRow - 0.5);   

	  sprintf(id, "hStrip%02d%d", iSM, igain);
	  sprintf(title, "LEDMon Amplitude: SM %d (%1s%d) gain %d", 
		  iSM, sideStr[iside], isector, igain);
	  hStrip[iSM][igain]= new TH1F(id, title, kNCol, -0.5, kNCol-0.5);

	  sprintf(id, "hStripCount%02d%d", iSM, igain);
	  sprintf(title, "LEDMon Entries: SM %d (%1s%d) gain %d", 
		  iSM, sideStr[iside], isector, igain);
	  hStripCount[iSM][igain]= new TH1F(id, title, kNCol, -0.5, kNCol-0.5);

	  sprintf(id, "hCount%02d%d", iSM, igain);
	  sprintf(title, "Tower Entries: SM %d (%1s%d) gain %d", 
		  iSM, sideStr[iside], isector, igain);
	  hCount[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5,
					kNRow, -0.5, kNRow - 0.5);

	  sprintf(id, "hRevStrip%02d%d", iSM, igain);
	  hRevStrip[iSM][igain]= new TH1F(id, title, kNStrip, -0.5, kNStrip-0.5);

	  sprintf(id, "h%02d%d", iSM, igain);
	  sprintf(title, "Tower Amplitude: SM %d (%1s%d) gain %d", 
		  iSM, sideStr[iside], isector, igain);
	  h[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5,
				   kNRow, -0.5, kNRow - 0.5);

	}
      }

      // fill the histograms; loop over trees
      int fChannelNum = 0; // for regular towers
      int fChannelN = 0; // for regular towers
      int fRefNum = 0; // for LED
      double fAvgAmp = 0;
      double fHour = 0;
      double fRMS=0;
      double totalChannelAmp=0.;
      double totalChannelRms=0.;
      int firstChannelNum=0, nentry=0, nsame=0;

      int mod = 0;
      int col = 0;
      int row = 0;
      int gain = 0;
      int strip = 0;
      int mod1 = 0;
      int col1 = 0;
      int row1 = 0;
      int gain1 = 0;
      int strip1 = 0;

      // associate variables 
      treeAmp->SetBranchAddress("fChannelNum",&fChannelNum);
      treeAmp->SetBranchAddress("fHour",&fHour);
      treeAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
      treeAmp->SetBranchAddress("fRMS",&fRMS);

      // associate variables  also for LED
      treeLEDAmp->SetBranchAddress("fRefNum",&fRefNum);
      treeLEDAmp->SetBranchAddress("fHour",&fHour);
      treeLEDAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
      treeLEDAmp->SetBranchAddress("fRMS",&fRMS);

      for (int i=0; i<treeLEDAmp->GetEntries(); i++) {
	nentry=nentry+nsame;
	if (nentry>treeLEDAmp->GetEntries()) continue;
	//treeLEDAmp->GetEntry(nentry);
	treeLEDAmp->GetEntry(i);
	nsame=0;
	emcCalibSignal->DecodeRefNum(fRefNum, &mod, &strip, &gain);

	Printf("fRefNum = %d, mod = %d, strip = %d, amp = %f, gain = %d", fRefNum, mod, strip, fAvgAmp, gain);
	//if(mod<10) continue;
	// if (fAvgAmp > kAmpCut) {
	//   totalChannelAmp=fAvgAmp;
	//   totalChannelRms=fRMS;
	//   firstChannelNum=fRefNum;
	// }
	//   else continue;
	// nsame=1;
	// for (int j=nentry+1; j<treeLEDAmp->GetEntries(); j++) {
	//   treeLEDAmp->GetEntry(j);
	//   emcCalibSignal->DecodeRefNum(fRefNum, &mod1, &strip1, &gain1);
	//   if (fRefNum>firstChannelNum) break;
	//   if ((fAvgAmp > kAmpCut)&&(fRefNum==firstChannelNum)&&(gain==gain1)) {
	// 	totalChannelAmp=totalChannelAmp+fAvgAmp;
	// 	totalChannelRms=totalChannelRms+fRMS;
	// 	nsame=nsame+1;
	//   }
	// }
	// hStrip[mod][gain]->Fill(2*strip, totalChannelAmp/nsame);
	// hStrip[mod][gain]->Fill(2*strip+1, totalChannelAmp/nsame);
	// hStripCount[mod][gain]->Fill(2*strip, totalChannelRms/nsame);
	// hStripCount[mod][gain]->Fill(2*strip+1, totalChannelRms/nsame);

	hStrip[mod][gain]->Fill(2*strip, fAvgAmp);
	hStrip[mod][gain]->Fill(2*strip+1, fAvgAmp);
	hStripCount[mod][gain]->Fill(2*strip, fRMS);
	hStripCount[mod][gain]->Fill(2*strip+1, fRMS);
      }

      nentry=0;
      nsame=0;
      totalChannelAmp=0.;
      totalChannelRms=0.;
      firstChannelNum=0;
      //Printf("AMP");
      for (int i=0; i<treeAmp->GetEntries(); i++) {
	nentry=nentry+nsame;
	//if (nentry>treeAmp->GetEntries()) continue;
	//treeAmp->GetEntry(nentry);
	treeAmp->GetEntry(i);
	// nsame=0;
	emcCalibSignal->DecodeChannelNum(fChannelNum, &mod, &col, &row, &gain);
	Printf("amp fChannelNum = %d, mod = %d, col = %d, amp = %f, gain = %d", fChannelNum, mod, col, fAvgAmp, gain);
	// if(mod>9) continue;
	// if (fAvgAmp > kAmpCut) {
	//   // decode channel num to more usual indices
	//   totalChannelAmp=fAvgAmp;
	//   totalChannelRms=fRMS;
	//   firstChannelNum=fChannelNum;
	// }
	// else continue;
	// nsame=1;
	Double_t weightf;
	// for (int j=nentry+1; j<treeAmp->GetEntries(); j++) {
	//   treeAmp->GetEntry(j);
	//   emcCalibSignal->DecodeChannelNum(fChannelNum, &mod1, &col1, &row1, &gain1);
	//   if (fChannelNum>firstChannelNum) break;
	//   if ((fAvgAmp > kAmpCut)&&(fChannelNum==firstChannelNum)&&(gain==gain1)) {
	// 	totalChannelAmp=totalChannelAmp+fAvgAmp;
	// 	totalChannelRms=totalChannelRms+fRMS*fRMS;
	// 	nsame=nsame+1;
	//   }
	// }

	// if (mod%2==0) {
	ledMonAmp=hStrip[mod][gain]->GetBinContent(col+1);
	errorLedMonAmp=hStripCount[mod][stripGain]->GetBinContent(col+1);
	// }
	// if (mod%2==1) {
	//   ledMonAmp=hStrip[mod][gain]->GetBinContent(kNCol-col);
	//   errorLedMonAmp=hStripCount[mod][stripGain]->GetBinContent(kNCol-col);
	// }

	//h[mod][gain]->Fill(col, row, (totalChannelAmp/nsame));  
	//weightf=(totalChannelAmp/nsame)/ledMonAmp;
	h[mod][gain]->Fill(col, row, fAvgAmp);  
	weightf=fAvgAmp/ledMonAmp;  

	if (ledMonAmp!=0) hAmpOverMon[mod][gain]->Fill(col, row, weightf);
	if (ledMonAmp==0) hAmpOverMon[mod][gain]->Fill(col, row, 0.0);

	Double_t ampInSM=fAvgAmp;
	Double_t ampError=fRMS;
	// error of the ratio
	Double_t error= TMath::Sqrt(((ampError/ledMonAmp)*(ampError/ledMonAmp))+(((errorLedMonAmp*ampInSM)/(ledMonAmp*ledMonAmp))*((errorLedMonAmp*ampInSM)/(ledMonAmp*ledMonAmp))));

	hCount[mod][gain]->Fill(col, row, fRMS);
      }

      for (int iSM=0; iSM<kNSM; iSM++) {
	nBadStripsInSM[iSM]=0;
      }

      for (int iSM=0; iSM<kNSM; iSM++) {
	int isector = iSM/2;
	int iside = iSM%2; 
	for (int ibinX=1; ibinX<=kNCol; ibinX++) {

	  if(iSM>11 && iSM<18 && ibinX > 32) continue;

	  //	printf("iSM=%d  ibinX=%d    content=%f  ", iSM, ibinX, hStrip[iSM][stripGain]->GetBinContent(ibinX));

	  if (hStrip[iSM][stripGain]->GetBinContent(ibinX)==0.) 
	    {
	      //	 printf("\n");
	      nBadStripsInSM[iSM]=nBadStripsInSM[iSM]+1;
	    }
	  else 
	    //	 printf(" -- not accepted\n");

	    //      if (mod%2==0) {
	    if (hStrip[iSM][stripGain]->GetBinContent(ibinX)!=0) hsmLEDAmpRmsOverMean[iSM]->Fill(ibinX-1, hStripCount[iSM][stripGain]->GetBinContent(ibinX)/hStrip[iSM][stripGain]->GetBinContent(ibinX));
	    else hsmLEDAmpRmsOverMean[iSM]->Fill(ibinX-1, 0);
	  //     }
	  // if (mod%2==1) {
	  // 	if (hStrip[iSM][stripGain]->GetBinContent(kNCol-ibinX+1)!=0) 
	  //        hsmLEDAmpRmsOverMean[iSM]->Fill(ibinX-1, hStripCount[iSM][stripGain]->GetBinContent(kNCol-ibinX+1)/hStrip[iSM][0]->GetBinContent(kNCol-ibinX+1));
	  // 	else hsmLEDAmpRmsOverMean[iSM]->Fill(ibinX-1, 0);
	  // }
	  for (int ibinY=1; ibinY<=kNRow; ibinY++) {
	    if (h[iSM][towerGain]->GetBinContent(ibinX, ibinY)!=0) {
	      hsmAmpRmsOverMean[iSM]->Fill(ibinX-1, ibinY-1, hCount[iSM][towerGain]->GetBinContent(ibinX, ibinY)/h[iSM][towerGain]->GetBinContent(ibinX, ibinY));  
	    }
	    else
	      hsmAmpRmsOverMean[iSM]->Fill(ibinX-1, ibinY-1, 0);
	  }
	}
      }

      for (int kSM=0; kSM<kNSM; kSM++) {
	Printf("nSM %d, number of bad strips = %d", kSM, nBadStripsInSM[kSM]);
	if(nBadStripsInSM[kSM]>0) numberOfSMsWithBadStrips=numberOfSMsWithBadStrips+1;
      }

      Int_t flagForLED=0;
      for (int iSM=0; iSM<kNSM; iSM++) {
	rmsValues[iSM]=hsmAmpRmsOverMean[iSM]->GetBinContent(hsmAmpRmsOverMean[iSM]->GetMaximumBin());
	// if (nBadStripsInSM[iSM]>4) {Printf ("runNumber is = %d, SM = %d", runno, iSM); flagForLED=1;}
      }

      // if(flagForLED==1) continue;  //This is to discard runs where a lot of strips are missing. However if the condition is true for all runs then all runs are discarded (like this year where SM7 has missing LED mon for all runs). Remove the condition. 

      Int_t FEEok=0, allFEEs=1;
      Int_t nFEE[kNSM][36];
      for (int iSM=0; iSM<kNSM; iSM++) {
	if (rmsValues[iSM]<previousRmsValues[iSM]) flagOfChangeRef=flagOfChangeRef+1;
	previousRmsValues[iSM]=rmsValues[iSM];
	//
	for (Int_t i=0; i<36; i++) {nFEE[iSM][i]=0;}

	Int_t iFEE=0;
	Printf("SM = %d", iSM);
	for (Int_t jcol=1; jcol<=kNCol; jcol=jcol+4) {
	  for (Int_t jrow=1; jrow<=kNRow; jrow=jrow+8) {
	    FEEok=0;
	    for (Int_t i=0; i<4; i++) {
	      for (Int_t j=0; j<8; j++) {
		if (h[iSM][1]->GetBinContent(jcol+i, jrow+j)>0) FEEok=1;
	      }
	    }
	    Printf("FEE card number = %d,  status = %d", iFEE, FEEok);
	nFEE[iSM][iFEE]=FEEok;
	iFEE=iFEE+1;
      }
    }
  }

  for (int iSM=0; iSM<kNSM; iSM++) {
    for (Int_t iFEE=0; iFEE<36; iFEE++) {
      if (nFEE[iSM][iFEE]==0) allFEEs=0;
    }
  }
      
  Printf("flag of change ref %d, previous = %d,  allFEEs = %d", flagOfChangeRef, prevflagOfChangeRef, allFEEs);
  Printf("numberOfSMsWithBadStrips = %d, previous = %d", numberOfSMsWithBadStrips, prevnumberOfSMsWithBadStrips);
  
  if ((allFEEs==1)&&(numberOfSMsWithBadStrips==0)) {referenceRun=runno; desiredFlag=1;}
  if (desiredFlag==0) {
    
    //   if (prevflagOfChangeRef<flagOfChangeRef) {referenceRun=runno;}
    // else
    if ((prevnumberOfSMsWithBadStrips>numberOfSMsWithBadStrips)) {referenceRun=runno;} }

  if (flagOfChangeRef>prevflagOfChangeRef) prevflagOfChangeRef=flagOfChangeRef;
  if (numberOfSMsWithBadStrips<prevnumberOfSMsWithBadStrips) prevnumberOfSMsWithBadStrips=numberOfSMsWithBadStrips;
  
  Printf("reference run = %d, desired flag = %d", referenceRun, desiredFlag);
 
  sprintf(filenameout, "hist5/LED_%09dNew.root",runno);

  TFile destFile(filenameout, "recreate");
  destFile.cd();
  

  for (int iSM=0; iSM<kNSM; iSM++) {

    hsmLEDAmpRmsOverMean[iSM]->Write();
    hsmAmpRmsOverMean[iSM]->Write();

    for (int igain=0; igain<2; igain++) {
      h[iSM][igain]->Write();
      hStrip[iSM][igain]->Write();
      hRevStrip[iSM][igain]->Write();
      hCount[iSM][igain]->Write();
      hAmpOverMon[iSM][igain]->Write();
      hStripCount[iSM][igain]->Write();
    }
  }

  destFile.Close();

  for (int iSM=0; iSM<kNSM; iSM++) {
    delete hsmAmpRmsOverMean[iSM];
    delete hsmLEDAmpRmsOverMean[iSM];
    for (int igain=0; igain<2; igain++) {
      delete h[iSM][igain];
      delete hCount[iSM][igain];
      delete hStrip[iSM][igain];
      delete hStripCount[iSM][igain];
      delete hAmpOverMon[iSM][igain];
    }
  }
 
  f->Close();
    }
  }
}

#endif
