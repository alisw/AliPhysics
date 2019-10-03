/*
 author: Ionut-Cristian Arsene
 email: i.c.arsene@cern.ch
 
 Macro to process the TRD QA output for a given run and obtain:
 1. Detailed QA plots
 2. QA trending values
 3. Trending values for calibration extracted from OCDB
 4. Other OCDB parameters (e.g. beam intensity)
 */

AliCDBEntry* GetCDBentry(TString path, Bool_t owner);

//_________________________________________________________________________________
void ProcessTRDRunQA(TString qaFile, Int_t runNumber, TString dataType, 
		     Int_t year, TString period, TString pass, 
		     TString ocdbStorage) {
  //
  // Process run level QA
  // Create standard QA plots and trending tree in the current directory
  const Char_t *friendsOpt="PID DET RES";
  
  //gStyle->SetTitleX(gStyle->GetPadLeftMargin());
	gStyle->SetGridColor(kBlack);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);

  // Load needed libraries
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT -I$ALICE_PHYSICS -I$ALICE_ROOT/TRD");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libSTAT");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGPP");
  
  // Initialize a tree streamer
  TTreeSRedirector *treeStreamer = new TTreeSRedirector("trending.root","RECREATE");
  (*treeStreamer)<< "trending"
          << "run=" << runNumber;
  
  // connect to grid if its the case
  if(qaFile.Contains("alien://") || ocdbStorage.Contains("alien://") || ocdbStorage[0]=='\0')
    TGrid::Connect("alien://");
  
  // trending values from the ESD task ------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TRD/macros/makeResults.C");
  Double_t esdTrendValues[100]; 
  for(Int_t i=0;i<100;i++) esdTrendValues[i]=0.0;
  makeSummaryESD(qaFile.Data(), esdTrendValues, 1);
  const Int_t kNESDtrends = 24;
  const TString kESDTrendNames[kNESDtrends] = {
    "TPCTRDmatchEffPosAll","TPCTRDmatchEffPosAllErr",
    "TPCTRDmatchEffNegAll","TPCTRDmatchEffNegAllErr",					      
    "TRDTOFmatchEffPosAll","TRDTOFmatchEffPosAllErr",					      
    "TRDTOFmatchEffNegAll","TRDTOFmatchEffNegAllErr",						     
    "AvTRDtrkltsPerTrackAll", "AvTRDtrkltsPerTrackAllErr",						      
    "AvNclsPerTrackAll", "AvNclsPerTrackAllErr",						      
    "PHplateauHeight", "PHplateauHeightErr",						
    "PHplateauSlope", "PHplateauSlopeErr",						
    "QtotLandauMPV1GeVAll", "QtotLandauWidth1GeVAll",
    "PHplateauHeightAbsolute", "PHplateauHeightErrAbsolute",
    "PHplateauSlopeAbsolute", "PHplateauSlopeErrAbsolute",						
    "QtotLandauMPV1GeVAllAbsolute", "QtotLandauWidth1GeVAllAbsolute"
  };
  for(Int_t i=0; i<kNESDtrends; ++i)
    (*treeStreamer)<< "trending" << Form("%s=",kESDTrendNames[i].Data()) << esdTrendValues[i];
  
  // process the QA output from the other tasks----------------------------------------
  if(dataType.Contains("sim"))
    makeResults(friendsOpt, qaFile.Data());
  else  
    makeResults(Form("NOMC %s", friendsOpt), qaFile.Data());
  
  TFile *trendFile = TFile::Open("TRD.Trend.root","READ");
  if(!trendFile || !trendFile->IsOpen() ){
    Warning("ProcessTRDRunQA.C", "Couldn't open the TRD.Trend.root file.");
    if(trendFile) delete trendFile; trendFile=0;
  }
  TFile *trendDB = TFile::Open("$ALICE_PHYSICS/PWGPP/TRD/data/TRD.Trend.root","READ");
  if(!trendDB || !trendDB->IsOpen() ){
    Error("ProcessTRDRunQA.C", "Couldn't open the Trending DB !");
    if(trendDB) delete trendDB; trendDB=0;
  }
  
  if(trendDB){
    TKey *tk(NULL); AliTRDtrendValue *tv(NULL); Int_t itv(0);
    const Int_t ntrends(5000);
    Double_t trendValues[ntrends]={0.0};
    TIterator *it(trendDB->GetListOfKeys()->MakeIterator());
    while((tk = (TKey*)it->Next()) && itv < ntrends){
      trendValues[itv] = -999.;
      if(trendFile && (tv = (AliTRDtrendValue*)trendFile->Get(tk->GetName()))){ 
        trendValues[itv] = tv->GetVal();
        printf("Found trend %03d \"%s\" %f\n", itv, tk->GetName(), trendValues[itv]);
      }
      (*treeStreamer)<< "trending" << Form("%s=", tk->GetName()) << trendValues[itv];
      itv++;
    }
  }
  
  // get OCDB information---------------------------------------------------------------
  // switch off grid infos to reduce output and logfilesize                                               
  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  AliCDBManager* man=AliCDBManager::Instance();
  if(ocdbStorage[0]=='\0')
    man->SetDefaultStorage(Form("alien://folder=/alice/data/%d/OCDB/", year));
  else
    man->SetDefaultStorage(ocdbStorage.Data());
  man->SetRun(runNumber);
  
  AliCDBEntry* entryExB = 0x0;
  AliCDBEntry* entryGainFactor = 0x0;
  AliCDBEntry* entryT0 = 0x0;
  AliCDBEntry* entryVdrift = 0x0;
  entryExB = man->Get("TRD/Calib/ChamberExB");
  entryGainFactor = man->Get("TRD/Calib/ChamberGainFactor");
  entryT0 = man->Get("TRD/Calib/ChamberT0");
  entryVdrift = man->Get("TRD/Calib/ChamberVdrift");
  AliTRDCalDet *caldetExB=0x0;
  AliTRDCalDet *caldetGainFactor=0x0;
  AliTRDCalDet *caldetT0=0x0;
  AliTRDCalDet *caldetVdrift=0x0;
  if(entryExB)        caldetExB        = (AliTRDCalDet*)entryExB->GetObject();
  if(entryGainFactor) caldetGainFactor = (AliTRDCalDet*)entryGainFactor->GetObject();
  if(entryT0)         caldetT0         = (AliTRDCalDet*)entryT0->GetObject();
  if(entryVdrift)     caldetVdrift     = (AliTRDCalDet*)entryVdrift->GetObject();
  // get the values
  Double_t meanExB        = (caldetExB ? caldetExB->CalcMean(1) : 0.0);
  Double_t rmsExB         = (caldetExB ? caldetExB->CalcRMS(1) : 0.0);
  Double_t meanGainFactor = (caldetGainFactor ? caldetGainFactor->CalcMean(1) : 0.0);
  Double_t rmsGainFactor  = (caldetGainFactor ? caldetGainFactor->CalcRMS(1) : 0.0);
  Double_t meanT0         = (caldetT0 ? caldetT0->CalcMean(1) : 0.0);
  Double_t rmsT0          = (caldetT0 ? caldetT0->CalcRMS(1) : 0.0);
  Double_t meanVdrift     = (caldetVdrift ? caldetVdrift->CalcMean(1) : 0.0);
  Double_t rmsVdrift      = (caldetVdrift ? caldetVdrift->CalcRMS(1) : 0.0);
  (*treeStreamer)<< "trending"
                 << "meanExB=" << meanExB
                 << "rmsExB=" << rmsExB
                 << "meanGainFactor=" << meanGainFactor
                 << "rmsGainFactor=" << rmsGainFactor
                 << "meanT0=" << meanT0
                 << "rmsT0=" << rmsT0
                 << "meanVdrift=" << meanVdrift
                 << "rmsVdrift=" << rmsVdrift;
  
  // Get the beam luminosity
  AliCDBEntry *entryLHCData = man->Get("GRP/GRP/LHCData");		 
  AliLHCData *lhcData = (entryLHCData ? (AliLHCData*)entryLHCData->GetObject() : 0x0);
  Double_t beamIntensityA=0.0;
  Double_t beamIntensityC=0.0;
  if(lhcData) {
    Int_t nLumiMeasA=lhcData->GetNLuminosityTotal(0); Int_t nA=0;
    Int_t nLumiMeasC=lhcData->GetNLuminosityTotal(1); Int_t nC=0;
    // Sum up the measurements
    AliLHCDipValF *dipVal0,*dipVal1;
    for(Int_t iLumiMeas=0;iLumiMeas<nLumiMeasA;iLumiMeas++){
      dipVal0 = lhcData->GetLuminosityTotal(0,iLumiMeas);
      if(dipVal0) {
        beamIntensityA += dipVal0->GetValue();
        ++nA;
      }
    }
    if(nA) beamIntensityA /= Double_t(nA);
    for(Int_t iLumiMeas=0;iLumiMeas<nLumiMeasC;iLumiMeas++){
      dipVal1 = lhcData->GetLuminosityTotal(1,iLumiMeas);
      if(dipVal1) {
        beamIntensityC += dipVal1->GetValue();
        ++nC;
      }
    }
    if(nC) beamIntensityC /= Double_t(nC);
  }
  (*treeStreamer)<< "trending"
                 << "beamIntensityA=" << beamIntensityA
                 << "beamIntensityC=" << beamIntensityC;
                 
  // Get the magnetic field polarity
  Double_t Bfield=-2.;
  AliGRPManager grpManager;
  if(grpManager.ReadGRPEntry() && grpManager.SetMagField()){
    AliMagF *f=TGeoGlobalMagField::Instance()->GetField();
    Bfield=f->Factor();
  }                  
  (*treeStreamer)<< "trending"
                 << "Bfield=" << Bfield;
		 
  (*treeStreamer)<< "trending"
                 << "\n";
  delete treeStreamer;		 
}
