#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TH2I.h>
  #include <TH1D.h>
  #include <TH2D.h>
  #include <TH1I.h>
  #include <TCanvas.h>
  #include <TMatrixD.h>
  #include <TSystem.h>
  #include <TROOT.h>
  #include <TClonesArray.h>
  #include <TTree.h>
  #include <TFile.h>  
  #include <TDirectoryFile.h>

  #include <AliRun.h>
  #include <AliStack.h>
  #include <AliLoader.h>
  #include <AliGeomManager.h>
  #include <AliITSUGeomTGeo.h>
  #include <AliITSMFTDigitPix.h>
  #include <AliITSMFTSegmentationPix.h>
 
  #include "AliITSUSuze02.h"

#endif

//v4 - SuZe limits version can be chosen. Requires MakeSuze_v2.h
//v5 - saving average number of digits per encoding window for each module/event
//v6 - saving the lowest remaining number of windows for each of three limits for eache module/event. Requires MakeSuze_v3.h
//v7 - saving distribution of number of digits per encoding window per layer and in total. Requires at least MakeSuze_v4.h
//v8 - changed the number of rows for the small and big sensor to 320 and 640 rescpectively. Used for the "final" production.
//v9 - compatibility changes to MakeSuze_v7 that returns the data size as well.       
//v10 - modules with signal hits are processed by SUZE
//v11 - version to use with MakeSuze_v8 where data size calculation has been corrected and SuzeLimitsVersion=99 has been added
//v12 - switch added to process modules with only signal hits        
//v13 - AliITSUSuze02.h is used instead of MakeSuze_v8.h
//v14 - QED digits are added from the separate file  
//v15 - switch added to add QED digits

Int_t GetSuzeLimits(Int_t DataSizePerWindowInRow, Int_t Version, Int_t& Limit32, Int_t& LimitHalfFSBB, Int_t& LimitFSBB){
  
  switch (Version){
    case 1:
      Limit32=6;
      LimitHalfFSBB=12;
      LimitFSBB=19;
      return 1;
    case 2:
      Limit32=6;
      LimitHalfFSBB=9;
      LimitFSBB=18;
      return 1;
    case 3:
      Limit32=5;
      LimitHalfFSBB=6;
      LimitFSBB=7;
      return 1;
    case 4:
      Limit32=4;
      LimitHalfFSBB=5;
      LimitFSBB=6;
      return 1;
    case 99:
      Limit32=(Int_t)180/DataSizePerWindowInRow;
      LimitHalfFSBB=(Int_t)360/DataSizePerWindowInRow;
      LimitFSBB=(Int_t)570/DataSizePerWindowInRow;
      return 1;
    default:
      return 0;
  }
}

void ScanDigitsSuze02_v15(Int_t Cycle=0, Int_t CollectMode=0, Bool_t ProcessOnlyChipsWithSignal=0, Bool_t AddQED=0, Int_t nEvents=-1, Int_t NRowsEncodingWindow=4,  Int_t NColsEncodingWindow=5, Int_t SuzeLimitsVersion=99, Bool_t SaveResults=kTRUE){ 
//CollectMode - defines which digits are added to SUZE matrix (-1 - Noise, +1 - Signal, 0 - Noise+Signal)
//ProcessOnlyChipsWithSignal - defines if only modules with signal hits are processed, if 0 all the modules are processed
 
  Int_t DataSizePerWindowInRow = 8 + NRowsEncodingWindow*NColsEncodingWindow+TMath::Ceil(TMath::Log2(NRowsEncodingWindow));

  gROOT->SetStyle("Plain");
  if(nEvents==-1) cout<<"Analysing all events in the run";
  else cout<<"Analysing first "<<nEvents<<" events of the run";
  cout<<" (cycle #"<<Cycle<<")";
  if(CollectMode==-1) cout<<" taking only noise pixels";
  else if(CollectMode==1) cout<<" taking only signal pixels";
  cout<<endl;    
  if(ProcessOnlyChipsWithSignal) cout<<"Only modules with signal hits will be processed";
  else cout<<"All modules will be processed";
  cout<<endl;
  //Int_t debugOn=0; //DEBUG
  Int_t nWindowsPerFSBBMax=19;
  Int_t nWindowsPerHalfFSBBMax=12;
  Int_t nWindowsPer32colsMax=6;

  if(!GetSuzeLimits(DataSizePerWindowInRow,SuzeLimitsVersion,nWindowsPer32colsMax,nWindowsPerHalfFSBBMax,nWindowsPerFSBBMax)){
   cout<<"Wrong suze limits version!"<<endl;
   return;
  }
  else{
    cout<<"Suze will use the follwing limits:"<<endl;
    cout<<nWindowsPer32colsMax<<" windows per 32 columns,"<<endl;
    cout<<nWindowsPerHalfFSBBMax<< " windows per half FSBB,"<<endl;
    cout<<nWindowsPerFSBBMax<<" windows per FSBB"<<endl;
  }

  Int_t nWindowsPerFSBBMin=nWindowsPerFSBBMax;
  Int_t nWindowsPerHalfFSBBMin=nWindowsPerHalfFSBBMax;
  Int_t nWindowsPer32colsMin=nWindowsPer32colsMax;

  Char_t logfile_name[100];
  sprintf(logfile_name,"ScanDigits_v15_log_Cycle_%d_nEvents_%d_EncWindow_%dx%d_SuzeLimitsVersion_%d_Mode_%d-%d_QED_%d.log",Cycle,nEvents,NRowsEncodingWindow,NColsEncodingWindow,SuzeLimitsVersion,CollectMode,ProcessOnlyChipsWithSignal,AddQED);
  FILE *logfile = fopen (logfile_name,"w");

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadSDigits();
  runLoader->LoadDigits();

  AliGeomManager::LoadGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE,kTRUE);
  //
  Int_t nLayers = gm->GetNLayers();
  Int_t nChips = gm->GetNChips();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  //DIGITS INIT
  TTree * digTree = 0x0; 
  TClonesArray *digArr = new TClonesArray("AliITSMFTDigitPix"); 
  
  TTree * digTreeQED = 0x0;
  TClonesArray *digArrQED = new TClonesArray("AliITSMFTDigitPix");
  TDirectoryFile* EventDirectoryQED = 0x0;
  TFile* DigitsFileQED = 0x0;
  
  if(AddQED) {
      
    DigitsFileQED = new TFile("QED/ITS.Digits.root");
    
  }
  if(nEvents==-1) nEvents=runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",nEvents);

  //Chip sizes
  Int_t Chip_Ncols=1362;
  Int_t Chip_Nrows_small=320;
  Int_t Chip_Nrows_big=640;
  Int_t ColAddress=0;
  Int_t RowAddress=0;

  Int_t DataSize=0;
  //Int_t ChipSum=0;

  TH1F* OverflowCodesPerChip = new TH1F("OverflowCodesPerChip","OverflowCodesPerChip",8,0,8);
  TH1F* OverflowCodes = new TH1F("OverflowCodes","Overflow codes",8,0,8);
  TH1F* OverflowCodesPerLayer[nLayers];

  TH1F* nDigitsPerEncodingWindowPerChip = new TH1F("nDigitsPerEncodingWindowPerChip","nDigitsPerEncodingWindowPerChip",NRowsEncodingWindow*NColsEncodingWindow,1,NRowsEncodingWindow*NColsEncodingWindow+1);
  TH1F* nDigitsPerEncodingWindow = new TH1F("nDigitsPerEncodingWindow","nDigitsPerEncodingWindow",NRowsEncodingWindow*NColsEncodingWindow,1,NRowsEncodingWindow*NColsEncodingWindow+1);
  TH1F* nDigitsPerEncodingWindowPerLayer[nLayers];

  Int_t nDigitsLostPerChip=0;
  Int_t nDigitsEncodedPerChip=0;
  Int_t nDigitsLostPerEvent=0;
  Int_t nDigitsPerEvent=0;
  Double_t FractionDigitsLostPerEvent=0;

  Int_t nDigitsLostPerEventPerLayer[nLayers];
  Int_t nDigitsPerEventPerLayer[nLayers];
  Double_t FractionDigitsLostPerEventPerLayer[nLayers];
  Int_t MaxNWindowsPerStavePerLayerPerEvent[nLayers];

  TH1I* NtracksPerEvent_hist = new TH1I("NtracksPerEvent","Ntracks per event",nEvents,0,nEvents);
  TH1I* NdigitsPerEvent_hist = new TH1I("NdigitsPerEvent","Ndigits per event",nEvents,0,nEvents);
  TH1I* NdigitsPerEventPerLayer_hist[nLayers];

  TH1I* nDigitsLostPerEvent_hist = new TH1I("nDigitsLostPerEvent","Digits lost per event",nEvents,0,nEvents);

  TH1D* FractionDigitsLostPerEvent_hist = new TH1D("FractionDigitsLostPerEvent","Fraction of Digits lost per event",nEvents,0,nEvents);
  TH1D* FractionDigitsLostPerEventPerLayer_hist[nLayers];

  TH1I* MaxNWindowsPerStavePerLayerPerEvent_hist[nLayers];

  Int_t nWindows=0;

  //complete maps
  TH2I* nDigitsPerChipPerEvent = new TH2I("nDigitsPerChipPerEvent","nDigits per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2I* nDigitsEncodedPerChipPerEvent = new TH2I("nDigitsEncodedPerChipPerEvent","nDigits encoded per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2I* nDigitsLostPerChipPerEvent = new TH2I("nDigitsLostPerChipPerEvent","nDigits lost per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2F* FractionDigitsLostPerChipPerEvent = new TH2F("FractionDigitsLostPerChipPerEvent","Fraction of digits lost per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2I* nEncodingWindowsPerChipPerEvent = new TH2I("nEncodingWindowsPerChipPerEvent","Encoding windows per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2F* nDigitsPerEncodingWindowPerChipPerEvent = new TH2F("nDigitsPerEncodingWindowPerChipPerEvent","Average digits per encoding window per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);

  TH2I* nWindowsPerFSBBMinPerChipPerEvent = new TH2I("nWindowsPerFSBBMinPerChipPerEvent","nWindowsPerFSBBMin per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2I* nWindowsPerHalfFSBBMinPerChipPerEvent = new TH2I("nWindowsPerHalfFSBBMinPerChipPerEvent","nWindowsPerHalfFSBBMin per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);
  TH2I* nWindowsPer32colsMinPerChipPerEvent = new TH2I("nWindowsPer32colsMinPerChipPerEvent","nWindowsPer32colsMin per Chip per Event",nChips,0,nChips,nEvents,0,nEvents);

  TH2F* DataSizePerChipPerEvent = new TH2F("DataSizePerChipPerEvent","DataSizePerChipPerEvent",nChips,0,nChips,nEvents,0,nEvents);

  Int_t current_stave=-1;
  Int_t current_layer=-1;
  Double_t NWindowsPerStave=0;

  for(Int_t i=0; i<nLayers; i++){
    nDigitsLostPerEventPerLayer[i]=0;
    nDigitsPerEventPerLayer[i]=0;
    FractionDigitsLostPerEventPerLayer[i]=0;

    NdigitsPerEventPerLayer_hist[i] = new TH1I(Form("NdigitsPerEventPerLayer_%d",i),Form("Ndigits at layer %d",i),nEvents,0,nEvents);
    FractionDigitsLostPerEventPerLayer_hist[i] = new TH1D(Form("FractionDigitsLostPerEventPerLayer_%d",i),Form("Fraction of digits lost per event at layer %d",i),nEvents,0,nEvents);

    MaxNWindowsPerStavePerLayerPerEvent[i]=0;
    MaxNWindowsPerStavePerLayerPerEvent_hist[i] = new TH1I(Form("MaxNWindowsPerStavePerLayerPerEvent_%d",i),Form("Max number of windows per stave per event at layer %d",i),nEvents,0,nEvents);

    OverflowCodesPerLayer[i] = new TH1F(Form("OverflowCodesPerLayer_%d",i),Form("Overflow codes at layer %d",i),8,0,8);
    nDigitsPerEncodingWindowPerLayer[i] = new TH1F(Form("nDigitsPerEncodingWindowPerLayer_%d",i),Form("nDigitsPerEncodingWindowPerLayer_%d",i),NRowsEncodingWindow*NColsEncodingWindow,1,NRowsEncodingWindow*NColsEncodingWindow+1);
  }

  Int_t NtracksKine=0;
  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
  //for (Int_t iEvent = 7; iEvent < nEvents; iEvent++) {    //!!!DEBUG
    printf("\n Event %i \n",iEvent);
    runLoader->GetEvent(iEvent);
    AliStack *stack = runLoader->Stack();
    NtracksKine=stack->TreeK()->GetEntries();
    digTree=dl->TreeD(); 
    
    //manually load digits tree for QED
    if(AddQED){
      EventDirectoryQED=(TDirectoryFile*)DigitsFileQED->Get(Form("Event%d",iEvent));          
      digTreeQED=(TTree*)EventDirectoryQED->Get("TreeD"); 
      digTreeQED->SetBranchAddress("ITSDigitsPix",&digArrQED);
    }
    //
    digTree->SetBranchAddress("ITSDigitsPix",&digArr);
    nDigitsLostPerEvent=0;
    nDigitsPerEvent=0;
    for(Int_t i=0; i<nLayers; i++){
      nDigitsPerEventPerLayer[i]=0;
      nDigitsLostPerEventPerLayer[i]=0;
      MaxNWindowsPerStavePerLayerPerEvent[i]=0;
    }
    current_stave=-1;
    current_layer=-1;
    NWindowsPerStave=0;
    
    Int_t ndigQED=0;
    
    for (Int_t imod=0;imod<nChips;imod++) {
      AliITSUSuze02* Chip;
      nDigitsLostPerChip=0;
      nDigitsEncodedPerChip=0;
      digTree->GetEntry(imod); 
      
      if(AddQED){
        digTreeQED->GetEntry(imod);
        ndigQED = digArrQED->GetEntries();
      }
      //Int_t detType = gm->GetChipChipTypeID(imod);
      //AliITSMFTSegmentationPix* segm = (AliITSMFTSegmentationPix*)gm->GetSegmentationByID(detType);
      Int_t lay,sta,det;
      Int_t ndig  = digArr->GetEntries();
      
      Int_t ndig_in_cycle=0; 
      Int_t ndig_signal_in_cycle=0;
      if (ndig<1) continue;
      gm->GetChipId(imod, lay,sta,det);
//       printf("\nChip %3d: (det %2d in stave %2d of Layer %d) | NDigits: %4d\n",imod,det,sta,lay,ndig);
      //    
      //if(iEvent==7 && lay==2) cout<<"Chip #"<<imod<<endl;     //!!!DEBUG   
      //if(iEvent==7 && lay==2 && imod==363) debugOn=1; //continue;  //!!!DEBUG
      //else debugOn=0; 
      if(lay>=0 && lay <=2){
        Chip = new AliITSUSuze02(Chip_Nrows_small,Chip_Ncols);
      }
      else{
        Chip = new AliITSUSuze02(Chip_Nrows_big,Chip_Ncols);
      }
      
      Chip->SetEncodingWindowSize(NRowsEncodingWindow,NColsEncodingWindow);
      Chip->SetQuotas(nWindowsPer32colsMax,nWindowsPerHalfFSBBMax,nWindowsPerFSBBMax);
      
      OverflowCodesPerChip->Reset();
      nDigitsPerEncodingWindowPerChip->Reset();

      nWindows=0;
      for(Int_t idig=0;idig<ndig;idig++) {
        AliITSMFTDigitPix *pDig = (AliITSMFTDigitPix*)digArr->At(idig);
        if(pDig->GetROCycle()!=Cycle) continue; //selection of the hits from a given RO cycle
        if((CollectMode==-1 && pDig->GetHit(0)!=-1) || (CollectMode==1 && pDig->GetHit(0)==-1)) continue; // selection between noise and signal or both
        ndig_in_cycle++; 
        if(pDig->GetHit(0)!=-1) ndig_signal_in_cycle++;  //counts signal hits in a given RO cycle
        ColAddress=pDig->GetCoord1();
        RowAddress=pDig->GetCoord2();
        Chip->AddDigit(RowAddress,ColAddress);
      }//diglist  
      
      if(AddQED){
        for(Int_t idig=0;idig<ndigQED;idig++) {
          AliITSMFTDigitPix *pDig = (AliITSMFTDigitPix*)digArrQED->At(idig);
          if(pDig->GetROCycle()!=Cycle) continue; //selection of the hits from a given RO cycle
          ndig_in_cycle++; 
          ColAddress=pDig->GetCoord1();
          RowAddress=pDig->GetCoord2();
          Chip->AddDigit(RowAddress,ColAddress);
        }//diglist
      }
      
      if(ndig_in_cycle<1){  
        delete Chip; 
        continue; //rejects when no hits in a given cycle
      }
      nDigitsPerEvent+=ndig_in_cycle;
      nDigitsPerEventPerLayer[lay]+=ndig_in_cycle;
      if(ProcessOnlyChipsWithSignal){       //if ProcessOnlyChipsWithSignal==1
        if(ndig_signal_in_cycle<1){   
          delete Chip;
          continue;  //rejects when only noise is the present
        }
      }
      
      Chip->Process(OverflowCodesPerChip,nDigitsPerEncodingWindowPerChip);
      DataSize=Chip->GetDataSize();
      nDigitsEncodedPerChip=Chip->GetNDigitsEncoded();
      nDigitsLostPerChip=Chip->GetNDigitsLost();
      nWindows=Chip->GetNEncodedWindows();                           
      nWindowsPer32colsMin=Chip->GetNWindowsPer32colsMin();
      nWindowsPerHalfFSBBMin=Chip->GetNWindowsPerHalfFSBBMin();
      nWindowsPerFSBBMin=Chip->GetNWindowsPerFSBBMin();
      
     //DataSize=MakeSuze(Chip_matrix,OverflowCodesPerChip,nDigitsEncodedPerChip,nDigitsLostPerChip,nWindows,nWindowsPer32colsMax,nWindowsPerHalfFSBBMax,nWindowsPerFSBBMax,nWindowsPer32colsMin,nWindowsPerHalfFSBBMin,nWindowsPerFSBBMin,nDigitsPerEncodingWindowPerChip);
//       cout<<"SUZE encoded "<<SuzeReturn<<" digits in "<<nWindows<<" windows"<<endl;
      OverflowCodes->Add(OverflowCodesPerChip);
      OverflowCodesPerLayer[lay]->Add(OverflowCodesPerChip);

      nDigitsPerEncodingWindow->Add(nDigitsPerEncodingWindowPerChip);
      nDigitsPerEncodingWindowPerLayer[lay]->Add(nDigitsPerEncodingWindowPerChip);

//       if(nDigitsLostPerChip){
//         cout<<"------- Some digits were lost. Check overflow errors"<<endl;
//         cout<<"Chip has "<<ChipSum<<" digits"<<endl;
//         cout<<"SUZE reported "<<SuzeReturn<<" encoded and "<<nDigitsLost<<" lost digits"<<endl;
//         cout<<"Chip n."<<imod<<" has lost "<<nDigitsLostPerChip<<" digits due to the overflow"<<endl;

//       }
      if(nDigitsLostPerChip){
	      printf("Event #%d mod. #%d (layer %d) has %d digits lost\n",iEvent,imod,lay,nDigitsLostPerChip);
	      fprintf(logfile,"Event #%d mod. #%d (layer %d) has %d (%f) digits lost\n",iEvent,imod,lay,nDigitsLostPerChip, (Float_t)nDigitsLostPerChip/ndig_in_cycle);
      }
      nDigitsLostPerEvent+=nDigitsLostPerChip;
      nDigitsLostPerEventPerLayer[lay]+=nDigitsLostPerChip;
//       cout<<"Lay:"<<lay<<" Sta:"<<sta<<" current_stave:"<<current_stave<<" current layer:"<<current_layer<<" Digits:"<<SuzeReturn<<endl;
      if(lay!=current_layer){
        cout<<"Layer #"<<lay<<endl;
      }
      if(sta!=current_stave || lay!=current_layer){
        current_stave=sta;
	      current_layer=lay;
	      NWindowsPerStave=0;
      }
      NWindowsPerStave+=nWindows;
      if(NWindowsPerStave>MaxNWindowsPerStavePerLayerPerEvent[current_layer]){
	      MaxNWindowsPerStavePerLayerPerEvent[current_layer]=NWindowsPerStave;
// 	cout<<"---- MaxNWindowsPerStavePerLayerPerEvent:"<<MaxNWindowsPerStavePerLayerPerEvent[current_layer]<<" at layer:"<<current_layer<<endl;
      }
      nDigitsPerChipPerEvent->SetBinContent(imod+1,iEvent+1,ndig_in_cycle);
      nDigitsEncodedPerChipPerEvent->SetBinContent(imod+1,iEvent+1,nDigitsEncodedPerChip);
      nDigitsLostPerChipPerEvent->SetBinContent(imod+1,iEvent+1,nDigitsLostPerChip);
      FractionDigitsLostPerChipPerEvent->SetBinContent(imod+1,iEvent+1,(Double_t)nDigitsLostPerChip/ndig_in_cycle);
      nEncodingWindowsPerChipPerEvent->SetBinContent(imod+1,iEvent+1,nWindows);
      nDigitsPerEncodingWindowPerChipPerEvent->SetBinContent(imod+1,iEvent+1,(Double_t)nDigitsEncodedPerChip/nWindows); 

      nWindowsPerFSBBMinPerChipPerEvent->SetBinContent(imod+1,iEvent+1,nWindowsPerFSBBMax-nWindowsPerFSBBMin);
      nWindowsPerHalfFSBBMinPerChipPerEvent->SetBinContent(imod+1,iEvent+1,nWindowsPerHalfFSBBMax-nWindowsPerHalfFSBBMin);
      nWindowsPer32colsMinPerChipPerEvent->SetBinContent(imod+1,iEvent+1,nWindowsPer32colsMax-nWindowsPer32colsMin);

      DataSizePerChipPerEvent->SetBinContent(imod+1,iEvent+1,DataSize);   
      
      delete Chip;
//      break;
    }//mod
    NtracksPerEvent_hist->SetBinContent(iEvent+1,NtracksKine);
    NdigitsPerEvent_hist->SetBinContent(iEvent+1,nDigitsPerEvent);
    for(Int_t i=0; i<nLayers; i++){
      NdigitsPerEventPerLayer_hist[i]->SetBinContent(iEvent+1,nDigitsPerEventPerLayer[i]);
    }

    nDigitsLostPerEvent_hist->SetBinContent(iEvent+1,nDigitsLostPerEvent);

    if(nDigitsPerEvent){
      FractionDigitsLostPerEvent=(Double_t)nDigitsLostPerEvent/nDigitsPerEvent;
      if(FractionDigitsLostPerEvent) cout<<"Fraction lost = "<<FractionDigitsLostPerEvent<<endl;
      FractionDigitsLostPerEvent_hist->SetBinContent(iEvent+1,FractionDigitsLostPerEvent);
      for(Int_t i=0; i<nLayers; i++){
        if(nDigitsPerEventPerLayer[i]){
          FractionDigitsLostPerEventPerLayer[i]=(Double_t)nDigitsLostPerEventPerLayer[i]/nDigitsPerEventPerLayer[i];
          if(FractionDigitsLostPerEventPerLayer[i]) cout<<"Fraction lost at layer:"<<i<<" is "<<FractionDigitsLostPerEventPerLayer[i]<<endl;
          FractionDigitsLostPerEventPerLayer_hist[i]->SetBinContent(iEvent+1,FractionDigitsLostPerEventPerLayer[i]);
        }
      }
    }
    for(Int_t i=0; i<nLayers; i++){
      MaxNWindowsPerStavePerLayerPerEvent_hist[i]->SetBinContent(iEvent+1,MaxNWindowsPerStavePerLayerPerEvent[i]);
    }
//     break;
  }//event loop
//   TCanvas* c1=new TCanvas();
//   c1->SetLogy();
//   OverflowCodes->DrawNormalized();
//   new TCanvas();
//   nMultipleEncodingsPerEvent_hist->Draw();
//   new TCanvas();
//   nDigitsLostPerEvent_hist->Draw();

  TFile* ResultsFile;
  Char_t ResultsFileName[100];
  if(SaveResults){
    sprintf(ResultsFileName,"ScanDigits_v15_results_cycle_%d_EncWindow_%dx%d_SuzeLimitsVersion_%d_mode_%d-%d_QED_%d.root",Cycle,NRowsEncodingWindow,NColsEncodingWindow,SuzeLimitsVersion,CollectMode,ProcessOnlyChipsWithSignal,AddQED);
    ResultsFile = new TFile(ResultsFileName,"RECREATE");

    OverflowCodes->Write();
    nDigitsPerEncodingWindow->Write();

    NtracksPerEvent_hist->Write();
    NdigitsPerEvent_hist->Write();

    nDigitsLostPerEvent_hist->Write();

    FractionDigitsLostPerEvent_hist->Write();

    for(Int_t i=0; i<nLayers; i++){
      NdigitsPerEventPerLayer_hist[i]->Write();
      FractionDigitsLostPerEventPerLayer_hist[i]->Write();

      MaxNWindowsPerStavePerLayerPerEvent_hist[i]->Write();

      OverflowCodesPerLayer[i]->Write();
      nDigitsPerEncodingWindowPerLayer[i]->Write();
    }

    nDigitsPerChipPerEvent->Write();
    nDigitsEncodedPerChipPerEvent->Write();
    nDigitsLostPerChipPerEvent->Write();
    FractionDigitsLostPerChipPerEvent->Write();
    nEncodingWindowsPerChipPerEvent->Write();
    nDigitsPerEncodingWindowPerChipPerEvent->Write();

    nWindowsPerFSBBMinPerChipPerEvent->Write();
    nWindowsPerHalfFSBBMinPerChipPerEvent->Write();
    nWindowsPer32colsMinPerChipPerEvent->Write();

    DataSizePerChipPerEvent->Write();

    ResultsFile->Close();
  }
  fclose (logfile);
//   cout<<"Multiple Encodings:"<<nMultipleEncodings<<endl;
//   cout<<"Lost digits:"<<nDigitsLostTotal<<endl;

}
