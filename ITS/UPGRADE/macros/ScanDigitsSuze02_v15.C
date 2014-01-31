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
  #include <AliITSUDigitPix.h>
  #include <AliITSUSegmentationPix.h>
 
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

void ScanDigitsSuze02_v15(Int_t Cycle=0, Int_t CollectMode=0, Bool_t ProcessOnlyModulesWithSignal=0, Bool_t AddQED=0, Int_t nEvents=-1, Int_t NRowsEncodingWindow=4,  Int_t NColsEncodingWindow=5, Int_t SuzeLimitsVersion=99, Bool_t SaveResults=kTRUE){ 
//CollectMode - defines which digits are added to SUZE matrix (-1 - Noise, +1 - Signal, 0 - Noise+Signal)
//ProcessOnlyModulesWithSignal - defines if only modules with signal hits are processed, if 0 all the modules are processed
 
  Int_t DataSizePerWindowInRow = 8 + NRowsEncodingWindow*NColsEncodingWindow+TMath::Ceil(TMath::Log2(NRowsEncodingWindow));

  gROOT->SetStyle("Plain");
  if(nEvents==-1) cout<<"Analysing all events in the run";
  else cout<<"Analysing first "<<nEvents<<" events of the run";
  cout<<" (cycle #"<<Cycle<<")";
  if(CollectMode==-1) cout<<" taking only noise pixels";
  else if(CollectMode==1) cout<<" taking only signal pixels";
  cout<<endl;    
  if(ProcessOnlyModulesWithSignal) cout<<"Only modules with signal hits will be processed";
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
  sprintf(logfile_name,"ScanDigits_v15_log_Cycle_%d_nEvents_%d_EncWindow_%dx%d_SuzeLimitsVersion_%d_Mode_%d-%d_QED_%d.log",Cycle,nEvents,NRowsEncodingWindow,NColsEncodingWindow,SuzeLimitsVersion,CollectMode,ProcessOnlyModulesWithSignal,AddQED);
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
  Int_t nModules = gm->GetNModules();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  //DIGITS INIT
  TTree * digTree = 0x0; 
  TClonesArray *digArr = new TClonesArray("AliITSUDigitPix"); 
  
  TTree * digTreeQED = 0x0;
  TClonesArray *digArrQED = new TClonesArray("AliITSUDigitPix");
  TDirectoryFile* EventDirectoryQED = 0x0;
  TFile* DigitsFileQED = 0x0;
  
  if(AddQED) {
      
    DigitsFileQED = new TFile("QED/ITS.Digits.root");
    
  }
  if(nEvents==-1) nEvents=runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",nEvents);

  //Module sizes
  Int_t Module_Ncols=1362;
  Int_t Module_Nrows_small=320;
  Int_t Module_Nrows_big=640;
  Int_t ColAddress=0;
  Int_t RowAddress=0;

  Int_t DataSize=0;
  //Int_t ModuleSum=0;

  TH1F* OverflowCodesPerModule = new TH1F("OverflowCodesPerModule","OverflowCodesPerModule",8,0,8);
  TH1F* OverflowCodes = new TH1F("OverflowCodes","Overflow codes",8,0,8);
  TH1F* OverflowCodesPerLayer[nLayers];

  TH1F* nDigitsPerEncodingWindowPerModule = new TH1F("nDigitsPerEncodingWindowPerModule","nDigitsPerEncodingWindowPerModule",NRowsEncodingWindow*NColsEncodingWindow,1,NRowsEncodingWindow*NColsEncodingWindow+1);
  TH1F* nDigitsPerEncodingWindow = new TH1F("nDigitsPerEncodingWindow","nDigitsPerEncodingWindow",NRowsEncodingWindow*NColsEncodingWindow,1,NRowsEncodingWindow*NColsEncodingWindow+1);
  TH1F* nDigitsPerEncodingWindowPerLayer[nLayers];

  Int_t nDigitsLostPerModule=0;
  Int_t nDigitsEncodedPerModule=0;
  Int_t nDigitsLostPerEvent=0;
  Int_t nDigitsPerEvent=0;
  Double_t FractionDigitsLostPerEvent=0;

  Int_t nDigitsLostPerEventPerLayer[nLayers];
  Int_t nDigitsPerEventPerLayer[nLayers];
  Double_t FractionDigitsLostPerEventPerLayer[nLayers];
  Int_t MaxNWindowsPerLadderPerLayerPerEvent[nLayers];

  TH1I* NtracksPerEvent_hist = new TH1I("NtracksPerEvent","Ntracks per event",nEvents,0,nEvents);
  TH1I* NdigitsPerEvent_hist = new TH1I("NdigitsPerEvent","Ndigits per event",nEvents,0,nEvents);
  TH1I* NdigitsPerEventPerLayer_hist[nLayers];

  TH1I* nDigitsLostPerEvent_hist = new TH1I("nDigitsLostPerEvent","Digits lost per event",nEvents,0,nEvents);

  TH1D* FractionDigitsLostPerEvent_hist = new TH1D("FractionDigitsLostPerEvent","Fraction of Digits lost per event",nEvents,0,nEvents);
  TH1D* FractionDigitsLostPerEventPerLayer_hist[nLayers];

  TH1I* MaxNWindowsPerLadderPerLayerPerEvent_hist[nLayers];

  Int_t nWindows=0;

  //complete maps
  TH2I* nDigitsPerModulePerEvent = new TH2I("nDigitsPerModulePerEvent","nDigits per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2I* nDigitsEncodedPerModulePerEvent = new TH2I("nDigitsEncodedPerModulePerEvent","nDigits encoded per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2I* nDigitsLostPerModulePerEvent = new TH2I("nDigitsLostPerModulePerEvent","nDigits lost per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2F* FractionDigitsLostPerModulePerEvent = new TH2F("FractionDigitsLostPerModulePerEvent","Fraction of digits lost per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2I* nEncodingWindowsPerModulePerEvent = new TH2I("nEncodingWindowsPerModulePerEvent","Encoding windows per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2F* nDigitsPerEncodingWindowPerModulePerEvent = new TH2F("nDigitsPerEncodingWindowPerModulePerEvent","Average digits per encoding window per Module per Event",nModules,0,nModules,nEvents,0,nEvents);

  TH2I* nWindowsPerFSBBMinPerModulePerEvent = new TH2I("nWindowsPerFSBBMinPerModulePerEvent","nWindowsPerFSBBMin per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2I* nWindowsPerHalfFSBBMinPerModulePerEvent = new TH2I("nWindowsPerHalfFSBBMinPerModulePerEvent","nWindowsPerHalfFSBBMin per Module per Event",nModules,0,nModules,nEvents,0,nEvents);
  TH2I* nWindowsPer32colsMinPerModulePerEvent = new TH2I("nWindowsPer32colsMinPerModulePerEvent","nWindowsPer32colsMin per Module per Event",nModules,0,nModules,nEvents,0,nEvents);

  TH2F* DataSizePerModulePerEvent = new TH2F("DataSizePerModulePerEvent","DataSizePerModulePerEvent",nModules,0,nModules,nEvents,0,nEvents);

  Int_t current_ladder=-1;
  Int_t current_layer=-1;
  Double_t NWindowsPerLadder=0;

  for(Int_t i=0; i<nLayers; i++){
    nDigitsLostPerEventPerLayer[i]=0;
    nDigitsPerEventPerLayer[i]=0;
    FractionDigitsLostPerEventPerLayer[i]=0;

    NdigitsPerEventPerLayer_hist[i] = new TH1I(Form("NdigitsPerEventPerLayer_%d",i),Form("Ndigits at layer %d",i),nEvents,0,nEvents);
    FractionDigitsLostPerEventPerLayer_hist[i] = new TH1D(Form("FractionDigitsLostPerEventPerLayer_%d",i),Form("Fraction of digits lost per event at layer %d",i),nEvents,0,nEvents);

    MaxNWindowsPerLadderPerLayerPerEvent[i]=0;
    MaxNWindowsPerLadderPerLayerPerEvent_hist[i] = new TH1I(Form("MaxNWindowsPerLadderPerLayerPerEvent_%d",i),Form("Max number of windows per ladder per event at layer %d",i),nEvents,0,nEvents);

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
      MaxNWindowsPerLadderPerLayerPerEvent[i]=0;
    }
    current_ladder=-1;
    current_layer=-1;
    NWindowsPerLadder=0;
    
    Int_t ndigQED=0;
    
    for (Int_t imod=0;imod<nModules;imod++) {
      AliITSUSuze02* Module;
      nDigitsLostPerModule=0;
      nDigitsEncodedPerModule=0;
      digTree->GetEntry(imod); 
      
      if(AddQED){
        digTreeQED->GetEntry(imod);
        ndigQED = digArrQED->GetEntries();
      }
      //Int_t detType = gm->GetModuleDetTypeID(imod);
      //AliITSUSegmentationPix* segm = (AliITSUSegmentationPix*)gm->GetSegmentationByID(detType);
      Int_t lay,lad,det;
      Int_t ndig  = digArr->GetEntries();
      
      Int_t ndig_in_cycle=0; 
      Int_t ndig_signal_in_cycle=0;
      if (ndig<1) continue;
      gm->GetModuleId(imod, lay,lad,det);
//       printf("\nModule %3d: (det %2d in ladder %2d of Layer %d) | NDigits: %4d\n",imod,det,lad,lay,ndig);
      //    
      //if(iEvent==7 && lay==2) cout<<"Module #"<<imod<<endl;     //!!!DEBUG   
      //if(iEvent==7 && lay==2 && imod==363) debugOn=1; //continue;  //!!!DEBUG
      //else debugOn=0; 
      if(lay>=0 && lay <=2){
        Module = new AliITSUSuze02(Module_Nrows_small,Module_Ncols);
      }
      else{
        Module = new AliITSUSuze02(Module_Nrows_big,Module_Ncols);
      }
      
      Module->SetEncodingWindowSize(NRowsEncodingWindow,NColsEncodingWindow);
      Module->SetQuotas(nWindowsPer32colsMax,nWindowsPerHalfFSBBMax,nWindowsPerFSBBMax);
      
      OverflowCodesPerModule->Reset();
      nDigitsPerEncodingWindowPerModule->Reset();

      nWindows=0;
      for(Int_t idig=0;idig<ndig;idig++) {
        AliITSUDigitPix *pDig = (AliITSUDigitPix*)digArr->At(idig);
        if(pDig->GetROCycle()!=Cycle) continue; //selection of the hits from a given RO cycle
        if((CollectMode==-1 && pDig->GetHit(0)!=-1) || (CollectMode==1 && pDig->GetHit(0)==-1)) continue; // selection between noise and signal or both
        ndig_in_cycle++; 
        if(pDig->GetHit(0)!=-1) ndig_signal_in_cycle++;  //counts signal hits in a given RO cycle
        ColAddress=pDig->GetCoord1();
        RowAddress=pDig->GetCoord2();
        Module->AddDigit(RowAddress,ColAddress);
      }//diglist  
      
      if(AddQED){
        for(Int_t idig=0;idig<ndigQED;idig++) {
          AliITSUDigitPix *pDig = (AliITSUDigitPix*)digArrQED->At(idig);
          if(pDig->GetROCycle()!=Cycle) continue; //selection of the hits from a given RO cycle
          ndig_in_cycle++; 
          ColAddress=pDig->GetCoord1();
          RowAddress=pDig->GetCoord2();
          Module->AddDigit(RowAddress,ColAddress);
        }//diglist
      }
      
      if(ndig_in_cycle<1){  
        delete Module; 
        continue; //rejects when no hits in a given cycle
      }
      nDigitsPerEvent+=ndig_in_cycle;
      nDigitsPerEventPerLayer[lay]+=ndig_in_cycle;
      if(ProcessOnlyModulesWithSignal){       //if ProcessOnlyModulesWithSignal==1
        if(ndig_signal_in_cycle<1){   
          delete Module;
          continue;  //rejects when only noise is the present
        }
      }
      
      Module->Process(OverflowCodesPerModule,nDigitsPerEncodingWindowPerModule);
      DataSize=Module->GetDataSize();
      nDigitsEncodedPerModule=Module->GetNDigitsEncoded();
      nDigitsLostPerModule=Module->GetNDigitsLost();
      nWindows=Module->GetNEncodedWindows();                           
      nWindowsPer32colsMin=Module->GetNWindowsPer32colsMin();
      nWindowsPerHalfFSBBMin=Module->GetNWindowsPerHalfFSBBMin();
      nWindowsPerFSBBMin=Module->GetNWindowsPerFSBBMin();
      
     //DataSize=MakeSuze(Module_matrix,OverflowCodesPerModule,nDigitsEncodedPerModule,nDigitsLostPerModule,nWindows,nWindowsPer32colsMax,nWindowsPerHalfFSBBMax,nWindowsPerFSBBMax,nWindowsPer32colsMin,nWindowsPerHalfFSBBMin,nWindowsPerFSBBMin,nDigitsPerEncodingWindowPerModule);
//       cout<<"SUZE encoded "<<SuzeReturn<<" digits in "<<nWindows<<" windows"<<endl;
      OverflowCodes->Add(OverflowCodesPerModule);
      OverflowCodesPerLayer[lay]->Add(OverflowCodesPerModule);

      nDigitsPerEncodingWindow->Add(nDigitsPerEncodingWindowPerModule);
      nDigitsPerEncodingWindowPerLayer[lay]->Add(nDigitsPerEncodingWindowPerModule);

//       if(nDigitsLostPerModule){
//         cout<<"------- Some digits were lost. Check overflow errors"<<endl;
//         cout<<"Module has "<<ModuleSum<<" digits"<<endl;
//         cout<<"SUZE reported "<<SuzeReturn<<" encoded and "<<nDigitsLost<<" lost digits"<<endl;
//         cout<<"Module n."<<imod<<" has lost "<<nDigitsLostPerModule<<" digits due to the overflow"<<endl;

//       }
      if(nDigitsLostPerModule){
	      printf("Event #%d mod. #%d (layer %d) has %d digits lost\n",iEvent,imod,lay,nDigitsLostPerModule);
	      fprintf(logfile,"Event #%d mod. #%d (layer %d) has %d (%f) digits lost\n",iEvent,imod,lay,nDigitsLostPerModule, (Float_t)nDigitsLostPerModule/ndig_in_cycle);
      }
      nDigitsLostPerEvent+=nDigitsLostPerModule;
      nDigitsLostPerEventPerLayer[lay]+=nDigitsLostPerModule;
//       cout<<"Lay:"<<lay<<" Lad:"<<lad<<" current_ladder:"<<current_ladder<<" current layer:"<<current_layer<<" Digits:"<<SuzeReturn<<endl;
      if(lay!=current_layer){
        cout<<"Layer #"<<lay<<endl;
      }
      if(lad!=current_ladder || lay!=current_layer){
        current_ladder=lad;
	      current_layer=lay;
	      NWindowsPerLadder=0;
      }
      NWindowsPerLadder+=nWindows;
      if(NWindowsPerLadder>MaxNWindowsPerLadderPerLayerPerEvent[current_layer]){
	      MaxNWindowsPerLadderPerLayerPerEvent[current_layer]=NWindowsPerLadder;
// 	cout<<"---- MaxNWindowsPerLadderPerLayerPerEvent:"<<MaxNWindowsPerLadderPerLayerPerEvent[current_layer]<<" at layer:"<<current_layer<<endl;
      }
      nDigitsPerModulePerEvent->SetBinContent(imod+1,iEvent+1,ndig_in_cycle);
      nDigitsEncodedPerModulePerEvent->SetBinContent(imod+1,iEvent+1,nDigitsEncodedPerModule);
      nDigitsLostPerModulePerEvent->SetBinContent(imod+1,iEvent+1,nDigitsLostPerModule);
      FractionDigitsLostPerModulePerEvent->SetBinContent(imod+1,iEvent+1,(Double_t)nDigitsLostPerModule/ndig_in_cycle);
      nEncodingWindowsPerModulePerEvent->SetBinContent(imod+1,iEvent+1,nWindows);
      nDigitsPerEncodingWindowPerModulePerEvent->SetBinContent(imod+1,iEvent+1,(Double_t)nDigitsEncodedPerModule/nWindows); 

      nWindowsPerFSBBMinPerModulePerEvent->SetBinContent(imod+1,iEvent+1,nWindowsPerFSBBMax-nWindowsPerFSBBMin);
      nWindowsPerHalfFSBBMinPerModulePerEvent->SetBinContent(imod+1,iEvent+1,nWindowsPerHalfFSBBMax-nWindowsPerHalfFSBBMin);
      nWindowsPer32colsMinPerModulePerEvent->SetBinContent(imod+1,iEvent+1,nWindowsPer32colsMax-nWindowsPer32colsMin);

      DataSizePerModulePerEvent->SetBinContent(imod+1,iEvent+1,DataSize);   
      
      delete Module;
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
      MaxNWindowsPerLadderPerLayerPerEvent_hist[i]->SetBinContent(iEvent+1,MaxNWindowsPerLadderPerLayerPerEvent[i]);
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
    sprintf(ResultsFileName,"ScanDigits_v15_results_cycle_%d_EncWindow_%dx%d_SuzeLimitsVersion_%d_mode_%d-%d_QED_%d.root",Cycle,NRowsEncodingWindow,NColsEncodingWindow,SuzeLimitsVersion,CollectMode,ProcessOnlyModulesWithSignal,AddQED);
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

      MaxNWindowsPerLadderPerLayerPerEvent_hist[i]->Write();

      OverflowCodesPerLayer[i]->Write();
      nDigitsPerEncodingWindowPerLayer[i]->Write();
    }

    nDigitsPerModulePerEvent->Write();
    nDigitsEncodedPerModulePerEvent->Write();
    nDigitsLostPerModulePerEvent->Write();
    FractionDigitsLostPerModulePerEvent->Write();
    nEncodingWindowsPerModulePerEvent->Write();
    nDigitsPerEncodingWindowPerModulePerEvent->Write();

    nWindowsPerFSBBMinPerModulePerEvent->Write();
    nWindowsPerHalfFSBBMinPerModulePerEvent->Write();
    nWindowsPer32colsMinPerModulePerEvent->Write();

    DataSizePerModulePerEvent->Write();

    ResultsFile->Close();
  }
  fclose (logfile);
//   cout<<"Multiple Encodings:"<<nMultipleEncodings<<endl;
//   cout<<"Lost digits:"<<nDigitsLostTotal<<endl;

}
