void ITSSDDQAMaker(char *iFile, Int_t MaxEvts=1000000, Int_t FirstEvt=0) {

  //To have Baseline Histos uncomment parts with " // BL!!! " comment

cout << "SDD Quality Assurance Prototype Macro" << endl; 

const Int_t nSDDmodules= 260;
const Int_t imodoffset = 240;
const Int_t modtotSDD  = nSDDmodules*2;
const Int_t anode = 256;

Float_t xi = -0.5;
Float_t xf = xi + nSDDmodules;
TH1F *hModulePattern = new TH1F("hModulePattern","Modules pattern",nSDDmodules,xi,xf);
xf = xi + modtotSDD;
TH1F *hModuleSidePattern = new TH1F("hModuleSidePattern","Modules/Side pattern",modtotSDD,xi,xf);

TH2F *hModuleChargeMap[modtotSDD];  //260 dx e 260 sx  with A, T, Q
TH2F *hModuleCountsMap[modtotSDD];  //260 dx e 260 sx  with A, T, Ncounts
TH2I *hModuleCarlos = new TH2I("hModuleCarlos","hModuleCarlos",modtotSDD,xi,xf,101,-0.5,100.5);
/*
  TH1F *hModuleSideBL[modtotSDD][anode];                                // BL!!!
*/ 


//-------histograms definition 
Char_t *hisnam = new Char_t[50];
Char_t *histit = new Char_t[50];  
Char_t *hisnam2 = new Char_t[50];
Char_t *histit2 = new Char_t[50];    
Char_t *hisnam3 = new Char_t[50];
 for(Int_t imod=0; imod<nSDDmodules;imod++){
   for(Int_t isid=0;isid<2;isid++){
     Int_t index=2*imod+isid;       //260*2 position

     sprintf(hisnam,"chargeMap%d",index);
     sprintf(histit,"Total Charge, module number %d",index);
     hModuleChargeMap[index]=new TH2F(hisnam,histit,256,-0.5,255.5,256,-0.5,255.5);  

     sprintf(hisnam2,"countsMap%d",index);
     sprintf(histit2,"Number of Counts, module number %d",index);
     hModuleCountsMap[index] = new TH2F(hisnam2,histit2,256,-0.5,255.5,256,-0.5,255.5);
     /*
       for(Int_t ianode=0; ianode<anode; ianode++){                         // BL!!!
       sprintf(hisnam3,"BL_module_%d_%d",index,ianode); 
       //cout<<hisnam3 <<endl;
       hModuleSideBL[index][ianode] = new TH1F(hisnam3,hisnam3,256,0.,1024.);
     }
     */
   }
 }

  AliRawReader *rd = new AliRawReaderDate(iFile,FirstEvt);  // open run
  Int_t evCounter = 0;
  Int_t eqOffset = 256;
  Int_t DDLid_range = 24;
  do{                       // start loop on events
    if(++evCounter > MaxEvts) { cout << MaxEvts << " events read, stop" << endl; evCounter--; break; }  
    cout << "Read Event: " << evCounter+FirstEvt-1 << endl;

    rd->RequireHeader(kFALSE);             
    rd->SelectEvents(7);                  
    rd->SelectEquipment(17,eqOffset+1,eqOffset+DDLid_range);  //17 states for "DRorc acquisition"
    rd->Reset();                           // reset the current position to the beginning of the event
 
    Int_t nSkip = 0;                     // number of skipped signals
    AliITSRawStreamSDD s(rd);            //This class provides access to ITS SDD digits in raw data.
    Int_t iddl;
    Int_t isddmod;
    Int_t moduleSDD;
    gStyle->SetPalette(1);
    while(s.Next()){                     //read the next raw digit; returns kFALSE if there is no digit left

      iddl=rd->GetDDLID(); // careful: it starts from 1!!!
	isddmod=s.GetModuleNumber(iddl-1,s.GetCarlosId());        //this is the FEE Carlos
       	//cout<<"DDLID= "<<iddl <<"; Module number= " <<isddmod  <<endl;
	if(isddmod >= imodoffset) { 
	  hModulePattern->Fill(isddmod-imodoffset);             // 0 to 259    so  240 to 499
	  moduleSDD=2*(isddmod-imodoffset)+s.GetChannel();
          hModuleSidePattern->Fill(moduleSDD);                  // 0 to 519
	  hModuleCarlos->Fill(isddmod-imodoffset,s.GetCarlosId());
          //cout << "anode " << s.GetCoord1() << ", time bin: " << s.GetCoord2() << ", charge: " << s.GetSignal() << endl;	  
	  Int_t coord1 = s.GetCoord1();
 	  Int_t coord2 = s.GetCoord2();
 	  Int_t signal = s.GetSignal();
	  hModuleChargeMap[moduleSDD]->Fill(coord2, coord1,signal);
          hModuleCountsMap[moduleSDD]->Fill(coord2, coord1 );   
	  //hModuleSideBL[moduleSDD][coord1]->Fill(signal);             // BL  !!!
	} else {
	  nSkip++;
	}
    }    
    cout << "End of Event " << evCounter+FirstEvt-1 << ", " << nSkip << " wrong module numbers" << endl;
  } while(rd->NextEvent()); // end loop on events
  delete rd;

  cout << "end after " << evCounter << " events" << endl;
  /*
  TNtuple *Baseline = new TNtuple("Baseline","Baseline","HalfModule:Anode:Mean:RMS");      // BL!!!
  Float_t meanBL;
  Float_t rmsBL;  
  */
  for(Int_t i=0; i<modtotSDD; i++){   
    if(hModuleSidePattern->GetBinContent(i+1)){              //check if they're not empty
      hModuleChargeMap[i]->GetXaxis()->SetTitle("Time Bin");
      hModuleChargeMap[i]->GetYaxis()->SetTitle("Anode"); 
      hModuleCountsMap[i]->GetXaxis()->SetTitle("Time Bin");
      hModuleCountsMap[i]->GetYaxis()->SetTitle("Anode");  
      /*
      for(Int_t ianode=0; ianode<anode; ianode++ ){                                      // BL!!!
	hModuleSideBL[i][ianode]->GetXaxis()->SetTitle("ADC counts");
	hModuleSideBL[i][ianode]->GetYaxis()->SetTitle("#"); 
	meanBL = hModuleSideBL[i][ianode]->GetMean();
	rmsBL = hModuleSideBL[i][ianode]->GetRMS();
	gaussfitBL = hModuleSideBL[i][ianode]->Fit("gaus");	
	Baseline->Fill(i,ianode,meanBL,rmsBL);
      }
      */
    }
  }  
  
  hModuleSidePattern->GetXaxis()->SetTitle("2*(Module Number-1)+Side");
  hModuleSidePattern->GetYaxis()->SetTitle("Counts");  
  hModulePattern->GetXaxis()->SetTitle("Module Number");  
  hModulePattern->GetYaxis()->SetTitle("Counts");  


  //-------store Histograms
  cout << "Store Histograms" << endl;
  TString oFileName(iFile);
  oFileName.Append(".root");
  TFile *oFile = TFile::Open(oFileName,"recreate");
  hModulePattern->Write();
  hModuleSidePattern->Write();
  hModuleCarlos->Write();
  for(Int_t i=0; i<modtotSDD; i++){ 
      if(hModuleSidePattern->GetBinContent(i+1)){            //check if they're not empty
	hModuleChargeMap[i]->Write();
	hModuleCountsMap[i]->Write();     
	/* 
	for(Int_t ianode=0; ianode<anode; ianode++ ){                           // BL!!!
 	  hModuleSideBL[i][ianode]->Write();
	  Baseline->Write();
 	}
	*/
      }
  }
  
  oFile->Close();
  cout << "Clear memory" << endl;
  for(Int_t imod=0; imod<nSDDmodules;imod++){
    for(Int_t isid=0;isid<2;isid++){
      Int_t index=2*imod+isid;       //260*2 position
      delete hModuleChargeMap[index];
      delete hModuleCountsMap[index];       
      /*
      for(Int_t ianode=0; ianode<anode; ianode++ ){                              // BL!!!
	delete hModuleSideBL[index][ianode]; 
	delete Baseline;
      }
      */
    }
  }
  delete hModulePattern;
  delete hModuleSidePattern;
  delete hModuleCarlos;

}
