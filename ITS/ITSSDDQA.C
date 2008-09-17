void ITSSDDQA(char *iFile, Int_t MaxEvts=1000000, Int_t FirstEvt=0) {

cout << "SDD Quality Assurance Prototype Macro" << endl; 

const Int_t nSDDmodules= 260;
const Int_t imodoffset = 240;
const Int_t modtotSDD  = nSDDmodules*2;

Float_t xi = 0.5;
Float_t xf = xi + nSDDmodules;
TH1F *modulePattern = new TH1F("patternModule","Modules pattern",nSDDmodules,xi,xf);
xf = xi + modtotSDD;
TH1F *moduleSidePattern = new TH1F("patternSide","Modules/Side pattern",modtotSDD,xi,xf);

TH2F *hismap[modtotSDD];  //260 dx e 260 sx  with A, T, Q
TH2F *hispop[modtotSDD];  //260 dx e 260 sx  with A, T, Ncounts
 
Char_t *cindex = new Char_t[5];
for(Int_t imod=0; imod<nSDDmodules;imod++){
  for(Int_t isid=0;isid<2;isid++){
    Int_t index=2*imod+isid;       //260*2 position

    sprintf(cindex,"%d",index+1); // imod,isid);
    TString sindex((const char *) cindex);

    TString histnam("chargeMap");
    TString histit("Total Charge, module number ");
    histnam.Append(sindex);
    histit.Append(sindex);
    hismap[index]=new TH2F(histnam.Data(),histit.Data(),256,-0.5,255.5,256,-0.5,255.5);
      
    TString hisnam2("countsMap");
    TString histit2("Number of Counts, module number ");
    hisnam2.Append(sindex);
    histit2.Append(sindex);
    hispop[index]=new TH2F(hisnam2.Data(),histit2.Data(),256,-0.5,255.5,256,-0.5,255.5);
  
    /*
      sprintf(hisnam,"hisprojX%03ds%d",imod,isid);
      sprintf(histitle,"layer , ladder, module position, channel, %d, %d", imod, isid);
      hisprojX[index]=new TH2F(hisnam,histitle,256,-0.5,255.5,256,-0.5,255.5);

      sprintf(hisnam,"hisprojY%03ds%d",imod,isid);
      sprintf(histitle,"layer , ladder, module position, channel, %d, %d", imod, isid);
      hisprojY[index]=new TH2F(hisnam,histitle,256,-0.5,255.5,256,-0.5,255.5);

      sprintf(hisnam,"hisprofX%03ds%d",imod,isid);
      sprintf(histitle,"layer , ladder, module position, channel, %d, %d", imod, isid);
      hisprofX[index]=new TH2F(hisnam,histitle,256,-0.5,255.5,256,-0.5,255.5);

      sprintf(hisnam,"hisprofY%03ds%d",imod,isid);
      sprintf(histitle,"layer , ladder, module position, channel, %d, %d", imod, isid);
      hisprofY[index]=new TH2F(hisnam,histitle,256,-0.5,255.5,256,-0.5,255.5);
      */
 
  }
}
delete [] cindex;

  TString strFile = iFile;
  strFile += "?EventType=7";
  AliRawReader *rd = new AliRawReaderDate(strFile.Data(),FirstEvt);  // open run
  Int_t evCounter = 0;

  //AliITS *itsRun = new AliITS();
  //  TGeoManager::Import("$ALICE_ROOT/EVE/alice-data/alice_fullgeo.root");
  /*
    AliITSInitGeometry *initgeom = new AliITSInitGeometry("AliITSvPPRasymmFMD",2);
    geom = initgeom->CreateAliITSgeom();
    delete initgeom;
  */

  do{       // start loop on events
    if(++evCounter > MaxEvts) { cout << MaxEvts << " events read, stop" << endl; evCounter--; break; }  
     cout << "Read Event: " << evCounter+FirstEvt-1 << endl;

    rd->RequireHeader(kFALSE);             


    rd->Reset();    // reset the current position to the beginning of the event
    AliITSRawStreamSDD s(rd);    //This class provides access to ITS SDD digits in raw data.
    Int_t iddl;
    Int_t isddmod;
    Int_t moduleSDD;
    while(s.Next()){                       //read the next raw digit; returns kFALSE if there is no digit left

	iddl=rd->GetDDLID();
	isddmod=s.GetModuleNumber(iddl,s.GetCarlosId());  //this is the FEE Carlos
	//cout<<"DDLID= "<<iddl <<"; Module number= " <<isddmod  <<endl;
	modulePattern->Fill(isddmod-imodoffset+1);     // 1 to 260
	moduleSDD=2*(isddmod-imodoffset)+s.GetChannel();
        moduleSidePattern->Fill(moduleSDD+1);          // 1 to 520
        //cout << "anode " << s.GetCoord1() << ", time bin: " << s.GetCoord2() << ", charge: " << s.GetSignal() << endl;
        hismap[moduleSDD]->Fill(s.GetCoord1(), s.GetCoord2(),s.GetSignal() );
        hispop[moduleSDD]->Fill(s.GetCoord1(), s.GetCoord2() );	 
    }    

  } while(rd->NextEvent()); // end loop on events
  delete rd;

  cout << "end after " << evCounter << " events" << endl;

   TString oFileName(iFile);
   oFileName.Append(".root");
   TFile *oFile = TFile::Open(oFileName,"recreate");
   modulePattern->Write();
   moduleSidePattern->Write();
   for(Int_t i=0; i<modtotSDD; i++){
     hismap[i]->Write();
     hispop[i]->Write();
   }
   oFile->Close();

   for(Int_t imod=0; imod<nSDDmodules;imod++){
     for(Int_t isid=0;isid<2;isid++){
       Int_t index=2*imod+isid;       //260*2 position
       delete hismap[index];
       delete hispop[index];
     }
   }
   delete modulePattern;
   delete moduleSidePattern;

}
