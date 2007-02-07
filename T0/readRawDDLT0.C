void readRawDDLT0()
{

  TH1F * hCFD[8]; TH1F *hLED[8]; TH1F*hQT01[4];TH1F*hQT02[4];
  TH1F*hQTD[4]; TH1F*hADC[8]; TH2*hQTCFD[4];
  Char_t buf1[10], buf2[10], buf3[10], buf4[10], buf5[10], buf6[10],buf7[10];
  for (Int_t ic=0; ic<8; ic++)
    {
      sprintf(buf1,"CFD%i",ic+1);
      hCFD[ic]= new TH1F(buf1,"CFD",500,0.5,6000.5);
      sprintf(buf2,"LED%i",ic+1);
      hLED[ic]= new TH1F(buf2,"LED",500,-0.5,2000.5);
      //LED-CFD
      sprintf(buf6,"ADC_%i",ic+1);
      hADC[ic]= new TH1F(buf6,"QT1",500,0,6000);
    }
  
  for (Int_t iq=0; iq<4; iq++)
    {
      //QT01 - QT04
      sprintf(buf3,"QT0%i",iq+1);
      hQT01[iq]= new TH1F(buf3,"QT01",500,0,6000);
      sprintf(buf4,"QT1_%i",iq+1);
      //QT11 - QT14 
      hQT02[iq]= new TH1F(buf4,"QT02",500,0,6000);
      sprintf(buf5,"QTD_%i",iq+1);
      //QT11-QT01 ....
      hQTD[iq]= new TH1F(buf5,"QTdiff",5000,1500,6000);

      sprintf(buf7,"QTCFD_%i",iq+1);
      hQTCFD[iq]= new TH2F(buf7,"QT vs CFD",500,0,6000,500,0,5000);
  }

  TH1F *hORA= new TH1F("hORA"," T0 A ",500, 1000,2000);
  TH1F *hORC= new TH1F("hORC"," T0 C ",500, 1000,2000);
  TH1F*hEffCFD= new TH1F("hEffCFD","Effeciency",8,0.25,4.25);

  //  TH2F*hQTCFD= new TH2F("hQTCFD","QT vs CFD",500,0.5,6000.5,500,0.5,5000.5);
  TH2F*hQTLED= new TH2F("hQTLED","QT vs LED",500,0.5,6000.5,500,0.5,5000.5);
  TH2F*hLEDCFD= new TH2F("hLEDCFD","LEd vs CFD",500,-0.5,10000.5,500,-0.5,10000.5);

  AliCDBManager* cdb      = AliCDBManager::Instance();
  AliCDBStorage *stor = cdb->GetStorage("local://$ALICE_ROOT");

  // AliRawReader *reader = new AliRawReaderDate("raw0/T0_3328.ddl");
   AliRawReader *reader = new AliRawReaderFile();
   //  reader->LoadEquipmentIdsMap("T0map.txt");
   // reader->RequireHeader(kFALSE);
  AliT0RawReader *start = new AliT0RawReader(reader);
   Int_t allData[110][5];

   while (reader->NextEvent()) {
     start->Next();
     for (Int_t i=0; i<110; i++) {
       allData[i][0]= start->GetData(i,0);
       cout<<" i "<< i<<" "<<allData[i][0] <<endl;
     }
   }
   
   /*
  Char_t filedata[20];

  sprintf(filedata,"START%i.root",runNumber);
  TFile *f1 = new TFile(filedata);
  //  f1->ls();
  TTree *tree = (TTree *)f1->Get("raw");
  //  t->Dump();
  AliSTARTdigit *digit = 0;
  tree->SetBranchAddress("START",&digit);
  Int_t nevent=tree->GetEntries();
  for(Int_t i=0; i<nevent; i++)
    {
      tree->GetEntry(i);
      Int_t time = digit->CFD01();
      hCFD[0]->Fill(time);
      time = digit->CFD02();
      hCFD[1]->Fill(time);
      time = digit->CFD03();
      hCFD[2]->Fill(time);
      time = digit->CFD04();
      hCFD[3]->Fill(time);
      time = digit->CFD11();
      hCFD[4]->Fill(time);
      time = digit->CFD12();
      hCFD[5]->Fill(time);
      time = digit->CFD13();
      hCFD[6]->Fill(time);
      time = digit->CFD14();
      hCFD[7]->Fill(time);
      time = digit->LED01();
      hLED[0]->Fill(time);
      time = digit->LED02();
      hLED[1]->Fill(time);
      time = digit->LED03();
      hLED[2]->Fill(time);
      time = digit->LED04();
      hLED[3]->Fill(time);
      time = digit->LED05();
      hLED[4]->Fill(time);
      time = digit->LED06();
      hLED[5]->Fill(time);
      time = digit->LED07();
      hLED[6]->Fill(time);
      time = digit->LED08();
      hLED[7]->Fill(time);

      Int_t led=(digit->LED01()-digit->CFD01());
      hLEDCFD->Fill(time,digit->CFD01());
      
      time = digit->QT01();
      hQT01[0]->Fill(time);
      time = digit->QT02();
      hQT01[1]->Fill(time);
      time = digit->QT03();
      hQT01[2]->Fill(time);
      time = digit->QT04();
      hQT01[3]->Fill(time);

      time = digit->QT11();
      hQT02[0]->Fill(time);
      time = digit->QT12();
      hQT02[1]->Fill(time);
      time = digit->QT13();
      hQT02[2]->Fill(time);
      time = digit->QT04();
      hQT02[3]->Fill(time);
      
      
      time = (digit->QT11() - digit->QT01());
      hQTD[0]->Fill(time);
      hQTLED->Fill(led,time);
      hQTCFD[0]->Fill(time,digit->CFD01());

      time = (digit->QT12() - digit->QT02());
      hQTD[1]->Fill(time);
      hQTCFD[1]->Fill(time,digit->CFD02());

      time = (digit->QT13() - digit->QT03());
      hQTD[2]->Fill(time);
      hQTCFD[2]->Fill(time,digit->CFD03());
      
      time = (digit->QT14() - digit->QT04());
      hQTD[3]->Fill(time);
      hQTCFD[3]->Fill(time,digit->CFD04());
      
      time = (digit->LED01() - digit->CFD01());
      hADC[0]->Fill(time);
      time = (digit->LED02() - digit->CFD02());
      hADC[1]->Fill(time);
      time = (digit->QT03() - digit->CFD03());
      hADC[2]->Fill(time);
      time = (digit->QT04() - digit->CFD04());
      hADC[3]->Fill(time);
      time = (digit->LED05() - digit->CFD11());
      hADC[4]->Fill(time);
      time = (digit->LED06() - digit->CFD12());
      hADC[5]->Fill(time);
      time = (digit->LED07() - digit->CFD13());
      hADC[6]->Fill(time);
      time = (digit->LED08() - digit->CFD14());
      hADC[7]->Fill(time);
     
      time = digit->ORA();
      hORA->Fill(time);
      time = digit->ORC();
      hORC->Fill(time);
      
    }
  cout<<nevent<<endl;
  for (Int_t iq=0; iq<4; iq++)
    {  
      Float_t NCFD = 0; 
      
      Float_t NCFD =  hCFD[iq]->Integral(hCFD[iq]->GetXaxis()->FindBin(200.5),
					 hCFD[iq]->GetXaxis()->FindBin(600.5));
      cout<<NCFD<<endl;
      
      hEffCFD->Fill(iq+1,NCFD/nevent);
      
    }

  Char_t filehist[20]; 
  sprintf(filehist,"hist%i.root",runNumber);
  TFile *hist = new TFile(filehist,"RECREATE");

  hLEDCFD->Write();
  for (i=0; i<8; i++)
    {
      hCFD[i]->Write();
      hLED[i]->Write();
      hADC[i]->Write();
    }
  for (i=0; i<4; i++)
    {
      hQT01[i]->Write();
      hQT02[i]->Write();
      hQTD[i]->Write();
      hQTCFD[i]->Write();
    }
  hEffCFD->Write();
  hORA->Write();
  hORC->Write();
  hQTLED->Write();
   */

}

  
