///////////////////////////////////////////////////////////
// Test macro for AliITStracksV1Pid.root file            //
// JINR Dubna Jan 2002                                   //
///////////////////////////////////////////////////////////
void
AliITSScanPIDV1(Int_t evNumber1=0,Int_t evNumber2=0) {
  //................. Prepare histogramms ................
     TH2F *qplot =  new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13);
     TH2F *qplotP=  new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13); 
     TH2F *qplotKa= new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13);
     TH2F *qplotPi= new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13);
     TH2F *qplotE=  new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13);
     qplotP.SetMarkerStyle(8); qplotP.SetMarkerColor(kBlack); qplotP.SetMarkerSize(.3);
     qplotKa.SetMarkerStyle(8); qplotKa.SetMarkerColor(kRed); qplotKa.SetMarkerSize(.3);
     qplotPi.SetMarkerStyle(8); qplotPi.SetMarkerColor(kBlue); qplotPi.SetMarkerSize(.3);
     qplotE.SetMarkerStyle(8); qplotE.SetMarkerColor(kGreen); qplotE.SetMarkerSize(.3);
  //......................................................
  TH1F *signal_mip = new TH1F("signal_mip","Signal (mips) for track",100,0.,15.);

//*****************************************************************************************************************************************

  const char *filename="itstracks.root";
  
  ///////////////// Dynamically link some shared libs ////////////////////////////////
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);

//
//   Loop over events 
//
   char tname[30];
   for (int nev=evNumber1; nev<= evNumber2; nev++) {
 

  // for (int nev=0; nev<= evNumber2; nev++) {

   sprintf(tname,"TreeT%d",nev);
   TTree *tracktree=(TTree*)file->Get(tname);
   TBranch *tbranch=tracktree->GetBranch("ITStracks");
   cout<<" nev = "<<nev<<"\n";
	//cout<<" open the file \n"; 
	
   Int_t nentr=tracktree->GetEntries();

   TObjArray tarray(nentr);
  // AliITSIOTrack *iotrack=0;
   printf("nentr %d\n",nentr);
	
   for (Int_t i=0; i<nentr; i++) {
      AliITSIOTrack *iotrack=new AliITSIOTrack;
      // tarray.AddAt(new AliITSIOTrack,i);
      // iotrack=(AliITSiotrack*)tarray.UncheckedAt(i);
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       tarray.AddLast(iotrack);
   }
   //file->Close();		 
	
	  AliITSIOTrack *iotrack;
   for (Int_t i=0; i<nentr; i++) {
         AliITSIOTrack *iotrack=new AliITSIOTrack;   	
	 iotrack=(AliITSIOTrack*)tarray.UncheckedAt(i);
	 if(!iotrack) continue;
	  
	 Float_t Px=iotrack->GetPx();
	 Float_t Py=iotrack->GetPy();
	 Float_t Pz=iotrack->GetPz();
	  
	 Float_t momentum = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
	 Float_t dEdx = iotrack->GetdEdx();
	 Int_t pcode = 211;//iotrack->GetPDG();

	 signal_mip->Fill(dEdx);
	 
         if(pcode == 2212) qplotP.Fill(1000*momentum,dEdx);
	 if(pcode ==  321) qplotKa.Fill(1000*momentum,dEdx);
	 if(pcode ==  211) qplotPi.Fill(1000*momentum,dEdx);
	 if(pcode ==   11) qplotE.Fill(1000*momentum,dEdx);
	 
    delete iotrack;		 
   }  

   }   // event loop 
   file->Close();   

//*****************************************************************************************************************************************
  
  //...................... Draw histogramms .................
   TCanvas *c1 = new TCanvas("PID_test","Scan PID ",200,10,900,700);
   c1->Divide(2,1);
  //.........................................................
   c1->cd(1); gPad->SetFillColor(33);
   signal_mip->Draw();

   c1->cd(2); //gPad->SetFillColor(33);
   qplot->Draw();
   qplotP.Draw("same"); qplotKa.Draw("same"); qplotPi.Draw("same"); qplotE.Draw("same");
   AliITSPid *pid =new AliITSPid(1000);
    c1->Range(0,0,1300,10);
    gStyle->SetLineColor(kRed);
    gStyle->SetLineWidth(2);
	TLine *lj[3],*lk[3]; 
	for(Int_t j=0;j<3;j++){
		Float_t x1,x2,y1,y2,xx1,xx2,yy1,yy2;
		x1=pid->cut[j+1][0]; x2=pid->cut[j+2][0];
		y1=y2=pid->cut[j+2][2];
	    lj[j]=new TLine(1000*x1,y1,1000*x2,y2); lj[j]->Draw();
	    if(j==0){yy1=10.;}else{yy1=lj[j-1]->GetY1();}
	    yy2=lj[j]->GetY1();
	    xx1=xx2=x1;
	    lk[j]=new TLine(1000*xx1,yy1,1000*xx2,yy2); lk[j]->Draw();
	}
	//Draw pions-kaons cuts.
	TLine *mj[7],*mk[7]; 
	for(Int_t j=0;j<7;j++){
		Float_t x1,x2,y1,y2,xx1,xx2,yy1,yy2;
		x1=pid->cut[j+2][0]; x2=pid->cut[j+3][0];
		y1=y2=pid->cut[j+3][5];
	    mj[j]=new TLine(1000*x1,y1,1000*x2,y2); mj[j]->Draw();
	    if(j==0){yy1=10.;}else{yy1=mj[j-1]->GetY1();}
	    yy2=mj[j]->GetY1();
	    xx1=xx2=x1;
	    mk[j]=new TLine(1000*xx1,yy1,1000*xx2,yy2); mk[j]->Draw();
	}
  cout<<"End of file AliITStracksV1Pid.root "<<endl; 
  return;
}

