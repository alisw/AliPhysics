///////////////////////////////////////////////////////////
// Test macro for AliITStracksV1Pid.root file            //
// JINR Dubna Jan 2002                                   //
///////////////////////////////////////////////////////////
void
AliITSScanPIDV2(Int_t evNumber1=0,Int_t evNumber2=0) {
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

TFile *fpid = new TFile("AliITStracksV2Pid.root","read");
fpid->ls();
//
//   Loop over events 
//
for (int nev=0; nev<= evNumber2; nev++) {
  char tpidname[30];
  sprintf(tpidname,"TreeT%d",nev);
  TTree *tracktree=(TTree*)fpid->Get(tpidname);
  TBranch *tbranch=tracktree->GetBranch("pids");
	
   Int_t nentr=tracktree->GetEntries();
   cout<<"Found PID for "<<nentr<<" ITS V1 tracks on "<<tpidname<<endl;

   AliITStrackV2Pid *iopid=0;
for(Int_t ii=0;ii<nentr;ii++)
  {
      AliITStrackV2Pid *iopid=new AliITStrackV2Pid;
      tbranch->SetAddress(&iopid);
      tracktree->GetEvent(ii);

      signal_mip->Fill(iopid->fSignal);

        if(iopid->fPcode ==2212)qplotP.Fill(1000*iopid->fMom,iopid->fSignal);
	if(iopid->fPcode == 321)qplotKa.Fill(1000*iopid->fMom,iopid->fSignal  );
	if(iopid->fPcode == 211)qplotPi.Fill(1000*iopid->fMom,iopid->fSignal  );
	if(iopid->fPcode ==  11)qplotE.Fill(1000*iopid->fMom,iopid->fSignal   );

      cout<<"PID pcode,fsignal,fmom= "
	  <<iopid->fPcode<<","<<iopid->fSignal<<","<<iopid->fMom<<endl;
      delete iopid;
  }// Enf for ii (tracks)
 }// End for nev (events)
 fpid->Close();
  //...................... Draw histogramms .................
   TCanvas *c1 = new TCanvas("PID_test","Scan PID ",200,10,900,700);
   c1->Divide(2,1);
  //.........................................................
   c1->cd(1); gPad->SetFillColor(33);
   signal_mip->Draw();

   c1->cd(2); //gPad->SetFillColor(33);
   qplot->Draw();
   qplotP.Draw("same"); qplotKa.Draw("same"); qplotPi.Draw("same"); qplotE.Draw("same");
   AliITSPid *pid =new AliITSPid(100);
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

