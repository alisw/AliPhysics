///////////////////////////////////////////////////////////
// Test macro for AliITStracksV2Pid.root file            //
// JINR Dubna Jan 2002                                   //
///////////////////////////////////////////////////////////
void
AliITSScanPIDV2(Int_t evNumber1=0,Int_t evNumber2=0) {
  //................. Prepare histogramms ................
     TH2F *qplot =  new TH2F("Qtrm","Qtrm vs Pmom",100,0,1.300,100,0,13);
     TH2F *qplotP=  new TH2F("QtrmP","Qtrm vs Pmom",100,0,1.300,100,0,13); 
     TH2F *qplotKa= new TH2F("QtrmKa","Qtrm vs Pmom",100,0,1.300,100,0,13);
     TH2F *qplotPi= new TH2F("QtrmPi","Qtrm vs Pmom",100,0,1.300,100,0,13);
     TH2F *qplotE=  new TH2F("QtrmE","Qtrm vs Pmom",100,0,1.300,100,0,13);
     qplotP.SetMarkerStyle(8); qplotP.SetMarkerColor(kBlue); qplotP.SetMarkerSize(.3);
     qplotKa.SetMarkerStyle(8); qplotKa.SetMarkerColor(kRed); qplotKa.SetMarkerSize(.3);
     qplotPi.SetMarkerStyle(8); qplotPi.SetMarkerColor(kBlack); qplotPi.SetMarkerSize(.3);
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
   cout<<"Found PID for "<<nentr<<" ITS V2 tracks on "<<tpidname<<endl;

   AliITStrackV2Pid *iopid=0;
for(Int_t ii=0;ii<nentr;ii++)
  {
      AliITStrackV2Pid *iopid=new AliITStrackV2Pid;
      tbranch->SetAddress(&iopid);
      tracktree->GetEvent(ii);

      signal_mip->Fill(iopid->fSignal);

        if(iopid->fPcode ==2212)qplotP.Fill(iopid->fMom,iopid->fSignal);
	if(iopid->fPcode == 321)qplotKa.Fill(iopid->fMom,iopid->fSignal  );
	if(iopid->fPcode == 211)qplotPi.Fill(iopid->fMom,iopid->fSignal  );
	if(iopid->fPcode ==  11)qplotE.Fill(iopid->fMom,iopid->fSignal   );
	/*
	
if(  (iopid->fWp<0.10)||(iopid->fWk<0.0)||(iopid->fWpi<0.0) ){
          cout<<"PID pcode,fsignal,fmom= "<<iopid->fPcode<<","<<iopid->fSignal<<","<<iopid->fMom<<endl;
	  cout<<"wpi,wka,wp="<<iopid->fWpi<<" "<<iopid->fWk<<" "<<iopid->fWp<<endl;
      }
	*/
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
   fcutka.Draw("same"); fcutpr.Draw("same");
   c1->Print("ITSPIDplot.ps");

  cout<<"End of file AliITStracksV2Pid.root "<<endl; 
  return;
}

