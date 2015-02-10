/// \file AliTPCScanPID.C
/// \brief Test macro for AliTPCtracksPid.root file
///
/// \author JINR Dubna
/// \date Aug 2002

void
AliTPCScanPID(Int_t evNumber2=3) {
  //................. Prepare histogramms ................

     TH2F *qplot =  new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13);
     TH2F *qplotP=  new TH2F("Qtrm1","Qtrm vs Pmom",100,0,1300,100,0,13); 
     TH2F *qplotKa= new TH2F("Qtrm2","Qtrm vs Pmom",100,0,1300,100,0,13);
     TH2F *qplotPi= new TH2F("Qtrm3","Qtrm vs Pmom",100,0,1300,100,0,13);
     TH2F *qplotE=  new TH2F("Qtrm4","Qtrm vs Pmom",100,0,1300,100,0,13);
     qplotP.SetMarkerStyle(8); qplotP.SetMarkerColor(kBlack); qplotP.SetMarkerSize(.3);
     qplotKa.SetMarkerStyle(8); qplotKa.SetMarkerColor(kRed); qplotKa.SetMarkerSize(.3);
     qplotPi.SetMarkerStyle(8); qplotPi.SetMarkerColor(kBlue); qplotPi.SetMarkerSize(.3);
     qplotE.SetMarkerStyle(8); qplotE.SetMarkerColor(kGreen); qplotE.SetMarkerSize(.3);
  //......................................................

TFile *fpid = new TFile("AliTPCtracksPid.root","read");
//
//   Loop over events 
//
for (int nev=0; nev< evNumber2; nev++) {
  char tpidname[30];
  sprintf(tpidname,"TreeT%d",nev);
  TTree *tracktree=(TTree*)fpid->Get(tpidname);
  TBranch *tbranch=tracktree->GetBranch("pids");
	
   Int_t nentr=tbranch->GetEntries(); 
   cout<<"Found PID for "<<nentr<<"  TPC tracks on "<<tpidname<<endl;

   AliTPCtrackPid *iopid=0;
for(Int_t ii=0;ii<nentr;ii++)
  {
      AliTPCtrackPid *iopid=new AliTPCtrackPid;
      tbranch->SetAddress(&iopid);
      tracktree->GetEvent(ii);

      Float_t xsignal=iopid->fSignal/50.;

        if(iopid->fPcode ==2212)qplotP.Fill(1000*iopid->fMom,xsignal);
	if(iopid->fPcode == 321)qplotKa.Fill(1000*iopid->fMom,xsignal  );
	if(iopid->fPcode == 211)qplotPi.Fill(1000*iopid->fMom,xsignal  );
	if(iopid->fPcode ==  11)qplotE.Fill(1000*iopid->fMom,xsignal   );

      //cout<<"PID pcode,fsignal,fmom= "
      //  <<iopid->fPcode<<","<<iopid->fSignal<<","<<iopid->fMom<<endl;
      delete iopid;
  }// Enf for ii (tracks)
 }// End for nev (events)
 fpid->Close();
  //...................... Draw histogramms .................
   TCanvas *c1 = new TCanvas("PID_test","Scan PID ",200,10,900,700);
   c1->Divide(1,1);
  //.........................................................
   c1->cd(1); gPad->SetFillColor(33);
   qplot->Draw();
   qplotP.Draw("same");qplotKa.Draw("same");qplotPi.Draw("same");qplotE.Draw("same");
   AliTPCPid *pid =new AliTPCPid(100);
    pid.fCutKa->SetLineColor(kRed);
    pid.fCutKa->Draw("same");
	
  c1->Print("PIDplot.ps");
  cout<<"End of file AliTPCtracksPid.root "<<endl; 

  return;
}


