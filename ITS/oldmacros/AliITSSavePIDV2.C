////////////////////////////////////////////////////////////
// Test macro for AliITSPid class and tracking version V2 // 
// Rev. 25 July 2002 v3-08-03 + Root 3.02/07 + RedHat 6.2 //
////////////////////////////////////////////////////////////
void
AliITSSavePIDV2(Int_t evNumber1=0,Int_t evNumber2=0) {
  const char *filename="AliITStracksV2.root";

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);

TFile *fpid = new TFile("AliITStracksV2Pid.root","recreate");


 Float_t factor=0;
 if(gAlice!=0)factor=gAlice->Field()->SolenoidField();
 if(factor==0.){
   cout<<" ================  WARNING ====================================\n";
   cout<<"  The default magnetic field value of 0.4 T will be used\n";
   cout<<" ==============================================================\n";
   factor=4.;  // Default  mag. field = 0.4T
 }
 else {
   cout<<"AliITSSavePIDV2.C:  Magnetic field is "<<factor/10.<<" T\n";
 }
 AliKalmanTrack::SetConvConst(1000/0.299792458/factor);



 AliITSPid *pid=new AliITSPid(100);
//
//   Loop over events 
//
 for (int nev=0; nev<= evNumber2; nev++) {
   char tname[30];
   sprintf(tname,"TreeT_ITS_%d;1",nev);
   TTree *tracktree=(TTree*)file->Get(tname);
   TBranch *tbranch=tracktree->GetBranch("tracks");
	
   Int_t nentr=tracktree->GetEntries();
   cout<<"Found "<<nentr<<" ITS tracks in event No "<<nev<<endl;

   char tpidname[30];
   sprintf(tpidname,"TreeT%d",nev);
   AliITStrackV2Pid pidtmp;
   TTree itspidTree(tpidname,"Tree with PID");
   AliITStrackV2Pid *outpid=&pidtmp;
   itspidTree.Branch("pids","AliITStrackV2Pid",&outpid,32000,1);
   AliITStrackV2 *iotrack=0;
   for (Int_t i=0; i<nentr; i++) {
      AliITStrackV2 *iotrack=new AliITStrackV2;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
              iotrack->CookdEdx();
              Float_t  signal=iotrack->GetdEdx();
	      iotrack->PropagateTo(3.,0.0028,65.19);
	      iotrack->PropagateToVertex();
	      Double_t xk,par[5]; iotrack->GetExternalParameters(xk,par);
	      Float_t lam=TMath::ATan(par[3]);
	      Float_t pt_1=TMath::Abs(par[4]);
	      Float_t mom=0.;
	      if( (pt_1*TMath::Cos(lam))!=0. ){ mom=1./(pt_1*TMath::Cos(lam)); }else{mom=0.;};
	      Float_t phi=TMath::ASin(par[2]) + iotrack->GetAlpha();
	      if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
	      if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
	      signal=signal/35.;
       pidtmp.fSignal=signal;
       pidtmp.fMom=mom;
       pidtmp.fPhi=phi;
       pidtmp.fLam=lam;
       pidtmp.fGlab=TMath::Abs(iotrack->GetLabel());
       Int_t pcode=pid->GetPcode(signal,mom);
       pidtmp.fPcode=pcode;
       pidtmp.fWpi=pid->GetWpi();
       pidtmp.fWk=pid->GetWk();
       pidtmp.fWp=pid->GetWp();
       //cout<<" pcode,sig,mom="<<pcode<<" "<<signal<<" "<<mom<<endl;
       //cout<<" wpi,wka,wp="<<pidtmp.fWpi<<" "<<pidtmp.fWk<<" "<<pidtmp.fWp<<endl;
       itspidTree.Fill();
       delete iotrack;
   }// End for i (tracks)
   cout<<"n ev="<<nev<<endl;
 fpid->Write();
 }// End for nev (events)
 file->Close();		 
 cout<<"File AliITStracksV2Pid.root written"<<endl;
 delete file;
 cout<<"end of AliITSSavePIDV2.C "<<endl;
 //return;
///////////////////////////////////////////////////////
}








