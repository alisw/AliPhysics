////////////////////////////////////////////////////////////
// Test macro for AliITSPid class and tracking version V2 // 
////////////////////////////////////////////////////////////
void
AliITSSavePIDV2(Int_t evNumber1=0,Int_t evNumber2=0) {
  const char *filename="AliITStracksV2.root";
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");    loadlibs();
      } else {    delete gAlice;    gAlice=0;
  }

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);

TFile *fpid = new TFile("AliITStracksV2Pid.root","recreate");
 Float_t factor=0;
 if(gAlice!=0)factor=gAlice->Field()->Factor();
 if(factor==0.)factor=1.;
 AliKalmanTrack::SetConvConst(100/0.299792458/0.2/factor);

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
   AliITSPid pid;

   AliITStrackV2 *iotrack=0;
   for (Int_t i=0; i<nentr; i++) {
      AliITStrackV2 *iotrack=new AliITStrackV2;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
              Int_t    pcode=pid->GetPcode(iotrack);
              Float_t  signal=iotrack->GetdEdx();
	      //iotrack->Propagate(iotrack->GetAlpha(),3.,0.1/65.19*1.848,0.1*1.848);
	      iotrack->PropagateTo(3.,0.0028,65.19);

	      iotrack->PropagateToVertex();
	      Double_t xk,par[5]; iotrack->GetExternalParameters(xk,par);
	      Float_t lam=TMath::ATan(par[3]);
	      Float_t pt_1=TMath::Abs(par[4]);
	      Float_t mom=0.;
	      if( (pt_1*TMath::Cos(lam))!=0. ){ mom=1./(pt_1*TMath::Cos(lam)); }else{mom=0.;};
       pidtmp.fPcode=pcode;
       pidtmp.fSignal=signal;
       pidtmp.fMom=mom;
       itspidTree.Fill();
       delete iotrack;
   }// End for i (tracks)
 fpid->Write();
 }// End for nev (events)
 fpid->Close();
 file->Close();		 
 cout<<"File AliITStracksV2Pid.root written"<<endl;
 return;
///////////////////////////////////////////////////////
}








