{
   TObjArray tarray(2000);
   TObjArray parray(2000);
//{
  TFile *cf=TFile::Open("AliITSclustersV2.root");      assert(cf);
  AliITSgeom *geom=(AliITSgeom*)cf->Get("AliITSgeom"); assert(geom);   
  AliITStrackerV2 tracker(geom);                             
  cf->Close();	          
   TFile *tf=TFile::Open("AliITStracksV2.root");
   if (!tf->IsOpen()) {cerr<<"Can't open AliITStracksV2.root !\n"; return 3;}
//   TObjArray tarray(2000);
   TTree *tracktree=(TTree*)tf->Get("TreeT_ITS_0;1");
   if (!tracktree) {cerr<<"Can't get a tree with ITS tracks !\n"; return 4;}
   TBranch *tbranch=tracktree->GetBranch("tracks");
   Int_t nentr=(Int_t)tracktree->GetEntries(),i;
   for (i=0; i<nentr; i++) {
       AliITStrackV2 *iotrack=new AliITStrackV2;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);

       Int_t tpcLabel=iotrack->GetLabel();
       //       tracker.CookLabel(iotrack,0.);
       Int_t itsLabel=iotrack->GetLabel();
       if (itsLabel != tpcLabel) iotrack->SetLabel(-TMath::Abs(itsLabel));
       if (tpcLabel < 0)         iotrack->SetLabel(-TMath::Abs(itsLabel));
      tarray.AddLast(iotrack);
   }   
    cout<<" N rec. tracks="<<nentr<<endl;
   tf->Close();
//----------------------------------------
TFile *fpid = new TFile("AliITStracksV2Pid.root","recreate");
AliITStrackV2Pid pidtmp;
AliITSPid* pid = new AliITSPid(1000);
  TTree itsTree("ITSf","Tree with PID");
  AliITStrackV2Pid *outpid=&pidtmp;
  itsTree.Branch("pids","AliITStrackV2Pid",&outpid,32000,1);
Float_t signal,pmom;
Double_t xv,par[5];
AliITStrackV2* track;
Float_t lam,pt_1;
Int_t pcode;
for(Int_t ii=0;ii<nentr;ii++)
{
  track = (AliITStrackV2*)tarray[ii];
  track->Propagate(track->GetAlpha(),3.,0.1/65.19*1.848,0.1*1.848);
  track->PropagateToVertex();
    signal=track->GetdEdx();
    track->GetExternalParameters(xv,par);  
    lam=TMath::ATan(par[3]);
    pt_1=TMath::Abs(par[4]);
    //cout<<"lam,pt_1="<<lam<<" " <<pt_1<<endl;
    if(TMath::Abs(pt_1)>0)pmom=1.*(1./(pt_1*TMath::Cos(lam)));
  pidtmp.fSignal=signal;
  pcode=pid->GetPcode(signal,pmom);
  pidtmp.fPcode=pcode;
  pidtmp.fWpi=pidtmp.fWk=pidtmp.fWp=-1;
  if(pcode==211)pidtmp.fWpi=1;
  if(pcode==321)pidtmp.fWk=1;
  if(pcode==2212)pidtmp.fWp=1;

  pidtmp.fMom=pmom;
  itsTree.Fill();
}
itsTree.Write();
fpid->Close();
cout<<"File AliITStracksV2Pid.root written"<<endl; 
//----------------------------------------

TFile *fpid = new TFile("AliITStracksV2Pid.root");
fpid->ls();
fpid->Close();
//-------------------------------------------------
   TFile *tfpid=TFile::Open("AliITStracksV2Pid.root");
   if (!tfpid->IsOpen()) {cerr<<"Can't open AliITStracksV2Pid.root !\n"; return 3;}
//   TObjArray parray(2000);
   TTree *pidtree=(TTree*)tfpid->Get("ITSf");
   if (!pidtree) {cerr<<"Can't get a tree with ITS pid !\n"; return 4;}
   TBranch *tbranch=pidtree->GetBranch("pids");
   Int_t nentr=(Int_t)pidtree->GetEntries(),i;
   for (i=0; i<nentr; i++) {
       AliITStrackV2Pid *iopid=new AliITStrackV2Pid;
       tbranch->SetAddress(&iopid);
       pidtree->GetEvent(i);
      parray.AddLast(iopid);
      //cout<<" fWpi,k,p, Signal = "<<iopid->fWpi<<","<<iopid->fWk<<","<<iopid->fWp<<","<<iopid->fSignal  <<endl;
      //cout<<" Pcode,pmom="<<iopid->fPcode<<" "<<iopid->fMom<<endl;
   }   
    cout<<" N rec. pids="<<nentr<<endl;
   tfpid->Close();

//-------------------------------------------------
}























