/// \file AliTPCSavePID.C
///
/// \author Dubna
/// \date 2003, 22 Jan

Int_t AliTPCSavePID(Int_t emax=3) {
///

#include<fstream.h>
TFile *inkin = TFile::Open("galice.root");
if (!inkin->IsOpen()) {
  cerr<<"Can't open galice.root !\n";
}                                                                                        
gAlice = (AliRun*)inkin->Get("gAlice");
cout<<"AliRun object found on file "<<gAlice<<endl;
cout<<"!!!! field ="<<gAlice->Field()->SolenoidField()<<endl;
AliKalmanTrack::SetFieldMap(gAlice->Field());
inkin->Close();                                  
 cout<<" Convconst="<<AliKalmanTrack::GetConvConst()<<endl;
/////////////////////////////////////// 
AliTPCtrack *iotrack=0; 
TFile *tf=TFile::Open("AliTPCtracks.root");
   if (!tf->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return 3;}
                                        
AliTPCPid pid;
AliTPCtrackPid pidtmp;                                                                
AliTPCtrackPid *outpid=&pidtmp;
TFile fpid("AliTPCtracksPid.root","recreate");

 Int_t nentr=0,i=0;
 for(int je=0;je<emax;je++){
   char tname[100]; sprintf(tname,"TreeT_TPC_%d",je);
   TTree *tracktree=(TTree*)tf->Get(tname);
   if (!tracktree) {cerr<<"Can't get a tree with TPC tracks !\n"; return 4;}

   TBranch *tbranch=tracktree->GetBranch("tracks");
   nentr=(Int_t)tracktree->GetEntries();
   cout<<nentr<<" "<<" tracks in track tree "<<tname<<"."<<endl;

    char tpname[100]; sprintf(tpname,"TreeT%d",je);
    TTree *ptree = new  TTree(tpname,"Tree with PID"); 
    ptree->Branch("pids","AliTPCtrackPid",&outpid,32000,1); 

   for (i=0; i<nentr; i++) {
       iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
      Double_t par[5],xx;xx=85.1951;
      iotrack->PropagateTo(xx); 
      iotrack->GetExternalParameters(xx,par);
      //cout<<" par="<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<endl;
      Float_t phi=TMath::ASin(par[2]) + iotrack->GetAlpha();
      if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
      if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
      Float_t lam=TMath::ATan(par[3]); 
      Float_t pt_1=TMath::Abs(par[4]);
      //cout<<"pt_1,lam="<<pt_1<<" "<<lam<<endl;
    if(pt_1!=0.){
      Float_t mom=1./(pt_1*TMath::Cos(lam));
      Float_t dedx=iotrack->GetdEdx();
      Int_t pcode=pid.GetPcode(dedx/50.,mom);
       pidtmp.fPcode=pcode;
       pidtmp.fSignal=dedx;
       pidtmp.fMom=mom;
       pidtmp.fPhi=phi;
       pidtmp.fLam=lam;
       pidtmp.fWpi=pid.fWpi;
       pidtmp.fWk=pid.fWk;
       pidtmp.fWp=pid.fWp;
       pidtmp.fLabel=TMath::Abs(iotrack->GetLabel());
       cout<<"tlab,dedx,mom,pcode="<<pidtmp.fLabel<<" "<<dedx<<" "<<mom<<" "<<pcode<<endl;      
       ptree->Fill();
     }//if(pt...
   }//for(i...
  cout<<"Event "<<je+1<<" ok."<<endl;      
   delete tracktree; 
 }//End for(je...)   
tf->Close();
fpid.Write();
cout<<"File AliTPCtracksPid.root written"<<endl;
return 0;
}
