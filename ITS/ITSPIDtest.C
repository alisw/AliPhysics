{
 AliKalmanTrack::SetConvConst(100/.299792458/.2);
   cout<<" Starting..."<<endl;                                                     
if(gClassTable->GetID("AliRun")<0){
    gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
    loadlibs();                                                                 
}

if(gROOT->IsBatch()){
    cout<<"Start tracking..."<<endl;
    gSystem->Exec("make all");
    cout<<"AliITSv2PID.root written."<<endl;
}else{

   fpid = new TFile("pidhit.root","RECREATE");

   gROOT->Macro("$ALICE_ROOT/ITS/load_particles.C");
   AliITSPid *pid = new AliITSPid(npart);          
   gROOT->LoadMacro("$ALICE_ROOT/ITS/dEdXgeant.C"); 
   gROOT->LoadMacro("$ALICE_ROOT/ITS/dedxanal.C");    
//----------------------------------------------------
NStat=pid.trs->GetEntries();      
//----------------------------------------------------                          
    TControlBar menu("vertical","PID menu",920,5);
    
    menu.AddButton("dEdX.C","pid->Reset();totpid=0;dEdXyy(0,0,pid);pid->Tab();"," Create new PID table ");               
    menu.AddButton("Save TAB",
	"pid->Tab();fpid->cd();pid->Write();fpid->Close();"," ");
   menu.AddButton("Load TAB","loadpid();"," ");    

    menu.AddButton("EFFALL","effall(); "," Efficiency if PID  ");
    menu.AddButton("dEdX spectra","qhisall(); "," dEdX for PI,K and P  ");
    menu.AddButton("dEdX-P plot","dedxhis(0); "," dEdX-P plot for PI,K,P  ");
    menu.AddButton("dEdX-P pions","dedxhis(211); "," dEdX-P plot for PI  ");
    menu.AddButton("dEdX-P kaons","dedxhis(321); "," dEdX-P plot for K  ");
    menu.AddButton("dEdX-P elect","dedxhis(11); "," dEdX-P plot for e+  ");
    menu.AddButton("dEdX-P prot ","dedxhis(2212); "," dEdX-P plot for P  ");
    
    menu.AddButton("Fit Kaons","fitkall(); "," Gaus Fit for Kaons   ");
    menu.AddButton("Fit Pions","fitpiall(); "," Gaus Fit for Kaons   ");
    menu.AddButton("Fit Protons","fitpall(); "," Gaus Fit for Protons   ");
    menu.AddButton("New cuts","newcuts(); "," Corrected cuts for PID object   ");
    menu.AddButton("pcode","pcode(); "," ...  ");              
    menu.AddButton("signal (mip)","signal(); "," ...  ");    
    menu.AddButton("pmom (MeV)","pmom(); "," ...  ");        
    menu.AddButton("tracks","tracks(); "," Track number histogram  ");   
    menu.AddButton("1 track","track(); "," Print next track  ");                                  
    menu.AddButton("test module","dEdXxx(0,0,pid,1); "," ...  ");
    menu.AddButton("fill tab test","filltab(); ","Fill track table with test data ");
    menu.AddButton("fill tab_tr","filltab_tracks(); ","Fill track table with reconstr. tracks ");

    menu.AddButton("Config.C","gSystem->Exec(\"make conf\");","Edit Config.C");
    menu.AddButton("Do tracking",
	"gSystem->Exec(\"make all\");pid->Reset();totpid=0;filltab_tracks();dedxhis(0)",
	"Start tracking");
    menu.AddButton("Exit","quit();","Quit");
  menu.Show();

}//if batch

}








