void run_all(int id=0, int end=50, int pro=0, int icut=0, int itrig=0){

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111111);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  //gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadBottomMargin(0.125);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetTitleYOffset(1.3);
  //gStyle->SetPadLeftMargin(0.1);
  cout<<"physics is always fun! "<<endl; 


  gSystem->Load("libCore");// no
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");// no
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");


  gSystem->Load("libana_sgl.so");
  ana_sgl *ana=new ana_sgl();

  char outputname[100];
  sprintf(outputname,"ntpair_rev2_%d_seg%d_%d_cut%d_trig%d_v2.root", pro, id,end, icut, itrig);


  ana->ana_init(outputname);

 //ana->set_tof_cuts(-3,3);
  //ana->enable_pair_emc_cut(0.7, 1.3);
  //ana->enable_pait_pt_cut(0.4, 20);

  if(pro==1){
    ana->set_veto_for_proton(-2, 3.5);
    ana->set_veto_for_kaon(-2, 3.5);
  }

  //0->MB, 1->SemiCent, 2-->Cent
  if(itrig==0){
    ana->select_trigger(2);
  }else if(itrig==1){
    ana->select_trigger(1);
  }

  if(icut==0){
    ana->reject_conversion(true);
    //ana->set_tpc_dedx_cuts(75, 90);
    ana->set_tpc_dedx_cuts(65, 90);
    ana->enable_pair_phiv_cut(0.6);
  }else if(icut==1){
    //ana->reject_conversion(true);
    //ana->set_tpc_dedx_cuts(80, 85);
    //ana->set_tpc_dedx_cuts(68, 75);
    ana->set_tpc_dedx_cuts(65, 90);
    ana->enable_pair_phiv_cut(0.6);
  }else if(icut==2){
    ana->reject_conversion(true);
    ana->set_tpc_dedx_cuts(75, 90);
    ana->enable_pair_phiv_cut(0.6);
  }else if(icut==3){
    ana->reject_conversion(true);
    ana->set_tpc_dedx_cuts(80, 85);
    ana->enable_pair_phiv_cut(0.6);
  }

  ana->print_cuts();

  if(pro==0){
    ifstream f("list_rev1_0.all");
    int nline = 1418;
  }else if(pro==1){
    ifstream f("list_rev2_6.all");
    int nline = 50;
  }else if(pro==2){
    ifstream f("list_rev2_7.all");
    int nline = 50;
  }


  char inputfilename[100];
  for(int i=0;i<nline;i++){
    f >> inputfilename; 
    if(i>=id && i<end){
      sprintf(inputfilename,"%s",inputfilename);
      ana->loop_a_file(inputfilename);
      cout<<" end of run : "<<i<<endl;
    }
  }

  ana->ana_end();

}

