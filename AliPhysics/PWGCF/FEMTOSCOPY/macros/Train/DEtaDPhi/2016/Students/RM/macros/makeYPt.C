/* by Rafal Maselek */
/* THIS IS A MACRO FOR MAKING YPt PLOTS FROM MC TRUTH OUTPUT */

void makeYPt(const char *inFile)
{
  //select particles for which you want to make YPt plots
 const int part_no = 8;
 const char* particles[part_no]={"pi+", "pi-","K+", "K-", "p+", "p-", "lambda", "lambdab"};
 const int* pid[part_no]={211, -211, 321, -321, 2212, -2212, 3122, -3122};
 const double rap_shift = 1.5344;

 int* selected_particles[part_no]={1,1,1,1,1 /* p+ */,1,1,1};

 TFile *file = new TFile(inFile);
 TFile *output_file = new TFile("YPt.root", "RECREATE");


  TH2D *h_YPt = new TH2D();
  TTree *particlesTree; //creating a tree
  file -> GetObject("events", particlesTree); //going to a certain tree in the input file
  TCut pidCut;

 //loop over all particle types
   for(int ii=0; ii<part_no; ii++)
   {
    if(selected_particles[ii]==1)
    {
     cout<<"processing: "<<particles[ii]<<endl;
     // h_YPt = new TH2D(Form("h_YPt_%s", particles[ii]), Form("Y vs Pt, %s",particles[ii]), 100, 0.0, 6.0, 100, 0.0, 3.0);
     pidCut = Form("fParticles.fPdg == %d", pid[ii]); //cut on PID
     //here histograms are added to tree
     particlesTree -> Draw("TMath::Sqrt(fParticles.fPx*fParticles.fPx+fParticles.fPy*fParticles.fPy):1.5344+0.5*TMath::Log((fParticles.fE+fParticles.fPz)/(fParticles.fE-fParticles.fPz))>>hist(100, 0.0, 6.0, 100, 0.0, 3.0)",pidCut, "goff");
     // particlesTree->Refresh();
     h_YPt = (TH2D*)gDirectory->Get("hist");
     h_YPt -> SetName(Form("y-p_{t}_%s", particles[ii]));
     h_YPt -> SetTitle(Form("Y vs Pt, %s",particles[ii]));
     output_file->cd();
     h_YPt->Write(); //writing histograms into the output file
    }
    
   }
  
}