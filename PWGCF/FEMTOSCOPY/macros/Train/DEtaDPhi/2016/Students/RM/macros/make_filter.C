/* by Rafal Maselek */
/* THIS IS A MACRO FOR MAKING FILTERS FOR ANALYSIS BY SIMPLY TAKING A RECONSTRUCTED Y VS PT PLOT
FOR RECONSTRUCTED TRACKS AND THEN DIVIDING BY MC TRUTH PLOT */

void make_filter(const char *inReco, const char *inTruth)
{
  gStyle -> SetOptStat(111111);
  
  //select particles for which you want to make filter
 const int part_no = 8;
 const char* particles[part_no]={"pi+", "pi-","K+", "K-", "p+", "p-", "lambda", "lambdab"};
 int* selected_particles[part_no]={1,1,1,1,1 /* p+ */,1,1,1};

 TFile *recoFile = new TFile(inReco);
 TFile *truthFile = new TFile(inTruth);
 TH2D *h_reco = new TH2D();
 TH2D *h_truth = new TH2D();
 TFile *output_file = new TFile("filters_10GeV.root", "RECREATE");
 TH2D *filter_hist = new TH2D();
 TH2D *filter_hist_corr = new TH2D();
 
 //loop over all particles
 for(int ii=0; ii<part_no; ii++)
 {
  if(selected_particles[ii]==1)
  {

   cout<<"processing: "<<particles[ii]<<endl;
   recoFile->cd(Form("KFTopoReconstructor/KFParticlesFinder/Particles/%s/Parameters/Signal", particles[ii]));
   h_reco=(TH2D*)gDirectory->Get("y-p_{t}");
   int reco_x_bins = h_reco ->GetNbinsX();
   int reco_y_bins = h_reco ->GetNbinsY();
   double reco_integral = h_reco->Integral();

   truthFile->cd();
   h_truth =(TH2D*)gDirectory->Get(Form("y-p_{t}_%s", particles[ii]));
   int truth_x_bins = h_truth ->GetNbinsX();
   int truth_y_bins = h_truth ->GetNbinsY();
   // double truth_integral = h_truth->Integral();
   h_truth->Scale(5000);

   //checking the bins of two reco and truth histograms
   const int filter_x_bins = reco_x_bins >= truth_x_bins ? truth_x_bins : reco_x_bins;
   const int filter_y_bins = reco_y_bins >= truth_y_bins ? truth_y_bins : reco_y_bins;
   filter_hist = new TH2D(Form("YPt_filter_%s", particles[ii]), Form("filter for %s", particles[ii]), filter_x_bins, 0.0, 6.0, filter_y_bins, 0.0, 3.0);
   filter_hist_corr = new TH2D(Form("YPt_filter_corr_%s", particles[ii]), Form("Corrected filter for %s", particles[ii]), filter_x_bins, 0.0, 6.0, filter_y_bins, 0.0, 3.0);
   // TH2D *filter_hist = new TH2D();
   double reco_bin_cont = 0.0;
   double truth_bin_cont = 0.0;

    h_truth->SetName(Form("truth_%s", particles[ii]));
    h_reco->SetName(Form("reco_%s", particles[ii]));
    output_file->cd();
    h_truth->Write();
    h_reco->Write();
   if(reco_integral != 0.0)
   {
	   //filling the filter
	   for(int xx=0; xx<filter_x_bins; xx++)
	   {
	      for(int yy=0; yy<filter_y_bins; yy++)
	      {
	        reco_bin_cont = h_reco->GetBinContent(xx, yy);
	        truth_bin_cont = h_truth->GetBinContent(xx,yy);
	        if(truth_bin_cont == 0)
	        {
	            filter_hist -> SetBinContent(xx, yy, 0.0);
	            filter_hist_corr ->SetBinContent(xx, yy, 0.0);
	        }
	        else
	        {
	            filter_hist -> SetBinContent(xx, yy, reco_bin_cont/truth_bin_cont);
	            if(reco_bin_cont/truth_bin_cont > 1.0)
	            	filter_hist_corr -> SetBinContent(xx, yy, 1.0);
	            else
	            	filter_hist_corr -> SetBinContent(xx, yy, reco_bin_cont/truth_bin_cont);
	        }
	      }
	   }
	   // filter_hist->Scale(truth_integral/reco_integral);
	  //  TCanvas *c1 = new TCanvas("c1","phi distributions",10,10,1200,750);
	  // c1 ->Divide(3,1);
	  //  c1->cd(3*ii+1);
	  //  h_reco->Draw();
	  //  c1->cd(3*ii+2);
	  //  h_truth->Draw();
	  //  c1->cd(3*ii+3);
	  //  filter_hist->Draw();
	  //  c1->Update();

	    filter_hist->Write();
	    filter_hist_corr->Write();
	}
   
  }
  
 }
    delete recoFile;
    delete truthFile;
    delete output_file;
}
