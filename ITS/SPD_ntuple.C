void SPD_ntuple()
{

TFile *f = new TFile("SPD_his.root");
//TH1F *Pix = (TH1F*)f->Get("Pix");
   
//gStyle->SetOptStat(1111111);
gStyle->SetOptLogy();
//TCanvas *c1 = new TCanvas("c1","SPD clusters",400,10,600,700);
TCanvas *c2 = new TCanvas("c2","SPD clusters",400,10,600,700);
//c1->Divide(2,2);
c2->Divide(2,2);
//c2->SetLogy();

/////////////////////////  Ntuple analysis ///////////////////////////////

// ntuple is created inside the hit loop for the hits in the cluster region;

// ntuple1 is created after a finish of the hit loop if one or more hits
// are in the cluster region;

// ntuple2 is created befor the hit loop for all clusters;

// -----------------------------------------------------------------------
// lay       - number of ITS layer;
// nx        - cluster size in the r*phi(x) direction;
// nz        - cluster size in the z direction;
// hitprim   - primary particle(hit) flag ( = 1 for primery particle);     
// dx        - difference of hit(mediate) and reconstructed (from cluster)
//             coordinates in r*phi(x) direction;
// dz        - difference of hit(mediate) and reconstructed (from cluster)
//             coordinates in z direction;
// noverlaps - number of particles overlapping in one cluster finded in this
//             macros in the cluster region;
// ntrover   - number of particles overlapping in one cluster finded in the
//             AliITSClusterFinderSPD class; 
// noverprim - number of primary particles overlapping in one cluster;
// qcl       - cluster charge in electron number  
// -------------------------------------------------------------------------

/*
c2->cd(1);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay == 1 && hitprim == 1");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay == 1 && hitprim == 1");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay == 2 && hitprim == 1");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay == 2 && hitprim == 1");
*/

/*      
c2->cd(1);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay == 1 && hitprim == 1");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dx","lay == 1 && hitprim == 1 && nx < 15");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay == 1 && hitprim == 1 && nx < 5");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dx","lay == 1 && hitprim == 1 && nx < 4");
*/
    

      /*                  
c2->cd(1);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay==1&&hitprim==1&&dx>-150&&dx<150");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay==1&&hitprim==1&&dz>-500&&dz<500");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay==2&&hitprim==1&&dx>-150&&dx<150");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay==2&&hitprim==1&&dz>-500&&dz<500");
      */

/*
c2->cd(1);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay==1&&hitprim==1&&nx<4");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay==1&&hitprim==1&&nz<4");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay==2&&hitprim==1&&nx<4");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay==2&&hitprim==1&&nz<4");
*/      

                              
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("noverlaps","lay == 1");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("noverlaps","lay == 2");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("ntrover","lay == 1");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("ntrover","lay == 2");

/*        
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("noverprim","lay == 1");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("noverprim","lay == 2");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("noverprim","lay == 1 && noverprim > 1");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("noverprim","lay == 2 && noverprim > 1");
*/       

/*
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("qcl","lay == 1&&noverprim>0");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("qcl","lay == 2&&noverprim>0");
*/

/*      
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("nx","lay == 1&&noverprim>0");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("nx","lay == 2&&noverprim>0");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("nz","lay == 1 && noverprim>0");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("nz","lay == 2 && noverprim>0");
*/
            
/*
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("dx","lay == 1&&noverprim>0");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("dx","lay == 2&&noverprim>0");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("dz","lay == 1 && noverprim>0");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("dz","lay == 2 && noverprim>0");
*/

      /*            
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("dx","lay == 1");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("dx","lay == 2");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("dz","lay == 1");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("dz","lay == 2");
      */

      /*                
c2->cd(1);
gPad->SetFillColor(33);
      ntuple2->SetFillColor(42);
      ntuple2->Draw("nx","lay == 1");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple2->SetFillColor(46);
      ntuple2->Draw("nz","lay == 1");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple2->SetFillColor(42);
      ntuple2->Draw("nx","lay == 2");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple2->SetFillColor(46);
      ntuple2->Draw("nz","lay == 2");
      */
      
/////////////////////   Histogramm analysis  ////////////////////////

      /*            
c2->cd(1);
gPad->SetFillColor(33);
      Xres1->SetFillColor(42);
      Xres1->Draw();
c2->cd(2);
gPad->SetFillColor(33);
      Zres1->SetFillColor(42);
      Zres1->Draw();
c2->cd(3);
gPad->SetFillColor(33);
      Xres2->SetFillColor(46);
      Xres2->Draw();
c2->cd(4);
gPad->SetFillColor(33);
      Zres2->SetFillColor(46);
      Zres2->Draw();
      */
      /*
c2->cd(1);
gPad->SetFillColor(33);
      Nzpix1->SetFillColor(42);
      Nzpix1->Draw();
c2->cd(2);
gPad->SetFillColor(33);
      Nxpix1->SetFillColor(42);
      Nxpix1->Draw();
c2->cd(3);
gPad->SetFillColor(33);
      Nzpix2->SetFillColor(46);
      Nzpix2->Draw();
c2->cd(4);
gPad->SetFillColor(33);
      Nxpix2->SetFillColor(46);
      Nxpix2->Draw();
      */                         
      
      /*                          
c2->cd(1);
gPad->SetFillColor(33);
      Zpix1->SetFillColor(42);
      Zpix1->Draw();
c2->cd(2);
gPad->SetFillColor(33);
      Xpix1->SetFillColor(46);
      Xpix1->Draw();
c2->cd(3);
gPad->SetFillColor(33);
      Zpix2->SetFillColor(42);
      Zpix2->Draw();
c2->cd(4);
gPad->SetFillColor(33);
      Xpix2->SetFillColor(46);
      Xpix2->Draw();
      */
     
      /*                       
c2->Draw();   
c2->Print("spd_res.ps");
      */
      /*                        
c2->Draw();   
c2->Print("spd_clsize.ps");
      */    

}








