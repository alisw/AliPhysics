void SSD_ntuple()
{

TFile *f = new TFile("SSD_his.root");
   
//gStyle->SetOptStat(1111111);
//gStyle->SetOptLogy();
//TCanvas *c1 = new TCanvas("c1","SPD clusters",400,10,600,700);
TCanvas *c2 = new TCanvas("c2","SPD clusters",400,10,600,700);
//c1->Divide(2,2);
c2->Divide(2,2);

/////////////////////////  Ntuple analysis ///////////////////////////////

// ntuple is created inside the hit loop for the hits corresponding to the
// recpoint;

// ntuple1 is created after a finish of the hit loop if one or more hits
// correspond to the recpoint;

// ntuple2 is created befor the hit loop for all recpoints;

// -----------------------------------------------------------------------
// lay       - number of ITS layer;
// lad       - number of ITS ladder;
// det       - number of ITS detector;
// nxP/N        - cluster size in the r*phi(x) direction for P/N sides;
// hitprim   - primary particle(hit) flag ( = 1 for primery particle);     
// x         - x local coordinate in mm; 
// z         - z local coordinate in mm; 
// dx        - difference of hit(mediate) and reconstructed (from cluster)
//             coordinates in r*phi(x) direction in microns;
// dz        - difference of hit(mediate) and reconstructed (from cluster)
//             coordinates in z direction in microns;
// noverlaps - number of particles overlapping in one cluster; 
// noverprim - number of primary particles overlapping in one cluster;
// qclP/N    - cluster signals in ADC normalized to the path length in Si  
// qrec      - recpoint signal (maximum from qclP and qclN)  
// pmod      - particle momentum at the vertex in MeV/c
// partcode  - particle code
// -------------------------------------------------------------------------



c2->cd(1);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay == 5 && hitprim == 1&&abs(dx)<200");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay == 5 && hitprim == 1&&abs(dz)<5000");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple->SetFillColor(42);
      ntuple->Draw("dx","lay == 6 && hitprim == 1&&abs(dx)<200");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple->SetFillColor(46);
      ntuple->Draw("dz","lay == 6 && hitprim == 1&&abs(dz)<5000");


      /*
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("nxP","lay == 5&&noverprim>=0");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("nxN","lay == 5&&noverprim>=0");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("nxP","lay == 6 && noverprim>=0");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("nxN","lay == 6 && noverprim>=0");
      */

/*
c2->cd(1);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("qclP","lay == 5&&noverprim>=0");
c2->cd(2);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(42);
      ntuple1->Draw("qclN","lay == 5&&noverprim>=0");
c2->cd(3);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("qclP","lay == 6&&noverprim>=0");
c2->cd(4);
gPad->SetFillColor(33);
      ntuple1->SetFillColor(46);
      ntuple1->Draw("qclN","lay == 6&&noverprim>=0");
*/


/////////////////////   Histogramm/ntuple  analysis  ////////////////////////

      /*
c2->cd(1);
gPad->SetFillColor(33);
      adcPadcN5all->SetFillColor(42);
      adcPadcN5all->Draw();
c2->cd(2);
gPad->SetFillColor(33);
      adcPadcN6all->SetFillColor(46);
      adcPadcN6all->Draw();
c2->cd(3);
gPad->SetFillColor(33);
      adcPadcN5cut->SetFillColor(42);
      adcPadcN5cut->Draw();
c2->cd(4);
gPad->SetFillColor(33);
      adcPadcN6cut->SetFillColor(46);
      adcPadcN6cut->Draw();
      */
      

      /*
c2->Draw();   
c2->Print("ssd_res.ps");
      */

      /*                        
c2->Draw();   
c2->Print("spd_clsize.ps");
      */    

}









