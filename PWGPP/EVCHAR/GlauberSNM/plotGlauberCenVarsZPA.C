///==========================================================================
///
///    macro to plot centrality bin values 
///==========================================================================
///
#include <cstdlib>

const Int_t nbins;

double centPercent[]={5.,10.,20.,40.,60.,80.,99.9999999}; 
//double centPercent[]={10.,20.,40.,60.,80.,99.9999999}; 


nbins = sizeof(centPercent)/sizeof(double);

TArrayI* binUp = new TArrayI(nbins);
TArrayF* Multbin = new TArrayF(nbins);

void plotGlauberCenVarsZPA()    
{
 TGraphErrors *gNpart=new TGraphErrors(0);
 gNpart->SetName("gNpart"); 
 TGraphErrors *gNcoll=new TGraphErrors(0);
 gNcoll->SetName("gNcoll"); 
 TGraphErrors *gtAA=new TGraphErrors(0);
 gtAA->SetName("gtAA"); 


  TFile *f=new TFile("ZPA_ntuple_188359.root");
  
  TNtuple *nt=(TNtuple*)f->Get("gnt");
  
  Float_t B;
  Float_t Npart;
  Float_t Ncoll;
  Float_t tAA;
  Float_t ntot;
  
  nt->SetBranchAddress("B",&B);
  nt->SetBranchAddress("Npart",&Npart);
  nt->SetBranchAddress("Ncoll",&Ncoll);
  nt->SetBranchAddress("tAA",&tAA);
  nt->SetBranchAddress("ntot",&ntot);


  Int_t totev=nt->GetEntries();
 
 
  TH1F*hB=new TH1F("hB","hB",28000,0.,14.);
  TH1F*hNpart=new TH1F("hNpart","hNpart",420,-0.5,419.5);
  TH1F*hNcoll=new TH1F("hNcoll","hNcoll",3500,0,35);
  TH1F*htAA=new TH1F("htAA","htAA",100000,-0.5,39.5);
  TH1F*hntot=new TH1F("hntot","hntot",100000,-0.5,39.5);
 
  
  for (Int_t iEvent = 0; iEvent < totev; iEvent++) {
    nt->GetEvent(iEvent); 
    
    Float_t nu = (Float_t) ntot; 
    nu = gRandom->Gaus(nu,0.5); 

    if(nu>0) {
  
    hNpart->Fill(Npart);
    hNcoll->Fill(Ncoll);
    hB->Fill(B);
    htAA->Fill(tAA);
    hntot->Fill(nu);
    }
  } 

   new TCanvas();
   hNcoll->Draw();

   //---------------------------------------------------
   getCentrality(hntot);
   //---------------------------------------------------


 TH1F* hnpartcutb[nbins];
 char histtitp[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitp,"npartcutb%d",i);
   hnpartcutb[i] = new TH1F(histtitp,histtitp,10000,-0.5,9999.5);
   hnpartcutb[i]->SetLineWidth(1);
   hnpartcutb[i]->SetStats(1);

 }

 TH1F* hncollcutb[nbins];
 char histtitc[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitc,"ncollcutb%d",i);
   hncollcutb[i] = new TH1F(histtitc,histtitc,10000,-0.5,9999.5);
   hncollcutb[i]->SetLineWidth(1);
   hncollcutb[i]->SetStats(1);

 }

 TH1F* htaacutb[nbins];
 char histtitt[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitt,"taacutb%d",i);
   htaacutb[i] = new TH1F(histtitt,histtitt,100000,-0.5,39.5);
   htaacutb[i]->SetLineWidth(1);
   htaacutb[i]->SetStats(1);

 }

 TH1F* hbcutb[nbins];
 char histtitt[100];
 for (Int_t i=0; i<binUp->GetSize(); i++) {
   sprintf(histtitt,"bcutb%d",i);
   hbcutb[i] = new TH1F(histtitt,histtitt,100000,-0.5,39.5);
   hbcutb[i]->SetLineWidth(1);
   hbcutb[i]->SetStats(1);

 }



 for (Int_t iEvent = 0; iEvent < totev; iEvent++) {
   nt->GetEvent(iEvent);
   for (int ibin=0; ibin<nbins; ibin++) {
     //     if (B<Multbin->At(ibin)) {

     Float_t nu = (Float_t) Ncoll; 
     nu = nu+1.*gRandom->Rndm();
  
     if (nu>Multbin->At(ibin)) {
       hnpartcutb[ibin]->Fill(Npart); 
       hncollcutb[ibin]->Fill(Ncoll); 
       htaacutb[ibin]->Fill(tAA); 
       hbcutb[ibin]->Fill(B); 
       break;
     } 
   } 
 } 
 

   


 // superimpose cut histograms
// new TCanvas();
// hnpartcutb[0]->Draw("");
 for (Int_t icentr=0; icentr<nbins;icentr++) {
   hnpartcutb[icentr]->SetLineColor(icentr);
   hnpartcutb[icentr]->Draw("same");
 } 

 for (Int_t icentr=0; icentr<nbins;icentr++) {
   new TCanvas();
   hnpartcutb[icentr]->SetLineColor(icentr);
   hnpartcutb[icentr]->SetStats(1); 
   hnpartcutb[icentr]->Draw("");
   gNpart->SetPoint(icentr,Float_t(icentr),hnpartcutb[icentr]->GetMean());
   gNpart->SetPointError(icentr,0,hnpartcutb[icentr]->GetRMS());
 }
 cout << endl;

 for (Int_t icentr=0; icentr<nbins;icentr++) {
   new TCanvas();
   hncollcutb[icentr]->SetLineColor(icentr);
   hncollcutb[icentr]->SetStats(1); 
   hncollcutb[icentr]->Draw("");
   gNcoll->SetPoint(icentr,Float_t(icentr),hncollcutb[icentr]->GetMean());
   gNcoll->SetPointError(icentr,0,hncollcutb[icentr]->GetRMS());
 }
 cout << endl;

 for (Int_t icentr=0; icentr<nbins;icentr++) {
   new TCanvas();
   htaacutb[icentr]->SetLineColor(icentr);
   htaacutb[icentr]->SetStats(1); 
   htaacutb[icentr]->Draw("");
   gtAA->SetPoint(icentr,Float_t(icentr),htaacutb[icentr]->GetMean());
   gtAA->SetPointError(icentr,0,htaacutb[icentr]->GetRMS());
 }

 for (Int_t icentr=0; icentr<nbins;icentr++) 
   cout<<icentr<<" | "<<setprecision(3)<<centPercent[icentr]<<" | "
       <<Multbin->At(icentr-1)<<" | "<<Multbin->At(icentr)<<" | "
       <<hnpartcutb[icentr]->GetMean()<<" | " <<hnpartcutb[icentr]->GetRMS()<< " | " 
       <<hncollcutb[icentr]->GetMean()<<" | " <<hncollcutb[icentr]->GetRMS()<< " | " 
       <<hbcutb[icentr]->GetMean()    <<" | " <<hbcutb[icentr]->GetRMS()    << " | "
       <<htaacutb[icentr]->GetMean()  <<" | " <<htaacutb[icentr]->GetRMS()  << " | "<<endl;


 //TString suffixhisto=Form("/home/atoia/GlauberNtuple/62/GlauberMC_PbPb_histoB_sigma%d_mind%d_r%d_a%d_%d.root",mysignn,mymind,myr,mya,nbins);
 //TString suffixhisto=Form("/home/atoia/GlauberNtuple/GlauberMC_PbPb_histoB_sigma%d_mind%d_r%d_a%d_%d.root",mysignn,mymind,myr,mya,nbins);
 //const Char_t* filehistoname=suffixhisto.Data();
 //TFile*filefinal=new TFile(filehistoname,"recreate");
 // --FOR TEST --
 //TFile*filefinal=new TFile("GlauberInfo.root","recreate");
 
 // gNpart->Write();
 // gNcoll->Write();
 // gtAA->Write();
 
 // hNpart->Write();
 // hNcoll->Write();
 // hB    ->Write();
 // htAA  ->Write();
  
 // for (int i=0; i<nbins; i++) {
 //   hnpartcutb[i]->Write();
 //   hncollcutb[i]->Write();
 //   htaacutb[i]  ->Write();
 // }
 
}


void getCentrality(TH1 *histNch, Float_t ff=1.0)
// histNch - histo of multiplicity distribution (better with binsize=1)x
// ff fraction of accepted events. All losses are assumed to occur in most
// peripheral bin
{

 //double sum= histNch->GetEntries() - histNch->GetBinContent(1);
 double sum= histNch->Integral(); 
 int nbinsx=histNch->GetNbinsX();
 double frac=0.;
 int ic=0;
 //for (int ib=1;ib<=nbinsx;ib++){
 for (int ib=nbinsx;ib>0;ib--){
   frac += histNch->GetBinContent(ib)/sum*100.*ff;
   if(frac > centPercent[ic]){
     binUp->SetAt(ib,ic);
     Multbin->SetAt(histNch->GetBinCenter(ib),ic);
     cout<<" centrality="<<centPercent[ic]<<"   ncoll <="<< histNch->GetBinCenter(ib) <<endl;
     //cout<<" centrality="<<centPercent[ic]<<" impact parameter <="<< histNch->GetBinCenter(ib) <<endl;
     ic++;
   }
   if(ic==nbins) break;
 }
 
 printf(" \n float binUp[%i] = {",nbins);  
 // cout <<" \n float multCent[nbins] = {";

 for (int ic=nbins-1; ic>-1; ic--){
   cout<< binUp->At(ic);
   if (ic!=0) cout<<", ";
 }
 cout<<"};\n"<<endl;


 printf(" \n float multCent[%i] = {",nbins);  
 // cout <<" \n float multCent[nbins] = {";

 for (int ic=nbins-1; ic>-1; ic--){
   cout<< histNch->GetBinCenter(binUp->At(ic));
   if (ic!=0) cout<<", ";
 }
 cout<<"};\n"<<endl;
}
