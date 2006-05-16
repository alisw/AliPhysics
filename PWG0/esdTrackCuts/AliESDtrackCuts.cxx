#include "AliESDtrackCuts.h"

#include <Riostream.h>

//____________________________________________________________________
ClassImp(AliESDtrackCuts);

//____________________________________________________________________
AliESDtrackCuts::AliESDtrackCuts() {

  //##############################################################################
  // setting default cuts

  SetMinNClustersTPC();
  SetMinNClustersITS();	    
  SetMaxChi2PerClusterTPC();
  SetMaxChi2PerClusterITS();  				    
  SetMaxCovDiagonalElements();  				    
  SetRequireTPCRefit();
  SetRequireITSRefit();
  SetAcceptKingDaughters();
  SetMinNsigmaToVertex();
  SetRequireSigmaToVertex();
  SetPRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();
  SetEtaRange();
  SetRapRange();

  SetHistogramsOn();

  // set the cut names
  fCutNames[0]  = "require TPC refit";
  fCutNames[1]  = "require ITS refit";  
  fCutNames[2]  = "n clusters TPC";
  fCutNames[3]  = "n clusters ITS";
  fCutNames[4]  = "#Chi^{2}/clusters TPC";
  fCutNames[5]  = "#Chi^{2}/clusters ITS";
  fCutNames[6]  = "cov 11";
  fCutNames[7]  = "cov 22";
  fCutNames[8]  = "cov 33";
  fCutNames[9]  = "cov 44";
  fCutNames[10] = "cov 55";
  fCutNames[11] = "trk-to-vtx";
  fCutNames[12] = "trk-to-vtx failed";
  fCutNames[13] = "kink daughters";

  fCutNames[14] = "p";
  fCutNames[15] = "p_{T}";
  fCutNames[16] = "p_{x}";
  fCutNames[17] = "p_{y}";
  fCutNames[18] = "p_{z}";
  fCutNames[19] = "y";
  fCutNames[20] = "eta";

}

//____________________________________________________________________
Bool_t 
AliESDtrackCuts::AcceptTrack(AliESDtrack* esdTrack) {
  // 
  // figure out if the tracks survives all the track cuts defined
  //

  UInt_t status = esdTrack->GetStatus();
  
  // getting quality parameters from the ESD track
  Int_t nClustersITS = esdTrack->GetITSclusters(fIdxInt);
  Int_t nClustersTPC = esdTrack->GetTPCclusters(fIdxInt);
  
  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = esdTrack->GetITSchi2()/Float_t(nClustersITS);
  if (nClustersTPC!=0)
    chi2PerClusterTPC = esdTrack->GetTPCchi2()/Float_t(nClustersTPC);  

  Double_t extCov[15];
  esdTrack->GetExternalCovariance(extCov);  

  // getting the track to vertex parameters
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);    
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution zero!");
    bCov[0]=0; bCov[1]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //
  // this is not correct - it will not give n sigma!!!
  // 
  Float_t nSigmaToVertex = -1;
  if (bRes[0]!=0 && bRes[1]!=0)
    nSigmaToVertex = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));  

  // getting the kinematic variables of the track 
  // (assuming the mass is known)
  Double_t p[3];
  esdTrack->GetPxPyPz(p);
  Float_t momentum = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2) + TMath::Power(p[2],2));
  Float_t pt       = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2));
  Float_t energy   = TMath::Sqrt(TMath::Power(esdTrack->GetMass(),2) + TMath::Power(momentum,2));

  //y-eta related calculations
  Float_t eta = -100.;
  Float_t y   = -100.;
  if((momentum != TMath::Abs(p[2]))&&(momentum != 0))
    eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  if((energy != TMath::Abs(p[2]))&&(momentum != 0))
    y = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));

  
  //########################################################################
  // cut the track?
  
  Bool_t cuts[fNCuts];
  for (Int_t i=0; i<fNCuts; i++) cuts[i]=kFALSE;
  
  // track quality cuts
  if (fCut_RequireTPCRefit && (status&AliESDtrack::kTPCrefit)==0)
    cuts[0]=kTRUE;
  if (fCut_RequireITSRefit && (status&AliESDtrack::kITSrefit)==0)
    cuts[1]=kTRUE;
  if (nClustersTPC<fCut_MinNClusterTPC) 
    cuts[2]=kTRUE;
  if (nClustersITS<fCut_MinNClusterITS) 
    cuts[3]=kTRUE;
  if (chi2PerClusterTPC>fCut_MaxChi2PerClusterTPC) 
    cuts[4]=kTRUE; 
  if (chi2PerClusterITS>fCut_MaxChi2PerClusterITS) 
    cuts[5]=kTRUE;
  if (extCov[0]  > fCut_MaxC11) 
    cuts[6]=kTRUE;  
  if (extCov[2]  > fCut_MaxC22) 
    cuts[7]=kTRUE;  
  if (extCov[5]  > fCut_MaxC33) 
    cuts[8]=kTRUE;  
  if (extCov[9]  > fCut_MaxC44) 
    cuts[9]=kTRUE;  
  if (extCov[14]  > fCut_MaxC55) 
    cuts[10]=kTRUE;  
  if (nSigmaToVertex > fCut_NsigmaToVertex) 
    cuts[11] = kTRUE;
  // if n sigma could not be calculated
  if (nSigmaToVertex<0 && fCut_SigmaToVertexRequired)   
    cuts[12]=kTRUE;
  if (!fCut_AcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0) 
    cuts[13]=kTRUE;
  // track kinematics cut
  if((momentum < fPMin) || (momentum > fPMax)) 
    cuts[14]=kTRUE;
  if((pt < fPtMin) || (pt > fPtMax)) 
    cuts[15] = kTRUE;
  if((p[0] < fPxMin) || (p[0] > fPxMax)) 
    cuts[16] = kTRUE;
  if((p[1] < fPyMin) || (p[1] > fPyMax)) 
    cuts[17] = kTRUE;
  if((p[2] < fPzMin) || (p[2] > fPzMax)) 
    cuts[18] = kTRUE;
  if((eta < fEtaMin) || (eta > fEtaMax)) 
    cuts[19] = kTRUE;
  if((y < fRapMin) || (y > fRapMax)) 
    cuts[20] = kTRUE;

  Bool_t cut=kFALSE;
  for (Int_t i=0; i<fNCuts; i++) 
    if (cuts[i]) cut = kTRUE;
  
  //########################################################################
  // filling histograms
  if (fHistogramsOn) {
    hCutStatistics->Fill(hCutStatistics->GetBinCenter(hCutStatistics->GetXaxis()->FindBin("n tracks")));
    
    if (cut)
      hCutStatistics->Fill(hCutStatistics->GetBinCenter(hCutStatistics->GetXaxis()->FindBin("n cut tracks")));
    
    for (Int_t i=0; i<fNCuts; i++) {
      if (cuts[i])
 	hCutStatistics->Fill(hCutStatistics->GetBinCenter(hCutStatistics->GetXaxis()->FindBin(fCutNames[i])));
      
      for (Int_t j=i; j<fNCuts; j++) {
 	if (cuts[i] && cuts[j]) {
 	  Float_t x = hCutCorrelation->GetXaxis()->GetBinCenter(hCutCorrelation->GetXaxis()->FindBin(fCutNames[i]));
 	  Float_t y = hCutCorrelation->GetYaxis()->GetBinCenter(hCutCorrelation->GetYaxis()->FindBin(fCutNames[j]));
 	  hCutCorrelation->Fill(x,y);
 	}
      }
    }
    

    hNClustersITS[0]->Fill(nClustersITS);        
    hNClustersTPC[0]->Fill(nClustersTPC);        
    hChi2PerClusterITS[0]->Fill(chi2PerClusterITS);
    hChi2PerClusterTPC[0]->Fill(chi2PerClusterTPC);   
    
    hC11[0]->Fill(extCov[0]);                 
    hC22[0]->Fill(extCov[2]);                 
    hC33[0]->Fill(extCov[5]);                 
    hC44[0]->Fill(extCov[9]);                                  
    hC55[0]->Fill(extCov[14]);                                  
    
    hDZ[0]->Fill(b[1]);     
    hDXY[0]->Fill(b[0]);    
    hDXYvsDZ[0]->Fill(b[1],b[0]);

    if (bRes[0]!=0 && bRes[1]!=0) {
      hDZNormalized[0]->Fill(b[1]/bRes[1]);     
      hDXYNormalized[0]->Fill(b[0]/bRes[0]);    
      hDXYvsDZNormalized[0]->Fill(b[1]/bRes[1], b[0]/bRes[0]);
    }
  }

  //########################################################################  
  // cut the track!
  if (cut) return kFALSE;

  //########################################################################  
  // filling histograms after cut
  if (fHistogramsOn) {
    hNClustersITS[1]->Fill(nClustersITS);        
    hNClustersTPC[1]->Fill(nClustersTPC);        
    hChi2PerClusterITS[1]->Fill(chi2PerClusterITS);
    hChi2PerClusterTPC[1]->Fill(chi2PerClusterTPC);   
    
    hC11[1]->Fill(extCov[0]);                 
    hC22[1]->Fill(extCov[2]);                 
    hC33[1]->Fill(extCov[5]);                 
    hC44[1]->Fill(extCov[9]);                                  
    hC55[1]->Fill(extCov[14]);                                  
    
    hDZ[1]->Fill(b[1]);     
    hDXY[1]->Fill(b[0]);    
    hDXYvsDZ[1]->Fill(b[1],b[0]);

    hDZNormalized[1]->Fill(b[1]/bRes[1]);     
    hDXYNormalized[1]->Fill(b[0]/bRes[0]);    
    hDXYvsDZNormalized[1]->Fill(b[1]/bRes[1], b[0]/bRes[0]);
  }
  
  return kTRUE;
}

//____________________________________________________________________
TObjArray*
AliESDtrackCuts::GetAcceptedTracks(AliESD* esd) {
  
  // returns an array of all tracks that pass the cuts
  fAcceptedTracks->Clear();
  
  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = esd->GetTrack(iTrack);
    
    if(AcceptTrack(track)) fAcceptedTracks->Add(track);
  }

  return fAcceptedTracks;
}

//____________________________________________________________________
void 
AliESDtrackCuts::DefineHistograms(Int_t color) {

  fHistogramsOn=kTRUE;

  //###################################################################################
  // defining histograms

  hCutStatistics = new TH1F("cut_statistics","cut statistics",fNCuts+4,-0.5,fNCuts+3.5);

  hCutStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
  hCutStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");

  hCutCorrelation = new TH2F("cut_correlation","cut correlation",fNCuts,-0.5,fNCuts-0.5,fNCuts,-0.5,fNCuts-0.5);;
  
  for (Int_t i=0; i<fNCuts; i++) {
    hCutStatistics->GetXaxis()->SetBinLabel(i+4,fCutNames[i]);
    hCutCorrelation->GetXaxis()->SetBinLabel(i+1,fCutNames[i]);
    hCutCorrelation->GetYaxis()->SetBinLabel(i+1,fCutNames[i]);
  } 

  hCutStatistics  ->SetLineColor(color);
  hCutCorrelation ->SetLineColor(color);
  hCutStatistics  ->SetLineWidth(2);
  hCutCorrelation ->SetLineWidth(2);


  hNClustersITS        = new TH1F*[2];
  hNClustersTPC        = new TH1F*[2];
  hChi2PerClusterITS   = new TH1F*[2];
  hChi2PerClusterTPC   = new TH1F*[2];
  		       
  hC11                 = new TH1F*[2];
  hC22                 = new TH1F*[2];
  hC33                 = new TH1F*[2];
  hC44                 = new TH1F*[2];
  hC55                 = new TH1F*[2];
  
  hDXY                 = new TH1F*[2];
  hDZ		       = new TH1F*[2];
  hDXYvsDZ	       = new TH2F*[2];

  hDXYNormalized       = new TH1F*[2];
  hDZNormalized        = new TH1F*[2];
  hDXYvsDZNormalized   = new TH2F*[2];


  Char_t str[256];
  for (Int_t i=0; i<2; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");

    hNClustersITS[i]        = new TH1F(Form("nClustersITS%s",str),"",8,-0.5,7.5);
    hNClustersTPC[i]        = new TH1F(Form("nClustersTPC%s",str),"",165,-0.5,164.5);
    hChi2PerClusterITS[i]   = new TH1F(Form("chi2PerClusterITS%s",str),"",500,0,10);
    hChi2PerClusterTPC[i]   = new TH1F(Form("chi2PerClusterTPC%s",str),"",500,0,10);
    
    hC11[i]                 = new TH1F(Form("covMatrixDiagonal11%s",str),"",1000,0,5);
    hC22[i]                 = new TH1F(Form("covMatrixDiagonal22%s",str),"",1000,0,5);
    hC33[i]                 = new TH1F(Form("covMatrixDiagonal33%s",str),"",1000,0,0.5);
    hC44[i]                 = new TH1F(Form("covMatrixDiagonal44%s",str),"",1000,0,5);
    hC55[i]                 = new TH1F(Form("covMatrixDiagonal55%s",str),"",1000,0,5);
    
    hDXY[i]                 = new TH1F(Form("dXY%s",str),"",500,-10,10);
    hDZ[i]                  = new TH1F(Form("dZ%s",str),"",500,-10,10);
    hDXYvsDZ[i]             = new TH2F(Form("dXYvsDZ%s",str),"",200,-10,10,200,-10,10);

    hDXYNormalized[i]       = new TH1F(Form("dXYNormalized%s",str),"",500,-10,10);
    hDZNormalized[i]        = new TH1F(Form("dZNormalized%s",str),"",500,-10,10);
    hDXYvsDZNormalized[i]   = new TH2F(Form("dXYvsDZNormalized%s",str),"",200,-10,10,200,-10,10);


    hNClustersITS[i]        ->SetXTitle("n ITS clusters");  
    hNClustersTPC[i]        ->SetXTitle("n TPC clusters"); 
    hChi2PerClusterITS[i]   ->SetXTitle("#Chi^{2} per ITS cluster"); 
    hChi2PerClusterTPC[i]   ->SetXTitle("#Chi^{2} per TPC cluster"); 
    
    hC11[i]                 ->SetXTitle("cov 11 : #sigma_{y}^{2} [cm^{2}]"); 
    hC22[i]                 ->SetXTitle("cov 22 : #sigma_{z}^{2} [cm^{2}]"); 
    hC33[i]                 ->SetXTitle("cov 33 : #sigma_{sin(#phi)}^{2}"); 
    hC44[i]                 ->SetXTitle("cov 44 : #sigma_{tan(#theta_{dip})}^{2}"); 
    hC55[i]                 ->SetXTitle("cov 55 : #sigma_{1/p_{T}}^{2} [(c/GeV)^2]"); 
   
    hDXY[i]                 ->SetXTitle("transverse impact parameter"); 
    hDZ[i]                  ->SetXTitle("longitudinal impact parameter"); 
    hDXYvsDZ[i]             ->SetXTitle("longitudinal impact parameter"); 
    hDXYvsDZ[i]             ->SetYTitle("transverse impact parameter"); 

    hDXYNormalized[i]       ->SetXTitle("normalized trans impact par"); 
    hDZNormalized[i]        ->SetXTitle("normalized long impact par"); 
    hDXYvsDZNormalized[i]   ->SetXTitle("normalized long impact par"); 
    hDXYvsDZNormalized[i]   ->SetYTitle("normalized trans impact par"); 

    hNClustersITS[i]        ->SetLineColor(color);   hNClustersITS[i]        ->SetLineWidth(2);
    hNClustersTPC[i]        ->SetLineColor(color);   hNClustersTPC[i]        ->SetLineWidth(2);
    hChi2PerClusterITS[i]   ->SetLineColor(color);   hChi2PerClusterITS[i]   ->SetLineWidth(2);
    hChi2PerClusterTPC[i]   ->SetLineColor(color);   hChi2PerClusterTPC[i]   ->SetLineWidth(2);
    						     					      
    hC11[i]                 ->SetLineColor(color);   hC11[i]                 ->SetLineWidth(2);
    hC22[i]                 ->SetLineColor(color);   hC22[i]                 ->SetLineWidth(2);
    hC33[i]                 ->SetLineColor(color);   hC33[i]                 ->SetLineWidth(2);
    hC44[i]                 ->SetLineColor(color);   hC44[i]                 ->SetLineWidth(2);
    hC55[i]                 ->SetLineColor(color);   hC55[i]                 ->SetLineWidth(2);
   						     					      
    hDXY[i]                 ->SetLineColor(color);   hDXY[i]                 ->SetLineWidth(2);
    hDZ[i]                  ->SetLineColor(color);   hDZ[i]                  ->SetLineWidth(2);
						     
    hDXYNormalized[i]       ->SetLineColor(color);   hDXYNormalized[i]       ->SetLineWidth(2);
    hDZNormalized[i]        ->SetLineColor(color);   hDZNormalized[i]        ->SetLineWidth(2);

  }
}

//____________________________________________________________________
void 
AliESDtrackCuts::Print(const Option_t*) const {

  AliInfo("AliESDtrackCuts...");
}


//____________________________________________________________________
void 
AliESDtrackCuts::SaveHistograms(Char_t* dir) {
  
  if (!fHistogramsOn) {
    AliDebug(0, "Histograms not on - cannot save histograms!!!");
    return;
  }

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");
 
  hCutStatistics->Write();
  hCutCorrelation->Write();

  for (Int_t i=0; i<2; i++) {
    if (i==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");
    
    hNClustersITS[i]        ->Write();
    hNClustersTPC[i]        ->Write();
    hChi2PerClusterITS[i]   ->Write();
    hChi2PerClusterTPC[i]   ->Write();
    
    hC11[i]                 ->Write();
    hC22[i]                 ->Write();
    hC33[i]                 ->Write();
    hC44[i]                 ->Write();
    hC55[i]                 ->Write();

    hDXY[i]                 ->Write();
    hDZ[i]                  ->Write();
    hDXYvsDZ[i]             ->Write();
    
    hDXYNormalized[i]       ->Write();
    hDZNormalized[i]        ->Write();
    hDXYvsDZNormalized[i]   ->Write();

    gDirectory->cd("../");
  }

  gDirectory->cd("../");
}



