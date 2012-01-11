// $Id$
/*
 * Plotting macro for comparing offline- and HLT- ESD trees from  
 * HLT-OFFLINE-GLOBAL-comparison.root produced using pwgpp-task:
 * compare-HLT-offline-local.C'("./AliESDs.root","pwgpp")' 
 * 
 * It allows you to choose from a detailed list of cuts or a combination of cuts.
 * 
 * Apply requested cuts in L300
 * Int_t  nCutsMin = 2;
 * Int_t  nCutsMax = 5;
 * Int_t  nCuts=nCutsMax+1;
 *
 * The example set above will set the cuts from 2 to 5, defined in function SetCuts().
 *
 * The current version of this macro needs modifications in order 
 * to set options here: (L312) 
 * Int_t  nCanMin = 0;
 * Int_t  nCanMax = 0 ;
 * Int_t  nCans=nCanMax+1;
 *
 *
 * Usage:
 * Running requires that you have the .root-files produced from the pwgpp-task 
 * in your local folder. 
 *
 * Run options:
 * 1) Run with script ./draw.sh This will create folders and organize the output files for you.  
 * 2) Run as aliroot -q drawPerformanceTPCQAofflineHLT.C'("./")'
 * When using the script you will also have to apply changes in the drax-script when applying cuts: 
 * Here example 2 to 5:
 * "for ii in {2..5} ; do"
 * NB! This line occurs 2 places
 *
 *
 * @ingroup alihlt_qa
 * @author jochen@thaeder.de, Camilla.Stokkevag@cern.ch, Kalliopi.Kanaki@ift.uib.no 
 */


// ----------------------------------------------------------
// ----------------------------------------------------------

void setCuts(THnSparse *htrack, Int_t cuts ) {
  
  if (cuts == 0) {       // no_cuts --------------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-150,149.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-150,149.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }

  else if (cuts == 1) { // nClust_[50,160]_+_p_t_[0.3,9999]_+_DCA_[-7,6.99]_+_Eta_[-1.0,1.0] -----
    htrack->GetAxis(0)->SetRangeUser(50,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-7,6.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-7,6.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.3,9999);     // Pt
    htrack->GetAxis(5)->SetRangeUser(-1, 1);    // Eta 

  }

  // ===========================================================================
  // == Single Cuts
  // ===========================================================================

  else if (cuts == 2) { // nClust_[70,160] -------------------------------------
    htrack->GetAxis(0)->SetRangeUser(70,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 3) { // nClust_[60,160] -------------------------------------
    htrack->GetAxis(0)->SetRangeUser(60,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 4) { // nClust_[50,160] -------------------------------------
    htrack->GetAxis(0)->SetRangeUser(50,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 5) { // p_t_[0.4,10] ----------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.4,10);       // Pt
  }
  else if (cuts == 6) { // p_t_[0.3,10] ----------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.3,10);       // Pt
  }
  else if (cuts == 7) { // p_t_[0.2,10] ----------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.2,10);       // Pt
  }
  else if (cuts == 8) { // DCA_[-3,2.99] ---------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-3,2.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-3,2.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 9) { // DCA_[-5,4.99] ---------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-5,4.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-5,4.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 10) { // DCA_[-7,6.99] ---------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-7,6.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-7,6.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }



  // ===========================================================================
  // == Combined Cuts
  // ===========================================================================

  else if (cuts == 11) { // nClust_[70,160]_+_p_t_[0.4,10] ---------------------
    htrack->GetAxis(0)->SetRangeUser(70,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.4,10);       // Pt
  }
  else if (cuts == 12) { // nClust_[70,160]_+_DCA_[-3,2.99] --------------------
    htrack->GetAxis(0)->SetRangeUser(70,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-3,2.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-3,2.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 13) { // DCA_[-3,2.99]_+_p_t_[0.4,10] -----------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-3,2.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-3,2.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.4,10);       // Pt
  }
  else if (cuts == 14) { // nClust_[60,160]_+_p_t_[0.3,10] ---------------------
    htrack->GetAxis(0)->SetRangeUser(60,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.3,10);       // Pt
  }
  else if (cuts == 15) { // nClust_[60,160]_+_DCA_[-5,4.99] --------------------
    htrack->GetAxis(0)->SetRangeUser(60,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-5,4.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-5,4.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 16) { // DCA_[-5,4.99]_+_p_t_[0.3,10] -----------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-5,4.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-5,4.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.3,10);       // Pt
  }

  else if (cuts == 17) { // nClust_[50,160]_+_p_t_[0.2,10] ---------------------
    htrack->GetAxis(0)->SetRangeUser(50,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.2,10);       // Pt
  }
  else if (cuts == 18) { // nClust_[50,160]_+_DCA_[-7,6.99] --------------------
    htrack->GetAxis(0)->SetRangeUser(50,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-7,6.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-7,6.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 19) { // DCA_[-7,6.99]_+_p_t_[0.2,10] -----------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-7,6.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-7,6.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.2,10);       // Pt
  }

  // ===========================================================================

  else if (cuts == 20) { // nClust_[70,160]_+_p_t_[0.4,10]_+_DCA_[-3,2.99] -----
    htrack->GetAxis(0)->SetRangeUser(70,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-3,2.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-3,2.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.4,10);       // Pt
  }
  else if (cuts == 21) { // nClust_[60,160]_+_p_t_[0.3,10]_+_DCA_[-6,4.99] -----
    htrack->GetAxis(0)->SetRangeUser(60,160);       // nClust
    htrack->GetAxis(3)->SetRangeUser(-5,4.99);      // DCAr
    htrack->GetAxis(4)->SetRangeUser(-5,4.99);      // DCAz
    htrack->GetAxis(7)->SetRangeUser(0.3,10);       // Pt
  }


  // ===========================================================================
  else if (cuts == 22) { // p_t_[1.0,10] ---------------------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(1.0,10);       // Pt
  }
  else if (cuts == 23) { // ControlCut_p_t_[0,1.0] -----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,1.0);        // Pt
  }
  // ===========================================================================
  // == Control Cuts
  // ===========================================================================

  else if (cuts == 24) { // ControlCut_nClust_[0,70] ---------------------------
    htrack->GetAxis(0)->SetRangeUser(0,70);         // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 25) { // ControlCut_nClust_[0,60] ---------------------------
    htrack->GetAxis(0)->SetRangeUser(0,60);         // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 26) { // ControlCut_nClust_[0,50] ---------------------------
    htrack->GetAxis(0)->SetRangeUser(0,50);         // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 27) { // ControlCut_p_t_[0,0.4] -----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,0.4);        // Pt
  }
  else if (cuts == 28) { // ControlCut_p_t_[0,0.3] -----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,0.3);        // Pt
  }
  else if (cuts == 29) { // ControlCut_p_t_[0,0.2] -----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,49.99);    // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,49.99);    // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,0.2);        // Pt
  }
  else if (cuts == 30) { // ControlCut_DCA_[2.99,49.99] ------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(2.99,49.99);   // DCAr
    htrack->GetAxis(4)->SetRangeUser(2.99,49.99);   // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 31) { // ControlCut_DCA_[-50,-3] ----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,-3);       // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,-3);       // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 32) { // ControlCut_DCA_[4.99,49.99] ------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(4.99,49.99);   // DCAr
    htrack->GetAxis(4)->SetRangeUser(4.99,49.99);   // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 33) { // ControlCut_DCA_[-50,-5] ----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,-5);       // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,-5);       // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 34) { // ControlCut_DCA_[6.99,49.99] ------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(6.99,49.99);   // DCAr
    htrack->GetAxis(4)->SetRangeUser(6.99,49.99);   // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }
  else if (cuts == 35) { // ControlCut_DCA_[-50,-7] ----------------------------
    htrack->GetAxis(0)->SetRangeUser(0,160);        // nClust
    htrack->GetAxis(3)->SetRangeUser(-50,-7);       // DCAr
    htrack->GetAxis(4)->SetRangeUser(-50,-7);       // DCAz
    htrack->GetAxis(7)->SetRangeUser(0,9999);       // Pt
  }

  // ===========================================================================
  if(cuts==0)
    htrack->GetAxis(5)->SetRangeUser(-1.5, 1.5);    // eta -1.0 - 1.0
  else
    htrack->GetAxis(5)->SetRangeUser(-1, 1);        // eta -1.0 - 1.0
  htrack->GetAxis(8)->SetRangeUser(-1, 1);          // Charge

  // ===========================================================================
}

// ----------------------------------------------------------
// ----------------------------------------------------------


drawPerformanceTPCQAofflineHLT(const Char_t* folder = "../..") {



  // ----------------------------------------------------
  // ------ APPLY YOUR CUTS -----------------------------
  // ----------------------------------------------------

  Int_t  nCutsMin = 0;
  Int_t  nCutsMax = 5;
  Int_t  nCuts=nCutsMax+1;

  //----------------------------------------------------


  //----------------------------------------------------
  //----CURRENT VERSION OF MACRO DOES NOT WORK ---------
  //-------IF CUTS ARE APPLIED IN LINES UNDER----------
  //----------------------------------------------------

  Int_t  nCanMin = 0;
  Int_t  nCanMax = 2;
  Int_t  nCans=nCanMax+1;

  // ----------------------------------------------------
  // ----------------------------------------------------

  gSystem->Load("libANALYSIS");                        
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGPP");

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(10);
  TH1::AddDirectory(kFALSE);

  // ----------------------------------------------------

  TString inFile(folder);
  inFile += "/TPC.Performance.root";
  TFile *file = TFile::Open(inFile.Data());
  if (!file) {
    printf("Error : No File\n");
    return;
  }

  TList* list = static_cast<TList*>(file->Get("TPC"));
  if (!list) {
    printf("Error : No List\n");
    return;
  }

  TString inFileH(folder);
  inFileH += "/HLTTPC.Performance.root";
  TFile *fileH = TFile::Open(inFileH.Data());
  if (!fileH) {
    printf("Error : No HLT File\n");
    return;
  }

  TList* listH = static_cast<TList*>(fileH->Get("HLTQA"));  
  if (!listH) {
    printf("Error : No HLT List\n");
    return;
  }

  // ----------------------------------------------------

  AliPerformanceTPC* obj = static_cast<AliPerformanceTPC*>(list->FindObject("AliPerformanceTPC"));
  if(obj==NULL) return;

  AliPerformanceTPC* objH = static_cast<AliPerformanceTPC*>(listH->FindObject("AliPerformanceTPC"));
  if(objH==NULL) return;

  TObjArray *aCan = new TObjArray();
  
  // ----------------------------------------------------

  /*  Xv:Yv:Zv:mult:multP:multN:vertStatus
      ------------------------------------
      0 Xv
      1 Yv
      2 Zv
      3 mult
      4 multP
      5 multN
      6 VertStatus
  */

  THnSparse *hevent = obj->GetTPCEventHisto(); 
  THnSparse *heventH = objH->GetTPCEventHisto(); 


  //  return;

  /*  nClust:chi2PerClust:nClust/nFindableClust:DCAr:DCAz:eta:phi:pt:charge 
      ---------------------------------------------------------------------
      0 nClust: 
      1 chi2PerClust:
      2 nClust/nFindableClust:
      3 DCAr:
      4 DCAz:
      5 eta:
      6 phi:
      7 pt:
      8 charge :
  */

  THnSparse *htrack  = obj->GetTPCTrackHisto(); 
  THnSparse *htrackH = objH->GetTPCTrackHisto(); 

  // ----------------------------------------------------

  TCanvas *can0 = new TCanvas("can0","TPC offline-HLT event information",550,750); 
  can0->Divide(2,3);
  trackInfo0(can0,hevent,heventH);
  can0->SaveAs("perfImg/qa/event/TPC_offline-HLT_event_info.png");
  can0->SaveAs("perfRoot/TPC_offline_HLT_event_info.root");


  // ----------------------------------------------------


  for (Int_t idxCut = nCutsMin; idxCut <= nCutsMax; ++idxCut){
    for ( Int_t idxCan = 2*nCanMin; idxCan < 2*(nCans); ++idxCan) {

      Int_t idx = idxCan/2;
      
      cout << "idx: " << idx << "    idxCut: " << idxCut <<endl;
  
      if(idx==0){    
	aCan->Add(new TCanvas(Form("can%d_%d", idxCan, idxCut), 
			      Form("TPC HLT tracks1 %d - cuts %d", idx, idxCut), 1600, 1200) );  //1200, 750) );
	(dynamic_cast<TCanvas*>(aCan->Last()))->Divide(3,3);  

	fillCanvas(dynamic_cast<TCanvas*>(aCan->Last()), htrack, htrackH, idxCut, idx, 1.0);

	(dynamic_cast<TCanvas*>(aCan->Last()))->SaveAs(Form("perfImg/qa/%d/TPC_offline_HLT_tracks_%d-%d.png", idxCut, idxCut, idx));
	(dynamic_cast<TCanvas*>(aCan->Last()))->SaveAs(Form("perfRoot/TPC_HLT_tracks_%d-%d.root", idxCut, idx));
      }
      else{
	// -- -- -- -- -- -- -- --
      
	//if(idx==0) continue;

	aCan->Add(new TCanvas(Form("canS%d_%d", idxCan, idxCut), 
			      Form("TPC tracks %d - cuts %d", idx, idxCut), 1200, 750) );
	(dynamic_cast<TCanvas*>(aCan->Last()))->Divide(4,3);  

	fillCanvas(dynamic_cast<TCanvas*>(aCan->Last()), htrack, idxCut, idx, 1.0);

	(dynamic_cast<TCanvas*>(aCan->Last()))->SaveAs(Form("perfImg/qa/%d/TPC_tracks_%d-%d.png", idxCut, idxCut, idx));
	(dynamic_cast<TCanvas*>(aCan->Last()))->SaveAs(Form("perfRoot/TPC_tracks_%d-%d.root", idxCut, idx));
      
	// -- -- -- -- -- -- -- -- 
	++idxCan;
	// -- -- -- -- -- -- -- -- 
      
	aCan->Add(new TCanvas(Form("can%d_%d", idxCan, idxCut), 
			      Form("HLTTPC tracks %d - cuts %d", idx, idxCut), 1200, 750) );
	(dynamic_cast<TCanvas*>(aCan->Last()))->Divide(4,3);
     

	fillCanvas(dynamic_cast<TCanvas*>(aCan->Last()), htrackH, idxCut, idx, 1.0);

	     
	(dynamic_cast<TCanvas*>(aCan->Last()))->SaveAs(Form("perfImg/qa/%d/HLTTPC_tracks_%d-%d.png", idxCut, idxCut, idx));
	(dynamic_cast<TCanvas*>(aCan->Last()))->SaveAs(Form("perfRoot/HLTTPC_tracks_%d-%d.root", idxCut, idx));
      }
    }
  }

  return;
}

//void fillCanvas(TCanvas* can, THnSparse *htrack, Int_t cuts, Int_t idx, Double_t scale ) {
void fillCanvas(TCanvas* can, THnSparse *htrack, THnSparse *htrackH, Int_t cuts, Int_t idx, Double_t scale ) { 
 
  if ( idx == 0 )
    trackInfo1(can, htrack, htrackH, cuts ,scale);

  return;
}

void fillCanvas(TCanvas* can, THnSparse *htrack, Int_t cuts, Int_t idx, Double_t scale ) {
 

  if ( idx == 1 )
    trackInfo2(can, htrack, cuts, scale);
  else if ( idx == 2 )
    trackInfo3(can, htrack, cuts, scale);

  return;
}


//--------------------------------------------------------------------------------------------
//----------- OFFLINE HLT EVENT INFO ---------------------------------------------------------
//--------------------------------------------------------------------------------------------


void trackInfo0(TCanvas* can, THnSparse *hevent,  THnSparse *heventH) {

  TH1D* histe = NULL;
  TH1D* histeH = NULL;
 
  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
  leg1->SetFillColor(10);
  leg1->SetLineColor(10);


 //======================

  histeH = heventH->Projection(6);                 
  histe = hevent->Projection(6);                   // vertexStatus
  histeH->SetTitle("VertexStatus");
  heventH->GetAxis(6)->SetRange(2,2);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  leg1->AddEntry(histeH,"HLT", "l");
  leg1->AddEntry(histe,"OFF", "l");

  can->cd(1);
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("sames");

  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)histeH->FindObject("stats");
  st1->SetLineColor(0);

  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)histe->FindObject("stats");
  st2->SetY1NDC(st1->GetY1NDC()-0.05);
  st2->SetY2NDC(st1->GetY2NDC()-0.05);
  st2->SetLineColor(0);
  st2->SetTextColor(histe->GetLineColor());
  st2->SetFillStyle(0);
  st2->Draw();
  
  //======================

  histeH  = heventH->Projection(0);                
  histe = hevent->Projection(0);             // Xv
  histeH->SetTitle("Vertex X");
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(2);
  gPad->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");  
  leg1->Draw("sames");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();

 //======================

  histeH = heventH->Projection(1);
  histe  = hevent->Projection(1);        // Yv
  histeH->SetTitle("Vertex Y");
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(3);
  gPad->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("sames");

  gPad->Update();
  TPaveStats *st5 = (TPaveStats*)histeH->FindObject("stats");
  st5->SetLineColor(0);

  gPad->Update();
  TPaveStats *st6 = (TPaveStats*)histe->FindObject("stats");
  st6->SetY1NDC(st5->GetY1NDC()-0.05);
  st6->SetY2NDC(st5->GetY2NDC()-0.05);
  st6->SetLineColor(0);
  st6->SetTextColor(histe->GetLineColor());
  st6->SetFillStyle(0);
  st6->Draw();

 //======================


  histeH = heventH->Projection(2); 
  histe  = hevent->Projection(2);                 // Zv
  histeH->SetTitle("Vertex Z");
  histe->SetLineColor(2); 
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(4);
  gPad->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("sames");

  gPad->Update();
  TPaveStats *st7 = (TPaveStats*)histeH->FindObject("stats");
  st7->SetLineColor(0);

  gPad->Update();
  TPaveStats *st8 = (TPaveStats*)histe->FindObject("stats");
  st8->SetY1NDC(st7->GetY1NDC()-0.05);
  st8->SetY2NDC(st7->GetY2NDC()-0.05);
  st8->SetLineColor(0);
  st8->SetTextColor(histe->GetLineColor());
  st8->SetFillStyle(0);
  st8->Draw();

 //======================

  histeH = heventH->Projection(3);
  histe  = hevent->Projection(3);                  // mult
  histeH->SetTitle("mult");
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(5);
  gPad->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");  
  leg1->Draw("sames");

  gPad->Update();
  TPaveStats *st9 = (TPaveStats*)histeH->FindObject("stats");
  st9->SetLineColor(0);

  gPad->Update();
  TPaveStats *st10 = (TPaveStats*)histe->FindObject("stats");
  st10->SetY1NDC(st9->GetY1NDC()-0.05);
  st10->SetY2NDC(st9->GetY2NDC()-0.05);
  st10->SetLineColor(0);
  st10->SetTextColor(histe->GetLineColor());
  st10->SetFillStyle(0);
  st10->Draw();

  return;

}

//--------------------------------------------------------------------------------------------
//----------- OFFLINE HLT TRACK INFO ---------------------------------------------------------
//--------------------------------------------------------------------------------------------

void trackInfo1(TCanvas* can, THnSparse *htrack, THnSparse *htrackH, Int_t cuts, Double_t scale) { 

  setCuts(htrack, cuts);
  setCuts(htrackH, cuts);

  TH1D* histe = NULL;
  TH1D* histeH = NULL; //HLT
 
  TH2D* hist = NULL;
  TPad* pad = NULL;

  TLegend *leg1 = new TLegend(0.6,0.6,0.8,0.8);
  leg1->SetFillColor(10);
  leg1->SetLineColor(10);

  TLegend *leg2 = new TLegend(0.6,0.2,0.8,0.4);
  leg2->SetFillColor(10);
  leg2->SetLineColor(10);

 //======================

  can->cd(1);
  histeH = htrackH->Projection(5);
  histe  = htrack->Projection(5);                // Eta
  histeH->SetTitle("Eta");
  histeH->Scale(1./scale); 
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  leg1->AddEntry(histeH,"HLT", "l");
  leg1->AddEntry(histe,"Offline", "l");
  leg2->AddEntry(histeH,"HLT", "l");
  leg2->AddEntry(histe,"Offline", "l");

  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("same");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();

 //======================
 
  histeH = htrackH->Projection(6);             
  histe  = htrack->Projection(6);                 // Phi
  histeH->SetTitle("Phi");
  histeH->Scale(1./scale);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(2);
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("same");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();

  
 //======================

  histeH  = htrackH->Projection(0); 
  histe  = htrack->Projection(0);                  // nClust
  histeH->SetTitle("nCluster");
  histeH->Scale(1./scale);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(3);
  //  gPad->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");  
  leg2->Draw("sames");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();


 //======================

  histeH  = htrackH->Projection(8);
  histe  = htrack->Projection(8);                  // Charge
  histeH->SetTitle("Charge");
  histeH->Scale(1./scale);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(4);
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("sames");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();


 //======================

  //pad = can->cd(5);
  histeH  = htrackH->Projection(7); 
  histe  = htrack->Projection(7);                  // Pt
  histeH->SetTitle("Pt");
  histeH->Scale(1./scale);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());  

  can->cd(5)->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("same");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();

 //======================

  histeH  = htrackH->Projection(3);
  histe  = htrack->Projection(3);                  // DCAr
  histeH->SetTitle("DCAr");
  histeH->Scale(1./scale);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(6)->SetLogy(); 
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("same");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();

 //======================

  histeH = htrackH->Projection(4); 
  histe  = htrack->Projection(4);                  // DCAz
  histeH->SetTitle("DCAz");
  histeH->Scale(1./scale);
  histe->SetLineColor(2);
  if(histeH->GetMaximum() > histe->GetMaximum()) histe->SetMaximum(1.1*histeH->GetMaximum());
  else histeH->SetMaximum(1.1*histe->GetMaximum());

  can->cd(7)->SetLogy();
  histeH->Draw("histeH");
  histe->Draw("sames");
  leg1->Draw("same");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)histeH->FindObject("stats");
  st3->SetLineColor(0);

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)histe->FindObject("stats");
  st4->SetY1NDC(st3->GetY1NDC()-0.05);
  st4->SetY2NDC(st3->GetY2NDC()-0.05);
  st4->SetLineColor(0);
  st4->SetTextColor(histe->GetLineColor());
  st4->SetFillStyle(0);
  st4->Draw();
}

 //===================================================================

void trackInfo2(TCanvas* can, THnSparse *htrack, Int_t cuts, Double_t scale) {

  TH1D* histe = NULL;
  TH2D* hist = NULL;
  TPad* pad = NULL;

  setCuts(htrack, cuts);                           // reset Cuts

  pad = can->cd(5);
  pad->SetLogz();
  hist  = htrack->Projection(5,6);                 // eta/phi
  hist->SetTitle("Eta-Phi");
  hist->Scale(1./scale);
  hist->Draw("colz");
  
  // -- -- -- -- -- -- -- -- -- -- -- -- --

  can->cd(6);
  histe  = htrack->Projection(0);                  // nClust
  histe->SetTitle("nCluster");
  histe->Scale(1./scale);
  histe->Draw("histe");

  can->cd(7);
  histe  = htrack->Projection(8);                  // Charge
  histe->SetTitle("Charge");
  histe->Scale(1./scale);
  histe->Draw("histe");

  pad = can->cd(8);
  pad->SetLogy();
  histe  = htrack->Projection(7);                  // Pt
  histe->SetTitle("Pt");
  histe->Scale(1./scale);
  histe->Draw("histe");

  can->cd(9);
  histe  = htrack->Projection(3);                  // DCAr
  histe->SetTitle("DCAr");
  histe->Scale(1./scale);
  histe->Draw("histe");

  can->cd(10);
  histe  = htrack->Projection(4);                  // DCAz
  histe->SetTitle("DCAz");
  histe->Scale(1./scale);
  histe->Draw("histe");

  pad = can->cd(11);
  pad->SetLogz();
  hist  = htrack->Projection(3,4);                 // DCAz/DCAr
  hist->SetTitle("DCAr-DCAz");
  hist->Scale(1./scale);
  hist->Draw("colz");

}


void trackInfo2(TCanvas* can, THnSparse *htrack, Int_t cuts, Double_t scale) {

  TH1D* histe = NULL;
  TH2D* hist = NULL;
  TPad* pad = NULL;

  setCuts(htrack, cuts);


  pad = can->cd(1);
  pad->SetLogz();
  hist  = htrack->Projection(3,5);                 // DCAr/Eta
  hist->SetTitle("DCAr-Eta");
  hist->Scale(1./scale);
  hist->Draw("colz");

  pad = can->cd(1);
  pad->SetLogz();
  hist  = htrack->Projection(3,5);                 // DCAr/Eta
  hist->SetTitle("DCAr-Eta");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(2);
  hist  = htrack->Projection(4,5);                 // DCAz/Eta
  hist->SetTitle("DCAz-Eta");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(3);
  hist  = htrack->Projection(3,6);                 // DCAr/Phi
  hist->SetTitle("DCAr-Phi");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(4);
  hist  = htrack->Projection(4,6);                 // DCAz/Phi
  hist->SetTitle("DCAz-Phi");
  hist->Scale(1./scale);
  hist->Draw("colz");

  // -- -- -- -- -- -- -- -- -- -- -- -- --

  pad = can->cd(5);
  pad->SetLogx();
  hist  = htrack->Projection(5,7);                 // Eta/pt
  hist->SetTitle("Eta-Pt");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  pad = can->cd(6);
  pad->SetLogx();
  hist  = htrack->Projection(6,7);                 // Phi/pt
  hist->SetTitle("Phi-Pt");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  pad = can->cd(7);
  pad->SetLogx();
  hist  = htrack->Projection(3,7);                 // DCAr/pt
  hist->SetTitle("DCAr-Pt");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  pad = can->cd(8);
  pad->SetLogx();
  hist  = htrack->Projection(4,7);                 // DCAz/pt
  hist->SetTitle("DCAz-Pt");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  pad = can->cd(9);
  pad->SetLogx();
  hist  = htrack->Projection(0,7);                 // nClust/pt
  hist->SetTitle("NClust-Pt");
  hist->Scale(1./scale);
  hist->Draw("colz"); 
  
  // -- -- -- -- -- -- -- -- -- -- -- -- --

  pad = can->cd(10);
  pad->SetLogy();
  htrack->GetAxis(5)->SetRangeUser(0.,0.799);      // Eta > 0 | A side
  histe  = htrack->Projection(7);                  // Pt
  histe->SetTitle("Pt A vs C side");
  histe->Draw("histe");
  htrack->GetAxis(5)->SetRangeUser(-0.8,-0.0001);  // Eta < 0 | C side
  TH1* h1 = (TH1*)htrack->Projection(7);
  h1->SetLineColor(kRed);
  histe->Scale(1./scale);
  h1->Scale(1./scale);
  h1->Draw("histesame");

  setCuts(htrack, cuts);

  pad = can->cd(11);
  pad->SetLogy();
  htrack->GetAxis(8)->SetRangeUser(0.,1);          // charge > 0
  histe  = htrack->Projection(7);                  // Pt
  histe->SetTitle("Pt Charge > 0 vs Charge < 0");
  histe->Draw("histe");
  htrack->GetAxis(8)->SetRangeUser(-1.,0.);        // charge < 0
  TH1* h2 = (TH1*)htrack->Projection(7);
  h2->SetLineColor(kGreen+2);
  histe->Scale(1./scale);
  h2->Scale(1./scale);
  h2->Draw("histesame");
}


void trackInfo3(TCanvas* can, THnSparse *htrack, Int_t cuts, Double_t scale) {

  setCuts(htrack, cuts);

  TH1D* histe = NULL;
  TH2D* hist = NULL;

  can->cd(1);
  hist  = htrack->Projection(0,5);                 // nClust/eta
  hist->SetTitle("NClust-Eta");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(2);
  hist  = htrack->Projection(0,6);                 // nClust/phi
  hist->SetTitle("NClust-Phi");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(3);
  htrack->GetAxis(5)->SetRangeUser(0.,0.799);      // Eta > 0 | A side
  hist  = htrack->Projection(0,6);                 // nClust/phi
  hist->SetTitle("NClust-Phi A Side");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(4);
  htrack->GetAxis(5)->SetRangeUser(-0.8,-0.0001);  // Eta < 0 | C side
  hist  = htrack->Projection(0,6);                 // nClust/phi
  hist->SetTitle("NClust-Phi - C side");
  hist->Scale(1./scale);
  hist->Draw("colz");

  setCuts(htrack, cuts);

  can->cd(5);
  hist  = htrack->Projection(0,3);                 // nClust/DCAr
  hist->SetTitle("NClust-DCAr");
  hist->Scale(1./scale);
  hist->Draw("colz");

  can->cd(6);
  hist  = htrack->Projection(0,4);                 // nClust/DCAz
  hist->SetTitle("NClust-DCAz");
  hist->Scale(1./scale);
  hist->Draw("colz");

  // -- -- -- -- -- -- -- -- -- -- -- -- -- 

  can->cd(7);
  hist  = htrack->Projection(0,8);                 // NClust/Charge
  hist->SetTitle("NClust-Charge");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  can->cd(8);
  hist  = htrack->Projection(7,8);                 // Pt/Charge
  hist->SetTitle("Pt-Charge");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  can->cd(9);
  hist  = htrack->Projection(5,8);                 // Eta/Charge
  hist->SetTitle("Eta-Charge");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  can->cd(10);
  hist  = htrack->Projection(6,8);                 // Phi/Charge
  hist->SetTitle("Phi-Charge");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  can->cd(11);
  hist  = htrack->Projection(3,8);                 // DCAr/Charge
  hist->SetTitle("DCAr-Charge");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  can->cd(12);
  hist  = htrack->Projection(4,8);                 // DCAz/Charge
  hist->SetTitle("DCAz-Charge");
  hist->Scale(1./scale);
  hist->Draw("colz"); 

  return;


}
