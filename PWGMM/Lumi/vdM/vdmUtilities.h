//-------------------------------------------------------
// this header gathers functions that are used somewhere
// else in the code
//-------------------------------------------------------

//-------------------------------------------------------
// Computation of cross sections
//-------------------------------------------------------

//-------------------------------------------------------
// cross section

Double_t GetXS( Double_t Ax, Double_t Ay, Double_t R00x, Double_t R00y, Double_t N1, Double_t N2)
// compute the visible cross sections given the areas, the rates and the intensities
{
  const Double_t kF = 1e+25;   // conversion factor from mm^2 to mili-barns

  Double_t Qx = R00x/Ax;
  Double_t Qy = R00y/Ay;
  Double_t R00 = 0.5*(R00x+R00y);
  return kF*R00/(g_kLHCFreq*N1*N2*Qx*Qy); 
}

//-------------------------------------------------------
// cross section error

Double_t GetXSerr(Double_t Ax, Double_t Axe, Double_t Ay, Double_t Aye, 
		   Double_t R00x, Double_t R00xe, Double_t R00y, Double_t R00ye, 
		   Double_t N1, Double_t N2)
// compute the error on the cross sections given the areas,
// the rates, and the intensities,
// as well as the corresponding errors
  
{
  // cross section
  Double_t xs = GetXS(Ax,Ay,R00x,R00y,N1,N2);
  // uncertainty from rates
  Double_t F = (1.0/R00x)+(1.0/R00y);
  Double_t DeltaF_2 = (((R00xe*R00xe)/(R00x*R00x*R00x*R00x)) +
		       ((R00ye*R00ye)/(R00y*R00y*R00y*R00y)))/(F*F);
  // uncertainty from areas
  Double_t DeltaAx_2 = (Axe*Axe)/(Ax*Ax);
  Double_t DeltaAy_2 = (Aye*Aye)/(Ay*Ay);
  
  // total uncertainty
  return xs*TMath::Sqrt(DeltaAx_2+DeltaAy_2+DeltaF_2);
}


//-------------------------------------------------------
// Set global pointers to input files and trees
//-------------------------------------------------------

void Set_pointers_to_input_files_and_trees()
{
  // set pointers to files
  g_vdm_File = new TFile(g_Input_vdm_File);
  g_vdm_DDL2_File = new TFile(g_Input_vdm_DDL2_File);
  g_vdm_BPTX_File = new TFile(g_Input_vdm_BPTX_File);

  // set pointers to trees
  if (g_vdm_File) g_vdm_Tree = (TTree *) g_vdm_File->Get("VdM");
  if (g_vdm_File) g_vdm_DDL2_Tree = (TTree *) g_vdm_DDL2_File->Get("DDL2");
  if (g_vdm_File) g_vdm_BPTX_Tree = (TTree *) g_vdm_BPTX_File->Get("VdM-BPTX");

  // add trees as friends of main tree
  g_vdm_Tree->AddFriend(g_vdm_BPTX_Tree);
  g_vdm_Tree->AddFriend(g_vdm_DDL2_Tree);
  
}

//-------------------------------------------------------
// Find start and end of x and y scans
//-------------------------------------------------------

void Find_start_and_end_of_scans()
{
  // cout << " Finding start and end of scans. It will take a moment."<< endl;
  
  // set  branches
  Int_t plane;
  g_vdm_Tree->ResetBranchAddresses();
  g_vdm_Tree->SetBranchAddress("plane",&plane);  

  // get storage space and initialize it
  g_Idx_Start_Scan_x = new Int_t[g_n_Scans_in_Fill];
  g_Idx_Start_Scan_y = new Int_t[g_n_Scans_in_Fill];
  g_Idx_End_Scan_x = new Int_t[g_n_Scans_in_Fill];
  g_Idx_End_Scan_y = new Int_t[g_n_Scans_in_Fill];
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++)
    g_Idx_Start_Scan_x[i]=g_Idx_Start_Scan_y[i]=g_Idx_End_Scan_x[i]=g_Idx_End_Scan_y[i]=-1;

  // loop over tree to find x-scans
  Int_t n_vdm = g_vdm_Tree->GetEntries();
  Int_t found_scans = 0;
  Int_t idx=0;
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    g_vdm_Tree->GetEntry(idx);
    while(plane!=1 && idx<n_vdm) { // find start of x-scan 
      idx++;
      g_vdm_Tree->GetEntry(idx);
    }
    g_Idx_Start_Scan_x[i]=idx;
    while(plane==1 && idx<n_vdm) { // find end of x-scan 
      idx++;
      g_vdm_Tree->GetEntry(idx);
    }
    g_Idx_End_Scan_x[i]=idx;
    idx++;
    found_scans++;
  }

  if (found_scans == g_n_Scans_in_Fill) {
    /*
      cout << found_scans << " x-scans found: " << endl;
    for(Int_t i=0;i<found_scans;i++)
      cout << "   indices are " << g_Idx_Start_Scan_x[i] << " and " << g_Idx_End_Scan_x[i] << endl;
    */
  } else {
    cout << " Number of found x-scans ("<< found_scans <<") different from expectations: " << g_n_Scans_in_Fill << endl;
    exit(-101);
  }
  
  // loop over tree to find y scans
  found_scans = 0;
  idx=0;
  for(Int_t i=0;i<g_n_Scans_in_Fill;i++) {
    g_vdm_Tree->GetEntry(idx);
    while(plane!=2 && idx<n_vdm) { // find start of y-scan 
      idx++;
      g_vdm_Tree->GetEntry(idx);
    }
    g_Idx_Start_Scan_y[i]=idx;
    while(plane==2 && idx<n_vdm) { // find end of y-scan 
      idx++;
      g_vdm_Tree->GetEntry(idx);
    }
    g_Idx_End_Scan_y[i]=idx;
    idx++;
    found_scans++;
  }

  if (found_scans == g_n_Scans_in_Fill) {
    /*
    cout << found_scans << " y-scans found: " << endl;
    for(Int_t i=0;i<found_scans;i++)
      cout << "   indices are " << g_Idx_Start_Scan_y[i] << " and " << g_Idx_End_Scan_y[i] << endl;
    */
  } else {
    cout << " Number of found y-scans ("<< found_scans <<") different from expectations: " << g_n_Scans_in_Fill << endl;
    exit(-102);
  }


}

//-------------------------------------------------------
// Find indices of start and end of scan
//-------------------------------------------------------

void FindIdicesOfScanStartEnd(Int_t scan, Int_t *indices)

{
   // get separation files
  char file_name[kg_string_size];
  sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
  TFile *ScanFileX = new TFile(file_name);
  if (ScanFileX == NULL) {
    cout << " file " << file_name << " not found " << endl;
    exit(-103);
  }
  sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);
  TFile *ScanFileY = new TFile(file_name);
  if (ScanFileY == NULL) {
    cout << " file " << file_name << " not found " << endl;
    exit(-103);
  }
  
 // find first index of x scan
  TTree *scan_tree_x = (TTree *) ScanFileX->Get("SepInfo");
  Int_t idx_separation_start = -1;  
  scan_tree_x->ResetBranchAddresses();
  scan_tree_x->SetBranchAddress("idx_separation_start",&idx_separation_start);
  scan_tree_x->GetEntry(0); // get first entry

  // find last index of y scan
  TTree *scan_tree_y = (TTree *) ScanFileY->Get("SepInfo");
  Int_t n_separations_y = scan_tree_y->GetEntries();
  Int_t idx_separation_end = -1;  
  scan_tree_y->ResetBranchAddresses();
  scan_tree_y->SetBranchAddress("idx_separation_end",&idx_separation_end);
  scan_tree_y->GetEntry(n_separations_y-1); // get last entry

  // return values
  indices[0] = idx_separation_start;
  indices[1] = idx_separation_end;  
}

//-------------------------------------------------------
// Find index between x and y scan
//-------------------------------------------------------

Int_t FindIdxBetweenScans(Int_t scan)

{
  // get separation files
  char file_name[kg_string_size];
  sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
  TFile *ScanFileX = new TFile(file_name);
  if (ScanFileX == NULL) {
    cout << " file " << file_name << " not found " << endl;
    exit(-103);
  }
  sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);
  TFile *ScanFileY = new TFile(file_name);
  if (ScanFileY == NULL) {
    cout << " file " << file_name << " not found " << endl;
    exit(-103);
  }
  

  // find last index of x scan
  TTree *scan_tree_x = (TTree *) ScanFileX->Get("SepInfo");
  Int_t n_separations_x = scan_tree_x->GetEntries();
  Int_t idx_separation_end = -1;  
  scan_tree_x->ResetBranchAddresses();
  scan_tree_x->SetBranchAddress("idx_separation_end",&idx_separation_end);
  scan_tree_x->GetEntry(n_separations_x-1); // get last entry

  // find first index of y scan
  TTree *scan_tree_y = (TTree *) ScanFileY->Get("SepInfo");
  Int_t idx_separation_start = -1;  
  scan_tree_y->ResetBranchAddresses();
  scan_tree_y->SetBranchAddress("idx_separation_start",&idx_separation_start);
  scan_tree_y->GetEntry(0); // get first entry

  Int_t idx = (Int_t) ((((Double_t) idx_separation_start) + ((Double_t) idx_separation_end))*0.5);
  return idx;
}

//-------------------------------------------------------
// Find the the number of steps in a given scan
//-------------------------------------------------------

Int_t FindNumberSeparations(Int_t scan_type, Int_t scan)
// scan_type: 1 => x-scan; 2 => y-scan
{
  // get correct separation file
  char file_name[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);
  TFile *ScanFile = new TFile(file_name);
  if (ScanFile == NULL) {
    cout << " file " << file_name << " not found " << endl;
    exit(-103);
  }

  TTree *scan_tree = (TTree *) ScanFile->Get("SepInfo");
  Int_t n_separations = scan_tree->GetEntries();
  return n_separations;

}

//-------------------------------------------------------
// Find the start and end of each step in the scan
//-------------------------------------------------------

void FindStepStartAndEnd(Int_t scan_type, Int_t scan, Int_t n_separations, Int_t *idx_start,Int_t *idx_end)
// scan_type: 1 => x-scan; 2 => y-scan
{
  // get correct separation file
  char file_name[kg_string_size];
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",g_vdm_Fill,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",g_vdm_Fill,scan);
  TFile *ScanFile = new TFile(file_name);
  if (ScanFile == NULL) {
    cout << " file " << file_name << " not found " << endl;
    exit(-103);
  }
  // set up separation tree
  Int_t idx_separation_start;
  Int_t idx_separation_end;  
  TTree *scan_tree = (TTree *) ScanFile->Get("SepInfo");
  scan_tree->ResetBranchAddresses();
  scan_tree->SetBranchAddress("idx_separation_start",&idx_separation_start);
  scan_tree->SetBranchAddress("idx_separation_end",&idx_separation_end);

  // fill separations
  for(Int_t i=0;i<n_separations;i++) {
    scan_tree->GetEntry(i);
    idx_start[i]=idx_separation_start;
    idx_end[i]=idx_separation_end;    
  }
}


//-------------------------------------------------------
// Find index in histogram corresponding to time
// (used to match DCCT histogram index to relative time in tree
//-------------------------------------------------------

Int_t GetHistogramIndex(TH1 *histo, Double_t time)
{
  Double_t diff = 1000000.0; // large number
  Int_t nbx = histo->GetNbinsX();
  Int_t idx = -1;
  for(Int_t j=0;j<=nbx;j++) {
    Double_t th = histo->GetBinCenter(j);
    if (TMath::Abs(th-time)<diff) {
      idx=j;
      diff = TMath::Abs(th-time);
    }
  }
  return idx;
}


//-------------------------------------------------------
// compute the raw rate and the corresponding error
//-------------------------------------------------------

Double_t RateRaw(Double_t counts, Double_t orbits)
// compute raw rate
{
  return g_kLHCFreq*counts/orbits;
}

Double_t RateRawErr(Double_t counts, Double_t orbits)
// compute error on raw rate
{
  return g_kLHCFreq*TMath::Sqrt(counts)/orbits;
}

//-------------------------------------------------------
// compute the bkdg correction and the corresponding error
//-------------------------------------------------------

Double_t BkgdCorrection(Double_t acc, Double_t tot)
// compute background correction factor
{
  if (tot == 0.0) return 1.0;
  return acc/tot;
}

//-------------------------------------------------------

Double_t BkgdCorrectionError(Double_t acc, Double_t tot)
// compute error on background correction factor
{
  if (tot == 0.0) return 0.0;
  return TMath::Sqrt(acc*(1-(acc/tot)))/tot;
}

//-------------------------------------------------------
// compute the error of the bkdg corrected rate
//-------------------------------------------------------

Double_t BkgdCorrectedRateError(Double_t r, Double_t re,Double_t f, Double_t fe)
// compute error on the background corrected rate
{
  return TMath::Sqrt(fe*fe*r*r+re*re*f*f);
}


//-------------------------------------------------------
// tools to compute the pile up as proposed by Martino
//-------------------------------------------------------

// this function defines the relationship between ratePerBC and mu
// Code by Martino
Double_t GetPileUp(Double_t * x, Double_t * par)
{

  //printf("using %f and %f\n",par[0],par[1]);
  
  Double_t eMinMu = TMath::Exp(-x[0]);
  Double_t eMinMuA = TMath::Exp(-x[0]*par[0]);
  Double_t eMinMuC = TMath::Exp(-x[0]*par[1]);
  
  return (1-eMinMu + eMinMu * (1-eMinMuA) * (1-eMinMuC));
}

// this function is called to invert numerically the GetPileUp function
// Code by Martino

Double_t trova(Double_t y, TF1* fun)
{
  //
  Double_t xmin=fun->GetXmin();
  Double_t xmax=fun->GetXmax();
  Double_t xmed,val;
  for(Int_t i=0; i<100000; i++){
    xmed=(xmax+xmin)/2.;
    val=fun->Eval(xmed);
    if(val>=y)
      xmax=xmed;
    else
      xmin=xmed;
  }
  xmed=(xmax+xmin)/2.;
  return xmed;
}

//-------------------------------------------------------
// this function should be called once per rate measurement, to get the corrected mu value.
// It takes as argument the total  counts divided by the LHC frequency
// Code by Martino
Double_t CorrectRateForPileUp(Double_t ratePerBC, TF1 *fPU)
{  
  Double_t Prob = trova(ratePerBC, fPU);
  return Prob;
}

//-------------------------------------------------------

Double_t RatePileUp(Double_t r, TF1 *fPU)
// compute pile-up rate
// r = rate
{
  Double_t rPU = g_kLHCFreq*CorrectRateForPileUp(r/g_kLHCFreq,fPU);
  return rPU;
}

//-------------------------------------------------------

Double_t RatePileUpErr(Double_t r, Double_t re , TF1 *fPU)
// compute error on pile-up rate
// r = rate, re = rate error
{
  Double_t rPUp = g_kLHCFreq*CorrectRateForPileUp((r+re)/g_kLHCFreq,fPU);
  Double_t rPUm = g_kLHCFreq*CorrectRateForPileUp((r-re)/g_kLHCFreq,fPU);
  Double_t err = 0.5*(rPUp-rPUm);
  return err;
}


//-------------------------------------------------------
// get bunch crossing information
//-------------------------------------------------------

Int_t GetNumberInteractingBunchCrossings()
{
  TTree *Tree = (TTree *) g_vdm_File->Get("BCinteracting");
  return Tree->GetEntries();
}

//-------------------------------------------------------
void GetBucketInfo(Int_t *BucketA, Int_t *BucketC)
{
  // get tree and set branches
  Int_t lhca;
  Int_t lhcc;
  TTree *Tree = (TTree *) g_vdm_File->Get("BCinteracting");
  TBranch *Branch2 = Tree->GetBranch("lhcA");
  Branch2->SetAddress(&lhca);
  TBranch *Branch3 = Tree->GetBranch("lhcC");
  Branch3->SetAddress(&lhcc);

  // loop over tree to fill in info
  Int_t nBC = Tree->GetEntries();
  for (Int_t i=0;i<nBC;i++) {
    Tree->GetEntry(i);
    BucketA[i]=lhca;
    BucketC[i]=lhcc;
  }
}

//-------------------------------------------------------
void GetBunchIndices(Int_t *bunches)
{
  // get tree and set branches
  Int_t bc;
  TTree *Tree = (TTree *) g_vdm_File->Get("BCinteracting");
  TBranch *Branch = Tree->GetBranch("BC");
  Branch->SetAddress(&bc);

  // loop over tree to fill in info
  Int_t nBC = Tree->GetEntries();
  for (Int_t i=0;i<nBC;i++) {
    Tree->GetEntry(i);
    bunches[i]=bc;
  }
}
