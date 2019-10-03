/**
 * @file   PlotSysInfo.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Oct 15 13:21:47 2014
 * 
 * @brief  A script to plot information from system watch files
 */
/** 
 * Make a canvas 
 * 
 * @param title Title on canvas 
 * 
 * @return The canvas 
 */
TVirtualPad* MakeCanvas(const char* title="")
{
  static Int_t cId = 0;
  TCanvas* c = new TCanvas(Form("c%02d", ++cId), title);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->cd();
  return c;
}
/** 
 * Summarize the usage 
 * 
 * @param tree tree to look in 
 * @param exp  expression to draw 
 * @param cut  Cut to use 
 * @param draw Whether to draw or not 
 * 
 * @return The summed usage 
 */
Double_t SumUsage(TTree*      tree, 
		  const char* exp, 
		  const char* cut, 
		  bool        draw=false)
{
  if (draw) MakeCanvas(Form("%s [%s]", exp, cut));
  
  Int_t  entries = tree->Draw(exp, cut, (draw ? "" : "goff"));
  if (entries==0) return 0;

  Double_t mean = TMath::Mean(entries, tree->GetV1());
  return entries * mean;
}
/** 
 * Return the most heavy load 
 * 
 * @param tree  tree to look in 
 * @param exp   expression to draw 
 * @param cut   Cut to use 
 * @param order Off set
 * 
 * @return Most heavy usage 
 */
Double_t TopUsage(TTree*      tree, 
		  const char* exp, 
		  const char* cut, 
		  Int_t       order)
{
  Int_t entries = tree->Draw(exp, cut, "goff");
  if (entries <= 1 || !tree->GetV1()) return -10000;
  
  TArrayI index(entries);
  TMath::Sort(entries, tree->GetV1(), index.fArray);

  Int_t    oindex = TMath::Min(order, entries);
  Double_t value  = tree->GetV1()[index[oindex-1]];
  
  return value;
}
/** 
 * Extract a histogram from a resource spec 
 * 
 * @param tree    tree to look in 
 * @param exp     expression to draw 
 * @param cut     Cut to use 
 * @param name    Name of histogram
 * @param xtitle  X axis title
 * @param ytitle  Y axis title 
 * @param draw    If true, draw 
 * 
 * @return The extract histogram or null
 */
  
TH1* ExtractHist(TTree*      tree,
		 const char* exp, 
		 const char* cut,
		 const char* name,
		 const char* xtitle="", 
		 const char* ytitle="",
		 Bool_t      draw=false)
{
  tree->Draw(Form("%s>>tmpa", exp), cut, "GOFF");
  if (!tree->GetHistogram()) return 0;

  TH1* ret = static_cast<TH1*>(tree->GetHistogram()->Clone(name));
  delete tree->GetHistogram();
  ret->SetXTitle(xtitle);
  ret->SetYTitle(ytitle);
  ret->SetMarkerStyle(22);
  ret->SetMarkerSize(1);
  if (draw) {
    if (draw) MakeCanvas(name);
    ret->Draw();
  }
  return ret;
}
		 
/** 
 * Print some information to output stream
 *  
 * @param o             Stream
 * @param name          Step
 * @param dT            Change in time
 * @param dVM           Change in VM usage
 * @param alldT         Sum of time
 * @param alldVM        Sum of VM usage 
 */		 
void Print(std::ostream& o,
	   const char*   name, 
	   Double_t      dT, 
	   Double_t      dVM, 
	   Double_t      alldT, 
	   Double_t      alldVM)
{
  o << name << "\t" 
    << dT   << "\t" 
    << dVM  << "\t" 
    << 100*(alldT  > 0 ? dT  / alldT  : 0) << "\t" 
    << 100*(alldVM > 0 ? dVM / alldVM : 0) << std::endl;
}
/**
 * Detectors
 * 
 */
const char* dets[] = {"ITS", 
		      "TPC", 
		      "TRD", 
		      "TOF", 
		      "PHOS", 
		      "HMPID", 
		      "EMCAL", 
		      "MUON", 
		      "FMD", 
		      "ZDC", 
		      "PMD", 
		      "T0", 
		      "VZERO", 
		      "ACORDE", 
		      "HLT",
		      0 };

/** 
 * Plot information from one file 
 * 
 * @param file  File to plot from 
 * @param draw  Drawing flags
 */
void
Plot1SysInfo(const char* file, UShort_t draw=0x1)
{
  // --- Create output file and tree ---------------------------------
  TString rootOut(file); 
  rootOut.ReplaceAll(".log", ".root");
  rootOut.ReplaceAll(".foo", ".root");
  Info("", "Writing to ROOT file %s", rootOut.Data());
  TFile* out  = TFile::Open(rootOut, "RECREATE");
  TTree* tree = AliSysInfo::MakeTree(file);
  tree->SetName("T");
  
  
  // --- Create ASCII output ------------------------------------------
  TString sumOut(rootOut); 
  sumOut.ReplaceAll(".root", ".sum");
  Info("", "Writing to ASCII file %s", sumOut.Data());
  std::ofstream ascii(sumOut.Data());
  ascii << "Det/C:sumDt/F:sumDvm/F:fracDt/F:fracDvm/F" << std::endl;

  // --- Get global stuff ---------------------------------------------
  const char* all     = "id0>=0&&id2>=0";
  Double_t sumdTAll   = SumUsage(tree, "deltaT",  all, draw & 0x1);
  Double_t sumdVMAll  = SumUsage(tree, "deltaVM", all, draw & 0x1);
  Double_t topdT      = TopUsage(tree, "deltaT",  "id2<3", 20);
  Double_t topdVM     = TopUsage(tree, "deltaVM", "", 20);
  TCut     cutT("cutDT", Form("deltaT > %f", topdT));
  TCut     cutVM("cutVM", Form("deltaVM > %f", topdVM));
  
  ExtractHist(tree, "deltaVM:sname", "1"+cutVM,
	      "DVMvsName","","#DeltaVM [MB]", draw&0x4);
  ExtractHist(tree, "VM:sname", "id2<3"+cutVM, 
	      "VMvsName", "", "VM [MB]", draw&0x4);
  ExtractHist(tree, "VM:T", "deltaVM>1", 
	      "VMvsTime", "Time [sec]", "VM [MB]", draw&0x4);
  ExtractHist(tree, "deltaT:sname","id2<3"+cutT,
	      "CPUvsName","","#DeltaT [sec]", draw&0x4);
  


  Print(ascii, "all", sumdTAll, sumdVMAll, sumdTAll, sumdVMAll);

  

  // --- Loop over detetors ------------------------------------------
  const char** pdet = dets;
  Int_t        idet = 0;
  while (*pdet) { 
    TString  cut    = Form("id0==%d && id2 >= 0", idet);
    Double_t sumdT  = SumUsage(tree, "deltaT",  cut, draw & 0x2);
    Double_t sumdVM = SumUsage(tree, "deltaVM", cut, draw & 0x2);
    Print(ascii, *pdet, sumdT, sumdVM, sumdTAll, sumdVMAll);
    
#if 0
    TString cut2    = Form("id0==%d",idet);
    ExtractHist(tree, "deltaVM:sname", cut2.Data()+cutVM,
		Form("DVMvsName_%02d", idet), "", "#DeltaVM [MB]", draw&0x8);
    ExtractHist(tree, "VM:sname",      cut2.Data()+cutVM,
		Form("VMvsName_%02d", idet), "", "VM [MB]", draw&0x8);
    ExtractHist(tree, "deltaT:sname",  cut2.Data()+cutT, 
		Form("CPUvsName_%02d", idet),"", "#DeltaT [sec]", draw&0x8);
#endif

    pdet++;
    idet++;
  }
  ascii.close();
  tree->Write();
  out->Write();
  out->ls();
  new TBrowser;
}

/** 
 * Plot for one job both simulation and reconstruction usage 
 * 
 * @param pid Job identifier 
 */
void
PlotSysInfo(ULong_t pid=431808952)
{
  Plot1SysInfo(Form("%d_simwatch.log", pid));
  Plot1SysInfo(Form("%d_recowatch.log", pid));
}
// 
// EOF
// 

  
