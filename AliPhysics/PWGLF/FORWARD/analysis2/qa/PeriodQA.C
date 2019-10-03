/**
 * @file   PeriodQA.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:35:08 2011
 * 
 * @brief  Script to run the QAPlotter over a tree with runs 
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
class TTree;
class TCanvas;

/** 
 * Run the QAPlotter
 * 
 * The QAPlotter is then run over the merged <tt>trending.root</tt>
 * file and produces two files
 * 
 * - <tt>index.root</tt> which contains TCanvas objects of the
 *   finished plots.  It also contains the TMultiGraph objects painted
 *   in the canvases.
 *
 * - <tt>index.pdf</tt> which is a PDF of the TCanvases mentioned
 *   above.
 * 
 * The QAPlotter will also produce PNGs of each canvas. 
 *
 * @param input  Input file
 * @param type   Data type (data or sim)
 * @param year   Year 
 * @param period Period (e.g., LHC10h)
 * @param pass   Pass (e.g., pass1, cpass1, passMC)
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void
PeriodQA(const char* input, 
	 const char* type, 
	 Int_t       year, 
	 const char* period,
	 const char* pass) 
{
  Bool_t useVar = true;
  TString fwd(gSystem->Getenv("QA_FWD"));
  TString mac(gROOT->GetMacroPath());
  if (!fwd.IsNull()) {
    mac.Prepend(Form(".:%s:",fwd.Data()));
    gSystem->AddIncludePath(Form("-I%s", fwd.Data()));
  }
  else { 
    fwd = gSystem->Getenv("ANA_SRC");
    if (fwd.IsNull()) 
      fwd = "${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2";
    mac.Prepend(Form(".:%s/qa:%s/corrs:",fwd.Data(), fwd.Data()));
    gSystem->AddIncludePath(Form("-I%s/qa", fwd.Data()));
  }
  gROOT->SetMacroPath(mac);
  gSystem->Load("libGpad");
  gSystem->Load("libTree");

  gROOT->LoadMacro("QAPlotter.C+g");
  QAPlotter p(type, year, period, pass, useVar);
  p.AddFile(input);
  // t.SetOutputName("trending.root");
  if (!p.Run()) exit(1);
}
//
// EOF
//
