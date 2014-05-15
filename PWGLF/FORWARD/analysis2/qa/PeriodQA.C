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
  gROOT->SetMacroPath(Form(".:$(ANA_SRC)/qa:$(ANA_SRC)/corrs:"
			   "$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/qa:"
			   "$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/corrs:"
			   "%s",
			   gROOT->GetMacroPath()));
  gSystem->AddIncludePath("-I${ALICE_ROOT}/PWGLF/FORWARD/analysis2/qa");
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
