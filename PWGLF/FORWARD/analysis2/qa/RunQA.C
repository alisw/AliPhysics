/**
 * @file   RunQA.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:35:08 2011
 * 
 * @brief  Script to run the QATrender over a single run
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
class TTree;
class TCanvas;

/** 
 * Run the QATrender
 * 
 * The QATrender is run over the list of files (runs) to produce the
 * file <tt>trending.root</tt> which contains a TTree of QA
 * information - one entry per run.  
 * 
 * The QATrender will also produce two files per run: 
 * 
 * - <tt>index.root</tt> which contains TCanvas objects of
 *   the finished plots.
 *
 * - <tt>index.pdf</tt> which is a PDF of the TCanvases
 *   mentioned above.
 *
 * @param input  Input file
 * @param type   Data type (data or sim)
 * @param year   Year 
 * @param period Period (e.g., LHC10h)
 * @param pass   Pass (e.g., pass1, cpass1, passMC)
 * @param runNo  Run number 
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void
RunQA(const char* input, 
      const char* type, 
      Int_t       year, 
      const char* period,
      const char* pass, 
      Long_t      runNo) 
{
  Bool_t keep = true;
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

  TString inp(input);
  if (inp.Contains("#QA_results")) inp.ReplaceAll("#QA_results", "#trending");
  
  gROOT->LoadMacro("QATrender.C+g");
  QATrender t(keep, type, year, period, pass, runNo);
  t.AddFile(inp);
  // t.SetOutputName("trending.root");
  if (!t.Run()) exit(1);
}
//
// EOF
//
