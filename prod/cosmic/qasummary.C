void qasummary()
{
  // checks all QA result and issues a summary

  TString cmd = Form(".! ls QA.*.root > tempo.txt") ;
  gROOT->ProcessLine(cmd.Data()) ; 		   
  ifstream in("tempo.txt") ; 	
  TString file ;
  while ( 1 ) {
    in >> file ; 
    if ( !in.good() ) 
      break ; 
    TFile qaf(file) ;
    Bool_t flag = kFALSE ; 
    AliQA * qa = static_cast<AliQA *>(qaf.Get(AliQA::GetQAName())) ; 
    for (Int_t det = 0; det < AliQA::kNDET; det++) {
      for (Int_t task = 0; task < AliQA::kNTASK; task++) {
	for (Int_t bit = 0; bit < AliQA::kNBIT; bit++) {
	  if (qa->IsSet(det, task, bit)) {
	    flag = kTRUE ; 
	    break ;
	  }
	  if (flag)
	    break ; 
	}
	if (flag) break ; 
      }
      if (flag) {
	printf(" ***************** QA error found in %s :\n", file.Data()) ; 
	qa->Show(det) ;
	flag = kFALSE ; 
      } 
    }
  }  
}
