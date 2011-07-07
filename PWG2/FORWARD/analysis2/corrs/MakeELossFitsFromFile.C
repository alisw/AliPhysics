/**
 * @file   MakeELossFitsFromFile.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Thu Jul  7 11:23:50 2011
 * 
 * @brief This script takes the output list from a file and applies
 * the energy fitter directly. In principle this should be done by the
 * train via AliFMDEnergyFitterTask::Terminate but if this fails or
 * impatience takes over this script can help...
 *
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
/** 
 * 
 * 
 * @param filename 
 * @param sys 
 * @param sNN 
 * @param field 
 * @param mc 
 *
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
void MakeELossFitsFromFile(const Char_t* filename="forward.root",
			   UShort_t sys=1, UShort_t sNN=900, Short_t field=5, 
			   Bool_t mc=false) 
{
  
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
  
  TFile::Open(filename);
  TList* list = (TList*)gFile->Get("Forward");
  
  AliForwardCorrectionManager::Instance().Init(sys, sNN, field, mc, 0);
  
  AliFMDEnergyFitter fitter("Fitter");
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  fitter.SetEtaAxis(200, -4, 6);
  // Set maximum energy loss to consider 
  fitter.SetMaxE(15); 
  // Set number of energy loss bins 
  fitter.SetNEbins(450);
  // Set whether to use increasing bin sizes 
  fitter.SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  fitter.SetDoFits(kTRUE);
  // Set whether to make the correction object 
  fitter.SetDoMakeObject(kTRUE);
  // Set the low cut used for energy
  fitter.SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  fitter.SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  fitter.SetNParticles(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  fitter.SetMinEntries(1000);
  fitter.SetMaxRelativeParameterError(0.12);
  fitter.SetMaxChi2PerNDF(20);
  fitter.SetMinWeight(1e-5);
  // --- Set limits on fits the energy -------------------------------
  // Maximum relative error on parameters 
  AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 20;

  TList* dummy = new TList;
  fitter.DefineOutput(dummy);
  TAxis a(200, -4, 6); // ignored, but needed
  fitter.Init(a); 

  fitter.Fit(list);

  list->ls();
  TDirectory* savdir = gDirectory;
  TFile* tmp = TFile::Open("elossfits.root", "RECREATE");
  list->Write(list->GetName(), TObject::kSingleKey);
  dummy->Write("dummy",  TObject::kSingleKey);
  tmp->Write();
  tmp->Close();
  savdir->cd();
}
//
// EOF
//
