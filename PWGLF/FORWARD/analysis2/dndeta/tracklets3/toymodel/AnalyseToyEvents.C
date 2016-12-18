
/**
 * @file   AnalyseToyEvents.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Aug 24 22:27:43 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets_toy
 * 
 */

/** 
 * Run the @f$ \mathrm{d}N/\mathrm{d}\vartheta@f$ analysis 
 * 
 * @param realFileName The measured tracklets 
 * @param simFileName  The tracklets used for corrections
 * @param weight       Reweighing factor 
 *
 * @relates ToyModeldNdTheta
 */
void
AnalyseToyEvents(const char* realFileName="real.root",
		 const char* simFileName="reduced.root",
		 Double_t    weight=1)
{
  TString dir("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
  if (gSystem->Getenv("ANA_SRC")) dir="$ANA_SRC";

  if (!gROOT->GetClass("ToyModeldNdTheta"))
    gROOT->LoadMacro(Form("%s/dndeta/tracklets3/toymodel/ToyModeldNdTheta.C+g",
			  dir.Data()));

  ToyModeldNdTheta* ana = new ToyModeldNdTheta;
  ana->Run(realFileName, simFileName, weight);
  ana->Visualize();
}

//
// EOF
//
