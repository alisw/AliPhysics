///
/// \file V0LamAnalysis/AliAnalysisV0LamCutProcessing.cxx
///

#include "AliAnalysisV0LamCutProcessing.h"

using namespace std;

AliAnalysisV0LamCut::AliAnalysisV0LamCut():
  fNumberOfCutValues(0),
  fIsAnUpperLimit(false),
  fCutValues(0)
{
}


AliAnalysisV0LamCut::AliAnalysisV0LamCut(std::vector<Double_t> variableCutValues, bool isUpperBound):
  fIsAnUpperLimit(isUpperBound)
{
  //Contructor for the cut class
  fNumberOfCutValues = variableCutValues.size();
  fCutValues = variableCutValues;
}

AliAnalysisV0LamCut::~AliAnalysisV0LamCut()
{
  //default destructor
}



AliAnalysisV0LamCutProcessing::AliAnalysisV0LamCutProcessing(TList *const outputList, Int_t varCutIndex):
  fOutputList(outputList) // use main task's output list
{
  //Constructor for the cut processing class.
  //Sets up the cuts that will be used for V0 reconstruction.  Also initializes
  //output histograms that will store information about the V0s.


  fVariableCutType = varCutIndex;
  if(fVariableCutType > 0) {
    fNumberOfVariableCutValues = 3;//needs to be manually set if values change
  }
  else fNumberOfVariableCutValues = 1;

  //New recontruction cuts should be put in here.
  // fOutputList = outputList; //Get the output list from the main analysis task

  vector<Double_t> valueDCAPrimProton;
  vector<Double_t> valueDCAPrimPion;
  vector<Double_t> valueDCADaughters;
  vector<Double_t> valueProperDecayLength;
  vector<Double_t> valueEta;
  vector<Double_t> valueCosPointing;
  vector<Double_t> valueV0DCA;
  vector<Double_t> valuePt;
  vector<Double_t> valueMassLamDiff;
  vector<Double_t> valueMassALamDiff;

  // Set all variable cut values, or set just the nominal value
  if(kDCAProton == fVariableCutType) {
    valueDCAPrimProton.push_back(0.09);
    valueDCAPrimProton.push_back(0.1);
    valueDCAPrimProton.push_back(0.11);
  } else {
    valueDCAPrimProton.push_back(0.1);
  }
  if(kDCAPion == fVariableCutType) {
    valueDCAPrimPion.push_back(0.27);
    valueDCAPrimPion.push_back(0.3);
    valueDCAPrimPion.push_back(0.33);
  } else {
    valueDCAPrimPion.push_back(0.3);
  }
  if(kDCADaughter == fVariableCutType) {
    valueDCADaughters.push_back(0.36);
    valueDCADaughters.push_back(0.4);
    valueDCADaughters.push_back(0.44);
  } else {
    valueDCADaughters.push_back(0.4);
  }
  if(kProperDecayLength == fVariableCutType) {
    valueProperDecayLength.push_back(54);
    valueProperDecayLength.push_back(60);
    valueProperDecayLength.push_back(66);
  } else {
    valueProperDecayLength.push_back(60);
  }
  if(kEta == fVariableCutType) {
    valueEta.push_back(0.72);
    valueEta.push_back(0.8);
    valueEta.push_back(0.88);
  } else {
    valueEta.push_back(0.8);
  }
  if(kCosPointing == fVariableCutType) {
    valueCosPointing.push_back(0.99923);
    valueCosPointing.push_back(0.9993);
    valueCosPointing.push_back(0.99937);
  } else {
    valueCosPointing.push_back(0.9993);
  }
  if(kDCAV0 == fVariableCutType) {
    valueV0DCA.push_back(0.45);
    valueV0DCA.push_back(0.5);
    valueV0DCA.push_back(0.55);
  } else {
    valueV0DCA.push_back(0.5);
  }
  if(kPt == fVariableCutType) {
    valuePt.push_back(0.36);
    valuePt.push_back(0.4);
    valuePt.push_back(0.44);
  } else {
    valuePt.push_back(0.4);
  }
  if(kMassLam == fVariableCutType) {
    valueMassLamDiff.push_back(0.00342);
    valueMassLamDiff.push_back(0.0038);
    valueMassLamDiff.push_back(0.00418);
  } else {
    valueMassLamDiff.push_back(0.0038);
  }
  if(kMassALam == fVariableCutType) {
    valueMassALamDiff.push_back(0.00342);
    valueMassALamDiff.push_back(0.0038);
    valueMassALamDiff.push_back(0.00418);
  } else {
    valueMassALamDiff.push_back(0.0038);
  }

  // //Make cut objects for each cut type
  fCuts.push_back(AliAnalysisV0LamCut(valueDCAPrimProton, false)); //v0->daughterPosDCAPrimaryVertex;
  fCuts.push_back(AliAnalysisV0LamCut(valueDCAPrimPion, false)); //v0->daughterNegDCAPrimaryVertex;
  fCuts.push_back(AliAnalysisV0LamCut(valueDCADaughters, true)); //v0->daughtersDCA;
  fCuts.push_back(AliAnalysisV0LamCut(valueProperDecayLength, true)); //v0->decayLength;
  fCuts.push_back(AliAnalysisV0LamCut(valueEta, true)); //v0->v0Eta;
  fCuts.push_back(AliAnalysisV0LamCut(valueCosPointing, false)); //v0->cosPointing;
  fCuts.push_back(AliAnalysisV0LamCut(valueV0DCA, true)); //v0->v0DCA;
  fCuts.push_back(AliAnalysisV0LamCut(valuePt, false)); //v0->Pt;
  fCuts.push_back(AliAnalysisV0LamCut(valueMassLamDiff, true)); //v0->massLamDifference;
  fCuts.push_back(AliAnalysisV0LamCut(valueMassALamDiff, true)); //v0->massALamDifference;

  fNumberOfCutTypes = fCuts.size(); //Manually set this to the number of recon cuts

  InitHistograms();
}

AliAnalysisV0LamCutProcessing::AliAnalysisV0LamCutProcessing(const AliAnalysisV0LamCutProcessing &processor)
{ // Not implemented
}

AliAnalysisV0LamCutProcessing AliAnalysisV0LamCutProcessing::operator=(const AliAnalysisV0LamCutProcessing &processor)
{ // Not implemented
  return *this;
}

AliAnalysisV0LamCutProcessing::~AliAnalysisV0LamCutProcessing()
{
}

void AliAnalysisV0LamCutProcessing::ProcessCut(AliReconstructedV0 *v0, int index, bool isLambdaCandidate)
{
  //This function processes a single reconstruction cut of a V0.  If the cut in
  //question is a variable cut, the function loops over each value of that cut
  //and stores whether the V0 passes or fails each particular value.
  //If isLambdaCandidate is false, the V0 is an Antilambda candidate

  //First we get the relevant V0 information that will be tested by this cut.
  Double_t v0Value = 1000;

  switch (index) {
    case 0: //Get the proton DCA to primary
      v0Value = (isLambdaCandidate)
              ? v0->daughterPosDCAPrimaryVertex
              : v0->daughterNegDCAPrimaryVertex;
      break;
    case 1: //Get the pion DCA to primary
      v0Value = (isLambdaCandidate)
              ? v0->daughterNegDCAPrimaryVertex
              : v0->daughterPosDCAPrimaryVertex;
      break;
    case 2:
      v0Value = v0->daughtersDCA;
      break;
    case 3:
      // something is wrong if gamma <=0, so set to a fail val
      if (v0->lorentzGammaLam > 0) {
        v0Value = v0->decayLength / v0->lorentzGammaLam;
      }
      break;
    case 4:
      v0Value = v0->v0Eta;
      break;
    case 5:
      v0Value = v0->cosPointing;
      break;
    case 6:
      v0Value = v0->v0DCA;
      break;
    case 7:
      v0Value = v0->v0Pt;
      break;
    case 8:
      v0Value = v0->massLamDifference;
      break;
    case 9:
      v0Value = v0->massALamDifference;
      break;
    default:
      cerr << "ERROR: No cut for index value " << index << "\n";
  }

  //Find how many cut values are associated with this cut
  UInt_t numberOfCutValues = fCuts[index].fNumberOfCutValues;
  if (0 == numberOfCutValues) {
    cerr<<"ERROR: Must have at least one cut value \n";
    return;
  }

  for (UInt_t i = 0; i < numberOfCutValues; i++) {
    //check if v0 passes each cut value and set v0 bools accordingly
    //Cut can either be an upper bound or a lower bound.
    if(fCuts[index].fIsAnUpperLimit){ //upper bound cuts
      if (fCuts[index].fCutValues[i] > v0Value) {
        //v0 passed cut.  Set v0 bool accordingly
        v0->hasPassedCut[index][i] = kTRUE;
      }
    }
    else { //lower bound cuts
      if(fCuts[index].fCutValues[i] < v0Value) {
        //v0 passed cut.  Set v0 bool accordingly
        v0->hasPassedCut[index][i] = kTRUE;
      }
    }
  }

}

void AliAnalysisV0LamCutProcessing::CheckIfV0PassesCuts(AliReconstructedV0 *v0)
{
  // Called by the Analysis Task
  // Function which checks if V0 passes cuts.
  // First the code processes each cut.
  // If the candidate passes all the cuts, DetermineIfTrueV0 sets
  // isLamCenter or isALamCenter to true

  for(int cutTypeIndex = 0; cutTypeIndex < fNumberOfCutTypes; cutTypeIndex++)
  {
    //First default these to false
    v0->isLamCenter[cutTypeIndex]      = kFALSE;
    v0->isALamCenter[cutTypeIndex]     = kFALSE;
    v0->isDeemedUnworthy[cutTypeIndex] = kFALSE;
    for(int varCutIndex = 0; varCutIndex < fNumberOfVariableCutValues; varCutIndex++){
      v0->hasPassedCut[cutTypeIndex][varCutIndex]=kFALSE;
    }
  }
  //At this point, a V0 could simultaneously be a lambda candidate and an
  //antilambda candidates (based only on the PID of the daughters).  It is
  //necessary to check both cases separately.

  //Process cuts for lambda candidates
  if(v0->hasProtonDaughter && v0->hasPiMinusDaughter){
    bool isLambdaCandidate = kTRUE;
    for(int cutIndex = 0; cutIndex < 8; cutIndex++) ProcessCut(v0,cutIndex,isLambdaCandidate); //General cuts
    ProcessCut(v0,8,isLambdaCandidate); //Lambda specific cut
    AliAnalysisV0LamCutProcessing::DetermineIfTrueV0(v0,isLambdaCandidate);
  }
  //Process cuts for antilambda candidates
  if(v0->hasAntiProtonDaughter && v0->hasPiPlusDaughter){
    bool isLambdaCandidate = kFALSE;
    for(int cutIndex = 0; cutIndex < 8; cutIndex++) ProcessCut(v0,cutIndex,isLambdaCandidate); //General cuts
    ProcessCut(v0,9,isLambdaCandidate); //AntiLambda specific cut
    AliAnalysisV0LamCutProcessing::DetermineIfTrueV0(v0,isLambdaCandidate);
  }
}

void AliAnalysisV0LamCutProcessing::DetermineIfTrueV0(AliReconstructedV0 *v0, bool isLambda)
{
  //Checks to see if the V0 passed all the standard cuts and any of the
  //variable cuts.  The function either checks if the v0 look like a Lambda
  //or looks like an antiLambda, as set by the isLambda input variable.
  vector<bool> passedCut(fNumberOfVariableCutValues,kTRUE);
  bool hasFailedStandardCut = kFALSE;
  //Check which cuts the V0 passed and failed
  for(int cutIndex = 0; cutIndex < fNumberOfCutTypes; cutIndex++){
    int numberOfCutValues = fCuts[cutIndex].fNumberOfCutValues; //find how many cut values
    for(int i = 0; i < numberOfCutValues; i++){
      if(8 > cutIndex /* update this if number of cuts changes*/){ //check generic V0 cuts
	if(!v0->hasPassedCut[cutIndex][i]){
	  if(1 == numberOfCutValues) hasFailedStandardCut = kTRUE;
	  else passedCut[i] = kFALSE;
	}
      }
      if((8 == cutIndex) && isLambda){ //Only check for Lambda case
	if(!v0->hasPassedCut[cutIndex][i]){
	  if(1 == numberOfCutValues) hasFailedStandardCut = kTRUE;
	  else passedCut[i] = kFALSE;
	}
      }
      if((9 == cutIndex) && !isLambda){ //Only check for AntiLambda case
	if(!v0->hasPassedCut[cutIndex][i]){
	  if(1 == numberOfCutValues) hasFailedStandardCut = kTRUE;
	  else passedCut[i] = kFALSE;
	}
      }
    }
  }
  if(hasFailedStandardCut){
    //if it failed any standard (i.e. non-variable cuts) cuts, set all
    //passedCut values to failed
    for(int i = 0; i < fNumberOfVariableCutValues; i++){
      passedCut[i] = kFALSE;
    }
  }
  for(int i = 0; i < fNumberOfVariableCutValues; i++){
    //Now, for each candidate that passed cuts, determine if Lam or antiLam
    if(passedCut[i]){ //if this is true at this point, the v0 has passed ALL the
                      //standard cuts AND the ith variable cut
      if(isLambda) v0->isLamCenter[i]=kTRUE;
      else v0->isALamCenter[i]=kTRUE;
    } //Otherwise isLamCenter and isALamCenter are false by default
  }
}




void AliAnalysisV0LamCutProcessing::DoV0Histogramming(AliReconstructedV0 *v0)
{
  //Called by the analysis task.  Histograms the V0's reconstruction parameters.
  //Does separate histogramming for lambda and antilambda cases
  if(v0->hasProtonDaughter && v0->hasPiMinusDaughter){
    SortAndFillCutHistograms(v0, true);
  }
  if(v0->hasAntiProtonDaughter && v0->hasPiPlusDaughter){
    SortAndFillCutHistograms(v0,false);
  }
}

void AliAnalysisV0LamCutProcessing::SortAndFillCutHistograms(AliReconstructedV0 *v0, bool isLambda)
{
  //Loops through each cut type to attempt to histogram it.  If all the *other*
  //cuts pass, then the V0's value for this cut is histogrammed.  FillHist is
  //called to handle the filling of the histogram
  //e.g. Only fill Minv hist if V0 passed decaylength, DCA, etc.
  for(int cutTypeIndex = 0; cutTypeIndex < fNumberOfCutTypes; cutTypeIndex++)
  { //Loop over each cut type
    bool passesOtherCuts = true;
    //Ignore irrelevant species specific cuts
    if(isLambda && (9 == cutTypeIndex)) continue; //ignore antilambda cut
    if(!isLambda && (8 == cutTypeIndex)) continue; //ignore lambda cut
    vector<bool> passesOtherVariableCuts(fNumberOfVariableCutValues, true);
    for(int variableCutIndex = 0; variableCutIndex < fNumberOfVariableCutValues; variableCutIndex++)
    { //Loop over number of variable cut values
      if(!passesOtherCuts) break; //If it fails a cut, don't need to histogram it
      for(int otherCutTypeIndex =0; otherCutTypeIndex < fNumberOfCutTypes; otherCutTypeIndex++)
      {
	if(cutTypeIndex == otherCutTypeIndex) continue; //only need to pass the OTHER cuts
	//Again, ignore irrelevant species specific cuts
	if(isLambda  && (9 == otherCutTypeIndex))  continue;
	if(!isLambda && (8 == otherCutTypeIndex))  continue;
	if(fCuts[otherCutTypeIndex].fNumberOfCutValues > 1){
	  if(!v0->hasPassedCut[otherCutTypeIndex][variableCutIndex]) {
	    passesOtherVariableCuts[variableCutIndex] = false;
	    break;
	  }
	}
	else if(!v0->hasPassedCut[otherCutTypeIndex][0]){
	  passesOtherCuts = false;
	  break;
	}
      }
      //Finally, if it passes all standard cuts and the current variable cut,
      //fill the histogram for that cut type.
      if(passesOtherCuts && passesOtherVariableCuts[variableCutIndex]){
	FillHist(v0,cutTypeIndex,variableCutIndex,isLambda);
      }
    }
  }
}

void AliAnalysisV0LamCutProcessing::InitHistograms()
{
  //Initilizes histograms for V0 characteristics
  TString nameDaughterPosDcaToPrimLam = "fHistDaughterPosDcaToPrimLam";
  TString nameDaughterPosDcaToPrimALam = "fHistDaughterPosDcaToPrimALam";
  TString nameDaughterNegDcaToPrimLam = "fHistDaughterNegDcaToPrimLam";
  TString nameDaughterNegDcaToPrimALam = "fHistDaughterNegDcaToPrimALam";
  TString nameDaughtersDcaLam = "fHistDaughtersDcaLam";
  TString nameDaughtersDcaALam = "fHistDaughtersDcaALam";
  TString nameDecayLengthLam = "fHistDecayLengthLam";
  TString nameDecayLengthALam = "fHistDecayLengthALam";
  TString nameProperDecayLengthLam = "fHistProperDecayLengthLam";
  TString nameProperDecayLengthALam = "fHistProperDecayLengthALam";
  TString nameEtaLam = "fHistEtaLam";
  TString nameEtaALam = "fHistEtaALam";
  TString nameCosPointingLam = "fHistCosPointingLam";
  TString nameCosPointingALam = "fHistCosPointingALam";
  TString nameDcaLam = "fHistDcaLam";
  TString nameDcaALam = "fHistDcaALam";
  TString namePtLam = "fHistPtLam";
  TString namePtALam = "fHistPtALam";
  TString nameMassLam = "fHistMassLam";
  TString nameMassALam = "fHistMassALam";
  TString nameMassCentralityLam = "fHistMassCentralityLam";
  TString nameMassCentralityALam = "fHistMassCentralityALam";
  fHistDaughterPosDcaToPrimLam = new TH2F((nameDaughterPosDcaToPrimLam).Data(),"DCA of proton daughter to primary vertex;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 10.);
  fHistDaughterPosDcaToPrimALam = new TH2F((nameDaughterPosDcaToPrimALam).Data(),"DCA of piplus daughter to primary vertex;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 10.);
  fHistDaughterNegDcaToPrimLam = new TH2F((nameDaughterNegDcaToPrimLam).Data(),"DCA of piminus daughter to primary vertex;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 10.);
  fHistDaughterNegDcaToPrimALam = new TH2F((nameDaughterNegDcaToPrimALam).Data(),"DCA of antiproton daughter to primary vertex;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 10.);
  fHistDaughtersDcaLam = new TH2F((nameDaughtersDcaLam).Data(), "DCA of Lambda daughters to each other;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 1.5);
  fHistDaughtersDcaALam = new TH2F((nameDaughtersDcaALam).Data(), "DCA of AntiLambda daughters to each other;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 1.5);
  fHistDecayLengthLam = new TH2F((nameDecayLengthLam).Data(), "Lambda decay length;Cut Bin;Decay length (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 500, 0., 100.);
  fHistDecayLengthALam = new TH2F((nameDecayLengthALam).Data(), "AntiLambda decay length;Cut Bin;Decay length (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 500, 0., 100.);
  fHistProperDecayLengthLam = new TH2F((nameProperDecayLengthLam).Data(), "Lambda proper decay length;Cut Bin;Proper decay length (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 500, 0., 100.);
  fHistProperDecayLengthALam = new TH2F((nameProperDecayLengthALam).Data(), "AntiLambda proper decay length;Cut Bin;Proper decay length (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 500, 0., 100.);
  fHistEtaLam = new TH2F((nameEtaLam).Data(), "Eta distribution of Lambdas;Cut Bin;Eta;Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 200, -1., 1.);
  fHistEtaALam = new TH2F((nameEtaALam).Data(), "Eta distribution of AntiLambdas;Cut Bin;Eta;Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 200, -1., 1.);
  fHistCosPointingLam = new TH2F((nameCosPointingLam).Data(), "Cosine(pointing angle) of Lambdas;Cut Bin;Cos(#Theta_{p});Counts",  fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,100, .97, 1.);
  fHistCosPointingALam = new TH2F((nameCosPointingALam).Data(), "Cosine(pointing angle) of AntiLambdas;Cut Bin;Cos(#Theta_{p});Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, .97, 1.);
  fHistDcaLam = new TH2F((nameDcaLam).Data(), "DCA of Lambdas to primary vertex;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 5.);
  fHistDcaALam = new TH2F((nameDcaALam).Data(), "DCA of AntiLambdas to primary vertex;Cut Bin;DCA (cm);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 100, 0., 5.);
  fHistPtLam = new TH2F((namePtLam).Data(), "Lambda pT distribution;Cut Bin;p_{T} (GeV/c);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 500, 0., 10.);
  fHistPtALam = new TH2F((namePtALam).Data(), "AntiLambda pT distribution;Cut Bin;p_{T} (GeV/c);Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 500, 0., 10.);
  fHistMassLam = new TH2F((nameMassLam).Data(), "Lambda Minv;Cut Bin;m_{inv} (GeV/c^{2});Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 600,1.04,1.19);
  fHistMassALam = new TH2F((nameMassALam).Data(), "AntiLambda Minv;Cut Bin;m_{inv} (GeV/c^{2});Counts", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 600,1.04,1.19);
  fHistMassCentralityLam = new TH3F((nameMassCentralityLam).Data(), "Lambda Minv vs Centrality;Cut Bin;Centrality Bin;m_{inv} (GeV/c^{2})", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 20, .5, 20+.5, 600, 1.04,1.19);
  fHistMassCentralityALam = new TH3F((nameMassCentralityALam).Data(), "AntiLambda Minv vs Centrality;Cut Bin;Centrality Bin;m_{inv} (GeV/c^{2})", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5, 20, .5, 20+.5, 600, 1.04,1.19);
  fOutputList->Add(fHistDaughterPosDcaToPrimLam);
  fOutputList->Add(fHistDaughterPosDcaToPrimALam);
  fOutputList->Add(fHistDaughterNegDcaToPrimLam);
  fOutputList->Add(fHistDaughterNegDcaToPrimALam);
  fOutputList->Add(fHistDaughtersDcaLam);
  fOutputList->Add(fHistDaughtersDcaALam);
  fOutputList->Add(fHistDecayLengthLam);
  fOutputList->Add(fHistDecayLengthALam);
  fOutputList->Add(fHistProperDecayLengthLam);
  fOutputList->Add(fHistProperDecayLengthALam);
  fOutputList->Add(fHistEtaLam);
  fOutputList->Add(fHistEtaALam);
  fOutputList->Add(fHistCosPointingLam);
  fOutputList->Add(fHistCosPointingALam);
  fOutputList->Add(fHistDcaLam);
  fOutputList->Add(fHistDcaALam);
  fOutputList->Add(fHistPtLam);
  fOutputList->Add(fHistPtALam);
  fOutputList->Add(fHistMassLam);
  fOutputList->Add(fHistMassALam);
  fOutputList->Add(fHistMassCentralityLam);
  fOutputList->Add(fHistMassCentralityALam);
  //Don't need to call PostData here, since PostData will be called by
  //the main analysis task
}

void AliAnalysisV0LamCutProcessing::FillHist(AliReconstructedV0 *v0, int cutTypeIndex, int variableCutIndex, bool isLambda)
{
  //Maps each cutTypeIndex to a different histogram and fills the associated value from the v0.
  //If isLambda is false, the V0 is an antilambda
  if(0 == cutTypeIndex){
    if(isLambda) fHistDaughterPosDcaToPrimLam->Fill(variableCutIndex, v0->daughterPosDCAPrimaryVertex);
    else fHistDaughterNegDcaToPrimALam->Fill(variableCutIndex, v0->daughterNegDCAPrimaryVertex);
  }
  if(1 == cutTypeIndex){
    if(isLambda) fHistDaughterNegDcaToPrimLam->Fill(variableCutIndex, v0->daughterNegDCAPrimaryVertex);
    else fHistDaughterPosDcaToPrimALam->Fill(variableCutIndex, v0->daughterPosDCAPrimaryVertex);
  }
  if(2 == cutTypeIndex){
    if(isLambda) fHistDaughtersDcaLam->Fill(variableCutIndex, v0->daughtersDCA);
    else fHistDaughtersDcaALam->Fill(variableCutIndex, v0->daughtersDCA);
  }
  if(3 == cutTypeIndex){
    if(isLambda){
      fHistDecayLengthLam->Fill(variableCutIndex, v0->decayLength);
      if(v0->lorentzGammaLam > 0){
	fHistProperDecayLengthLam->Fill(variableCutIndex, v0->decayLength/v0->lorentzGammaLam);
      }
    }
    else{
      fHistDecayLengthALam->Fill(variableCutIndex, v0->decayLength);
      if(v0->lorentzGammaLam > 0){
	fHistProperDecayLengthALam->Fill(variableCutIndex, v0->decayLength/v0->lorentzGammaLam);
      }
    }
  }
  if(4 == cutTypeIndex){
    if(isLambda) fHistEtaLam->Fill(variableCutIndex, v0->v0Eta);
    else fHistEtaALam->Fill(variableCutIndex, v0->v0Eta);
  }
  if(5 == cutTypeIndex){
    if(isLambda) fHistCosPointingLam->Fill(variableCutIndex, v0->cosPointing);
    else fHistCosPointingALam->Fill(variableCutIndex, v0->cosPointing);
  }
  if(6 == cutTypeIndex){
    if(isLambda) fHistDcaLam->Fill(variableCutIndex, v0->v0DCA);
    else fHistDcaALam->Fill(variableCutIndex, v0->v0DCA);
  }
  if(7 == cutTypeIndex){
    if(isLambda) fHistPtLam->Fill(variableCutIndex, v0->v0Pt);
    else fHistPtALam->Fill(variableCutIndex, v0->v0Pt);
  }
  if(8 == cutTypeIndex){
    if(isLambda){ //only do this histogramming for Lambdas
      fHistMassLam->Fill(variableCutIndex, v0->massLam);
      fHistMassCentralityLam->Fill(variableCutIndex, fCurrentCentralityBin+1,v0->massLam);
    }
  }
  if(9 == cutTypeIndex){
    if(!isLambda){ //Only do this histogramming for AntiLambdas
      fHistMassALam->Fill(variableCutIndex, v0->massALam);
      fHistMassCentralityALam->Fill(variableCutIndex, fCurrentCentralityBin+1,v0->massALam);
    }
  }
}
