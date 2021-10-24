#include "TString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "/home/martin2/Documents/macros/SystematicErrorEstimation.C"

enum modes { kPt = 0, kZ = 1, kXi = 2};


TString getStringFromTObjStrArray(TObjArray *arr, Int_t position)
{
  if (position > arr->GetEntriesFast()-1)
    return "";
  
  TObjString* objStr = (TObjString*)(arr->At(position));
  if (!objStr)
    return "";
  
  return objStr->GetString();
}

//_____________________________________________________________________________________________________________________________
Bool_t CentralityHasDecimalsPlaces(Double_t cent)
{
  // Check, whether the centrality has decimals places or is an integer
  const Double_t temp1 = ((Int_t)cent)*1e6;
  const Double_t temp2 = cent*1e6;
  
  return TMath::Abs(temp1 - temp2) > 1;
}


//_____________________________________________________________________________________________________________________________
void setMeanThreshold(Double_t& setMeanLowerThreshold, Double_t& setMeanUpperThreshold, Int_t mode, Double_t lowerJetPt, Double_t upperJetPt)
{
  const Double_t setMeanLowerThresholdPt = 4;
  const Double_t setMeanUpperThresholdPt = 999.;
  
  // Take effective jetPt for z and xi
  Double_t effectiveJetPt = 0.5 * (lowerJetPt + upperJetPt);
  
  if (mode == kPt) {
    setMeanLowerThreshold = setMeanLowerThresholdPt;
    setMeanUpperThreshold = setMeanUpperThresholdPt;
  }
  else if (mode == kZ) {
    setMeanLowerThreshold = setMeanLowerThresholdPt / effectiveJetPt;
    setMeanUpperThreshold = setMeanUpperThresholdPt / effectiveJetPt;
  }
  else if (mode == kXi) {
    // Thresholds are swapped!
    setMeanLowerThreshold = TMath::Log(effectiveJetPt / setMeanUpperThresholdPt);
    setMeanUpperThreshold = TMath::Log(effectiveJetPt / setMeanLowerThresholdPt);
    if (setMeanUpperThresholdPt > 900)
      setMeanLowerThreshold = 0.;
  }
  else {
    printf("ERROR: Unknown mode!\n");
    exit(-1);
  }
  
}
/* "", "_negCharge", "_posCharge"*/
void runSystematicErrorEstimation_new(TString jetString = "", TString chargeString = "", TString referencepath = "", TString outfilepath = "", TString referenceFileSchema = "", TString stringForSystematicInformation = "", TString centStepsString = "0;100;0;10;60;100", TString jetPtStepsString = "5;10;10;15;15;20;20;30", TString modesInputString = "Pt", Int_t nSigma = 0)
{
  // NOTE:
  // If mean instead of default is to be used, just set setMean to kTRUE in the following.
  // ALSO: In addUpSystematicErrors the reference must be adapted to the corresponding mean (not the default).
  // Currently, ONLY the graph and the plot with the single outputs for the sys errors are changed.
  // TODO for average instead of default for spline systetematics at high pT > 4 GeV/c:
  // - Seems to make no difference more or less (at least for the old "symmetric" systetematics)
  // - Questionable to set the threshold as function of momentum! Could mean that K and p are move away from the best estimated
  //   because only the different models only affect pions in the rel. rise but for quite some pT not the other species
  // - Also a problem:
  // ATTENTION: Normally, weighted average should be done either for the yields or should be WEIGHTED with the total yield
  // (must be done anyway for the to-pion-ratios!). BUT: The whole approach relies on the assumption that the overall statistics
  // is the same for all compared data points. This must be checked anyway and is true within percent typically!
  //
  // => Looks like much to do and many problems with most likely little gain. Better don't use the average instead of the default!
	printf("Start systematic error estimation\n");
  
  const Bool_t ignoreSigmaErrors = kTRUE;
  
  TObjArray* jetPtBins = jetPtStepsString.Tokenize(";");
	Int_t numJetPtBins = jetPtBins->GetEntriesFast();
	Int_t* jetPt = new Int_t[numJetPtBins];
	for (Int_t i=0;i<numJetPtBins;i++) {
		jetPt[i] = getStringFromTObjStrArray(jetPtBins, i).Atoi();
	}
	
  TString jetPtString = "";
	
  TObjArray* centBins = centStepsString.Tokenize(";");
	Int_t numCentralities = centBins->GetEntriesFast();
	Int_t* centralities = new Int_t[numCentralities];
	for (Int_t i=0;i<numCentralities;i++) {
		centralities[i] = getStringFromTObjStrArray(centBins, i).Atoi();
	}
  
  TString centralityString = "";
  
  TObjArray* modeArray = modesInputString.Tokenize(";");
	Int_t maxModes = modeArray->GetEntriesFast();
	TString* modeString = new TString[maxModes];
	for (Int_t i=0;i<maxModes;i++) {
    modeString[i] = getStringFromTObjStrArray(modeArray, i);
	}
	
  TObjArray* systematicsStringsArray = stringForSystematicInformation.Tokenize("|");
  
	Int_t nOfSystematics = systematicsStringsArray->GetEntriesFast();
  
  TObjArray* systematicsArray = new TObjArray(nOfSystematics);
	for (Int_t i=0;i<nOfSystematics;i++) {
    TObjArray* singleSystematicsStringArray = getStringFromTObjStrArray(systematicsStringsArray, i).Tokenize("§");
    
    TObjArray* singleSystematicsArray = new TObjArray(singleSystematicsStringArray->GetEntriesFast());

    TObjArray* singleSystematicsSetting = getStringFromTObjStrArray(singleSystematicsStringArray, 0).Tokenize("$");
    
    singleSystematicsArray->AddAt(singleSystematicsSetting, 0); //Adding reference histogram title, outfile and input path at position 0
    
    //Looping through systematic files
    for (Int_t j=1;j<singleSystematicsStringArray->GetEntriesFast();++j) {
      TObjArray* singleSysEntryArray = getStringFromTObjStrArray(singleSystematicsStringArray, j).Tokenize("$");
      singleSystematicsArray->AddAt(singleSysEntryArray, j);
    }
    
    systematicsArray->AddAt(singleSystematicsArray, i);
	}
  
  Double_t setMeanLowerThreshold = 999.;
  Double_t setMeanUpperThreshold = 999.;

  Bool_t useCentralities = kTRUE;
  Bool_t useJetPt = !(jetString.Contains("Inclusive")) && jetString.Contains("Jets");
  if (!useJetPt)
    numJetPtBins = 2;
  
  Bool_t setMean = kFALSE;

  // Inclusive means: Only pT // In case of taking the average of the systematics: Only do this within certain range, HERE: pT only
  if (useJetPt) 
    setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, -1, -1);
  
  for (Int_t iCent = 0; iCent < numCentralities/2; iCent++) {
    centralityString = useCentralities ? Form("_centrality%d_%d", centralities[2*iCent], centralities[2*iCent + 1]) : ""; 

    for (Int_t iJetPt = 0; iJetPt < numJetPtBins/2; iJetPt++)  {
      jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[2*iJetPt], jetPt[2*iJetPt + 1]) : "";
      jetPtString.Append(chargeString);        
      for (Int_t mode = 0; mode < maxModes; mode++) {  
        TString referenceFile = Form(referencepath + "/" + referenceFileSchema, jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
        for (Int_t sysEntry=0;sysEntry<systematicsArray->GetEntriesFast();sysEntry++) {
          TObjArray* singleSystematicsArray = (TObjArray*)systematicsArray->At(sysEntry); //Contains all systematic entries
          TObjArray* generalProperties = (TObjArray*)singleSystematicsArray->At(0);       //Contains output Pattern and reference histTitle
          
          Int_t numFiles = singleSystematicsArray->GetEntriesFast();
          
          TString* fileNames = new TString[numFiles];
          TString* histTitles = new TString[numFiles];
          
          fileNames[0] = referenceFile;
          histTitles[0] = getStringFromTObjStrArray(generalProperties, 1);   //Reference hist Title
          TString inputPath = getStringFromTObjStrArray(generalProperties, 2);
          
          for (Int_t singleSysEntry=1;singleSysEntry<numFiles;singleSysEntry++) {
            TObjArray* singleSysEntryArray = (TObjArray*)singleSystematicsArray->At(singleSysEntry);
            TString fileString = getStringFromTObjStrArray(singleSysEntryArray, 0);
            fileNames[singleSysEntry] = inputPath + "/" + Form(fileString, jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
            histTitles[singleSysEntry] = getStringFromTObjStrArray(singleSysEntryArray, 1);
          }
          
          TString outFilePattern = getStringFromTObjStrArray(generalProperties, 0);
          TString outFileTitle = Form(outFilePattern, jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          
          //TODO: Last parameter should be mode, this has to be adjusted (needed for titles)
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold);  
        }        
      }
    }
  }
}
