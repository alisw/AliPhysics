#include "TString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "/home/martin2/Documents/macros/AddUpSystematicErrors.C"

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
void runAddUpSystematicErrors_new(
  TString jetString,
  TString chargeString,
  TString date,
  TString referenceFile,
  TString centStepsString,
  TString jetPtStepsString,
  TString modesInputString,
  TString outPutFilesToAddUp,
  TString outPath,
  Int_t nSigma
)
{
  const TString sigmaString = Form("_nSigma%.1f", (Double_t)nSigma);// nSigma > 0 ? Form("_nSigma%.1f", nSigma) : "";

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
  
  Bool_t useCentralities = kTRUE;
  Bool_t useJetPt = !(jetString.Contains("Inclusive")) && jetString.Contains("Jets");
  if (!useJetPt)
    numJetPtBins = 2;
  
  TObjArray* modeArray = modesInputString.Tokenize(";");
	Int_t maxModes = modeArray->GetEntriesFast();
	TString* modeString = new TString[maxModes];
	for (Int_t i=0;i<maxModes;i++) {
    modeString[i] = getStringFromTObjStrArray(modeArray, i);
	}
	
  TObjArray* outPutFilesArray = outPutFilesToAddUp.Tokenize("|");
	Int_t numFiles = outPutFilesArray->GetEntriesFast();
	TString* filePatterns = new TString[numFiles];
	for (Int_t i=0;i<numFiles;i++) {
    filePatterns[i] = getStringFromTObjStrArray(outPutFilesArray, i);
	}
	
	TString* fileNames = new TString[numFiles];

  for (Int_t iCent = 0; iCent <= numCentralities/2 - 1; iCent++) {
    centralityString = useCentralities ? Form("_centrality%d_%d", centralities[2*iCent], centralities[2*iCent + 1]) : ""; 

    for (Int_t iJetPt = 0; iJetPt <= numJetPtBins/2 - 1; iJetPt++)  {
      jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[2*iJetPt], jetPt[2*iJetPt + 1]) : "";
      jetPtString.Append(chargeString);        
      for (Int_t mode = 0; mode <= maxModes - 1; mode++) {
        TString outFileTitle = Form("SummedSystematicErrors_%s_%s__%s%s__%s", jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), date.Data());
        
        for (Int_t iSyst = 0;iSyst < numFiles; ++iSyst) {
          fileNames[iSyst] = Form(filePatterns[iSyst] + sigmaString + "__" + date + ".root", jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
        }

        TString fileNameReference = Form(referenceFile, jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
        
        AddUpSystematicErrors(outPath, outFileTitle, fileNames, numFiles, fileNameReference);
      }
    }
  }
}
