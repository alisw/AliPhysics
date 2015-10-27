/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliMuonAnalysisOutput.h"

// ROOT includes
#include "TAxis.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "THashList.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TMath.h"

// STEER inlcudes
#include "AliLog.h"

// PWG includes
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"

/// \cond CLASSIMP
ClassImp(AliMuonAnalysisOutput) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliMuonAnalysisOutput::AliMuonAnalysisOutput() :
  TObject(),
  fCounterCollection(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0)
//,
//  fOutputPrototypeList(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliMuonAnalysisOutput::AliMuonAnalysisOutput ( TObjArray* outputList ) :
  TObject(),
  fCounterCollection(0x0),
  fMergeableCollection(0x0),
  fOutputList(outputList)
//,
//  fOutputPrototypeList(0x0)
{
  //
  /// Constructor
  //
  Init();
}

//________________________________________________________________________
AliMuonAnalysisOutput::AliMuonAnalysisOutput ( const char *filename, const char* outputName ) :
  TObject(),
  fCounterCollection(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0)
//,
//  fOutputPrototypeList(0x0)
{
  //
  /// Constructor
  //
  TFile* file = TFile::Open(filename);
  if ( file ) {
    fOutputList = static_cast<TObjArray*>(file->FindObjectAny(outputName));
    delete file;
  }
  else AliError(Form("Cannot open %s",filename));

  Init();
}

//________________________________________________________________________
AliMuonAnalysisOutput::~AliMuonAnalysisOutput()
{
  //
  /// Destructor
  //
//  delete fOutputPrototypeList;
}

//________________________________________________________________________
Bool_t AliMuonAnalysisOutput::Init()
{
  //
  /// Initialize objects
  //
  if ( ! fOutputList ) {
    AliError("No output list provided or found!");
    return kTRUE;
  }
  fCounterCollection = static_cast<AliCounterCollection*>(fOutputList->At(0));
  fMergeableCollection = static_cast<AliMergeableCollection*>(fOutputList->At(1));

  return kTRUE;
}


////________________________________________________________________________
//Bool_t AliMuonAnalysisOutput::AddObjectToCollection ( TObject* object )
//{
//  //
//  /// Add object to collection
//  //
//  
//  if ( ! fOutputPrototypeList ) {
//    fOutputPrototypeList = new THashList();
//    fOutputPrototypeList->SetOwner();
//  }
//  if ( fOutputPrototypeList->FindObject(object->GetName() ) ) {
//    AliWarning(Form("Object with name %s already in the list", object->GetName()));
//    return kFALSE;
//  }
//  fOutputPrototypeList->Add(object);
//  
//  return kTRUE;
//}

//________________________________________________________________________
TObject* AliMuonAnalysisOutput::GetMergeableObject ( TString physSel, TString trigClassName, TString centrality, TString objectName )
{
  //
  /// Get mergeable object
  /// (create collection if necessary)
  //

  if ( ! fMergeableCollection ) return 0x0;

  TString identifier = Form("/%s/%s/%s/", physSel.Data(), trigClassName.Data(), centrality.Data());
  
  return fMergeableCollection->GetObject(identifier.Data(), objectName.Data());
//  if ( ! obj ) {
//    CreateMergeableObjects(physSel, trigClassName, centrality);
//    obj = fMergeableCollection->GetObject(identifier.Data(), objectName.Data());
//    AliInfo(Form("Mergeable object collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
//  }
//  return obj;
}

//________________________________________________________________________
TObject* AliMuonAnalysisOutput::GetSum ( TString physSel, TString trigClassNames, TString centrality, TString objectPattern )
{
  //
  /// Sum objects
  /// - physSel, trigClassNames must be in the form: key1,key2
  /// - centrality must be in the form minValue_maxValue
  /// - objectPattern must be in the form match1,match2
  ///   meaning that the object name must contain match1 or match2
  ///   wildcard * is allowed
  
  if ( ! fMergeableCollection ) return 0x0;
  
  TString sumCentralityString = GetCentralities(centrality);

  TObjArray objectNameInCollection;
  objectNameInCollection.SetOwner();
  TObjArray* physSelList = physSel.Tokenize(",");
  TObjArray* trigClassList = trigClassNames.Tokenize(",");
  TObjArray* centralityList = sumCentralityString.Tokenize(",");
  for ( Int_t isel=0; isel<physSelList->GetEntries(); isel++ ) {
    for ( Int_t itrig = 0; itrig<trigClassList->GetEntries(); itrig++ ) {
      for ( Int_t icent=0; icent<centralityList->GetEntries(); icent++ ) {
        TString currId = Form("/%s/%s/%s/", physSelList->At(isel)->GetName(), trigClassList->At(itrig)->GetName(),centralityList->At(icent)->GetName());
        TList* objNameList = fMergeableCollection->CreateListOfObjectNames(currId.Data());
        for ( Int_t iobj=0; iobj<objNameList->GetEntries(); iobj++ ) {
          TString objName = objNameList->At(iobj)->GetName();
          if ( ! objectNameInCollection.FindObject(objName.Data()) ) objectNameInCollection.Add(new TObjString(objName.Data()));
        }
        delete objNameList;
      }
    }
  }
  delete physSelList;
  delete trigClassList;
  delete centralityList;

  TObjArray* objPatternList = objectPattern.Tokenize(",");

  TString matchingObjectNames = "";
  for ( Int_t iobj=0; iobj<objectNameInCollection.GetEntries(); iobj++ ) {
    TString objName = objectNameInCollection.At(iobj)->GetName();
    for ( Int_t ipat=0; ipat<objPatternList->GetEntries(); ipat++ ) {
      TString currPattern = objPatternList->At(ipat)->GetName();
      if ( currPattern.Contains("*") ) {
        if ( ! objName.Contains(TRegexp(currPattern.Data(),kTRUE)) ) continue;
      }
      else if ( objName != currPattern ) continue;

      if ( ! matchingObjectNames.IsNull() ) matchingObjectNames.Append(",");
      matchingObjectNames += objName;
    }
  }
  delete objPatternList;

  TString idPattern = Form("/%s/%s/%s/%s", physSel.Data(), trigClassNames.Data(), sumCentralityString.Data(), matchingObjectNames.Data());
  idPattern.ReplaceAll(" ","");
  
  AliDebug(1,Form("Sum pattern %s", idPattern.Data()));
  
  return fMergeableCollection->GetSum(idPattern.Data());
}

////___________________________________________________________________________
//void AliMuonAnalysisOutput::CreateMergeableObjects(TString physSel, TString trigClassName, TString centrality)
//{
//  //
//  /// Create mergeable objects
//  //
//  TObject* obj = 0x0;
//  TString objectName = "";
//  TString identifier = Form("/%s/%s/%s/", physSel.Data(), trigClassName.Data(), centrality.Data());
//  for ( Int_t iobj=0; iobj<fOutputPrototypeList->GetEntries(); ++iobj ) {
//    objectName = fOutputPrototypeList->At(iobj)->GetName();
//    obj = fOutputPrototypeList->At(iobj)->Clone(objectName.Data());
//    fMergeableCollection->Adopt(identifier, obj);
//  } // loop on histos
//}

//______________________________________________
TString AliMuonAnalysisOutput::GetCentralities ( const TString centralityRange ) const
{
  /// Get list of centrality bins in mergeable collection
  /// matching the chosen centrality range
  /// Centrality range is in the form minVal1_maxVal1,minVal2_maxVal2,etc...

  if ( centralityRange.Contains(",") ) return centralityRange;

  TString matchingCentralities = "";
  if ( ! centralityRange.Contains("_") ) {
    AliError("Centrality range must be in the form: minVal_maxVal");
    return matchingCentralities;
  }

  TObjArray* centralityRangeArr = centralityRange.Tokenize("_");
  Double_t minCentr = static_cast<TObjString*>(centralityRangeArr->At(0))->GetString().Atof();
  Double_t maxCentr = static_cast<TObjString*>(centralityRangeArr->At(1))->GetString().Atof();
  delete centralityRangeArr;


  TList* centrKeyList = fMergeableCollection->CreateListOfKeys(2);
  for ( Int_t ikey=0; ikey<centrKeyList->GetEntries(); ikey++ ) {
    TString centr = centrKeyList->At(ikey)->GetName();
    TObjArray* array = centr.Tokenize("_");
    Double_t currMinCentr = static_cast<TObjString*>(array->At(0))->GetString().Atof();
    Double_t currMaxCentr = static_cast<TObjString*>(array->At(1))->GetString().Atof();
    delete array;
    if ( currMinCentr >= minCentr && currMaxCentr <= maxCentr ) {
      if ( ! matchingCentralities.IsNull() ) matchingCentralities.Append(",");
      matchingCentralities += centr;
    }
  }

  AliDebug(2,Form("Centralities in %s: %s",centralityRange.Data(),matchingCentralities.Data()));

  return matchingCentralities;
}
