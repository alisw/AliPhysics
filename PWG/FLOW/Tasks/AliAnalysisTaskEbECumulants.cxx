/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/************************************** 
* template class for student projects * 
**************************************/ 
  
#include <Riostream.h>
#include <AliAnalysisTaskEbECumulants.h>
#include <AliLog.h>
#include <AliAODEvent.h>
#include <AliAODInputHandler.h>
#include <TFile.h>
#include <TRandom3.h>

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

ClassImp(AliAnalysisTaskEbECumulants)

//================================================================================================================

AliAnalysisTaskEbECumulants::AliAnalysisTaskEbECumulants(const char *name): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 fUseFisherYates(kFALSE),
 fRandomIndices(NULL),
 fVerbose(kFALSE),
 // Event histograms:
 fEventHistogramsList(NULL),
 fEventHistogramsPro(NULL),
 // Final results:
 fFinalResultsList(NULL)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskEbECumulants::AliAnalysisTaskEbECumulants(const char *name)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("outputEbECumulantsAnalysis");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();

  // Define input and output slots here. 
  // Below is the old example, now obsolete.

  // Input slot #0 works with an AliFlowEventSimple
  //DefineInput(0, AliFlowEventSimple::Class());  
  // Input slot #1 is needed for the weights input file:
  //if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());   
  //}  
  // Output slot #0 is reserved              
  // Output slot #1 writes into a TList container

  DefineOutput(1, TList::Class());  

} // AliAnalysisTaskEbECumulants::AliAnalysisTaskEbECumulants(): 

//================================================================================================================

AliAnalysisTaskEbECumulants::AliAnalysisTaskEbECumulants(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 fUseFisherYates(kFALSE),
 fRandomIndices(NULL),
 fVerbose(kFALSE),
 // Event histograms:
 fEventHistogramsList(NULL),
 fEventHistogramsPro(NULL),
 // Final results:
 fFinalResultsList(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskEbECumulants::AliAnalysisTaskEbECumulants()");

} // AliAnalysisTaskEbECumulants::AliAnalysisTaskEbECumulants():

//================================================================================================================

AliAnalysisTaskEbECumulants::~AliAnalysisTaskEbECumulants()
{
 // Destructor.

 if(fHistList) delete fHistList;

 if(fUseFisherYates) { delete fRandomIndices; fRandomIndices = NULL; }
  
} // AliAnalysisTaskEbECumulants::~AliAnalysisTaskEbECumulants()

//================================================================================================================

void AliAnalysisTaskEbECumulants::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // *) Trick to avoid name clashes, part 2.
  
 // a) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // *) Book random generator:
 delete gRandom;
 gRandom = new TRandom3(0); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID 
 //gRandom = new TRandom3(fRandomSeed); // TBI 20211115 add setter for arbitrary seed 

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 this->BookEventHistograms();
 this->BookFinalResultsHistograms();

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskEbECumulants::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskEbECumulants::UserExec(Option_t *) 
{
 // Main loop (called for each event).
 // a) Get pointer to AOD event:
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.

 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}

 // b) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to a track
  if(!aTrack){continue;} // protection against NULL pointers
  if(!aTrack->TestFilterBit(128)){continue;} // filter bit 128 denotes TPC-only tracks

  // example variables for each track: (for more options, please see class STEER/AOD/AliAODTrack.h)  
  Double_t pt = aTrack->Pt(); // Pt
  //Double_t px = aTrack->Px(); // x-component of momenta
  //Double_t py = aTrack->Py(); // y-component of momenta
  //Double_t pz = aTrack->Pz(); // z-component of momenta
  //Double_t e = aTrack->E();  // energy
  //Double_t phi = aTrack->Phi(); // azimuthal angle
  //Double_t eta = aTrack->Eta(); // pseudorapidity
  //Double_t charge = aTrack->Charge(); // charge
 
  // apply some cuts: e.g. take for the analysis only particles in 0.2 < pT < 5.0 GeV
  if( (0.2 < pt) && (pt < 5.0) ) // example cuts
  {
   // fill some control histograms:
   //fPtHist->Fill(pt); // filling pt distribution
 
   // do some analysis only with the particles which passed the cuts
   // ... your analysis code ... 

  } // if( (0.2 < pT) && (pT < 5.0) )

 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
  
 // c) Reset event-by-event objects:
 // ...

 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskEbECumulants::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskEbECumulants::Terminate(Option_t *)
{
 // Accessing the merged output list.

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:

 // ... your code for offline calculations ...


} // end of void AliAnalysisTaskEbECumulants::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskEbECumulants::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

 // *) Event histograms;
 
 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // *) Event histograms:
 for(Int_t t=0;t<gEventHistogramsEbE;t++) // type, see enum eEventHistograms
 {
  for(Int_t rs=0;rs<2;rs++) // reco/sim
  {
   for(Int_t ba=0;ba<2;ba++) // before/after cuts
   {
    fEventHistograms[t][rs][ba] = NULL;
   } // for(Int_t ba=0;ba<2;ba++)
  } // for(Int_t rs=0;rs<2;rs++) // reco/sim
 } // for(Int_t t=0;t<gEventHistogramsEbE;t++) // type, see enum eEventHistograms

} // void AliAnalysisTaskEbECumulants::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskEbECumulants::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskEbECumulants::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 /*
 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("ControlHistograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);
 */

 // b) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fHistList->Add(fFinalResultsList);

} // void AliAnalysisTaskEbECumulants::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskEbECumulants::BookEventHistograms()
{
 // Book all event histograms.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 TString stype[gEventHistogramsEbE] = {"NumberOfEvents","TotalMultiplicity","SelectedParticles","Centrality","Vertex_x","Vertex_y","Vertex_z"}; // keep in sync. with enum eEventHistograms
 TString srs[2] = {"rec","sim"};
 TString sba[2] = {"before cuts","after cuts"};

 for(Int_t t=0;t<gEventHistogramsEbE;t++) // type, see enum eEventHistograms
 {
  if(!fBookEventHistograms[t]){continue;}
  for(Int_t rs=0;rs<2;rs++) // reco/sim
  {
   for(Int_t ba=0;ba<2;ba++) // before/after cuts
   {
    // Skip exceptional cases:
    if(eSelectedParticles == t && eBefore == ba){continue;} // Number of selected particles makes sense only after cuts
    // ...
    // Book the rest:
    fEventHistograms[t][rs][ba] = new TH1D(Form("fEventHistograms[%d][%d][%d]",t,rs,ba),Form("%s, %s, %s",stype[t].Data(),srs[rs].Data(),sba[ba].Data()),(Int_t)fEventHistogramsBins[t][0],fEventHistogramsBins[t][1],fEventHistogramsBins[t][2]); 
    //fEventHistograms[t][rs][ba]->SetLineColor(fBeforeAfterColor[ba]);  
    //fEventHistograms[t][rs][ba]->SetFillColor(fBeforeAfterColor[ba]-10);
    fEventHistogramsList->Add(fEventHistograms[t][rs][ba]);
   } // for(Int_t ba=0;ba<2;ba++)
  } // for(Int_t rs=0;rs<2;rs++) // reco/sim
 } // for(Int_t t=0;t<gEventHistogramsEbE;t++) // type, see enum 'eEvent'

} // void AliAnalysisTaskEbECumulants::BookEventHistograms()

//=======================================================================================================================

void AliAnalysisTaskEbECumulants::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

} // void AliAnalysisTaskEbECumulants::BookFinalResultsHistograms()

//=======================================================================================================================

void AliAnalysisTaskEbECumulants::Red(const char* text)
{ 
 cout<<"\n\033[1;31m"<<text<<"\033[0m\n"<<endl;
} 

//=======================================================================================

void AliAnalysisTaskEbECumulants::Green(const char* text)
{ 
 cout<<"\n\033[1;32m"<<text<<"\033[0m\n"<<endl;
}
//=======================================================================================

void AliAnalysisTaskEbECumulants::Yellow(const char* text)
{ 
 cout<<"\n\033[1;33m"<<text<<"\033[0m\n"<<endl;
} 

//=======================================================================================

void AliAnalysisTaskEbECumulants::Blue(const char* text)
{ 
 cout<<"\n\033[1;34m"<<text<<"\033[0m\n"<<endl;
} 

//=======================================================================================

TObject* AliAnalysisTaskEbECumulants::GetObjectFromList(TList *list, Char_t *objectName)
{
 // Get TObject pointer from TList, even if it's in some nested TList. Foreseen to be used to fetch histograms or profiles from files directly. 
 // Some ideas taken from TCollection::ls() 
 // If you have added histograms directly to files (without TList's), then you can fetch them directly with file->Get("hist-name"). 
 
 // Usage: TH1D *hist = (TH1D*) GetObjectFromList("some-valid-TList-pointer","some-object-name");

 // Example: GetObjectFromList("some-valid-TList-pointer","some-object-name")->Draw(); // yes, for histograms and profiles this is just fine

 // Last update: 20210911

 // To do: 
 // a) If I have objects with same name, nested in different TLists, what then?
 
 // Insanity checks:  
 if(!list){cout<<__LINE__<<endl;exit(1);}
 if(!objectName){cout<<__LINE__<<endl;exit(1);}
 if(0 == list->GetEntries()){return NULL;}

 // The object is in the current base list:
 TObject *objectFinal = list->FindObject(objectName); // final object I am after
 if(objectFinal) return objectFinal;

 // Search for object recursively in the nested lists:
 TObject *objectIter; // iterator object in the loop below
 TIter next(list);
 while((objectIter = next())) // double round braces are to silent the warnings
 {
  if(TString(objectIter->ClassName()).EqualTo("TList"))
  {
   objectFinal = GetObjectFromList((TList*)objectIter,objectName);
   if(objectFinal) return objectFinal;
  }
 } // while(objectIter = next()) 

 return NULL;

} // TObject* AliAnalysisTaskEbECumulants::GetObjectFromList(TList *list, Char_t *objectName)

//=======================================================================================

Int_t AliAnalysisTaskEbECumulants::NumberOfNonEmptyLines(const char *externalFile)
{
 // Count number of non-empty lines in some external file.

 if(gSystem->AccessPathName(externalFile,kFileExists))
 {
  Red(Form("if(gSystem->AccessPathName(externalFile,kFileExists)), externalFile = %s",externalFile)); 
  cout<<__LINE__<<endl;
  exit(1);
 }

 string line;
 ifstream myfile;
 myfile.open(externalFile);
 Int_t nLines = 0;
 while (getline(myfile,line))
 { 
  if(TString(line).EqualTo("")){continue;}
  nLines++;
 }
 myfile.close();

 return nLines;

} // Int_t AliAnalysisTaskEbECumulants::NumberOfNonEmptyLines(const char *externalFile)

//=======================================================================================

void AliAnalysisTaskEbECumulants::ResetEventByEventQuantities()
{
 // Reset all global event-by-event quantities here:

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // *) Fisher-Yates algorithm:
 if(fUseFisherYates)
 {
  delete fRandomIndices; fRandomIndices = NULL; 
 }

} // void AliAnalysisTaskEbECumulants::ResetEventByEventQuantities()

//=======================================================================================

void AliAnalysisTaskEbECumulants::RandomIndices(AliVEvent *ave)
{
 // Randomize indices using Fisher-Yates algorithm. 

 // a) Determine Ali{MC,ESD,AOD}Event;
 // b) Get total number of tracks;
 // c) Fisher-Yates algorithm.

 if(fVerbose){Green(__PRETTY_FUNCTION__);}

 // a) Determine Ali{MC,ESD,AOD}Event:
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 //AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 // b) Get total number of tracks:
 Int_t nTracks = 0;
 if(aAOD)
 {
  nTracks = aAOD->GetNumberOfTracks();
 }
 else if(aMC)
 {
  nTracks = aMC->GetNumberOfTracks();
 }

 if(nTracks<1){return;}

 // c) Fisher-Yates algorithm:
 fRandomIndices = new TArrayI(nTracks);
 fRandomIndices->Reset(); 
 for(Int_t i=0;i<nTracks;i++)
 {
  fRandomIndices->AddAt(i,i);
 }
 for(Int_t i=nTracks-1;i>=1;i--)
 {
  Int_t j = gRandom->Integer(i+1);
  Int_t temp = fRandomIndices->GetAt(j);
  fRandomIndices->AddAt(fRandomIndices->GetAt(i),j);
  fRandomIndices->AddAt(temp,i);
 } // end of for(Int_t i=nTracks-1;i>=1;i--) 

} // void AliAnalysisTaskEbECumulants::RandomIndices(AliVEvent *ave)

//=======================================================================================



