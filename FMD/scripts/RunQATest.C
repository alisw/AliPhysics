#ifdef _CINT__
class AliFMDInput;
class AliQADataMaker;
#endif
#ifdef BUILD 
# ifndef _CINT__
#  include <AliCDBManager.h>
#  include <AliQADataMaker.h>
#  include "AliFMDQADataMakerRec.h"
#  include "AliFMDQADataMakerSim.h"
#  include <AliQAv1.h>
#  include <TObjArray.h>
#  include <AliRecoParam.h>
#  include <AliRawReader.h>
#  include <AliQACheckerBase.h>
#  include <AliFMDQAChecker.h>
#  include <AliQAChecker.h>
#  include <TCanvas.h>
#  include <TMath.h>
#  include <AliFMDInput.h>
# endif

/** 
 * Class to test the QA code 
 * 
 * 
 */
class QATest : public AliFMDInput
{
public:
  /** 
   * Constructor 
   */
  QATest()
    : AliFMDInput("galice.root"),
      fMaker(0), 
      fArray(0),
      fSpecie(AliRecoParam::kLowMult), 
      fCycleLength(-1)
  {
    for (Int_t i = 0; i < AliQAv1::kNTASKINDEX; i++)  
      fTasks[i] = AliQAv1::kNULLTASKINDEX;
    Int_t nArray = AliQAv1::kNTASKINDEX * AliRecoParam::kNSpecies;
    Info("QAtest", "Allocating %dx%d=%d TObjArrays",
	 AliQAv1::kNTASKINDEX, AliRecoParam::kNSpecies, nArray);
    fArray = new TObjArray*[nArray];
    for (Int_t i = 0; i < nArray; i++) fArray[i] = 0;
  }
  /** 
   * Set the location of the QA reference storage 
   * 
   * @param refloc Location of QA storage (e.g., local://${ALICE_ROOT}/QAref)
   */
  void SetQARefStorage(const char* refloc)
  {
    AliQAv1::SetQARefStorage(refloc);
  }
  /** 
   * Calculate index of TObjArray* in fArray
   * 
   * @param specie Event species
   * @param task   Task index 
   * 
   * @return Index. 
   */
  Int_t CalcArrayIndex(AliQAv1::TASKINDEX_t task,
		       AliRecoParam::EventSpecie_t specie) const
  {
    Int_t es = AliRecoParam::AConvert(specie);
    if (es >= AliRecoParam::kNSpecies) return -1;
    
    Int_t base = CalcSpeciesIndex(task);
    if (base < 0) return base;
    return base + es;
  }
  /** 
   * Get index to species array 
   * 
   * @param task Task 
   * 
   * @return 
   */
  Int_t CalcSpeciesIndex(AliQAv1::TASKINDEX_t task) const
  {
    if (task >= AliQAv1::kNTASKINDEX) return -1;
    // Species are consequtive 
    return task * AliRecoParam::kNSpecies;
  }
    
  /** 
   * Initialise the code
   * 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Init()
  {
    // --- Create the maker ------------------------------------------
    if (IsLoaded(kHits)       || 
	IsLoaded(kDigits)     || 
	IsLoaded(kKinematics) ||
	IsLoaded(kSDigits)) {
      AddLoad(kHeader);
      fMaker = new AliFMDQADataMakerSim();
    }
    else 
      fMaker = new AliFMDQADataMakerRec();

    // --- Figure out tasks ------------------------------------------
    Int_t j = 0;
    if (IsLoaded(kHits))      fTasks[j++]	  = AliQAv1::kHITS;
    if (IsLoaded(kDigits))    fTasks[j++]	  = AliQAv1::kDIGITS;
    if (IsLoaded(kSDigits))   fTasks[j++]	  = AliQAv1::kSDIGITS;
    if (IsLoaded(kRecPoints)) fTasks[j++]         = AliQAv1::kRECPOINTS;
    if (IsLoaded(kESD))       fTasks[j++]         = AliQAv1::kESDS;
    if (IsLoaded(kRaw))       fTasks[j++]         = AliQAv1::kRAWS;
    if (j == 0) { 
      AliError(Form("Loaded trees (%s) cannot be processed by QA", 
		    LoadedString(true)));
      return false;
    }
    // Int_t nTasks = j;

    // --- Data maker ------------------------------------------------
    Info("TestQA", "Creating data maker");
    fMaker = new AliFMDQADataMakerRec();
    
    // --- Init all species histograms -------------------------------
    Info("TestQA", "Setup data species");
    AliQAv1* qa = AliQAv1::Instance();
    for (unsigned int es = 0; es < AliRecoParam::kNSpecies; es++) {
      AliRecoParam::EventSpecie_t specie = AliRecoParam::ConvertIndex(es);
      fMaker->SetEventSpecie(specie);
      qa->SetEventSpecie(specie);
      for (Int_t i = 0; i < AliQAv1::kNTASKINDEX; i++) {
	if (fTasks[i] == AliQAv1::kNULLTASKINDEX) continue;
	Int_t k = CalcArrayIndex(fTasks[i], specie);
	Info("Init", "Array for task %d (%s), specie %d (%s) @ %d/%d", 
	     fTasks[i], AliQAv1::GetTaskName(fTasks[i]).Data(), 
	     specie, AliRecoParam::GetEventSpecieName(specie), k,
	     AliQAv1::kNTASKINDEX * AliRecoParam::kNSpecies);
	fArray[k] = fMaker->Init(fTasks[i], specie);
	fArray[k]->ls();
      }
    }

    // --- Start of cycle --------------------------------------------
    Int_t  run  = AliCDBManager::Instance()->GetRun();
    Bool_t same = false;
    for (Int_t i = 0; i < AliQAv1::kNTASKINDEX; i++) {
      if (fTasks[i] == AliQAv1::kNULLTASKINDEX) continue;
      fMaker->StartOfCycle(fTasks[i], run, same);
      same = !same;
    }
    
    
    return AliFMDInput::Init();
  }
  /** 
   * Process the hits 
   * 
   * @return true on success 
   */
  virtual Bool_t ProcessHits()
  {
    fMaker->SetEventSpecie(fSpecie);
    fMaker->MakeHits(fTreeH);
    return true;
  }
  /** 
   * Process the digits
   * 
   * @return true on success 
   */
  virtual Bool_t ProcessDigits()
  {
    fMaker->SetEventSpecie(fSpecie);
    fMaker->MakeDigits(fTreeD);
    return true;
  }
  /** 
   * Process the summable digits
   * 
   * @return true on success 
   */
  virtual Bool_t ProcessSDigits()
  {
    fMaker->SetEventSpecie(fSpecie);
    fMaker->MakeSDigits(fTreeS);
    return true;
  }
  /** 
   * Process the reconstructed points 
   * 
   * @return true on success 
   */
  virtual Bool_t ProcessRecPoints()
  {
    fMaker->SetEventSpecie(fSpecie);
    fMaker->MakeRecPoints(fTreeR);
    return true;
  }
  /** 
   * Process the event summary data
   * 
   * @return true on success 
   */
  virtual Bool_t ProcesssESDs()
  {
    fMaker->SetEventSpecie(fSpecie);
    fMaker->MakeESDs(fESDEvent);
    return true;
  }
  /** 
   * Process the raw data 
   * 
   * @return true on success 
   */
  virtual Bool_t ProcessRawDigits()
  {
    fMaker->SetEventSpecie(fSpecie);
    fMaker->MakeRaws(fReader);
    return true;
  }
  virtual Bool_t End()
  {
    Bool_t ret = AliFMDInput::End();
    if (fCycleLength < 0) return ret;

    if (!(fEventCount != 0 && (fEventCount % fCycleLength) == 0))
      return ret;

    // --- End of cycle - this calls the ecker ---------------------
    Info("End", "End of cycle");
    for (Int_t i = 0; i < AliQAv1::kNTASKINDEX; i++) {
      if (fTasks[i] == AliQAv1::kNULLTASKINDEX) continue;
      fMaker->EndOfCycle(fTasks[i]);
    }

    // --- Get the checker -------------------------------------------
    Info("End", "Running checker");
    AliQACheckerBase * checker = AliQAChecker::Instance()->
      GetDetQAChecker(AliQAv1::GetDetIndex("FMD"));
    ((AliFMDQAChecker*)checker)->SetDoScale();
    
    // --- Test: Remake plots ----------------------------------------
    for (unsigned int idx = 0; idx < AliQAv1::kNTASKINDEX; idx++) {
      // AliRecoParam::EventSpecie_t specie = AliRecoParam::ConvertIndex(es);
      // if (!qa->IsEventSpecieSet(specie)) continue;
      if (fTasks[idx] == AliQAv1::kNTASKINDEX) continue;
      AliQAv1::TASKINDEX_t task = AliQAv1::TASKINDEX_t(idx); // AliQAv1::kRAWS;
      AliQAv1::MODE_t      mode = AliQAv1::kRECMODE;
      Int_t k = CalcSpeciesIndex(task);
      Info("End", "Array for task %d (%s) @ %d: %p", 
	   task, AliQAv1::GetTaskName(task).Data(), k, fArray[k]);
      if (!fArray[k]) continue;
      fArray[k]->ls();
      checker->MakeImage(&(fArray[k]), task, mode);
    }
    return ret;
  }
  /** 
   * Called at the end of the job.  Runs the checkers
   * 
   * @return true on success 
   */
  virtual Bool_t Finish()
  {
    // --- Finish maker ----------------------------------------------
    fMaker->Finish();
  
    // --- Get the checker ------------------------------------------- 
    AliQACheckerBase * checker = AliQAChecker::Instance()->
      GetDetQAChecker(AliQAv1::GetDetIndex("FMD")); 

    // --- Get images from checker -----------------------------------
    AliQAv1* qa = AliQAv1::Instance();
    TObjArray* canvases = new TObjArray();
    for (unsigned int es = 0; es < AliRecoParam::kNSpecies; es++) {
      AliRecoParam::EventSpecie_t specie = AliRecoParam::ConvertIndex(es);
      if (!qa->IsEventSpecieSet(specie)) continue;
    
      TCanvas* c = checker->GetImage(specie);
      if (!c) continue; 

      canvases->Add(c);
    }
  
    // --- Create summary image --------------------------------------
    TCanvas* out = new TCanvas("summary", "Summary", 1024, 1024);
    out->SetFillColor(kWhite);
    out->SetBorderSize(0);
    out->SetBorderMode(0);
  
    Int_t nImgs = canvases->GetEntriesFast();
    Int_t nX    = Int_t(TMath::Sqrt(nImgs) + .5);
    Int_t nY    = nX;				
    out->Divide(nX, nY);
    for (Int_t i = 0; i < nImgs; i++) { 
      TVirtualPad* p = out->cd(i + 1);
      p->SetRightMargin(0.001);
      p->SetTopMargin(0.001);
      if (!p) { 
	Warning("TestQA", "No pad at index %d / %d", i+1, nImgs);
	continue;
      }

      TCanvas* c = static_cast<TCanvas*>(canvases->At(i));
      c->DrawClonePad();
    }
    out->Print("summary.png");    

    return true;
  }
  void SetSpecie(AliRecoParam::EventSpecie_t s) { fSpecie = s; }
protected:
  AliQADataMaker*             fMaker; // Data maker 
  AliQAv1::TASKINDEX_t        fTasks[AliQAv1::kNTASKINDEX]; // Tasks to do
  TObjArray**                 fArray;
  AliRecoParam::EventSpecie_t fSpecie;
  Int_t                       fCycleLength;
  ClassDef(QATest, 0); 
};

#else
void
RunQATest(const char* src, const char* specie="low", Int_t runno=0)
{
  gROOT->LoadMacro("$ALICE_ROOT/FMD/scripts/Compile.C");
  gSystem->AddIncludePath("-DBUILD=1");
  Compile("$ALICE_ROOT/../trunk/FMD/scripts/RunQATest.C", "+g");

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(runno);
  
  QATest*                     qaTest = new QATest;
  TString                     what(src);
  TString                     spec(specie);
  Int_t                       colon = what.Index(":");
  TString                     type = "";
  if (colon != TString::kNPOS) { 
    type = what(0, colon);
    what = what(colon+1, what.Length()-colon-1);
  }
  else { 
    if      (what.Contains("AliESD")) type = "esd";
    else if (what.Contains("galice")) type = "sim";
    else if (what.Contains(".raw"))   type = "raw";
    else if (what.Contains(".root"))  type = "raw";
  }
  Info("RunQATest", "type=%s, what=%s specie=%s", 
       type.Data(), what.Data(), spec.Data());

  spec.ToLower();
  if      (spec.Contains("low"))    qaTest->SetSpecie(AliRecoParam::kLowMult);
  else if (spec.Contains("high"))   qaTest->SetSpecie(AliRecoParam::kHighMult);
  else if (spec.Contains("cosmic")) qaTest->SetSpecie(AliRecoParam::kCosmic);
  else if (spec.Contains("calib"))  qaTest->SetSpecie(AliRecoParam::kCalib);

  type.ToLower();
  if (type.CompareTo("esd") == 0) { 
    qaTest->AddLoad(AliFMDInput::kESD);
    qaTest->SetInputDir(what);
  }
  else if (type.CompareTo("raw") == 0) {
    qaTest->AddLoad(AliFMDInput::kRaw);
    qaTest->SetRawFile(what);
  }
  else if (type.CompareTo("sim") == 0) { 
    qaTest->AddLoad(AliFMDInput::kHits);
    qaTest->AddLoad(AliFMDInput::kSDigits);
    qaTest->AddLoad(AliFMDInput::kDigits);
  }
  else if (type.CompareTo("rec") == 0) { 
    qaTest->AddLoad(AliFMDInput::kDigits);
    qaTest->AddLoad(AliFMDInput::kRecPoints);
    qaTest->AddLoad(AliFMDInput::kESD);
  }
  else { 
    Error("RunQATest", "Unknown type='%s' in '%s'", type.Data(), src);
    return;
  }
  qaTest->Run(1000);
}
#endif
  
  

//
// EOF
//
