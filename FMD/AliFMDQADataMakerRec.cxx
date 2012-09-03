/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
// --- ROOT system ---
#include <iostream>
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2I.h> 
#include <TGeoManager.h>

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliFMDQADataMakerRec.h"
#include "AliFMDDigit.h"
#include "AliFMDRecPoint.h"
#include "AliQAChecker.h"
#include "AliESDFMD.h"
#include "AliFMDParameters.h"
#include "AliFMDRawReader.h"
#include "AliFMDReconstructor.h"
#include "AliRawReader.h"
#include "AliFMDAltroMapping.h"
#include "AliFMDDebug.h"

namespace {
  Int_t colors[3] = {kRed,kGreen,kBlue};
}
//_____________________________________________________________________
// This is the class that collects the QA data for the FMD during
// reconstruction.  
//
// The following data types are picked up:
// - rec points
// - esd data
// - raws
// Author : Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
//_____________________________________________________________________

ClassImp(AliFMDQADataMakerRec)
#if 0
; // For Emacs - do not delete!
#endif
           
//_____________________________________________________________________
AliFMDQADataMakerRec::AliFMDQADataMakerRec() 
  : AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kFMD), 
		      "FMD Quality Assurance Data Maker"),
    fRecPointsArray("AliFMDRecPoint", 1000), 
    fReconstructor(0),
    fUseReconstructor(true)
{
  // ctor
 
}

//_____________________________________________________________________
AliFMDQADataMakerRec::AliFMDQADataMakerRec(const AliFMDQADataMakerRec& qadm) 
  : AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kFMD), 
		      "FMD Quality Assurance Data Maker"),
    fRecPointsArray(qadm.fRecPointsArray), 
    fReconstructor(qadm.fReconstructor),
    fUseReconstructor(qadm.fUseReconstructor)
{
  // copy ctor 
  // Parameters: 
  //    qadm    Object to copy from
  
}
//_____________________________________________________________________
AliFMDQADataMakerRec& 
AliFMDQADataMakerRec::operator = (const AliFMDQADataMakerRec& qadm ) 
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    qadm What to assign from 
  // 
  // Return:
  //    Reference to this
  //
  fRecPointsArray   = qadm.fRecPointsArray;
  fReconstructor    = qadm.fReconstructor;
  fUseReconstructor = qadm.fUseReconstructor;
  return *this;
}
//_____________________________________________________________________
AliFMDQADataMakerRec::~AliFMDQADataMakerRec()
{
  // 
  // Destrcutor 
  // 
}


//_____________________________________________________________________ 

void 
AliFMDQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, 
					 TObjArray ** list)
{
  // Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses(); // reset triggers list to select all histos
  AliLog::Message(5,"FMD: end of detector cycle",
		  "AliFMDQADataMakerRec","AliFMDQADataMakerRec",
		  "AliFMDQADataMakerRec::EndOfDetectorCycle",
		  "AliFMDQADataMakerRec.cxx",95);
  AliQAChecker::Instance()->Run(AliQAv1::kFMD, task, list);
}

//_____________________________________________________________________ 
TH1* AliFMDQADataMakerRec::MakeADCHist(UShort_t d, Char_t r, Short_t b)
{
  TString name("adc"); 
  TString title("ADC counts");
  Int_t   color = kRed+1; 
  if (d > 0) { 
    name.Append(Form("FMD%d%c", d, r));
    title.Append(Form(" in FMD%d%c", d, r));
    color = colors[d-1]+3+(r == 'I' || r == 'i' ? 0 : 1);
    if (b >= 0) {
      name.Append(Form("_0x%02x", b));
      title.Append(Form("[0x%02x]", b));
    }
  }
  TH1* hist = new TH1F(name, title,1024,0,1024);
  hist->SetXTitle("Amplitude [ADC counts]");
  hist->SetYTitle("Events [log]");
  hist->SetFillStyle(3001);
  hist->SetFillColor(color);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->GetXaxis()->SetNdivisions(408,false);
  hist->SetDirectory(0);
  // hist->SetStats(0);

  return hist;
}
//_____________________________________________________________________ 
TH1* AliFMDQADataMakerRec::MakeELossHist(UShort_t d, Char_t r, Short_t b)
{
  TString name("eloss"); 
  TString title("Energy loss");
  Int_t   color = kBlue+1;
  if (d > 0) { 
    name.Append(Form("FMD%d%c", d, r));
    title.Append(Form(" in FMD%d%c", d, r));
    color = colors[d-1]+3+(r == 'I' || r == 'i' ? 0 : 1);
    if (b >= 0) {
      name.Append(Form("_0x%02x", b));
      title.Append(Form("[0x%02x]", b));
    }
  }
  TH1* hist = new TH1F(name, title,600,0, 15);
  hist->SetXTitle("#Delta/#Delta_{mip}");
  hist->SetYTitle("Events [log]");
  hist->SetFillStyle(3001);
  hist->SetFillColor(color);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->Sumw2();
  hist->SetDirectory(0);
  // hist->SetStats(0);

  return hist;
}


//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitESDs()
{
  // create Digits histograms in Digits subdir
  Info("InitESDs", "Initializing ESDs");
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1* hist = MakeELossHist();
  Add2ESDsList(hist, 0, !expert, image);
  ClonePerTrigClass(AliQAv1::kESDS); // this should be the last line    
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  Info("InitDigits", "Initializing Digits");
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1* hist = MakeADCHist();
  Add2DigitsList(hist, 0, !expert, image);
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  Info("InitRecPoints", "Initializing RecPoints");
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1* hist = MakeELossHist();
  Add2RecPointsList(hist,0, !expert, image);
  ClonePerTrigClass(AliQAv1::kRECPOINTS); // this should be the last linea
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::InitRaws()
{
  // create Raws histograms in Raws subdir  
  Info("InitRaws", "Initializing Raws");
  const Bool_t expert   = kTRUE ; 
  const Bool_t saveCorr = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH2I* hErrors = new TH2I("readoutErrors", "Read out errors", 3, .5, 3.5,
			   160, -.5, 159.5); 
  hErrors->GetXaxis()->SetBinLabel(1, "FMD1");
  hErrors->GetXaxis()->SetBinLabel(2, "FMD2");
  hErrors->GetXaxis()->SetBinLabel(3, "FMD3");
  hErrors->SetYTitle("# errors [log]");
  hErrors->SetZTitle("Events [log]");
  hErrors->SetDirectory(0);
  Add2RawsList(hErrors, 1, !expert, image, !saveCorr);
  //AliInfo(Form("Adding %30s to raw list @ %2d", hErrors->GetName(), 1));

  if (fUseReconstructor && !fReconstructor) {
    // Int_t oldDbg = AliLog::GetDebugLevel("FMD","");
    // AliLog::SetModuleDebugLevel("FMD", 5);

    if (!gGeoManager) {
      Info("InitRaws", "Loading the geometry");
      AliGeomManager::LoadGeometry();
    }

    fReconstructor = new AliFMDReconstructor();
    fReconstructor->SetDiagnose(false);
    fReconstructor->Init();
    // AliLog::SetModuleDebugLevel("FMD", oldDbg);
  }

  TH2* status = new TH2D("status", "Status per cycle", 
			  5, .5, 5.5, 4, -.5, 3.5);
  status->SetDirectory(0);
  status->SetXTitle("Detector");
  status->SetYTitle("Status");
  status->SetZTitle("N_{cycles} [LOG]");
  status->GetXaxis()->SetBinLabel(1, "FMD1i");
  status->GetXaxis()->SetBinLabel(2, "FMD2i");
  status->GetXaxis()->SetBinLabel(3, "FMD2o");
  status->GetXaxis()->SetBinLabel(4, "FMD3i");
  status->GetXaxis()->SetBinLabel(5, "FMD3o");
  status->GetYaxis()->SetBinLabel(1, "OK");
  status->GetYaxis()->SetBinLabel(2, "Problem");
  status->GetYaxis()->SetBinLabel(3, "Bad");
  status->GetYaxis()->SetBinLabel(4, "What the ...?");
  status->SetStats(0);
  Add2RawsList(status, GetHalfringIndex(4, 'i', 0, 0), 
	       !expert, image, !saveCorr);
	       
  TH1* hist;
  Int_t idx = 0;
  for(UShort_t d = 1; d<=3; d++) {
    UShort_t nR = (d == 1 ? 1 : 2); 
    for(UShort_t q = 0; q < nR; q++) {
      Char_t r = (q == 1 ? 'O' : 'I');
      hist = (fUseReconstructor ? 
	      MakeELossHist(d, r, -1) : 
	      MakeADCHist(d, r, -1));

      Int_t index1 = GetHalfringIndex(d, r, 0, 1);
      idx          = TMath::Max(index1, idx);
      Add2RawsList(hist, index1, !expert, image, !saveCorr);
      //AliInfo(Form("Adding %30s to raw list @ %2d", hist->GetName(), index1));
      
      // If we're using the reconstructor, do not make expert histograms 
      if (fUseReconstructor) continue;
      
      for(UShort_t b = 0; b <= 1; b++) {
	//Hexadecimal board numbers 0x0, 0x1, 0x10, 0x11;
	UShort_t board = (q == 1 ? 0 : 1) + b*16;

	hist = MakeADCHist(d, r, board);
	Int_t index2 = GetHalfringIndex(d, r, board/16,0);
	idx          = TMath::Max(index2, idx);
	Add2RawsList(hist, index2, expert, !image, !saveCorr);
	//AliInfo(Form("Adding %30s to raw list @ %2d",hist->GetName(),index2));
      }
    }
  }
  //
  ClonePerTrigClass(AliQAv1::kRAWS); // this should be the last line
}

#if 0
struct FillESDHist : public AliESDFMD::ForOne
{
  FillESDHist(AliFMDQADataMakerRec* m) : fM(m) {}
  FillESDHist(const FillESDHist& o) : fM(o.fM) {}
  FillESDHist& operator=(const FillESDHist& o) { fM = o.fM; return *this;  }
  Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, Float_t m, Float_t) 
  {
    // Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
    if(m == AliESDFMD::kInvalidMult) return true;
    
    fM->FillESDsData(0,m);
    return true;
  }
  AliFMDQADataMakerRec* fM;
};
#endif

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // 
  // Analyse ESD event
  // 
  // Parameters:
  //    esd ESD event
  //
  if(!esd) {
    AliError("FMD ESD object not found!!") ; 
    return;
  }
  AliFMDDebug(2, ("Will loop over ESD data and fill histogram"));

  AliESDFMD* fmd = esd->GetFMDData();
  if (!fmd) return;

#if 0
  FillESDHist f(this);
  fmd->ForEach(f);
#else



  // FIXME - we should use AliESDFMD::ForOne subclass to do this!
  for(UShort_t det=1;det<=3;det++) {
    UShort_t nrng = (det == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nrng; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult) continue;
	  
	  FillESDsData(0,mult);
	}
      }
    }
  }
#endif
  IncEvCountCycleESDs();
  IncEvCountTotalESDs();
}


//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeDigits()
{
  // makes data from Digits  
  if(!fDigitsArray)  {
    AliError("FMD Digit object not found!!") ;
    return;
  }
  
  for(Int_t i=0;i<fDigitsArray->GetEntriesFast();i++) {
    //Raw ADC counts
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(fDigitsArray->At(i));
    FillDigitsData(0,digit->Counts());
  }
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeDigits(TTree * digitTree)
{
  // 
  // Analyse digits
  // 
  // Parameters:
  //    digitTree Tree of digits
  //
  if (fDigitsArray) 
    fDigitsArray->Clear();
  else 
    fDigitsArray = new TClonesArray("AliFMDDigit", 1000);

  TBranch*      branch = digitTree->GetBranch("FMD");
  if (!branch) {
    AliWarning("FMD branch in Digit Tree not found") ; 
    return;
  } 
  branch->SetAddress(&fDigitsArray);
  branch->GetEntry(0); 
  MakeDigits();
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeRaws(AliRawReader* rawReader)
{
  // 
  // Analyse raw 
  // 
  // Parameters:
  //    rawReader Raw reader
  //
  AliFMDRawReader fmdReader(rawReader,0);
  if (fDigitsArray) fDigitsArray->Clear();
  else              fDigitsArray = new TClonesArray("AliFMDDigit", 1000);
    
  TClonesArray* digitsAddress = fDigitsArray;
    
  rawReader->Reset();
    
  digitsAddress->Clear();
  fmdReader.ReadAdcs(digitsAddress);
  //
  FillRawsData(1,1, fmdReader.GetNErrors(0));
  FillRawsData(1,2, fmdReader.GetNErrors(1));
  FillRawsData(1,3, fmdReader.GetNErrors(2));

  if (fUseReconstructor) { 
    AliESDFMD* fmd = fReconstructor->GetESDObject();
    fmd->Clear();
    
    // AliLog::SetModuleDebugLevel("FMD", 15);
    fReconstructor->ProcessDigits(digitsAddress, fmdReader);

    if (!fmd) AliFatal("No ESD object from reconstructor");

    for(UShort_t det=1;det<=3;det++) {
      UShort_t nrng = (det == 1 ? 1 : 2);
      for (UShort_t ir = 0; ir < nrng; ir++) {
	Char_t   ring = (ir == 0 ? 'I' : 'O');
	UShort_t nsec = (ir == 0 ? 20  : 40);
	UShort_t nstr = (ir == 0 ? 512 : 256);
	for(UShort_t sec =0; sec < nsec;  sec++)  {
	  for(UShort_t strip = 0; strip < nstr; strip++) {
	    Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	    if(mult == AliESDFMD::kInvalidMult) continue;
	    
	    Int_t index1 = GetHalfringIndex(det, ring, 0, 1);
	    FillRawsData(index1,mult);
	  }
	}
      }
    }
  }
  else {
    for(Int_t i=0;i<digitsAddress->GetEntriesFast();i++) {
      //Raw ADC counts
      AliFMDDigit*      digit = static_cast<AliFMDDigit*>(digitsAddress->At(i));
      UShort_t          det   = digit->Detector();
      Char_t            ring  = digit->Ring();
      UShort_t          sec   = digit->Sector();
      // UShort_t strip = digit->Strip();
      AliFMDParameters* pars  = AliFMDParameters::Instance();
      Short_t           board = pars->GetAltroMap()->Sector2Board(ring, sec);
      
      Int_t index1 = GetHalfringIndex(det, ring, 0, 1);
      FillRawsData(index1,digit->Counts());
      Int_t index2 = GetHalfringIndex(det, ring, board/16,0);
      FillRawsData(index2,digit->Counts());
    }
  }
  //
  IncEvCountCycleRaws();
  IncEvCountTotalRaws();
}

//_____________________________________________________________________
void AliFMDQADataMakerRec::MakeRecPoints(TTree* clustersTree)
{
  // makes data from RecPoints
  
   AliFMDParameters* pars = AliFMDParameters::Instance();
  fRecPointsArray.Clear();
  TBranch *fmdbranch = clustersTree->GetBranch("FMD");
  if (!fmdbranch) { 
    AliError("can't get the branch with the FMD recpoints !");
    return;
  }
  
  TClonesArray* RecPointsAddress = &fRecPointsArray;
  
  fmdbranch->SetAddress(&RecPointsAddress);
  fmdbranch->GetEntry(0);
  TIter next(RecPointsAddress) ; 
  AliFMDRecPoint * rp ; 
  while ((rp = static_cast<AliFMDRecPoint*>(next()))) {
    FillRecPointsData(0,rp->Edep()/pars->GetEdepMip());
  }
  IncEvCountCycleRecPoints();
  IncEvCountTotalRecPoints();
  //
}

//_____________________________________________________________________ 
void AliFMDQADataMakerRec::StartOfDetectorCycle()
{
  // Do an init on the reconstructor.  If we have the
  // same run nothing happens, but if we have a new run, we update our
  // parameters.
  if (fUseReconstructor && fReconstructor) fReconstructor->Init();
  if (fRawsQAList) {
    for (Int_t index = 0 ; index < AliRecoParam::kNSpecies ; index++) {
      if (!fRawsQAList[index]) continue;
      AliRecoParam::EventSpecie_t specie = AliRecoParam::ConvertIndex(index);
      if (specie == AliRecoParam::kCalib || specie == AliRecoParam::kCosmic) 
	continue;

      TIter    nextObject(fRawsQAList[index]);
      TObject* object = 0;
      while ((object = nextObject())) { 
	if (!object->InheritsFrom(TH1::Class())) continue;
	TH1* hist = static_cast<TH1*>(object);
	if (!hist->TestBit(BIT(23))) continue;
	
	AliInfoF("Resetting histogram %s", hist->GetName());
	hist->Reset("M");
	hist->SetBit(BIT(23), false);
      }
    }
  }
}
//_____________________________________________________________________ 
Int_t AliFMDQADataMakerRec::GetHalfringIndex(UShort_t det, 
					     Char_t ring, 
					     UShort_t board, 
					     UShort_t monitor)
{
  // 
  // Get the half-ring index
  // 
  // Parameters:
  //    det      Detector
  //    ring     Ring
  //    board    Board number
  //    monitor  Monitor 
  // 
  // Return:
  //    Half ring index
  //  
  UShort_t iring = (ring == 'I' || ring == 'i' ? 1 : 0);
  Int_t    index = ((((det-1) & 0x3) << 3) | 
		    ((iring   & 0x1) << 2) |
		    ((board   & 0x1) << 1) | 
		    ((monitor & 0x1) << 0));
#if 0
  AliInfo(Form("d=%d, r=%c, b=%d, m=%d -> (%d<<3)|(%d<<2)|(%d<<1)|(%d<<0)=%2d",
	       det, ring, board, monitor, (det-1) & 0x3, iring & 0x1, 
	       board & 0x1, monitor & 0x1, index));
#endif
  return index-2;
}
//_____________________________________________________________________ 
void AliFMDQADataMakerRec::GetHalfringFromIndex(Int_t     idx, 
						UShort_t& det, 
						Char_t&   ring, 
						UShort_t& board, 
						UShort_t& monitor)
{
  det     = ((idx >> 3) & 0x3) + 1;
  ring    = ((idx >> 2) & 0x1) == 1 ? 'I' : 'O';
  board   = ((idx >> 1) & 0x1);
  monitor = ((idx >> 0) & 0x1);
}


//_____________________________________________________________________ 


//
// EOF
//
