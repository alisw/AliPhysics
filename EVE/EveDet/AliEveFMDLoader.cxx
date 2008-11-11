/**************************************************************************
 * Copyright(c) 2008, Christian Holm Christensen                          *
 *                                                                        *
 * Author: Christian Holm Christensen.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation for any purposes is hereby granted without fee,          * 
 * provided that the above copyright notice appears in all copies and     *
 * that both the copyright notice and this permission notice remains      *
 * intact.  The authors make no claims about the suitability of this      * 
 * software for any purpose. It is provided "as is" without express or    *
 * implied warranty.                                      		  * 
 **************************************************************************/
/* $Id$ */
/** @file    AliEveFMDLoader.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:59:18 2006
    @brief   Implementation of AliEveFMDLoader singleton class 
*/
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. This class
// is the loader for the event display. 
//
// This class is a singleton, meaning that there's only one instance
// of this.  This is done to speed up the processing by putting all
// things that are needed every time into the constructor. 
//
#include "AliRunLoader.h"
#include "EveBase/AliEveEventManager.h"
#include "AliEveFMDLoader.h"
#include "../FMD/AliFMDUShortMap.h"
#include "../FMD/AliFMDBoolMap.h"
#include "../FMD/AliFMDGeometry.h"
#include "../FMD/AliFMDParameters.h"
#include "../FMD/AliFMDDetector.h"
#include "../FMD/AliFMDRing.h"
#include "../FMD/AliFMDBaseDigit.h"
#include "../FMD/AliFMDDigit.h"
#include "../FMD/AliFMDRawReader.h"
#include "../FMD/AliFMDHit.h"
#include "AliESDEvent.h"
#include "AliESDFMD.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include <TClonesArray.h>
#include <TTree.h>
#include <TGeoShape.h>
#include <TGeoManager.h>
#include <TEveGeoNode.h>
#include <TEveBoxSet.h>
#include <TEveQuadSet.h>
#include <TEveManager.h>
#include <TEveUtil.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>

// Some very private variables. 
namespace 
{
  const Char_t* kDetector = "FMD%d";
  const Char_t* kRing     = "FMD%d%c";
  const Char_t* kModule   = "FMD%d%c[%02d-%02d]";
  const Char_t* kSector   = "FMD%d%c[%02d] %s";

  const Char_t* kHits     = "Hits";
  const Char_t* kDigits   = "Digits";
  const Char_t* kRaw      = "Raw";
  const Char_t* kESD      = "ESD";
  
}

//____________________________________________________________________
AliEveFMDLoader* AliEveFMDLoader::fgInstance = 0;

//____________________________________________________________________
AliEveFMDLoader* AliEveFMDLoader::Instance()
{
  // Get the singleton instance.  If the instance has not been
  // instantised yet, it will be after this call.
  if (!fgInstance) 
    fgInstance = new AliEveFMDLoader();
  return fgInstance;
}

//____________________________________________________________________
AliEveFMDLoader::AliEveFMDLoader(const char* name, Bool_t useBoxes, 
				 Bool_t old)
  : TEveElementList(name, 0), 
    fHitPalette(0, 1000),
    fDigitPalette(0, 1023), 
    fMultPalette(0, 20),
    fUseBoxDigits(useBoxes), 
    fHitCache("AliFMDHit",0),
    fDigitCache("AliFMDDigit", 0),
    fRawCache("AliFMDDigit", 0)
{
  // Constructor 
  // Parameters:
  // @param name     Name of the folder. 
  // @param useBoxes Whether to use boxes or Quads for the signals 
  // @param old      Whether to enable reading old RCU data format

  // increase reference count 
  IncDenyDestroy();
  
  // Increase reference counts on palettes 
  fHitPalette.IncRefCount();
  fDigitPalette.IncRefCount();
  fMultPalette.IncRefCount();
  
  
  // Initialize the FMD geometry manager 
  AliEveEventManager::AssertGeometry();
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();

  AliFMDParameters* pars = AliFMDParameters::Instance();
  // pars->UseRcuTrailer(!old);
  pars->UseCompleteHeader(old);
  pars->SetSampleRate(4);
  pars->Init(kFALSE, 0);

  // Get shapes
  TGeoShape* inner = static_cast<TGeoShape*>(gGeoManager->GetListOfShapes()
					     ->FindObject("FMDI_physical_sensor"));
  if (!inner) throw TEveException("Shape of inner type sensors not found");
  TGeoShape* outer = static_cast<TGeoShape*>(gGeoManager->GetListOfShapes()
					     ->FindObject("FMDO_physical_sensor"));
  if (!outer) throw TEveException("Shape of outer type sensors not found");

  // Emulate reference counting 
  inner->SetUniqueID(1000);
  outer->SetUniqueID(1000);
  
  // Loop over detectors 
  for (UShort_t d = 1; d <= 3; d++) { 
    AliFMDDetector*  detector = geom->GetDetector(d);
    if (!detector) continue;
    TEveElementList* ed       = new TEveElementList(Form(kDetector, 
							 detector->GetId()));
    AddElement(ed);
    ed->IncDenyDestroy();
    ed->SetUserData(detector);

    // Loop over rings 
    Char_t           rings[]  = { 'I', 'O', 0 };
    Char_t*          pr       = &(rings[0]);
    while (*pr) { 
      AliFMDRing* ring = detector->GetRing(*pr);
      pr++;
      if (!ring) continue;
      TEveElementList* er      = new TEveElementList(Form(kRing, 
							  detector->GetId(),
							  ring->GetId()));
      ed->AddElement(er);
      er->IncDenyDestroy();
      er->SetUserData(ring);
      
      // UShort_t      nsec    = ring->GetNSectors();
      // UShort_t      nstr    = ring->GetNStrips();
      UShort_t         nmod    = ring->GetNModules();
      // Loop over modules 
      for (UShort_t m = 0; m < nmod; m++) { 
	TEveGeoShape* em = new TEveGeoShape(Form(kModule, 
						 detector->GetId(),
						 ring->GetId(), 
						 2*m, 2*m+1));
	er->AddElement(em);
	em->SetTransMatrix(*(detector->FindTransform(ring->GetId(), 2*m)));
	em->SetShape(ring->GetId() == 'I' ? inner : outer);
	em->SetMainColor(Color_t(kRed));
	em->SetMainTransparency(32);
	em->IncDenyDestroy();

#if 0
	for (UShort_t s = 2*m; s < 2*m+2 && s < nsec; s++) { 
	  TEveDigitSet* eb = MakeDigitSet(Form(kSector,
					       detector->GetId(), 
					       ring->GetId(), s), nstr);
	  em->AddElement(eb);
	  eb->SetEmitSignals(kFALSE);
	  eb->SetPickable(kTRUE);
	  // eb->SetOwnIds(kTRUE);
	} // for (UShort_t s ...)
#endif
      }  // for (UShort_t m ...)
    }  // while (pr)
  }  // for (UShort_t d ...)
}

//____________________________________________________________________
AliEveFMDLoader::~AliEveFMDLoader()
{
  // Destructor
  AliWarning("AliEveFMDLoader being destroyed!");
}

//____________________________________________________________________
Int_t
AliEveFMDLoader::RemoveFromListTrees(TEveElement* el)
{
  // Called when the element should be removed from the list. We
  // overload this to allow clearing of signals.
  // Parameters:
  // @param el Tree to remove from.

  // Since we're most likely setting up for a new event, we clear all
  // signals here - a little tricky, but it works(tm). 
  ClearDigitSets("All");
  
  // Do normal TEveElement::RemoveElement
  return TEveElementList::RemoveFromListTrees(el);
}
//____________________________________________________________________
void
AliEveFMDLoader::RemoveParent(TEveElement* el)
{
  // Called when the element should be removed from the list. We
  // overload this to allow clearing of signals.
  // Parameters:
  // @param el Parent to remove from.
  TEveElementList::RemoveParent(el);
}

//____________________________________________________________________
TEveDigitSet* 
AliEveFMDLoader::MakeDigitSet(const char* name, UShort_t nstr)
{
  // Make a digit set.  The type of digit set depends on the setting
  // of fUseBoxDigits.  If this is true, we return a TEveBoxSet,
  // otherwise a TEveQuadSet   
  // Parameters: 
  //    name	The name 
  //	nstr	The number of strips 
  // Return 
  //    newly allocated digit set
  TEveDigitSet* ret = 0;
  if (fUseBoxDigits) { 
    TEveBoxSet* boxes = new TEveBoxSet(name);
    // boxes->Reset(TEveBoxSet::kBT_AABox, kFALSE, nstr);
    boxes->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, nstr);
    ret = boxes;
  }
  else { 
    TEveQuadSet* quads = new TEveQuadSet(name);
    quads->Reset(TEveQuadSet::kQT_RectangleXY, kFALSE, nstr);
    ret = quads;
  }
  return ret;
}

//____________________________________________________________________
void
AliEveFMDLoader::ClearDigitSets(const char* type)
{
  // Clear signals of some type.   
  // Parameters:
  // @param type Type of signals to clear 
  //     Type can be one of 
  //     - All    All signals 
  //     - Hits   Hits 
  //     - Digits Digits 
  //     - Raw    Raw 
  //     - ESD    ESD 
  TString stype(type);

  for (TEveElement::List_i di = BeginChildren(); 
       di != EndChildren(); ++di) { 
    for (TEveElement::List_i ri = (*di)->BeginChildren(); 
	 ri != (*di)->EndChildren(); ++ri) { 
      for (TEveElement::List_i mi = (*ri)->BeginChildren();
	   mi != (*ri)->EndChildren(); ++mi) { 
	if (stype == "All") {
	  (*mi)->RemoveElements();
	  continue;
	}
	for (TEveElement::List_i si = (*mi)->BeginChildren(); 
	     si != (*mi)->EndChildren(); ++si) { 
	  TEveDigitSet* signals = static_cast<TEveDigitSet*>((*si));
	  if (!signals) continue;
	  TString s(signals->GetName());
	  if (!s.Contains(type)) continue;
	  (*mi)->RemoveElement(signals);
	}
      }
    }
  }
}

//____________________________________________________________________
TEveDigitSet*
AliEveFMDLoader::FindDigitSet(const char* t, UShort_t d, Char_t r, UShort_t s)
{
  // Find a digit set corresponding to the passed parameters.  If it
  // is not found, one is created
  // Parameters:
  // @param type   Type of data 
  // @param d      Detector 
  // @param r      Ring 
  // @param s      Sector 
  // @return a digit set
  TEveElement* detector = FindChild(Form(kDetector, d));
  if (!detector) { 
    AliError(Form("Detector %s not found", Form(kDetector, d)));
    return 0;
  }
  
  TEveElement* ring = detector->FindChild(Form(kRing, d, r));
  if (!ring) { 
    AliError(Form("Ring %s not found", Form(kRing, d, r)));
    return 0;
  }
  
  Int_t mod = 2*(s/2);
  TEveElement* module = ring->FindChild(Form(kModule, d, r, mod, mod+1));
  if (!module) { 
    AliError(Form("Module %s not found", Form(kModule, d, r, s, s+1)));
    return 0;
  }
  
  TEveElement*  sector = module->FindChild(Form(kSector, d, r, s, t));
  TEveDigitSet* signal = static_cast<TEveDigitSet*>(sector);
  if (!sector) { 
    AliFMDRing* rng = AliFMDGeometry::Instance()->GetRing(r);
    signal = MakeDigitSet(Form(kSector, d, r, s, t), rng->GetNStrips());
    module->AddElement(signal);
    signal->SetEmitSignals(kFALSE);
    signal->SetPickable(kTRUE);
    TString st(t);
    if      (t == kHits)   signal->SetPalette(&fHitPalette);
    else if (t == kDigits) signal->SetPalette(&fDigitPalette);
    else if (t == kRaw)    signal->SetPalette(&fDigitPalette);
    else if (t == kESD) {
      signal->SetPalette(&fMultPalette);
      signal->SetOwnIds(kTRUE);
    }
  }
  return signal;
}

//____________________________________________________________________
void
AliEveFMDLoader::AddSignal(const char* t, 
			   UShort_t det, Char_t rng, UShort_t sec, 
			   UShort_t str, Float_t signal, Float_t min, 
			   Float_t  max, TObject* ref)
{
  // Add a signal to a digit set
  // Parameters:
  // @param type   Type of data 
  // @param det    Detector 
  // @param rng    Ring 
  // @param sec    Sector 
  // @param str    Strip
  // @param signal Signal value 
  // @param min    Minimum of this kind of signal 
  // @param max    Maximum of this kind of signal 
  // @param ref    Reference object 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  Double_t x, y, z;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);
  AddSignal(t, det, rng, sec, str, x, y, z, signal, min, max, ref);
}

//____________________________________________________________________
void
AliEveFMDLoader::AddSignal(const char* t, 
			   UShort_t det, Char_t rng, UShort_t sec, 
			   UShort_t str, Double_t x, Double_t y, Double_t z, 
			   Float_t signal, Float_t min, Float_t  max, 
			   TObject* ref)
{
  // Add a signal to a digit set, with known (x,y,z) coordinates
  // (this is for hits)
  // Parameters:
  // @param type   Type of data 
  // @param det    Detector 
  // @param rng    Ring 
  // @param sec    Sector 
  // @param str    Strip
  // @param x      X coordinate  
  // @param y      Y coordinate  
  // @param z      Z coordinate  
  // @param signal Signal value 
  // @param min    Minimum of this kind of signal 
  // @param max    Maximum of this kind of signal 
  // @param ref    Reference object 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  AliFMDRing*     ring = geom->GetRing(rng);
  if (!ring) return;

  TEveDigitSet* signals = FindDigitSet(t, det, rng, sec);
  if (!signals) { 
    AliWarning(Form("No signal (%s) found for FMD%d%c[%02d,%03d]", 
		    t, det, rng, sec, str));
    return;
  }
  
  Float_t  scaled = TMath::Min((signal - min) / (max - min) * 10., 10.);
  Double_t w      = 2*ring->GetPitch();
  Int_t    value  = int(TMath::Nint(signal));
  AliDebug(1, Form("New signal at FMD%d%c[%02d,%03d]=%f (v=%d, s=%f)", 
		   det, rng, sec, str, signal, value, scaled));
  AddDigit(signals, x, y, z, w, scaled, value, ref);
}

//____________________________________________________________________
void
AliEveFMDLoader::AddDigit(TEveDigitSet* signals, 
			  Double_t x, Double_t y, Double_t z, 
			  Double_t w, Float_t scaled, Int_t value, 
			  TObject* ref)
{
  // Add a digit to a digit set. 
  // Parameters:
  // @param signals Digit set. 
  // @param x      X coordinate  
  // @param y      Y coordinate  
  // @param z      Z coordinate  
  // @param w      strip pitch 
  // @param scaled Scaled value 
  // @param value  Signal value 
  // @param ref    Reference object
  if (fUseBoxDigits) { 
    TEveBoxSet* boxes = static_cast<TEveBoxSet*>(signals);
    Float_t zc   = (z > 0 ? -1 : 1) * scaled + z;
    Float_t vs[] = { -w, -5*w, zc-scaled,   // Lower back  left
		     +w, -5*w, zc-scaled,   // Lower back  right
		     +w, +5*w, zc-scaled,   // Lower front right
		     -w, +5*w, zc-scaled,   // Lower front left 
		     -w, -5*w, zc+scaled,   // Upper back  left
		     +w, -5*w, zc+scaled,   // Upper back  right
		     +w, +5*w, zc+scaled,   // Upper front right
		     -w, +5*w, zc+scaled }; // Upper front left
    Float_t ang  = TMath::ATan2(y,x);
    for (size_t i = 0; i < 8; i++) { 
      Float_t bx = vs[3*i+0];
      Float_t by = vs[3*i+1];
      Float_t ca = TMath::Cos(ang);
      Float_t sa = TMath::Sin(ang);
      vs[3*i+0]  = bx * ca - by * sa + x;
      vs[3*i+1]  = bx * sa + by * ca + y;
    }
    // boxes->AddBox(x, y, (z > 0 ? -scaled : 0) + z , 5*w, w, scaled);
    boxes->AddBox(vs);
    boxes->DigitValue(value);
    if (ref) boxes->DigitId(ref);
  }
  else { 
    TEveQuadSet* quads = static_cast<TEveQuadSet*>(signals);
    quads->AddQuad(x,y,z,w,w);
    quads->QuadValue(value);
    if (ref) quads->QuadId(ref);
  }
}


//____________________________________________________________________
void
AliEveFMDLoader::CheckAdd()
{
  // check if we shoul re-add ourselves to the current event node 
  TEveElement* event = gEve->GetCurrentEvent();
  if (event && event->FindChild(GetName())) return;
  gEve->AddElement(this);
}


//____________________________________________________________________
void
AliEveFMDLoader::LoadHits()
{
  // Load and display hits
  ClearDigitSets(kHits);

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  if (!rl) { 
    AliError("No run loader");
    return;
  }
  
  rl->LoadHits("FMD");
  TTree* ht = rl->GetTreeH("FMD", false);
  if (!ht) { 
    AliError("No FMD tree");
    return;
  }
  
  TClonesArray* hits = &fHitCache;
  fHitCache.Clear();
  ht->SetBranchAddress("FMD", &hits);

  Float_t min = fHitPalette.GetMinVal();
  Float_t max = fHitPalette.GetMaxVal();

  Int_t nTracks = ht->GetEntriesFast();
  for (Int_t i = 0; i < nTracks; i++) {
    Int_t hitRead  = ht->GetEntry(i);
    if (hitRead <= 0) continue;

    Int_t nHit = hits->GetEntriesFast();
    if (nHit <= 0) continue;
  
    for (Int_t j = 0; j < nHit; j++) {
      AliFMDHit* hit = static_cast<AliFMDHit*>(hits->At(j));
      if (!hit) continue;
  
      AddSignal(kHits, 
		hit->Detector(), hit->Ring(), hit->Sector(), hit->Strip(),
		hit->X(), hit->Y(), hit->Z(), int(hit->Edep()*1000), 
		min, max, hit);
    }
  }
  CheckAdd();
}

//____________________________________________________________________
void
AliEveFMDLoader::DoLoadDigits(const char* t, TClonesArray* digits)
{
  // Do the actual display of digits 
  // Parameters:
  // @param type What to show 
  // @param digits The digits
  Float_t min = fDigitPalette.GetMinVal();
  Float_t max = fDigitPalette.GetMaxVal();

  Int_t n = digits->GetEntriesFast();
  for (Int_t i = 0; i < n; i++) { 
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    if (!digit) return;
    AddSignal(t, digit->Detector(), digit->Ring(), digit->Sector(),
	      digit->Strip(), digit->Counts(), min, max, digit);
  }
  CheckAdd();
}

//____________________________________________________________________
void
AliEveFMDLoader::LoadDigits()
{
  // Load and display simulated digits 
  ClearDigitSets(kDigits);

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  if (!rl) { 
    AliError("No run-loader");
    return;
  }
  
  rl->LoadDigits("FMD");
  TTree* dt = rl->GetTreeD("FMD", false);
  if (!dt) { 
    AliError("No FMD tree");
    return;
  }
  
  TClonesArray* digits = &fDigitCache;
  fDigitCache.Clear();
  dt->SetBranchAddress("FMD", &digits);

  Int_t read = dt->GetEntry(0);
  if (read <= 0) { 
    AliWarning("Nothing read");
    return;
  }
  DoLoadDigits(kDigits, digits);
}


//____________________________________________________________________
void
AliEveFMDLoader::LoadRaw()
{
  // Load and display raw digits 
  ClearDigitSets(kRaw);

  AliRawReader* rr =  AliEveEventManager::AssertRawReader();
  if (!rr) { 
    AliError("No raw-reader");
    return;
  }
  rr->Reset();
  std::cout<<"Now in event # " << *(rr->GetEventId()) << std::endl;
  AliFMDRawReader* fr = new AliFMDRawReader(rr, 0);
  TClonesArray* digits = &fRawCache;
  fRawCache.Clear();
  
  fr->ReadAdcs(digits);
  
  DoLoadDigits(kRaw, digits);
}

//____________________________________________________________________
void
AliEveFMDLoader::LoadESD()
{
  // Load and display ESD information 
  ClearDigitSets(kESD);

  AliESDEvent* esd =  AliEveEventManager::AssertESD();
  if (!esd) { 
    AliError("No ESD");
    return;
  }

  AliESDFMD* fmd = esd->GetFMDData();
  if (!fmd) { 
    AliError("No FMD ESD data");
    return;
  }

  Float_t min = fMultPalette.GetMinVal();
  Float_t max = fMultPalette.GetMaxVal();
  
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* rng = rings; *rng != '\0'; rng++) {
      UShort_t nsec = (*rng == 'I' ?  20 :  40);
      UShort_t nstr = (*rng == 'I' ? 512 : 256);
      for (UShort_t sec = 0; sec < nsec; sec++) {
	for (UShort_t str = 0; str < nstr; str++) {
	  Float_t mult = fmd->Multiplicity(det,*rng,sec,str);
	  if (mult == AliESDFMD::kInvalidMult) continue;
	  Float_t eta  = fmd->Eta(det,*rng,sec,str);
	  AddSignal(kESD, det, *rng, sec, str, mult, min, max, 
		    new TNamed(Form("FMD%d%c[%02d,%03d]", det, *rng, sec, str), 
			       Form("Mch=%f, eta=%f", mult, eta)));
	}
      }
    }
  }
  CheckAdd();
}


	
      

  
//____________________________________________________________________
//
// EOF
//
