#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "AliParticleYield.h"
#include "TDatabasePDG.h"
#include "AliLog.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "AliPDG.h"
#include "TBranch.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TEventList.h"

using std::endl;
using std::left;
using std::setw;

ClassImp(AliParticleYield)

// set statics
const char * AliParticleYield::kStatusString[] = {"Published", "Preliminary", "Final, but not published", "May change"} ;
const char * AliParticleYield::kSystemString[] = {"pp", "p-Pb", "Pb-Pb"} ;
Int_t AliParticleYield::fSignificantDigits = 3;
Float_t AliParticleYield::fEpsilon = 0.000000000000000001;

AliParticleYield::AliParticleYield() :
TObject(),
fPdgCode(0),
fPdgCode2(0),
fPartName(""),
fCollisionSystem(0),
fSqrtS(0),
fYield(0),
fStatError(0),
fSystError(0),
fNormErrorPos(0),
fNormErrorNeg(0),
fYMin(0),
fYMax(0),
fStatus(0),
fMeasurementType(0),
fCentr(""),
fIsSum(0),
fTag("")
{
  // constructor
  AliPDG::AddParticlesToPdgDataBase(); // Make sure that ALICE-defined particles were added to the PDG DB
}

AliParticleYield::AliParticleYield(Int_t pdg, Int_t system, Float_t sqrts, Float_t value, Float_t stat, Float_t syst, Float_t norm, Float_t ymin, Float_t ymax, Int_t status, Int_t type, TString centr, Int_t isSum, TString tag):
TObject(),
fPdgCode(pdg),
fPdgCode2(0),
fPartName(""),
fCollisionSystem(system),
fSqrtS(sqrts),
fYield(value),
fStatError(stat),
fSystError(syst),
fNormErrorPos(norm),
fNormErrorNeg(0),
fYMin(ymin),
fYMax(ymax),
fStatus(status),
fMeasurementType(type),
fCentr(centr),
fIsSum(isSum),
fTag(tag)

{
  // Constructor
  fPartName = TDatabasePDG::Instance()->GetParticle(fPdgCode)->GetName();
  AliPDG::AddParticlesToPdgDataBase(); // Make sure that ALICE-defined particles were added to the PDG DB
}

AliParticleYield::AliParticleYield(Int_t pdg, Int_t system, Float_t sqrts, Float_t value, Float_t stat, Float_t syst, Float_t normPos, Float_t normNeg, Float_t ymin, Float_t ymax, Int_t status, Int_t type, TString centr, Int_t isSum, TString tag):
TObject(),
fPdgCode(pdg),
fPdgCode2(0),
fPartName(""),
fCollisionSystem(system),
fSqrtS(sqrts),
fYield(value),
fStatError(stat),
fSystError(syst),
fNormErrorPos(normPos),
fNormErrorNeg(normNeg),
fYMin(ymin),
fYMax(ymax),
fStatus(status),
fMeasurementType(type),
fCentr(centr),
fIsSum(isSum),
fTag(tag)

{
  // Constructor
  fPartName = TDatabasePDG::Instance()->GetParticle(fPdgCode)->GetName();
  AliPDG::AddParticlesToPdgDataBase(); // Make sure that ALICE-defined particles were added to the PDG DB
}


AliParticleYield::AliParticleYield(const AliParticleYield& part) : 
TObject(),
fPdgCode(part.fPdgCode),
fPdgCode2(part.fPdgCode2),
fPartName(part.fPartName),
fCollisionSystem(part.fCollisionSystem),
fSqrtS(part.fSqrtS),
fYield(part.fYield),
fStatError(part.fStatError),
fSystError(part.fSystError),
fNormErrorPos(part.fNormErrorPos),
fNormErrorNeg(part.fNormErrorNeg),
fYMin(part.fYMin),
fYMax(part.fYMax),
fStatus(part.fStatus),
fMeasurementType(part.fMeasurementType),
fCentr(part.fCentr),
fIsSum(part.fIsSum),
fTag(part.fTag){
  // Copy constructor
}

AliParticleYield::~AliParticleYield() {

}

TTree * AliParticleYield::ReadFromASCIIFileAsTree(const char * fileName, const char * separators){
  // Read the table from an ASCII File and returns a tree of particles. See ReadFromASCIIFile for detailed info on the format
  TClonesArray * arr = ReadFromASCIIFile(fileName, separators);
  AliParticleYield * part = 0;
  TTree * tree = new TTree ("treePart", "Particle Yields and Ratios");
  tree->Branch("particles", &part);
  TIter iterPart (arr);
  while ((part = (AliParticleYield*) iterPart.Next())){
    tree->Fill();
  }
  
  delete arr;
  delete part;
  return tree;
}

TClonesArray * AliParticleYield::GetEntriesMatchingSelection(TTree * tree, TString selection) {
  // Returns an array of particles from a tree created with ReadFromASCIIFileAsTree matching the selection. You can use the normal tree sintax for the selection, e.g. "fCentr == \"V0M0010\" && fStatus == 0".

  TClonesArray * arr = new TClonesArray("AliParticleYield");
  AliParticleYield * part = 0;
  tree->SetBranchAddress("particles", &part);
  // In order to get the array, we first create an entry list matching the selection using TTree::Draw, and them we loop over all entries in the tree.
  tree->Draw(">>particlelist", selection);// Produce selection list
  TEventList *elist = (TEventList*)gDirectory->Get("particlelist");
  Int_t npart = elist->GetN();
  for(Int_t ipart = 0; ipart < npart; ipart++){
    tree->GetEntry(elist->GetEntry(ipart));
    new((*arr)[ipart]) AliParticleYield(*part);// We need to clone part, because it is overwritten by the next read
  }
  elist->Delete();
  return arr;
}


TClonesArray * AliParticleYield::ReadFromASCIIFile(const char * fileName, const char * separators){
  // Read the table from an ASCII File with the format indicated
  // below. Returns a TClonesArray of AliParticleyields with the
  // content of the lines. Lines beginning by "#" are skipped.
  // The order of the columns is compulsory, but the separator can be set (by default whitespaces are assumed).

  // The columns should be:
  // PDG   NAME  SYSTEM  SQRTS  VALUE  SYST  STAT  NORM  YMIN  YMAX  STATUS  TYPE  CENTR     ISSUM TAG

  // PDG should either be an integher or the ratio of two integers (in case of particle ratios), with the following format:
  //  pdgcode/pdgcode2
  // NO SPACES ARE ALLOWED IN NAMES AND PDG CODE, unless you use a separator which is not a whitespace

  // A Header can be present (lines beginning with the word "PDG" are also skipped

  const Int_t kNCols = 14; // The lines are actually 15, but the last one (TAG) can be empty, so we put 14 here.

  TClonesArray * arr = new TClonesArray ("AliParticleYield");
  ifstream filestream (fileName);
  if(!filestream.is_open()) {
    Printf("Cannot open file %s\n", fileName);
    exit(1);
  }
  TString line;
  Int_t ipart = 0;
  std::cout << "Reading " << fileName << std::endl;
  
  while (line.ReadLine(filestream) ) {
    // Strip trailing and leading whitespaces
    line = line.Strip(TString::kLeading,  ' ');
    line = line.Strip(TString::kTrailing, ' ');

    // Skip commented lines and headers
    if (line.BeginsWith("#")) {
      //print comments. It if they look like warnings, print them such that they are really visible
      if(line.Contains("warn", TString::kIgnoreCase)) std::cout << std::endl << "********************************************************" <<std::endl ;
      std::cout << " " << line.Data() << std::endl;      
      if(line.Contains("warn", TString::kIgnoreCase)) std::cout << "********************************************************" <<std::endl << std::endl;

      continue;
    }
    if (line.BeginsWith("PDG")) continue;

    // Tokenize line using custom separator
    TObjArray * cols = line.Tokenize(separators);

    // Check the number of columns
    if(cols->GetEntries() < kNCols) {
      Printf("Wrong number of columns in table %d vs %d expected" , cols->GetEntries(), kNCols);
      delete arr;
      return NULL;
    }

    // Get Values
    // get type first, as some operations are type-specific
    UInt_t  type   = ((TObjString*)cols->At(11)) ->String().Atoi();

    // if it's a ratio, try to get the 2 pdg codes
    Int_t pdg =0, pdg2 = 0;
    
    if (type & kTypeParticleRatio) {      
      TString col0 = ((TObjString*)cols->At(0))  ->String();
      TObjArray * tokens = col0.Tokenize("/");
      if(tokens->GetEntries() != 2) {
	Printf("ERROR: Cannot get both PDGs for ratios");		
      } else {
	pdg  = ((TObjString*)tokens->At(0))  ->String().Atoi();
	pdg2 = ((TObjString*)tokens->At(1))  ->String().Atoi();
      }
    }
    else {
      pdg    = ((TObjString*)cols->At(0))  ->String().Atoi();
    }
    TString name   = ((TObjString*)cols->At(1))  ->String();
    Int_t   system = ((TObjString*)cols->At(2))  ->String().Atoi();
    Float_t sqrts  = ((TObjString*)cols->At(3))  ->String().Atof();
    Float_t yield  = ((TObjString*)cols->At(4))  ->String().Atof();
    // The "GetError" function can handle % errors. 
    Float_t stat   = GetError(((TObjString*)cols->At(5))  ->String(), yield);
    Float_t syst   = GetError(((TObjString*)cols->At(6))  ->String(), yield);
    TString normString(((TObjString*)cols->At(7))->String());

    Float_t normPos = 0;
    Float_t normNeg = 0;
    if (normString.Contains("+") && normString.Contains("-")) {
      
      // If the string for the normalization uncertainty contains a + and a -, it means it is asymmetric
      if(normString.First("+") < normString.First("-") ) {// the + error is quoted first
        normPos = GetError(normString(1,normString.First("-")-1)+normString(normString.First("e"),normString.Length()), yield); // start from 1 (skip + sign). The second bit is to propagate the scientific notation to the first part of the error
        normNeg = GetError(normString(normString.First("-")+1,normString.Length()), yield); // +1 -> skip sign
      } 
      else {
        // This is the opposite case
        normNeg = GetError(normString(1,normString.First("+")-1)+normString(normString.First("e"),normString.Length()), yield); // start from 1 (skip + sign). The second bit is to propagate the scientific notation to the first part of the error
        normPos = GetError(normString(normString.First("+")+1,normString.Length()), yield); // +1 -> skip sign
      }
      
    } else {
      // symmetric error: set only normpos
      normPos   = GetError(((TObjString*)cols->At(7))  ->String(), yield);
    }
    Float_t ymin   = ((TObjString*)cols->At(8))  ->String().Atof();
    Float_t ymax   = ((TObjString*)cols->At(9))  ->String().Atof();
    Int_t   status = ((TObjString*)cols->At(10)) ->String().Atoi();
    TString centr  = ((TObjString*)cols->At(12)) ->String();
    Int_t   issum  = ((TObjString*)cols->At(13)) ->String().Atoi();    
    TString tag    = cols->At(14) ? ((TObjString*)cols->At(14)) ->String() : ""; // tag can be empty

    // Cleanup strings
    name  = name.Strip(TString::kLeading,  ' ');
    name  = name.Strip(TString::kTrailing, ' ');
    centr = centr.Strip(TString::kLeading,  ' ');
    centr = centr.Strip(TString::kTrailing, ' ');
    tag   = tag.Strip(TString::kLeading,  ' ');
    tag   = tag.Strip(TString::kTrailing, ' ');
    
    // add to array
    AliParticleYield * part = new  AliParticleYield(pdg,system,sqrts,yield,stat,syst,normPos, normNeg,ymin,ymax,status,type,centr,issum,tag);
    part->SetPartName(name); // Check name and PDG code consistency   
    part->SetPdgCode2(pdg2); // Set second PDG code in case of ratios 
    part->CheckTypeConsistency();                                     
    if(!part->CheckForDuplicates(arr)) {
      new ((*arr)[ipart++]) AliParticleYield(*part); 
    }
      //    delete part;

  }
  std::cout << "<- File read" << std::endl;


  return arr;
}

const char * AliParticleYield::GetLatexName(Int_t pdg) const {
  
  // Returns a TLatex compatible name for the particle
  // if pdg == 0 uses fPdgcode;
  // We need the pdg argument for particle ratios

  if(!pdg && fMeasurementType & kTypeParticleRatio) {
    // If it's a ratio, we try to build the ratio name. To avoid an infinite loop we have to call GetLatexname with a non-zero argument.
    static TString name;
    name ="";
    name += GetLatexName(fPdgCode);
    name += " / ";
    name += GetLatexName(fPdgCode2);
    return name.Data();
  }

  if(!pdg) pdg = fPdgCode;

  switch (pdg) {
  case 211:
    if (fIsSum) return "(#pi^{+} + #pi^{-})";
    return "#pi^{+}";
  case -211:
    return "#pi^{-}";
  case 321:
    if (fIsSum) return "(K^{+} + K^{-})";
    return "K^{+}";
  case -321:
    return "K^{-}";
  case 2212:
    if (fIsSum) return "(p + #bar{p})";
    return "p";
  case -2212:
    return "#bar^{p}";
  case 3122:
    if (fIsSum) return "(#Omega^{-} + #bar{#Omega}^{+})";
    return "#Omega^{-}";
  case -3122:
    return "#bar{#Omega}^{+}";
  case 3312:
    if (fIsSum) return "(#Xi^{-} + #bar{#Xi}^{+})";
    return "#Xi^{-}";
  case -3312:
    return "#bar{#Xi}^{+}";
  case 3334:
    if (fIsSum) return "(#Omega^{-} + #bar{#Omega}^{+})";
    return "#Omega^{-}";
  case -3334:
    return "#bar{#Omega}^{+}";
  case 310:
    return "K^{0}_{S}";
  case 333:
    return "#phi";
  case 313:
    return "K^{*}";
  case 1000010020:
    if(fIsSum) return "(d + #bar{d})";
    return "d";// Deuteron
  case -1000010020:
    return "#bar{d}";// Deuteron
  case 1000020030:
    if(fIsSum) return "(^{3}He + #bar{^{3}He})";
    return "^{3}He";
  case -1000020030:
    return "#bar{^{3}He}";
  case 1010010030:
    if(fIsSum) return "^{3}_{#Lambda}H + #bar{^{3}_{#Lambda}H}";
    return "^{3}_{#Lambda}H";
  case -1010010030:
    return "#bar{^{3}_{#Lambda}H}";    
  }
  AliWarning("Latex Name not know for this particle");
  return fPartName.Data();

}

Float_t AliParticleYield::GetTotalError(Bool_t includeNormalization) const {
  // Returns the total error, including or not the normalization uncertainty
  // All uncertainties are supposed to be uncorrelated (e.g. summed in quadrature)
  // If stat and syst are stored separately, the total error is computed summing them in quadrature
  Float_t error = GetSystError();
  if (!(fMeasurementType & kTypeOnlyTotError)) error = TMath::Sqrt(error*error + GetStatError()*GetStatError());
  if(includeNormalization) error = TMath::Sqrt(error*error + GetNormError()*GetNormError());
  
  return error;


}


void AliParticleYield::SaveAsASCIIFile(TClonesArray * arr, const char * fileName, const char * separator, Int_t colWidth){
  // Saves the array as an ASCII File with the format indicated
  // below. 

  // The columns should be:
  // PDG   NAME  SYSTEM  SQRTS  VALUE  STAT  SYST  NORM  YMIN  YMAX  STATUS  TYPE  CENTR     ISSUM TAG
  if(!arr) {
    Printf("<AliParticleYield::SaveAsASCIIFile> Error: no array provided");
    return;
  }
  if(!fileName) {
    Printf("<AliParticleYield::SaveAsASCIIFile> Error: no filename provided");
  }


  ofstream fileOut(fileName);
  //print header
  fileOut << FormatCol("PDG", colWidth, separator) <<   FormatCol("NAME", colWidth, separator) <<  FormatCol("SYSTEM", colWidth, separator) <<  FormatCol("SQRTS", colWidth, separator) <<  FormatCol("VALUE", colWidth, separator) <<  FormatCol("STAT" , colWidth, separator)<<  FormatCol("SYST", colWidth, separator) <<  FormatCol("NORM", colWidth, separator) <<  FormatCol("YMIN", colWidth, separator) <<  FormatCol("YMAX", colWidth, separator) <<  FormatCol("STATUS", colWidth, separator) <<  FormatCol("TYPE", colWidth, separator) <<  FormatCol("CENTR", colWidth, separator) <<     FormatCol("ISSUM", colWidth, separator) <<  FormatCol("TAG", colWidth, separator) << endl;
  

  // This is used for float numbers in the table.
  // The "g" options switches between the normal or scientific notation, whathever is more appropriate.
  // We want to have up to fSignificantDigits digits after the .
  char format[20];
  snprintf(format,20,"%%%dg", fSignificantDigits);

  char formatA[30];// We have to rebuild the format for asymmetric uncertainties...
  snprintf(formatA,30,"+%%%dg-%%%dg", fSignificantDigits, fSignificantDigits);

  TIter iter(arr);
  AliParticleYield * part = 0;
  TString normError ;
  while ((part = (AliParticleYield*) iter.Next())){    
    if(part->GetNormErrorNeg()) {
      normError = FormatCol(Form(formatA, // Asymmetric error format  
                                 RoundToSignificantFigures(part->GetNormErrorPos(),fSignificantDigits), 
                                 RoundToSignificantFigures(part->GetNormErrorNeg(),fSignificantDigits)),
                            colWidth,
                            separator);
    }
    else {
      normError = FormatCol(Form(format, RoundToSignificantFigures(part->GetNormError(),fSignificantDigits)) , colWidth , separator);
    }
    fileOut 
      << FormatCol(Form("%d",part->GetPdgCode())                                                    , colWidth , separator) 
      << FormatCol(part->GetPartName()                                                              , colWidth , separator) 	    
      << FormatCol(Form("%d", part->GetCollisionSystem())                                           , colWidth , separator) 
      << FormatCol(Form(format, part->GetSqrtS())                                                   , colWidth , separator)	    
      << FormatCol(Form(format, RoundToSignificantFigures(part->GetYield(),    fSignificantDigits)) , colWidth , separator)
      << FormatCol(Form(format, RoundToSignificantFigures(part->GetStatError(),fSignificantDigits)) , colWidth , separator) 
      << FormatCol(Form(format, RoundToSignificantFigures(part->GetSystError(),fSignificantDigits)) , colWidth , separator)
      << normError.Data()	
      << FormatCol(Form(format, part->GetYMin())                                                    , colWidth , separator) 
      << FormatCol(Form(format, part->GetYMax())                                                    , colWidth , separator)	    
      << FormatCol(Form("%d",part->GetStatus()          )                                           , colWidth , separator) 
      << FormatCol(Form("%d",part->GetMeasurementType() )                                           , colWidth , separator)       
      << FormatCol(part->GetCentr()                                                                 , colWidth , separator) 
      << FormatCol(Form("%d",part->GetIsSum())                                                      , colWidth , separator) 
      << FormatCol(part->GetTag()                                                                   , colWidth , separator) 
      << endl;
  }


}

void AliParticleYield::WriteThermusFile(TClonesArray * arr, const char * filename, Int_t colwidth) {
  // Writes a txt file which can we used as input in therums fits

  if(!arr) {
    Printf("<AliParticleYield::WriteThermusFile> Error: no array provided");
    return;
  }
  if(!filename) {
    Printf("<AliParticleYield::WriteThermusFile> Error: no filename provided");
    return;
  }

  ofstream fileOut(filename);

  TIter iter(arr);
  AliParticleYield * part = 0;
  char format[20];
  // This is used for float numbers in the table.
  // The "g" options switches between the normal or scientific notation, whathever is more appropriate.
  // We want to have up to fSignificantDigits digits after the .
  snprintf(format,20,"%%%dg", fSignificantDigits);
  
  //  snprintf(format, 20, "%d.%d%%f", fSignificantDigits, fSignificantDigits);
  while ((part = (AliParticleYield*) iter.Next())){    
    
    if(part->IsTypeRatio()) { 
      // If it's a ratio we have to write the 2 pdg codes
      fileOut << FormatCol(Form("%d %d ",part->GetPdgCode(), part->GetPdgCode2())                              , colwidth) 
	      << FormatCol(part->GetTag()                                                                      , colwidth)
	      << FormatCol(Form(format, RoundToSignificantFigures(part->GetYield()      , fSignificantDigits)) , colwidth)
	      << FormatCol(Form(format, RoundToSignificantFigures(part->GetTotalError() , fSignificantDigits)) , colwidth)
	      << endl;
    }
    else {
      fileOut << FormatCol(Form("%d",part->GetPdgCode())                                                       , colwidth) 
	      << FormatCol(part->GetTag()                                                                      , colwidth)
	      << FormatCol(Form(format, RoundToSignificantFigures(part->GetYield()      , fSignificantDigits)) , colwidth)
	      << FormatCol(Form(format, RoundToSignificantFigures(part->GetTotalError() , fSignificantDigits)) , colwidth)
	      << endl;      
    }
  
  }
  
}


const char * AliParticleYield::FormatCol(const char * text, Int_t width,  const char * sep) {
  
  TString format(Form("%%-%ds %s", width, sep));
  return Form(format.Data(), text);

}

Double_t AliParticleYield::RoundToSignificantFigures(double num, int n) {
  // Rounds num to n significant digits.
  // Recipe from :http://stackoverflow.com/questions/202302/rounding-to-an-arbitrary-number-of-significant-digits
  // Basically the log is used to determine the number of leading 0s, than convert to an integer by multipliing by the expo, 
  // round the integer and shift back.
  if(num == 0) {
    return 0;
  }

  Double_t d = TMath::Ceil(TMath::Log10(num < 0 ? -num: num));
  Int_t power = n - (int) d;

  Double_t magnitude = TMath::Power(10, power);
  Long_t shifted = TMath::Nint(num*magnitude);
  return shifted/magnitude;

}


Float_t AliParticleYield::GetError(TString error, Float_t yield) {
  // The "GetError" function can handle % errors. 
  if(error.Contains("%")) {
    return yield * error.Atof()/100;
  }
  return error.Atof();
}

void AliParticleYield::SetPdgCode(TString partName) {
  // Set pdg code from part name, if there was a previous name, check if it is consistent
  TParticlePDG * part  = TDatabasePDG::Instance()->GetParticle(partName);
  if(IsTypeRatio() || fIsSum) return; // Name check does not make sense for ratios and sums
  if(!part) {
    AliError(Form("No particle %s in TDatabasePDG", partName.Data()));
    return;
  }
  if(fPdgCode != part->PdgCode() &&  !(fMeasurementType&kTypeParticleRatio)) { // The consistency check on PDG codes is disabled case of ratios
    AliError(Form("Name and pdg code are not consistent! fPartName: %s partName %s, pdgcode %d fPdgCode %d", fPartName.Data(), partName.Data(), part->PdgCode(), fPdgCode));
  }
  fPdgCode = part->PdgCode();

}

void AliParticleYield::SetPartName(Int_t pdgCode) {
  // Set part name from pdg code, if there was a previous code, check if it is consistent
  TParticlePDG * part  = TDatabasePDG::Instance()->GetParticle(pdgCode);
  if(IsTypeRatio() || fIsSum) return; // Name check does not make sense for ratios and sums
  if(!part) {
    AliError(Form("No particle with code %d in TDatabasePDG", pdgCode));
    return;
  }
  if(fPdgCode != part->PdgCode() && !(fMeasurementType&kTypeParticleRatio)) { // The consistency check on particle names is disabled case of ratios
    AliError(Form("Name and pdg code are not consistent! fPartName: %s partName %s, pdgcode %d fPdgCode %d", fPartName.Data(), part->GetName(), part->PdgCode(), fPdgCode));
  }
  fPartName = part->GetName();

}

Bool_t AliParticleYield::CheckTypeConsistency() const {
  // Check for any inconsistency with the type mask. Returns true if the object is fully consistent.
  Bool_t isOK = kTRUE;

  if((fMeasurementType & kTypeOnlyTotError) && GetStatError()) {
    AliError(Form("Error flagged as total only, but stat error is not 0 (%f)!",GetStatError()));
    isOK= kFALSE;
  } else if (!(fMeasurementType & kTypeOnlyTotError) && (!GetStatError() || !GetSystError())) {
    AliError("Stat or syst errors are null");
  } 

  return isOK;
}

void AliParticleYield::Print (Option_t *) const {
  // Print
  Printf("-------------------------------");
  CheckTypeConsistency();
  Printf("%s [%s] (%d) %s %s", fPartName.Data(), GetLatexName(), fPdgCode, fIsSum ? "particle + antiparticle" : "", fTag.Length() ? Form("[%s]", fTag.Data()) : "" );
  Printf("Status: %s, %s", kStatusString[fStatus], fMeasurementType & kTypeLinearInterpolation ? "Interpolated" : "Measured");
  Printf("%s , sqrt(s) = %2.2f GeV, %2.2f < y < %2.2f %s", kSystemString[fCollisionSystem], fSqrtS, fYMin, fYMax, fCentr.Data()); 
  if(fMeasurementType & kTypeOnlyTotError) {
    Printf("%f +- %f (total error)", fYield, fSystError);
  }
  else {
    Printf("%f +- %f (stat) +- %f (syst)", fYield, fStatError, fSystError);
  }
  if(fNormErrorNeg) {
    Printf("Normalizaion uncertainty: +%f-%f", fNormErrorPos, fNormErrorNeg);    
  }
  else {
    Printf("Normalizaion uncertainty: %f", fNormErrorPos);
  }
}

Bool_t AliParticleYield::operator==(const AliParticleYield& rhs) {
  // Check if the two particles are identical

  Bool_t isEqual = 
    (fPdgCode         == rhs.fPdgCode              ) &&
    (fPdgCode2        == rhs.fPdgCode2             ) &&
    (fPartName        == rhs.fPartName             ) &&
    (fCollisionSystem == rhs.fCollisionSystem      ) &&
    Compare2Floats(fSqrtS,rhs.fSqrtS               ) &&
    Compare2Floats(fYield,rhs.fYield               ) &&
    Compare2Floats(fStatError,rhs.fStatError       ) &&
    Compare2Floats(fSystError,rhs.fSystError       ) &&
    Compare2Floats(fNormErrorPos,rhs.fNormErrorPos ) &&
    Compare2Floats(fNormErrorNeg,rhs.fNormErrorNeg ) &&
    Compare2Floats(fYMin,rhs.fYMin                 ) &&
    Compare2Floats(fYMax,rhs.fYMax                 ) &&
    (fStatus          == rhs.fStatus               ) &&
    (fMeasurementType == rhs.fMeasurementType      ) &&
    (fCentr           == rhs.fCentr                ) &&
    (fIsSum           == rhs.fIsSum                ) &&
    (fTag             == rhs.fTag                  ) ;
  
  return isEqual;
  
}
Bool_t AliParticleYield::IsTheSameMeasurement(AliParticleYield &rhs) {

  // Check the two particles represent the same measurement (independently of the values)
  Bool_t isEqual = 
    (fPdgCode         == rhs.fPdgCode         ) &&
    (fPdgCode2        == rhs.fPdgCode2        ) &&
    (fCollisionSystem == rhs.fCollisionSystem ) &&
    Compare2Floats(fSqrtS,rhs.fSqrtS          ) &&
    Compare2Floats(fYMin,rhs.fYMin            ) &&
    Compare2Floats(fYMax,rhs.fYMax            ) &&
    (fStatus          == rhs.fStatus          ) &&
    (fCentr           == rhs.fCentr           ) &&
    (fIsSum           == rhs.fIsSum           ) &&
    (fTag             == rhs.fTag             ) ;
  
  return isEqual;


}

Bool_t AliParticleYield::CheckForDuplicates(TClonesArray * arr) {

  // loop over all elements on the array and check for duplicates
  TIter iter(arr);
  AliParticleYield * part = 0;
  Bool_t isDuplicate = kFALSE;

  while ((part = (AliParticleYield*) iter.Next())) {
    if (IsTheSameMeasurement(*part)){
      AliWarning("Duplicated measurement found");
      isDuplicate = kTRUE;
      if (!((*this) == (*part))) {
	part->Print();
	Print();
	AliFatal("The 2 particles are different!");
      }
    }
  }
  return isDuplicate;

} 
 
Bool_t AliParticleYield::Compare2Floats(Float_t a, Float_t b) {
  // just a simple helper for the comparison methods
  if(!a) {
    if(!b) return kTRUE; // They are both 0;
    return kFALSE;// return here to avoid division by 0
  }
  Bool_t areEqual = (TMath::Abs((a - b)/a) < fEpsilon); // If the relative difference is < epsilon, returns true
  if(!areEqual) {
    Printf("Warning: %f and %f are different", a,b); 
  }
  return areEqual;
}


Float_t AliParticleYield::GetNormError() const {
  // Returs a symmetrized error in case the normalizatione Error is asymmetric
  if(fNormErrorNeg) {
    AliWarning("Error is asymmetric, returining symmetrized uncertainty");
    return (TMath::Abs(fNormErrorNeg)+TMath::Abs(fNormErrorPos))/2;
  }
  else return fNormErrorPos; // If the uncertainty is not asymmetric, fNormErrorPos stores it.

}
