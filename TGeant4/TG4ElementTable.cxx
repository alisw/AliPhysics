// $Id$
// Category: geometry
//
// See the class description in the header file.

#include "TG4ElementTable.h"
#include "TG4Globals.h"

#include <G4Element.hh>

// static data members

TG4ElementTable* TG4ElementTable::fgInstance = 0;

// lifecycle

TG4ElementTable::TG4ElementTable() {
// 
  Construct();
}

TG4ElementTable::TG4ElementTable(const TG4ElementTable& right) { 
//
  TG4Globals::Exception(
    "Attempt to copy TG4ElementTable singleton.");
}

TG4ElementTable::~TG4ElementTable(){
//
}

// operators

TG4ElementTable& TG4ElementTable::operator=(const TG4ElementTable& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4ElementTable singleton.");
    
  return *this;  
}    
          
// static methods
  
TG4ElementTable* TG4ElementTable::Instance() 
{
// singleton access method
// ---

  if (fgInstance == 0 ) {
    fgInstance = new TG4ElementTable();
  }
  
  return fgInstance;
}

// private methods

void TG4ElementTable::Construct() 
{ 
// construct element table
// ---

  // new G4Element(name, symbol, z, a);
  // check names
  new G4Element("Hydrogen",  "H",   1.,  1.00797*g/mole);
  new G4Element("Helium",    "He",  2.,  4.00260*g/mole);
  new G4Element("Lithium",   "Li",  3.,  6.941*g/mole);
  new G4Element("Berylium",  "Be",  4.,  9.01218*g/mole);
  new G4Element("Bohr",      "B",   5.,  10.811*g/mole);
  new G4Element("Carbon",    "C",   6.,  12.01115*g/mole);
  new G4Element("Nitrogen",  "N",   7.,  14.0067*g/mole);
  //in periodic table
  //new G4Element("Oxygen",    "O",   8.,  15.9994*g/mole);
  new G4Element("Oxygen",    "O",   8.,  16.00*g/mole);
  new G4Element("Fluorine",  "F",   9.,  18.99840*g/mole);
  new G4Element("Neon",      "Ne", 10.,  20.179*g/mole);
  new G4Element("Sodium",    "Na", 11.,  22.98977*g/mole);
  new G4Element("Magnesium", "Mg", 12.,  24.305 *g/mole);
  new G4Element("Aluminium", "Al", 13.,  26.98154*g/mole);
  new G4Element("Silicon",   "Si", 14.,  28.086*g/mole);
  new G4Element("Phosphorus","P",  15.,  30.97376*g/mole);
  new G4Element("Sulfur",    "S",  16.,  32.064*g/mole);
  new G4Element("Chlorine",  "Cl", 17.,  35.453*g/mole);
  new G4Element("Argon",     "Ar", 18,   39.948*g/mole);
  new G4Element("Pottassium","K",  19.,  39.098*g/mole);
  new G4Element("Calcium",   "Ca", 20.,  40.08*g/mole);
  new G4Element("Scandium",  "Sc", 21.,  44.9559*g/mole);
  new G4Element("Titanium",  "Ti", 22.,  47.90*g/mole);
  new G4Element("Vanadium",  "V",  23.,  50.9414*g/mole);
  new G4Element("Chromium",  "Cr", 24.,  51.996*g/mole);
  new G4Element("Manganese", "Mn", 25.,  54.9380*g/mole);
  new G4Element("Iron",      "Fe", 26.,  55.847*g/mole);
  new G4Element("Cobalt",    "Co", 27.,  58.9332*g/mole);
  new G4Element("Nickel",    "Ni", 28.,  58.70*g/mole);
  new G4Element("Copper",    "Cu", 29.,  63.546*g/mole);
  new G4Element("Zinc",      "Zn", 30.,  65.38*g/mole);
  new G4Element("Gallium",   "Ga", 31.,  69.72*g/mole);
  new G4Element("Germanium", "Ge", 32.,  72.59*g/mole);
  new G4Element("Arsenic",   "As", 33.,  74.9216*g/mole);
  new G4Element("Selenium",  "Se", 34.,  78.96*g/mole);
  new G4Element("Bromine",   "Br", 35.,  79.904*g/mole);
  new G4Element("Krypton",   "Kr", 36.,  83.80*g/mole);
  new G4Element("Rubidium",  "Rb", 37.,  85.4678*g/mole);
  new G4Element("Strontium", "Sr", 38.,  87.62*g/mole);
  new G4Element("Yttrium",   "Y",  39.,  88.9059*g/mole);
  new G4Element("Zirconium", "Zr", 40.,  91.22*g/mole);
  new G4Element("Niobium",   "Nb", 41.,  92.9064*g/mole);
  new G4Element("Molybdenum","Mo", 42.,  95.94*g/mole);
  new G4Element("Technetium","Tc", 43.,  97.*g/mole);
  new G4Element("Ruthenium", "Ru", 44.,  101.07*g/mole);
  new G4Element("Rhodium",   "Rh", 45.,  102.9055*g/mole);
  new G4Element("Palladium", "Pd", 46.,  106.4*g/mole);
  new G4Element("Silver",    "Ag", 47.,  107.868*g/mole);
/*  
  new G4Element("Cadmium",   "", ,  *g/mole);
  new G4Element("Indium",   "", ,  *g/mole);
  new G4Element("Tin",   "", ,  *g/mole);
  new G4Element("Antimony",   "", ,  *g/mole);
  new G4Element("Tellurium",   "", ,  *g/mole);
  new G4Element("Iodine",   "", ,  *g/mole);
  new G4Element("Xenon",   "", ,  *g/mole);
  new G4Element("Cesium",   "", ,  *g/mole);
  new G4Element("Ba??",   "", ,  *g/mole);
  new G4Element("Lanthanum",   "", ,  *g/mole);
  new G4Element("Hafnium",   "", ,  *g/mole);
  new G4Element("Tantalum",   "", ,  *g/mole);
  new G4Element("Tungsten",   "", ,  *g/mole);
  new G4Element("Rhenium",   "", ,  *g/mole);
  new G4Element("Osmium",   "", ,  *g/mole);
  new G4Element("Iridium",   "", ,  *g/mole);
  new G4Element("Platinum",   "", ,  *g/mole);
  new G4Element("Gold",   "", ,  *g/mole);
  new G4Element("Mercury",   "", ,  *g/mole);
  new G4Element("Thallium",   "", ,  *g/mole);
  new G4Element("Lead",      "Pb", 82., 207.19*g/mole);
  new G4Element("Bismuth",   "", ,  *g/mole);
  new G4Element("Polonium",   "", ,  *g/mole);
  new G4Element("Astatine",   "", ,  *g/mole);
  new G4Element("Radon",   "", ,  *g/mole);
  new G4Element("Franicum",   "", ,  *g/mole);
  new G4Element("Radium",   "", ,  *g/mole);
  new G4Element("Actinium",   "", ,  *g/mole);
*/
}

