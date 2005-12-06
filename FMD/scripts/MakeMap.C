//
// $Id$
//
// Script to make a class derived from AliFMDMap. 
//
#ifndef __CINT__
#include <fstream>
#include <TString.h>
#include <TDatime.h>
#include <TSystem.h>
#include <iostream>
using namespace std;
#endif

void MakeMap(const Char_t* type="Int_t", const Char_t* name=0) 
{
  TString base;
  TString ttype(type);
  if (ttype.EndsWith("_t")) {
    Ssiz_t undert = ttype.Index("_t");
    ttype.Remove(undert);
  }
  if (!name) 
    base = Form("AliFMD%sMap", ttype.Data());
  else
    base = name;

  cout << "Base name is " << base << endl;
  
  TString decl_name(Form("%s.h", base.Data()));
  TString impl_name(Form("%s.cxx", base.Data()));
  ofstream decl(decl_name.Data());
  ofstream impl(impl_name.Data());

  if (!decl) {
    cerr << "Cannot open declaration file " << decl_name << endl;
    return;
  }
  if (!impl) {
    cerr << "Cannot open implementation file " << impl_name << endl;
    return;
  }
  
  TDatime now;
  cout << "The time is now " << now.AsString() << endl;
  UserGroup_t* uinfo = gSystem->GetUserInfo();
  TString     uname(uinfo->fRealName);
  Ssiz_t      comma = uname.Index(",");
  if (comma != kNPOS) uname.Remove(comma);
  cout << "User's name is " << uname << endl;
  TString guard(base);
  guard.Append("_h");
  guard.ToUpper();

  cout << "Writing declaration file " << decl_name << " ... " 
       << flush;
  decl << "#ifndef " << guard << "\n";
  decl << "#define " << guard << "\n";
  decl << "/* Copyright (c) " << now.GetYear() ;
  decl << ", ALICE Experiment @ CERN.\n" ;
  decl << " * All rights reserved\n";
  decl << " * See " << impl_name << " for full copyright notice\n";
  decl << " * \n" ;
  decl << " * Created " << now.AsString() << " by " << uname << "\n";
  decl << " */\n";
  decl << "/* $Id$ */\n";
  decl << "//__________________________________________________________\n";
  decl << "// \n";
  decl << "// Map of " << type << " for each FMD strip\n" ;
  decl << "// \n";
  decl << "#ifndef ALIFMDMAP_H\n";
  decl << "# include <AliFMDMap.h>\n";
  decl << "#endif\n\n";
  decl << "class " << base << " : public AliFMDMap\n";
  decl << "{\n";
  decl << "public:\n";
  decl << "  " << base << "(const " << base << "& other);\n";
  decl << "  " << base << "(size_t maxDet  = kMaxDetectors,\n";
  decl << "                 size_t maxRing = kMaxRings,\n";
  decl << "                 size_t maxSec  = kMaxSectors,\n";
  decl << "                 size_t maxStr  = kMaxStrips);\n";
  decl << "  virtual ~" << base << "() { delete [] fData; }\n";
  decl << "  " << base << "& operator=(const " << base << "& other);\n";
  decl << "  virtual void Clear(const " << type << "& v=" << type << "());\n";
  decl << "  virtual " << type << "& operator()(UShort_t det,\n";
  decl << "                                     Char_t   ring,\n";
  decl << "                                     UShort_t sec,\n";
  decl << "                                     UShort_t str);\n";
  decl << "  virtual const " << type << "& operator()(UShort_t det,\n";
  decl << "                                           Char_t   ring,\n";
  decl << "                                           UShort_t sec,\n";
  decl << "                                           UShort_t str) const;\n";
  decl << "protected:\n";
  decl << "  " << type << "* fData; // The Data\n";
  decl << "  ClassDef(" << base << ",1) // Map of " << type ;
  decl << " data per strip\n" ;
  decl << "};\n\n";
  decl << "#endif\n";
  decl << "//__________________________________________________________\n";
  decl << "// \n";
  decl << "// Local Variables:\n";
  decl << "//   mode: C++\n";
  decl << "// End:\n";
  decl << "//" << endl;;
  decl.close();
  cout << "done" << endl;

  cout << "Writing implementation file " << impl_name << " ... " 
       << flush;
  impl << "/**************************************************************\n";
  impl << " * Copyright(c) 1998-1999, ALICE Experiment at CERN.          *\n";
  impl << " * All rights reserved.                                       *\n";
  impl << " *                                                            *\n";
  impl << " * Author: The ALICE Off-line Project.                        *\n";
  impl << " * Contributors are mentioned in the code where appropriate.  *\n";
  impl << " *                                                            *\n";
  impl << " * Permission to use, copy, modify and distribute this        *\n";
  impl << " * software and its documentation strictly for non-commercial *\n";
  impl << " * purposes is hereby granted without fee, provided that the  *\n";
  impl << " * above copyright notice appears in all copies and that both *\n";
  impl << " * the copyright notice and this permission notice appear in  *\n";
  impl << " * the supporting documentation. The authors make no claims   *\n";
  impl << " * about the suitability of this software for any purpose. It *\n";
  impl << " * is provided \"as is\" without express or implied warranty.   *\n";
  impl << " **************************************************************/\n";
  impl << "/* $Id$ */\n";
  impl << "//__________________________________________________________\n";
  impl << "// \n";
  impl << "// Map of per strip " << type << " information\n";
  impl << "// \n";
  impl << "// Created " << now.AsString() << " by " << uname << "\n";
  impl << "// \n";
  impl << "#include \"" << decl_name << "\"\t//" << guard << "\n";
  impl << "//__________________________________________________________\n";
  impl << "ClassImp(" << base << ");\n";
  impl << "//__________________________________________________________\n";
  impl << base << "::" << base << "(const " << base << "& other)\n";
  impl << "  : AliFMDMap(other.fMaxDetectors,\n";
  impl << "              other.fMaxRings,\n";
  impl << "              other.fMaxSectors,\n";
  impl << "              other.fMaxStrips),\n";
  impl << "    fData(0)\n";
  impl << "{\n";
  impl << "  // Copy constructor\n";
  impl << "  fData = new " << type << "[fMaxDetectors * fMaxRings ";
  impl << "* fMaxSectors * fMaxStrips];\n" ;
  impl << "  for (size_t i = 0; i < fMaxDetectors * fMaxRings ";
  impl << "* fMaxSectors * fMaxStrips; i++)\n";
  impl << "    fData[i] = other.fData[i];\n";
  impl << "}\n\n";
  impl << "//__________________________________________________________\n";
  impl << base << "::" << base << "(size_t maxDet,\n";
  impl << "                         size_t maxRing,\n";
  impl << "                         size_t maxSec,\n";
  impl << "                         size_t maxStr)\n";
  impl << "  : AliFMDMap(maxDet, maxRing, maxSec, maxStr),\n";
  impl << "    fData(0)\n";
  impl << "{\n";
  impl << "  // Constructor.\n";
  impl << "  // Parameters:\n";
  impl << "  //\tmaxDet\tMaximum number of detectors\n";
  impl << "  //\tmaxRing\tMaximum number of rings per detector\n";
  impl << "  //\tmaxSec\tMaximum number of sectors per ring\n";
  impl << "  //\tmaxStr\tMaximum number of strips per sector\n";
  impl << "  fData = new " << type << "[fMaxDetectors * fMaxRings ";
  impl << "* fMaxSectors * fMaxStrips];\n" ;
  impl << "  Clear();\n";
  impl << "}\n\n";
  impl << "//__________________________________________________________\n";
  impl << base << "&\n";
  impl << base << "::operator=(const " << base << "& other)\n";
  impl << "{\n";
  impl << "  // Assignment operator \n";
  impl << "  fMaxDetectors = other.fMaxDetectors;\n";
  impl << "  fMaxRings     = other.fMaxRings;\n";
  impl << "  fMaxSectors   = other.fMaxSectors;\n";
  impl << "  fMaxStrips    = other.fMaxStrips;\n";
  impl << "  if (fData) delete [] fData;\n";
  impl << "  fData = new " << type << "[fMaxDetectors * fMaxRings ";
  impl << "* fMaxSectors * fMaxStrips];\n" ;
  impl << "  for (size_t i = 0; i < fMaxDetectors * fMaxRings ";
  impl << "* fMaxSectors * fMaxStrips; i++)\n";
  impl << "    fData[i] = other.fData[i];\n";
  impl << "}\n\n"    ;
  impl << "//__________________________________________________________\n";
  impl << "void\n";
  impl << base << "::Clear(const " << type << "& val)\n";
  impl << "{\n";
  impl << "  // Reset map to val\n";
  impl << "  for (size_t i = 0; i < fMaxDetectors * fMaxRings ";
  impl << "* fMaxSectors * fMaxStrips; i++)\n";
  impl << "    fData[i] = val;\n";
  impl << "}\n\n"    ;
  impl << "//__________________________________________________________\n";
  impl << type << "&\n";
  impl << base << "::operator()(UShort_t det, Char_t ring, UShort_t sec, " ;
  impl << "UShort_t str)\n" ;
  impl << "{\n" ;
  impl << "  // Get data\n";
  impl << "  // Parameters:\n";
  impl << "  //\tdet\tDetector #\n";
  impl << "  //\tring\tRing ID\n";
  impl << "  //\tsec\tSector #\n";
  impl << "  //\tstr\tStrip #\n" ;
  impl << "  // Returns appropriate data\n";
  impl << "  return fData[CalcIndex(det, ring, sec, str)];\n";
  impl << "}\n\n";
  impl << "//__________________________________________________________\n";
  impl << "const " << type << "&\n";
  impl << base << "::operator()(UShort_t det, Char_t ring, UShort_t sec, " ;
  impl << "UShort_t str) const\n" ;
  impl << "{\n" ;
  impl << "  // Get data\n";
  impl << "  // Parameters:\n";
  impl << "  //\tdet\tDetector #\n";
  impl << "  //\tring\tRing ID\n";
  impl << "  //\tsec\tSector #\n";
  impl << "  //\tstr\tStrip #\n" ;
  impl << "  // Returns appropriate data\n";
  impl << "  return fData[CalcIndex(det, ring, sec, str)];\n";
  impl << "}\n\n";
  impl << "//__________________________________________________________\n";
  impl << "// \n";
  impl << "// EOF\n";
  impl << "// \n";
  impl << endl;
  impl.close();
  cout << "done" << endl;
}

  


#ifndef __CINT__
int main() 
{
  makemap();
  return 0;
}
#endif

//____________________________________________________________________
//
// EOF
//
