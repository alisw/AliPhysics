/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONIniReader
// ----------------------
// General class to read data in ASCII file format,
// similar to the Windows ".ini" files (a set of sections tagged by a 
//                      [ sectionName ]
//  and values defined in the way:
//                    parameterName = value
//
// comment lines can be introduced if the first non-blank character
// is either ';' or '#'
// Included in AliRoot 2003/01/28

#if !defined(__HP_aCC) && !defined(__alpha)
  #include <sstream>
#endif

#include <Riostream.h>
#include <Rstrstream.h>

#include "AliMUONSt1IniReader.h"

//______________________________________________________________________
AliMUONSt1IniReader::AliMUONSt1IniReader()
  :fFile(),fCurrentType(kUndef),fEndOfFile(true)
{
// default constructor
// ---
}

//______________________________________________________________________
AliMUONSt1IniReader::AliMUONSt1IniReader(string fileName)
{
// normal constructor
// ---

  fFile.open(fileName.c_str());
  if (!fFile) {cerr<<"Unable to open file "<<fileName<<endl;}
  fEndOfFile = !fFile.good();
  fCurrentType=kUndef;
}

//______________________________________________________________________
AliMUONSt1IniReader::~AliMUONSt1IniReader()
{
  //destructor
 fFile.close();
}

//______________________________________________________________________
void AliMUONSt1IniReader::Reset()
{
// Reset the input stream. The file can be re-read after calling this function
// ---

  fFile.clear();
  fFile.seekg(0,ios::beg);
  fCurrentType=kUndef;
  fEndOfFile=!fFile.good();
}

//______________________________________________________________________
bool AliMUONSt1IniReader::ReadNextLine()
{
// The main function of this class.
// Read next line in the file and set CurrentType(), CurrentName() and
// CurrentValue() with the line's content
// ---

  if ( (!fFile) || (fFile.eof()) || (!fFile.good()) ) 
    {fEndOfFile=true; fCurrentType=kUndef; return false;}

  string line;
  getline(fFile,line);
  if ( line.empty()) {             // this is a blank line
    return ReadNextLine();
  }
#if defined (__HP_aCC) || (__alpha)
  strstream l;
  l << line;
#else
  istringstream l(line); 
#endif    
  
  char c;

  l>>c;
  if ( (c==';') || (c=='#') ) {    // this is a comment
    return ReadNextLine();
  }
  
  if (c=='[') {       	     // this is a chapter name
    getline(l,fCurrentName,']');
    fCurrentName=Trail(fCurrentName);
    fCurrentType=kChapter;
    return true;
  } else {
    if (line.find_first_of("=") != string::npos ) {
      l.putback(c);
      getline(l,fCurrentName,'=');
      fCurrentName = Trail(fCurrentName);
       
      getline(l,fCurrentValue);
      fCurrentValue = Trail(fCurrentValue);
      fCurrentType=kValue;
      return true;
    } else {
      cerr<<"Warning, badly formated line..."<<line<<endl;
      fCurrentType=kUndef;
      return false;
    }
  }
  // fCurrentType=kUndef;
  // return false;
       // unreachable
}

//______________________________________________________________________
AliMUONSt1IniReader::ValueList AliMUONSt1IniReader::MakeCurrentValueList()
{
// Read the next lines in the file
// until eof() or a new section is found. 
// Return the list of (name,value) pairs read.
// ---

  ValueList ans;
  while (true){
    if (fCurrentType==kValue){
      ans.push_back( ValuePair(fCurrentName,fCurrentValue));
    } else break;
    ReadNextLine();
  }
  return ans;
}

//______________________________________________________________________
AliMUONSt1IniReader::Chapter AliMUONSt1IniReader::MakeCurrentChapter()
{
// Searches in the rest file for a new section
// and return it's name and the list of (name,value) pairs in it
// ---

  while ((!Eof()) && (fCurrentType != kChapter)) ReadNextLine();
  if (Eof()) return Chapter();
  string name = fCurrentName;
  ReadNextLine();
  return Chapter(name,MakeCurrentValueList());
}

//______________________________________________________________________
AliMUONSt1IniReader::ChapterList AliMUONSt1IniReader::MakeChapterList()
{
// Read the rest of the file and return all the chapter names and
// (name,value) pair lists found after the current position
// ---

  ChapterList ans;
  while (true) {
    if (fCurrentType==kChapter) {
      string s= fCurrentName;
      ReadNextLine();
      //ans.insert(Chapter(s,MakeCurrentValueList()));
                   // does not compile on SunOS
      ans.insert(ChapterList::value_type(s,MakeCurrentValueList()));
    } else ReadNextLine();
    if (fEndOfFile) break;
  }
  return ans;
}

//______________________________________________________________________
string AliMUONSt1IniReader::Trail(const string& s) const
{
// Utility function: clear the blanks before and after the string <s>
// ---

  string::size_type p1=s.find_first_not_of(" ");
  if (p1==string::npos) return "";
  string::size_type p2=s.find_last_not_of(" ");
  return s.substr(p1,p2-p1+1);
}
