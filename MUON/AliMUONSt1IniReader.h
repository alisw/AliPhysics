#ifndef ALI_MUON_ST1_INI_READER_H
#define ALI_MUON_ST1_INI_READER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <fstream>

#include "AliMUONSt1Types.h"

class AliMUONSt1IniReader
{
  public:
    enum TType {kUndef,kChapter,kValue};

    typedef pair<string,string> TValuePair;
    typedef vector<TValuePair> TValueList;
    typedef pair<string,TValueList> TChapter;
    typedef multimap <string,TValueList> TChapterList;

  public:
    AliMUONSt1IniReader();
    AliMUONSt1IniReader(string fileName);
    virtual ~AliMUONSt1IniReader();
  
    bool   ReadNextLine();
    TType  GetCurrentType()  const  {return fCurrentType; }
    string GetCurrentName()  const  {return fCurrentName; }
    string GetCurrentValue() const  {return fCurrentValue;}
    TChapter     MakeCurrentChapter();
    TValueList   MakeCurrentValueList();
    TChapterList MakeChapterList();
    bool Eof() const {return fEndOfFile;}
    void Reset() ;

  private:
    string Trail(const string& s) const;

    ifstream fFile;        // the file to be read
    TType    fCurrentType; // current type of line (either kChapter or kValue)
    string   fCurrentName; // name of chapter / name of parameter pair
    string   fCurrentValue;// value of the parameter pair if the type is kValue
    bool     fEndOfFile;   // true if the file is entirely read
};

#endif //ALI_MUON_ST1_INI_READER_H 
