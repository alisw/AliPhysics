#ifndef ALI_MUON_ST1_INI_READER_H
#define ALI_MUON_ST1_INI_READER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

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

#ifndef __HP_aCC
  using std::string;
  using std::pair;
  using std::vector;
  using std::multimap;
#endif  

class AliMUONSt1IniReader
{
  public:
    enum IniType {kUndef, kChapter, kValue};

    typedef pair<string, string> ValuePair;
    typedef vector<ValuePair>  ValueList;
    typedef pair<string, ValueList> Chapter;
    typedef multimap <string, ValueList> ChapterList;

  public:
    AliMUONSt1IniReader();
    AliMUONSt1IniReader(string fileName);
    virtual ~AliMUONSt1IniReader();
  
    bool     ReadNextLine();
    IniType  GetCurrentType()  const  {return fCurrentType; }
    string   GetCurrentName()  const  {return fCurrentName; }
    string   GetCurrentValue() const  {return fCurrentValue;}
    Chapter     MakeCurrentChapter();
    ValueList   MakeCurrentValueList();
    ChapterList MakeChapterList();
    bool Eof() const {return fEndOfFile;}
    void Reset() ;

  private:
    string Trail(const string& s) const;

    ifstream fFile;        // the file to be read
    IniType  fCurrentType; // current type of line (either kChapter or kValue)
    string   fCurrentName; // name of chapter / name of parameter pair
    string   fCurrentValue;// value of the parameter pair if the type is kValue
    bool     fEndOfFile;   // true if the file is entirely read
};

#endif //ALI_MUON_ST1_INI_READER_H 
