//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtParserXml.cc
//
// Description: Reading the decay XML file.
//
// Modification history:
//
//    DCC     24 October, 2011        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include "EvtGenBase/EvtParserXml.hh"
#include "EvtGenBase/EvtReport.hh"
using namespace std;

EvtParserXml::EvtParserXml(){
  _line = "";
  _lineNo=0;
  _tag = "";
  _tagTitle = "";
}

EvtParserXml::~EvtParserXml(){


}


bool EvtParserXml::open(std::string filename){
  
  if(!expandEnvVars(filename)) {
    report(Severity::Error,"EvtGen") << "Error while expanding environment variables in file name '"<<filename.c_str()<<"'"<<endl;
    return false;
  }

  _fin.open(filename.c_str());
  if (!_fin) {
    report(Severity::Error,"EvtGen") << "Could not open file '"<<filename.c_str()<<"'"<<endl;
    return false;
  }

  return true;
  
}

bool EvtParserXml::close() {
  _fin.close();
  return true;
}

bool EvtParserXml::readNextTag() {
    if(!processTagTree()) {
      report(Severity::Error,"EvtGen")
                  << "Unexpected end tag "<<_tagTitle<<" found near line "<<_lineNo<<endl;
      report(Severity::Error,"EvtGen")
                  << "Will terminate execution!"<<endl;
      return false;
    }//first process the previous tag to find out where we are in the tag tree

    while(_line.find("<") == std::string::npos) {//add lines until we find start of a tag
      std::string addLine;
      if(!std::getline(_fin, addLine)) return false;
      _lineNo++;
      _line += " ";
      _line += addLine;
    }

    unsigned int startTag;
    unsigned int endTag;
    unsigned int endTagTitle;

    startTag = _line.find("<");

    if(_line[startTag+1] =='?') { //XML header tag - ignore then read the next tag
      while(_line.find("?>", startTag) == std::string::npos) {
        std::string addLine;
        if(!std::getline(_fin, addLine)) return false;
        _lineNo++;
        _line += " ";
        _line += addLine;
      }
      endTag = _line.find("?>", startTag);
      _line = _line.substr(endTag + 2);
      return readNextTag();
    } else if(_line[startTag+1] == '!') { //XML comment tag - ignore then read the next tag
      while(_line.find("-->", startTag) == std::string::npos) {
        std::string addLine;
        if(!std::getline(_fin, addLine)) return false;
        _lineNo++;
        _line += " ";
        _line += addLine;
      }
      endTag = _line.find("-->", startTag);
      _line = _line.substr(endTag + 3);
      _tagTitle = "";
      _tag = "";
      return readNextTag();
    } else { //parsable

      while(_line.find(">", startTag) == std::string::npos) {//find end of a tag
        std::string addLine;
        if(!std::getline(_fin, addLine)) return false;
        _lineNo++;
        _line += " ";
        _line += addLine;
      }
      endTag = _line.find(">", startTag);
      _inLineTag = false;
      if(_line.find("/>", startTag) < endTag) {
        endTag--;
        _inLineTag = true;
      }

      if(_line.find(" ", startTag) != std::string::npos && _line.find(" ", startTag) < endTag) {//find end of the first word in the tag
        endTagTitle = _line.find(" ", startTag);
      } else {
        endTagTitle = endTag;
      }

      _tagTitle = _line.substr(startTag + 1, endTagTitle - startTag - 1);
      _tag = _line.substr(startTag + 1, endTag - startTag - 1);

      //now we have the tag lets remove it from the line
      if(_inLineTag) {
        _line = _line.substr(endTag+2);
      } else {
        _line = _line.substr(endTag+1);
      }
      return true;
    }
}

std::string EvtParserXml::getParentTagTitle() {
  if(_tagTree.empty()) return "";
  else return _tagTree.back();
}

std::string EvtParserXml::readAttribute(std::string attribute, std::string defaultValue) {
  std::string whitespace = " \t\n\v\f\r";
  for(unsigned int i=0; i<whitespace.size(); i++) {
  //find any whitespace followed by the attribute name followed by an '='
    std::string attName = whitespace[i] + attribute + "=";
    if (_tag.find(attName) != std::string::npos) {
      int startAttri = _tag.find(attName);
      int startQuote = _tag.find("\"", startAttri + 1);
      int endQuote = _tag.find("\"", startQuote + 1);
      return _tag.substr(startQuote + 1, endQuote - startQuote - 1);
    }
  }
  return defaultValue;
}

bool EvtParserXml::readAttributeBool(std::string attribute, bool defaultValue) {
  std::string valStr = readAttribute(attribute);
  if(!defaultValue) return (valStr == "true" || valStr == "1" || valStr == "on" || valStr == "yes");
  else return (valStr != "false" && valStr != "0" && valStr != "off" && valStr != "no");
}

int EvtParserXml::readAttributeInt(std::string attribute, int defaultValue) {
  std::string valStr = readAttribute(attribute);
  if (valStr == "") return defaultValue;
  std::istringstream valStream(valStr);
  int retVal;
  valStream >> retVal;
  return retVal;
}

double EvtParserXml::readAttributeDouble(std::string attribute, double defaultValue) {
  std::string valStr = readAttribute(attribute);
  if (valStr == "") return defaultValue;
  std::istringstream valStream(valStr);
  double retVal;
  valStream >> retVal;
  return retVal;
}

bool EvtParserXml::processTagTree() {
  if(_tagTitle == "") return true;
  if(_tagTitle[0] == '/') {
    if(_tagTitle.substr(1) == _tagTree.back()) {
      _tagTree.pop_back();
    } else {
      return false;
    }
  } else if(!_inLineTag) {
    _tagTree.push_back(_tagTitle);
  }
  return true;
}

bool EvtParserXml::expandEnvVars(std::string& str) {
  while(str.find('$') != std::string::npos) {
    size_t varStart = str.find('$');
    size_t varNameLength;
    std::string varName;
    
    //if this is the last character then just remove the $
    if(varStart == str.length()-1) {
      str.erase(varStart);
      return true;
    }
    
    if(str[varStart+1] == '{') {
      //deal with environment variables in {}s
      size_t braceStart = varStart+1;
      size_t braceEnd = str.find('}',braceStart);
      
      if(braceEnd == std::string::npos) {
        report(Severity::Error,"EvtGen")
          << "Incomplete environment variable found in text: "<<str<<endl;
        report(Severity::Error,"EvtGen")
          << "Will terminate execution!"<<endl;
          return false;
      }

      varName = str.substr(braceStart+1,braceEnd-braceStart-1);
      varNameLength = braceEnd-varStart;

    } else {
      //deal with everything else
      varNameLength=0;

      while(varNameLength + varStart + 1 < str.length() && isAlphaNum(str[varStart+varNameLength+1])) {
        ++varNameLength;
      }

      varName = str.substr(varStart+1,varNameLength);
    }

    char* envVar = getenv(varName.c_str());

    if(envVar) str.replace(varStart,varNameLength+1,envVar);
    else {
      report(Severity::Warning,"EvtGen")
        << "Undefined environment variable found in text: "<<varName<<endl;
      str.replace(varStart,varNameLength+1,"");
    }
  }
  return true;
}

bool EvtParserXml::isAlphaNum(char c) {
  if(c>='0' && c<='9') return true;
  if(c>='A' && c<='Z') return true;
  if(c>='a' && c<='z') return true;
  if(c=='_') return true;
  return false;
}
