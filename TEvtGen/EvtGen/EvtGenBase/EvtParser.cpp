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
// Module: EvtParser.cc
//
// Description: Reading the decay table and produce a list of tokens.
//
// Modification history:
//
//    RYD     Febuary 11, 1998        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include <fstream>
#include <sstream>
#include <string.h>
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtReport.hh"
using namespace std;

#define MAXBUF 1024

EvtParser::EvtParser(){
  _ntoken=0;
  _lengthoftokenlist=0;
  _tokenlist=0;
  _linelist=0;
}

EvtParser::~EvtParser(){

  delete [] _tokenlist;
  delete [] _linelist;

}


int EvtParser::getNToken(){

  return _ntoken;

}

const std::string& EvtParser::getToken(int i){

  return _tokenlist[i];

}

int EvtParser::getLineofToken(int i){

  return _linelist[i];

}

int EvtParser::read(const std::string filename){
  ifstream fin;
  
  fin.open(filename.c_str());
  if (!fin) {
    report(Severity::Error,"EvtGen") << "Could not open file '"<<filename.c_str()<<"'"<<endl;
    return -1;
  }

  char buf[MAXBUF];
  char buf2[MAXBUF];
  char c;

  int line=0;
  int i;

  while(fin.peek() != EOF){
    line++;
    
    i=0;
    while((c=fin.get()) != '\n' && i<MAXBUF) {
      buf[i]=c;
      i++;
    }
    if(i==MAXBUF) {
      report(Severity::Error,"EvtGen") << "Error in EvtParser: line:"
			     <<line<<" to long"<<endl;
    }
    else {
      buf[i] = '\0';
    }
    
    //search for '#' which indicates comment for rest of line!
    i=0;
    do{
      if (buf[i]=='#') buf[i]=0;
      i++;
    }while(buf[i-1]!=0);

    string tmp(buf,strlen(buf));

    //read each token
    istringstream ist(tmp);
    while(ist>>buf2){
      i=0;
      int semicolon=0;
      do{
	if (buf2[i]==';') {
	  buf2[i]=0;
	  semicolon=1;
	}
      }while(buf2[i++]!=0);
      if (buf2[0]!=0){
	addToken(line,buf2);
      }
      if (semicolon) addToken(line,";");
    }
  }

  fin.close();

  return 0;
  
}



void EvtParser::addToken(int line,const std::string& string){

  //report(Severity::Info,"EvtGen") <<_ntoken<<" "<<line<<" "<<string<<endl;  

  if (_ntoken==_lengthoftokenlist) {

    int new_length=1000+4*_lengthoftokenlist;


    
    int*     newlinelist= new int[new_length];
    std::string* newtokenlist= new std::string[new_length];
  
    int i;

    for(i=0;i<_ntoken;i++){
     newlinelist[i]=_linelist[i];
     newtokenlist[i]=_tokenlist[i];
    }

    delete [] _tokenlist;
    delete [] _linelist;

    _tokenlist=newtokenlist;
    _linelist=newlinelist;    

    _lengthoftokenlist=new_length;

  }


  _tokenlist[_ntoken]=string;

  _linelist[_ntoken]=line;
 
  _ntoken++;  

  //report(Severity::Info,"EvtGen") << "First:"<<_tokenlist[0]<<" last:"<<_tokenlist[_ntoken-1]<<endl;

}
   
