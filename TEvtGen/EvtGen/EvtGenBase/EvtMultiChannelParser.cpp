//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtMultiChannelParser.cpp,v 1.4 2009-03-16 15:46:01 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998 Caltech, UCSB
//
// Module creator:
//      Alexei Dvoretskii, Caltech, 2001-2002.
//-----------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtMultiChannelParser.hh"
#include "EvtGenBase/EvtDecayMode.hh"
#include "EvtGenBase/EvtPDL.hh"

using std::string;
using std::vector;

EvtDecayMode EvtMultiChannelParser::getDecayMode(const char* file)
{
  // Open file, read tokens

  EvtParser parser;
  parser.read(file);

  // Seek Decay

  int i = 0;
  int N = parser.getNToken();
  while(i<N) {
    
    std::string tok = parser.getToken(i++);
    if(tok == std::string("Decay")) break;
  }
  
  // Get mother

  string mother = string(parser.getToken(i++).c_str());
  std::string bf = parser.getToken(i++);

  vector<string> dauV;
  // Get daughters
    
  while(1) {

    std::string d = parser.getToken(i++);

    if(EvtPDL::getStdHep(EvtPDL::getId(d.c_str())) == 0) break;
    
    dauV.push_back(string(d.c_str()));
  }
  
  EvtDecayMode mode(mother,dauV);
  printf("Decay File defines mode %s\n",mode.mode().c_str());

  return mode;
}



void EvtMultiChannelParser::parse(const char* file, const char* model)
{
  // Open file, read tokens

  EvtParser parser;
  parser.read(file);


  // Get parameters (tokens between the model name and ;)

  int i = 0;
  int N = parser.getNToken();

  // Seek the model name

  while(i<N) {
   
    std::string tok = parser.getToken(i++);
    if(tok == std::string(model)) break;
  }
  if(i == N) {

    printf("No model %s found in decay file %s",model,file);
    exit(0);
  }  

  
  // Add all tokens up to a semicolon to vector 

  std::vector<std::string> v;
  while(i<N) {

    std::string tok = parser.getToken(i++);
    if(tok == std::string(";")) break;
    else v.push_back(tok);
  }
  
  if(i == N) {

    printf("No terminating ; found in decay file %s",file);
    assert(0);
  }

  parse(v);
}


void EvtMultiChannelParser::parse(const std::vector<std::string>& v)
{
  // place holder for strtod
  char** tc = 0;


  // Get PDF maximum or number of points to 
  // use in the scan.

  if(v[0] == std::string("MAXPDF")) {
    
    _pdfMax = strtod(v[1].c_str(),tc);
    if(_pdfMax <= 0) { printf("Bad pdfMax=%f\n",_pdfMax); assert(0); }
  }
  else 
    if(v[0] == std::string("SCANPDF")) {

      _nScan = atoi(v[1].c_str());
    }
    else {

      printf("Error parsing decay file\n");
      assert(0);
    }


  // Now parse the rest of file for amplitude specifications.

  bool conjugate = false;
  size_t i = 2;
  assert(isKeyword(v[2]));
  
  while(i  < v.size()) {
  
    size_t i0 = i;
  
    // Switch to conjugate amplitudes after keyword
    if(v[i] == std::string("CONJUGATE")) {
      
      assert(conjugate == false);
      conjugate = true;
      i++;
      _dm =  strtod(v[i++].c_str(),tc);
      _mixAmpli = strtod(v[i++].c_str(),tc);
      _mixPhase = strtod(v[i++].c_str(),tc);
    }

    if (i >= v.size()) break;
    std::vector<std::string> params;
    EvtComplex c;
    int format;

    if(!conjugate && v[i] == std::string("AMPLITUDE")) {

      while(!isKeyword(v[++i])) params.push_back(v[i]);
      _amp.push_back(params);
      
      parseComplexCoef(i,v,c,format);
      _ampCoef.push_back(c);
      _coefFormat.push_back(format);
      continue;
    }
    else 
      if(conjugate && v[i] == std::string("AMPLITUDE")) {

	while(!isKeyword(v[++i])) params.push_back(v[i]);
	_ampConj.push_back(params);
	parseComplexCoef(i,v,c,format);
	_ampConjCoef.push_back(c);
	_coefConjFormat.push_back(format);
	continue;
      }
      else {
	
	printf("Expect keyword, found parameter %s\n",v[i].c_str());
	assert(0);
      }
    

    assert(i > i0); _unused( i0 );
  }

  printf("PARSING SUCCESSFUL\n");
  printf("%d amplitude terms\n",(int)_amp.size());
  printf("%d conj amplitude terms\n",(int)_ampConj.size());
}



void EvtMultiChannelParser::parseComplexCoef(size_t& i, const std::vector<std::string>& v,
					     EvtComplex& c, int& format) 
{
  // place holder for strtod
  char** tc = 0;

  std::string coefString = v[i++];
  assert(coefString == std::string("COEFFICIENT"));

  if(v[i] == std::string("POLAR_DEG")) {

      double mag = strtod(v[i+1].c_str(),tc);
      double phaseRad = strtod(v[i+2].c_str(),tc)*EvtConst::pi/180.0;
      i += 3;
      c = EvtComplex(mag*cos(phaseRad),mag*sin(phaseRad));  
      format = POLAR_DEG;
  }
  else if(v[i] == std::string("POLAR_RAD")) {
    
    double mag = strtod(v[i+1].c_str(),tc);
    double phaseRad = strtod(v[i+2].c_str(),tc);
    i += 3;
    c = EvtComplex(mag*cos(phaseRad),mag*sin(phaseRad));  
    format = POLAR_RAD;    
  }
  else if(v[i] == std::string("CARTESIAN")) {

    double re = strtod(v[i+1].c_str(),tc);
    double im = strtod(v[i+2].c_str(),tc);
    i += 3;
    c = EvtComplex(re,im);  
    format = CARTESIAN;
  }
  else {
    
    printf("Invalid format %s for complex coefficient\n",v[i].c_str());
    exit(0);
  }
}


double EvtMultiChannelParser::parseRealCoef(int& i, const std::vector<std::string>& v)
{
  // place holder for strtod
  char** tc = 0;
  double value = 0;

  if(v[i] == std::string("COEFFICIENT")) {

    value = strtod(v[i+1].c_str(),tc);
  }
  else assert(0);

  i += 2;

  assert(value > 0.);
  return value;
}


bool EvtMultiChannelParser::isKeyword(const std::string& s)
{
  if(s == std::string("AMPLITUDE")) return true;
  if(s == std::string("CONJUGATE")) return true;
  if(s == std::string("COEFFICIENT")) return true;
  return false;
}



