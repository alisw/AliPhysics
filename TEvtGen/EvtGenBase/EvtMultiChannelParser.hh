//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtMultiChannelParser.hh,v 1.2 2009-03-16 16:41:09 robbep Exp $
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

// Parse decay file tokens from a file or a vector

#ifndef EVT_MULTI_CHANNEL_PARSER
#define EVT_MULTI_CHANNEL_PARSER


#include <vector>
#include <string>
#include "EvtGenBase/EvtComplex.hh"

class EvtDecayMode;

enum {POLAR_RAD,POLAR_DEG,CARTESIAN};  

class EvtMultiChannelParser {

public:

  EvtMultiChannelParser() : _pdfMax(-1.), _nScan(0), _dm(0.) {}
  ~EvtMultiChannelParser() {}

  static EvtDecayMode getDecayMode(const char* file);

  void parse(const char* file, const char* model);
  void parse(const std::vector<std::string>& v);

  static void parseComplexCoef(size_t& i, const std::vector<std::string>& v,
			       EvtComplex& c, int& format);
  static double parseRealCoef(int& i, const std::vector<std::string>& v);
  static bool isKeyword(const std::string& s);

  inline double pdfMax() const { return _pdfMax; }
  inline int nScan() const { return _nScan; }
  inline double dm() const { return _dm; }
  inline double mixPhase() const { return _mixPhase; }
  inline double mixAmpli() const { return _mixAmpli; }

  inline std::vector<std::string> amp(int i) const { return _amp[i]; }
  inline std::vector<std::string> ampConj(int i) const { return _ampConj[i]; }
  inline EvtComplex ampCoef(int i) const { return _ampCoef[i]; }
  inline EvtComplex ampConjCoef(int i) const { return _ampConjCoef[i]; }

  inline int coefFormat(int i) const { return _coefFormat[i]; }
  inline int coefConjFormat(int i) const { return _coefConjFormat[i]; }

  inline int getNAmp() const { return _amp.size(); } 
  inline int getNAmpConj() const { return _ampConj.size(); } 

private:

  double _pdfMax;
  int _nScan;
  double _dm;
  double _mixPhase;
  double _mixAmpli;

  std::vector<std::vector<std::string> >  _amp;
  std::vector<std::vector<std::string> >  _ampConj;
  std::vector<EvtComplex> _ampCoef;
  std::vector<int> _coefFormat;
  std::vector<EvtComplex> _ampConjCoef;
  std::vector<int> _coefConjFormat;
};

#endif





