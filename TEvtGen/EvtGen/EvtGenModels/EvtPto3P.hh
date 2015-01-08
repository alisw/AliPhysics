//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtPto3P.hh,v 1.2 2009-03-16 16:31:05 robbep Exp $
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

#ifndef EVT_PTO3P_HH
#define EVT_PTO3P_HH

#include <vector>
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenModels/EvtIntervalDecayAmp.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"


class EvtPto3P : public  EvtIntervalDecayAmp<EvtDalitzPoint> {
  
public:
  
  EvtPto3P() {}
  virtual ~EvtPto3P() {}
  std::string getName() { return "PTO3P"; }
  EvtDecayBase* clone() { return new EvtPto3P(); } 
  
  
  virtual EvtAmpFactory<EvtDalitzPoint>* createFactory(const EvtMultiChannelParser& parser);
  virtual std::vector<EvtVector4R> initDaughters(const EvtDalitzPoint& p) const;
  
  EvtDalitzPlot dp();
  
};

#endif

