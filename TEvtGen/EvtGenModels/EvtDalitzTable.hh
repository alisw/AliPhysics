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
// Module: EvtGen/EvtGenericDalitz.hh
//
// Description: Model to describe a generic dalitz decay
//
// Modification history:
//
//    DCC     16 December, 2011         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDALITZTABLE_HPP
#define EVTDALITZTABLE_HPP

#include "EvtGenModels/EvtDalitzDecayInfo.hh"
#include "EvtGenBase/EvtId.hh"

#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <map>
#include <string>
#include <vector>

class EvtDalitzTable {
public:

  static EvtDalitzTable* getInstance(const std::string dec_name="", bool verbose=true);

  bool fileHasBeenRead(const std::string dec_name);
  void readXMLDecayFile(const std::string dec_name, bool verbose=true);
  void checkParticle(std::string particle);

  void addDecay(EvtId parent, const EvtDalitzDecayInfo& dec);
  void copyDecay(EvtId parent, EvtId* daughters, EvtId copy, EvtId* copyd);

  std::vector<EvtDalitzDecayInfo> getDalitzTable(const EvtId& parent);

protected:

  EvtDalitzTable();
  ~EvtDalitzTable();

private:

  EvtDalitzReso getResonance(std::string shape, EvtDalitzPlot dp, EvtCyclic3::Pair angPair, EvtCyclic3::Pair resPair,
                             EvtSpinType::spintype spinType, double mass, double width, double FFp, double FFr, double alpha,
                             double aLass, double rLass, double BLass, double phiBLass, double RLass, double phiRLass, double cutoffLass);
  int getDaughterPairs(EvtId* resDaughter, EvtId* daughter, std::vector< std::pair<EvtCyclic3::Pair,EvtCyclic3::Pair> >& angAndResPairs);

  std::map<EvtId, std::vector<EvtDalitzDecayInfo> > _dalitztable;
  std::vector<std::string> _readFiles;

  EvtDalitzTable(const EvtDalitzTable&);
  EvtDalitzTable& operator=(const EvtDalitzTable&);

  //to calculate probMax
  double calcProbMax(EvtDalitzPlot dp, EvtDalitzDecayInfo* model);
  double calcProb(EvtDalitzPoint point, EvtDalitzDecayInfo* model);
};

#endif
