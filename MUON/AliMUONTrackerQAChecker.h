#ifndef ALIMUONTRACKERQACHECKER_H
#define ALIMUONTRACKERQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec 
/// \class AliMUONTrackerQAChecker
/// \brief Implementation of QAChecker for MCH
///
//  Author: Laurent Aphecetche

#include "AliMUONVQAChecker.h"

class TH1;

class AliMUONTrackerQAChecker: public AliMUONVQAChecker {

public:
  AliMUONTrackerQAChecker();
  AliMUONTrackerQAChecker(const AliMUONTrackerQAChecker& qac);
  virtual ~AliMUONTrackerQAChecker();

  virtual AliMUONVQAChecker::ECheckCode* CheckRaws(TObjArray** list, const AliMUONRecoParam* recoParam);
  virtual AliMUONVQAChecker::ECheckCode* CheckRecPoints(TObjArray** list, const AliMUONRecoParam* recoParam);
  virtual AliMUONVQAChecker::ECheckCode* CheckESD(TObjArray** list, const AliMUONRecoParam* recoParam);

private:
  
  AliMUONVQAChecker::ECheckCode MarkHisto(TH1& histo, AliMUONVQAChecker::ECheckCode value) const;
  
  AliMUONVQAChecker::ECheckCode BeautifyOccupancyHistograms(TH1& hddl,
                                                            TH1& hbp, 
                                                            const TH1* hbuspatchconfig, 
                                                            Int_t neventsseen,
                                                            Int_t neventsused,
                                                            const AliMUONRecoParam& recoParam);

  AliMUONVQAChecker::ECheckCode BeautifyReadoutHistograms(TH1& hroe,
                                                          TH1& hroenorm,
                                                          const TH1& hbuspatchtokenerrors,
                                                          Int_t neventsseen,
                                                          Int_t neventsused,
                                                          const AliMUONRecoParam& recoParam);

  AliMUONVQAChecker::ECheckCode BeautifyEventsizeHistograms(TH1& heventsize,
                                                            TH1& heventsizeperevent,
                                                            Int_t neventsseen,
                                                            Int_t neventsused,
                                                            const AliMUONRecoParam& recoParam);
  
  ClassDef(AliMUONTrackerQAChecker,1)  // MUON quality assurance checker

};

#endif 
