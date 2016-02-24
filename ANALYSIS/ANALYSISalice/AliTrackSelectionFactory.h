/**
 * \file AliTrackSelectionFactory.h
 * \brief Declaration of AliTrackSelectionFactory
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \since Feb. 24, 2015
 */
#ifndef ALITRACKSELECTIONFACTORY_H_
#define ALITRACKSELECTIONFACTORY_H_
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class AliVTrackSelection;

/**
 * \class AliTrackSelectionFactory
 * \brief Base class for track selection generators
 */
class AliTrackSelectionFactory : public TObject {
public:
  /**
   * Switch for different data types
   */
  enum DataType_t{
    kESD,//!< ESD analysis
    kAOD //!< AOD analysis
  };
  AliTrackSelectionFactory();
  virtual ~AliTrackSelectionFactory();

  /**
   * Method creating track selection, to be implemented by the user
   * @return virtual track selection object
   */
  virtual AliVTrackSelection *CreateTrackCuts(DataType_t datatype) const = 0;

  ClassDef(AliTrackSelectionFactory, 1);
};

#endif /* ALITRACKSELECTIONFACTORY_H_ */
