#ifndef ALIMUONSIMLOADER_H
#define ALIMUONSIMLOADER_H

/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup sim
/// \class AliMUONSimLoader
/// \brief Implements AliLoader for MUON subsystem
///
/// \author Gines Martinez

#include "AliLoader.h"

class AliMUONSimData;


class AliMUONSimLoader : public AliLoader 
{
  public:
    AliMUONSimLoader();
    AliMUONSimLoader(const Char_t *detname,const Char_t *eventfoldername); //contructor with name of the top folder of the tree
    AliMUONSimLoader(const Char_t *detname,TFolder* eventfolder);
    virtual ~AliMUONSimLoader();

    void              SetMUONData(AliMUONSimData * MUONData);
    AliMUONSimData *  GetMUONData();
 
  protected:
    /// Not implemented
    AliMUONSimLoader(const AliMUONSimLoader& rhs);
    /// Not implemented
    AliMUONSimLoader& operator=(const AliMUONSimLoader& rhs);

    AliMUONSimData * fMUONData; ///< data for MUON subsystem 

  ClassDef(AliMUONSimLoader,1)
};

#endif
