#ifndef ALIMUONRECLOADER_H
#define ALIMUONRECLOADER_H

/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup rec
/// \class AliMUONRecLoader
/// \brief Implements AliLoader for MUON subsystem
///
/// \author Gines Martinez

#include "AliLoader.h"

class AliMUONRecData;


class AliMUONRecLoader : public AliLoader 
{
  public:
    AliMUONRecLoader();
    AliMUONRecLoader(const Char_t *detname,const Char_t *eventfoldername); //contructor with name of the top folder of the tree
    AliMUONRecLoader(const Char_t *detname,TFolder* eventfolder);
    virtual ~AliMUONRecLoader();

    void              SetMUONData(AliMUONRecData * MUONData);
    AliMUONRecData *  GetMUONData();
 
  protected:
    /// Not implemented
    AliMUONRecLoader(const AliMUONRecLoader& rhs);
    /// Not implemented
    AliMUONRecLoader& operator=(const AliMUONRecLoader& rhs);

    AliMUONRecData * fMUONData; ///< data for MUON subsystem 

  ClassDef(AliMUONRecLoader,1)
};

#endif
