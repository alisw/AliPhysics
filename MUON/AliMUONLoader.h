#ifndef ALIMUONLOADER_H
#define ALIMUONLOADER_H

/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONLoader
/// \brief Implements AliLoader for MUON subsystem
///
/// \author Gines Martinez

#include "AliLoader.h"

class AliMUONData;


class AliMUONLoader : public AliLoader 
{
  public:
    AliMUONLoader();
    AliMUONLoader(const Char_t *detname,const Char_t *eventfoldername); //contructor with name of the top folder of the tree
    AliMUONLoader(const Char_t *detname,TFolder* eventfolder);
    virtual ~AliMUONLoader();

    void           SetMUONData(AliMUONData * MUONData);
    AliMUONData *  GetMUONData();
 
  protected:
    AliMUONLoader(const AliMUONLoader& rhs);
    AliMUONLoader& operator=(const AliMUONLoader& rhs);

    AliMUONData * fMUONData; ///< data for MUON subsystem 

  private:
    //descendant classes should
    //use protected interface methods to access these folders

    /**********************************************/
    /***********     P U B L I C     **************/
    /*********       S T A T I C       ************/
    /*********         METHODS         ************/
    /*********     They are used by    ************/
    /*********** AliRunLoader as well**************/
    /**********************************************/

  ClassDef(AliMUONLoader,1)
};

#endif
