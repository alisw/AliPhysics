#ifndef ALIMUONLOADER_H
#define ALIMUONLOADER_H

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include "AliLoader.h"

//__________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliMUONLoader                                            //
//                                                                 //
/////////////////////////////////////////////////////////////////////

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
    AliMUONData * fMUONData; // data for MUON subsystem 

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
