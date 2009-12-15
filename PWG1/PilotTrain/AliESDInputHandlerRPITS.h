#ifndef ALIESDINPUTHANDLERRPITS_H
#define ALIESDINPUTHANDLERRPITS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliESDInputHandler.h 24521 2008-03-14 16:43:54Z morsch $ */

//-------------------------------------------------------------------------
//     ESD Input Handler realisation of the AliVEventHandler interface
//     Automatic loading of RecPoint Trees
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliESDInputHandlerRP.h"
#include "AliESDEvent.h"
class TList;
class TTree;
class TDirectoryFile;
class TString;


class AliESDInputHandlerRPITS : public AliESDInputHandlerRP {

 public:
    AliESDInputHandlerRPITS();
    AliESDInputHandlerRPITS(const char* name, const char* title);
    virtual ~AliESDInputHandlerRPITS();
    virtual Bool_t       Notify(const char* path);
    ClassDef(AliESDInputHandlerRPITS, 1);
};

#endif
