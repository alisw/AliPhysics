#ifndef ALIMUONLOADER_H
#define ALIMUONLOADER_H

#include "AliLoader.h"
#include "AliDataLoader.h"

//__________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliMUONLoader                                            //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliMUONLoader : public AliLoader {
 public:
    AliMUONLoader();
    AliMUONLoader(const Char_t *detname,const Char_t *eventfoldername); //contructor with name of the top folder of the tree
    AliMUONLoader(const Char_t *detname,TFolder* eventfolder);

    virtual ~AliMUONLoader();
   
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
 public:
 
    ClassDef(AliMUONLoader,1)
 };

#endif
