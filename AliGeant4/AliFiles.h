// $Id$
// Category: global
//
// Class for file names and paths.
// It is protected from instantiating (only static data members
// and static methods are defined).

#ifndef ALI_FILES_H
#define ALI_FILES_H

#include <globals.hh>

#ifdef G4USE_STL
#include <string>
#endif

class AliFiles
{
  public:
    // --> protected
    // AliFiles();
    virtual ~AliFiles();

    // static get methods
    static G4String Config();   
    static G4String DetConfig1();
    static G4String DetConfig2();
    static G4String DetConfig3();
    static G4String DetConfig4();
    static G4String DetConfigName1();
    static G4String DetConfigName2();
    static G4String DetData1();
    static G4String DetData2();
    static G4String DetData3();
    static G4String STRUCT();

  protected:
    AliFiles();  
       // only static data members and methods
    
  private:       
    // static data members  
    static const G4String  fgTop;        //top directory
    static const G4String  fgConfig;     //path to general Config.C
    static const G4String  fgDetConfig1; //path (part 1) to module Config.C/in
    static const G4String  fgDetConfig2; //path (part 2) to module Config.C/in
    static const G4String  fgDetConfig3; //path (part 3) to module Config.C/in
    static const G4String  fgDetConfig4; //path (part 2) to module Config.C/in
    static const G4String  fgDetConfigName1;  //config macro name (part 1)
    static const G4String  fgDetConfigName2;  //config macro name (part 2)
    static const G4String  fgDetData1;   //path (part 1) to module g3calls.dat
    static const G4String  fgDetData2;   //path (part 2) to module g3calls.dat
    static const G4String  fgDetData3;   //path (part 3) to module g3calls.dat
    static const G4String  fgSTRUCT;     //structure directory name
};  

// inline methods

inline G4String AliFiles::Config()
{ return fgConfig; }

inline G4String AliFiles::DetConfig1()
{ return fgDetConfig1; }

inline G4String AliFiles::DetConfig2()
{ return fgDetConfig2; }

inline G4String AliFiles::DetConfig3()
{ return fgDetConfig3; }

inline G4String AliFiles::DetConfig4()
{ return fgDetConfig4; }

inline G4String AliFiles::DetConfigName1()
{ return fgDetConfigName1; }

inline G4String AliFiles::DetConfigName2()
{ return fgDetConfigName2; }

inline G4String AliFiles::DetData1()
{ return fgDetData1; }

inline G4String AliFiles::DetData2()
{ return fgDetData2; }

inline G4String AliFiles::DetData3()
{ return fgDetData3; }

inline G4String AliFiles::STRUCT()
{ return fgSTRUCT; }

#endif //ALI_FILES_H
