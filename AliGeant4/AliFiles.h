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
    static const G4String  fgkTop;        //top directory
    static const G4String  fgkConfig;     //path to general Config.C
    static const G4String  fgkDetConfig1; //path (part 1) to module Config.C/in
    static const G4String  fgkDetConfig2; //path (part 2) to module Config.C/in
    static const G4String  fgkDetConfig3; //path (part 3) to module Config.C/in
    static const G4String  fgkDetConfig4; //path (part 2) to module Config.C/in
    static const G4String  fgkDetConfigName1;  //config macro name (part 1)
    static const G4String  fgkDetConfigName2;  //config macro name (part 2)
    static const G4String  fgkDetData1;   //path (part 1) to module g3calls.dat
    static const G4String  fgkDetData2;   //path (part 2) to module g3calls.dat
    static const G4String  fgkDetData3;   //path (part 3) to module g3calls.dat
    static const G4String  fgkSTRUCT;     //structure directory name
};  

// inline methods

inline G4String AliFiles::Config()
{ return fgkConfig; }

inline G4String AliFiles::DetConfig1()
{ return fgkDetConfig1; }

inline G4String AliFiles::DetConfig2()
{ return fgkDetConfig2; }

inline G4String AliFiles::DetConfig3()
{ return fgkDetConfig3; }

inline G4String AliFiles::DetConfig4()
{ return fgkDetConfig4; }

inline G4String AliFiles::DetConfigName1()
{ return fgkDetConfigName1; }

inline G4String AliFiles::DetConfigName2()
{ return fgkDetConfigName2; }

inline G4String AliFiles::DetData1()
{ return fgkDetData1; }

inline G4String AliFiles::DetData2()
{ return fgkDetData2; }

inline G4String AliFiles::DetData3()
{ return fgkDetData3; }

inline G4String AliFiles::STRUCT()
{ return fgkSTRUCT; }

#endif //ALI_FILES_H
