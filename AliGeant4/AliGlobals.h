// $Id$
// Category: global
//
// Class for generally used basic functions.
// It is protected from instantiating (only static data members
// and static methods are defined).


#ifndef ALI_GLOBALS_H
#define ALI_GLOBALS_H

#include <globals.hh>

#ifdef G4USE_STL
#include <string>
#endif

class AliGlobals
{
  public:
    // --> protected
    // AliGlobals();
    virtual ~AliGlobals();

    // static methods
    static void Exception(const char* s=0);
      // Global error function prints string to cerr, and aborts
      // program - according to G4Exception.cc
    static void Warning(const char* s=0);
      // Global warning function prints string to cerr
#ifdef G4USE_STL
    static void Exception(G4std::string s);
    static void Exception(G4String s);
    static void Warning(G4std::string s);
    static void Warning(G4String s);
#endif

    static void AppendNumberToString(G4String& string, G4int number);
    static G4int StringToInt(G4String string);
    
    // static get methods
    static G4double DefaultCut();

  protected:
    AliGlobals();  
       // only static data members and methods
    
  private:       
    // static data members  
    static const G4double  fgkDefaultCut; //default cut value
};  

// inline methods

inline G4double AliGlobals::DefaultCut()
{ return fgkDefaultCut; }

#endif //ALI_GLOBALS_H
