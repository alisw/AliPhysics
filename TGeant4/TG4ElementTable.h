// $Id$
// Category: geometry
//
// Singleton class that builds the default element table
// Currently defined elements: up to Z=47
// to be continued

#ifndef TG4_ELEMENT_TABLE_H
#define TG4_ELEMENT_TABLE_H

class TG4ElementTable 
{
  public:
    // --> protected
    // TG4ElementTable();
    // TG4ElementTable(const TG4ElementTable& right);
    virtual ~TG4ElementTable();
    
    // static methods
    static TG4ElementTable* Instance();

  protected:
    TG4ElementTable();    
    TG4ElementTable(const TG4ElementTable& right); 
    
    // operators
    TG4ElementTable& operator=(const TG4ElementTable& right);
          
  private:
    void Construct();

    // data members
    static TG4ElementTable*  fgInstance; //this instance
};   

#endif //G4_ELEMENT_TABLE_H
