#ifndef FLUGGNAVIGATOR_H
#define FLUGGNAVIGATOR_H 1

//To have access to the private data members, let's play a dirty trick
#define private protected
#include "G4Navigator.hh"
#undef private


class FluggNavigator : public G4Navigator
{
  public:

  friend std::ostream& operator << (std::ostream &os, const FluggNavigator &n);

  FluggNavigator();
    // Constructor - initialisers and setup.

  virtual ~FluggNavigator() {}
    // Destructor. No actions.

  // flugg member function: reinitialization of navigator history with 
  // secondary particle history banked on fluka side
  void UpdateNavigatorHistory(const G4NavigationHistory* newNavHistory);

};

#endif
