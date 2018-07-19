///
/// \file AliFemtoEventCut.h
///

#ifndef AliFemtoEventCut_hh
#define AliFemtoEventCut_hh

class AliFemtoEvent;
class AliFemtoAnalysis;

#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoString.h"

/// \class AliFemtoEventCut
/// \brief The pure virtual base class for the event cut
///
/// All event cuts must inherit from this one and implement the ::Pass and
/// ::Report methods. The ::Clone() function simply returns NULL, so if users
/// want their cuts to behave as expected, they should also write their own.
///
class AliFemtoEventCut : public AliFemtoCutMonitorHandler {

  friend class AliFemtoAnalysis;

public:

  /// Default constructor
  ///
  /// Parent analysis set to NULL.
  AliFemtoEventCut();

  /// Simply copies pointer to parent analysis
  AliFemtoEventCut(const AliFemtoEventCut& c);

  /// no-op destructor
  virtual ~AliFemtoEventCut();

  /// Simply copies pointer to parent analysis
  AliFemtoEventCut& operator=(const AliFemtoEventCut& aCut);

  /// true if event passes, false if not
  ///
  /// This is an abstract method and MUST be implemented by a subclass.
  virtual bool Pass(const AliFemtoEvent* event) = 0;

  /// Return new settings list.
  ///
  /// This method creates a new list of TObjStrings describing cut parameters.
  /// The default implementation automatically calls the AppendSettings method
  /// to fill the list, so users only need to overload that method.
  virtual TList* ListSettings() const;


  /// Appends cut settings to a TList
  ///
  /// This method should be overloaded by the user to add any relevent settings
  /// of the cut to the list
  ///
  /// No settings are added by this class. Simply returns the incoming TList.
  ///
  /// \param A list to append settings to.
  /// \param prefix An optional prefix to prepend to the beginning of each setting
  /// \return The same pointer as the parameter
  virtual TList* AppendSettings(TList*, const TString& prefix="") const;

  virtual AliFemtoString Report() = 0; ///< A user-written method to return a string describing cuts
  virtual AliFemtoEventCut* Clone() const;   ///< Returns NULL - users should overload.

  /// Returns the analysis this cut belongs to
  AliFemtoAnalysis* HbtAnalysis();

  /// Sets the analysis this cut will belong to
  void SetAnalysis(AliFemtoAnalysis* aAnalysis);

protected:

  AliFemtoAnalysis* fyAnalysis; //!<! Pointer to the 'parent' analysis

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoEventCut, 2);
  /// \endcond
#endif
};

inline AliFemtoEventCut::AliFemtoEventCut():
  AliFemtoCutMonitorHandler(),
  fyAnalysis(NULL)
{
}

inline AliFemtoEventCut::AliFemtoEventCut(const AliFemtoEventCut& c):
  AliFemtoCutMonitorHandler(),
  fyAnalysis(c.fyAnalysis)
{
}

inline AliFemtoEventCut& AliFemtoEventCut::operator=(const AliFemtoEventCut& aCut)
{
  if (this == &aCut) {
    return *this;
  }
  fyAnalysis = aCut.fyAnalysis;
  return *this;
}

inline AliFemtoEventCut::~AliFemtoEventCut()
{ // no-op
}

inline AliFemtoAnalysis* AliFemtoEventCut::HbtAnalysis()
{
  return fyAnalysis;
}


inline void AliFemtoEventCut::SetAnalysis(AliFemtoAnalysis* analysis)
{
  fyAnalysis = analysis;
}

inline AliFemtoEventCut* AliFemtoEventCut::Clone() const
{
  return NULL;
}

inline TList* AliFemtoEventCut::ListSettings() const
{
  return AppendSettings(new TList());
}

inline TList* AliFemtoEventCut::AppendSettings(TList *setting_list,
                                               const TString &prefix) const
{
  return setting_list;
}


#endif
