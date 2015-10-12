///
/// \file AliFemtoHiddenInfo.h
///

#ifndef AliFemtoHiddenInfo_hh
#define AliFemtoHiddenInfo_hh

/// \class AliFemtoHiddenInfo
/// \brief A pure virtual base class for the hidden (Monte Carlo) data.
///
/// Hidden info stores additional information, which is not in a standard
/// track.
///
class AliFemtoHiddenInfo {
public:

  /// Trivial Constructor
  AliFemtoHiddenInfo();

  /// Trivial Destructor
  virtual ~AliFemtoHiddenInfo();

  /// !!! MANDATORY !!!
  /// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* Clone() const;

protected:
  // This is called by Clone to do the actual copying and return hidden info
  virtual AliFemtoHiddenInfo* GetParticleHiddenInfo() const = 0;

};

inline AliFemtoHiddenInfo::AliFemtoHiddenInfo()
{ // no-op
}

inline AliFemtoHiddenInfo::~AliFemtoHiddenInfo()
{ // no-op
}

//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoHiddenInfo::Clone() const
{
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
