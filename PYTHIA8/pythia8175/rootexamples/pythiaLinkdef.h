#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
 
#pragma link C++ namespace Pythia8;
#pragma link C++ class Pythia8::Event+;
#pragma link C++ class Pythia8::Particle+;
#pragma link C++ class Pythia8::Junction+;
#pragma link C++ class Pythia8::Vec4+;

#ifdef PYTHIA8_COMPLETE_ROOT_DICTIONARY
// For tree.cc, none of the following (generator-internal)
// classes is needed. Thus excluded from the default Linkdef.h.
#pragma link C++ class Pythia8::ParticleData+;
#pragma link C++ class Pythia8::ParticleDataEntry+;
#pragma link C++ class Pythia8::DecayChannel+;
#pragma link C++ class Pythia8::Pythia+;
#pragma link C++ class Pythia8::CoupSM+;
#pragma link C++ class Pythia8::Couplings+;
#pragma link C++ class Pythia8::InBeam+;
#pragma link C++ class Pythia8::InPair+;
#pragma link C++ class Pythia8::AlphaStrong+;
#pragma link C++ class Pythia8::AlphaEM+;
#pragma link C++ class Pythia8::Info+;
#pragma link C++ class Pythia8::Rndm+;
#pragma link C++ class Pythia8::RndmEngine+;
#pragma link C++ class Pythia8::Settings+;
#pragma link C++ class Pythia8::ResonanceTop+;
#pragma link C++ class Pythia8::ResonanceGeneric+;
#pragma link C++ class Pythia8::ResonanceGmZ+;
#pragma link C++ class Pythia8::ResonanceW+;
#pragma link C++ class Pythia8::ResonanceFour+;
#pragma link C++ class Pythia8::ResonanceH+;
#pragma link C++ class Pythia8::ResonanceHchg+;
#pragma link C++ class Pythia8::ResonanceZprime+;
#pragma link C++ class Pythia8::ResonanceWprime+;
#pragma link C++ class Pythia8::ResonanceRhorizontal+;
#pragma link C++ class Pythia8::ResonanceExcited+;
#pragma link C++ class Pythia8::ResonanceGraviton+;
#pragma link C++ class Pythia8::ResonanceLeptoquark+;
#pragma link C++ class Pythia8::ResonanceNuRight+;
#pragma link C++ class Pythia8::ResonanceZRight+;
#pragma link C++ class Pythia8::ResonanceWRight+;
#pragma link C++ class Pythia8::ResonanceHchgchgLeft+;
#pragma link C++ class Pythia8::ResonanceHchgchgRight+;
#pragma link C++ class Pythia8::ResonanceWidths+;
#pragma link C++ class Pythia8::SigmaProcess+;
#pragma link C++ class std::vector<ResonanceWidths>+;
#pragma link C++ class std::pair<int,ParticleDataEntry>+;
#endif
#endif
