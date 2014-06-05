#ifndef HEPMC_DEFS_H
#define HEPMC_DEFS_H
// ----------------------------------------------------------------------
//
// HepMCDefs.h
// Author:  Lynn Garren
//
//  Define various useful macros
//  Most allow users to check for various code features
//
// ----------------------------------------------------------------------

// the HeavyIon class is available in HepMC  
#ifndef HEPMC_HAS_HEAVY_ION
#define HEPMC_HAS_HEAVY_ION
#endif

// the IO_Ascii class is NOT available in HepMC   
#ifndef HEPMC_IO_ASCII_REMOVED
#define HEPMC_IO_ASCII_REMOVED
#endif

// the ParticleData class is NOT available in HepMC   
#ifndef HEPMC_PARTICLE_DATA_REMOVED
#define HEPMC_PARTICLE_DATA_REMOVED
#endif

// the IO_GenEvent class is available in HepMC   
#ifndef HEPMC_HAS_IO_GENEVENT
#define HEPMC_HAS_IO_GENEVENT
#endif

// the PdfInfo class is available in HepMC  
#ifndef HEPMC_HAS_PDF_INFO
#define HEPMC_HAS_PDF_INFO
#endif

// HepMC uses SimpleVector (FourVector) to store momentum and position  
#ifndef HEPMC_HAS_SIMPLE_VECTOR
#define HEPMC_HAS_SIMPLE_VECTOR
#endif

// units are defined in HepMC  
#ifndef HEPMC_HAS_UNITS
#define HEPMC_HAS_UNITS
#endif

// the GenCrossSection class is available in HepMC  
#ifndef HEPMC_HAS_CROSS_SECTION
#define HEPMC_HAS_CROSS_SECTION
#endif

// the iterator range classes are available in HepMC  
#ifndef HEPMC_HAS_ITERATOR_RANGES
#define HEPMC_HAS_ITERATOR_RANGES
#endif

// the HepMC::WeightContainer class allows named weights
#ifndef HEPMC_HAS_NAMED_WEIGHTS
#define HEPMC_HAS_NAMED_WEIGHTS
#endif

// define the version of HepMC. 
#ifndef HEPMC_VERSION
#define HEPMC_VERSION "2.06.08"
#endif

#endif  // HEPMC_DEFS_H
