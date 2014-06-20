includedir=$(prefix)/include/Photos

lib_LTLIBRARIES = libPhotosFortran.la libPhotosCxxInterface.la

libPhotosFortran_la_SOURCES = photos-fortran/photosC.cxx \
photos-fortran/forW-MEc.cxx                     \
photos-fortran/forZ-MEc.cxx

################################################################################

libPhotosCxxInterface_la_SOURCES =              \
eventRecordInterfaces/PhotosHEPEVTEvent.cxx     \
eventRecordInterfaces/PhotosHEPEVTParticle.cxx  \
photosCInterfaces/Photos.cxx                    \
photosCInterfaces/PhotosParticle.cxx            \
photosCInterfaces/PhotosBranch.cxx              \
photosCInterfaces/PhotosEvent.cxx               \
photosFortranInterfaces/PH_HEPEVT_Interface.cxx \
utilities/Log.cxx                               \
utilities/PhotosUtilities.cxx                   \
utilities/PhotosRandom.cxx

include_HEADERS =                             \
eventRecordInterfaces/PhotosHEPEVTParticle.h  \
eventRecordInterfaces/PhotosHEPEVTEvent.h     \
photosCInterfaces/Photos.h                    \
photosCInterfaces/PhotosParticle.h            \
photosCInterfaces/PhotosBranch.h              \
photosCInterfaces/PhotosEvent.h               \
photosFortranInterfaces/PH_HEPEVT_Interface.h \
photosFortranInterfaces/f_Init.h              \
utilities/Log.h                               \
utilities/PhotosRandom.h

if HAS_HEPMC
  libPhotosCxxInterface_la_SOURCES += \
  eventRecordInterfaces/PhotosHepMCEvent.cxx    \
  eventRecordInterfaces/PhotosHepMCParticle.cxx

  include_HEADERS += \
  eventRecordInterfaces/PhotosHepMCParticle.h  \
  eventRecordInterfaces/PhotosHepMCEvent.h
endif

INCLUDES = -I$(top_srcdir)/src/photosFortranInterfaces \
           -I$(top_srcdir)/src/photosCInterfaces \
           -I$(top_srcdir)/src/utilities \
           -I$(HEPMC_DIR)/include

libPhotosCxxInterface_la_LIBADD = libPhotosFortran.la
