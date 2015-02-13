# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************

# Check for daqDA library in order to build detectors DA

set(daqDA_FOUND false)

if(daqDA)
    # Check for header
    find_path(DAQDAH daqDA.h PATHS ${daqDA})
    
    if(DAQDAH-NOTFOUND)
        message(FATAL_ERROR "daqDA enabled but daqDA.h not found. Please check that daqDA points to your installation")
    endif(DAQDAH-NOTFOUND)
    mark_as_advanced(DAQDAH)

    find_path(DAQDALIB libdaqDA.a PATHS ${daqDA})
    
    if(DAQDALIB-NOTFOUND)
        message(FATAL_ERROR "daqDA enabled but libdaqDA.a not found. Please check that daqDA points to your installation")
    endif(DAQDALIB-NOTFOUND)
    mark_as_advanced(DAQDALIB)

    set(daqDA_FOUND TRUE)
endif(daqDA)