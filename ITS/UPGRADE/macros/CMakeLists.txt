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

# Macros in alphabetical order
set(MACRO_SRCS
    Config_pileup.C
    efficiency.C
    FastVsSlowSim.C
    PlotDoseExtrapolat.C
    PlotOptimizeLenght.C
    recGloCooked.C
    ScanDigitsSuze02_v15.C
    testTrackerCA.C
    testTrackerCooked.C
)

install(FILES ${MACRO_SRCS} DESTINATION macros/ITS/UPGRADE)