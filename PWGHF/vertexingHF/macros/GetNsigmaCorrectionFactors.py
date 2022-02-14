'''
Python script for the computation of the Nsigma post-calibration factors
Gets as input the output file of https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/EstimateSingleTrackPIDsyst.C

run: python GetNsigmaCorrectionFactors.py inFile.root
'''

import argparse
from ROOT import TFile

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileName', metavar='text', default='inFile.root')
args = parser.parse_args()

species = ['Pion', 'Kaon', 'Proton']
tags = {'Pion': 'V0', 'Kaon': 'TOF', 'Proton': 'V0'}
hMeanVsEtaVsP, hSigmaVsEtaVsP = {}, {}
inFile = TFile.Open(args.inFileName)
for specie in species:
    hMeanVsEtaVsP[specie] = inFile.Get(f'Data/hMean{specie}TPCData{tags[specie]}tagVsEtaVsP')
    hSigmaVsEtaVsP[specie] = inFile.Get(f'Data/hSigma{specie}TPCData{tags[specie]}tagVsEtaVsP')

nPBins = hMeanVsEtaVsP["Pion"].GetXaxis().GetNbins()
nEtaBins = hMeanVsEtaVsP["Pion"].GetYaxis().GetNbins()

# pylint: disable=anomalous-backslash-in-string,line-too-long
outFileTxt = open(args.inFileName.replace('.root', '.txt'), 'w')
outFileTxt.write(f'    nPbins = {nPBins};\n')
outFileTxt.write('    vector<Float_t> pTPClims = {')
for iP in range(1, nPBins+1):
    outFileTxt.write(f'{hMeanVsEtaVsP["Pion"].GetXaxis().GetBinLowEdge(iP)}, ')
outFileTxt.write(f'{hMeanVsEtaVsP["Pion"].GetXaxis().GetBinUpEdge(iP)}')
outFileTxt.write('};\n')
outFileTxt.write(f'    nEtabins = {nEtaBins};\n')
outFileTxt.write('    vector<Float_t> absetalims = {')
for iEta in range(1, nEtaBins+1):
    outFileTxt.write(f'{hMeanVsEtaVsP["Pion"].GetYaxis().GetBinLowEdge(iEta)}, ')
outFileTxt.write(f'{hMeanVsEtaVsP["Pion"].GetYaxis().GetBinUpEdge(iEta)}')
outFileTxt.write('};\n\n')

# means
for specie in species:
    outFileTxt.write(f'    vector<vector<Float_t> > mean{specie} = ')
    outFileTxt.write('{\n')
    for iEta in range(1, nEtaBins+1):
        if specie != 'Proton':
            outFileTxt.write('                                         {')
        else:
            outFileTxt.write('                                           {')
        for iP in range(1, nPBins+1):
            if iP != nPBins:
                outFileTxt.write(f'{hMeanVsEtaVsP[specie].GetBinContent(iP, iEta):.6f}, ')
            else:
                outFileTxt.write(f'{hMeanVsEtaVsP[specie].GetBinContent(iP, iEta):.6f}')
        if iEta != nEtaBins:
            outFileTxt.write('},\n')
        else:
            outFileTxt.write('}\n')
    if specie != 'Proton':
        outFileTxt.write('                                        };\n')
    else:
        outFileTxt.write('                                          };\n')

# widths
for specie in species:
    outFileTxt.write(f'    vector<vector<Float_t> > sigma{specie} = ')
    outFileTxt.write('{\n')
    for iEta in range(1, nEtaBins+1):
        if specie != 'Proton':
            outFileTxt.write('                                          {')
        else:
            outFileTxt.write('                                            {')
        for iP in range(1, nPBins+1):
            if iP != nPBins:
                outFileTxt.write(f'{hSigmaVsEtaVsP[specie].GetBinContent(iP, iEta):.6f}, ')
            else:
                outFileTxt.write(f'{hSigmaVsEtaVsP[specie].GetBinContent(iP, iEta):.6f}')
        if iEta != nEtaBins:
            outFileTxt.write('},\n')
        else:
            outFileTxt.write('}\n')
    if specie != 'Proton':
        outFileTxt.write('                                         };\n')
    else:
        outFileTxt.write('                                           };\n')

outFileTxt.close()

print(f'File saved in {args.inFileName.replace(".root", ".txt")}')
