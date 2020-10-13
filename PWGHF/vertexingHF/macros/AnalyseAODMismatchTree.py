import uproot
import numpy as np
import pandas as pd
from ROOT import TH1F, TH2F, TCanvas, TLegend
from ROOT import kRed, kAzure, gStyle, kIsland

def GetMaskOfBits(bits):
    '''
    Helper method to get bit mask from bits

    Arguments
    ----------

    - list of bits

    Returns
    ----------

    - mask corresponding to the input bits
    '''
    mask = 0
    for bit in bits:
        mask += 2**bit

    return mask


def FilterBitDf(dfToFilter, column, bitsToTest, logic='or'):
    '''
    Method to apply selection testing one or more bits

    Arguments
    ----------

    - pandas dataframe to filter
    - colum with bitmap
    - list of bits to test
    - logic to combine the bits (and, or)

    Returns
    ----------
    - filtered pandas dataframe
    '''
    maskOfBits = GetMaskOfBits(bitsToTest)
    flags = dfToFilter[column].astype(int) & maskOfBits
    if logic == 'or':
        flags = flags.astype('bool')
    elif logic == 'and':
        flags -= maskOfBits
        flags = ~flags.astype('bool')
    elif logic == 'not':
        flags = ~flags.astype('bool')
    else:
        print('Error: only and, or, and not logics are supported for bitwise operations')
        return None

    dfFilt = dfToFilter[flags.values]

    return dfFilt


def main():
    '''
    Main function
    '''

    prod = 'LHC20g11a'

    tree = uproot.open('AnalysisResults.root')['AOD_dAOD_Matching/fTreeMismatch']
    df = tree.pandas.df()
    df = df.sort_values(by=['file_name'])
    pd.set_option('display.max_colwidth', None)

    nFiles = len(df)
    nEvents = sum(df['n_events'].values)

    dfSel = {'good_files': df.query('mismatch_status == 0'),
             'mism_ev': FilterBitDf(df, 'mismatch_status', [0]),
             'mism_TProcessID': FilterBitDf(df, 'mismatch_status', [1]),
             'mism_cand': FilterBitDf(df, 'mismatch_status', [2]),
             'mism_ev_and_TProcessID': FilterBitDf(df, 'mismatch_status', [0, 1], logic='and'),
             'mism_ev_and_cand': FilterBitDf(df, 'mismatch_status', [0, 2], logic='and'),
             'mism_cand_and_TProcessID': FilterBitDf(df, 'mismatch_status', [1, 2], logic='and'),
             'mism_all': FilterBitDf(df, 'mismatch_status', [1, 2, 3], logic='and')}

    fracFiles, fracEv = {}, {}
    for mism in dfSel:
        fracFiles[mism] = len(dfSel[mism]) / nFiles
        fracEv[mism] = sum(dfSel[mism]['n_events'].values) / nEvents
        print(f'\nfraction of files with flag \"{mism}\": {fracFiles[mism]}')
        print(f'fraction of events with flag \"{mism}\": {fracEv[mism]}')
    
    gStyle.SetTitleSize(0.045, 'xy')
    gStyle.SetLabelSize(0.04, 'xy')
    gStyle.SetPadTopMargin(0.035)
    gStyle.SetPadRightMargin(0.035)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetOptStat(0)
    gStyle.SetPalette(kIsland)
    
    hAODMism = TH1F('hAODMism', ';;fraction', 8, 0.5, 8.5)
    hAODMism.SetLineWidth(2)
    hAODMism.SetLineColor(kRed+1)
    hAODMism.GetYaxis().SetRangeUser(1.e-5, 1.)
    hEventMism = TH1F('hEventMism', ';;fraction', 8, 0.5, 8.5)
    hEventMism.SetLineWidth(2)
    hEventMism.SetLineColor(kAzure+4)
    hEventMism.GetYaxis().SetRangeUser(1.e-5, 1.)
    for iMism, mism in enumerate(dfSel):
        hAODMism.GetXaxis().SetBinLabel(iMism+1, mism)
        hEventMism.GetXaxis().SetBinLabel(iMism+1, mism)
        hAODMism.SetBinContent(iMism+1, fracFiles[mism])
        hEventMism.SetBinContent(iMism+1, fracEv[mism])
    
    leg = TLegend(0.6, 0.7, 0.8, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hAODMism, 'AOD files', 'l')
    leg.AddEntry(hEventMism, 'events', 'l')

    cMismFrac = TCanvas('cMismFrac', '', 1920, 1080)
    cMismFrac.SetLogy()
    hAODMism.Draw()
    hEventMism.Draw('same')
    leg.Draw()
    cMismFrac.Modified()
    cMismFrac.Update()

    dfSel['mism_cand'][['file_name']].to_csv(f'AOD_mismatch_{prod}_cand.txt', header=False, index=False)
    dfSel['mism_ev'][['file_name']].to_csv(f'AOD_mismatch_{prod}_nevents.txt', header=False, index=False)
    dfSel['mism_TProcessID'][['file_name']].to_csv(f'AOD_mismatch_{prod}_TProcessID.txt', header=False, index=False)

    cMismFrac.SaveAs(f'AODMismatch_fractions_{prod}.pdf')

    # check for files not tested (jobs failed)
    runs = np.unique(df['run_number'].values)
    nRuns = len(runs)
    for iRun, run in enumerate(runs):
        dfRunSel = df.query(f'run_number == {run}')
        lastProcessedFile = list(dfRunSel['file_name'].values)[-1]
        numLastProcFile = int(lastProcessedFile.decode().rpartition('AOD/')[2].rpartition('/')[0])
        hFilesTested = TH2F(f'hFilesTested{run}', f'run {run};AOD number;', numLastProcFile, 0.5, numLastProcFile+0.5, 1, 0., 1.)
        hFilesTested.GetZaxis().SetRangeUser(-0.001, 1.)
        cFilesTested = TCanvas(f'cFilesTested{run}', '', 1920, 1080)
        cFilesTested.SetTopMargin(0.12)
        cFilesTested.SetRightMargin(0.12)
        for fileName in dfRunSel['file_name']:
            numProcFile = int(fileName.decode().rpartition('AOD/')[2].rpartition('/')[0])
            hFilesTested.Fill(numProcFile, 0.5)
        hFilesTested.Draw('colz')
        cFilesTested.Modified()
        cFilesTested.Update()
        if iRun == 0:
            cFilesTested.SaveAs(f'FilesTested_{prod}.pdf[')
        cFilesTested.SaveAs(f'FilesTested_{prod}.pdf')
        if iRun == nRuns-1:
            cFilesTested.SaveAs(f'FilesTested_{prod}.pdf]')
    input()

# call main function
main()