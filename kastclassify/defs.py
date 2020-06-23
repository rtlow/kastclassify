'''
Contains the dicts for the indices, EWs,
and chi factors.
'''

import numpy as np


#global strings
line_keyword = 'line'
wave_range_keyword = 'wave_range'
feat_type_keyword = 'feat_type'
name_keyword = 'name'

indices_keyword = 'indices'
continuum_keyword = 'continuum'
feature_keyword = 'feature'
method_keyword = 'method'
sample_keyword = 'sample'
ref_keyword = 'ref'


#index dict

index_dict = {
'lepine2007':
    
     #leipine (2007)
    
     {'indices': {\
        'CaH3': {'feature': [[6960,6990]],\
                 'continuum': [[7042,7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaH2': {'feature': [[6814, 6846]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO5': {'feature': [[7126,7135]],\
                 'continuum': [[7042,7046]],\
                 'method': 'ratio',\
                 'sample': 'average'}
        }},\
 
 'lepine2003':
    
     #lepine (2003)

    {'indices': {\
        'CaH1': {'feature': [[6380,6390]],\
                 'continuum': [[6410,6420], [6345, 6355]],\
                 'method': 'avdenom',\
                 'sample': 'average'},\
        'CaH2': {'feature': [[6814, 6846]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaH3': {'feature': [[6960,6990]],\
                 'continuum': [[7042,7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO5': {'feature': [[7126,7135]],\
                 'continuum': [[7042,7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'VO1': {'feature': [[7430, 7470]],\
                 'continuum': [[7550,7570]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO6': {'feature': [[7550,7570]],\
                 'continuum': [[7745,7765]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'VO2': {'feature': [[7920, 7960]],\
                 'continuum': [[8440, 8470]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO7': {'feature': [[8440, 8470]],\
                 'continuum': [[8400, 8420]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'Color-M': {'feature': [[8105, 8155]],\
                 'continuum': [[6510, 6560]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        }},\

 'gizis':
    
    #gizis (1997)

    {'indices': {\
        'TiO5': {'feature': [[7126, 7135]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaH1': {'feature': [[6380, 6390]],\
                 'continuum': [[6345, 6355], [6410, 6420]],\
                 'method': 'avdenom',\
                 'sample': 'average'},\
        'CaH2': {'feature': [[6814, 6846]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaH3': {'feature': [[6960, 6990]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'}
        }},\
    
 'reid':
     
    #reid (1995)

    {'indices': {\
        'TiO1': {'feature': [[6718, 6723]],\
                 'continuum': [[6703, 6708]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO2': {'feature': [[7058, 7061]],\
                 'continuum': [[7043, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO3': {'feature': [[7092, 7097]],\
                 'continuum': [[7079, 7084]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'TiO4': {'feature': [[7130, 7135]],\
                 'continuum': [[7115, 7120]],\
                 'method': 'ratio',
                 'sample': 'average'},\
        'TiO5': {'feature': [[7126, 7135]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaH1': {'feature': [[6380, 6390]],\
                 'continuum': [[6345, 6355], [6410, 6420]],\
                 'method': 'avdenom',
                 'sample': 'average'},\
        'CaH2': {'feature': [[6814, 6846]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaH3': {'feature': [[6960, 6990]],\
                 'continuum': [[7042, 7046]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        'CaOH': {'feature': [[6230, 6240]],\
                 'continuum': [[6345, 6354]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        
        }},\
 
 
 'kirkpatrick':
        
    #kirkpatrick (1999)

    {'indices': {\
        'Rb-a': {'feature': [[7775.2, 7785.2], [7815.2, 7825.2]],\
                 'continuum': [[7795.2, 7805.2]],\
                 'method': 'avnum',\
                 'sample': 'integrate'},\
        'Rb-b': {'feature': [[7922.6, 7932.6], [7962.6, 7972.6]],\
                 'continuum': [[7942.6, 7952.6]],\
                 'method': 'avnum',\
                 'sample': 'integrate'},\
        'Na-a': {'feature': [[8153.3, 8163.3]],\
                 'continuum': [[8178.3, 8188.3]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Na-b': {'feature': [[8153.3, 8183.3]],\
                 'continuum': [[8189.8, 8199.8]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Cs-a': {'feature': [[8496.1, 8506.1], [8536.1, 8546.1]],\
                 'continuum': [[8516.1, 8526.1]],\
                 'method': 'avnum',\
                 'sample': 'integrate'},\
        'Cs-b': {'feature': [[8918.5, 8928.5], [8958.3, 8968.3]],\
                 'continuum': [[8938.5, 8948.3]],\
                 'method': 'avnum',\
                 'sample': 'integrate'},\
        'TiO-a': {'feature': [[7033.0, 7048.0]],\
                 'continuum': [[7958.0, 7973.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'TiO-b': {'feature': [[8400.0, 8415.0]],\
                 'continuum': [[8435.0, 8470.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'VO-a': {'feature': [[7350.0, 7370.0], [7550.0, 7570.0]],\
                 'continuum': [[7430.0, 7470.0]],\
                 'method': 'sumnum',\
                 'sample': 'integrate'},\
        'VO-b': {'feature': [[7860.0, 7880.0], [8080.0, 8100.0]],\
                 'continuum': [[7960.0, 8000.0]],\
                 'method': 'sumnum',\
                 'sample': 'integrate'},\
        'CrH-a': {'feature': [[8580.0, 8600.0]],\
                 'continuum': [[8621.0, 8641.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'CrH-b': {'feature': [[9940.0, 9960.0]],\
                 'continuum': [[9970.0, 9990.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'FeH-a': {'feature': [[8660.0, 8680.0]],\
                 'continuum': [[8700.0, 8720.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'FeH-b': {'feature': [[9863.0, 9883.0]],\
                 'continuum': [[9908.0, 9928.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Color-a': {'feature': [[9800.0, 9850.0]],\
                 'continuum': [[7300.0, 7350.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Color-b': {'feature': [[9800.0, 9850.0]],\
                 'continuum': [[7000.0, 7050.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Color-c': {'feature': [[9800.0, 9850.0]],\
                 'continuum': [[8100.0, 8150.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Color-d': {'feature': [[9675.0, 9875.0]],\
                 'continuum': [[7350.0, 7550.0]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        
        }},\
 
 'burgasser':
        
    #burgasser (2003)

    {'indices': {\
        'CsI(A)': {'feature': [[8496.1, 8506.1], [8536.1, 8546.1]],\
                 'continuum': [[8516.1, 8626.1]],\
                 'method': 'sumnum_twicedenom',\
                 'sample': 'average'},\
        'CsI(B)': {'feature': [[8918.5, 8928.5], [8958.3, 8968.3]],\
                 'continuum': [[8938.5, 8948.3]],\
                 'method': 'sumnum_twicedenom',\
                 'sample': 'average'},\
        'H2O': {'feature': [[9220, 9240]],\
                 'continuum': [[9280, 9300]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'CrH(A)': {'feature': [[8560, 8600]],\
                 'continuum': [[8610, 8650]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'CrH(B)': {'feature': [[9855, 9885]],\
                 'continuum': [[9970, 10000]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'FeH(A)': {'feature': [[8560, 8600]],\
                 'continuum': [[8685, 8725]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'FeH(B)': {'feature': [[9855, 9885]],\
                 'continuum': [[9905, 9935]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'Color-e': {'feature': [[9140, 9240]],\
                 'continuum': [[8400, 8500]],\
                 'method': 'ratio',\
                 'sample': 'average'},\
        }},\
 
 'martin':
        
    #martin (1999)

    {'indices': {\
        'PC3': {'feature': [[8230, 8270]],\
                 'continuum': [[7540, 7580]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'PC6': {'feature': [[9090, 9130]],\
                 'continuum': [[6500, 6540]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'CrH1': {'feature': [[8560, 8600]],\
                 'continuum': [[8610, 8650]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'CrH2': {'feature': [[9840, 9880]],\
                 'continuum': [[9970, 10010]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'FeH1': {'feature': [[8560, 8600]],\
                 'continuum': [[8685, 8725]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'FeH2': {'feature': [[9840, 9880]],\
                 'continuum': [[9900, 9940]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'H2O1': {'feature': [[9190, 9230]],\
                 'continuum': [[9280, 9320]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'TiO1': {'feature': [[7000, 7040]],\
                 'continuum': [[7060, 7100]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'TiO2': {'feature': [[8380, 8420]],\
                 'continuum': [[8440, 8480]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'VO1': {'feature': [[7540, 7580]],\
                 'continuum': [[7420, 7460]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        'VO2': {'feature': [[7990, 8030]],\
                 'continuum': [[7900, 7940]],\
                 'method': 'ratio',\
                 'sample': 'integrate'},\
        }}
 
 
}
 
    
    
#EW info dict
EW_dict = {
    
    'RbIa': {
        'name': 'RbIa',\
        'line': 7800,\
        'wave_range': [7770, 7830],\
        'feat_type': 'absorption'},\
        
    'RbIb': {
        'name': 'RbIb',\
        'line': 7948,\
        'wave_range': [7918, 7978],\
        'feat_type': 'absorption'},\
        
    'NaIa': {
        'name': 'NaIa',\
        'line': 8183,\
        'wave_range': [8153, 8213],\
        'feat_type': 'absorption'},\
        
    'NaIb': {
        'name': 'NaIb',\
        'line': 8195,\
        'wave_range': [8165, 8225],\
        'feat_type': 'absorption'},\
    
    'CsIa': {
        'name': 'CsIa',\
        'line': 8521,\
        'wave_range': [8491, 8551],\
        'feat_type': 'absorption'},\
        
    'CsIb': {
        'name': 'CsIb',\
        'line': 8943,\
        'wave_range': [8913, 8973],\
        'feat_type': 'absorption'},\
    
    'Halpha': {
        'name': 'Halpha',\
        'line': 6563,\
        'wave_range': [6533, 6593],\
        'feat_type': 'emission'},\
        
    'LiI': {
        'name': 'LiI',\
        'line': 6708,\
        'wave_range': [6678, 6738],\
        'feat_type': 'absorption'},\
        
    'KIa': {
        'name': 'KIa',\
        'line': 7665,\
        'wave_range': [7635, 7695],\
        'feat_type': 'absorption'},\
        
    'KIb': {
        'name': 'KIb',\
        'line': 7699,\
        'wave_range': [7669, 7729],\
        'feat_type': 'absorption'},\
        
    'TiIa': {
        'name': 'TiIa',\
        'line': 6573,\
        'wave_range': [6543, 6603],\
        'feat_type': 'absorption'},\
        
    'TiIb': {
        'name': 'TiIb',\
        'line': 8542,\
        'wave_range': [8512, 8572],\
        'feat_type': 'absorption'},\
        
    'CaII': {
        'name': 'CaII',\
        'line': 8436,\
        'wave_range': [8406, 8466],\
        'feat_type': 'absorption'},\
        
    'CaIa': {
        'name': 'CaIa',\
        'line': 7209,\
        'wave_range': [7179, 7239],\
        'feat_type': 'absorption'},\
        
    'CaIb': {
        'name': 'CaIb',\
        'line': 7213,\
        'wave_range': [7183, 7243],\
        'feat_type': 'absorption'},\
        
    'CaIc': {
        'name': 'CaIc',\
        'line': 7326,\
        'wave_range': [7296, 7356],\
        'feat_type': 'absorption'},\
    
}


#chi factor dict

#calibrated chi values from Douglas, 2014 2014ApJ...795..161D for M dwarfs
# from Schmidt, 2014 2014PASP..126..642S for L dwarfs

chi_factor_types = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]

chi_factor_values = [6.6453e-5, 6.0334e-5, 5.2658e-5, 4.4872e-5, 3.5926e-5, 2.4768e-5,\
                    1.7365e-5, 1.2057e-5, 0.6122e-5, 0.3522e-5, 1.98e-6, 2.25e-6, 2.11e-6,\
                    1.67e-6, 1.16e-6, 1.46e-6, 1.23e-6, 0.73e-6]

chi_factor_values_err = [0.6207e-5, 0.5326e-5, 0.5963e-5, 0.4967e-5, 0.5297e-5, 0.4860e-5,\
                        0.3475e-5, 0.3267e-5, 0.2053e-5, 0.1432e-5, 0.27e-6, 0.11e-6, 0.36e-6,\
                        0.22e-6, np.nan, 0.28e-6, np.nan, np.nan]

