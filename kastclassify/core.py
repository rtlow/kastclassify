'''
DEPENDENCIES:

kastredux
numpy
matplotlib
pandas
scipy
astropy
specutils
splat
'''

#Internal Imports
import glob
import os
import datetime

#Standard External Imports
import numpy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import copy

import astropy.nddata
import astropy.io.fits as fits
from astropy import units as u
from astropy.visualization import quantity_support

quantity_support()

from scipy.interpolate import interp1d
from scipy.integrate import trapz


#Unique External Imports
import kastredux

import specutils
import specutils.fitting

from splat.utilities import typeToNum


####################################################
### END IMPORTS ####################################
####################################################


####################################################
### BEGIN ANALYSIS FUNCTIONS #######################
####################################################

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


#chi2 classification functions
def spt_from_name(name):
    
    spt = name.split('_')[0]
    
    spt = spt[:-1] + '.' + spt[-1:]
    
    return spt

def classify_by_standard(spec, ref='lepine', plot=False, diag_path=None):
    
    lepine_str = 'lepine'
    kirkpatrick_str = 'kirkpatrick'
    burgasser_str = 'burgasser'
    sdss_str = 'SDSS'
    all_str = 'all'
    
    STANDARD_SPECTRA = []
    LEPINE_STANDARDS = []
    KIRKPATRICK_STANDARDS = []
    BURGASSER_STANDARDS = []
    SDSS_STANDARDS = []
    
    STANDARD_PATHS = glob.glob(kastredux.SPTSTDFOLDER + '*.txt')

    
    for path in STANDARD_PATHS:
        name = os.path.basename(os.path.splitext(path)[0])
        stand_spec = kastredux.readSpectrum(path)
        stand_spec.name = name
        
        STANDARD_SPECTRA.append(stand_spec)
        
        if (lepine_str in name):
            LEPINE_STANDARDS.append(stand_spec)
        elif (kirkpatrick_str in name):
            KIRKPATRICK_STANDARDS.append(stand_spec)
        elif (burgasser_str in name):
            BURGASSER_STANDARDS.append(stand_spec)
        elif sdss_str in name:
            SDSS_STANDARDS.append(stand_spec)
            


    
    minimal = 1e90
    minimal_scale = None
    minimal_standard = None
    minimal_standard_name = None
    
    COMPARISON = []
    
    if (ref in lepine_str):
        COMPARISON = LEPINE_STANDARDS
    elif (ref in kirkpatrick_str):
        COMPARISON = KIRKPATRICK_STANDARDS
    elif (ref in burgasser_str):
        COMPARISON = BURGASSER_STANDARDS
    elif ref in sdss_str:
        COMPARISON = SDSS_STANDARDS
    else:
        COMPARISON = STANDARD_SPECTRA
        
    for stand in COMPARISON:
        
        chi2, scalefactor = kastredux.compareSpectra_simple(spec, stand)

        
        if chi2 < minimal:
            minimal = chi2
            minimal_scale = scalefactor
            minimal_standard = stand
            minimal_standard_name = stand.name
            
    if plot or (diag_path is not None):
        placeholder = kastredux.compareSpectra_simple(spec, minimal_standard, plot=True)
        
        if diag_path is not None:
            
            plt.savefig(diag_path + spec.name + '_standardComparison.png' )
        
        if plot:
            plt.show()
        plt.close()
     
    spt = spt_from_name(minimal_standard_name)

    return spt, minimal_standard_name, minimal, minimal_scale


#index measurement functions

#measure a single index
def measure_index(waves, fluxes, uncs, continuum_ranges, feature_ranges, method='ratio', sample='integrate', plot=False, name=None, diag_path=None):
    
    #number of samples in interpolation
    no_samples = 10000
    
    #creating flux interpolation
    fluxinterp = interp1d(waves, fluxes, bounds_error=False, fill_value=0.)
    uinterp = interp1d(waves, uncs, bounds_error=False, fill_value=0.)
    
    continuum_waves = []
    continuum_fluxes = []
    continuum_uncs = []
    
    #getting the continuum flux
    for continuum_range in continuum_ranges:
            
        min_wavelength = continuum_range[0]
        max_wavelength = continuum_range[1]
            
        #no_divisions = (max_wavelength - min_wavelength) * no_samples
            
        w = numpy.linspace(min_wavelength, max_wavelength, no_samples)
        
        continuum_waves.append(w)
        
        f = fluxinterp(w)
        
        u = uinterp(w)
            
        continuum_flux = f[ (w >= min_wavelength) & (w <= max_wavelength) ]
        
        continuum_unc = u[ (w >= min_wavelength) & (w <= max_wavelength) ]
            
        continuum_fluxes.append(continuum_flux)
        
        continuum_uncs.append(continuum_unc)
    
    #check if there is any flux in the continuum regions
    for continuum_flux, continuum_unc in zip(continuum_fluxes, continuum_uncs):
        
        if np.all(np.abs(continuum_flux) < continuum_unc):
            print('Warning! No flux in continuum region for {}!'.format(name))
            return 'No Flux in Continuum'
    
    feature_waves = []
    feature_fluxes = []
    feature_uncs = []
    
    #getting the feature flux
    for feature_range in feature_ranges:
            
        min_wavelength = feature_range[0]
        max_wavelength = feature_range[1]
            
        #no_divisions = (max_wavelength - min_wavelength) * no_samples
            
        w = numpy.linspace(min_wavelength, max_wavelength, no_samples)
        
        feature_waves.append(w)
        
        f = fluxinterp(w)
        u = uinterp(w)
            
        feature_flux = f[ (w >= min_wavelength) & (w <= max_wavelength) ]
        feature_unc = u[ (w >= min_wavelength) & (w <= max_wavelength) ]
            
        feature_fluxes.append(feature_flux)
        feature_uncs.append(feature_unc)
     
    #check if there is no flux in the feature regions
    for feature_flux, feature_unc in zip(feature_fluxes, feature_uncs):
        
        if np.all( np.abs(feature_flux) < feature_unc):
            print('Warning! No flux in feature region for {}!'.format(name))
            return 'No Flux in Feature'

    #plot the regions
    if plot or (diag_path is not None):
        
        for waves, fluxes in zip(continuum_waves, continuum_fluxes):
            plt.figure(figsize=[12,7])
            plt.plot(waves, fluxes)
            plt.xlabel('Wavelength (Angstrom)')
            plt.ylabel('Flux (erg/(Angstrom cm2 s))')
            plt.title('{} Continuum'.format(name))
            
            if diag_path is not None:
                plt.savefig(diag_path + '{}_index_Continuum.png'.format(name))
            
            if plot:
                
                plt.show()
            plt.close()
            
        for waves, fluxes in zip(feature_waves, feature_fluxes):
            plt.figure(figsize=[12,7])
            plt.plot(waves, fluxes)
            plt.xlabel('Wavelength (Angstrom)')
            plt.ylabel('Flux (erg/(Angstrom cm2 s))')
            plt.title('{} Feature'.format(name))
            if diag_path is not None:
                plt.savefig(diag_path + '{}_index_Feature.png'.format(name))
            
            if plot:
                
                plt.show()
            plt.close()
    
    #calculate samples
    continuum_values = []
    feature_values = []
    
    #only average and integrate are implemented
    if 'average' in sample:
        
        for waves, fluxes in zip(continuum_waves, continuum_fluxes):
            
            continuum_value = np.mean(fluxes)
            continuum_values.append(continuum_value)
            
        for waves, fluxes in zip(feature_waves, feature_fluxes):
            
            feature_value = np.mean(fluxes)
            feature_values.append(feature_value)
            
    elif 'integrate' in sample:
        
        for waves, fluxes in zip(continuum_waves, continuum_fluxes):
            
            continuum_value = trapz(fluxes, waves)
            continuum_values.append(continuum_value)
            
        for waves, fluxes in zip(feature_waves, feature_fluxes):
            
            feature_value = trapz(fluxes, waves)
            feature_values.append(feature_value)
    
    elif 'sum' in sample:
        
        for waves, fluxes in zip(continuum_waves, continuum_fluxes):
            
            continuum_value = np.sum(fluxes)
            continuum_values.append(continuum_value)
            
        for waves, fluxes in zip(feature_waves, feature_fluxes):
            
            feature_value = np.sum(fluxes)
            feature_values.append(feature_value)
        
            
    
    #calculate the index
    if 'ratio' in method:
        
        index = feature_values[0] / continuum_values[0]
        
        return index
        
    elif 'sumnum' in method:
        
        numerator = np.sum(feature_values)
        index = numerator / continuum_values[0]
        
        return index
        
    elif 'sumdenom' in method:
        
        denominator = np.sum(continuum_values)
        index = feature_values[0] / denominator
        
        return index
    
    elif 'avnum' in method:
        
        numerator = np.mean(feature_values)
        index = numerator / continuum_values[0]
        
        return index
    
    elif 'avdenom' in method:
        
        denominator = np.mean(continuum_values)
        index = feature_values[0] / denominator
        
        return index
    
    elif 'sumnum_twicedenom' in method:
        
        denominator = 2 * continuum_values[0]
        
        numerator = np.sum(feature_values)
        
        index = numerator / denominator
        
        return index
    
    else:
        return np.nan
            
            

#measure a set of indices from a  reference
def measure_index_set(spec, ref='lepine2007', index_info=None, plot=False, no_trials=1000, diag_path=None):
    
    index_info = None
    
    for reference in list(kastclassify.defs.index_dict):
        
        if ref in reference:
            index_info = kastclassify.defs.index_dict[reference]
            
    
    names = list(index_info[indices_keyword])
    
    result = {}
    
    waves = spec.wave.value
    flx = spec.flux.value
    unc = spec.unc.value
    
    for n in names:

        
        #measuring the value here

        index = measure_index(waves, flx, unc,\
                 index_info[indices_keyword][n][continuum_keyword],\
                 index_info[indices_keyword][n][feature_keyword],\
                 index_info[indices_keyword][n][method_keyword],\
                 index_info[indices_keyword][n][sample_keyword], plot=plot, name=n, diag_path=diag_path)
            
        
        #doing MC here for error
        
        #if we get the no flux string, pass that back as the result
        #and don't do MC
        if isinstance(index, str):
            result[n] = (index, index)
            continue
        
        index_samples = []
        
        
        for i in range(no_trials):
            flxdist = np.random.normal(flx, unc)
            
            index_sample = measure_index(waves, flxdist, unc,\
             index_info[indices_keyword][n][continuum_keyword],\
             index_info[indices_keyword][n][feature_keyword],\
             index_info[indices_keyword][n][method_keyword],\
             index_info[indices_keyword][n][sample_keyword], name=n)
            
            if isinstance(index_sample, str):
                print('Warning! Flux may not be sufficiently above noise level!')
                continue
            
            index_samples.append(index_sample)
            

        err = np.std(index_samples)
        
        result[n] = (index, err)
        
        
    return result


#calculations using indices

#lepine solar TiO5 (Lepine 2007)
def measure_solar_TiO5(index_dict):
    
    
    indices = index_dict
    #unpacking values
    
    CaH3 = indices['CaH3'][0]
    CaH3_err = indices['CaH3'][1]
    
    CaH2 = indices['CaH2'][0]
    CaH2_err = indices['CaH2'][1]
    
    #getting the measurement here
    solar_TiO5_measurement = solar_TiO5(CaH2, CaH3)
    
    #doing MC here
    no_trials = 1000
    
    CaH3_samps = np.random.normal(CaH3, CaH3_err, no_trials)
    
    CaH2_samps = np.random.normal(CaH2, CaH2_err, no_trials)
    
    solar_TiO5_samps = solar_TiO5(CaH2_samps, CaH3_samps)
    

    
    solar_TiO5_err = np.std(solar_TiO5_samps)
    
    return solar_TiO5_measurement, solar_TiO5_err


def measure_lepine_zeta(index_dict):
    
    indices = index_dict
    
    
    #unpacking values
    
    
    TiO5 = indices['TiO5'][0]
    TiO5_err = indices['TiO5'][1]
    
    solar_TiO5, solar_TiO5_err = measure_solar_TiO5(indices)
    
    #getting the measurement
    lepine_zeta_measurement = lepine_zeta(TiO5, solar_TiO5)
    
    #doing MC for error here
    no_trials = 1000
    
    TiO5_samps = np.random.normal(TiO5, TiO5_err, no_trials)
    
    solar_TiO5_samps = np.random.normal(solar_TiO5, solar_TiO5_err, no_trials)
    
    lepine_zeta_samps = lepine_zeta(TiO5_samps, solar_TiO5_samps)
    
    
    lepine_zeta_err = np.std(lepine_zeta_samps)
    
    return lepine_zeta_measurement, lepine_zeta_err
    

def solar_TiO5(CaH2, CaH3):
    
    CaH = CaH2 + CaH3
    
    solar_TiO5 = -0.164 * CaH**3 + 0.670 * CaH**2 - 0.118 * CaH - 0.050
    
    return solar_TiO5

def lepine_zeta(TiO5, solar_TiO5):
    
    zeta = (1 - TiO5)/(1 - solar_TiO5)
    
    return zeta

    
def determine_metallicity_class(zeta):
    
    
    metallicity_class = ''
    
    
    if zeta > 0.825:
        
        metallicity_class = 'Dwarf'
        
    elif ((zeta > 0.500) and (zeta < 0.825)):
        
        metallicity_class = 'sd'
        
    elif ((zeta > 0.200) and (zeta < 0.500)):
        
        metallicity_class = 'esd'
        
    elif zeta < 0.200:
        
        metallicity_class = 'usd'
        
    else:
        metallicity_class = 'Unknown'
          
    return metallicity_class

#TODO

def classify_by_index(spec, index_dict=None, no_trials=100, plot=False):
    
    if index_dict is not None:
        lepine_indices = index_dict['lepine2007']
        gizis_indices = index_dict['gizis']
        
    else:  
        lepine_indices = measure_index_set(spec, ref='lepine2007', no_trials=no_trials, plot=plot)

        gizis_indices = measure_index_set(spec, ref='gizis', no_trials=no_trials, plot=plot)
    
    zeta, zeta_err = measure_lepine_zeta(lepine_indices)
    
    metallicity_class = determine_metallicity_class(zeta)

    indices = gizis_indices

    TiO5_index = indices['TiO5'][0]
    TiO5_err = indices['TiO5'][1]

    CaH2_index = indices['CaH2'][0]
    CaH2_err = indices['CaH2'][1]

    CaH3_index = indices['CaH3'][0]
    CaH3_err = indices['CaH3'][1]

    TiO5_spt, CaH2_spt, CaH3_spt = gizis_spt(TiO5_index, CaH2_index, CaH3_index)

    #doing MC

    TiO5_samps = np.random.normal(TiO5_index, TiO5_err, no_trials)
    CaH2_samps = np.random.normal(CaH2_index, CaH2_err, no_trials)
    CaH3_samps = np.random.normal(CaH3_index, CaH3_err, no_trials)

    TiO5_spt_samps, CaH2_spt_samps, CaH3_spt_samps = gizis_spt(TiO5_samps, CaH2_samps, CaH3_samps)

    TiO5_spt_err = np.std(TiO5_spt_samps)
    CaH2_spt_err = np.std(CaH2_spt_samps)
    CaH3_spt_err = np.std(CaH3_spt_samps)

    result = {}

    result['gizis_TiO5_spt'] = TiO5_spt
    result['gizis_TiO5_spt_err'] = TiO5_spt_err
    result['gizis_CaH2_spt'] = CaH2_spt
    result['gizis_CaH2_spt_err'] = CaH2_spt_err
    result['gizis_CaH3_spt'] = CaH3_spt
    result['gizis_CaH3_spt_err'] = CaH3_spt_err


    indices = lepine_indices
    #unpacking values

    CaH3 = indices['CaH3'][0]
    CaH3_err = indices['CaH3'][1]

    CaH2 = indices['CaH2'][0]
    CaH2_err = indices['CaH2'][1]

    #getting the value here
    lepine_spt_measurement = lepine_spt(CaH2, CaH3)

    #doing MC here

    CaH3_samps = np.random.normal(CaH3, CaH3_err, no_trials)

    CaH2_samps = np.random.normal(CaH2, CaH2_err, no_trials)

    lepine_spt_samps = lepine_spt(CaH2_samps, CaH3_samps)


    lepine_spt_err = np.std(lepine_spt_samps)

    result['lepine_spt'] = lepine_spt_measurement
    result['lepine_spt_err'] = lepine_spt_err

    return result, (zeta, zeta_err), metallicity_class    

    

def lepine_spt(CaH2, CaH3):
    CaH = CaH2 + CaH3
    
    spt = 1.4 * CaH**2 - 10. * CaH + 12.4
    
    return spt

def gizis_spt(TiO5, CaH2, CaH3):
    
    TiO5_spt = -9.64 * TiO5 + 7.76
    CaH2_spt = 7.91 * CaH2**2 - 20.63 * CaH2 + 10.71
    CaH3_spt = -18. * CaH3 + 15.8
    
    return TiO5_spt, CaH2_spt, CaH3_spt


#EW functions

def measure_EW_simple(waves, fluxes, unc, line, wave_range, feat_type, line_width=[0.,0.], plot=False, file=None, verbose=False, name='Unknown Line', diag_path=None):


    feat_type = feat_type.lower()
    
    #check to see if there's flux
    if np.all(np.abs(fluxes) <= 2 * unc):
            print('Warning! No flux in region for {}!'.format(name))
            
            res = 'No Flux'
            
            return res, res
    
    feature_waves = waves * u.AA
    
    feature_flux = fluxes * u.erg / (u.AA * u.s * (u.cm)**2)
    
    uncertainties = unc
    
    feature_unc = astropy.nddata.StdDevUncertainty( unc * u.erg / (u.AA * u.s * (u.cm)**2))
    
    feature_spec = specutils.Spectrum1D(spectral_axis=feature_waves, flux=feature_flux, uncertainty=feature_unc)
    
    continuum_fit = specutils.fitting.fit_generic_continuum(feature_spec)(feature_spec.spectral_axis)
    

    #looping to refine the continuum fit
    
    #initial mask around the feature
    mask = np.ones(len(feature_flux), dtype='bool')
    mask[ (feature_waves.value > line_width[0]) & (feature_waves.value < line_width[1])] = False
    
    subtr_spec = feature_spec - continuum_fit

    
    for i in range(5):
        
        masked_spec = specutils.Spectrum1D(spectral_axis=feature_waves[mask], flux=feature_flux[mask], uncertainty=feature_unc[mask])
        
        continuum_model = specutils.fitting.fit_generic_continuum(masked_spec)
        continuum_fit = continuum_model(feature_spec.spectral_axis)
        
        subtr_spec = feature_spec - continuum_fit
        
        
        mask = np.abs(subtr_spec.flux.value) < 3 * uncertainties
        
        if verbose:
            plt.plot(feature_waves, continuum_fit)
            plt.plot(feature_waves, feature_spec.flux)
            plt.title('refine {}'.format(i))
            plt.show()
            
     
    continuum_fit = continuum_model(feature_spec.spectral_axis)
    
    subtr_spec = feature_spec - continuum_fit
        
    if plot or (diag_path is not None):
        
        plt.figure(figsize=[12,7])
        
        plt.plot(feature_spec.wavelength, feature_spec.flux)
        plt.plot(feature_spec.wavelength, continuum_fit)
        
        #plt.xlabel('Wavelength')
        #plt.ylabel('Flux Density')
        
        plt.title('Continuum fit of {}'.format(name))
        
        if isinstance(file, str):
            
            plt.savefig('kastclassify_EW_continuum_' + file)
            
        if diag_path is not None:
            
            plt.savefig(diag_path + '{}_EW_continuum_fit.png'.format(name))
        
        if plot:
            
            plt.show()
            
        plt.close()
    
    
    #finding line center
    
    near_line_waves = feature_waves[(feature_waves.value >= (line - 10)) & (feature_waves.value <= (line + 10))]
    near_line_fluxes = subtr_spec.flux.value[(feature_waves.value >= (line - 10)) & (feature_waves.value <= (line + 10))]
    

    if feat_type == 'absorption':
        
        line = near_line_waves.value[near_line_fluxes == np.amin(near_line_fluxes)]
        
    elif feat_type == 'emission':
        
        line == near_line_waves.value[near_line_fluxes == np.amax(near_line_fluxes)]
        
    else:
        print('WARNING! Feature type {} is unknown'.format(feat_type))
        
        res = 'Unknown Feature Type'
        
        return res, res
    
    
    
    #create continuum normalized spectrum and continuum subtracted spectrum
    norm_spec = feature_spec / continuum_fit
    
    #re-centering using subtracted spectrum

    
    #find all features in this spectrum
    line_info = specutils.fitting.find_lines_threshold(subtr_spec)

    
    if len(line_info) == 0:
        print('Warning! No features detected for {}!'.format(name))
        
        res = 'No Feature Detected'
        
        return res, res
    
    try:
        #the actual line center is the closest feature to our expected center
        nearest = np.argmin(np.abs(line_info[ line_info['line_type'] == feat_type]['line_center'].value - line))

    except:
        print('Warning! No features detected for {}!'.format(name))

        res = 'No Feature Detected'

        return res, res
        
    line_center = line_info[ line_info['line_type'] == feat_type][nearest]['line_center'].value

    
    
    #if the line width is not defined, calculate it
    #calculated as 1.5 * FWHM
    if line_width[0] == 0.:
        
        
        near_line_waves = feature_waves[(feature_waves.value >= (line - 15)) & (feature_waves.value <= (line + 15))].value
        near_line_fluxes = subtr_spec.flux.value[(feature_waves.value >= (line - 15)) & (feature_waves.value <= (line + 15))]
        
        if feat_type == 'absorption':
            
            near_line_fluxes = np.abs(near_line_fluxes) 
        
        elif feat_type == 'emission':
            
            pass
            
        else:
            print('WARNING! Feature type {} is unknown'.format(feat_type))

            res = 'Unknown Feature Type'

            return res, res
        
        max_flx = np.amax(near_line_fluxes)
        
        half_max = max_flx / 2
        
        outside = near_line_fluxes < half_max
        
        try:
            lower_bound = max( near_line_waves[(near_line_waves < line_center) & outside])
            upper_bound = min( near_line_waves[(near_line_waves > line_center) & outside])
        except:
            print('Warning! No features detected for {}!'.format(name))
        
            res = 'No Feature Detected'

            return res, res
        
        FWHM = upper_bound - lower_bound
        
        line_width[0] = 1.5 * FWHM
        line_width[1] = 1.5 * FWHM
        
    
    
    #pick out the region containing the feature
    region = specutils.SpectralRegion( (line_center - line_width[0]) * u.AA, (line_center + line_width[1]) * u.AA )

    
    #calculate the width
    EW = specutils.analysis.equivalent_width(norm_spec, regions=region)

    
    if plot or (diag_path is not None):
        
        fig, ax = plt.subplots(1, figsize=[12,7])
        
        ax.plot(norm_spec.wavelength, norm_spec.flux, label=EW)
        
        ax.axvline(line_center - line_width[0])
        ax.axvline(line_center + line_width[1])
        
        rect = matplotlib.patches.Rectangle( (line_center - EW.value/2, 0), EW.value, 1, hatch='x', fill=False)
        
        ax.add_patch(rect)
        
        ax.set_xlabel('Wavelength ({})'.format(u.AA))
        ax.set_ylabel('Normalized Flux Density')
        
        ax.set_title('Equivalent Width of {}'.format(name))

        fig.legend()
        
        if isinstance(file, str):
            
            fig.savefig('kastclassify_EW_width_' + file)
            
        if diag_path is not None:
            
            fig.savefig(diag_path + '{}_EW_measurement.png'.format(name))
            
        if plot:
            
            fig.show()
        
        plt.close(fig)
        
    return EW, line_center

def measure_EW_element(spec, ref, no_trials=100, plot=False, file=None, verbose=False, diag_path=None):
    
    line_info = None
    
    for el in list(kastclassify.defs.EW_dict):
        
        if ref == el:
            
            line_info = kastclassify.defs.EW_dict[el]
            
    
    result = {}
    
    wave = spec.wave.value
    flux = spec.flux.value
    uncs = spec.unc.value
    
    
    no_samples = 10000
    
    finterp = interp1d(wave, flux, bounds_error=False, fill_value=0.)
    uinterp = interp1d(wave, uncs, bounds_error=False, fill_value=0.)
    
    waves = np.linspace(line_info[wave_range_keyword][0], line_info[wave_range_keyword][1], no_samples)
    
    fluxes = finterp(waves)
    unc = uinterp(waves)
    
    EW_info = measure_EW_simple(waves, fluxes, unc,\
                           line_info[line_keyword], line_info[wave_range_keyword],\
                           line_info[feat_type_keyword], plot=plot, verbose=verbose,\
                           name=line_info[name_keyword], file=file, diag_path=diag_path)
    
    if plot:
        plt.pause(1)
    
    EW = EW_info[0]
    line_center = EW_info[1]

    #if we get any error string, pass that back as the result
    #and don't do MC
    if isinstance(EW, str):
        result['EW'] = EW
        result['EW_err'] = EW
        result['line_center'] = EW
        result['line_center_err'] = EW
        return result
    
    #otherwise, break out the value
    else:
        EW = EW.value
    
    #doing MC here
    EW_samples = []
    line_center_samples = []
    for i in range(no_trials):
        
        
        flxdist = np.random.normal(fluxes, unc)
        
        EW_trial = measure_EW_simple(waves, flxdist, unc,\
                           line_info[line_keyword], line_info[wave_range_keyword],\
                           line_info[feat_type_keyword], plot=False, name=line_info[name_keyword])
        
        EW_sample = EW_trial[0].value
        line_center_sample = EW_trial[1]
        
        EW_samples.append(EW_sample)
        line_center_samples.append(line_center_sample)

        
    err_EW = np.std(EW_samples)
    err_line_center = np.std(line_center_samples)
    
    result['EW'] = EW
    result['EW_err'] = err_EW
    result['line_center'] = line_center
    result['line_center_err'] = err_line_center
    
    return result


#TODO
#interpolate for chi factors
def calculate_L_Lbol(Halpha_EW, Halpha_EW_err, spt):
    
    chi_interp = interp1d(kastclassify.defs.chi_factor_types, kastclassify.defs.chi_factor_values, bounds_error=False, fill_value=np.nan)

    chi_err_interp = interp1d(kastclssify.defs.chi_factor_types, kastclassify.defs.chi_factor_values_err, bounds_error=False, fill_value=np.nan)
    
    stype = typeToNum(spt)
    
    chi = chi_interp(stype)
    
    chi_err = chi_err_interp(stype)

    Lratio = chi * Halpha_EW
    Lratio_err = chi_err * Halpha_EW_err
    
    return Lratio, Lratio_err


####################################
### END ANALYSIS FUNCTIONS #########
####################################


####################################
### BEGIN WORKFLOW FUNCTIONS #######
####################################

#calculating information for one spectrum
def calc_spec_info(filepath, plot=False, no_index_trials=100, no_EW_trials=30, diag_path=None):
    
    
    t1 = datetime.datetime.now()
    
    #collecting header info
    print('Reading in {}'.format(filepath))
    
    spec = kastredux.readSpectrum(filepath)
    spec.name = fits.getheader(filepath)['OBJECT']
    
    ra = fits.getheader(filepath)['RA']
    
    dec = fits.getheader(filepath)['DEC']
    
    date = fits.getheader(filepath)['DATE']
    
    filename = os.path.basename(filepath)
    
    if plot:
        spec.plot()
        plt.show()
    if diag_path is not None:
        
        spec.plot()
        plt.savefig(diag_path + spec.name + '_spec.png')
        
    
    print('Classifying by standard comparison...')
    spt, stand_name, spt_chi2, spt_scale = classify_by_standard(spec, ref='all', plot=plot, diag_path=diag_path)
    
    print('Done!')
    
    print('Measuring indices...')
    indices = {}
    
    for ref in list(kastclassify_index_dict):
        
        print('Mesuring index set {}...'.format(ref))
        
        indices[ref] = measure_index_set(spec, ref=ref, no_trials = no_index_trials)
        
        print('Done!')
    
    print('Done!')
    
    print('Measuring Equivalent Widths...')
    EWs = {}
    
    for ref in list(kastclassify_EW_dict):
        
        print('Measuring EW for {}...'.format(ref))
        
        EWs[ref] = measure_EW_element(spec, ref=ref, plot=plot, no_trials=no_EW_trials, diag_path=diag_path)
        
        print('Done!')
     
    
    print('Done!')
    #calculating with the measurements
    
    print('Calculating spectral types by index...')
    index_spt, zeta_info, metallicity_class = classify_by_index(spec, index_dict=indices)
    
    print('Done!')
    
    print('Calculating L_Halpha/L_bol...')
    
    Lratio, Lratio_err = calculate_L_Lbol(EWs['Halpha']['EW'], EWs['Halpha']['EW_err'], spt)
    
    print('Done!')
    
    plt.close('all')
    #unpacking results
    
    print('Packaging results...')
    
    results = {}
    
    results['name'] = spec.name
    results['filename'] = filename
    results['RA'] = ra
    results['DEC'] = dec
    results['OBSDATE'] = date
    
    results['spt'] = spt
    results['subtype'] = typeToNum(spt) - 10
    results['spt_chi2'] = spt_chi2
    results['spt_scale'] = spt_scale
    
    
    for ref in list(indices):
        
        for name in list(indices[ref]):
            
            results[ref + '_' + name] = indices[ref][name][0]
            results[ref + '_' + name + '_err'] = indices[ref][name][1]
            
    for ref in list(EWs):
        
        for name in list(EWs[ref]):
            
            results[ref + '_' + name + '_(Angstrom)'] = EWs[ref][name]
            
    for ref in list(index_spt):
        
        results[ref] = index_spt[ref]
        
    results['zeta'] = zeta_info[0]
    
    results['zeta_err'] = zeta_info[1]
    
    results['metallicity_class'] = metallicity_class
    
    results['L_Halpha/L_bol'] = Lratio
    
    results['L_Halpha/L_bol_err'] = Lratio_err
    
    print('Done!')
    
    t2 = datetime.datetime.now()
    dt = t2 - t1
    
    print('All done in {}! All done at {}'.format(dt, t2))
    
    return results


#calculate info for many spectra
def classify(file_list, outfile, infile=None, no_index_trials=100, no_EW_trials=30, plot=False, diag_path=None):
    
    filepaths = copy.deepcopy(file_list)
    
    final_results = pd.DataFrame()
    
    if diag_path is not None:
        
        if os.path.exists(diag_path):
            print('Path {} already exists'.format(diag_path))
        else:
            
            try:
                os.makedirs(diag_path)
            except OSError:
                print('Failed to create directory {}'.format(diag_path))
                return None
        
        print('Saving diagnostic plots to {}'.format(diag_path))
    
    if infile is not None:
        
        final_results = pd.read_csv(infile)
        
        known_filenames = final_results['filename'].values
        
        known_removed = [path for path in filepaths if not any(name in path for name in known_filenames)]
        
        filepaths = known_removed
    
    for filepath in filepaths:
        
        if diag_path is not None:
            diag_dir = diag_path + '/{}/'.format(os.path.splitext(os.path.basename(filepath))[0])
            if os.path.exists(diag_dir):
                print('Path {} already exists'.format(diag_dir))
            else:

                try:
                    os.makedirs(diag_dir)
                except OSError:
                    print('Failed to create directory {}'.format(diag_dir))
                    diag_dir = None
        else:
            diag_dir = diag_path
        
        spec_info = calc_spec_info(filepath, plot=plot, no_index_trials=no_index_trials, no_EW_trials=no_EW_trials, diag_path=diag_dir)
        
        s = pd.Series(spec_info)
        
        s = s.fillna('Invalid Data!')
        
        df = s.to_frame().transpose()
        
        final_results = final_results.append(df)
        
        final_results.to_csv(outfile, index=False)
        
##################################
### END WORKFLOW FUNCTIONS #######
##################################