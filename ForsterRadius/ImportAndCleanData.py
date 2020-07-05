# %%
import requests
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import warnings
from fuzzywuzzy import process, fuzz

PI = np.pi
AVOGADRO = 6.022e23 # molecules / mole
MIN_WAVELENGTH = 100.0
MAX_WAVELENGTH = 901.0
STEP_SIZE = 1
WAVELENGTHS_USED = list(np.arange(MIN_WAVELENGTH, MAX_WAVELENGTH, STEP_SIZE, dtype=float))

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

# Takes [[wavelengths],[values]] and returns an interpolated Series; drops duplicate entries.
def interpolate(spectrum):
    orig_idx = spectrum[0]
    combined_idx = sorted(list(set().union(WAVELENGTHS_USED, orig_idx))) # set structure is unordered
    df = pd.Series(spectrum[1], index=orig_idx)
    df = df.loc[~df.index.duplicated()]
    df = df.reindex(index=combined_idx) # don't fill with 0s yet, that's for interpolate to do
    df = df.interpolate(method='index', limit_area='inside')
    df = df.reindex(WAVELENGTHS_USED)
    df = df.fillna(0)
    return df

# Returns new DataFrame with old columns merged into one new col.
# If override, the first value will be kept; else, values are concatenated as strings.
def merge_columns(df, col_labels, new_col_name, override=False, separator=None):
    merged_col = df.loc[:,col_labels[0]]
    if override:
        for col in col_labels[1:]:
            merged_col = merged_col.combine_first(df.loc[:,col])
    else:
        merged_col = merged_col.apply(str).str.cat(df.loc[:,col_labels[1:]].applymap(str),sep=separator,na_rep='')
    df = df.drop(columns=col_labels)
    df.loc[:,new_col_name] = merged_col
    return df

# Stores data into FPBase Spectra.json
def import_fpbase_spectra():
    r = requests.get('https://www.fpbase.org/api/proteins/spectra/?format=json').json()
    with open('FPBase Spectra.json', 'w', encoding='utf-8') as f:
        json.dump(r, f, ensure_ascii=False, indent=4)
    return True

# Retrieves data from FPBase Spectra.json
def get_fpbase_spectra():
    with open('FPBase Spectra.json', 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

# Verifies that each "state" value in the FPBase spectra contains either 'ex' or 'em'
def verify_separation(): 
    data = get_fpbase_spectra()
    for fluor in data:
        for spectrum in fluor["spectra"]:
            state = spectrum["state"]
            if (state.find("em") >= 0) is (state.find("ex") >= 0): # if both or neither of 'ex', 'em' are found
                print (fluor["name"], state)
                return False
    return True

def fpbase_spectra_to_df():
    data = get_fpbase_spectra()
    df = pd.json_normalize(data,'spectra',['name'])
    df.loc[:,'spectrum'] = df.loc[:,'data'].apply(lambda x: list(zip(*x)))
    df = df.drop(columns='data')
    df.loc[:,'role'] = np.where(df.loc[:,'state'].map(lambda x: 'ex' in x), 'Acceptor', 'Donor')
    df = df.set_index(['name','state'])
    df.loc[:,'max'] = df.loc[:,'max'].fillna(df.loc[:,'spectrum'].map(lambda x: x[0][np.argmax(x[1])]))
    
    df_s = pd.DataFrame(index=WAVELENGTHS_USED, columns=df.index)
    for i in df.index:
        df_s[i] = interpolate(df.loc[i,'spectrum'])
    df = pd.concat([df.drop(columns='spectrum'), df_s.T], axis=1)
    # TODO total
    df = df.reset_index().rename(columns={'state': 'Other names'})
    df = df.rename(index=lambda x: x + 2000)
    df = df.reindex(columns=['name', 'role', 'Other names', 'qy', 'ec', 'max', *WAVELENGTHS_USED])
    is_donor = df['role'].map(lambda x: x == 'Donor') & df['qy'] > 0
    donors = df[is_donor]
    is_acceptor = df['role'].map(lambda x: x == 'Acceptor') & df['ec'] > 0
    acceptors = df[is_acceptor]
    donors['total'] = donors.index.to_series().map(lambda i: np.trapz(list(donors.loc[i,np.arange(100,901)]), WAVELENGTHS_USED))
    df['total'] = donors['total']

    df.to_excel('fp_data.xlsx')
    donors.to_excel('fp_donors.xlsx')
    acceptors.to_excel('fp_acceptors.xlsx')

    return df, donors, acceptors

# NOTE Important changes to Utzinger data that were MADE MANUALLY:
    # DOWNLOAD http://spectra.arizona.edu/utzinger_spectra.sql, save as 'Utzinger Spectra.sql'

    # INSERT a backslash before all (14) double quotes. " -> \"
    # REPLACE all (9) instances of four single quotes with an escaped double quote. '''' -> \"
    # REPLACE all instances of a number followed by two double quotes with the number followed by an escaped single quote. 2'' -> 2\'
    #   2'' (1), 3'' (5), 4'' (3), 4'' (5), 7'' (1)

    # Lines 43-15748 are 'attributes', records in the form (id, attribute, value):
    #   43      (1, 'Long Dye Name', '1-anilinonaphthalene-8-sulfonic acid in EtOH'),
    #   ...
    #   15748   (1112, 'Source', 'DRBIO');

    # Lines 15770-17538 are 'items' of the form (`id`, `name`, `type`, `subtype`, `published`):
    #   15770   (1, '1-ANS', 'dye', 'AB', 'Y'),
    #   ...
    #   17538   (1831, 'SP 0532 - BSP01-532R', 'filter', 'SP', 'Y');

    # Lines 17557-1104997 are spectra, flattened into the form (`id`, `wavelength`, `value`):
    #   17557 (1, '250.00', 1.1315900000),
    #   ...
    #   1104997 (1831, '1210.00', 0.0002410000);

# Stores data into Utzinger Spectra.sql
def get_ut_spectra():
    with open('Utzinger Spectra.sql', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    return lines

def ut_spectra_to_df():
 
    lines = get_ut_spectra()    
    # I manually checked that any line that contains NULL is irrelevant (none of them contain data about dyes)
    # so I just ignore any lines with NULL in any value
    attributes = [eval(x.strip()[:-1]) for x in lines[42:15748] if (x[0] == '(' and 'NULL' not in x)]
    items = [eval(x.strip()[:-1]) for x in lines[15769:17538] if (x[0] == '(' and 'NULL' not in x)]
    spectra = [eval(x.strip()[:-1]) for x in lines[17556:1104997] if (x[0] == '(' and 'NULL' not in x)]
    del lines

    df_a = pd.DataFrame(attributes).rename(columns={0:'id',1:'attribute',2:'value'})
    df_a_pivoted = df_a.drop_duplicates(subset=['id','attribute']).pivot(index='id',columns='attribute',values='value')
    df_i = pd.DataFrame(items).set_index(0).drop(4, axis=1).rename(columns={1:'name',2:'type',3:'role'})
    df = pd.concat([df_i,df_a_pivoted],axis=1).applymap(lambda x: x.strip() if type(x) == str else x) # strip() all the string elements
    
    # Incorporate the spectral data
    spectra.sort(key=lambda x: (x[0],float(x[1])))
    spectra_arr = np.array(list(zip(*spectra)),dtype=float).T
    del spectra
    split_by_id = np.split(spectra_arr, np.where(np.diff(spectra_arr[:,0]))[0]+1) # split whenever the id changes; this is why we needed to sort spectra by id
    spectra_dict = {x[0][0]: x.T[1:] for x in split_by_id}
    df_s = pd.DataFrame(index=WAVELENGTHS_USED, columns=df.index)
    for key in spectra_dict:
         df_s.loc[:,key] = interpolate(spectra_dict[key])
    df_s = df_s.drop(columns=1271)
    df = pd.concat([df, df_s.T],axis=1)
    
    # DROP all the completely useless columns and rows
    # we only want dyes. not filters, lights, mirrors, dyecats, filtersets, and ESPECIALLY not NaN's >:(
    type_is_dye = df.loc[:,'type'].map(lambda x: x == 'dye')
    role_in_emabex = df.loc[:,'role'].apply(str).str.upper().map(lambda x: x in 'EM,AB,EX' and len(x)>1)
    df = df[type_is_dye & role_in_emabex]
    has_spectrum = df.index.to_series().map(lambda x: x in spectra_dict.keys())
    df = df[has_spectrum]

    df.loc[:,'role'] = df.loc[:,'role'].map({'EM': 'Donor', 'AB': 'Acceptor', 'EX': 'Acceptor'})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df.loc[:,'acronym or common name'] = df.loc[:,'acronym or common name'].fillna('',inplace=True)
    other_names = ['acronym or common name','alternative name', 'Acronym','Alternative','acronyms and alternative names','acronyms or common names','alternate name','Full name', 'Long Dye Name', 'Long dye name', 'Dye name']
    df = merge_columns(df, other_names, 'Other names', separator=';')
    df = df.drop(columns=['Common name', 'Name', 'acronmym or common name', 'common name','LONG NAME','Dye','dye name'])
    df = df.drop(columns=['type'])
    df = df.drop(columns=['', '2P', '2P brightness'])
    df = df.drop(columns=['1P primary source', '2P primary source', 'Data source', 'SOURCE', 'Secondary Source', 'Secondary source', 
        'Source', 'Source for Extinction Coefficient', 'Source`', 'PRIMARY SOURCE', 'Primary Source', 'Primary source', 
        'primary source', 'reference', 'secondary source', 'source', 'Reference', 'Full Reference', 'Data Source'])
    df = df.drop(columns=['Additional Info', 'Comment', 'Comment 2nd', 'Comments', 'NOTES', 'Note', 'Notes', 'note', 'notes'])
    df = df.drop(columns=['BANDCENTER BANDPASS', 'BANDPASS', 'band edges', 'bandcenter/bandpass', 'Number bands', 'NUMBER BANDS'])
    df = df.drop(columns=['Brightness', 'Brightness (% of DsRed)', 
        'CF', 'CF260', 'CF280', 
        'Category', 'Category of dye', 
        'Company', 'Contact', 'Data entry'])
    env = ['Environment', 'environment', 'solvent', 'Solvent', 'solvent   alternative name']
    df = merge_columns(df, env, 'Environment', override=True)
    df = df.drop(columns=['pH', 'pKa'])
    df = df.drop(columns=['Ex  Coeff', 'Visible light Extinction coefficient value', 'molar extinction M^-1 cm^-1', 
        'Extinction Coeficient', 'Extinction', 'extinction coefficient or quantum yield']) # Contains only nonunique or incorrect information.
     
    ex_coeff = ['Ex. Coeff', 'Exctinction coefficient', 'Extinction Coefficient', 'Extinction Coefficient:', 
        'Molar Ex  Coeff', 'Molar Extinction', 'extinction coefficient', 'molar extinction coefficient']
    df = merge_columns(df, ex_coeff, 'ec', override=True)
    df.loc[508,'ec'] = 13400 # Filled from df['extinction coefficient or quantum yield']
    df.loc[933,'ec'] = 50000 # Filled from df['Ex  Coeff']
    df.loc[375,'ec'] = 148000 # Value taken from catalog no. D3911 at https://www.thermofisher.com/us/en/home/references/molecular-probes-the-handbook/probes-for-lipids-and-membranes/dialkylcarbocyanine-and-dialkylaminostyryl-probes.html#datatable
    df.loc[377,'ec'] = 138000 # Catalog no. D3898

    df = df.drop(columns=['FLUOR LIFETIME', 'Lifetime', 'Itensity max', 'Medium', 'Module Configurator',
        'Location', 'Field19', 'Fluorophore (typical)', 'MANUFACTURER', 'Manufacturer', 'manufacturer',
        'Have Item', 'intensity max', 'INTENSITY MAX', 'Max intensity', 'Imaging quality', 'Intensity @ 200mm (mW/cm2', 'Intensity max'])
    df = df.drop(columns=['Mol  Wt', 'Mol Wt.', 'Mol. Wt.', 'Molar Weight', 'Molar Weigth', 'mol. wt.'])
    df = df.drop(columns=['PART NAME'])
    df = df.drop(columns=['R <0> Alexa 488', 'R<0> Val. with AF488 acceptor', 'R<0> Val. with AF546 acceptor', 'R<0> Val. with AF555 acceptor', 'R<0> Val. with AF568 acceptor', 'R<0> Val. with AF594 acceptor', 'R<0> Val. with AF647 acceptor', 'R<0> Value source', 'R<0> Values', 
        'Test', 'Total Power @ 200 mm (mW)', 'URL', 'URL Quantum Yield', 'URL have item', 'Units'])
    df = df.drop(columns=['block Ti:Sa', 'category', 'category of dye', 'component of set', 'filename', 'filter lamda', 'filter pass', 'filter wavelength', 
        'fliter pass', 'full reference', 'have doc', 'have item', 'have spectr', 'have spectrum', 'how many', 'in database', 'in db', 'in system', 
        'include', 'lot #',  'max value 381 to 700\r\n', 'miscellaneous', 'position', ])
    df = df.drop(columns=['t0 5 bleach', 't0 5 bleach (min)', 't0 5 maturatin', 't0 5 maturation (hr)', 't0.5 bleach', 't0.5 bleach (min)', 't0.5 maturatin', 't0.5 maturation (hr)', 
        'thick (mm)', 'thick(mm)', 'thickness','dia (mm)'])
    df = df.drop(columns=['Spectrum Type', 'Spectrum type', 'spectrum type', 'spectrum type (1P,2p)', 'spectrum type (a,x,m)', 'spectrum type (a,x,m, 2P)', 'spectrum type (a,x,m,2p)', 'type (ex, bs, em)'])
    df = df.drop(columns=['Original Data Rescaled By', 'rescaled by (max value)', 'rescaling factor', 'Rescale 2P by:', 'Rescale Em by:', 'Rescale Ex by:', 'Rescaled 2P by:'])

    quant_yield = ['Quant  Yield', 'Quant. Yield', 'Quantum Yield']
    df = merge_columns(df, quant_yield, 'qy', override=True)
    df.loc[573,'qy'] = 0.4300 # Filled from df['extinction coefficient or quantum yield']
    # The wavelength data is WAY too messy... see test_lam.xlsx for the raw data. 
    # I decided to determine the wavelengths myself using the spectrum data.
    df = df.drop(columns=['Wavelength of Max Excitation', 'Wavelength Used To Determine Ex. Coeff.', 'ABS', 'Abs/Em maxima', 'EX EM', 'Em max', 'Em max (nm)', 'Ex max', 
        'Max Wavelength',  'Max wavelength', 'Molar Ex Coeff. Wavelength', 'Visible Light at Ec Wavelength', 'Visible light Extinction coefficient wavelength', 'WAVELENGTH', 
        'Wavelength Used To Determine Ex  Coeff', 'Wavelength at Ec', 'Wavelength at Ec (Table)', 'Wavelength at Ec (from source)', 'Wavelength at Ec (source)', 'Wavelength at Ec (table)', 'Wavelength at Max Absorption', 
        'Wavelength at Max Emission', 'Wavelength of Max Absorption', 'Wavelength used to determine Ex. Coeff.', 'Wavelength used to measure Ext  Coefficient', 'Wavelength used to measure extinction', 
        'Wavlength used to measure extinction coefficient', 'max wavelength', 'max wavlength', 'visble light Ec wavelength', 'visible light EC wavelength', 'visible light Ec wavelength', 
        'visible light Ec wavlength', 'visible light at Ec wavelength', 'visible light at Ec wavelenth', 'visible light wavelength', 'visiblle light Ec wavelength', 'wavelength at EC (from table)', 'wavelength at Ec', 'wavelength at Ec (from source', 
        'wavelength at Ec (from source)', 'wavelength at Ec (from table)', 'wavelength at Ec or Em (table below)', 'wavlelngth at Ec (from source)', 'wavlength at Ec', 'wavlength at Ec (from source)'])
    
    df = df.reindex(columns=['name', 'role', 'Other names', 'qy', 'ec', 'max', *WAVELENGTHS_USED])
    # only accept donors/acceptors with valid qy/ec values respectively
    df['max'] = df.index.to_series().map(lambda i: spectra_dict[i][0][np.argmax(spectra_dict[i][1])])    
    is_donor = df['role'].map(lambda x: x == 'Donor') & df['qy'].map(lambda x: isfloat(x) and float(x) > 0)
    donors = df[is_donor]
    is_acceptor = df['role'].map(lambda x: x == 'Acceptor') & df['ec'].map(lambda x: isfloat(x) and float(x) > 0)
    acceptors = df[is_acceptor]
    df = df[is_donor | is_acceptor]
    
    donors['qy'] = donors['qy'].apply(float)
    acceptors['ec'] = acceptors['ec'].apply(float)
    donors['total'] = donors.index.to_series().map(lambda i: np.trapz(list(donors.loc[i,np.arange(100,901)]), WAVELENGTHS_USED))
    df['total'] = donors['total']

    df.to_excel('ut_data.xlsx')
    donors.to_excel('ut_donors.xlsx')
    acceptors.to_excel('ut_acceptors.xlsx')
    return df, donors, acceptors

# donor and acceptor should be records from the dataframes generated in the main methods.
def calc_forster_radius(donor, acceptor, quant_yield=None, kappa2 = 2./3, refractive_index = 1.4): 
    donor_em = donor.loc[np.arange(100,901)].div(donor.total) 
    # nm^-1
    acceptor_ex = acceptor.loc[np.arange(100,901)].mul(acceptor.ec / max(acceptor.loc[np.arange(100,901)].values)) 
    # M^-1 cm^-1 = mol^-1 dm^3 cm^-1
    if quant_yield is None:
        quant_yield = donor.qy
    # unitless

    coeff = 9. * np.log(10) * kappa2 * float(quant_yield) / (128. * PI**5 * AVOGADRO * refractive_index**4) 
    # mol

    integrand = [donor_em.loc[lam] * acceptor_ex.loc[lam] * lam ** 4 for lam in WAVELENGTHS_USED]
    integral = si.simps(integrand, WAVELENGTHS_USED) * STEP_SIZE 
    # nm^-1 mol^-1 dm^3 cm^-1 nm^4 nm = mol^-1 nm^6 10^17
    # dm = nm * 10^8
    # cm = nm * 10^7
    
    return (1e17 * coeff * integral)**(1./6) # nm

# df: sub-DataFrame of data
def plot(df, ids, integer_index = False):
    plt.ion()
    if integer_index:
        for i in ids:
            plt.plot(WAVELENGTHS_USED, df.iloc[i].loc[WAVELENGTHS_USED],
            label=df.iloc[i].loc['name'] + ' ' + ('EM' if df.iloc[i].loc['role'] == 'Donor' else 'EX'))
    else:
        for i in ids:
            plt.plot(WAVELENGTHS_USED, df.loc[i,WAVELENGTHS_USED],
            label=df.loc[i,'name'] + ' ' + ('EM' if df.loc[i,'role'] == 'Donor' else 'EX'))
    plt.legend()
    plt.show()
# %% Main methods. Cleans data into pandas DataFrame,
#    combines the donors and acceptors from each data set
fp, fp_donors,fp_acceptors = fpbase_spectra_to_df()
ut, ut_donors,ut_acceptors = ut_spectra_to_df()
data = pd.concat([ut,fp])
donors = pd.concat([ut_donors, fp_donors])
acceptors = pd.concat([ut_acceptors, fp_acceptors])
# %% Calculates Forster Radius for each donor/acceptor pair.
# WARNING takes a long time (hours)

fr = pd.DataFrame(index=[list(donors.index),list(donors.name)],
    columns=[list(acceptors.index), list(acceptors.name)]).rename_axis(['id','name']).rename_axis(['id','name'],axis=1)
# could have used MultiIndex.from_product, then pivot/stack after calculating values.

fr = fr.apply(lambda a: a.index.to_series().map(lambda d: calc_forster_radius(donors.loc[d[0]], acceptors.loc[a.name[0]])))
fr.to_excel('Forster Radii.xlsx')
# %%
# Assumes df has columns 'name' and 'Other names'
def search_by_name(df, query, results=10):
    names_dict = df['name'].str.lower().to_dict()
    other_dict = df['Other names'].str.lower().to_dict()
    name_matches = np.array(process.extract(query.lower(),names_dict,limit=results,scorer=fuzz.ratio))
    other_matches = np.array(process.extract(query.lower(),other_dict,limit=results,scorer=fuzz.partial_ratio))
    other_matches[:,1] = other_matches[:,1].astype(float)/2.
    matches = np.concatenate([name_matches, other_matches])
    unique_matches=matches[np.unique(matches[:,2],return_index=True)[1],]
    sorted_unique_matches = unique_matches[unique_matches[:,1].astype(float).argsort()][:-(results+1):-1]
    match_ids = sorted_unique_matches[:,2]
    return df.loc[match_ids.astype(float),:].drop(columns=WAVELENGTHS_USED)


# given a network of fluorophores with known rate constants, we need a mathematical frmaework to 
# predict interactions. we want to build this up for more complex networks, using more complex analysis or simulation.

# expanding functionality using photoswitching

# algorithm for weight training given a FRET system