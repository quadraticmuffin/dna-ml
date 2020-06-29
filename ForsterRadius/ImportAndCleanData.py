import requests
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si

def import_fpbase_spectra():
    r = requests.get('https://www.fpbase.org/api/proteins/spectra/?format=json').json()
    with open('FPBase Spectra.json', 'w', encoding='utf-8') as f:
        json.dump(r, f, ensure_ascii=False, indent=4)
    return True
import_fpbase_spectra()

def get_fpbase_spectra():
    with open('FPBase Spectra.json', 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

MIN_WAVELENGTH = 100.0
MAX_WAVELENGTH = 901.0
STEP_SIZE = 1
WAVELENGTHS_USED = list(np.arange(MIN_WAVELENGTH, MAX_WAVELENGTH, STEP_SIZE, dtype=float))
def interpolate(spectrum): # takes a [[wavelengths],[values]] and returns a Series returns the same format
    orig_idx = spectrum[1]
    combined_idx = sorted(list(set().union(WAVELENGTHS_USED, orig_idx))) # set structure is unordered
    df = pd.Series(orig_idx, index=spectrum[0])
    df = df.reindex(combined_idx) # don't fill with 0s yet, that's for interpolate to do
    df = df.interpolate(method='index', limit_area='inside')
    df = df.reindex(WAVELENGTHS_USED)
    df = df.fillna(0)
    return df
def fpbase_spectra_to_df():
    data = get_fpbase_spectra()
    df = pd.json_normalize(data,'spectra',['name'])
    df['data'] = df['data'].apply(lambda x: list(zip(*x))) # 
    df['role'] = np.where(df['state'].map(lambda x: 'ex' in x), 'Acceptor', 'Donor')
    df = df.set_index(['name','state'])
    df = df.reindex(columns=['role', 'qy', 'ec', 'max', 'data'])
    df['max'] = df['max'].fillna(df['data'].map(lambda x: x[0][np.argmax(x[1])]))
    df['data'] = df['data'].map(lambda x: interpolate(x))
    df['total'] = df['data'].map(lambda x: np.trapz(list(x), list(x.index))) # turns out si.trapz is copied directly from np.trapz!
    return df
df = fpbase_spectra_to_df()

def verify_separation(): # Verifies that each "state" value in the FPBase spectra contains either 'ex' or 'em'
    data = get_fpbase_spectra()
    for fluor in data:
        for spectrum in fluor["spectra"]:
            state = spectrum["state"]
            if (state.find("em") >= 0) is (state.find("ex") >= 0): # if both or neither of 'ex', 'em' are found
                print (fluor["name"], state)
                return False
    return True
verify_separation()

PI = np.pi
AVOGADRO = 6.022e23 # molecules / mole

# donor_em and acceptor_ex should be records from fpbase_spectra_to_df()
def calc_forster_radius(donor, acceptor, kappa2 = 2./3, refractive_index = 1.4): 
    donor_em = donor.data.divide(donor.total) # nm^-1
    acceptor_ex = acceptor.data.multiply(acceptor.ec / max(acceptor.data)) # M^-1 cm^-1 = mol^-1 dm^3 cm^-1
    quant_yield = donor.qy # unitless

    coeff = 9. * np.log(10) * kappa2 * float(quant_yield) / (128. * PI**5 * AVOGADRO * refractive_index**4) # mol

    integrand = [donor_em.loc[lam] * acceptor_ex.loc[lam] * lam ** 4 for lam in WAVELENGTHS_USED]
    integral = si.simps(integrand, WAVELENGTHS_USED) * STEP_SIZE # nm^-1 mol^-1 dm^3 cm^-1 nm^4 nm = mol^-1 nm^6 10^17

    return (1e17 * coeff * integral)**(1./6) # nm
donor = df.loc['mTFP1','default_em']
acceptor = df.loc['mNeonGreen','default_ex']
plt.plot(list(donor.data.index),list(donor.data))
plt.plot(list(acceptor.data.index),list(acceptor.data))
plt.xlim(400, 700)
plt.ion()
calc_forster_radius(donor,acceptor)

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

def get_ut_spectra():
    with open('Utzinger Spectra.sql', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    return lines

def ut_spectra_to_pd():
    lines = get_ut_spectra()
    # I manually checked that any line that contains NULL is irrelevant (none of them contain data about dyes)
    attributes = [eval(x.strip()[:-1]) for x in lines[42:15748] if (x[0] == '(' and 'NULL' not in x)]
    items = [eval(x.strip()[:-1]) for x in lines[15769:17538] if (x[0] == '(' and 'NULL' not in x)]
    spectra = [eval(x.strip()[:-1]) for x in lines[17556:1104997] if (x[0] == '(' and 'NULL' not in x)]
    spectra.sort(key=lambda x: (x[0],x[1]))

    df_a = pd.DataFrame(attributes).rename(columns={0:'id',1:'attribute',2:'value'})
    df_a_pivoted = df_a.drop_duplicates(subset=['id','attribute']).pivot(index='id',columns='attribute',values='value')
    df_i = pd.DataFrame(items).set_index(0).drop(4, axis=1).rename(columns={1:'name',2:'type',3:'role'})
    df_all = pd.concat([df_i,df_a_pivoted],axis=1).applymap(lambda x: x.strip() if type(x) == str else x) # strip() all the string elements
    # df_s = pd.DataFrame(spectra).set_index(0).rename(columns={1:'wavelength',2:'value'})

    # DROP all the completely useless columns and rows
    df_all = df_all[(df_all.type == 'dye')] # we only want dyes. not filters, lights, mirrors, dyecats, filtersets, and ESPECIALLY not NaN's >:(
    df_all = df_all[df_all['role'].map(lambda x: x in 'EM,AB,EX' and x != '')]
    df_all['role'] = df_all['role'].map({'EM': 'Donor', 'AB': 'Acceptor', 'EX': 'Acceptor'})
    df_all = df_all.drop(columns=['type', '', '1P primary source', '2P', '2P brightness', '2P primary source', 'Additional Info', 'BANDCENTER BANDPASS', 'BANDPASS', 'Brightness', 'Brightness (% of DsRed)', 'CF', 
    'CF260', 'CF280', 'Category', 'Category of dye', 'Comment', 'Comment 2nd', 'Comments', 'Company', 'Contact', 'Data Source', 'Data entry', 'Data source', 'Field19', 'Fluorophore (typical)', 'Full Reference',
    'Have Item', 'Imaging quality', 'Intensity @ 200mm (mW/cm2', 'Location', 'MANUFACTURER', 'Manufacturer', 'Module Configurator', 'NOTES', 'NUMBER BANDS', 'Note', 'Notes', 'Number bands', 'PART NAME', 'PRIMARY SOURCE', 'Primary Source', 'Primary source',
    'R <0> Alexa 488', 'R<0> Val. with AF488 acceptor', 'R<0> Val. with AF546 acceptor', 'R<0> Val. with AF555 acceptor', 'R<0> Val. with AF568 acceptor', 'R<0> Val. with AF594 acceptor', 'R<0> Val. with AF647 acceptor', 'R<0> Value source', 'R<0> Values',
    'Reference', 'Rescale 2P by:', 'Rescaled 2P by:', 'SOURCE', 'Secondary Source', 'Secondary source', 'Source', 'Source for Extinction Coefficient', 'Source`', 'Total Power @ 200 mm (mW)', 'URL', 'URL Quantum Yield', 'URL have item', 'band edges', 
    'bandcenter/bandpass', 'block Ti:Sa', 'category', 'category of dye', 'component of set', 'dia (mm)', 'filename', 'filter lamda', 'filter pass', 'filter wavelength', 'fliter pass', 'full reference', 'have doc', 'have item', 'have spectr', 'have spectrum', 
    'how many', 'in database', 'in db', 'in system', 'include', 'lot #', 'manufacturer', 'max value 381 to 700\r\n', 'miscellaneous', 't0 5 bleach', 't0 5 bleach (min)', 't0 5 maturatin', 't0 5 maturation (hr)', 't0.5 bleach', 't0.5 bleach (min)', 't0.5 maturatin', 't0.5 maturation (hr)'])
    # Some R0 values are at http://www.invitrogen.com/site/us/en/home/References/Molecular-Probes-The-Handbook/tables/R0-values-for-some-Alexa-Fluor-dyes.html
    # AGGREGATE the columns that contain the same type of information but the labels are misspelled -_-
    ex_c = ['Ex  Coeff', 'Ex max', 'Ex. Coeff', 'Exctinction coefficient', 'Extinction', 'Extinction Coefficient', 'Extinction Coefficient:', 'Extinction Coeficient', 'extinction coefficient']
    df_all = merge_columns(df_all, ex_c, 'Extinction Coefficient')
    wav_at_ec = ['Molar Ex Coeff. Wavelength', 'Visible light Extinction coefficient wavelength', 'WAVELENGTH', 'Wavelength Used To Determine Ex  Coeff', 'Wavelength Used To Determine Ex. Coeff.', 'Wavelength at Ec', 
    'Wavelength at Ec (Table)', 'Wavelength at Ec (from source)', 'Wavelength at Ec (source)', 'Wavelength at Ec (table)', 'Wavelength at Max Absorption', 'Wavelength at Max Emission', 
    'Wavelength of Max Absorption', 'Wavelength of Max Excitation', 'Wavelength used to determine Ex. Coeff.', 'Wavelength used to measure Ext  Coefficient', 'Wavelength used to measure extinction', 
    'Wavlength used to measure extinction coefficient', 'visble light Ec wavelength', 'visible light EC wavelength', 'visible light Ec wavelength', 'visible light Ec wavlength', 'visible light at Ec wavelength', 'visible light at Ec wavelenth', 
    'visible light wavelength', 'visiblle light Ec wavelength', 'wavelength at EC (from table)', 'wavelength at Ec', 'wavelength at Ec (from source', 'wavelength at Ec (from source)', 
    'wavelength at Ec (from table)', 'wavelength at Ec or Em (table below)', 'wavlelngth at Ec (from source)', 'wavlength at Ec', 'wavlength at Ec (from source)']
    df_all = merge_columns(df_all, wav_at_ec, 'Wavelength at Extinction Coefficient')
    other_name = ['Common name', 'Acronym', 'Alternative', 'acronmym or common name', 'acronym or common name', 'acronyms and alternative names', 'acronyms or common names', 'alternate name', 'alternative name', 'common name']
    df_all = merge_columns(df_all, other_name, 'Acronym or Common Name')
    thicc = ['thick (mm)', 'thick(mm)', 'thickness']
    df_all = merge_columns(df_all, thicc, 'Thickness(mm)')

    arr = np.array(list(zip(*spectra)),dtype=float).T
    split_by_id = np.split(arr, np.where(np.diff(arr[:,0]))[0]+1) # split whenever the id changes; this is why we needed to sort spectra by id
    arr = [x.T for x in split_by_id]

    return df_all, arr

# df_all, arr = ut_spectra_to_pd()

def merge_columns(df, col_labels, new_col_name): # takes a DataFrame and returns a Series with the lowest value in each 
    merged_col = df[col_labels[0]]
    assert type(merged_col) == type(pd.Series(dtype=object))
    for col in col_labels[1:]:
        merged_col = merged_col.combine_first(df[col])
    df = df.drop(columns=col_labels)
    df[new_col_name] = merged_col
    return df
    

# TODO 
