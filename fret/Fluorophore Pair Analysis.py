import requests
r = requests.get('https://www.fpbase.org/api/proteins/spectra/?format=json')

import json
# with open('spectra.json', 'w', encoding='utf-8') as f:
#     json.dump(r.json(), f, ensure_ascii=False, indent=4)

# TODO ask API authors for their data on non-protein dyes that are not in API

# Take into account:
# - spectral overlap
# - distance between emission spectra
# - quantum yield?
# - extinction coeff?

# Here I'm only concerned with Forster radius. Ideally we also minimize the overlap
# between donor emission and acceptor emission, but that's not my goal with this code.

# Naive fret function, numerically integrates spectra then multiplies by qy and ec
def fret (em, ex):
    emArray = em["data"] # array of [wavelength, intensity] arrays
    exArray = ex["data"]

    minLam = int(max(emArray[0][0], exArray[0][0]))
    maxLam = int(min(emArray[-1][0], exArray[-1][0]))
    emIdx = 0
    while emArray[emIdx][0] < minLam:
        emIdx += 1
    exIdx = 0
    while exArray[exIdx][0] < minLam:
        exIdx += 1

    J = 0    
    for lam in range(minLam, maxLam+1):
        J += emArray[emIdx][1] * exArray[exIdx][1]
        emIdx += 1
        exIdx += 1
    J *= em["qy"] * ex["ec"]
    return J

with open('FPBase Data.json', 'r', encoding='utf-8') as f:
    fpbase_data = json.load(f)

def test1(): 
    fluor1 = data[0]["spectra"]
    fluor2 = data[1]["spectra"]
    # ex is 0, em is 1
    print ("Ex:" + data[0]["name"])
    print ("Em:" + data[1]["name"])
    print ("Relative efficiency: " + str(fret(fluor2[1], fluor1[0])))
    print ("Ex:" + data[1]["name"])
    print ("Em:" + data[0]["name"])
    print ("Relative efficiency: " + str(fret(fluor1[1], fluor2[0])))
    # expect one of these to be zero
# Realized that some fluors in data have more than two spectra,
# accounting for multiple states

# Checked that I can separate the spectra into excitation and emission 
# since no spectra states have both or neither of "ex", "em"
def verify_separation(data):
    for fluor in data:
        for spectrum in fluor["spectra"]:
            state = spectrum["state"]
            if state.find("em") >= 0 and state.find("ex") >= 0:
                print (fluor["name"])
                print (state)
            if state.find("em") < 0 and state.find("ex") < 0:
                print (fluor["name"])
                print (state)

def clean(spectrum):
    arr = [0] * 1000
    for item in spectrum:
        arr[int(item[0])] = item[1]
    return arr

# generate separate JSON lists for emission and excitation data from FPBase
# Mostly (if not all) fluorescent proteins, few to no organic dyes etc.
def separate(data):
    em = []
    ex = []
    for fluor in data:
        for spectrum in fluor["spectra"]:
            state = spectrum["state"]
            newItem = {}
            newItem["name"] = fluor["name"]
            newItem["state"] = state
            newItem["max"] = spectrum["max"]

            if state.find("em") >= 0:
                newItem["qy"] = spectrum["qy"]
                newItem["data"] = clean(spectrum["data"])
                em.append(newItem)

            else: # if not emission, must be excitation
                newItem["ec"] = spectrum["ec"]
                newItem["data"] = clean(spectrum["data"])
                ex.append(newItem)

    def getMax(spectrum):
        if spectrum["max"] is None:
            return 12345
        return spectrum["max"]

    em.sort(key=getMax)
    with open('FPBase Emission.json', 'w', encoding='utf-8') as f:
        json.dump(em, f, ensure_ascii=False, indent=4)
    print("FPBase Emission done")
    ex.sort(key=getMax)
    with open('FPBase Excitation.json', 'w', encoding='utf-8') as f:
        json.dump(ex, f, ensure_ascii=False, indent=4)
    print ("FPBase Excitation done") # TODO are em and ex flipped? 

# separate(fpbase_data)
del(fpbase_data)

# From Arizona 
# Cleaned files: turned into JSON format, removed some NULL entries

def arizona_compile_and_separate():
    arizona_data = [{} for i in range(2000)]

    with open('Arizona Items.json', 'r', encoding='utf-8') as f:
        arizona_items = json.load(f)
    for item in arizona_items:
        arizona_data[int(item[0])]["name"] = item[1]
        arizona_data[int(item[0])]["state"] = item[3]
    del(arizona_items)

    with open('Arizona Info.json', 'r', encoding='utf-8') as f:
        arizona_info = json.load(f)
    for item in arizona_info:
        arizona_data[int(item[0])][item[1]] = item[2]
    del(arizona_info)

    with open('Arizona Spectra.json', 'r', encoding='utf-8') as f:
        arizona_spectra = json.load(f)
    for fluor in arizona_data:
        fluor["data"] = [0]*1000
    for item in arizona_spectra:
        if "state" in arizona_data[int(item[0])] and any(arizona_data[int(item[0])]["state"] == i for i in ["EX", "EM", "AB"]):
            if item[1] >= 1000 or item[2] > 1: 
                # TODO fix this... do i really need to worry about any wavelengths over 1000 nm?
                # and I don't know what unit is being used when item[2] goes over 1.
                # It's definitely not relative intensity anymore.
                continue
            arizona_data[int(item[0])]["data"][int(item[1])] = item[2]
        # arizona_spectra contains non-integer wavelengths; I ignore and round.
        # jk this might be a bad idea, for the spectra that also skip wavelengths. 
        # I think I have to just make a smarter Forster radius function
        # which accounts for this with thicker/thinner Riemann sum.
    del(arizona_spectra)

    with open('Arizona Data.json', 'w', encoding='utf-8') as f:
        json.dump(arizona_data, f, ensure_ascii=False, indent=4)
    
    em = []
    ex = []
    for fluor in arizona_data:
        state = fluor["state"]
        newItem = {}
        newItem["name"] = fluor["name"]
        newItem["state"] = state
        newItem["max"] = fluor["max"] #TODO This goes by different names in Arizona

        if state == "EM":
            newItem["qy"] = fluor["qy"] #TODO this too
            newItem["data"] = fluor["data"]
            em.append(newItem)

        else: # if not emission, must be excitation
            newItem["ec"] = fluor["ec"] #TODO this too
            newItem["data"] = fluor["data"]
            ex.append(newItem)

    def getMax(fluor):
        if fluor["max"] is None:
            return 12345
        return fluor["max"]

    em.sort(key=getMax)
    with open('Arizona Emission.json', 'w', encoding='utf-8') as f:
        json.dump(em, f, ensure_ascii=False, indent=4)
    print("FPBase Emission done")
    ex.sort(key=getMax)
    with open('Arizona Excitation.json', 'w', encoding='utf-8') as f:
        json.dump(ex, f, ensure_ascii=False, indent=4)

    print ("Arizona done")
# arizona_compile_and_separate()

import scipy.integrate
import numpy as np
PI = np.pi
AVOGADRO = 6.022e23 # molecules / mole
LIGHT_SPEED = 3e8 # m/s
# Forster radius
#TODO brUH WHAT UNITS DOES THIS OUTPUT
def calc_forster_radius(donor_emission = None, acceptor_absorbance = None, wavelength_low = None, wavelength_high = None, kappa2 = 2./3, quant_yield = 1.0, refractive_index = 1.4):
    # donor_emission: 1D numpy array, units of inverse distance
    # acceptor_absorbance: 1D numpy array, in terms of wavelength, units of inverse concentration * inverse distance
    # wavelength_low, wavelength_high: float, units of distance
    # kappa2 (default 2/3): float, unitless
    # quant_yield (default 1.0): float, unitless
    # refractive_index: float, unitless

    assert(donor_emission is not None)
    assert(acceptor_absorbance is not None)
    assert(len(donor_emission) == len(acceptor_absorbance))

    if quant_yield is None: #TODO fix this abomination? why are some of the quantum yields 0
        quant_yield = 1.0
    coeff = 9. * np.log(10) * kappa2 * quant_yield / (128. * PI**5 * AVOGADRO * refractive_index**4)
    # step_size = (wavelength_high - wavelength_low) / (len(donor_emission)-1)
    # step_size = 1
    em, ab = interpolate(donor_emission, acceptor_absorbance)
    
    # integrate
    integral = 0
    for i in range(len(em) - 1):
        step_size = em[i+1][0] - em[i][0]
        lam = em[i][0]
        integral += step_size * lam**4 * em[i][1] * ab[i][1]
    # integral = scipy.integrate.simps(donor_emission * acceptor_absorbance * np.arange(wavelength_low, wavelength_high+step_size/2, step_size)**4) * step_size
    # print (coeff, integral)
    return (coeff * integral)**(1./6)

def interpolate(a, b):
    a.sort(key=lambda x: x[0])
    b.sort(key=lambda x: x[0])
    floor = max(a[0][0], b[0][0])
    ceil = min(a[-1][0], b[-1][0])
    ai = 0
    bi = 0
    while a[ai][0] < floor:
        ai += 1
    while b[bi][0] < floor:
        bi += 1
    
    while a[ai][0] < ceil or b[bi][0] < ceil:
        if a[ai][0] > b[bi][0]:
            a.append([b[bi][0], a[ai-1][1] + (a[ai][1] - a[ai-1][1])
            *(b[bi][0]-a[ai-1][0])
            /(a[ai][0]-a[ai-1][0])])
            bi += 1
        elif b[bi][0] > a[ai][0]:
            b.append([a[ai][0], b[bi-1][1] + (b[bi][1] - b[bi-1][1])
            *(a[ai][0]-b[bi-1][0])
            /(b[bi][0]-b[bi-1][0])])
            ai += 1
        else: 
            ai += 1
            bi += 1
            
    
    if a[ai][0] > b[bi][0]:
        a.append([b[bi][0], a[ai-1][1] + (a[ai][1] - a[ai-1][1])
        *(b[bi][0]-a[ai-1][0])
        /(a[ai][0]-a[ai-1][0])])
    elif b[bi][0] > a[ai][0]:
        b.append([a[ai][0], b[bi-1][1] + (b[bi][1] - b[bi-1][1])
        *(a[ai][0]-b[bi-1][0])
        /(b[bi][0]-b[bi-1][0])])

    a.sort(key=lambda x: x[0])
    b.sort(key=lambda x: x[0])
    while a[0][0] < floor:
        del(a[0])
    while a[-1][0] > ceil:
        del(a[-1])
    while b[0][0] < floor:
        del(b[0])
    while b[-1][0] > ceil:
        del(b[-1])
    return a,b

a = [[1, 1],[2,2],[3,3],[5,5],[6,6]]
b = [[0,2],[0.5,2],[1.25,5],[1.75,6],[2,5],[2.25,4],[2.5,2]]
print (interpolate(a,b)[0])
print (interpolate(a,b)[1])

def main():
    with open('FPBase Emission.json', 'r', encoding='utf-8') as f:
        fpbase_emission = json.load(f)
    with open('FPBase Excitation.json', 'r', encoding='utf-8') as f:
        fpbase_excitation = json.load(f)
    results = [0] * 75000
    i=0
    for donor in fpbase_emission:
        for acceptor in fpbase_excitation:
            result = ("donor: " + donor["name"],"acceptor: " + acceptor["name"], 
            # "Avg wavelength: " + str(donor["max"]/2 + acceptor["max"]/2), 
            "Forster radius: " + str(calc_forster_radius(np.array(donor["data"]), np.array(acceptor["data"]), wavelength_low=0, wavelength_high=999, quant_yield=donor["qy"])))
            results[i] = result
            i += 1
    def temp(x): #TODO deal with this...
        if isinstance(x, int):
            return -1000000
        return -x[2]
    results.sort(key=temp)
    results = results[results.count(0):]
    with open('Results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, ensure_ascii=False, indent=4)
    print("done ayyyyy")