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
    with open('fpbase_emission.json', 'w', encoding='utf-8') as f:
        json.dump(em, f, ensure_ascii=False, indent=4)
    print("FPBase Emission done")
    ex.sort(key=getMax)
    with open('fpbase_excitation.json', 'w', encoding='utf-8') as f:
        json.dump(ex, f, ensure_ascii=False, indent=4)
    print ("FPBase Excitation done")

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
            if item[1] >= 1000: 
                # TODO fix this... do i really need to worry about any wavelengths over 1000 nm?
                continue
            arizona_data[int(item[0])]["data"][int(item[1])] = item[2]
        # arizona_spectra contains non-integer wavelengths; I ignore and round.
        # jk this might be a bad idea, for the spectra that also skip wavelengths. 
        # I think I have to just make a smarter Forster radius function
        # which accounts for this with thicker/thinner Riemann sum.
    del(arizona_spectra)

    with open('Arizona Data.json', 'w', encoding='utf-8') as f:
        json.dump(arizona_data, f, ensure_ascii=False, indent=4)
    
    
    print ("Arizona done")
arizona_compile_and_separate()