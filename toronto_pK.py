from equilibrator_api import ComponentContribution, Q_
cc = ComponentContribution()
import cobra
import pandas as pd

# Those are your parameters and model name
model = cobra.io.read_sbml_model("resources/iML1515.xml")
cc.p_h = Q_(7.5)
cc.p_mg = Q_(3.0)
cc.ionic_strength = Q_("0.25M")

# Find in the database, (change bigg metabolite names slighly to fit eQ)
## Seperate known from unknown for further steps down
pka = []
Name_known = []
Name_unknown = []
ID_unknown = []
ID_known = []
Formula_known = []
Formula_unknown = []
Charge_known = []
Charge_unknown = []
inchikey_known = []
for mets in model.metabolites:
    new_mets = mets.id.replace("_c", "")
    new_mets = new_mets.replace("_e", "")
    new_mets = new_mets.replace("_p", "")
    MetsSearch = "bigg.metabolite:"+ new_mets
    compound = cc.get_compound(MetsSearch)
    if hasattr(compound, 'dissociation_constants'):
        pka.append(compound.dissociation_constants)
        Name_known.append(mets.name)
        ID_known.append(new_mets)
        Charge_known.append(mets.charge)
        Formula_known.append(mets.formula)
        inchikey_known.append(compound.inchi_key)
    else:
        Name_unknown.append(mets.name)
        ID_unknown.append(new_mets)
        Formula_unknown.append(mets.formula)
        Charge_unknown.append(mets.charge)

# Use search function only for the unknown (can be quite long)
## To ensure that those are the right metabolites, compare chemical formula

Success = []
newpKa = []
Formula = []
inchikey_unknown = []
l = (list(range(0,len(Name_unknown))))
for i in l:
    Unknown_compound = cc.search_compound(Name_unknown[i])
    newpKa.append(Unknown_compound.dissociation_constants)
    Formula.append(Unknown_compound.formula)
    inchikey_unknown.append(Unknown_compound.inchi_key)
    if Formula_unknown[i] == Unknown_compound.formula:
        Success.append('Yes')
    else:
        Success.append('No')

# Now visually inspect the compounds where there is no match
## most of the time it is only a different charged species and can be used anyway

k = [i for i,x in enumerate(Success) if x=='No']
for i in k:
    print(Name_unknown[i])
    print(Formula_unknown[i])
    print(Formula[i])
    print(newpKa[i])

# Create a class with all the compounds now that we know them all
## make sure that the size fits the model you are using

l = (list(range(0,len(Name_unknown))))
for i in l:
    Name_known.append(Name_unknown[i])
    pka.append(newpKa[i])
    ID_known.append(ID_unknown[i])
    Formula_known.append(Formula[i])
    Charge_known.append(Charge_unknown[i])
    inchikey_known.append(inchikey_unknown[i])
print(len(Name_known))
print(len(pka))
print(len(ID_known))
print(len(Formula_known))
print(len(Charge_known))
print(len(inchikey_known))

# Find minimal and maximal charges
## BiGG uses pH 7.2, thus it will be our reference

min_charge = []
max_charge = []
l = list(range(0,len(Name_known)))
for k in l:
    C = Charge_known[k]
    if pka[k] is None:
        max_charge.append(C)
        min_charge.append(C)
    else:
        m = list(range(0,len(pka[k])))
        for n in m:
            if pka[k][n] > 7.2:
                C = C - 1
        min_charge.append(C)

# Find effective charges
## I: ionic strength
## T: temperature in Kelvin
## a: Deviation from standard temperature
## P: Effective polynomial
### Find the pKa at I = 0.1

pka_I = []
I = 0.25
T = 37 + 273.15
a = 1.10708 - 1.54508*10**(-3)*T + 5.95584*10**(-6)*T**(2)
B = 1.6
f = (a*I**(1/2))/(1+B*I**(1/2))
K_I = []
l = list(range(0,len(Name_known)))
for k in l:
    if pka[k] is None:
        pka_I.append(None)
        K_I.append(None)
    elif pka[k] == []:
        pka_I.append(None)
        K_I.append(None)
    else:
        pka_kn = []
        K_kn = []
        m = list(range(0,len(pka[k])))
        for n in m:
            Temp = pka[k][n]+(1-2*min_charge[k]-2*(n+1))*f
            pka_kn.append(Temp)
            K_kn.append(10**-Temp)
        pka_I.append(pka_kn)
        K_I.append(K_kn)

### Find the binding polynomial with those new pKa(I)

import numpy
pH = 7.5
P = []
for k in l:
    if pka_I[k] is None:
        P.append(None)
    else:
        Temp = []
        P_kn = []
        m = len(pka_I[k])
        while m > 0:
            Temp = ((10**-pH)**m)/numpy.prod(K_I[k][0:m])
            P_kn.append(Temp)
            m = m - 1
        P.append(1+sum(P_kn))

### Sum all the fractions to get the effective charge

Fractions = []
Eff_charge = []
m = []
l = list(range(0,len(Name_known)))
for k in l:
    if pka_I[k] is None:
        Fractions.append(1)
        Eff_charge.append(min_charge[k])
    else:
        Temp = []
        frac_kn = []
        charge_kn = []
        m = len(pka_I[k])
        n = 0
        while n < m:
            Temp = ((10**-pH)**n)/(numpy.prod(K_I[k][0:n]))*(1/P[k])
            frac_kn.append(Temp)
            charge_kn.append(Temp*(min_charge[k]+n))
            n = n + 1
        Fractions.append(frac_kn)
        Eff_charge.append(sum(charge_kn))

# Create a class with all the compounds now that we know them all
## make sure that the size fits the model you are using

class Model:
    def __init__(self, Name, Formula, ID, Charge, pka, K, P, Fractions, Inchikey):
        self.N = Name
        self.F = Formula
        self.ID = ID
        self.C = Charge
        self.pka = pka
        self.K = K
        self.P = P
        self.Frac = Fractions
        self.Ich = Inchikey

unique_strings = []
indices = {}

for i, string in enumerate(ID_known):
    if string not in unique_strings:
        unique_strings.append(string)
    else:
        if string in indices:
            indices[string].append(i)
        else:
            indices[string] = [unique_strings.index(string), i]

repeated_indices = [index for indices_list in indices.values() for index in indices_list[1:]]
repeated_indices.sort(reverse=True)

for i in repeated_indices:
    del ID_known[i]
    del Name_known[i]
    del Formula_known[i]
    del Eff_charge[i]
    del pka_I[i]
    del K_I[i]
    del P[i]
    del Fractions[i]
    del inchikey_known[i]

model = Model(Name_known, Formula_known, ID_known, Eff_charge, pka_I, K_I, P, Fractions, inchikey_known)

# Save
import pickle
with open(f'iML1515_metabolites.pickle', 'wb') as file:
    pickle.dump(model, file)
