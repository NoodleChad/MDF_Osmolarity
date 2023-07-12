import pandas as pd
import math
import json
from functions import find_intercept, find_slope, find_r2

Structure = 1 # file type    
data = pd.read_csv('KCl.csv') # File name
interval = 0.25 # time interval between OD measurement
plate_row = ["A", "B", "C","D","E","F","G"] # non-empty row
col = range(2,12) # non-empty column
replicate = 6 # number of experimental replicates

## Blank information
blank = {} # initialize dictionary
# location of blank (A2) and treatments associated to it (B2,...,G7)
blank["A2"] = ["B2","C2","D2","E7","F7","G7"] 
blank["A3"] = ["B3","C3","D3","E8","F8","G8"]
blank["A4"] = ["B4","C4","D4","E9","F9","G9"]
blank["A5"] = ["B5","C5","D5","E10","F10","G10"]
blank["A6"] = ["B6","C6","D6","E11","F11","G11"]
blank["A7"] = ["B7","C7","D7","E2","F2","G2"]
blank["A8"] = ["B8","C8","D8","E3","F3","G3"]
blank["A9"] = ["B9","C9","D9","E4","F4","G4"]
blank["A10"] = ["B10","C10","D10","E5","F5","G5"]
blank["A11"] = ["B11","C11","D11","E6","F6","G6"]
## please enter blank information

# calculates size of plate
if Structure == 1:
    m = range(0,int((len(data['2'])-8)/9+1))
elif Structure == 2:
    m = range(0,int(len(data[list(data)[1]])))
n = range(0,len(plate_row))

# Will create a dictionary (raw) using all the data, will contain all raw reads
# over the entire experiment
raw = {}
if Structure == 1:
    for i in col:
        for x in plate_row:
            globals()["temp"+x] = []
        temp = []
        temp.append(data[str(i)])
        for j in m:
            for k in n:
                    globals()["temp"+plate_row[k]].append(temp[0][j*9+k])
        for k in n:
            pos = plate_row[k]+str(i)
            raw[pos]=[]    
            raw[list(raw)[(i-2)*(len(plate_row))+k]] = globals()["temp"+plate_row[k]]
elif Structure == 2:
    raw = data


# from this the we will create a dictionary with the actual OD
OD = {}
l = range(0,len(list(blank)))
o = range(0,replicate)
for i in l:
    for j in o:
        temp = []
        for k in m:
            temp.append(raw[blank[list(blank)[i]][j]][k] - raw[list(blank)[i]][j])
        OD[blank[list(blank)[i]][j]] = temp

# Now we can start calculating the slope between each points
slope = {}
u_max = {}
max_itt = 12
q = range(0,len(list(OD)))
for i in q:
    umax = []
    itt = 0
    while not umax and itt <= max_itt:
        slope_range = 20
        slope_range = slope_range - itt
        if Structure == 1:
            p = range(0,int((len(data[list(data)[0]])-8)/9+1)-slope_range)
        elif Structure == 2:
            p = range(0,int(len(data[list(data)[1]]))-slope_range)
        r = range(0,slope_range)
        interval = (0.25*slope_range)
        slope_val = []
        r2 = []
        for j in p:
            if OD[list(OD)[i]][j] > 0 and OD[list(OD)[i]][j+slope_range] > 0:
                slope_val.append(find_slope(math.log(OD[list(OD)[i]][j]),math.log(OD[list(OD)[i]][j+slope_range]),interval))
                a = slope_val[-1]
                b = find_intercept(a,(j+slope_range)*0.25,math.log(OD[list(OD)[i]][j+slope_range]))
                pred = []
                real = []
                for k in r:
                    pred.append(a*(0.25*(j+k))+b)
                    if OD[list(OD)[i]][j+k] > 0:
                        real.append(math.log(OD[list(OD)[i]][j+k]))
                    else:
                        real.append(0.001)
                r2.append(find_r2(real,pred))
            else:
                slope_val.append(0.01)
                r2.append(0)
        temp = slope_val.copy()
        temp.sort()
        itt_1 = 1
        while itt_1 <=3:
            if r2[slope_val.index(temp[-itt_1])] >= 0.99:
                umax = slope_val[slope_val.index(temp[-itt_1])]
                slope[list(OD)[i]] = slope_val
                slope[list(OD)[i]+"_r2"] = r2
                u_max[list(OD)[i]] = umax
                break
            else:
                itt_1 = itt_1+1
        else:
            itt = itt + 1
    if itt == max_itt+1:
        slope[list(OD)[i]] = []
        slope[list(OD)[i]+"_r2"] = []
        u_max[list(OD)[i]] = []
treat = {}
for i in l:
    temp = []
    for j in o:
        if not u_max[list(u_max)[i*replicate+j]]:
            temp.append(u_max[list(u_max)[i*replicate+j]])
        else:
            temp.append(round(u_max[list(u_max)[i*replicate+j]],4))
    treat["t"+str(i)] = temp
    
print("Results")
for i in l:
    print(f"Treatment {i}: { treat[list(treat)[i]]}")

#with open("./NaCl_low.json", "w") as f:
 #   json.dump(treat,f)
