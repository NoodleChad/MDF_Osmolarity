This data was used as an imput on https://escher.github.io/#/app?map=iJO1366.Central%20metabolism&tool=Builder&model=iJO1366
in the data load reaction data.
To make those, first the CSV was devided in files with different growth rates. These files are called z_and_flux
Then this csv was converted to JSON using convertCSVtoJSON.py (All codes are found in the Tools_Figure_3 folder).
Then this the FWD and REV data was paired together to find the direction of the fluxes using Pair_fwd.py.
Finally, the reactions from the imL1515 were converted to the iJO1366 model to prepare the Escher maps using 4map.py
