## Functions for the statistical analysis
# Slope calculation formula
def find_slope(x1,x2,interval):
    x = (x2-x1)/(interval)
    return x

# Max growth calculation
def max_growth(x1,x2,interval):
    x = math.log(x2/x1)/(interval)
    return x 

# Calculate the y-intercept
def find_intercept(a,x,y):
    b = y-a*x
    return b

# Calculate r2
def find_r2(y,y_pred):
    l = range(0,len(y))
    y_bar = sum(y)/len(y)
    SSres_temp= []
    SStot_temp = []
    for i in l:
        SSres_temp.append((y[i]-y_pred[i])**2)
        SStot_temp.append((y[i]-y_bar)**2)
    SSres = sum(SSres_temp)
    SStot = sum(SStot_temp)
    if SStot != 0:
        r2 = 1-(SSres/SStot)
    else:
        r2 = 0
    return r2