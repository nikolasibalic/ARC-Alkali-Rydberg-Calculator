import numpy as np

def si(datum,precision=4):
# format datum with SI abbreviation to specified precision (# digits)

    exponent = np.floor(np.log10(np.abs(datum)))
    expInt   = np.floor(exponent/3).astype('int')
    expRange = (expInt * 3).astype('double')

    digitsLeftOfDecimal  = exponent - expRange + 1
    digitsRightOfDecimal = np.max((precision - digitsLeftOfDecimal,0)) 
    
    newDatum = datum * 10**(-expRange);

    sisym = ('y','z','a','f','p','n','\mu','m','','k','M','G','T','P','E','Z','Y')
    if np.abs(expRange) <= 24:
        sym = " " + sisym[expInt + 8]
    else:
        sym = " x 10^{%d}"%expRange

    if digitsLeftOfDecimal == precision: # if the last significant figure is in the
                                         # ones place, add the decimal to indicate
                                         # it as such
        sym = "." + sym


    # Formally, if digitsLeftOfDecimal > precision, newDatum should be rounded off
    # to requested precision, but since we are showing no more than 3 digits left
    # of the decimal, it's probably better not to round off
        
    fmtString = "%%%d.%df%s"%(digitsLeftOfDecimal,digitsRightOfDecimal,sym);

    return fmtString%(newDatum)
