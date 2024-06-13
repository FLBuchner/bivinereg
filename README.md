# Y-vine regression

An R package for Y-vine regression.

## Installation

    remotes::install_github("FLBuchner/bivinereg", quiet = TRUE)

## Example

    library(bivinereg)
    data(data)

    # fit model
    (fit <- bivinereg(cbind(U1, U4) ~ U2 + U3 + U5 + U6, family_Set = "parametric", data = data))
    # Y-vine regression model: U1, U4 | U3, U2, U5, U6 
    # nobs = 500, edf = 7, acll = 1230.62, acaic = -2447.23, acbic = -2417.73 

    summary(fit)
    #      var edf        cll        caic        cbic       p_value
    # 1 U1, U4  NA         NA          NA          NA            NA
    # 2     U3   2 667.382084 -1330.76417 -1322.33495 1.444254e-290
    # 3     U2   3 544.317155 -1082.63431 -1069.99049 1.063757e-235
    # 4     U5   1  11.013633   -20.02727   -15.81266  2.688046e-06
    # 5     U6   1   7.902634   -13.80527    -9.59066  7.020677e-05
