#!/bin/bash

###########################
#Magnetic flux phase scan
phi=(
0.0000
0.0616
0.1232
0.1848
0.2464
0.3080
0.3696
0.4312
0.4928
0.5544
0.6160
0.6776
0.7392
0.8008
0.8624
0.9240
0.9856
1.0472
1.1088
1.1704
1.2320
1.2936
1.3552
1.4168
1.4784
1.5400
1.6016
1.6632
1.7248
1.7864
1.8480
1.9096
1.9712
2.0328
2.0944
2.1560
2.2176
2.2792
2.3408
2.4024
2.4640
2.5256
2.5872
2.6488
2.7104
2.7720
2.8336
2.8952
2.9568
3.0184
3.0800
3.1416
)

down=0;         # Initial-phase
upper=52;       # final-phase

dt=0.1;         # time-step
mt2=2.54;       # local potential
LCP=+01;        # ellipticity for LCP laser pulse
RCP=-01;        # ellipticity for RCP laser pulse

#################################
## LOOP ON magnetic flux phase
for((iparam=${down}; iparam<=${upper}; iparam++ ))
do
    echo "###################"
    echo "phi = ${phi[iparam]}"

    "./"launchingJobs.sh ${LCP} ${phi[iparam]} ${mt2} ${dt}
    "./"launchingJobs.sh ${RCP} ${phi[iparam]} ${mt2} ${dt}

done
echo "****************************  "
echo "End of launching jobs on phi  "
echo "----------------------------  "
