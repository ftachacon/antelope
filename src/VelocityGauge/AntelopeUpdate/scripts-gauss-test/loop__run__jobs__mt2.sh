#!/bin/bash

########################
#Local potential scan
mt2=(
0.0000
0.1765
0.3529
0.5294
0.7059
0.8824
1.0588
1.2353
1.4118
1.5882
1.7647
1.9412
2.1176
2.2941
2.4706
2.6471
2.8235
3.0000
3.1765
3.3529
3.5294
3.7059
3.8824
4.0588
4.2353
4.4118
4.5882
4.7647
4.9412
5.1176
5.2941
5.4706
5.6471
5.8235
6.0000
6.1765
6.3529
6.5294
6.7059
6.8824
7.0588
7.2353
7.4118
7.5882
7.7647
7.9412
8.1176
8.2941
8.4706
8.6471
8.8235
9.0000
)

down=0;         # Initial-phase
upper=52;       # Final-phase

dt=0.1;         # time-step
phi=1.571;      # magnetic flux phase
LCP=+01;        # ellipticity for LCP laser pulse
RCP=-01;        # ellipticity for RCP laser pulse

########################
#LOOP ON Local potential ratio M/t2
for((iparam=${down}; iparam<=${upper}; iparam++ ))
do
    echo "###################"
    echo "M/t2 = ${mt2[iparam]}"

    "./"launchingJobs.sh ${LCP} ${phi} ${mt2[iparam]} ${dt}
    "./"launchingJobs.sh ${RCP} ${phi} ${mt2[iparam]} ${dt}

done
echo "****************************  "
echo "End of launching jobs on m/t2 "
echo "----------------------------  "
