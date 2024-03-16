# A low-complexity Beamforming Design for Beyond-Diagonal RIS aided Multi-User Networks.

This is a code package related to the following paper:

T. Fang and Y. Mao, "A Low-Complexity Beamforming Design for Beyond-Diagonal RIS Aided Multi-User Networks," in IEEE Communications Letters, vol. 28, no. 1, pp. 203-207, Jan. 2024, doi: 10.1109/LCOMM.2023.3333411.

# Content of Code Package
Here is a detailed description of the package: 
- The code in all packages are implemented in Matlab environment, and prat of them is assisted by Manopt toolbox. 
- In the 'Figure 1' code package,  Fig. 1 of the above paper will be reproduced by running the Matlab script 'main.m'. By changing the variables 'N' (the number of BD-RIS elements), 'Type_BD_RIS'(the architecture of BD-RIS belong to), 'channel_number' (channel realizations), you can reproduce Fig. 1.
- In the 'Figure 2&3' code package,  Fig. 2-3 of the above paper will be reproduced by running the Matlab script 'main.m'. By changing the variable 'N','Type_BD_RIS' and the channel realizations, you can reproduce Fig. 2-3.


# Abstract of the Article
Beyond-diagonal reconfigurable intelligent surface (BD-RIS) has been proposed recently as a novel and generalized RIS architecture that offers enhanced wave manipulation flexibility and large coverage expansion. 

However, the beyond-diagonal mathematical model in BD-RIS inevitably introduces additional optimization challenges in beamforming design.

In this letter, we derive a closed-form solution for the BD-RIS passive beamforming matrix that maximizes the sum of the effective channel gains among users. 

We further propose a computationally efficient two-stage beamforming framework to jointly design the active beamforming at the base station and passive beamforming at the BD-RIS to enhance the sum-rate for a BD-RIS aided multi-user multi-antenna network.

Numerical results show that our proposed algorithm achieves a higher sum-rate while requiring less computation time compared to state-of-the-art algorithms. The proposed algorithm paves the way for practical beamforming design in BD-RIS aided wireless networks.


# License and Referencing
This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

# Acknowledgements
This work has been supported in part by the National Nature Science Foundation of China under Grant 62201347; and in part by Shanghai Sailing Program under Grant 22YF1428400.

