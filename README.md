# Robust control with an anti-windup technique based in relaxed LMI conditions for LTV system

![Figure](https://github.com/roscibely/robust-predictive-control-with-AW/blob/main/abstractpng.png)

Abstract: This paper proposes a new technique to address the anti-windup (AW) with model predictive control (MPC) scheme for linear time-varying (LTV) systems. The main advantage of this new approach is the reduced conservativeness compared with other well-known anti-windup techniques and to prevent integration windup in MPC controllers when the actuators are saturated. The control with AW is applied in a three-state switching cell (3SSC) DC-DC converter operating under saturation conditions. The MPC with proposed anti-windup is compared with the MPC technique and with MPC-AW without relaxation. The MPC-AW with relaxation improves the performance when the converter operated in the saturated mode and allows the rational use of the converter. The simulation results validated the efficiency of the proposed approach and showed that the proposed approach not only allows working better with the polytope modeling but also improves the response under LTV disturbance.

## Package required:
   
   📍[YALMIP](https://yalmip.github.io/)
   📍[SEDUMI](https://yalmip.github.io/solver/sedumi/) 
   
## How to run: 
   run the [main file](https://github.com/roscibely/robust-predictive-control-with-AW/blob/main/main_conversor.m) 

## How to cite:
    @article{rego2020robust,
     title={Robust control with an anti-windup technique based in relaxed LMI conditions for LTV system},
     author={Rego, Rosana and Costa, Marcus},
     journal={International Journal of Modelling, Identification and Control},
     volume={35},
     number={4},
     pages={298--304},
     year={2020},
     publisher={Inderscience Publishers (IEL)}
    }

Rego, Rosana, and Marcus Costa. "[Robust control with an anti-windup technique based in relaxed LMI conditions for LTV system.](https://www.inderscienceonline.com/doi/abs/10.1504/IJMIC.2020.114785)" International Journal of Modelling, Identification and Control 35.4 (2020): 298-304.
