function sum_rate = Two_stage_beamforming_FP(G,H,E,tolerance,SNR)
%TWO_STAGE_BEAMFORMING_FP 


Theta=symuni(H*G'*E');
F=G+E'*Theta'*H;

[~,sum_rate]=FP_algorithm(F,SNR,tolerance);

end

