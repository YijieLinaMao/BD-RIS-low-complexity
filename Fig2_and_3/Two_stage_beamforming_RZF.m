function sum_rate = Two_stage_beamforming_RZF(G,H,E,SNR)
%TWO_STAGE_BEAMFORMING_ZF 此处显示有关此函数的摘要
%   此处显示详细说明

Theta=symuni(H*G'*E');
F=G+E'*Theta'*H;

sum_rate=classical_RZF(F,SNR);
end

