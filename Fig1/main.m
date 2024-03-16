% This is a code package related to the following journal paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T. Fang and Y. Mao, "A Low-Complexity Beamforming Design for Beyond-Diagonal" + ...
%" RIS Aided Multi-User Networks," in IEEE Communications Letters, vol. 28, no. 1,
% pp. 203-207, Jan. 2024, doi: 10.1109/LCOMM.2023.3333411.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The code is written by Tianyu Fang
%
% The code is implemented in Matlab environment
%
% Fig. 1 of the above paper will be reproduced by running the Matlab script 'main.m'.
% By changing the variables 'N' (the number of BD-RIS elements), 'channel_number' (channel realizations), 
% you can reproduce Fig. 1.
%
% Notice that the PoO method may take a lot of time.

clc
clear

tolerance=1e-5; %accuracy of convergence

L=4; %number of transmit antennas

K=4; %number of users

[data_PoP_F_norm,data_PoO_F_norm,data_DD_F_norm]=deal(zeros(100,5,3));  % preallocate the data. 
Without_RIS=zeros(100,1);
for channel_number=1:100 % iteration for 100 channel realizations

    for N_RIS=1:5  % iteration for element number from 2^2 to 2^6
  
        disp([channel_number N_RIS]) %output the progress bar
        N=2^(N_RIS+1);
        [G,H,E]=channel(channel_number,1,K,L,N);
        
        data_PoP_F_norm(channel_number,N_RIS,:)=Projection_of_Proposed(G,H,E);  %Using the Projection of the proposed solutions

        data_PoO_F_norm(channel_number,N_RIS,:)=Projection_of_Optimal(G,H,E,tolerance); %using the projection of the optimal solution of the relaxte problem

        data_DD_F_norm(channel_number,N_RIS,:)=Directly_Design(G,H,E,tolerance); %using the direct design method proposed in [5]

        Without_RIS(channel_number)=sqrt(trace(G*G'));
    end
end


figure

T0=Without_RIS;

T1=[T0,data_DD_F_norm(:,:,1)];
T2=[T0,data_PoP_F_norm(:,:,1)];
T3=[T0,data_DD_F_norm(:,:,2)];
T4=[T0,data_PoP_F_norm(:,:,2)];
T5=[T0,data_DD_F_norm(:,:,3)];
T6=[T0,data_PoP_F_norm(:,:,3)];
T7=[T0,data_PoO_F_norm(:,:,1)];
T8=[T0,data_PoO_F_norm(:,:,2)];
T9=[T0,data_PoO_F_norm(:,:,3)];

x=[2,4,8,16,32,64];
y1=mean(T1);
y2=mean(T2);
y3=mean(T3);
y4=mean(T4);
y5=mean(T5);
y6=mean(T6);
y7=mean(T7);
y8=mean(T8);
y9=mean(T9);

slg=semilogx(x,y2,'-o',x,y4,'-v',x,y6,'-d',x,y1,'--o',x,y3,'--v'...
    ,x,y5,'--d',x,y7,'-.o',x,y9,'-.v',x,y8,'-.d');
slg(1).LineWidth=1.5;
slg(2).LineWidth=1.5;
slg(3).LineWidth=1.5;
slg(4).LineWidth=1.5;
slg(5).LineWidth=1.5;
slg(6).LineWidth=1.5;
slg(7).LineWidth=1.5;
slg(8).LineWidth=1.5;
slg(9).LineWidth=1.5;
slg(1).Color=color(1);
slg(2).Color=color(1);
slg(3).Color=color(1);
slg(4).Color=color(3);
slg(5).Color=color(3);
slg(6).Color=color(3);
slg(7).Color=color(2);
slg(8).Color=color(2);
slg(9).Color=color(2);

set(gca,'XTick',x,'XTickLabel',[0,4,8,16,32,64])
xlabel('RIS Elements $N$','interpreter','latex')
ylabel('$f(\mathbf\Theta)$','interpreter','latex')
xlim([2 64])
grid on
legend({'FC-PoP','SC-PoP','GC-PoP','FC-DD','SC-DD','GC-DD'},...
    'Location','northwest')
axesNew=axes('position',get(gca,'position'),'visible','off');

legend(axesNew,[slg(7),slg(8),slg(9)],{'FC-PoO','SC-PoO','GC-PoO'},'Location','northeast')


