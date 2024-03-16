function [out] = color(index)
%COLOR 此处显示有关此函数的摘要
% 1-blue
% 2-orage
% 3-cyan
% 4-pink
% 5-green
% 6-purple
% 7-yellow
% 8-red
% 9-gray
% 10-black
mycolors=[    0.0353    0.5176    0.8902
    0.8824    0.4392    0.3333
         0    0.8078    0.7882
    0.9098    0.2627    0.5765
         0    0.7216    0.5804
    0.4235    0.3608    0.9059
    0.9922    0.7961    0.4314
    0.8392    0.1882    0.1922
    0.6980    0.7451    0.7647
    0.1765    0.2039    0.2118
    0.4549    0.7255    1.0000
    0.9804    0.6941    0.6275
    0.5059    0.9255    0.9255
    0.9922    0.4745    0.6588
    0.3333    0.9373    0.7686
    0.6353    0.6078    0.9961
    1.0000    0.9176    0.6549
    1.0000    0.4627    0.4588
    0.8745    0.9020    0.9137
    0.3882    0.4314    0.4471];
switch isnumeric(index)
    case 1
        out=mycolors(index,:);
    case 0
        x=linspace(0,4*pi);
        figure
        hold on
        for i=1:20
            y=sin(x-0.1*i);

            plot(x,y,'Color',mycolors(i,:),'LineWidth',1.5)
        end


end
end

