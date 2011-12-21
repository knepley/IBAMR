clear all;
clc;

load 'Drag'
DF = Drag;


Time = DF(:,1);
START = 1;
END   = size(DF,1);
Trep = 0.002*5;
NSteps = floor(Time(end)/Trep);

Dx   = -2*DF(START:END,2)  + 2*DF(START:END,5);
Dy   = DF(START:END,3)  - DF(START:END,6);
Dz   = DF(START:END,4)  - DF(START:END,7);

%Do a running avg;
k = 1;
for j = 1:NSteps
    target_time(j) = j*Trep;
    Favg_x = 0;
    Favg_y = 0;
    Favg_z = 0;
    counter = 0;
    
    for i = k:length(Time)
        if (Time(i) <= target_time(j))
            counter = counter +1;
            Favg_x    = Favg_x + Dx(i);
            Favg_y    = Favg_y + Dy(i);
            Favg_z    = Favg_z + Dz(i);
            
        else
            k = i
            break;
        end
    end
    
    Dx_target(j) = Favg_x/counter;
    Dy_target(j) = Favg_y/counter;
    Dz_target(j) = Favg_z/counter;
end

Cd = 2.09*ones(length(target_time),1);
plot(target_time,Dx_target, target_time, Cd , '--')
axis([0 17 0 5])
xlabel('t')
ylabel('C_{d}')