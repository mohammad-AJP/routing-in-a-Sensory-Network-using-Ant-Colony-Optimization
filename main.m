clc;
clear all;
close all;
%% loading data %%
load('nodes.mat')
load('energy.mat')
old_energies = energies;
Communication_radius = 150;
distance = dist_calc(nods,Communication_radius);
pheromone = (distance<Inf) - (distance==0);
%% pheromone and eta initializing %%
for i=1:100
    pheromone(i,i)=0;
end
pheromone_old = pheromone;
eta = 1./distance;
for i=1:100
    eta(i,i)=0;
end
ant_num = 20;
Packet_Size = 4098; 
rho = 0.9;
E_elec = 50*10^(-9);
E_amp = 0.1*10^(-9);
E_recieve = Packet_Size*E_elec;
E_transfer = Packet_Size*E_elec + E_amp*Packet_Size*(distance.^2);
E_transfer = (Packet_Size*E_elec + E_amp*Packet_Size*(distance.^2)).*(E_transfer>E_recieve);
max_iteration = 1000;
max_round = 2;
Avg_Res = zeros(max_round,10);
Res_Stan_Dev = zeros(max_round,10);
Min_Res = zeros(max_round,10);
k = 0;
while k<max_round
    pheromone = pheromone_old;
    energies = old_energies;
    source_index = randi([1,100]);
    sink_index = randi([1,100]); 
    while source_index==sink_index
            sink_index = randi([1,100]);
    end
    source_nod = nods(source_index,:);
    sink_nod = nods(sink_index,:);
    neighbour_node = find_neighbour(sink_index,distance);
    for i=1:max_iteration
        alpha = i/max_iteration;
        dlta_pheromone = zeros(100,100);
        for j=1:ant_num
            delta = find_path(alpha,eta,pheromone,energies,E_recieve,E_transfer,1,source_index,sink_index,neighbour_node);
            dlta_pheromone = dlta_pheromone + delta;
        end
        pheromone = (1-rho)*pheromone + dlta_pheromone;
    end
    k = k+1;
    for i=3:12
        [~,energy_consumption] = energy_calc(alpha,eta,pheromone,E_recieve,E_transfer,i*100,source_index,sink_index,neighbour_node);
        Min_Res(k,i-2) = min(energies - energy_consumption);
        Avg_Res(k,i-2) = mean(energies - energy_consumption);
        Res_Stan_Dev(k,i-2) = std(energies - energy_consumption);
    end
end
Avg_Res = mean(Avg_Res);
Res_Stan_Dev = mean(Res_Stan_Dev);
Min_Res = mean(Min_Res);
[path,energy_consumption] = energy_calc(1,eta,pheromone,E_recieve,E_transfer,1200,source_index,sink_index,neighbour_node);
X=(nods(:,1))';
Y=(nods(:,2))';
X_length=length(X);
for r=1:X_length
    x=X(r);
    y=Y(r);
    hold on
    if r==source_index
        plot(x,y,'rs','LineWidth',5,'MarkerSize',8);
    elseif r==sink_index
        plot(x,y,'rp','LineWidth',5,'MarkerSize',8);
    else
        plot(x,y,'bo','LineWidth',5,'MarkerSize',2);
    end
end
axis([0 1100 0 1100])
ind = zeros(1,1201);
count = 1;
for i=1:length(path)
    if path(i)==sink_index
        count = count + 1;
        ind(count) = i;
    end
end
for i=1:1200
    p=[];
    q=[];
    pat = path(ind(i)+1:ind(i+1));
    for qq=1:length(pat)
        p(qq) = nods(pat(qq),1);
    end
    for qq=1:length(pat)
        q(qq) = nods(pat(qq),2);
    end
    plot(p,q,'r')
end
ite = 3:12;
ite = 100*ite;
figure(2)
plot(ite,Res_Stan_Dev,'--r','LineWidth',5)
xlabel('Packets')
ylabel('Residual Energy Standard Deviation')
figure(3)
plot(ite,Avg_Res,'--r','LineWidth',5)
xlabel('Packets')
ylabel('Average Residual Energy')
figure(4)
plot(ite,Min_Res,'--r','LineWidth',5)
xlabel('Packets')
ylabel('Minimum Residual Energy')
