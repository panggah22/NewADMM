%% Implementation part-by-part

%% Constraint 1. Active Power Balance
% Searching for bus connection
buscon10 = cell(data.num_bus,6);
for ii = 1:data.num_bus
    buscon10{ii,1} = data.bus(ii,1);
    buscon10{ii,2} = data.branch(data.branch(:,3)==ii,1);
    buscon10{ii,3} = data.branch(data.branch(:,2)==ii,1);
    buscon10{ii,4} = data.ficgen(data.ficgen(:,2)==ii,1);
    buscon10{ii,5} = data.evcs(data.evcs(:,2)==ii,1);
    buscon10{ii,6} = data.cb(data.cb(:,2)==ii,1);
end

Pg10 = zeros(data.num_bus,data.num_ficgen);
Pline10 = zeros(data.num_bus,data.num_branch);
Pevcs10 = zeros(data.num_bus,data.num_evcs);

% Input variable constant (-1 and 1) according to the signs of constraint
for ii = 1:data.num_bus
    Pline10(buscon10{ii,1},buscon10{ii,2}) = 1; % P_hi,t
    Pline10(buscon10{ii,1},buscon10{ii,3}) = -1; % P_ij,t
    Pg10(buscon10{ii,1},buscon10{ii,4}) = 1; % P_g,t
    Pevcs10(buscon10{ii,1},buscon10{ii,5}) = 1; % P_evcs
end

equ(10).Aeq = zeros(data.num_bus*ts,lent.total);
equ(10).Aeq(:,inp.Pline) = diagrep(Pline10,ts);
equ(10).Aeq(:,inp.Pg) = diagrep(Pg10,ts);
equ(10).Aeq(:,inp.Pevcs) = diagrep(Pevcs10,ts);

equ(10).beq = zeros(data.num_bus*ts,1); % Load's active power consumption
crapp = data.bus(:,3);
for ii = 1:ts
    equ(10).beq((ii-1)*data.num_bus+(1:data.num_bus)) = crapp .* data.loadcoeff(ii);
end

%% Constraint 2. Reactive Power Balance
equ(11).Aeq = zeros(data.num_bus*ts,lent.total);
equ(11).Aeq(:,inp.Qline) = diagrep(Pline10,ts);
equ(11).Aeq(:,inp.Qg) = diagrep(Pg10,ts);
equ(11).Aeq(:,inp.Qevcs) = diagrep(Pevcs10,ts);
Qcb11 = zeros(data.num_bus,data.num_cb);
for ii = 1:data.num_bus
    Qcb11(buscon10{ii,1},buscon10{ii,6}) = 1;
end
equ(11).Aeq(:,inp.Qcb) = diagrep(Qcb11,ts);

equ(11).beq = zeros(data.num_bus*ts,1); % Load's reactive power consumption
craqq = data.bus(:,4);
for ii = 1:ts
    equ(11).beq((ii-1)*data.num_bus+(1:data.num_bus)) = craqq .* data.loadcoeff(ii);
end
