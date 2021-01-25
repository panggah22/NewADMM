%% Implementation part-by-part
clear;

% Parameters
ts = 24;
unc_mean = 0;
unc_stdv = 1;
unc_val = 0.99;
maxiter = 500;
rho = 10;
isl = 1:33;
fic = 1;

data = data33(ts,unc_val,isl,fic);
evr = length(data.k2);
[len, lent] = lengthvar(data,evr,ts);
inp = inputvar(lent);

%% Objective Function : r * (Pline^2 + Qline^2)
%% Execute MIQP
H1 = zeros(lent.total);
H1(inp.Pline,inp.Pline) = diagrep(eye(len.Pline).*data.branch(:,6)*2,ts);
H1(inp.Qline,inp.Qline) = diagrep(eye(len.Qline).*data.branch(:,6)*2,ts);

H2 = zeros(lent.total);

H = H1+(rho/2*H2);

f = zeros(lent.total,1);

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

%% Constraint 3. Voltage Drop
V12 = zeros(data.num_branch,data.num_bus);
for ii = 1:data.num_branch
    V12(ii,data.branch(ii,2)) = -1;
    V12(ii,data.branch(ii,3)) = 1;
end
Pline12 = eye(len.Pline) .* data.branch(:,6);
Qline12 = eye(len.Qline) .* data.branch(:,7);

equ(12).Aeq = zeros(lent.Pline,lent.total);
equ(12).Aeq(:,[inp.V, inp.Pline, inp.Qline]) = [diagrep(V12,ts), diagrep(Pline12,ts), diagrep(Qline12,ts)];
equ(12).beq = zeros(size(equ(12).Aeq,1),1);

%% Constraint 4. Voltage Limit
% A
V13a = -eye(len.V);
b13a = -ones(len.V,1) .* 0.9;

% B
V13b = eye(len.V);
b13b = ones(len.V,1) .* 1.1;

ineq(13).A = zeros(lent.V*2,lent.total);
ineq(13).A(:,inp.V) = [diagrep(V13a,ts); diagrep(V13b,ts)];
ineq(13).b = [repmat(b13a,ts,1); repmat(b13b,ts,1)];

%% Constraint 5. Tap Changer
V14 = zeros(1,len.V);
V14(1) = 1;
Tap14 = -eye(size(V14,1)) .* 0.00625;

equ(14).Aeq = zeros(ts,lent.total);
equ(14).Aeq(:,inp.V) = diagrep(V14,ts);
equ(14).Aeq(:,inp.Tap) = diagrep(Tap14,ts);
equ(14).beq = ones(size(equ(14).Aeq,1),1);

%% Concatenate constraints
Cequ = cell(size(equ,2),1);
Cbequ = Cequ;
for i = 1:size(equ,2)
    Cequ{i} = equ(i).Aeq;
    Cbequ{i} = equ(i).beq;
end
Aeq = cell2mat(Cequ);
beq = cell2mat(Cbequ);

Cineq = cell(size(ineq,2),1);
Cbineq = Cineq;
for i = 1:size(ineq,2)
    Cineq{i} = ineq(i).A;
    Cbineq{i} = ineq(i).b;
end
Aineq = cell2mat(Cineq);
bineq = cell2mat(Cbineq);

%% Bounds
lb = -inf(lent.total,1);
lb(inp.PEV) = 0;
lb(inp.bin) = 0;
lb(inp.Tap) = -16;

ub = inf(lent.total,1);
ub(inp.bin) = 1;
ub(inp.Tap) = 16;
ub([inp.cmin inp.cmax]) = 20;

%% Continuous and integer
ctypenum = 67*ones(1,lent.total);
ctypenum(inp.bin) = 73;

ctype = char(ctypenum);

options = cplexoptimset('cplex');
options.display = 'on';
% options.mip.tolerances.mipgap = 1/100;
options.timelimit = 600;
%% Execute MIQP
[xx, fval, exitflag, output] = cplexmiqp (H, f, Aineq, bineq, Aeq, beq,...
    [], [], [], lb, ub, ctype, [], options);