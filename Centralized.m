%% Implementation part-by-part
clear; tic

% Parameters
ts = 24;
unc_mean = 0;
unc_stdv = 0.7;
uncertainty = [0.90];
maxiter = 500;
rho = 10;
isl = 1:33;
fic = 1;
loadmultiplier = 1;

%% Options for enabling/disabling QV curve and uncertainties
%% A. enable QV curve
enQV = true;
% enQV = false;
%% B. enable PV uncertainty?
% enPV = true;
enPV = false;
%% C. enable EV uncertainty?
% enEV = true;
enEV = false;
if (enPV || enEV) == false
    uncertainty = 0.8;
end
    
for uu = 1:length(uncertainty)
unc_val = uncertainty(uu);
data = data33(ts,unc_val,isl,fic);
evr = length(data.k2);
[len, lent] = lengthvar(data,evr,ts);
inp = inputvar(lent);
data.bus(:,[3 4]) = data.bus(:,[3 4])*loadmultiplier;

%% Objective Function : r * (Pline^2 + Qline^2)
H1 = zeros(lent.total);
H1(inp.Pline,inp.Pline) = diagrep(eye(len.Pline).*data.branch(:,6)*2,ts);
H1(inp.Qline,inp.Qline) = diagrep(eye(len.Qline).*data.branch(:,6)*2,ts);

H2 = zeros(lent.total);

H = H1+(rho/2*H2);

f = zeros(lent.total,1);

%% Constraint 1. Active Power Balance
% Searching for bus connection
buscon1 = cell(data.num_bus,6);
for ii = 1:data.num_bus
    buscon1{ii,1} = data.bus(ii,1);
    buscon1{ii,2} = data.branch(data.branch(:,3)==ii,1);
    buscon1{ii,3} = data.branch(data.branch(:,2)==ii,1);
    buscon1{ii,4} = data.ficgen(data.ficgen(:,2)==ii,1);
    buscon1{ii,5} = data.evcs(data.evcs(:,2)==ii,1);
    buscon1{ii,6} = data.cb(data.cb(:,2)==ii,1);
end

Pg1 = zeros(data.num_bus,data.num_ficgen);
Pline1 = zeros(data.num_bus,data.num_branch);
Pevcs1 = zeros(data.num_bus,data.num_evcs);

% Input variable constant (-1 and 1) according to the signs of constraint
for ii = 1:data.num_bus
    Pline1(buscon1{ii,1},buscon1{ii,2}) = 1; % P_hi,t
    Pline1(buscon1{ii,1},buscon1{ii,3}) = -1; % P_ij,t
    Pg1(buscon1{ii,1},buscon1{ii,4}) = 1; % P_g,t
    Pevcs1(buscon1{ii,1},buscon1{ii,5}) = 1; % P_evcs
end

equ(1).Aeq = zeros(data.num_bus*ts,lent.total);
equ(1).Aeq(:,inp.Pline) = diagrep(Pline1,ts);
equ(1).Aeq(:,inp.Pg) = diagrep(Pg1,ts);
equ(1).Aeq(:,inp.Pevcs) = diagrep(Pevcs1,ts);

equ(1).beq = zeros(data.num_bus*ts,1); % Load's active power consumption
crapp = data.bus(:,3);
for ii = 1:ts
    equ(1).beq((ii-1)*data.num_bus+(1:data.num_bus)) = crapp .* data.loadcoeff(ii);
end

%% Constraint 2. Reactive Power Balance
equ(2).Aeq = zeros(data.num_bus*ts,lent.total);
equ(2).Aeq(:,inp.Qline) = diagrep(Pline1,ts);
equ(2).Aeq(:,inp.Qg) = diagrep(Pg1,ts);
equ(2).Aeq(:,inp.Qevcs) = diagrep(Pevcs1,ts);
Qcb2 = zeros(data.num_bus,data.num_cb);
for ii = 1:data.num_bus
    Qcb2(buscon1{ii,1},buscon1{ii,6}) = 1;
end
% equ(2).Aeq(:,inp.Qcb) = diagrep(Qcb2,ts);

equ(2).beq = zeros(data.num_bus*ts,1); % Load's reactive power consumption
craqq = data.bus(:,4);
for ii = 1:ts
    equ(2).beq((ii-1)*data.num_bus+(1:data.num_bus)) = craqq .* data.loadcoeff(ii);
end

%% Constraint 3. Voltage Drop
V3 = zeros(data.num_branch,data.num_bus);
for ii = 1:data.num_branch
    V3(ii,data.branch(ii,2)) = -1;
    V3(ii,data.branch(ii,3)) = 1;
end
Pline3 = eye(len.Pline) .* data.branch(:,6);
Qline3 = eye(len.Qline) .* data.branch(:,7);

equ(3).Aeq = zeros(lent.Pline,lent.total);
equ(3).Aeq(:,[inp.V, inp.Pline, inp.Qline]) = [diagrep(V3,ts), diagrep(Pline3,ts), diagrep(Qline3,ts)];
equ(3).beq = zeros(size(equ(3).Aeq,1),1);

%% Constraint 4. Voltage Limit
% A
V4a = -eye(len.V);
b4a = -ones(len.V,1) .* 0.9;

% B
V4b = eye(len.V);
b4b = ones(len.V,1) .* 1.1;

ineq(4).A = zeros(lent.V*2,lent.total);
ineq(4).A(:,inp.V) = [diagrep(V4a,ts); diagrep(V4b,ts)];
ineq(4).b = [repmat(b4a,ts,1); repmat(b4b,ts,1)];

%% Constraint 5. Tap Changer
V5 = zeros(1,len.V);
V5(1) = 1;
Tap5 = -eye(size(V5,1)) .* 0.00625;

equ(5).Aeq = zeros(ts,lent.total);
equ(5).Aeq(:,inp.V) = diagrep(V5,ts);
equ(5).Aeq(:,inp.Tap) = diagrep(Tap5,ts);
equ(5).beq = ones(size(equ(5).Aeq,1),1);

%% Constraint 6. CB Reactive Power
Qcb6b = eye(len.Qcb);
bcb6b = -eye(len.bcb) .* data.cb(:,4);
b6b = zeros(len.Qcb,1);

equ(6).Aeq = zeros(lent.Qcb,lent.total);
equ(6).Aeq(:,inp.Qcb) = diagrep(Qcb6b,ts);
equ(6).Aeq(:,inp.bcb) = diagrep(bcb6b,ts);
equ(6).beq = repmat(b6b,ts,1);

%% Constraint 7. SOC of EV (Sigma Form)
PEV7 = -eye(len.PEV) * 0.95 ./ data.evv(1);
bEV7 = repmat(cumsum(data.dr{3}),len.PEV,1);

equ(7).Aeq = zeros(lent.SocEV,lent.total);
equ(7).Aeq(:,inp.SocEV) = eye(lent.SocEV);
equ(7).Aeq(:,inp.PEV) = tril(repmat(PEV7,ts));

equ(7).beq = 0.2 - bEV7(:) ./ data.evv(1);

%% Constraint 8. Bounds of EV SoC
ineq(8).A = zeros(lent.SocEV*2,lent.total);
ineq(8).A(:,inp.SocEV) = [-eye(lent.SocEV); eye(lent.SocEV)];
ineq(8).b = [-ones(lent.SocEV,1) .* 0.2; ones(lent.SocEV,1)]; % 20 ~ 100 %

%% Constraint 9. Binary of EV
bEV9 = double(repmat(data.dr{3},data.num_ev,1)==0) .* data.evv(2);
ineq(9).A = zeros(lent.PEV*2,lent.total);
ineq(9).A(:,inp.PEV) = [-eye(lent.PEV); eye(lent.PEV)];
ineq(9).b = [zeros(lent.PEV,1); bEV9(:)];

%% Constraint 10-17. Reserved for the EV uncertainty
%% Constraint 10 >> Sigma form
ineq(10).A = repmat(zeros(lent.PEV,lent.total),length(data.k1),1);
al = cell(ts);
for ii = 1:ts
    for jj = 1:ts
        if jj <= ii
            al{ii,jj} = eye(len.PEV) .* 0.95 ./ repmat(data.evv(1),data.num_ev,1);
        else
            al{ii,jj} = zeros(len.PEV);
        end
    end
end
ineq(10).A(:,inp.PEV) = -repmat(cell2mat(al),length(data.k1),1);

celldr = cell(length(data.k1),1);
for kk = 1:length(data.k1)
    ilh = repmat(data.k1{kk},data.num_ev,1);
    ap = zeros(size(ilh));
    for ii = 1:size(ilh,1)
        ap(ii,:) = cumsum(ilh(ii,:)) ./ data.evv(1);
    end
    celldr{kk} = ap(:);
end
ineq(10).b = 0.2 - 0.2 - cell2mat(celldr);

%% Constraint 11 >> Sigma form
ineq(11).A = repmat(zeros(lent.PEV,lent.total),length(data.k1),1);
al = cell(ts);
for ii = 1:ts
    for jj = 1:ts
        if jj <= ii
            al{ii,jj} = eye(len.PEV) .* 0.95 ./ repmat(data.evv(1),data.num_ev,1);
        else
            al{ii,jj} = zeros(len.PEV);
        end
    end
end
ineq(11).A(:,inp.PEV) = repmat(cell2mat(al),length(data.k1),1);

celldr = cell(length(data.k1),1);
for kk = 1:length(data.k1)
    ilh = repmat(data.k1{kk},data.num_ev,1);
    ap = zeros(size(ilh));
    for ii = 1:size(ilh,1)
        ap(ii,:) = cumsum(ilh(ii,:)) ./ data.evv(1);
    end
    celldr{kk} = ap(:);
end
ineq(11).b = 1 - 0.2 + cell2mat(celldr);

%% Constraint 12
cell12 = cell(length(data.k1),1);
for ii = 1:length(data.k1)
    coy = double(repmat(data.k1{ii},data.num_ev,1)==0) .* data.evv(2);
    cell12{ii} = coy(:);
end
bEV12 = cell2mat(cell12);

ineq(12).A = repmat(zeros(lent.PEV,lent.total),length(data.k1),1);
ineq(12).A(:,inp.PEV) = repmat(eye(lent.PEV),length(data.k1),1);
ineq(12).b = bEV12;

%% Constraint 13 >> Sigma form
ineq(13).A = repmat(zeros(lent.PEV,lent.total),length(data.k2),1);
am = cell(ts);
for ii = 1:ts
    for jj = 1:ts
        if jj <= ii
            am{ii,jj} = eye(len.PEV) .* 0.95 ./ repmat(data.evv(1),data.num_ev,1);
        else
            am{ii,jj} = zeros(len.PEV);
        end
    end
end
ineq(13).A(:,inp.PEV) = -repmat(cell2mat(am),length(data.k2),1);

celldr2 = cell(length(data.k2),1);
for kk = 1:length(data.k2)
    ilh2 = repmat(data.k2{kk},data.num_ev,1);
    aq = zeros(size(ilh2));
    for ii = 1:size(ilh2,1)
        aq(ii,:) = cumsum(ilh2(ii,:)) ./ data.evv(1);
    end
    celldr2{kk} = aq(:);
end
ineq(13).A(:,inp.z) = diagrep(repmat(-eye(len.z),ts,1),length(data.k2)) .* cell2mat(celldr2);
ineq(13).b = 0.2 - 0.2 - cell2mat(celldr2);

%% Constraint 14 >> Sigma form
ineq(14).A = repmat(zeros(lent.PEV,lent.total),length(data.k2),1);
am = cell(ts);
for ii = 1:ts
    for jj = 1:ts
        if jj <= ii
            am{ii,jj} = eye(len.PEV) .* 0.95 ./ repmat(data.evv(1),data.num_ev,1);
        else
            am{ii,jj} = zeros(len.PEV);
        end
    end
end
ineq(14).A(:,inp.PEV) = repmat(cell2mat(am),length(data.k2),1);

celldr2 = cell(length(data.k2),1);
for kk = 1:length(data.k2)
    ilh2 = repmat(data.k2{kk},data.num_ev,1);
    aq = zeros(size(ilh2));
    for ii = 1:size(ilh2,1)
        aq(ii,:) = cumsum(ilh2(ii,:)) ./ data.evv(1);
    end
    celldr2{kk} = aq(:);
end
ineq(14).A(:,inp.z) = diagrep(repmat(-eye(len.z),ts,1),length(data.k2)) .* (data.evv(2)/data.evv(1)-cell2mat(celldr2));
ineq(14).b = 1 - 0.2 + cell2mat(celldr2);

%% Constraint 15
cell15 = cell(length(data.k2),1);
for ii = 1:length(data.k2)
    coy = double(repmat(data.k2{ii},data.num_ev,1)==0) .* data.evv(2);
    cell15{ii} = coy(:);
end
bEV15 = cell2mat(cell15);

ineq(15).A = repmat(zeros(lent.PEV,lent.total),length(data.k2),1);
ineq(15).A(:,inp.z) = -diagrep(repmat(eye(len.z),ts,1),length(data.k2)) .* data.evv(2);
ineq(15).A(:,inp.PEV) = repmat(eye(lent.PEV),length(data.k2),1);
ineq(15).b = bEV15;

%% Constraint 16
ineq(16).A = zeros(length(data.k2),lent.total);
ineq(16).A(:,inp.z) = diagrep(ones(1,len.z),length(data.k2)) .* data.rea(data.rea(:,5)<=(1-unc_val),5);
ineq(16).b = (1-unc_val) .* ones(length(data.k2),1);

% %% Constraint 16
% ineq(16).A = zeros(len.z,lent.total);
% cell16 = cell(1,length(data.k2));
% dat16 = data.rea(data.rea(:,5)<=(1-unc_val),5);
% for ii = 1:length(dat16)
%     cell16{ii} = dat16(ii) .* ones(1,len.z);6
% end
%     
% ineq(16).A(:,inp.z) = repmat(eye(len.z),1,length(data.k2)) .* cell2mat(cell16);
% ineq(16).b = (1-unc_val) .* ones(len.z,1);

%% Constraint 18. SOC of BESS
Pbessc18 = eye(len.Pbessc) .* 0.95 ./ data.evcs(:,4);
Pbessd18 = eye(len.Pbessd) ./ (0.95 * data.evcs(:,4));

equ(18).Aeq = zeros(lent.SocBess,lent.total);
equ(18).Aeq(:,inp.SocBess) = eye(lent.SocBess);
equ(18).Aeq(:,inp.Pbessc) = -tril(repmat(Pbessc18,ts));
equ(18).Aeq(:,inp.Pbessd) = tril(repmat(Pbessd18,ts));

equ(18).beq = 0.2 * ones(lent.SocBess,1);

%% Constraint 19. Bounds of BESS SoC
ineq(19).A = zeros(lent.SocBess*2,lent.total);
ineq(19).A(:,inp.SocBess) = [-eye(lent.SocBess); eye(lent.SocBess)];
ineq(19).b = [-ones(lent.SocBess,1) .* 0.2; ones(lent.SocBess,1)]; % Limit of Soc from 20% to 100%

%% Constraint 20. Binary of Pbessc
ineq(20).A = zeros(lent.Pbessc*2,lent.total);
ineq(20).A(:,inp.Pbessc) = [-eye(lent.Pbessc); eye(lent.Pbessc)];
ineq(20).A(:,inp.bBess) = [zeros(lent.bBess); -eye(lent.bBess) .* repmat(data.evcs(:,5),ts,1)];

ineq(20).b = zeros(lent.Pbessc*2,1);

%% Constraint 21. Binary of Pbessd
ineq(21).A = zeros(lent.Pbessd*2,lent.total);
ineq(21).A(:,inp.Pbessd) = [-eye(lent.Pbessd); eye(lent.Pbessd)];
ineq(21).A(:,inp.bBess) = [zeros(lent.bBess); eye(lent.bBess) .* repmat(data.evcs(:,5),ts,1)];

ineq(21).b = [zeros(lent.Pbessd,1); ones(lent.bBess,1) .* repmat(data.evcs(:,5),ts,1)];

%% Constraint 22. Power balance of EVCS
% Create sum of PEVs according to their respective EVCS position
am = cell(len.Pevcs);
for ii = 1:len.Pevcs % trace the column
    for jj = 1:len.Pevcs % trace the row
        if ii == jj
            am{jj,ii} = ones(1,data.evcs(ii,8));
        else
            am{jj,ii} = zeros(1,data.evcs(ii,8));
        end
    end
end

equ(22).Aeq = zeros(lent.Pevcs,lent.total);
equ(22).Aeq(:,inp.Pevcs) = -eye(lent.Pevcs);
equ(22).Aeq(:,inp.PEV) = -diagrep(cell2mat(am),ts);
equ(22).Aeq(:,inp.Pbessc) = -eye(lent.Pbessc);
equ(22).Aeq(:,inp.Pbessd) = eye(lent.Pbessd);

PVco = cell(ts,1);
for ii = 1:length(PVco)
    PVco{ii} = data.evcs(:,7) .* data.PVcoeff(ii);
end
equ(22).beq = -cell2mat(PVco);

%% Replacement for PV uncertainty using ESS SOC limit
% EV SOC calculation is located at constraint 18
% EV SOC limit is at constraint 19
%% Constraint 23
inv_cdf = icdf('Normal',unc_val,unc_mean,unc_stdv);

Pbessc23 = eye(len.Pbessc) .* 0.95 ./ data.evcs(:,4);
Pbessd23 = eye(len.Pbessd) ./ (0.95 * data.evcs(:,4));
PV23 = cumsum(data.evcs(:,7)*data.PVcoeff,2) ./ data.evcs(:,4);


ineq(23).A = zeros(lent.Pbessc,lent.total);
ineq(23).A(:,inp.Pbessc) =  tril(repmat(Pbessc23,ts));
ineq(23).A(:,inp.Pbessd) = -tril(repmat(Pbessd23,ts));

% - Sigma PV/Ecap * (std+mean) + Socmax - Socini
ineq(23).b = PV23(:) * (-unc_stdv*inv_cdf - unc_mean) + 1 - 0.2;


%% Constraint 24
ineq(24).A = zeros(lent.Pbessc,lent.total);
ineq(24).A(:,inp.Pbessc) = -tril(repmat(Pbessc23,ts));
ineq(24).A(:,inp.Pbessd) =  tril(repmat(Pbessd23,ts));

ineq(24).b = PV23(:) * (-unc_stdv*inv_cdf + unc_mean) - 0.2 + 0.2;

%% Constraint 25~27
% For constraint 29, polygon-based linearization is required. Go to
% Constraint 70 !!
Sij = data.evcs(:,6) .* sqrt((2*pi/6)/sin(2*pi/6));


%% Constraint 25
ineq(25).A = zeros(lent.Qevcs*2,lent.total);
ineq(25).A(:,inp.Pevcs) = [-sqrt(3) .* eye(lent.Pevcs); sqrt(3) .* eye(lent.Pevcs)];
ineq(25).A(:,inp.Qevcs) = [-eye(lent.Qevcs); eye(lent.Qevcs)];

ineq(25).b = ones(lent.Qevcs*2,1) .* (sqrt(3)*repmat(Sij,ts*2,1));

%% Constraint 26
ineq(26).A = zeros(lent.Qevcs*2,lent.total);
ineq(26).A(:,inp.Qevcs) = [-eye(lent.Qevcs); eye(lent.Qevcs)];

ineq(26).b = ones(lent.Qevcs*2,1) .* (sqrt(3)/2*repmat(Sij,ts*2,1));

%% Constraint 27
ineq(27).A = zeros(lent.Qevcs*2,lent.total);
ineq(27).A(:,inp.Qevcs) = [-eye(lent.Qevcs); eye(lent.Qevcs)];
ineq(27).A(:,inp.Pevcs) = [sqrt(3) .* eye(lent.Pevcs); -sqrt(3) .* eye(lent.Pevcs)];

ineq(27).b = ones(lent.Qevcs*2,1) .* (sqrt(3)*repmat(Sij,ts*2,1));

%% Constraint 30. Relation of Qevcs and Qstar
equ(30).Aeq = zeros(lent.Qstar,lent.total);
equ(30).Aeq(:,inp.Qevcs) = eye(lent.Qevcs);
equ(30).Aeq(:,inp.Qstar) = -eye(lent.Qstar) .* repmat(data.evcs(:,6),ts,1);

equ(30).beq = zeros(lent.Qstar,1);

%% QV curve constraints // the previous one does not seem reliable
%% Reserved for constraints 31-39

%% Constraint 31 // cmin<=cmax
ineq(31).A = zeros(lent.cmin,lent.total);
ineq(31).A(:,inp.cmin) = eye(lent.cmin);
ineq(31).A(:,inp.cmax) = -eye(lent.cmax);

ineq(31).b = zeros(lent.cmax,1);

%% Constraint 32 // cmin = sigma 2^m.lmin
lmin32 = ones(1,len.lmin/data.num_evcs) .* [2^0 2^1 2^2 2^3 2^4];

equ(32).Aeq = zeros(data.num_evcs,lent.total);
equ(32).Aeq(:,inp.cmin) = eye(lent.cmin);
equ(32).Aeq(:,inp.lmin) = -diagrep(lmin32,data.num_evcs);

equ(32).beq = zeros(lent.cmin,1);

%% Constraint 33 // cmax = sigma 2^m.lmax ... lmax constant = lmin constant
lmax33 = lmin32;

equ(33).Aeq = zeros(data.num_evcs,lent.total);
equ(33).Aeq(:,inp.cmax) = eye(lent.cmax);
equ(33).Aeq(:,inp.lmax) = -diagrep(lmax33,data.num_evcs);

equ(33).beq = zeros(lent.cmax,1);

%% Constraint 34 // BigM method
M = 100;
a34 = zeros(1,len.a/data.num_evcs); a34(3) = 1;
a34 = diagrep(repmat(a34,5,1),data.num_evcs);
lmin34 = diagrep(eye(len.lmin/data.num_evcs),data.num_evcs);

ineq(34).A = zeros(lent.wmin*2,lent.total);
ineq(34).A(:,inp.wmin) = [-eye(lent.wmin); eye(lent.wmin)];
ineq(34).A(:,inp.lmin) = [M*repmat(lmin34,ts,1); zeros(lent.wmin,lent.lmin)];
ineq(34).A(:,inp.a) = [diagrep(a34,ts); -diagrep(a34,ts)];

ineq(34).b = [M*ones(lent.wmin,1); zeros(lent.wmin,1)];

%% Constraint 35
ineq(35).A = zeros(lent.wmin*2,lent.total);
ineq(35).A(:,inp.wmin) = [-eye(lent.wmin); eye(lent.wmin)];
ineq(35).A(:,inp.lmin) = [zeros(lent.wmin,lent.lmin); -M*repmat(lmin34,ts,1)];

ineq(35).b = [zeros(lent.wmin,1); zeros(lent.wmin,1)];

%% Constraint 36
a36 = zeros(1,len.a/data.num_evcs); a36(4) = 1;
a36 = diagrep(repmat(a36,5,1),data.num_evcs);
lmax36 = diagrep(eye(len.lmax/data.num_evcs),data.num_evcs);

ineq(36).A = zeros(lent.wmax*2,lent.total);
ineq(36).A(:,inp.wmax) = [-eye(lent.wmax); eye(lent.wmax)];
ineq(36).A(:,inp.lmax) = [M*repmat(lmax36,ts,1); zeros(lent.wmax,lent.lmax)];
ineq(36).A(:,inp.a) = [diagrep(a36,ts); -diagrep(a36,ts)];

ineq(36).b = [M*ones(lent.wmax,1); zeros(lent.wmax,1)];

%% Constraint 37
ineq(37).A = zeros(lent.wmax*2,lent.total);
ineq(37).A(:,inp.wmax) = [-eye(lent.wmax); eye(lent.wmax)];
ineq(37).A(:,inp.lmax) = [zeros(lent.wmax,lent.lmax); -M*repmat(lmax36,ts,1)];

ineq(37).b = [zeros(lent.wmax,1); zeros(lent.wmax,1)];

%% Constraint 38
a38 = zeros(1,len.a/data.num_evcs); a38(3) = 1;
a38 = diagrep(a38,data.num_evcs);
wmin38 = diagrep(ones(1,5) .* [2^0 2^1 2^2 2^3 2^4],data.num_evcs);

equ(38).Aeq = zeros(lent.smin,lent.total);
equ(38).Aeq(:,inp.smin) = eye(lent.smin);
equ(38).Aeq(:,inp.a) = - 0.9 * diagrep(a38,ts);
equ(38).Aeq(:,inp.wmin) = - 0.01 * diagrep(wmin38,ts);

equ(38).beq = zeros(lent.smin,1);

%% Constraint 39
a39 = zeros(1,len.a/data.num_evcs); a39(4) = 1;
a39 = diagrep(a39,data.num_evcs);
wmax39 = diagrep(ones(1,5) .* [2^0 2^1 2^2 2^3 2^4],data.num_evcs);

equ(39).Aeq = zeros(lent.smax,lent.total);
equ(39).Aeq(:,inp.smax) = eye(lent.smax);
equ(39).Aeq(:,inp.a) = - 0.9 * diagrep(a39,ts);
equ(39).Aeq(:,inp.wmax) = - 0.01 * diagrep(wmax39,ts);

equ(39).beq = zeros(lent.smax,1);
%% Constraint 40
V40 = zeros(data.num_evcs,len.V);
for ii = 1:data.num_evcs
    V40(ii,data.bus(:,1)==data.evcs(ii,2)) = 1;
end
a40 = diagrep(ones(1,6) .* [0.8 0.9 0 0 1.1 1.2],data.num_evcs);

equ(40).Aeq = zeros(lent.smin,lent.total);
equ(40).Aeq(:,inp.V) = -diagrep(V40,ts);
equ(40).Aeq(:,inp.a) = diagrep(a40,ts);
equ(40).Aeq(:,inp.smin) = eye(lent.smin);
equ(40).Aeq(:,inp.smax) = eye(lent.smax);

equ(40).beq = zeros(lent.smin,1);

%% Constraint 41
a41 = ones(1,len.a/data.num_evcs) .* [1 1 0 0 -1 -1];
a41 = diagrep(a41,data.num_evcs);

equ(41).Aeq = zeros(lent.Qstar,lent.total);
equ(41).Aeq(:,inp.Qstar) = -eye(lent.Qstar);
equ(41).Aeq(:,inp.a) = diagrep(a41,ts);

equ(41).beq = zeros(lent.Qstar,1);

%% Constraint 42
a42a = zeros(1,len.a/data.num_evcs); a42a(1) = 1;
a42a = diagrep(a42a,data.num_evcs);
d42a = zeros(1,len.d/data.num_evcs); d42a(1) = 1;
d42a = diagrep(d42a,data.num_evcs);

a42b = zeros(1,len.a/data.num_evcs); a42b(6) = 1;
a42b = diagrep(a42b,data.num_evcs);
d42b = zeros(1,len.d/data.num_evcs); d42b(5) = 1;
d42b = diagrep(d42b,data.num_evcs);

ineq(42).A = zeros(data.num_evcs*ts*2,lent.total);
ineq(42).A(:,inp.a) = [diagrep(a42a,ts); diagrep(a42b,ts)];
ineq(42).A(:,inp.d) = [-diagrep(d42a,ts); -diagrep(d42b,ts)];

ineq(42).b = zeros(data.num_evcs*ts*2,1);

%% Constraint 43
a43 = eye(len.a/data.num_evcs); a43 = a43([2:5],:);
a43 = diagrep(a43,data.num_evcs);
d43 = time_relate(1,1,5);
d43 = diagrep(d43,data.num_evcs);

ineq(43).A = zeros(4*data.num_evcs*ts,lent.total);
ineq(43).A(:,inp.a) = diagrep(a43,ts);
ineq(43).A(:,inp.d) = -diagrep(d43,ts);

ineq(43).b = zeros(4*data.num_evcs*ts,1);

%% Constraint 44
ineq(44).A = zeros(lent.a,lent.total);
ineq(44).A(:,inp.a) = -eye(lent.a);

ineq(44).b = zeros(lent.a,1);

%% Constraint 45
a45 = diagrep(ones(1,len.a/data.num_evcs),data.num_evcs);

equ(45).Aeq = zeros(data.num_evcs*ts,lent.total);
equ(45).Aeq(:,inp.a) = diagrep(a45,ts);

equ(45).beq = ones(data.num_evcs*ts,1);

%% Constraint 46
d46 = diagrep(ones(1,len.d/data.num_evcs),data.num_evcs);

equ(46).Aeq = zeros(data.num_evcs*ts,lent.total);
equ(46).Aeq(:,inp.d) = diagrep(d46,ts);

equ(46).beq = ones(data.num_evcs*ts,1);

%% Constraint 47 boundary of Qevcs
% ineq(47).A = zeros(lent.Qevcs*2,lent.total);
% ineq(47).A(:,inp.Qevcs) = [-eye(lent.Qevcs); eye(lent.Qevcs)];
% 
% ineq(47).b = [repmat(data.evcs(:,9),ts,1); repmat(data.evcs(:,9),ts,1)];

%% Constraint 48 // power factor of EVCS
setPF = 0.95;
ineq(48).A = zeros(lent.Pevcs,lent.total);
ineq(48).A(:,inp.Qevcs) = eye(lent.Qevcs);
ineq(48).A(:,inp.Pevcs) = -eye(lent.Pevcs) .* tan(acos(setPF));

ineq(48).b = zeros(lent.Pevcs,1);

%% Enable QV curve constraint
if enQV == false
    for ii = 31:44 
        ineq(ii).A = []; ineq(ii).b = [];
    end
end
if enQV == false
    for ii = 31:46
        equ(ii).Aeq = []; equ(ii).beq = [];
    end
end

%% Enable EV Uncertainty
if enEV == false % delete old constraints
    for ii = 10:16
        ineq(ii).A = []; ineq(ii).b = [];
    end
else
    for ii = 7:9
        equ(ii).Aeq = []; equ(ii).beq = [];
        ineq(ii).A = []; ineq(ii).b = [];
    end
end

%% Enable PV Uncertainty (NEW VER)
if enPV == false
    for ii = 23:24
        ineq(ii).A = []; ineq(ii).b = []; % delete the uncertainty constraint
    end
else
    equ(18).Aeq = []; equ(18).beq = []; % delete the original constraint
    ineq(19).A = []; ineq(19).b = []; % delete the original constraint
end

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
% options.mip.tolerances.mipgap = 0.5/100;
% options.timelimit = 300;

%% Execute MIQP
[xx, fval, exitflag, output] = cplexmiqp (H, f, Aineq, bineq, Aeq, beq,...
    [], [], [], lb, ub, ctype, [], options);
toc
%% Extract result
res = extractresult(xx,len,inp,ts);

% Plot the QV Curve
Vmins = 0.9+0.01*res.cmin;
Vmaxs = 0.9+0.01*res.cmax;
figure;
for ii = 1:length(Vmins)
    xx = [0.8 0.9 Vmins(ii) Vmaxs(ii) 1.1 1.2];
    yy = [1 1 0 0 -1 -1];
    pl(ii) = plot(xx,yy,'LineWidth',2);
    hold on
    ylim([-1.5 1.5])
end
% Scatter Qstar in every V
Vscat = res.V(:,data.evcs(:,2)');
Qscat = res.Qstar;
for ii = 1:size(res.Qstar,2)
    sc(ii) = scatter(Vscat(:,ii),Qscat(:,ii));
    hold on
end

title('QV control curve and reactive power dispatch of each EVCS')
xlabel('Voltage')
ylabel('Q^*')
VV = [Vmins Vmaxs];
for ii = 1:6
    xl(ii) = xline(VV(ii),':');
end
legend([pl sc],'EVCS1','EVCS2','EVCS3','Q^* EVCS1','Q^* EVCS2','Q^* EVCS3')
hold off
result{1,uu} = res;
result{2,uu} = fval;
result{3,uu} = output;
end

%!!! NOTICE : parameter in constraint 30 originally from evcs(:,6) changed
%into evcs(:,9) but has not been computed yet