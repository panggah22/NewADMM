function data = data33(tt,unc,pk,fic)
%% This provides data for the whole formulation
% tt  : time step for the simulation
% unc : value of uncertainty
% pk  : list of buses in an area
% fic : fictitious generator
% All conversions and data selection should be in consecutive steps:
% A. select buses in an area in each data section
% B. counting data size and internal numbering
% C. convert the data into per unit value
%
% 1. Bus data contains the bus ID and load demands
% 2. Line data contains the connection between two bus IDs as well as
% resistance and reactance


data.kVbase = 12.66;
data.MVAbase = 1;
data.zbase = (data.kVbase^2)/data.MVAbase;

if size(pk,1) == 1 % pk must be a column vector
    pk = pk';
end

%% Bus data 
data.bus = [
 % BusID    P(kW)   Q(kW)  
    1       0       0	
    2       100     60	
    3       90      40	
    4       120     80	
    5       60      30	
    6       60      20	
    7       200     100	
    8       200     100	
    9       60      20	
    10      60      20	
    11      45      30	
    12      60      35	
    13      60      35	
    14      120     80	
    15      60      10	
    16      60      20	
    17      60      20	
    18      90      40	
    19      90      40	
    20      90      40	
    21      90      40	
    22      90      40	
    23      90      50	
    24      420     200	
    25      420     200	
    26      60      25	
    27      60      25	
    28      60      20	
    29      120     70	
    30      200     600	
    31      150     70	
    32      210     100	
    33      60      40	

    ];

% A
data.bus = data.bus(ismember(data.bus(:,1),pk,'rows'),:);
% B
data.num_bus = size(data.bus,1);
data.bus = [(1:data.num_bus)', data.bus]; % for internal numbering
% C
data.bus(:,3:4) = data.bus(:,3:4)./(data.MVAbase*1000);
% Load variation in 24hr span
data.loadcoeff = [0.556 0.487 0.452 0.443 0.405 0.462 0.542 0.652 0.712 0.872 0.904 0.906 0.912 0.915 0.915 0.872 0.861 0.855 0.832 0.724 0.667 0.657 0.623 0.584]';
data.loadcoeff = data.loadcoeff * (1/max(data.loadcoeff));

%% Branch data
data.branch = [
% From  to    r       x     Cap
    1	2	0.0922	0.0470	5064	
    2	3	0.4930	0.2511	3798	
    3	4	0.3660	0.1864	3798	
    4	5	0.3811	0.1941	3798	
    5	6	0.8190	0.7070	3798	
    6	7	0.1872	0.6188	3798	
    7	8	0.7114	0.2351	3798	
    8	9	1.0300	0.7400	3798	
    9	10	1.0440	0.7400	3798	
    10	11	0.1966	0.0650	3798	
    11	12	0.3744	0.1238	3798    
    12	13	1.4680	1.1550	3798	
    13	14	0.5416	0.7129	3798	
    14	15	0.5910	0.5260	3798	
    15	16	0.7463	0.5450	3798	
    16	17	1.2890	1.7210	3798	
    17	18	0.7320	0.5740	3798	
    2	19	0.1640	0.1565	5064	
    19	20	1.5042	1.3554	3798	
    20	21	0.4095	0.4784	3798	
    21	22	0.7089	0.9373	3798	
    3	23	0.4512	0.3083	3798	
    23	24	0.8980	0.7091	2532	
    24	25	0.8960	0.7011	2532	
    6	26	0.2030	0.1034	2532	
    26	27	0.2842	0.1447	2532	
    27	28	1.0590	0.9337	2532	
    28	29	0.8042	0.7006	2532	
    29	30	0.5075	0.2585	2532	
    30	31	0.9744	0.9630	2532	
    31	32	0.3105	0.3619	2532	
    32	33	0.3410	0.5302	2532	
%     8	21	2.0000	2.0000	2532	
%     9	15	2.0000	2.0000	2532	
%     12	22	2.0000	2.0000	2532	
%     18	33	0.5000	0.5000	2532	
%     25	29	0.5000	0.5000	2532	
    ];

% A
data.branch = data.branch(ismember(data.branch(:,1),pk,'rows')&ismember(data.branch(:,2),pk,'rows'),:);
% B
data.num_branch = size(data.branch,1);
fromto = zeros(data.num_branch,2);
for i = 1:data.num_branch % branch numbering
    fromto(i,1) = data.bus((data.bus(:,2) == data.branch(i,1)),1);
    fromto(i,2) = data.bus((data.bus(:,2) == data.branch(i,2)),1);
end
data.branch = [(1:data.num_branch)', fromto, data.branch];
% C
data.branch(:,6:7) = data.branch(:,6:7)./data.zbase;
data.branch(:,8) = data.branch(:,8)./(data.MVAbase*1000);

%% CB data
data.cb = [
%   BusID   Q(kVar)
    11       100
    23      100
    28      100
];

% A
data.cb = data.cb(ismember(data.cb(:,1),pk,'rows'),:);
% B
data.num_cb = size(data.cb,1);
cbbus = zeros(data.num_cb,1);
for ii = 1:data.num_cb
    cbbus(ii) = data.bus(data.bus(:,2)==data.cb(ii,1),1);
end
data.cb = [(1:data.num_cb)', cbbus ,data.cb];
% C
data.cb(:,4) = data.cb(:,4)./(data.MVAbase*1000);

%% Fictitious generator
if size(fic,1) == 1
    data.ficgen = fic';
else
    data.ficgen = fic;
end

% A
data.ficgen = data.ficgen(ismember(data.ficgen(:,1),pk,'rows'),:);
% B
data.num_ficgen = size(data.ficgen,1);
ficbus = zeros(data.num_ficgen,1);
for ii = 1:data.num_ficgen
    ficbus(ii) = data.bus(data.bus(:,2)==data.ficgen(ii),1);
end
data.ficgen = [(1:data.num_ficgen)', ficbus, data.ficgen];

%% EVCS data
data.evcs = [
%  Node  	Capacity   kWmax   Smax     PVkW   EV_num
    18      1000         200      400    50     10
    3     1000         200      400    50     10
    32     1000         200      400    50     10
];
Qmax = sqrt(data.evcs(:,4).^2 / 2);
data.evcs = [data.evcs Qmax];
% A 
data.evcs = data.evcs(ismember(data.evcs(:,1),pk,'rows'),:);
% B 
data.num_evcs = size(data.evcs,1);
evcsbus = zeros(data.num_evcs,1);
for ii = 1:data.num_evcs
    evcsbus(ii) = data.bus(data.bus(:,2)==data.evcs(ii,1),1);
end
data.evcs = [(1:data.num_evcs)', evcsbus, data.evcs];
% C 
data.evcs(:,[4:7, 9]) = data.evcs(:,[4:7, 9])./(data.MVAbase*1000);
% PV coefficient
data.PVcoeff = 1.25*[0 0 0 0 0 0 0.05 0.2 0.35 0.50 0.65 0.75 0.8 0.75 0.65 0.50 0.35 0.2 0.05 0 0 0 0 0];

% A
data.num_ev = sum(data.evcs(:,8));
%% Driving pattern realization
%           kWh     Pcmax   kWh/km  Socmin  Socmax
data.evv = [60      10      0.15    0.2     1];
data.evv(:,1:3) = data.evv(:,1:3)./(data.MVAbase*1000);

data.rea = [
    1       4       23      140     0.00169
    2       4       24      140     0.00003
    3       5       23      140     0.71116
    4       5       24      140     0.03383
    5       4       23      150     0.00008
    6       5       23      150     0.00671
    7       5       24      150     0.00049
    8       4       23      160     0.00032
    9       5       23      160     0.04367
    10      5       24      160     0.00506
    ];
    
t_travel = data.rea(:,3) - data.rea(:,2) - 1;
consum = data.rea(:,4) .* data.evv(3);
kW_each_hour = consum./t_travel;

ddd = cell(size(data.rea,1),1);
for ii = 1:length(ddd)
    ak = zeros(1,tt);
    ak(data.rea(ii,2):(data.rea(ii,3)-1)) = kW_each_hour(ii);
    ddd{ii} = ak;
end

data.dr = ddd;

data.k1 = data.dr(data.rea(:,5) > (1-unc));
data.k2 = data.dr(data.rea(:,5) <= (1-unc));


