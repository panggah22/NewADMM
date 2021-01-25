function [len, lent] = lengthvar(data,evr,tt)

len.Pg = data.num_ficgen;
len.Qg = data.num_ficgen;
len.Pline = data.num_branch;
len.Qline = data.num_branch;
len.Pbessc = data.num_evcs;
len.Qbessc = data.num_evcs;
len.Pbessd = data.num_evcs;
len.Qbessd = data.num_evcs;
len.Pevcs = data.num_evcs;
len.Qevcs = data.num_evcs;
len.Ppv = data.num_evcs;
len.Qstar = data.num_evcs;
len.Qcb = data.num_cb;
len.PEV = data.num_ev;
len.SocEV = data.num_ev;
len.SocBess = data.num_evcs;
len.V = data.num_bus;
len.a = 6*data.num_evcs;

len.wmin = 5*data.num_evcs;
len.wmax = 5*data.num_evcs;
len.smin = data.num_evcs;
len.smax = data.num_evcs;

len.Tap = 1;
len.bcb = data.num_cb; 
len.bBess = data.num_evcs;
len.cmin = data.num_evcs; % integer 0~20
len.cmax = data.num_evcs; % integer 0~20
len.lmin = 5*len.cmin; % binary
len.lmax = 5*len.cmax; % binary
len.z = data.num_ev;
len.d = 5*data.num_evcs; % binary


lent.Pg = len.Pg*tt;
lent.Qg = len.Qg*tt;
lent.Pline = len.Pline*tt;
lent.Qline = len.Qline*tt;
lent.Pbessc = len.Pbessc*tt;
lent.Qbessc = len.Qbessc*tt;
lent.Pbessd = len.Qbessd*tt;
lent.Qbessd = len.Qbessd*tt;
lent.Pevcs = len.Pevcs*tt;
lent.Qevcs = len.Qevcs*tt;
lent.Ppv = len.Ppv*tt;
lent.Qstar = len.Qstar*tt;
lent.Qcb = len.Qcb*tt;
lent.PEV = len.PEV*tt;
lent.SocEV = len.SocEV*(tt+1);
lent.SocBess = len.SocBess*(tt+1);
lent.V = len.V*tt;
lent.a = len.a*tt;

lent.wmin = len.wmin*tt;
lent.wmax = len.wmax*tt;
lent.smin = len.smin*tt;
lent.smax = len.smax*tt;

lent.Tap = len.Tap*tt; % integer -8~8
lent.bcb = len.bcb*tt; % binary
lent.bBess = len.bBess*tt; % binary
lent.cmin = len.cmin; % integer 0~20
lent.cmax = len.cmax; % binary
lent.lmin = len.lmin; % binary
lent.lmax = len.lmax; % binary
lent.z = data.num_ev * evr; % binary
lent.d = len.d*tt; % binary


lent.total = sum(cell2mat(struct2cell(lent)));




end