function [inp] = inputvar(lent)

inp.Pg = 1:lent.Pg;
inp.Qg = inp.Pg(end) + (1:lent.Qg);
inp.Pline = inp.Qg(end) + (1:lent.Pline);
inp.Qline = inp.Pline(end) + (1:lent.Qline);
inp.Pbessc = inp.Qline(end) + (1:lent.Pbessc);
inp.Pbessd = inp.Pbessc(end) + (1:lent.Pbessd);
inp.Pevcs = inp.Pbessd(end) + (1:lent.Pevcs);
inp.Qevcs = inp.Pevcs(end) + (1:lent.Qevcs);
inp.Qstar = inp.Qevcs(end) + (1:lent.Qstar);
inp.Qcb = inp.Qstar(end) + (1:lent.Qcb);
inp.PEV = inp.Qcb(end) + (1:lent.PEV);
inp.SocEV = inp.PEV(end) + (1:lent.SocEV);
inp.SocBess = inp.SocEV(end) + (1:lent.SocBess);
inp.V = inp.SocBess(end) + (1:lent.V);
inp.a = inp.V(end) + (1:lent.a);
inp.wmin = inp.a(end) + (1:lent.wmin);
inp.wmax = inp.wmin(end) + (1:lent.wmax);
inp.smin = inp.wmax(end) + (1:lent.smin);
inp.smax = inp.smin(end) + (1:lent.smax);

inp.Tap = inp.smax(end) + (1:lent.Tap);
inp.bcb = inp.Tap(end) + (1:lent.bcb);
inp.bBess = inp.bcb(end) + (1:lent.bBess);
inp.cmin = inp.bBess(end) + (1:lent.cmin);
inp.cmax = inp.cmin(end) + (1:lent.cmax);
inp.lmin = inp.cmax(end) + (1:lent.lmin);
inp.lmax = inp.lmin(end) + (1:lent.lmax);
inp.z = inp.lmax(end) + (1:lent.z);
inp.d = inp.z(end) + (1:lent.d);

inp.cont = 1:inp.smax(end);
inp.bin = inp.Tap(1):inp.d(end);

end