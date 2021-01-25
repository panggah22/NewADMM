function res = extractresult(len,inp,ts)
%% Extract result
res.Pg = reshape(xx(inp.Pg),len.Pg,ts)';
res.Qg = reshape(xx(inp.Qg),len.Qg,ts)';
res.Pline = reshape(xx(inp.Pline),len.Pline,ts)';
res.Qline = reshape(xx(inp.Qline),len.Qline,ts)';
res.Pbessc = reshape(xx(inp.Pbessc),len.Pbessc,ts)';
res.Qbessc = reshape(xx(inp.Qbessc),len.Qbessc,ts)';
res.Pbessd = reshape(xx(inp.Pbessd),len.Pbessd,ts)';
res.Qbessd = reshape(xx(inp.Qbessd),len.Qbessd,ts)';
res.Pevcs = reshape(xx(inp.Pevcs),len.Pevcs,ts)';
res.Qevcs = reshape(xx(inp.Qevcs),len.Qevcs,ts)';
res.Qstar = reshape(xx(inp.Qstar),len.Qstar,ts)';
res.Qcb = reshape(xx(inp.Qcb),len.Qcb,ts)';
res.PEV = reshape(xx(inp.PEV),len.PEV,ts)';
res.SocEV = reshape(xx(inp.SocEV),len.SocEV,ts+1)';
res.SocBess = reshape(xx(inp.SocBess),len.SocBess,ts+1)';
res.V = reshape(xx(inp.V),len.V,ts)';
res.a = reshape(xx(inp.a),len.a,ts)';
res.wmin = reshape(xx(inp.wmin),len.wmin,ts)';
res.wmax = reshape(xx(inp.wmax),len.wmax,ts)';
res.smin = reshape(xx(inp.smin),len.smin,ts)';
res.smax = reshape(xx(inp.smax),len.smax,ts)';

res.Tap = reshape(xx(inp.Tap),len.Tap,ts)';
res.bcb = reshape(xx(inp.bcb),len.bcb,ts)';
res.bBess = reshape(xx(inp.bBess),len.bBess,ts)';
res.cmin = xx(inp.cmin)';
res.cmax = xx(inp.cmax)';
res.lmin = xx(inp.lmin)';
res.lmax = xx(inp.lmax)';
res.z = xx(inp.z)';
res.d = reshape(xx(inp.d),len.d,ts)';
end