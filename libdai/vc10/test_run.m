fg = dai_readfg('../tests/alarm.fg');
[logz, q, md, qv, qf, qmap] = dai(fg, 'BP', '[inference=SUMPROD,updates=SEQMAX,tol=1e-6,maxiter=100,logdomain=0]');
