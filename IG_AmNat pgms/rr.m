function profilex
tic
addpath /usr/local/packages/oclJacket/Lnx64/engine
profile on
IG_AmNat;
profile off
profsave(profile('info'),datestr(now))
toc