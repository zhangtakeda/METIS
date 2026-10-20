function tests = test_sol_candidate
% Returned SOL fields and residual must describe one evaluated candidate.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.prior_path=path;
addpath(fullfile(root,'zerod'));
addpath(fullfile(root,'solver'));
end

function teardownOnce(testCase)
path(testCase.TestData.prior_path);
end

function testReturnedResidualUsesTheReturnedFields(testCase)
[o,c,g,z,p]=inputs;
prior_rng=rng; cleanup=onCleanup(@() rng(prior_rng));
rng(7,'twister');
[tu,nt,tt,q,err,nb,bad,fmom,rad,neutral,incoming,power,zeff,gamma,mach]= ...
    z0convergence_2points_dic(o,c,g,z,p,4);
gamma_e=gamma-2.5*(z.tibord./z.tebord);
corrected=min(incoming-rad,(incoming-rad)./(1+o.Recycling*o.eioniz./gamma_e./max(.1,tt)));
verifyEqual(testCase,err,max(abs(q-corrected)./max(1,abs(incoming))),'RelTol',2e-13);
verifyEqual(testCase,neutral,incoming-min(incoming,corrected+rad),'RelTol',2e-13);
verifyEmpty(testCase,bad);
% Undo the returned Mach transform before checking the sheath equation.
factor=(gamma_e.*(1+mach.^2)./(2*abs(mach).*(gamma_e-1+mach.^2))).^2;
base_tt=tt./factor; base_nt=nt.*factor.*(1+mach.^2)/2;
e=1.602176462e-19; ua=1.66053873e-27;
sheath=gamma*e.*base_nt.*base_tt.*sqrt((z.zeff+z.tibord./z.tebord)*e.*base_tt./(ua*z.meff));
verifyEqual(testCase,sheath,q,'RelTol',2e-12);
end

function testInvalidPointFlagSurvivesSearchCompletion(testCase)
[o,c,g,z,p]=inputs;
o.mach_corr=0; o.Recycling=0;
z.nebord(2)=Inf;
[tu,nt,tt,q,err,nb,bad]=z0convergence_2points_dic(o,c,g,z,p,0);
verifyEqual(testCase,bad,2);
end

function [o,c,g,z,p]=inputs
% Construct a two-time, circular, radiation-free example from the equations.
x=linspace(0,1,21); v=ones(2,1);
g.R=2.5*v; g.a=.8*v;
c.temps=[0;1]; c.iso=zeros(2,1);
o.sol_model='2_points'; o.W_effect=0; o.configuration=1;
o.Recycling_target=0; o.Recycling=.8; o.eioniz=30; o.detach=0;
o.lcx=1; o.Sq=0; o.fpower=.5; o.imp_div='fixed';
o.zimp=6; o.zmax=18; o.rimp=0; o.fzmax_div=0;
o.carbonblow=0; o.gaz=2; o.frhe0=0; o.Sn_fraction=0;
o.de=.1; o.fmom=.5; o.fcond=1; o.alpha_e=.82; o.mach_corr=1;
z.pin=[3e6;7e6]; z.pioniz=0*v; z.pcyclo=0*v;
z.prad=0*v; z.pbrem=0*v; z.pradsol=0*v;
z.qa=3*v; z.peri=2*pi*g.a; z.xpoint=0*v; z.dsol=.01*v;
z.tebord=100*v; z.tibord=100*v; z.tite=v;
z.nebord=3e19*v; z.zeff=v; z.meff=2*v;
z.pel=.6*z.pin; z.pion=.4*z.pin;
p.xli=x; p.Raxe=g.R*ones(size(x)); p.rmx=g.a*x;
p.fdia=6.25*ones(2,numel(x)); p.psi=-.2*v*x.^2;
p.bpol=.2*ones(2,numel(x)); p.vpr_tor=10*v*x;
p.nep=3e19*ones(2,numel(x)); p.zeff=ones(2,numel(x));
p.nwp=zeros(2,numel(x)); p.tip=100*ones(2,numel(x));
end
