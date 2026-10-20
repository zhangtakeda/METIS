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

function testAutomaticMomentumMatchesTheSolvedPressureBalance(testCase)
[o,c,g,z,p]=inputs;
o.fmom=0; z.nebord=1e20*ones(2,1); z.pin=[3e5;7e5];
prior_rng=rng; cleanup=onCleanup(@() rng(prior_rng));
for conduction=[1,-1]
    for mach_enabled=[0,1]
        o.fcond=conduction; o.mach_corr=mach_enabled;
        rng(7,'twister');
        [tu,nt,tt,q,err,nb,bad,fmom,rad,neutral,incoming,power,zeff,gamma,mach]= ...
            z0convergence_2points_dic(o,c,g,z,p,4);
        % The Mach transform changes the target momentum by 2/(1+M^2).
        pressure=2*(1+mach.^2).*nt.*tt;
        expected=fmom.*(1+z.tibord./z.tebord).*z.nebord.*tu;
        verifyEqual(testCase,pressure,expected,'RelTol',2e-12);
        verifyEmpty(testCase,bad);
    end
end
end

function testCandidateRankingUsesItsOwnIncomingPower(testCase)
% Manufactured wall feedback changes incoming power at every trial flux.
[o,c,g,z,p]=wall_inputs;
directory=tempname; mkdir(directory); prior_path=path;
cleanup=onCleanup(@() remove_helpers(directory,prior_path));
write_helper(fullfile(directory,'z0acctungsten.m'),sprintf([ ...
    'function [nw,nwm,z1,z2]=z0acctungsten(o,g,c,z,p)\n' ...
    'nw=zeros(size(p.nwp)); nwm=zeros(size(z.nem)); z1=nw; z2=nw;\nend\n']));
write_helper(fullfile(directory,'zrad0.m'),sprintf([ ...
    'function [prad,pbrem,pcyclo,pradsol,p]=zrad0(varargin)\n' ...
    'persistent initial;\n' ...
    'tt=varargin{22}; nt=varargin{23}; e=1.602176462e-19; ua=1.66053873e-27;\n' ...
    'gamma=2.5+2/.9-.5*log((2*pi*9.10938188e-31/1.6726485e-27)*2/.9^2);\n' ...
    'q=gamma*e.*nt.*tt.*sqrt(2*e.*tt/(2*ua));\n' ...
    'if isempty(initial); initial=q; end\n' ...
    'fraction=.2*ones(size(q)); fraction(q>=.999*initial)=.99;\n' ...
    'prad=[3e6;7e6].*(1-fraction.*q./initial);\n' ...
    'pbrem=zeros(size(q)); pcyclo=pbrem; pradsol=pbrem; p=struct;\nend\n']));
addpath(directory,'-begin'); clear z0acctungsten zrad0;
prior_rng=rng; rng_cleanup=onCleanup(@() rng(prior_rng)); rng(7,'twister');
[tu,nt,tt,q,err,nb,bad,fmom,rad,neutral,incoming]= ...
    z0convergence_2points_dic(o,c,g,z,p,4);
% The initial candidate has residual 1/99; every lower-flux candidate has 4.
verifyEqual(testCase,incoming,.99*q,'RelTol',2e-12);
verifyEqual(testCase,err,1/99,'RelTol',2e-11);
verifyEmpty(testCase,bad);
end

function [o,c,g,z,p]=wall_inputs
[o,c,g,z,p]=inputs;
o.W_effect=1; o.mach_corr=0; o.Recycling=0; o.configuration=3;
o.modeh=1; o.frad=1; o.matthews=0; o.z_prad='zmax';
o.sol_rad='coupled'; o.gaunt=0; o.noncoronal=0; o.fnesol=1;
o.te_max=1e5; o.cor_rel_brem='Stott'; z.xpoint=ones(2,1);
g.K=ones(2,1); g.b0=2.5*ones(2,1);
z.nem=3e19*ones(2,1); z.tem=1e3*ones(2,1); z.zmszl=ones(2,1);
z.ane=ones(2,1); z.ate=ones(2,1); z.te0=1e3*ones(2,1);
z.nbar=z.nem; z.taue=ones(2,1); z.vp=ones(2,1); z.sext=ones(2,1);
end

function write_helper(filename,source)
fid=fopen(filename,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',source);
end

function remove_helpers(directory,prior_path)
path(prior_path); clear z0acctungsten zrad0;
rmdir(directory,'s');
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
