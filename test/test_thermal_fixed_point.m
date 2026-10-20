function tests = test_thermal_fixed_point
% A symmetric heat balance must reach its known electron/ion fixed point.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.prior_path = path;
addpath(fullfile(root,'zerod'));
addpath(fullfile(root,'solver'));
end

function teardownOnce(testCase)
path(testCase.TestData.prior_path);
end

function testSymmetricManufacturedTemperature(testCase)
x = linspace(0,1,21);
ne = 1e19*ones(size(x));
expected = 100+2000*(1-x.^2);
p.xli=x; p.vpr=100*x; p.grho2=ones(size(x)); p.rmx=x;
p.Raxe=3*ones(size(x)); p.ftrap=.6*ones(size(x));
p.nep=ne; p.nip=ne; p.n1p=ne;
p.nhep=0*ne; p.nzp=0*ne; p.nwp=0*ne;
p.tep=100+1.2*(expected-100); p.tip=100+.8*(expected-100);
p.source_el=1e5*ones(size(x)); p.source_ion=p.source_el;
p.xieshape=ones(size(x)); p.xieshape_itb=p.xieshape;
p.xie=ones(size(x)); p.xii=p.xie;
energy=3*1.602176462e-19*trapz(x,ne.*expected.*p.vpr);
[unused1,unused2,unused3,unused4,unused5,out]=zeqie0prof( ...
    trapz(x,p.vpr),1,1,p.tep(1),1e19,1,1,energy,1,1,1,2,0,1,0, ...
    1e19,0,1e19,0,0,6,0,18,0,100,1e19,1,0,1,.95,0,p,1, ...
    'bgb',0,0,3,0,1e5,'off',0,0,13.6);
% The solver stops at a 1e-3 residual normalized by the mean temperature.
% Compare the manufactured fields on their common temperature scale.
verifyLessThan(testCase,max(abs(out.tep-expected))/max(expected),2e-3);
verifyLessThan(testCase,max(abs(out.tip-expected))/max(expected),2e-3);
actual_energy=1.5*1.602176462e-19*trapz(x,ne.*(out.tep+out.tip).*p.vpr);
verifyEqual(testCase,actual_energy,energy,'RelTol',1e-3);
verifyEqual(testCase,out.tep(end),100,'AbsTol',1e-12);
verifyLessThan(testCase,abs(out.tip(end)-100)/max(expected),2e-3);
end
