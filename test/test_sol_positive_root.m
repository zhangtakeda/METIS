function tests = test_sol_positive_root
% Regression for cold SOL targets and the unchanged public call signature.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.prior_path = path;
addpath(fullfile(root,'zerod'));
end

function teardownOnce(testCase)
path(testCase.TestData.prior_path);
end

function testColdTargetAgainstIndependent80DigitRoot(testCase)
n = 31;
q = logspace(2,7,n)';
ne = logspace(18,21,n)';
lc = linspace(5,200,n)';
zeff = linspace(1,4,n)';
mass = linspace(1,4,n)';
fpe = linspace(.1,1,n)';
fei = linspace(1,3,n)';
tite = linspace(.5,3,n)';
gamma = linspace(5,9,n)';
[tu,tt,nt,bad] = z0twopoints(q,ne,lc,zeff,mass,fpe,fei,tite,gamma);
% Decimal Newton solution of tau^3+C*tau=(A*B)^(7/2), at 80 digits.
index = [23,29,30,31];
expected_tu = [36.44972752021934;83.9712836175745;96.18582072573463;110.09075936542457];
expected_tt = [6.29815757561871e-5;4.531552425450344e-5;4.344416535752357e-5;4.1780211974907296e-5];
expected_nt = [5.6562843970624005e25;8.37916974667876e26;1.289677203052111e27;1.9762482194599166e27];
verifyEqual(testCase,tu(index),expected_tu,'RelTol',2e-13);
verifyEqual(testCase,tt(index),expected_tt,'RelTol',2e-13);
verifyEqual(testCase,nt(index),expected_nt,'RelTol',2e-13);
verifyEmpty(testCase,bad);
end

function testManufacturedWarmTarget(testCase)
% Construct coefficients from positive target values and the three balances.
e = 1.602176462e-19; ua = 1.66053873e-27;
tt = [1;10;100]; tu = [50;100;500]; nt = [1e20;2e20;3e20];
zeff = [1;2;3]; mass = [2;2.5;3]; tite = [1;1;2];
gamma = 7*ones(3,1); fei = 2*ones(3,1); fpe = .5*ones(3,1);
ne = 4*nt.*tt./(fei.*tu);
q = gamma*e.*nt.*tt.*sqrt((zeff+tite)*e.*tt./(ua*mass));
lc = (tu.^(7/2)-tt.^(7/2)).*(2000./zeff)./((7/2)*fpe.*q);
[actual_tu,actual_tt,actual_nt,bad] = z0twopoints(q,ne,lc,zeff,mass,fpe,fei,tite,gamma);
verifyEqual(testCase,actual_tu,tu,'RelTol',2e-13);
verifyEqual(testCase,actual_tt,tt,'RelTol',2e-13);
verifyEqual(testCase,actual_nt,nt,'RelTol',2e-13);
verifyEmpty(testCase,bad);
end

function testZeroConductiveCoefficient(testCase)
q = [1e5;1e7]; ne = [1e19;1e20]; v = ones(2,1);
[tu,tt,nt,bad] = z0twopoints(q,ne,10*v,v,2*v,0*v,2*v,v,7*v);
verifyEqual(testCase,tu,tt,'RelTol',2e-13);
verifyTrue(testCase,all(isfinite(nt) & nt>0));
verifyEmpty(testCase,bad);
end

function testManufacturedColdBalanceWithoutDimensionalEpsilonFloor(testCase)
% A cold positive root still satisfies pressure, sheath and conduction.
% These algebraic limits do not establish a cold-plasma validity domain.
e=1.602176462e-19; ua=1.66053873e-27;
q=ones(3,1); ne=[1e19;1e20;1e21]; tu=200*ones(3,1);
lc=100*ones(3,1); v=ones(3,1); gamma=7*v; fei=2*v;
B=q./gamma/e./sqrt(2*e/(2*ua));
tt=(4*B./(fei.*ne.*tu)).^2;
nt=fei.*ne.*tu./(4*tt);
fpe=(tu.^(7/2)-tt.^(7/2))*2000./((7/2)*q.*lc);
[actual_tu,actual_tt,actual_nt,bad]=z0twopoints(q,ne,lc,v,2*v,fpe,fei,v,gamma);
verifyEqual(testCase,actual_tu,tu,'RelTol',2e-13);
verifyEqual(testCase,actual_tt,tt,'RelTol',2e-13);
verifyEqual(testCase,actual_nt,nt,'RelTol',2e-13);
verifyEmpty(testCase,bad);
verifyEqual(testCase,4*actual_nt.*actual_tt,fei.*ne.*actual_tu,'RelTol',2e-13);
verifyEqual(testCase,gamma*e.*actual_nt.*actual_tt.*sqrt(2*e*actual_tt/(2*ua)),q,'RelTol',2e-13);
verifyEqual(testCase,actual_tu.^(7/2)-actual_tt.^(7/2), ...
    (7/2)*fpe.*q.*lc/2000,'RelTol',2e-13);
end

function testMismatchedCoefficientShapesAreRejected(testCase)
% Preserve the equal-size coefficient contract of the original polynomial solve.
verifyError(testCase,@() z0twopoints(1e6,[1e19;2e19],10,1,2,.5,2,1,7), ...
    'METIS:twoPointShape');
verifyError(testCase,@() z0twopoints(1e6,1e19,[10;20],1,2,.5,2,1,7), ...
    'METIS:twoPointShape');
verifyError(testCase,@() z0twopoints(1e6,[1e19;2e19],[10;20],1,2,.5,2,1,7), ...
    'METIS:twoPointShape');
end

function testExplicitVectorInputsSolveEveryPoint(testCase)
v=ones(2,1);
[tu,tt,nt,bad]=z0twopoints(1e6*v,[1e19;2e19],10*v,v,2*v,.5*v,2*v,v,7*v);
verifySize(testCase,tu,[2,1]);
verifySize(testCase,tt,[2,1]);
verifySize(testCase,nt,[2,1]);
% Solving a batch must give the same physical values as each scalar call.
for k=1:2
    [single_tu,single_tt,single_nt,single_bad]= ...
        z0twopoints(1e6,k*1e19,10,1,2,.5,2,1,7);
    verifyEqual(testCase,[tu(k),tt(k),nt(k)], ...
        [single_tu,single_tt,single_nt],'RelTol',2e-13);
    verifyEmpty(testCase,single_bad);
end
verifyEmpty(testCase,bad);
end

function testInvalidPointsDoNotAbortOrChangeValidPoints(testCase)
v=ones(2,1);
[good_tu,good_tt,good_nt,good_bad]=z0twopoints(1e6,1e19,10,1,2,.5,2,1,7);
verifyEmpty(testCase,good_bad);
for invalid_value=[NaN,Inf,-.5]
    % Exercise invalid conductive and forcing coefficients independently.
    for invalid_position=1:2
        fraction=.5*v; pressure=2*v;
        if invalid_position==1
            fraction(2)=invalid_value;
        else
            pressure(2)=invalid_value;
        end
        [tu,tt,nt,bad]=z0twopoints(1e6*v,1e19*v,10*v,v,2*v,fraction,pressure,v,7*v);
        verifyEqual(testCase,[tu(1),tt(1),nt(1)], ...
            [good_tu,good_tt,good_nt],'RelTol',2e-13);
        verifyEqual(testCase,bad,2);
        verifySize(testCase,tu,[2,1]);
        verifySize(testCase,tt,[2,1]);
        verifySize(testCase,nt,[2,1]);
    end
end
end
