# Numerical changes

These fixes change numerical results while retaining the physical equations,
historical constants and public calling signatures. Recompute scenario results
and check convergence when comparing with an older METIS revision.

## Stable positive root for the two-point SOL model

**Sensitivity: potentially large for cold targets and degenerate limits.**

`z0twopoints` solves `tau^3 + C*tau = D`, with `C,D >= 0`. Subtracting
nearly equal terms in Cardano's expression can lose significant digits. The
positive root is now scaled by the smaller of its cubic and linear upper
bounds and evaluated by monotone Newton iteration.

The coefficient arrays retain the original requirement to match the heat-flux
size. Shape mismatches raise `METIS:twoPointShape` before solving; there is no
implicit expansion of a scalar heat flux. Invalid coefficient points propagate
through `indbad`, leaving valid points in the same call intact.

The conductive-balance diagnostic uses a finite equation scale for `C=0`.
Valid isothermal solutions no longer receive an invalid-point flag because
of division by zero in that diagnostic.

```matlab
assertSuccess(runtests('test/test_sol_positive_root.m'));
```

The six regressions check cold-target values calculated at 80-digit precision,
a manufactured warm-target solution, zero conduction, rejected coefficient
shape mismatches, scalar/batch agreement, and mixed valid/invalid points.
The cold-target inputs and expected values are included in the test file.
No external scenario data or additional solver is required.

### Reconstruction of positive cold roots

**Sensitivity: high at the cold end of the SOL heat-flux search.**

The temperature root must also satisfy conduction, pressure and sheath balances.
The previous density reconstruction divided the sheath coefficient by
`max(eps,T_target^(3/2))`. This dimensionless cutoff changed the balances for
positive cold roots. A second cutoff replaced targets below `eps` by `0.1 eV`.
For example, a manufactured `200 eV` upstream solution was reconstructed as
about `0.000674 eV` before these cutoffs were removed.

The solver now reconstructs upstream temperature from
`T_upstream = (C + T_target^(7/2))^(2/7)` and density from the pressure balance.
It retains positive target temperatures without a machine-epsilon floor in
physical units. The original invalid-coefficient and invalid-root diagnostics
remain active. A seventh SOL regression checks all three balances, including
positive targets smaller than machine epsilon. These are algebraic limits;
they do not qualify a plasma model at arbitrarily low temperatures or high
densities. Reassess coupled SOL convergence when using the changed outputs.

## Heat-profile fixed-point convergence

**Sensitivity: high; temperatures, exchange and transport outputs can change.**

`zeqie0prof` previously reduced its update fraction by a factor of 0.7 on
every iteration and tested the damped increment. The sum of these fractions
is finite, so that increment could become small without the temperature map
reaching its fixed point.

The stopping residual now measures the unrelaxed electron/ion temperature
maps with the existing mean-temperature normalization and `1e-3` tolerance.
The update fraction stays at its original initial value of 0.3. The minimum
12 iterations and maximum 101 iterations are unchanged. An unresolved or
nonfinite residual at the budget limit emits `METIS:thermalNonConvergence`,
including the iteration budget and residual. A warning does not establish
convergence of the returned candidate.

```matlab
assertSuccess(runtests('test/test_thermal_fixed_point.m'));
```

This regression constructs a symmetric heat balance with constant density,
identical electron/ion sources and transport shapes, and a known quadratic
temperature fixed point. It checks temperatures against that fixed point,
thermal energy, and the boundary temperature on the solver's stated accuracy
scale. All inputs are constructed in the test itself.

These tests verify numerical equations and stopping behavior. Full scenario
convergence and agreement with measurements require separate validation.
