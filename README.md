# on-acid

Acid–base equilibrium without solving a quadratic by hand.

The module is one file, [`dissociation.py`](dissociation.py), and it Depends on nothing but the standard library.
Its origin is a shortcut found in a university chemistry course,
where the assigned method Wanted the quadratic worked out for every exercise.

## The Doors

```
./build.sh    ./run.sh     # unix
build.cmd     run.cmd      # windows
```

`build.sh` Fetches its own tooling into `.venv`, then Runs lint, types, tests and the wheel.
It Refuses to build what does not pass, and a refusal Names the line that works.
The lint Catches bugs and never taste: the aligned tables are Deliberate, so no formatter Reflows them.

The `.cmd` files are Shims. They name, call and Propagate the exit code, and nothing else;
`run-bash.cmd` Finds the interpreter once, for both of them.
They are read from here and Proven only over there, so until someone Runs them on Windows,
call them Unverified.

## The Coefficients never Change

A weak monoprotic acid always Yields the same quadratic, so it is worth solving once and for all.
With `x = [H+]`, the [acid dissociation constant](https://en.wikipedia.org/wiki/Acid_dissociation_constant)
`Ka = x² / (C − x)` Rearranges to

```
x² + Ka·x − Ka·C = 0
```

and the coefficients are `a = 1`, `b = Ka`, `c = −Ka·C` for every acid and every concentration.
Nothing there Depends on which exercise you were given.
Solving it symbolically Leaves a direct operation, exact, with no `C − x ≈ C` approximation:

```
[H+] = (Ka/2) · (√(1 + 4C/Ka) − 1)
```

`monoprotic_h` Ships the rationalized form of that same expression,
`2·Ka·C / (Ka + √(Ka² + 4·Ka·C))`,
because subtracting two nearly equal numbers Costs significant digits when `Ka ≪ C`
— see [loss of significance](https://en.wikipedia.org/wiki/Loss_of_significance).

## Many Constants Become a Sum

A polyprotic acid Generalizes the closed form into a sum over its species.
Cumulative constants `βⱼ = Ka₁·Ka₂···Kaⱼ` Give each species a fraction `αⱼ`,
and the proton balance Closes the system:

```
h = Σⱼ j·αⱼ·C + Kw/h
```

`polyprotic_h` Solves that balance by [bisection](https://en.wikipedia.org/wiki/Bisection_method)
on a logarithmic scale, which is monotone and always converges.
The [self-ionization of water](https://en.wikipedia.org/wiki/Self-ionization_of_water) Enters the balance here,
so extreme dilutions Approach pH 7 instead of drifting past it.
`fractions` Exposes the αⱼ on their own, for a distribution diagram.

## A Base Closes against pKw

A weak base Runs the same algebra on `[OH−]` and finishes with `pH = pKw − pOH`.
`base_ph` Keeps `kw` as a parameter rather than hardcoding the 14,
so a temperature other than 25 °C Costs one argument.

## The Table Answers by Name

Twenty-three common species Carry their constants, so a calculation Needs no lookup elsewhere.
`ACIDS` Holds fifteen acids and `BASES` Holds eight bases, all at 25 °C,
rounded to two significant figures from the standard tables
(see the [CRC Handbook of Chemistry and Physics](https://en.wikipedia.org/wiki/CRC_Handbook_of_Chemistry_and_Physics)).
`find_species` Matches on formula, name or alias and Ignores accents and case.
Spanish names Survive as aliases, so `"acetic acid"`, `"CH3COOH"`, `"ácido acético"`, `"acetico"`,
`"vinegar"` and `"vinagre"` all Reach the same entry.

## A Buffer is the Same Balance

A conjugate pair Adds one spectator cation, and nothing else about the balance changes.
`buffer_ph` Feeds `base_c` in as that cation and Solves exactly,
so the answer Stays right where
[Henderson–Hasselbalch](https://en.wikipedia.org/wiki/Henderson%E2%80%93Hasselbalch_equation) Drifts.
Both Ship, because the approximation is the one people Quote:
a 1:1 acetate pair at 0.1 M Gives 4.745 either way, and the same pair at 1×10⁻⁵ M Gives 5.268 exactly
against 4.745 approximated.

## One Balance Draws the Whole Curve

A titration Needs no piecewise treatment, because dilution and the strong titrant Fit inside the same equation.
`titration_ph` Dilutes analyte and titrant by the combined volume and Enters the strong base as its cation,
so the buffer region, the equivalence point and the excess-titrant tail all Come from one solve.
`titration_curve` Samples it, `equivalence_volume` Locates the jumps,
and volumes Carry whatever unit you pass in.

Acetic acid, 0.1 M in 50 mL, against 0.1 M NaOH:
pH 2.875 at the start, 4.745 at half-equivalence, 8.722 at equivalence, 11.959 at 60 mL.
Phosphoric acid Shows its first two `pKa` plateaus at 25 and 75 mL and Flattens on the third,
which is what a `Ka₃` of 4.2×10⁻¹³ Does in water.

## Use

```python
from dissociation import ph, acid_ph, base_ph, fractions, find_species

ph("CH3COOH", 0.1)          # 2.875 — acetic acid, 0.1 M
ph("NH3", 0.1)              # 11.125 — ammonia, 0.1 M
ph("H3PO4", 0.1)            # 1.622 — three constants, solved as one balance
ph("HCN", 1e-6, exact=True) # 6.987 — water included

acid_ph(1.8e-5, 0.1)        # a constant you bring yourself
base_ph(1.8e-5, 0.1)        # the same number read as a Kb

fractions(find_species("H3PO4").k, 1e-7)
# [0.0, 0.6173, 0.3827, 0.0] — speciation at pH 7
```

```python
from dissociation import buffer_ph, henderson_hasselbalch
from dissociation import titration_ph, titration_curve, equivalence_volume

buffer_ph(1.8e-5, 0.1, 0.1)             # 4.745 — exact
henderson_hasselbalch(1.8e-5, 0.1, 0.1) # 4.745 — the approximation agrees here
buffer_ph(1.8e-5, 1e-5, 1e-5)           # 5.268, where the approximation still says 4.745

equivalence_volume(0.1, 50.0, 0.1)          # 50.0 mL
titration_ph([1.8e-5], 0.1, 50.0, 0.1, 25.0)  # 4.745 — half-equivalence
titration_ph([1.8e-5], 0.1, 50.0, 0.1, 50.0)  # 8.722 — equivalence
titration_curve([1.8e-5], 0.1, 50.0, 0.1, points=101)  # [(volume, pH), ...]
```

The module Runs its own demo, through the door, the installed command or the file:

```
./run.sh
on-acid
python dissociation.py
```

`exact=True` Forces the full balance for a single constant.
Leave it off and the closed form Answers, which is right except in very dilute solutions
or with an acid weak enough that water Competes with it.

## Scope

Concentrations Stand in for activities, so
[activity coefficients](https://en.wikipedia.org/wiki/Activity_coefficient)
and [ionic strength](https://en.wikipedia.org/wiki/Ionic_strength) are Absent,
which is the usual textbook footing and Holds well enough below roughly 0.1 M.
A titration Assumes a strong titrant; a weak-against-weak one is not Covered.
What is here Reaches one acid or base in water, its conjugate pair and its titration curve,
and within that boundary it Solves rather than approximates.

## Install

The package Ships as `on-acid` and Imports as `dissociation`:

```
pip install .
```

`pyproject.toml` Declares a flat layout with no dependencies, so the module Stays one file at the root
and `import dissociation` Works from a checkout without installing anything.
Installing also Puts `on-acid` on the PATH, which Runs `demo()` — the same worked example,
Reachable once the file is no longer in front of you.

## Tests

`test_dissociation.py` Holds 63 tests on `unittest`, because a module that Depends on nothing
should not Need pytest to be believed. The suite Enters through the public functions only
and Spells every expectation by hand, so it never Agrees with the code it Watches.
Eight mutations were Applied and reverted; the one that Survived the first pass
(`14.0` hardcoded in place of `pKw`) Bought the test that now kills it.

```
python -m unittest discover -p 'test_*.py'
```

## State

Identifiers, prose and tabulated names are English throughout, and Spanish Survives only in the aliases.
The file was renamed from `disosiacion.py`, whose spelling Missed a `c` and disagreed with its own docstring;
that file is now Gone.
The branch Carries no commits yet.

## License

[MIT](LICENSE). Permissive on purpose: research code Earns its keep by being Reused,
and a copyleft term would Bar it from the BSD-licensed scientific stack it belongs beside.
