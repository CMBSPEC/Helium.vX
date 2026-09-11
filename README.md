# Helium.vX

This repository contains the neutral helium atom module used in the CMBSPEC
atomic-model setup. The code was developed for detailed cosmological
recombination calculations, where helium bound-bound, bound-free, and
recombination processes have to be treated consistently.

The repository is intended to hold the Helium atom class and the data tables
that belong directly to that class, without copying the more general Tools or
Development libraries from the larger codes.

The main interface is provided by `HeI_Atom.h` / `HeI_Atom.cpp`. The module
sets up singlet and triplet helium levels, transition data, photoionization
cross sections, and recombination-related rates. Runtime data are stored in
`Helium.Data/`, including transition tables and TOPbase photoionization data.

## Data Path

By default, the code looks for data through the `HEIDATADIR` build definition
and then `./Helium.Data/`. Recent versions also allow the data path to be set
explicitly through the helium atom setup, which is useful when the module is
used outside its original CosmoSpec/CosmoRec directory layout.

## Dependencies

This repository is not a standalone application. It depends on common helper
code from the surrounding CMBSPEC/CosmoSpec/CosmoTherm toolchains, including
physical constants, file I/O helpers, interpolation/integration routines, Voigt
profiles, hydrogenic bound-bound and photoionization routines.

A small external demonstration project,
[Helium-demo](https://github.com/CMBSPEC/Helium-demo), will show the minimal
set of required Tools and a concrete build/use example.

## Related Literature And Data Sources

The code and tables are connected to several atomic-physics and recombination
references, including:

- G. W. F. Drake and D. C. Morton, "A Multiplet Table for Neutral Helium
  (^4He I) with Transition Rates", ApJS 170, 251, 2007.
- R. A. Benjamin, E. D. Skillman, and D. P. Smits, "Improving Predictions for
  Helium Emission Lines", ApJ 514, 307, 1999.
- G. Lach and K. Pachucki, "Forbidden transitions in the helium atom",
  Phys. Rev. A 64, 042510, 2001.
- W. Cunto, C. Mendoza, F. Ochsenbein, and C. J. Zeippen, "TOPbase at the CDS",
  A&A 275, L5, 1993.
- J. A. Rubino-Martin, J. Chluba, and R. A. Sunyaev, "Lines in the cosmic
  microwave background spectrum from the epoch of cosmological helium
  recombination", A&A 485, 377, 2008.
- J. Chluba and R. M. Thomas, "Towards a complete treatment of the
  cosmological recombination problem", MNRAS 412, 748, 2011.

These references are listed to document the origin and scientific context of
the data and methods used here; the repository itself only provides the helium
atom setup.

**Acknowledgements:** This repository was made available and documented with the help of Codex. The related release work was supported in part by a grant of access to OpenAI models through the ChatGPT for Academic Researchers program.
