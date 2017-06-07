/// \file
/// Contains constants for unit conversions.
///
/// The units for this code are:
///     - Time in femtoseconds (fs)
///     - Length in Angstroms (Angs)
///     - Energy in electron Volts (eV)
///     - Mass read in as Atomic Mass Units (amu) and then converted for
///       consistency (energy*time^2/length^2)
/// Values are taken from NIST, http://physics.nist.gov/cuu/Constants/

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

/// 1 amu in kilograms
static constexpr const double  amuInKilograms = 1.660538921e-27;

/// 1 fs in seconds
static constexpr const double  fsInSeconds  =   1.0e-15;

/// 1 Ang in meters
static constexpr const double  AngsInMeters  =  1.0e-10;

/// 1 eV in Joules
static constexpr const double  eVInJoules   =   1.602176565e-19;

/// Internal mass units are eV * fs^2 / Ang^2
static constexpr const double amuToInternalMass =
         amuInKilograms * AngsInMeters * AngsInMeters
         / (fsInSeconds * fsInSeconds  * eVInJoules);

/// Boltmann constant in eV's
static constexpr const double kB_eV = 8.6173324e-5;  // eV/K

/// Hartrees to eVs
static constexpr const double hartreeToEv = 27.21138505;

/// Bohrs to Angstroms
static constexpr const double bohrToAngs = 0.52917721092;

#endif
