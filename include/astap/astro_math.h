// Small astronomical / numerical helpers used by the solver.
// Ported from unit_command_line_solving.pas and unit_command_line_general.pas.

#pragma once

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace astap {
  constexpr double kPi = 3.14159265358979323846;

  // FPC's Round() rounds half to even (the default FPU rounding mode). std::round
  // rounds half away from zero, which would drift from the original on ties, so
  // the port uses nearbyint() with the default rounding mode instead.
  inline long pround(double x) { return static_cast<long>(std::nearbyint(x)); }

  // Pascal trunc(): rounds towards zero.
  inline long ptrunc(double x) { return static_cast<long>(std::trunc(x)); }

  // Pascal frac(): fractional part, keeps the sign of the argument.
  inline double pfrac(double x) { return x - std::trunc(x); }

  inline double sqr(double x) { return x * x; }

  // Pascal: fnmodulo. Always returns a value in [0, range).
  double fnmodulo(double x, double range);

  // Median of the first `leng` values. The list is sorted in place, exactly like
  // the Pascal SMedian, so callers that still need the original order must pass a
  // copy. The pointer overload lets hot callers use a stack buffer.
  double smedian(double *list, size_t leng);

  inline double smedian(std::vector<double> &list, size_t leng) {
    return smedian(list.data(), leng);
  }

  // Angular separation, formula 9.1 (old Meeus) / 16.1 (new Meeus).
  void ang_sep(double ra1, double dec1, double ra2, double dec2, double &sep);

  // Rotate a vector point, angle seen from the y-axis, counter clockwise.
  void rotate(double rot, double x, double y, double &x2, double &y2);

  // Position angle of a body at ra1,dec1 as seen at ra0,dec0. Rigorous method,
  // see Meeus, Astronomical Algorithms, formula 48.5 (1998 edition).
  double position_angle(double ra1, double dec1, double ra0, double dec0);

  // Transformation of equatorial coordinates into CCD pixel coordinates for
  // optical (gnomonic) projection, rigid method.
  //   ra0, dec0 : right ascension and declination of the optical axis [rad]
  //   ra, dec   : star position [rad]
  //   cdelt     : CCD scale in arcsec per pixel
  //   xx, yy    : resulting CCD coordinates
  void equatorial_standard(double ra0, double dec0, double ra, double dec,
                           double cdelt, double &xx, double &yy);

  // Inverse of equatorial_standard. The arctan form is used because it stays
  // accurate close to the poles.
  void standard_equatorial(double ra0, double dec0, double x, double y,
                           double cdelt, double &ra, double &dec);

  // Assume the optical system can be modelled by a simple arctan function like a
  // standard rectilinear (pinhole) lens.
  double apply_arctan(double fov);

  // Angular distance to string, the unit is chosen from `dist`.
  std::string distance_to_string(double dist, double inp);

  // Radians to text, format "hh: mm  ss.s".
  std::string prepare_ra(double rax, const std::string &sep);

  // Radians to text, format "+dd° mm  ss".
  std::string prepare_dec(double decx, const std::string &sep);

  // Sexagesimal (or decimal) text to radians. Accepts ':', 'd', 'h', 'm', 's' and
  // ',' as separators. Returns false when the text could not be parsed.
  bool ra_text_to_radians(const std::string &inp, double &ra);

  bool dec_text_to_radians(const std::string &inp, double &dec);

  // Simple precession correction, new Meeus chapter precession formula 20.1.
  void precession_jnow_to_j2000(double equinox, double &ra, double &dec);

  std::string float_to_str(double x, int decimals);
} // namespace astap
