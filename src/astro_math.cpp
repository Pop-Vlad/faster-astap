#include "astap/astro_math.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace astap {
  double fnmodulo(double x, double range) {
    double result = x - range * std::floor(x / range);
    if (result >= range) result = result - range;
    if (result < 0) result = result + range;
    return result;
  }

  double smedian(std::vector<double> &list, size_t leng) {
    if (leng == 0) return std::nan("");
    if (leng == 1) return list[0];

    std::sort(list.begin(), list.begin() + static_cast<long>(leng));
    size_t mid = (leng - 1) / 2;
    if (leng % 2 == 1) {
      if (leng <= 3) return list[mid];
      // The original averages three values around the centre to smooth the result.
      return (list[mid - 1] + list[mid] + list[mid + 1]) / 3.0;
    }
    return (list[mid] + list[mid + 1]) / 2.0;
  }

  void ang_sep(double ra1, double dec1, double ra2, double dec2, double &sep) {
    double sin_dec1 = std::sin(dec1), cos_dec1 = std::cos(dec1);
    double sin_dec2 = std::sin(dec2), cos_dec2 = std::cos(dec2);
    // Clamp to prevent a domain error on values such as 1.000000000002.
    double cos_sep = std::max(-1.0, std::min(1.0, sin_dec1 * sin_dec2 +
                                                  cos_dec1 * cos_dec2 * std::cos(ra1 - ra2)));
    sep = std::acos(cos_sep);
  }

  void rotate(double rot, double x, double y, double &x2, double &y2) {
    double sin_rot = std::sin(rot), cos_rot = std::cos(rot);
    x2 = x * +cos_rot + y * sin_rot;
    y2 = x * -sin_rot + y * cos_rot;
  }

  double position_angle(double ra1, double dec1, double ra0, double dec0) {
    double sinDeltaRa = std::sin(ra1 - ra0), cosDeltaRa = std::cos(ra1 - ra0);
    double sinDec0 = std::sin(dec0), cosDec0 = std::cos(dec0);
    double sinDec1 = std::sin(dec1), cosDec1 = std::cos(dec1);
    return std::atan2(cosDec1 * sinDeltaRa, sinDec1 * cosDec0 - cosDec1 * sinDec0 * cosDeltaRa);
  }

  void equatorial_standard(double ra0, double dec0, double ra, double dec, double cdelt,
                           double &xx, double &yy) {
    double sin_dec0 = std::sin(dec0), cos_dec0 = std::cos(dec0);
    double sin_dec = std::sin(dec), cos_dec = std::cos(dec);
    double sin_deltaRA = std::sin(ra - ra0), cos_deltaRA = std::cos(ra - ra0);
    // cdelt/(3600*180/pi) converts standard coordinates into CCD pixels.
    double dv = (cos_dec0 * cos_dec * cos_deltaRA + sin_dec0 * sin_dec) * cdelt / (3600 * 180 / kPi);
    xx = -cos_dec * sin_deltaRA / dv; // tangent of the angle in RA
    yy = -(sin_dec0 * cos_dec * cos_deltaRA - cos_dec0 * sin_dec) / dv; // tangent of the angle in DEC
  }

  void standard_equatorial(double ra0, double dec0, double x, double y, double cdelt,
                           double &ra, double &dec) {
    double sin_dec0 = std::sin(dec0), cos_dec0 = std::cos(dec0);
    x = x * cdelt / (3600 * 180 / kPi); // scale CCD pixels to standard coordinates
    y = y * cdelt / (3600 * 180 / kPi);

    double delta = cos_dec0 - y * sin_dec0;
    ra = ra0 + std::atan2(-x, delta); // atan2 is required for images containing the celestial pole
    dec = std::atan((sin_dec0 + y * cos_dec0) / std::sqrt(x * x + delta * delta));
    if (ra > kPi * 2) ra = ra - kPi * 2; // prevent values above 2*pi
    if (ra < 0) ra = ra + kPi * 2;
  }

  double apply_arctan(double fov) { return 2 * std::atan((fov / 2) * kPi / 180) * 180 / kPi; }

  std::string float_to_str(double x, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, x);
    return std::string(buf);
  }

  std::string distance_to_string(double dist, double inp) {
    if (std::fabs(dist) < kPi / (180 * 60)) // unit seconds
      return float_to_str(inp * 3600 * 180 / kPi, 1) + "\"";
    if (std::fabs(dist) < kPi / 180) // unit minutes
      return float_to_str(inp * 60 * 180 / kPi, 1) + "'";
    return float_to_str(inp * 180 / kPi, 1) + "d";
  }

  static std::string leading_zero(int v) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%02d", v);
    return std::string(buf);
  }

  std::string prepare_ra(double rax, const std::string &sep) {
    // Add 1/10 of half a second to round correctly and avoid 7:60 results.
    rax = rax + kPi * 0.1 / (24 * 60 * 60);
    rax = rax * 12 / kPi; // make hours
    int h = static_cast<int>(ptrunc(rax));
    int m = static_cast<int>(ptrunc((rax - h) * 60));
    int s = static_cast<int>(ptrunc((rax - h - m / 60.0) * 3600));
    int ds = static_cast<int>(ptrunc((rax - h - m / 60.0 - s / 3600.0) * 36000));
    return leading_zero(h) + sep + leading_zero(m) + "  " + leading_zero(s) + "." +
           static_cast<char>(ds + '0');
  }

  std::string prepare_dec(double decx, const std::string &sep) {
    char sign = decx < 0 ? '-' : '+';
    // Add half a second to round correctly and avoid 7:60 results.
    decx = std::fabs(decx) + kPi / (360 * 60 * 60);
    decx = decx * 180 / kPi; // make degrees
    int g = static_cast<int>(ptrunc(decx));
    int m = static_cast<int>(ptrunc((decx - g) * 60));
    int s = static_cast<int>(ptrunc((decx - g - m / 60.0) * 3600));
    return std::string(1, sign) + leading_zero(g) + sep + leading_zero(m) + "  " + leading_zero(s);
  }

  // Shared helper for the two sexagesimal parsers. `separators` lists the unit
  // letters that have to be treated as white space.
  static bool parse_sexagesimal(std::string inp, const std::string &separators, double &v1,
                                double &v2, double &v3, int &plusmin) {
    for (char &c: inp) {
      if (c == ',') c = '.';
      if (c == ':') c = ' ';
      if (separators.find(c) != std::string::npos) c = ' ';
    }
    plusmin = inp.find('-') != std::string::npos ? -1 : 1;

    std::istringstream is(inp);
    std::vector<double> parts;
    std::string token;
    while (is >> token) {
      try {
        size_t used = 0;
        double d = std::stod(token, &used);
        if (used != token.size()) return false;
        parts.push_back(d);
      } catch (...) {
        return false;
      }
    }
    if (parts.empty()) return false;
    v1 = parts[0];
    v2 = parts.size() > 1 ? parts[1] : 0.0;
    v3 = parts.size() > 2 ? parts[2] : 0.0;
    return true;
  }

  bool ra_text_to_radians(const std::string &inp, double &ra) {
    double rah, ram, ras;
    int plusmin;
    if (!parse_sexagesimal(inp, "hms", rah, ram, ras, plusmin)) return false;
    ra = plusmin * (std::fabs(rah) + ram / 60 + ras / 3600) * kPi / 12;
    return ra <= 2 * kPi;
  }

  bool dec_text_to_radians(const std::string &inp, double &dec) {
    double decd, decm, decs;
    int plusmin;
    if (!parse_sexagesimal(inp, "dms\xb0", decd, decm, decs, plusmin)) return false;
    dec = plusmin * (std::fabs(decd) + decm / 60 + decs / 3600) * kPi / 180;
    return true;
  }

  void precession_jnow_to_j2000(double equinox, double &ra, double &dec) {
    double t = (equinox - 2000) / 100; // time in julian centuries since J2000
    double m = 3.07496 + 0.00186 * t; // seconds
    double n = 1.33621 - 0.00057 * t; // seconds
    double n2 = 20.0431 - 0.0085 * t; // arcsec
    double dra = (m + n * std::sin(ra) * std::tan(dec)) * kPi / (3600 * 12);
    double ddec = n2 * std::cos(ra) * kPi / (3600 * 180);
    ra = ra - (dra * t * 100); // t*100 is the number of years, subtract to go back to J2000
    dec = dec - (ddec * t * 100);
  }
} // namespace astap
