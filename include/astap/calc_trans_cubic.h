// Third order (cubic) coordinate transformation between two matched star lists,
// used to derive the SIP distortion coefficients.
// Ported from unit_command_line_calc_trans_cubic.pas, which itself is derived
// from the "match" package by Michael Richmond.

#pragma once

#include <string>
#include <vector>

namespace astap {
  struct SStar {
    double x = 0;
    double y = 0;
  };

  // x' = x00 + x10*x + x01*y + x20*x² + x11*x*y + x02*y² + x30*x³ + x21*x²y +
  //      x12*x*y² + x03*y³, and the same for y' with the y.. coefficients.
  struct TTrans {
    double x00 = 0, x10 = 0, x01 = 0, x20 = 0, x11 = 0, x02 = 0, x30 = 0, x21 = 0, x12 = 0, x03 = 0;
    double y00 = 0, y10 = 0, y01 = 0, y20 = 0, y11 = 0, y02 = 0, y30 = 0, y21 = 0, y12 = 0, y03 = 0;
  };

  // Given a set of matched pairs, find the transformation that takes the
  // coordinates of `stars_reference` into those of `stars_distorted`.
  // At least 10 pairs are required. Returns false with a message in `err_mess`
  // when the normal equations cannot be solved.
  bool calc_trans_cubic(const std::vector<SStar> &stars_reference,
                        const std::vector<SStar> &stars_distorted, TTrans &trans,
                        std::string &err_mess);
} // namespace astap
