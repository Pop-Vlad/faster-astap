// Step 6 to 8 of the ASTAP method: find quad hash code matches, remove the
// outliers on the longest side ratio and solve the overdetermined system of
// linear equations for the six plate constants.
//
// Ported from unit_command_line_solving.pas (lsq_fit, find_fit,
// find_fit_using_hash, find_offset_and_rotation).

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "astap/types.h"

namespace astap {
  // The global matching state of the Pascal unit, bundled so that it can be
  // passed around explicitly.
  struct MatchState {
    RowList quad_star_distances1; // database / reference quads
    RowList quad_star_distances2; // image quads

    RowList a_xy_positions; // A matrix, quad centres of the image
    std::vector<double> b_xrefpositions; // X_ref, quad centres of the database
    std::vector<double> b_yrefpositions; // Y_ref

    int nr_references = 0; // matches left after the outlier filter
    int nr_references2 = 0; // matches found on the hash codes

    SolutionVector solution_vector_x; // X_ref := a1*x + b1*y + c1
    SolutionVector solution_vector_y; // Y_ref := a2*x + b2*y + c2
  };

  // Find the solution vector of an overdetermined system of linear equations
  // with the method of least squares, using Givens rotations.
  // See Montenbruck & Pfleger, Astronomy on the personal computer.
  bool lsq_fit(const RowList &a_matrix, std::vector<double> b_matrix, SolutionVector &x_matrix);

  // Brute force hash code matching, used for fewer than 120 database quads.
  bool find_fit(MatchState &st, int minimum_count, double quad_tolerance,
                const std::function<void(const std::string &)> &log = nullptr);

  // Hash based matching, used from 120 database quads onwards.
  bool find_fit_using_hash(MatchState &st, int minimum_count, double quad_tolerance,
                           const std::function<void(const std::string &)> &log = nullptr);

  // Resets the solution vectors to a scaled identity transformation.
  void reset_solution_vectors(MatchState &st, double factor);

  // Find the difference between the reference image / database and the new image.
  // Returns true when a solution was found.
  bool find_offset_and_rotation(MatchState &st, int minimum_quads, double tolerance,
                                const std::function<void(const std::string &)> &log = nullptr);
} // namespace astap
