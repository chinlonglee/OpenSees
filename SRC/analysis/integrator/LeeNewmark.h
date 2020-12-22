/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

#ifndef LeeNewmark_h
#define LeeNewmark_h

#include "Newmark.h"
#include <vector>
#include <sparseGEN/sp_sparse/triplet_form.hpp>

class LeeNewmark final : public Newmark {
	enum class MatType { None, InitialStiffness, CurrentStiffness, Stiffness, Mass, Damping };

	MatType which_matrix = MatType::None;

	std::vector<double> mass_coef, stiffness_coef;

	double omega = 0.;

	index_t n_damping = 0;
	index_t n_block = 0;

	double C6 = 0., C7 = 0.;

	bool first_iteration = true;

	const bool use_initial;

	Vector current_internal, trial_internal;
public:
	LeeNewmark(double, double, const std::vector<double>&, const std::vector<double>&, bool = false);

	int formTangent(int) override;

	int formEleTangent(FE_Element*) override;

	int formNodTangent(DOF_Group*) override;

	int formUnbalance() override;

	int commit() override;
	int revertToLastStep() override;
	int revertToStart() override;

	int domainChanged() override;

	int update(const Vector&) override;

	int newStep(double) override;
};

void* OPS_LeeNewmark(bool);

#endif
