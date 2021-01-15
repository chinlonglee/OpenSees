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

#ifndef LeeNewmarkFull_h
#define LeeNewmarkFull_h

#include "Newmark.h"
#include <vector>
#include <sparseGEN/sp_sparse/triplet_form.hpp>

class LeeNewmarkFull final : public Newmark {
public:
	enum class Type { T0, T1, T2, T3 };

	struct Basis {
		Type t;
		std::vector<double> p;
		double zeta, omega;
	};

private:
	std::vector<Basis> damping_basis;

	triplet_form<double> rabbit;
	const triplet_form<double> current_stiffness;
	const triplet_form<double> current_mass;

	// ---------------------------------------------------------------------------------------------------

	unsigned get_amplifier() const;
	unsigned get_total_size() const;
	void update_residual() const;

	void assemble_by_type_zero(triplet_form<double>&, unsigned&, double, double) const;
	void assemble_by_type_one(triplet_form<double>&, unsigned&, double, double, double) const;
	void assemble_by_type_two(triplet_form<double>&, unsigned&, double, double, double, double) const;
	void assemble_by_type_three(triplet_form<double>&, unsigned&, double, double, double) const;
	void assemble_by_type_four(triplet_form<double>&, unsigned&, double, double, double) const;

	enum class MatType : unsigned { None, InitialStiffness, CurrentStiffness, TangentStiffness, Mass, Damping };

	MatType which_matrix = MatType::None;

	const MatType stiffness_type = MatType::CurrentStiffness;

	unsigned n_block = 0;

	bool first_iteration = true;

	Vector current_internal, trial_internal;
public:
	LeeNewmarkFull(double,               // alpha
	               double,               // beta
	               std::vector<Basis>&&, // damping basis vector
	               unsigned              // flag to indicate which stiffness to be used
	);

	int formTangent(int) override;

	int formEleTangent(FE_Element*) override;

	int formNodTangent(DOF_Group*) override;

	int commit() override;
	int revertToLastStep() override;
	int revertToStart() override;

	int domainChanged() override;

	int update(const Vector&) override;
};

void* OPS_LeeNewmarkFull(unsigned);

#endif
