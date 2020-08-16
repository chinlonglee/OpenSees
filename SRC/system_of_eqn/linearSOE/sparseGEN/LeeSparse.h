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

#ifndef LeeSparse_h
#define LeeSparse_h

#include "SparseGenColLinSOE.h"
#include "sp_sparse/triplet_form.hpp"

class LeeSparse final : public SparseGenColLinSOE {
public:
	using SparseGenColLinSOE::SparseGenColLinSOE;

	int setSize(Graph&) override;

	int solve() override;

	int add_current_stiffness(const Matrix&, const ID&, double = 1.);
	int add_stiffness(const Matrix&, const ID&, double = 1.);
	int add_damping(const Matrix&, const ID&, double = 1.);
	int add_mass(const Matrix&, const ID&, double = 1.);

	void zero_current_stiffness() const;
	void zero_stiffness() const;
	void zero_damping() const;
	void zero_mass() const;

	triplet_form<double> global_stiffness, stiffness, damping, mass, current_stiffness;

	Vector increment, residual;
};

#endif
