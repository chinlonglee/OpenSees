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

// ReSharper disable CppClangTidyClangDiagnosticShorten64To32
// ReSharper disable CppClangTidyClangDiagnosticSignConversion
#include "LeeNewmarkFull.h"
#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <elementAPI.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <iostream>
#include <LinearSOE.h>
#include <sparseGEN/LeeSparse.h>

void* OPS_LeeNewmarkFull() {
	// TODO: need implementation

	return nullptr;
}

unsigned LeeNewmarkFull::get_amplifier() const {
	auto n_size = 2u;

	for(auto& I : damping_basis)
		if(Type::T0 == I.t) n_size += 5u;
		else if(Type::T1 == I.t) n_size += 5u + 6u * static_cast<unsigned>(I.p[0]);
		else if(Type::T2 == I.t) n_size += 4u + 5u * static_cast<unsigned>(.5 * (I.p[0] + I.p[1] - 1.));
		else if(Type::T3 == I.t) n_size += 9u;

	return n_size;
}

unsigned LeeNewmarkFull::get_total_size() const {
	auto n_size = 1u;

	for(auto& I : damping_basis)
		if(Type::T0 == I.t) n_size += 1u;
		else if(Type::T1 == I.t) n_size += 2u * static_cast<unsigned>(I.p[0]) + 1u;
		else if(Type::T2 == I.t) n_size += static_cast<unsigned>(I.p[0] + I.p[1]) + 1u;
		else if(Type::T3 == I.t) n_size += 2u;

	return n_size * n_block;
}

void LeeNewmarkFull::update_residual() const {
	// TODO: to be implemented
	
}

void LeeNewmarkFull::assemble_by_type_zero(triplet_form<double>&, unsigned&, double, double) const {}

void LeeNewmarkFull::assemble_by_type_one(triplet_form<double>& stiffness, unsigned& current_pos, const double mass_coef, const double stiffness_coef, double order) const {
	const auto mass_coefs = .5 * mass_coef;           // eq. 10
	const auto stiffness_coefs = .5 * stiffness_coef; // eq. 10

	// current masss and stiffness are now stored as objects so use dot (.) instead of ->
	const auto& m_row = current_mass.row_idx;
	const auto& m_col = current_mass.col_idx;
	const auto& m_val = current_mass.val_idx;
	const auto& s_row = current_stiffness.row_idx;
	const auto& s_col = current_stiffness.col_idx;
	const auto& s_val = current_stiffness.val_idx;

	index_t I = 0llu;
	index_t J = current_pos;
	index_t K = current_pos += n_block;
	index_t L = current_pos += n_block;
	index_t M = current_pos += n_block;

	while(order > 1.5) {
		// eq. 61
		for(index_t N = 0; N < current_mass.c_size; ++N) {
			const auto O = m_row[N], P = m_col[N], Q = O + J, R = P + J;
			stiffness.at(O + I, R) = stiffness.at(Q, P + I) = mass_coef * m_val[N];
			stiffness.at(Q, P + K) = stiffness.at(O + K, R) = stiffness.at(O + L, P + M) = stiffness.at(O + M, P + L) = mass_coefs * m_val[N];
		}

		for(index_t N = 0; N < current_stiffness.c_size; ++N) {
			const auto O = s_row[N], P = s_col[N], Q = O + K, R = O + L, S = P + K, T = P + L;
			stiffness.at(Q, T) = stiffness.at(R, S) = stiffness_coef * s_val[N];
			stiffness.at(O + J, S) = stiffness.at(Q, P + J) = stiffness.at(R, P + M) = stiffness.at(O + M, T) = stiffness_coefs * s_val[N];
		}

		I = current_pos;
		J = current_pos += n_block;
		K = current_pos += n_block;
		L = current_pos += n_block;
		M = current_pos += n_block;
		order -= 2.;
	}

	if(order < .5) assemble_by_type_zero(stiffness, current_pos = J, mass_coef, stiffness_coef);
	else {
		// eq. 53 
		for(index_t N = 0; N < current_mass.c_size; ++N) {
			const auto O = m_row[N], P = m_col[N], Q = O + J, R = P + J;
			stiffness.at(O + L, P + L) = -(stiffness.at(O + I, R) = stiffness.at(Q, P + I) = mass_coef * m_val[N]);
			stiffness.at(Q, P + K) = stiffness.at(O + K, R) = mass_coefs * m_val[N];
		}

		for(index_t N = 0; N < current_stiffness.c_size; ++N) {
			const auto O = s_row[N], P = s_col[N], Q = O + K, R = P + K, S = O + L, T = P + L;
			stiffness.at(S, T) = -(stiffness.at(Q, T) = stiffness.at(S, R) = stiffness_coef * s_val[N]);
			stiffness.at(O + J, R) = stiffness.at(Q, P + J) = stiffness_coefs * s_val[N];
		}
	}
}

void LeeNewmarkFull::assemble_by_type_two(triplet_form<double>&, unsigned&, double, double, double, double) const {}

void LeeNewmarkFull::assemble_by_type_three(triplet_form<double>&, unsigned&, double, double, double) const {}

void LeeNewmarkFull::assemble_by_type_four(triplet_form<double>&, unsigned&, double, double, double) const {}

LeeNewmarkFull::LeeNewmarkFull(const double _gamma, const double _beta, const std::vector<double>& X, const std::vector<double>& F, const bool _use_initial)
	: Newmark(_gamma, _beta) {}

int LeeNewmarkFull::formTangent(const int F) {
	statusFlag = F;

	auto result = 0;

	auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());
	auto t_model = getAnalysisModel();

	if(t_soe == nullptr || t_model == nullptr) {
		opserr << "WARNING LeeNewmarkFull::formTangent() no LeeSparse or AnalysisModel has been set\n";
		return -1;
	}

	if(first_iteration) {
		t_soe->zero_current_stiffness();
		t_soe->zero_mass();
	}
	t_soe->zero_stiffness();
	t_soe->zero_damping();

	auto& t_dof = t_model->getDOFs();
	DOF_Group* t_dofptr;

	while((t_dofptr = t_dof()) != nullptr) {
		if(first_iteration) {
			which_matrix = stiffness_type;
			if(t_soe->add_current_stiffness(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:dof\n";
				result = -1;
			}
			which_matrix = MatType::Mass;
			if(t_soe->add_mass(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_mass:dof\n";
				result = -1;
			}
		}
		which_matrix = MatType::TangentStiffness;
		if(t_soe->add_stiffness(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:dof\n";
			result = -1;
		}
		which_matrix = MatType::Damping;
		if(t_soe->add_damping(t_dofptr->getTangent(this), t_dofptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_damping:dof\n";
			result = -1;
		}
	}

	auto& t_ele = t_model->getFEs();
	FE_Element* t_eleptr;
	while((t_eleptr = t_ele()) != nullptr) {
		if(first_iteration) {
			which_matrix = stiffness_type;
			if(t_soe->add_current_stiffness(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:ele\n";
				result = -2;
			}
			which_matrix = MatType::Mass;
			if(t_soe->add_mass(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
				opserr << "LeeNewmarkFull::formTangent() - failed to add_mass:ele\n";
				result = -2;
			}
		}
		which_matrix = MatType::TangentStiffness;
		if(t_soe->add_stiffness(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_stiffness:ele\n";
			result = -2;
		}
		which_matrix = MatType::Damping;
		if(t_soe->add_damping(t_eleptr->getTangent(this), t_eleptr->getID()) < 0) {
			opserr << "LeeNewmarkFull::formTangent() - failed to add_damping:ele\n";
			result = -2;
		}
	}

	/* After this line, all matrices (of original size) are assembled. Elemental work is done.
	 * Now we deal with the formulation of the big matrix.
	 */

	auto& stiffness = t_soe->global_stiffness;

	stiffness.zeros();

	// global stiffness is now populated with global damping matrix
	// need to add K+M+C to top left corner

	if(first_iteration) {
		// preallocate memory
		stiffness.resize(get_amplifier() * t_soe->stiffness.c_size);
		stiffness.zeros();

		t_soe->mass.csc_condense();
		t_soe->current_stiffness.csc_condense();

		access::rw(current_mass) = t_soe->mass;
		access::rw(current_stiffness) = t_soe->current_stiffness;

		// now current mass and stiffness are formulated
		// assemble unrolled damping matrix and the corresponding damping force
		// assemble stiffness

		// populating global stiffness matrix with global damping matrix
		auto IDX = n_block;

		for(auto& I : damping_basis) {
			const auto mass_coef = 4. * I.zeta * I.omega * c2;
			const auto stiffness_coef = 4. * I.zeta / I.omega * c2;
			if(Type::T0 == I.t) assemble_by_type_zero(stiffness, IDX, mass_coef, stiffness_coef);
			else if(Type::T1 == I.t) assemble_by_type_one(stiffness, IDX, mass_coef, stiffness_coef, I.p[0]);
			else if(Type::T2 == I.t) assemble_by_type_two(stiffness, IDX, mass_coef, stiffness_coef, I.p[0], I.p[1]);
			else if(Type::T3 == I.t) assemble_by_type_three(stiffness, IDX, mass_coef, stiffness_coef, I.p[0]);
		}

		// now stiffness holds global damping matrix
		stiffness.csc_condense();

		// update residual according to global damping matrix
		update_residual();

		const auto& row = stiffness.row_idx;
		const auto& col = stiffness.col_idx;
		const auto& val = stiffness.val_idx;

		rabbit = triplet_form<double>(n_block, n_block, t_soe->stiffness.c_size);
		for(unsigned I = 0; I < stiffness.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// check in left top block of unrolled damping matrix to be used in subsequent iterations
			if(row[I] < n_block) rabbit.at(row[I], col[I]) = val[I];
		}

		// after this line, rabbit holds the top left corner of global damping matrix may or may not be zeros
		// depending on types

		// now add K, M, C to global effective stiffness matrix
		stiffness += t_soe->stiffness + t_soe->mass * c3 + t_soe->damping * c2;

		// no need to condense again
		// it will be done before solving automatically

		access::rw(first_iteration) = false;
	} else {
		// if not first iteration
		// erase the tangent stiffness entries
		if(!stiffness.csc_sort()) return -1;

		const auto& row = stiffness.row_idx;
		const auto& col = stiffness.col_idx;
		const auto& val = stiffness.val_idx;

		for(size_t I = 0; I < stiffness.c_size; ++I) {
			// quit if current column is beyond the original size of matrix
			if(col[I] >= n_block) break;
			// erase existing entries if fall in intact stiffness matrix
			if(row[I] < n_block) val[I] = 0.;
		}

		// check in original nonzero entries in unrolled damping matrix
		stiffness += rabbit;

		// ! now global damping matrix is recovered and stored in stiffness
		update_residual();

		// formulate global effective stiffness matrix
		stiffness += t_soe->stiffness + t_soe->mass * c3 + t_soe->damping * c2;
	}

	return result;
}

int LeeNewmarkFull::formEleTangent(FE_Element* theEle) {
	theEle->zeroTangent();

	if(MatType::InitialStiffness == which_matrix) theEle->addKiToTang();
	else if(MatType::TangentStiffness == which_matrix || MatType::CurrentStiffness == which_matrix)
		if(statusFlag == CURRENT_TANGENT) theEle->addKtToTang();
		else if(statusFlag == INITIAL_TANGENT) theEle->addKiToTang();
		else if(statusFlag == HALL_TANGENT) {
			theEle->addKtToTang(cFactor);
			theEle->addKiToTang(iFactor);
		} else { opserr << "Newmark::formEleTangent - unknown FLAG\n"; }
	else if(MatType::Mass == which_matrix) theEle->addMtoTang();
	else if(MatType::Damping == which_matrix) theEle->addCtoTang();

	return 0;
}

int LeeNewmarkFull::formNodTangent(DOF_Group* theDof) {
	theDof->zeroTangent();

	if(MatType::Mass == which_matrix) theDof->addMtoTang();
	else if(MatType::Damping == which_matrix) theDof->addCtoTang();

	return 0;
}

int LeeNewmarkFull::formUnbalance() {
	const auto flag = TransientIntegrator::formUnbalance();
	if(0 != flag) return flag;

	const auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());

	if(nullptr == t_soe) {
		opserr << "WARNING LeeNewmarkFull::formUnbalance() need a LeeSparse system\n";
		return -1;
	}

	auto& residual = t_soe->residual;

	if(0 == residual.Size()) residual.resize(static_cast<int>(get_total_size()));

	// TODO: formulate residual based on global effective stiffness matrix
	//
	// 
	// vec trial_vel = -trial_internal;
	// trial_vel.head(n_block) = factory->get_trial_velocity() / -C1;
	// access::rw(residual) = stiffness * trial_vel;
	// // ? may potentially improve performance
	// // access::rw(residual) = csr_form<double>(stiffness->triplet_mat) * trial_vel;
	//
	// // ! check in damping force
	// get_trial_damping_force(factory) -= residual.head(n_block);
	// get_incre_damping_force(factory) -= residual.head(n_block);
	// // ! update left hand side
	// get_sushi(factory) -= residual.head(n_block);

	return 0;
}

int LeeNewmarkFull::commit() {
	first_iteration = true;

	current_internal = trial_internal;

	return IncrementalIntegrator::commit();
}

int LeeNewmarkFull::revertToLastStep() {
	first_iteration = true;

	trial_internal = current_internal;

	return Newmark::revertToLastStep();
}

int LeeNewmarkFull::revertToStart() {
	first_iteration = true;

	trial_internal.Zero();
	current_internal.Zero();

	return Newmark::revertToStart();
}

int LeeNewmarkFull::domainChanged() {
	const auto flag = Newmark::domainChanged();

	if(flag != 0) return flag;

	n_block = getLinearSOE()->getX().Size();

	const auto n_size = get_total_size();

	first_iteration = true;

	current_internal = Vector(n_size);
	trial_internal = Vector(n_size);

	return 0;
}

int LeeNewmarkFull::update(const Vector& DU) {
	const auto t_soe = dynamic_cast<LeeSparse*>(getLinearSOE());

	trial_internal += t_soe->increment;

	return Newmark::update(DU);
}
