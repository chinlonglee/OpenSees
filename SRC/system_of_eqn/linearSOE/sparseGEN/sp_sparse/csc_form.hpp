/*******************************************************************************
 * Copyright (C) 2017-2020 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#ifndef CSC_FORM_HPP
#define CSC_FORM_HPP

#include "sparse_form.hpp"
#include <cstring>

template<typename T> class triplet_form;
template<typename T> class csr_form;

template<typename T> class csc_form final : public sparse_form<T, csc_form<T>> {
	using sparse_form<T, csc_form<T>>::bin;

	void copy_memory(index_t, const index_t*, const index_t*, const T*) override;
public:
	using sparse_form<T, csc_form<T>>::n_rows;
	using sparse_form<T, csc_form<T>>::n_cols;
	using sparse_form<T, csc_form<T>>::n_elem;
	using sparse_form<T, csc_form<T>>::c_size;

	index_t* row_idx = nullptr; // index storage
	index_t* col_ptr = nullptr; // index storage
	T* val_idx = nullptr;       // value storage

	csc_form() = default;

	~csc_form() override;

	csc_form(const csc_form&);                // copy ctor
	csc_form(csc_form&&) noexcept;            // move ctor
	csc_form& operator=(const csc_form&);     // copy assignment
	csc_form& operator=(csc_form&&) noexcept; // move assignment

	void reset() const override;
	void zeros() const override;

	T max() const override;

	bool init() override;
	bool init(index_t) override;
	bool init(index_t, index_t, index_t) override;
	bool resize() override;
	bool resize(index_t) override;
	bool resize(index_t, index_t, index_t) override;

	csc_form<T> transpose() const;

	template<typename T2> csc_form<T> operator*(T2);
	template<typename T2> csc_form<T> operator/(T2);
	template<typename T2> csc_form<T>& operator*=(T2);
	template<typename T2> csc_form<T>& operator/=(T2);

	explicit csc_form(triplet_form<T>&);
	explicit csc_form(triplet_form<T>&&);
	csc_form& operator=(const triplet_form<T>&);

	explicit csc_form(const csr_form<T>&);
	csc_form& operator=(const csr_form<T>&);

	const T& operator()(index_t, index_t) const;

	Vector operator*(const Vector&) override;
};

template<typename T> void csc_form<T>::copy_memory(const index_t size, const index_t* const in_row_idx, const index_t* const in_col_ptr, const T* const in_val_idx) {
	if(size > n_elem) resize(size);

	auto bytes = size * sizeof(index_t);
	memcpy(this->row_idx, in_row_idx, bytes);
	bytes = (n_cols + 1l) * sizeof(index_t);
	memcpy(this->col_ptr, in_col_ptr, bytes);
	bytes = size * sizeof(T);
	memcpy(this->val_idx, in_val_idx, bytes);

	access::rw(c_size) = size;
}

template<typename T> csc_form<T>::~csc_form() { csc_form<T>::reset(); }

template<typename T> csc_form<T>::csc_form(const csc_form& in_mat)
	: sparse_form<T, csc_form<T>>(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem) {
	csc_form<T>::init();
	csc_form<T>::copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_ptr, in_mat.val_idx);
}

template<typename T> csc_form<T>::csc_form(csc_form&& in_mat) noexcept {
	csc_form<T>::reset();
	n_rows = in_mat.n_rows;
	n_cols = in_mat.n_cols;
	n_elem = in_mat.n_elem;
	c_size = in_mat.c_size;
	row_idx = in_mat.row_idx;
	col_ptr = in_mat.col_ptr;
	val_idx = in_mat.val_idx;
	in_mat.n_rows = in_mat.n_cols = in_mat.n_elem = in_mat.c_size = 0;
	in_mat.row_idx = in_mat.col_ptr = nullptr;
	in_mat.val_idx = nullptr;
}

template<typename T> csc_form<T>& csc_form<T>::operator=(const csc_form& in_mat) {
	if(this != &in_mat) {
		init(in_mat.n_rows, in_mat.n_cols, in_mat.n_elem);
		copy_memory(in_mat.c_size, in_mat.row_idx, in_mat.col_ptr, in_mat.val_idx);
	}

	return *this;
}

template<typename T> csc_form<T>& csc_form<T>::operator=(csc_form&& in_mat) noexcept {
	reset();
	access::rw(n_rows) = in_mat.n_rows;
	access::rw(n_cols) = in_mat.n_cols;
	access::rw(n_elem) = in_mat.n_elem;
	access::rw(c_size) = in_mat.c_size;
	row_idx = in_mat.row_idx;
	col_ptr = in_mat.col_ptr;
	val_idx = in_mat.val_idx;
	access::rw(in_mat.n_rows) = access::rw(in_mat.n_cols) = access::rw(in_mat.n_elem) = access::rw(in_mat.c_size) = 0;
	in_mat.row_idx = in_mat.col_ptr = nullptr;
	in_mat.val_idx = nullptr;
	return *this;
}

template<typename T> void csc_form<T>::reset() const {
	zeros();
	delete[] row_idx;
	delete[] col_ptr;
	delete[] val_idx;
}

template<typename T> void csc_form<T>::zeros() const { access::rw(c_size) = 0; }

template<typename T> T csc_form<T>::max() const { return *std::max_element(val_idx, val_idx + c_size); }

template<typename T> bool csc_form<T>::init() {
	reset();
	row_idx = new(std::nothrow) index_t[n_elem];
	col_ptr = new(std::nothrow) index_t[n_cols + 1];
	val_idx = new(std::nothrow) T[n_elem];

	if(row_idx == nullptr || col_ptr == nullptr || val_idx == nullptr) {
		reset();
		return false;
	}
	return true;
}

template<typename T> bool csc_form<T>::init(const index_t in_elem) {
	if(in_elem <= n_elem) {
		zeros();
		return true;
	}
	access::rw(n_elem) = in_elem;
	return init();
}

template<typename T> bool csc_form<T>::init(const index_t in_row, const index_t in_col, const index_t in_elem) {
	if(n_rows != in_row) access::rw(n_rows) = in_row;
	if(n_cols != in_col) access::rw(n_cols) = in_col;

	return init(in_elem);
}

template<typename T> bool csc_form<T>::resize() {
	const auto copy = *this;

	if(!init(n_elem == 0 ? 1 : 2 * n_elem)) return false;

	copy_memory(copy.c_size, copy.row_idx, copy.col_ptr, copy.val_idx);

	return true;
}

template<typename T> bool csc_form<T>::resize(const index_t in_elem) {
	const auto copy = *this;

	if(in_elem <= c_size || !init(in_elem)) return false;

	copy_memory(copy.c_size, copy.row_idx, copy.col_ptr, copy.val_idx);

	return true;
}

template<typename T> bool csc_form<T>::resize(const index_t in_row, const index_t in_col, const index_t in_elem) {
	const auto copy = *this;

	if(in_row < n_rows || in_col < n_cols || in_elem < c_size || !init(in_row, in_col, in_elem)) return false;

	copy_memory(copy.c_size, copy.row_idx, copy.col_ptr, copy.val_idx);

	return true;
}

template<typename T> csc_form<T> csc_form<T>::transpose() const {
	csr_form<T> copy(*this);
	csc_form<T> out;

	out.n_cols = n_rows;
	out.n_rows = n_cols;
	out.n_elem = n_elem;
	out.c_size = c_size;

	out.row_idx = copy.col_idx;
	out.col_ptr = copy.row_ptr;
	out.val_idx = copy.val_idx;

	copy.row_ptr = copy.col_idx = nullptr;
	copy.val_idx = nullptr;

	return out;
}

template<typename T> template<typename T2> csc_form<T> csc_form<T>::operator*(const T2 scalar) {
	csc_form<T> copy = *this;

#ifdef SUANPAN_MT
	tbb::parallel_for(index_t(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] *= T(scalar); });
#else
	for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] *= T(scalar);
#endif

	return copy;
}

template<typename T> template<typename T2> csc_form<T> csc_form<T>::operator/(const T2 scalar) {
	csc_form<T> copy = *this;

#ifdef SUANPAN_MT
	tbb::parallel_for(index_t(0), copy.c_size, [&](const index_t I) { copy.val_idx[I] /= T(scalar); });
#else
	for(auto I = 0; I < copy.c_size; ++I) copy.val_idx[I] /= T(scalar);
#endif

	return copy;
}

template<typename T> template<typename T2> csc_form<T>& csc_form<T>::operator*=(const T2 scalar) {
#ifdef SUANPAN_MT
	tbb::parallel_for(index_t(0), c_size, [&](const index_t I) { val_idx[I] *= T(scalar); });
#else
	for(auto I = 0; I < c_size; ++I) val_idx[I] *= T(scalar);
#endif

	return *this;
}

template<typename T> template<typename T2> csc_form<T>& csc_form<T>::operator/=(const T2 scalar) {
#ifdef SUANPAN_MT
	tbb::parallel_for(index_t(0), c_size, [&](const index_t I) { val_idx[I] /= T(scalar); });
#else
	for(auto I = 0; I < c_size; ++I) val_idx[I] /= T(scalar);
#endif

	return *this;
}

template<typename T> csc_form<T>::csc_form(triplet_form<T>& old_mat) { *this = old_mat; }

template<typename T> csc_form<T>::csc_form(triplet_form<T>&& old_mat) { *this = old_mat; }

template<typename T> csc_form<T>& csc_form<T>::operator=(const triplet_form<T>& in_mat) {
	in_mat.csc_condense();

	init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

	if(in_mat.c_size == 0) return *this;

	access::rw(c_size) = in_mat.c_size;

	auto bytes = in_mat.c_size * sizeof(index_t);
	memcpy(this->row_idx, in_mat.row_idx, bytes);
	bytes = in_mat.c_size * sizeof(T);
	memcpy(this->val_idx, in_mat.val_idx, bytes);

	index_t current_pos = 0, current_col = 0;
	while(current_pos < in_mat.c_size)
		if(in_mat.col_idx[current_pos] < current_col) ++current_pos;
		else col_ptr[current_col++] = current_pos;

	col_ptr[n_cols] = c_size;

	return *this;
}

template<typename T> csc_form<T>::csc_form(const csr_form<T>& in_mat) { *this = in_mat; }

template<typename T> csc_form<T>& csc_form<T>::operator=(const csr_form<T>& in_mat) {
	csc_form<T>::init(in_mat.n_rows, in_mat.n_cols, in_mat.c_size);

#ifdef SUANPAN_MT
	tbb::parallel_for(index_t(0), n_cols + 1, [&](const index_t I) { col_ptr[I] = 0; });
#else
	for(index_t I = 0; I <= n_cols; ++I) col_ptr[I] = 0;
#endif

	for(index_t I = 0; I < in_mat.c_size; ++I) ++col_ptr[in_mat.col_idx[I] + 1];
	for(index_t I = 2; I <= n_cols; ++I) col_ptr[I] += col_ptr[I - 1];

	std::vector<index_t> counter(n_cols, 0);

	index_t r_idx = 1;
	for(index_t I = 0; I < in_mat.c_size; ++I) {
		if(I >= in_mat.row_ptr[r_idx]) ++r_idx;
		const auto& c_idx = in_mat.col_idx[I];
		const auto c_pos = counter[c_idx]++ + col_ptr[c_idx];
		row_idx[c_pos] = r_idx - 1;
		val_idx[c_pos] = in_mat.val_idx[I];
	}

	c_size = in_mat.c_size;

	return *this;
}

template<typename T> const T& csc_form<T>::operator()(const index_t in_row, const index_t in_col) const {
	if(in_row < n_rows && in_col < n_cols) for(auto I = col_ptr[in_col]; I < col_ptr[in_col + 1]; ++I) if(row_idx[I] == in_row) return val_idx[I];

	access::rw(bin) = 0.;
	return bin;
}

template<typename T> Vector csc_form<T>::operator*(const Vector& in_mat) {
	Vector out_mat(in_mat.Size());

	for(index_t I = 0; I < n_cols; ++I) for(auto J = col_ptr[I]; J < col_ptr[I + 1]; ++J) out_mat(row_idx[J]) += val_idx[J] * in_mat(I);

	return out_mat;
}

#endif
