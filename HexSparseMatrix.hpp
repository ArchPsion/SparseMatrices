#ifndef __HEX_SPARSE_MATRIX_HPP__
#define __HEX_SPARSE_MATRIX_HPP__

// Standard Libraries
#include <vector>

// Qt Libraries
#include <QtGlobal>
#include <QtMath>

// Custom Libraries
#include "HexRandomGenerator.hpp"

struct HexColumnValuePair
{
	qreal	value = 0.;
	qint32	column = 0;
	
	HexColumnValuePair(qreal v = 0., qint32 c = 0) : value(v), column(c)
	{
	}
};

/* It's more efficient to tie these two values together
 * and make one vector, rather than use two separate vectors
 * that would always have the same length, and that would
 * always be parsed at the same time one after the other.
 */

struct HexDecomposition;

class HexSparseMatrix
{
	private:
	
		template<typename Type> inline static void				Change(std::vector<HexColumnValuePair>&, Type, Type, qreal);
		inline static qreal							Normalise(std::vector<HexColumnValuePair>&);
		template<typename Type1, typename Type2> inline static void		Rewrite(Type1&, Type2, Type2);
		template<typename Type1, typename Type2> inline static qreal		Scalar(Type1, Type1, Type2, Type2);
	
		std::vector<HexColumnValuePair>						pairs;
		std::vector<qint32>							rowOffsets;
		
		qint32									numberOfRows = 0;
		qint32									numberOfColumns = 0;
		
		inline void								addValue(qint32, qint32, qreal);
		inline qint32								indexOfFreeCell(qint32, qint32) const;
		inline bool								insertColumnValuePair(qint32, qint32, qint32, qreal);
		inline void								removeValue(qint32, qint32);
		inline void								updateNumberOfColumns(void);
		inline void								updateNumberOfRows(void);
	
	public:
	
		inline									HexSparseMatrix(void);
		inline									HexSparseMatrix(const std::vector<std::vector<HexColumnValuePair>>&, qint32);
	
		inline void								downsize(void);
		inline QString								getData(void) const;
		inline HexDecomposition							getDecomposition(void) const;
		inline std::vector<qreal>						getDenseMatrix(void) const;
		inline qreal								getDensity(void) const;
		inline qint32								getNumberOfColumns(void) const;
		inline qint32								getNumberOfRows(void) const;
		inline const std::vector<HexColumnValuePair>&				getPairs(void) const;
		inline qint32								getRank(void) const;
		inline const std::vector<qint32>&					getRowOffsets(void) const;
		inline qreal								getSparsity(void) const;
		inline bool								insertOne(HexRandomGenerator&);
		inline void								setValue(qint32, qint32, qreal);
		inline bool								shuffle(HexRandomGenerator&);
		inline void								swapColumns(qint32, qint32);
		inline void								swapRows(qint32, qint32);
		inline void								transpose(void);
		inline HexSparseMatrix							transposed(void) const;
};

/* This class is supposed to recalculate its numberOfRows and
 * numberOfColumns so that there is always a non-zero value in
 * its last row and a non-zero value in its last column.
 */

struct HexDecomposition // QR decomposition of matrix, where Q is unitary and R is upper triangular
{
	HexSparseMatrix unitary;		// Q
	HexSparseMatrix triangular;		// R
};

HexSparseMatrix::HexSparseMatrix(void)
{
}

HexSparseMatrix::HexSparseMatrix(const std::vector<std::vector<HexColumnValuePair>>& rows, qint32 noc)
{
	HexSparseMatrix::numberOfRows = static_cast<qint32>(rows.size());
	HexSparseMatrix::numberOfColumns = noc;
	
	HexSparseMatrix::rowOffsets.push_back(0);
	
	for (const auto& row : rows)
	{
		if (not row.empty())
			HexSparseMatrix::pairs.insert(HexSparseMatrix::pairs.end(), row.cbegin(), row.cend());
		
		const auto nextIndex = HexSparseMatrix::rowOffsets.back() + static_cast<qint32>(row.size());
		HexSparseMatrix::rowOffsets.push_back(nextIndex);
	}
}

void HexSparseMatrix::addValue(qint32 row, qint32 column, qreal value)
{
	if (HexSparseMatrix::pairs.empty())
	{
		HexSparseMatrix::rowOffsets = std::vector<qint32>(row + 1, 0);
		HexSparseMatrix::rowOffsets.push_back(1);
			
		HexSparseMatrix::pairs.emplace_back(value, column);
		
		HexSparseMatrix::numberOfRows = row + 1;
		HexSparseMatrix::numberOfColumns = column + 1;
		return;
	}
	
	if (column >= HexSparseMatrix::numberOfColumns)
		HexSparseMatrix::numberOfColumns = column + 1;
	
	if (HexSparseMatrix::numberOfRows <= row) // Here HexSparseMatrix::numberOfRows represents the first row that doesn't exist
	{
		const auto indexOfNewValue = static_cast<qint32>(HexSparseMatrix::pairs.size());
		
		HexSparseMatrix::rowOffsets.insert(HexSparseMatrix::rowOffsets.end(), row - HexSparseMatrix::numberOfRows, indexOfNewValue);
		HexSparseMatrix::rowOffsets.push_back(indexOfNewValue + 1);
		
		HexSparseMatrix::pairs.emplace_back(value, column);
		HexSparseMatrix::numberOfRows = row + 1;
		return;
	}
	
	if (HexSparseMatrix::numberOfRows - 1 == row) // Here HexSparseMatrix::numberOfRows - 1 represents the last row
	{	
		const auto firstPossibleIndex = HexSparseMatrix::rowOffsets.back() - 1;
		const auto lastPossibleIndex = HexSparseMatrix::rowOffsets[row];
		const auto newValueWasInserted = HexSparseMatrix::insertColumnValuePair(firstPossibleIndex, lastPossibleIndex, column, value);
		
		if (newValueWasInserted)
			++HexSparseMatrix::rowOffsets.back();
		
		return;
	}
	
	const auto startIndex = HexSparseMatrix::rowOffsets[row];
	
	if (HexSparseMatrix::rowOffsets[row + 1] == startIndex) // The row is full of zeroes
	{
		for (auto it = HexSparseMatrix::rowOffsets.begin() + row + 1; it != HexSparseMatrix::rowOffsets.end(); ++it)
			++(*it);
		
		HexSparseMatrix::pairs.emplace(HexSparseMatrix::pairs.begin() + startIndex, value, column);
		return;
	}
	
	const auto firstPossibleIndex = HexSparseMatrix::rowOffsets[row + 1] - 1;
	const auto lastPossibleIndex = HexSparseMatrix::rowOffsets[row];
	const auto newValueWasInserted = HexSparseMatrix::insertColumnValuePair(firstPossibleIndex, lastPossibleIndex, column, value);
	
	if (newValueWasInserted)
	{
		for (auto it = HexSparseMatrix::rowOffsets.begin() + row + 1; it != HexSparseMatrix::rowOffsets.end(); ++it)
			++(*it);
	}
}

template<typename Type>
void HexSparseMatrix::Change(std::vector<HexColumnValuePair>& vect, Type beg2, Type end2, qreal coeff)
{
	if (beg2 == end2 or coeff == 0.) // If row is full of zeroes OR if nothing to add
		return;
	
	auto newVect = std::vector<HexColumnValuePair>();
	
	const auto end1 = vect.cend();
	auto cit1 = vect.cbegin();
	
	for (auto cit2 = beg2; cit2 != end2; ++cit2)
	{
		while (cit1 != end1 and cit1->column < cit2->column)
		{
			newVect.emplace_back(cit1->value, cit1->column);
			++cit1;
		}
		
		if (cit1 == end1 or cit1->column != cit2->column)
			newVect.emplace_back(cit2->value*coeff, cit2->column);
		else
		{
			newVect.emplace_back(cit1->value + cit2->value*coeff, cit2->column);
			++cit1;
		}
	}
	
	newVect.insert(newVect.end(), cit1, end1);
	vect.swap(newVect);
}

void HexSparseMatrix::downsize(void)
{
	if (HexSparseMatrix::pairs.empty())
	{
		HexSparseMatrix::numberOfRows = 0;
		HexSparseMatrix::numberOfColumns = 0;
	}
	else
	{
		HexSparseMatrix::updateNumberOfRows();
		HexSparseMatrix::updateNumberOfColumns();
	}
}

/* The hard copy at the end of this function can't really be avoided
 * as you can't populate Q with a new vector until you know for sure
 * that this vector is independent from the previous ones. I chose to
 * populate Q with all the independent vectors that were found at the
 * end of the function, rather than do it one vector at a time at
 * each iteration. Same thing with R and its scalar coeffs.
 */
HexDecomposition HexSparseMatrix::getDecomposition(void) const
{
	auto decomp = HexDecomposition();
	
	if (HexSparseMatrix::pairs.empty())
		return decomp;
	
	const auto transposed = HexSparseMatrix::transposed();
	auto stopIndex = transposed.rowOffsets.cbegin() + 1u;
	
	std::vector<std::vector<HexColumnValuePair>> rowBase;
	std::vector<std::vector<HexColumnValuePair>> coeffs;
	
	for (const auto& startIndex : transposed.rowOffsets)
	{
		auto& rowCoeffs = coeffs.emplace_back();
		
		if (startIndex != *stopIndex) // Row is not full of zeroes, it might be a new independent vector
		{
			const auto start = transposed.pairs.cbegin() + startIndex;
			const auto stop = transposed.pairs.cbegin() + (*stopIndex);
			
			auto newCandidateForBase = std::vector<HexColumnValuePair>(start, stop);
			auto rowIndex = 0;
			
			for (const auto& vct : rowBase)
			{
				const auto scalar = HexSparseMatrix::Scalar(vct.cbegin(), vct.cend(), start, stop);
				
				rowCoeffs.emplace_back(scalar, rowIndex);
				++rowIndex;
			}
			
			auto cit = rowCoeffs.cbegin();
			
			for (const auto& vct : rowBase)
			{
				HexSparseMatrix::Change(newCandidateForBase, vct.cbegin(), vct.cend(), -cit->value);
				++cit;
			}
			
			const auto norm = HexSparseMatrix::Normalise(newCandidateForBase);
			
			if (norm != 0.)
			{
				rowCoeffs.emplace_back(norm, rowIndex);
				rowBase.emplace_back().swap(newCandidateForBase);
			}
		}
		
		++stopIndex;
		
		if (transposed.rowOffsets.cend() == stopIndex)
			break;
	}
	
	decomp.unitary = HexSparseMatrix(rowBase, HexSparseMatrix::numberOfRows).transposed();
	decomp.triangular = HexSparseMatrix(coeffs, HexSparseMatrix::numberOfRows).transposed();
	
	return decomp;
}

/* This function assumes the user knows how to read the vector,
 * meaning they know how many rows and columns this matrix has.
 */
std::vector<qreal> HexSparseMatrix::getDenseMatrix(void) const
{
	std::vector<qreal> matrix;
	matrix.reserve(HexSparseMatrix::numberOfRows*HexSparseMatrix::numberOfColumns);
	
	auto stopIndex = HexSparseMatrix::rowOffsets.cbegin() + 1u;
	
	for (const auto& startIndex : HexSparseMatrix::rowOffsets)
	{
		auto currentColumn = 0;
		
		for (auto index = startIndex; index < *stopIndex; ++index)
		{
			const auto& currentPair = HexSparseMatrix::pairs[index];
			
			while (currentColumn < currentPair.column)
			{
				matrix.push_back(0.);
				++currentColumn;
			}
			
			matrix.push_back(currentPair.value);
			++currentColumn;
		}
		
		while (currentColumn < HexSparseMatrix::numberOfColumns)
		{
			matrix.push_back(0.);
			++currentColumn;
		}
		
		++stopIndex;
		
		if (HexSparseMatrix::rowOffsets.cend() == stopIndex)
			break;
	}
	
	return matrix;
}

qreal HexSparseMatrix::getDensity(void) const
{
	return static_cast<qreal>(HexSparseMatrix::pairs.size())/static_cast<qreal>(HexSparseMatrix::numberOfColumns*HexSparseMatrix::numberOfRows);
}

qint32 HexSparseMatrix::getNumberOfColumns(void) const
{
	return HexSparseMatrix::numberOfColumns;
}

qint32 HexSparseMatrix::getNumberOfRows(void) const
{
	return HexSparseMatrix::numberOfRows;
}

const std::vector<HexColumnValuePair>& HexSparseMatrix::getPairs(void) const
{
	return HexSparseMatrix::pairs;
}

const std::vector<qint32>& HexSparseMatrix::getRowOffsets(void) const
{
	return HexSparseMatrix::rowOffsets;
}

qreal HexSparseMatrix::getSparsity(void) const
{
	return 1. - HexSparseMatrix::getDensity();
}

bool HexSparseMatrix::insertColumnValuePair(qint32 firstPossibleIndex, qint32 lastPossibleIndex, qint32 column, qreal value)
{
	const auto end = HexSparseMatrix::pairs.begin() + lastPossibleIndex - 1;
	auto it = HexSparseMatrix::pairs.begin() + firstPossibleIndex;
	
	while (it != end and it->column > column)
		--it;
	
	if (it == end or it->column != column)
	{
		HexSparseMatrix::pairs.emplace(it + 1u, value, column);
		return true;
	}
	
	it->value = value;
	return false;
}

qint32 HexSparseMatrix::indexOfFreeCell(qint32 row, qint32 column) const
{
	const auto startIndex = HexSparseMatrix::rowOffsets[row];
	const auto stopIndex = HexSparseMatrix::rowOffsets[row + 1];
	
	if (stopIndex - startIndex == HexSparseMatrix::numberOfColumns)
		return -1;
	
	if (startIndex == stopIndex)
		return startIndex;
	
	const auto end = HexSparseMatrix::pairs.cbegin() + stopIndex;
	auto index = startIndex;
	
	for (auto it = HexSparseMatrix::pairs.cbegin() + startIndex; it != end; ++it)
	{
		if (it->column >= column)
			return (it->column == column ? -1 : index);
		
		++index;
	}
	
	return index;
}

bool HexSparseMatrix::insertOne(HexRandomGenerator& generator)
{
	const auto numberOfCells = HexSparseMatrix::numberOfRows*HexSparseMatrix::numberOfColumns;
	const auto numberOfElements = static_cast<qint32>(HexSparseMatrix::pairs.size());
	
	if (numberOfElements >= numberOfCells)
		return false;
	
	auto cell = generator.getNumberWithinRange(numberOfCells);
	auto column = cell % HexSparseMatrix::numberOfColumns;
	auto row = cell/HexSparseMatrix::numberOfColumns;
	
	if (numberOfElements < 1)
		HexSparseMatrix::pairs.emplace_back(1., column);
	else // Very unoptimised for large, dense matrices, but then this project is all about sparse matrices...
	{
		auto index = HexSparseMatrix::indexOfFreeCell(row, column);
		
		while (index < 0)
		{
			cell = generator.getNumberWithinRange(numberOfCells);
			column = cell % HexSparseMatrix::numberOfColumns;
			row = cell/HexSparseMatrix::numberOfColumns;
			
			index = HexSparseMatrix::indexOfFreeCell(row, column);
		}
		
		const auto it = HexSparseMatrix::pairs.begin() + index;
		HexSparseMatrix::pairs.emplace(it, 1., column);
	}
	
	for (auto r = row + 1; r <= HexSparseMatrix::numberOfRows; ++r)
		++HexSparseMatrix::rowOffsets[r];
	
	return true;
}

/* Vectors with very short norms are considered null to avoid
 * numerical instability. I guess 0.001 is still too high though.
 */
qreal HexSparseMatrix::Normalise(std::vector<HexColumnValuePair>& vect)
{
	auto norm = 0.;
	
	for (const auto& pr : vect)
		norm += pr.value*pr.value;
	
	norm = qSqrt(norm);
	
	if (norm < 0.001)
		return 0.;
	
	for (auto& pr : vect)
		pr.value /= norm;

	return norm;
}

void HexSparseMatrix::removeValue(qint32 row, qint32 column)
{
	if (row >= HexSparseMatrix::numberOfRows or column >= HexSparseMatrix::numberOfColumns)
		return;
	
	const auto startIndex = HexSparseMatrix::rowOffsets[row];
	const auto stopIndex = HexSparseMatrix::rowOffsets[row + 1];
	
	if (startIndex == stopIndex)
		return;
	
	const auto iteratorToNextRow = HexSparseMatrix::pairs.begin() + stopIndex;
	
	for (auto it = HexSparseMatrix::pairs.begin() + startIndex; it != iteratorToNextRow; ++it)
	{
		if (it->column < column)
			continue;
		
		if (it->column > column) // It is now clear that there is no value to remove
			return;
		
		HexSparseMatrix::pairs.erase(it);
		break;
	}
	
	for (auto r = row + 1; r <= HexSparseMatrix::numberOfRows; ++r)
		--HexSparseMatrix::rowOffsets[r];
}

template<typename Type1, typename Type2>
void HexSparseMatrix::Rewrite(Type1& it, Type2 beg, Type2 end)
{
	for (auto itt = beg; itt != end; ++itt)
	{
		*it = *itt;
		++it;
	}
}

template<typename Type1, typename Type2>
qreal HexSparseMatrix::Scalar(Type1 beg1, Type1 end1, Type2 beg2, Type2 end2)
{
	if (beg1 == end1 or beg2 == end2) // If either row is full of zeroes
		return 0.;
	
	auto cit2 = beg2;
	auto result = 0.;
	
	for (auto cit1 = beg1; cit1 != end1; ++cit1)
	{
		while (cit2 != end2 and cit2->column < cit1->column)
			++cit2;
		
		if (cit2 == end2)
			break;
		
		if (cit2->column == cit1->column)
			result += cit1->value*cit2->value;
	}
	
	return result;
}

void HexSparseMatrix::setValue(qint32 row, qint32 column, qreal value)
{
	if (row < 0 or column < 0)
		return;
	
	if (value != 0.)
		return HexSparseMatrix::addValue(row, column, value);
	
	HexSparseMatrix::removeValue(row, column);
}

bool HexSparseMatrix::shuffle(HexRandomGenerator& generator)
{
	const auto numberOfElements = HexSparseMatrix::pairs.size();
	
	if (numberOfElements < 1u)
		return false;
	
	auto newRowOffsets = std::vector<qint32>(HexSparseMatrix::numberOfRows + 1, 0);
	auto rowCounts = std::vector<qint32>(HexSparseMatrix::numberOfRows, 0);
	
	for (auto i = 0u; i < numberOfElements; ++i)
	{
		auto row = generator.getNumberWithinRange(HexSparseMatrix::numberOfRows);
		
		while (rowCounts[row] >= HexSparseMatrix::numberOfColumns)
			row = generator.getNumberWithinRange(HexSparseMatrix::numberOfRows);
		
		++rowCounts[row];
	}
	
	for (auto row = 1; row <= HexSparseMatrix::numberOfRows; ++row)
		newRowOffsets[row] = newRowOffsets[row - 1] + rowCounts[row - 1];
	
	HexSparseMatrix::rowOffsets.swap(newRowOffsets);
	generator.shuffle(HexSparseMatrix::pairs);
	
	auto stopIndex = HexSparseMatrix::rowOffsets.cbegin() + 1u;
	
	for (const auto& startIndex : HexSparseMatrix::rowOffsets)
	{
		const auto numberOfElementsInRow = (*stopIndex) - startIndex;
		
		if (numberOfElementsInRow >= 1)
		{
			const auto newColumns = generator.getNumbersWithinRange(numberOfElementsInRow, HexSparseMatrix::numberOfColumns);
			auto it = HexSparseMatrix::pairs.begin() + startIndex;
			
			for (const auto& column : newColumns)
			{
				it->column = column;
				++it;
			}
		}
		
		++stopIndex;
		
		if (HexSparseMatrix::rowOffsets.cend() == stopIndex)
			break;
	}
	
	return true;
}

void HexSparseMatrix::swapColumns(qint32 column1, qint32 column2)
{	
	if (HexSparseMatrix::pairs.empty())
		return;
	
	if (column1 == column2)
		return;
	
	if (column1 < 0 or column2 < 0)
		return;
	
	if (column1 >= HexSparseMatrix::numberOfColumns or column2 >= HexSparseMatrix::numberOfColumns)
		return;
	
	if (column1 > column2)
		std::swap(column1, column2);
	
	auto stopIndex = HexSparseMatrix::rowOffsets.cbegin() + 1u;
	
	for (const auto& startIndex : HexSparseMatrix::rowOffsets) // We have to perform the swap row by row.
	{
		if (*stopIndex != startIndex)
		{
			auto cit = HexSparseMatrix::pairs.cbegin() + startIndex;
			
			auto index1 = -1;
			auto index2 = -1; 
			
			for (auto index = startIndex; index < *stopIndex; ++index)
			{
				if (cit->column == column1)
					index1 = index;
				
				if (cit->column == column2)
					index2 = index;
				
				++cit;
			}
			
			if (index1 >= 0 and index2 >= 0) // Both cells were non-zeroes.
				std::swap(HexSparseMatrix::pairs[index1].value, HexSparseMatrix::pairs[index2].value);
			else if (index1 >= 0) // One cell was non-zero but the other was zero.
			{
				const auto value = HexSparseMatrix::pairs[index1].value;
					
				while (index1 + 1 != *stopIndex and HexSparseMatrix::pairs[index1 + 1].column < column2)
				{
					HexSparseMatrix::pairs[index1] = HexSparseMatrix::pairs[index1 + 1];
					++index1;
				}
				
				HexSparseMatrix::pairs[index1] = HexColumnValuePair(value, column2);
			}
			else if (index2 >= 0) // One cell was non-zero but the other was zero.
			{
				const auto value = HexSparseMatrix::pairs[index2].value;
				
				while (index2 != startIndex and HexSparseMatrix::pairs[index2 - 1].column > column1)
				{
					HexSparseMatrix::pairs[index2] = HexSparseMatrix::pairs[index2 - 1];
					--index2;
				}
				
				HexSparseMatrix::pairs[index2] = HexColumnValuePair(value, column1);
			}
			// Else: both cells are zeroes so there's nothing to swap here.
		}
		
		++stopIndex;
		
		if (HexSparseMatrix::rowOffsets.cend() == stopIndex)
			break;
	}
}

void HexSparseMatrix::swapRows(qint32 row1, qint32 row2)
{
	if (HexSparseMatrix::pairs.empty())
		return;
	
	if (row1 == row2)
		return;
	
	if (row1 < 0 or row2 < 0)
		return;
	
	if (row1 >= HexSparseMatrix::numberOfRows or row2 >= HexSparseMatrix::numberOfRows)
		return;
	
	if (row1 > row2)
		std::swap(row1, row2);
	
	const auto numberOfElementsInRow1 = HexSparseMatrix::rowOffsets[row1 + 1] - HexSparseMatrix::rowOffsets[row1];
	const auto numberOfElementsInRow2 = HexSparseMatrix::rowOffsets[row2 + 1] - HexSparseMatrix::rowOffsets[row2];
	
	if (numberOfElementsInRow1 < 1 and numberOfElementsInRow2 < 1)
		return;
	
	if (numberOfElementsInRow2 <= numberOfElementsInRow1)
	{
		const auto beg1 = HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row1];
		const auto end1 = HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row1 + 1];
		
		const auto beg2 = HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row2];
		const auto end2 = HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row2 + 1];
		
		const auto values1 = std::vector<HexColumnValuePair>(beg1, end1);
		auto it = HexSparseMatrix::pairs.begin() + HexSparseMatrix::rowOffsets[row1];
		
		HexSparseMatrix::Rewrite(it, beg2, end2);
		HexSparseMatrix::Rewrite(it, end1, beg2);
		HexSparseMatrix::Rewrite(it, values1.cbegin(), values1.cend());
	}
	else
	{
		const auto beg1 = std::make_reverse_iterator(HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row1 + 1]);
		const auto end1 = std::make_reverse_iterator(HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row1]);
		
		const auto beg2 = std::make_reverse_iterator(HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row2 + 1]);
		const auto end2 = std::make_reverse_iterator(HexSparseMatrix::pairs.cbegin() + HexSparseMatrix::rowOffsets[row2]);
		
		const auto values2 = std::vector<HexColumnValuePair>(beg2, end2);
		auto it = std::make_reverse_iterator(HexSparseMatrix::pairs.begin() + HexSparseMatrix::rowOffsets[row2 + 1]);
		
		HexSparseMatrix::Rewrite(it, beg1, end1);
		HexSparseMatrix::Rewrite(it, end2, beg1);
		HexSparseMatrix::Rewrite(it, values2.cbegin(), values2.cend());
	}
	
	const auto differenceOfElements = numberOfElementsInRow2 - numberOfElementsInRow1;
	
	for (auto row = row1 + 1; row <= row2; ++row) // There is no difference beyond row2 as the sum of past elements is unchanged.
		HexSparseMatrix::rowOffsets[row] += differenceOfElements;
}

void HexSparseMatrix::transpose(void)
{
	auto newMatrix = HexSparseMatrix::transposed();
	
	HexSparseMatrix::rowOffsets.swap(newMatrix.rowOffsets);
	HexSparseMatrix::pairs.swap(newMatrix.pairs);
	
	std::swap(HexSparseMatrix::numberOfRows, HexSparseMatrix::numberOfColumns);
}

HexSparseMatrix HexSparseMatrix::transposed(void) const
{
	auto rowCounts = std::vector<qint32>(HexSparseMatrix::numberOfColumns, 0);
	
	for (const auto& pr : HexSparseMatrix::pairs)
		++rowCounts[pr.column];
	
	auto transposed = HexSparseMatrix();
	transposed.rowOffsets = std::vector<qint32>(HexSparseMatrix::numberOfColumns + 1, 0);
	
	for (auto row = 1; row <= HexSparseMatrix::numberOfColumns; ++row)
		transposed.rowOffsets[row] = transposed.rowOffsets[row - 1] + rowCounts[row - 1];
	
	transposed.pairs = std::vector<HexColumnValuePair>(HexSparseMatrix::pairs.size());
	auto rowIndexes = transposed.rowOffsets;
	
	auto stopIndex = HexSparseMatrix::rowOffsets.cbegin() + 1u;
	auto row = 0;
	
	for (const auto& startIndex : HexSparseMatrix::rowOffsets)
	{
		for (auto oldIndex = startIndex; oldIndex < *stopIndex; ++oldIndex)
		{
			const auto& pr = HexSparseMatrix::pairs[oldIndex];
			auto& newIndex = rowIndexes[pr.column];
			
			transposed.pairs[newIndex] = HexColumnValuePair(pr.value, row);
			++newIndex; // We need to increment that reference so that, next time another element of the same column is found, it is located right next to the previous one in the vector.
		}
		
		++row;
		++stopIndex;
		
		if (HexSparseMatrix::rowOffsets.cend() == stopIndex)
			break;
	}
	
	transposed.numberOfColumns = HexSparseMatrix::numberOfRows;
	transposed.numberOfRows = HexSparseMatrix::numberOfColumns;
	
	return transposed;
}

void HexSparseMatrix::updateNumberOfColumns(void)
{
	auto newNumberOfColumns = 0;
	
	for (const auto& pr : HexSparseMatrix::pairs)
	{
		if (pr.column + 1 == HexSparseMatrix::numberOfColumns)
			return;
		
		if (newNumberOfColumns <= pr.column)
			newNumberOfColumns = pr.column + 1;
	}
	
	HexSparseMatrix::numberOfColumns = newNumberOfColumns;
}

void HexSparseMatrix::updateNumberOfRows(void)
{
	const auto lastValue = HexSparseMatrix::rowOffsets.back();
	
	while (HexSparseMatrix::rowOffsets.back() == lastValue)
		HexSparseMatrix::rowOffsets.pop_back();
	
	HexSparseMatrix::numberOfRows = static_cast<qint32>(HexSparseMatrix::rowOffsets.size());
	HexSparseMatrix::rowOffsets.push_back(lastValue);
}

#endif
