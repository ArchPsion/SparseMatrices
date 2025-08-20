#ifndef __Q_SPARSE_MATRIX_WINDOW_HPP__
#define __Q_SPARSE_MATRIX_WINDOW_HPP__

// Qt Libraries
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QPushButton>
#include <QGroupBox>
#include <QRegularExpressionValidator>
#include <QTextEdit>

// Custom Libraries
#include "HexSparseMatrix.hpp"

class QSparseMatrixWindow : public QMainWindow
{
	Q_OBJECT
	
	private:
	
		inline static void		FillEntryWithMatrix(const HexSparseMatrix&, QTextEdit*, QString&&, qint32, bool);
		inline static QString		NumberString(qreal);
	
		QWidget* const			mainWidget = new QWidget();
		
		QTextEdit* const		matrixEdit = new QTextEdit("", mainWidget);
		QTextEdit* const		unitaryEdit = new QTextEdit("", mainWidget);
		QTextEdit* const		triangularEdit = new QTextEdit("", mainWidget);
		QTextEdit* const		dataEdit = new QTextEdit("", mainWidget);
		
		QLineEdit* const		rowSetEdit = new QLineEdit("", mainWidget);
		QLineEdit* const		columnSetEdit = new QLineEdit("", mainWidget);
		QLineEdit* const		valueSetEdit = new QLineEdit("", mainWidget);
		
		QLineEdit* const		swapEdit1 = new QLineEdit("", mainWidget);
		QLineEdit* const		swapEdit2 = new QLineEdit("", mainWidget);
		
		QLineEdit* const		dimensionEdit = new QLineEdit("0 × 0", mainWidget);
		QLineEdit* const		rankEdit = new QLineEdit("Rank 0", mainWidget);
		QLineEdit* const		densityEdit = new QLineEdit("-", mainWidget);
		QLineEdit* const		sparsityEdit = new QLineEdit("-", mainWidget);
		
		QLineEdit* const		rowOffsetsEdit = new QLineEdit("[ ]", mainWidget);
		QLineEdit* const		columnIndicesEdit = new QLineEdit("[ ]", mainWidget);
		QLineEdit* const		valuesEdit = new QLineEdit("[ ]", mainWidget);
		
		HexSparseMatrix			matrix;
		HexDecomposition		decomp;
		
		inline void			updateEntries(void) const;
	
	private slots:
	
		inline void			decompose(void);
		inline void			setValue(void);
		inline void			swapColumns(void);
		inline void			swapRows(void);
		inline void			transpose(void);
	
	public:
	
		inline				QSparseMatrixWindow(void);
};

QSparseMatrixWindow::QSparseMatrixWindow(void) : QMainWindow()
{
	QMainWindow::setCentralWidget(QSparseMatrixWindow::mainWidget);
	QMainWindow::setWindowTitle("Sparse Matrix CSR Interface");
	QMainWindow::setMinimumWidth(400);
	
	const auto setValueButton = new QPushButton("Set Value", QSparseMatrixWindow::mainWidget);
	const auto swapRowsButton = new QPushButton("Swap Rows", QSparseMatrixWindow::mainWidget);
	const auto swapColumnsButton = new QPushButton("Swap Columns", QSparseMatrixWindow::mainWidget);
	const auto transposeButton = new QPushButton("Transpose", QSparseMatrixWindow::mainWidget);
	const auto decomposeButton = new QPushButton("Decompose", QSparseMatrixWindow::mainWidget);
	
	const auto rowLabel = new QLabel("Row Offsets", QSparseMatrixWindow::mainWidget);
	const auto columnLabel = new QLabel("Column Indices", QSparseMatrixWindow::mainWidget);
	const auto valueLabel = new QLabel("Matrix Values", QSparseMatrixWindow::mainWidget);
	
	const auto box1 = new QGroupBox("Sparse Matrix A");
	const auto layout1 = new QHBoxLayout();
	
	layout1->addWidget(QSparseMatrixWindow::matrixEdit);
	QSparseMatrixWindow::matrixEdit->setFrameStyle(QFrame::NoFrame);
	box1->setLayout(layout1);
	
	const auto box2 = new QGroupBox("Unitary Matrix Q");
	const auto layout2 = new QHBoxLayout();
	
	layout2->addWidget(QSparseMatrixWindow::unitaryEdit);
	QSparseMatrixWindow::unitaryEdit->setFrameStyle(QFrame::NoFrame);
	box2->setLayout(layout2);
	
	const auto box3 = new QGroupBox("Triangular Matrix R");
	const auto layout3 = new QHBoxLayout();
	
	layout3->addWidget(QSparseMatrixWindow::triangularEdit);
	QSparseMatrixWindow::triangularEdit->setFrameStyle(QFrame::NoFrame);
	box3->setLayout(layout3);
	
	const auto box4 = new QGroupBox("Non-Zero Values");
	const auto layout4 = new QHBoxLayout();
	
	layout4->addWidget(QSparseMatrixWindow::dataEdit);
	QSparseMatrixWindow::dataEdit->setFrameStyle(QFrame::NoFrame);
	box4->setLayout(layout4);
	
	const auto grid = new QGridLayout();
	auto rowCount = 0;
	
	grid->addWidget(box1, rowCount, 0, 1, 3);
	grid->addWidget(box4, rowCount, 3, 3, 1);
	
	grid->addWidget(box2, ++rowCount, 0, 1, 3);
	grid->addWidget(box3, ++rowCount, 0, 1, 3);
	
	grid->addWidget(QSparseMatrixWindow::dimensionEdit, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::rankEdit, rowCount, 1, 1, 1);
	grid->addWidget(QSparseMatrixWindow::sparsityEdit, rowCount, 2, 1, 1);
	grid->addWidget(QSparseMatrixWindow::densityEdit, rowCount, 3, 1, 1);
	
	grid->addWidget(rowLabel, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::rowOffsetsEdit, rowCount, 1, 1, 3);
	grid->addWidget(columnLabel, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::columnIndicesEdit, rowCount, 1, 1, 3);
	grid->addWidget(valueLabel, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::valuesEdit, rowCount, 1, 1, 3);
	
	grid->addWidget(QSparseMatrixWindow::rowSetEdit, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::columnSetEdit, rowCount, 1, 1, 1);
	grid->addWidget(QSparseMatrixWindow::valueSetEdit, rowCount, 2, 1, 1);
	grid->addWidget(setValueButton, rowCount, 3, 1, 1);
	
	grid->addWidget(QSparseMatrixWindow::swapEdit1, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::swapEdit2, rowCount, 1, 1, 1);
	grid->addWidget(swapRowsButton, rowCount, 2, 1, 1);
	grid->addWidget(swapColumnsButton, rowCount, 3, 1, 1);
	
	grid->addWidget(transposeButton, ++rowCount, 0, 1, 2);
	grid->addWidget(decomposeButton, rowCount, 2, 1, 2);
	
	grid->setRowStretch(0, 100);
	grid->setRowStretch(1, 100);
	grid->setRowStretch(2, 100);
	grid->setColumnStretch(1, 100);
	grid->setColumnStretch(2, 100);
	
	QSparseMatrixWindow::mainWidget->setLayout(grid);
	
	const auto realValidator = new QRegularExpressionValidator(QRegularExpression("\\-?[0-9]*\\.?[0-9]*"), this);
	const auto integerValidator = new QRegularExpressionValidator(QRegularExpression("[0-9]*"), this);
	
	QSparseMatrixWindow::matrixEdit->setReadOnly(true);
	QSparseMatrixWindow::unitaryEdit->setReadOnly(true);
	QSparseMatrixWindow::triangularEdit->setReadOnly(true);
	QSparseMatrixWindow::dataEdit->setReadOnly(true);
	
	QSparseMatrixWindow::valueSetEdit->setAlignment(Qt::AlignCenter);
	QSparseMatrixWindow::valueSetEdit->setValidator(realValidator);
	
	const auto entries = { QSparseMatrixWindow::rowSetEdit, QSparseMatrixWindow::columnSetEdit, QSparseMatrixWindow::swapEdit1, QSparseMatrixWindow::swapEdit2 };
	const auto edits = { QSparseMatrixWindow::rowOffsetsEdit, QSparseMatrixWindow::columnIndicesEdit, QSparseMatrixWindow::valuesEdit, QSparseMatrixWindow::dimensionEdit, QSparseMatrixWindow::rankEdit, QSparseMatrixWindow::sparsityEdit, QSparseMatrixWindow::densityEdit };
	const auto labels = { rowLabel, columnLabel, valueLabel };
	
	for (const auto& entry : entries)
	{
		entry->setAlignment(Qt::AlignCenter);
		entry->setValidator(integerValidator);
	}
	
	for (const auto& edit : edits)
	{
		edit->setReadOnly(true);
		edit->setAlignment(Qt::AlignCenter);
		edit->setObjectName("A");
	}
	
	for (const auto& label : labels)
		label->setAlignment(Qt::AlignCenter);
	
	QObject::connect(setValueButton, SIGNAL(clicked(void)), this, SLOT(setValue(void)));
	QObject::connect(swapRowsButton, SIGNAL(clicked(void)), this, SLOT(swapRows(void)));
	QObject::connect(swapColumnsButton, SIGNAL(clicked(void)), this, SLOT(swapColumns(void)));
	QObject::connect(transposeButton, SIGNAL(clicked(void)), this, SLOT(transpose(void)));
	QObject::connect(decomposeButton, SIGNAL(clicked(void)), this, SLOT(decompose(void)));
}

void QSparseMatrixWindow::decompose(void)
{
	QSparseMatrixWindow::decomp = QSparseMatrixWindow::matrix.getDecomposition();
	
	QSparseMatrixWindow::FillEntryWithMatrix(decomp.unitary, QSparseMatrixWindow::unitaryEdit, "unused orthonormal column(s)", QSparseMatrixWindow::matrix.getNumberOfRows(), true);
	QSparseMatrixWindow::FillEntryWithMatrix(decomp.triangular, QSparseMatrixWindow::triangularEdit, "row(s) of zeroes", QSparseMatrixWindow::matrix.getNumberOfRows(), false);
	
	const auto rank = decomp.unitary.getNumberOfColumns();
	const auto rankString = "Rank " + QString::number(rank);
	
	QSparseMatrixWindow::rankEdit->setText(rankString);
}

void QSparseMatrixWindow::FillEntryWithMatrix(const HexSparseMatrix& matrix, QTextEdit* edit, QString&& str, qint32 biggerNumber, bool column)
{
	const auto numberOfColumns = matrix.getNumberOfColumns();
	const auto numberOfRows = matrix.getNumberOfRows();
	const auto denseMatrix = matrix.getDenseMatrix();
	
	auto matrixString = QString();
	auto index = 0;
	
	for (auto cit = denseMatrix.cbegin(); denseMatrix.cend() != cit; ++cit)
	{
		matrixString += QSparseMatrixWindow::NumberString(*cit);
		matrixString += ((index + 1) % numberOfColumns != 0 ? ' ' : '\n');
		
		++index;
	}
	
	const auto smallerNumber = (column ? numberOfColumns : numberOfRows);
	const auto diff = biggerNumber - smallerNumber;
	
	if (diff != 0)
		matrixString += "and " + QString::number(diff) + " more " + str;
	
	edit->setText(matrixString);
}

QString QSparseMatrixWindow::NumberString(qreal val)
{
	auto str = QString::number(val, 'f', 6);
	const auto pos = str.indexOf('.');
	
	if (pos > 3)
		return str.sliced(0, pos);
	
	while (str.size() > 4 or (str.size() > 0 and (str.back() == '.' or str.back() == '0')))
		str.chop(1u);
	
	switch (str.size())
	{
		case 0u:
			return " 0  ";
		
		case 1u:
			return ' ' + str + "  ";
		
		case 2u:
			return ' ' + str + ' ';
		
		case 3u:
			return str + ' ';
		
		default:
			break;
	}
	
	return str;
}

void QSparseMatrixWindow::setValue(void)
{
	if (QSparseMatrixWindow::rowSetEdit->text().isEmpty() or QSparseMatrixWindow::columnSetEdit->text().isEmpty() or QSparseMatrixWindow::valueSetEdit->text().isEmpty())
		return;
	
	const auto row = QSparseMatrixWindow::rowSetEdit->text().toInt();
	const auto column = QSparseMatrixWindow::columnSetEdit->text().toInt();
	const auto value = QSparseMatrixWindow::valueSetEdit->text().toDouble();
	
	QSparseMatrixWindow::matrix.setValue(row, column, value);
	QSparseMatrixWindow::updateEntries();
	
	QSparseMatrixWindow::rowSetEdit->clear();
	QSparseMatrixWindow::columnSetEdit->clear();
	QSparseMatrixWindow::valueSetEdit->clear();
}

void QSparseMatrixWindow::swapColumns(void)
{
	if (QSparseMatrixWindow::swapEdit1->text().isEmpty() or QSparseMatrixWindow::swapEdit2->text().isEmpty())
		return;
	
	const auto column1 = QSparseMatrixWindow::swapEdit1->text().toInt();
	const auto column2 = QSparseMatrixWindow::swapEdit2->text().toInt();
	
	QSparseMatrixWindow::matrix.swapColumns(column1, column2);
	QSparseMatrixWindow::updateEntries();
	
	QSparseMatrixWindow::swapEdit1->clear();
	QSparseMatrixWindow::swapEdit2->clear();
}

void QSparseMatrixWindow::swapRows(void)
{
	if (QSparseMatrixWindow::swapEdit1->text().isEmpty() or QSparseMatrixWindow::swapEdit2->text().isEmpty())
		return;
	
	const auto row1 = QSparseMatrixWindow::swapEdit1->text().toInt();
	const auto row2 = QSparseMatrixWindow::swapEdit2->text().toInt();
	
	QSparseMatrixWindow::matrix.swapRows(row1, row2);
	QSparseMatrixWindow::updateEntries();
	
	QSparseMatrixWindow::swapEdit1->clear();
	QSparseMatrixWindow::swapEdit2->clear();
}

void QSparseMatrixWindow::transpose(void)
{
	QSparseMatrixWindow::matrix.transpose();
	QSparseMatrixWindow::updateEntries();
}

void QSparseMatrixWindow::updateEntries(void) const
{
	const auto denseMatrix = QSparseMatrixWindow::matrix.getDenseMatrix();
	const auto numberOfRows = QSparseMatrixWindow::matrix.getNumberOfRows();
	const auto numberOfColumns = QSparseMatrixWindow::matrix.getNumberOfColumns();
	
	auto matrixString = QString();
	auto dataString = QString();
	auto index = 0;
	
	for (auto cit = denseMatrix.cbegin(); denseMatrix.cend() != cit; ++cit)
	{
		matrixString += QSparseMatrixWindow::NumberString(*cit);
		matrixString += ((index + 1) % numberOfColumns != 0 ? ' ' : '\n');
		
		if (*cit != 0.)
			dataString += "a(" + QString::number(index/numberOfColumns) + ", " + QString::number(index % numberOfColumns) + ") = " + QString::number(*cit) + '\n';
		
		++index;
	}
	
	dataString.chop(1u);
	
	QSparseMatrixWindow::dataEdit->setText(dataString);
	QSparseMatrixWindow::matrixEdit->setText(matrixString);
	QSparseMatrixWindow::unitaryEdit->setText("");
	QSparseMatrixWindow::triangularEdit->setText("");
	
	const auto& rowOffsets = QSparseMatrixWindow::matrix.getRowOffsets();
	const auto& pairs = QSparseMatrixWindow::matrix.getPairs();
	
	auto rowOffsetsString = QString("[");
	auto columnIndicesString = QString("[");
	auto valuesString = QString("[");
	
	for (const auto& offset : rowOffsets)
		rowOffsetsString += ' ' + QString::number(offset);
	
	for (const auto& pr : pairs)
	{
		columnIndicesString += ' ' + QString::number(pr.column);
		valuesString += ' ' + QString::number(pr.value);
	}
	
	QSparseMatrixWindow::rowOffsetsEdit->setText(rowOffsetsString + " ]");
	QSparseMatrixWindow::columnIndicesEdit->setText(columnIndicesString + " ]");
	QSparseMatrixWindow::valuesEdit->setText(valuesString + " ]");
	
	const auto rankString = QString(pairs.empty() ? "Rank 0" : "-");
	const auto dimensionString = QString::number(numberOfRows) + " × " + QString::number(numberOfColumns);
	const auto densityString = (pairs.empty() ? "-" : QString::number(QSparseMatrixWindow::matrix.getDensity()*100., 'f', 2));
	const auto sparsityString = (pairs.empty() ? "-" : QString::number(QSparseMatrixWindow::matrix.getSparsity()*100., 'f', 2));
	
	QSparseMatrixWindow::rankEdit->setText(rankString);
	QSparseMatrixWindow::dimensionEdit->setText(dimensionString);
	QSparseMatrixWindow::densityEdit->setText(densityString);
	QSparseMatrixWindow::sparsityEdit->setText(sparsityString);
}

#endif
