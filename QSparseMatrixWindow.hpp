#ifndef __Q_SPARSE_MATRIX_WINDOW_HPP__
#define __Q_SPARSE_MATRIX_WINDOW_HPP__

// Qt Libraries
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QMessageBox>
#include <QPushButton>
#include <QRegularExpressionValidator>
#include <QShortcut>
#include <QTextEdit>

// Custom Libraries
#include "HexRandomGenerator.hpp"
#include "HexSparseMatrix.hpp"

class QSparseMatrixWindow : public QMainWindow
{
	Q_OBJECT
	
	private:
	
		static constexpr qint32		CharactersPerField = 7;
		
		inline static QString		CellString(bool);
		inline static QString		NumberString(qreal);
	
		QWidget* const			mainWidget = new QWidget();
		
		QGroupBox* const		matrixBox = new QGroupBox("Sparse Matrix A", mainWidget);
		QTextEdit* const		matrixEdit = new QTextEdit("", mainWidget);
		
		QGroupBox* const		unitaryBox = new QGroupBox("Unitary Matrix Q", mainWidget);
		QTextEdit* const		unitaryEdit = new QTextEdit("", mainWidget);
		
		QGroupBox* const		triangularBox = new QGroupBox("Triangular Matrix R", mainWidget);
		QTextEdit* const		triangularEdit = new QTextEdit("", mainWidget);
		
		QGroupBox* const		dataBox = new QGroupBox("Non-Zero Values", mainWidget);
		QTextEdit* const		dataEdit = new QTextEdit("", mainWidget);
		
		QLineEdit* const		rowSetEdit = new QLineEdit("", mainWidget);
		QLineEdit* const		columnSetEdit = new QLineEdit("", mainWidget);
		QLineEdit* const		valueSetEdit = new QLineEdit("", mainWidget);
		
		QLineEdit* const		swapEdit1 = new QLineEdit("", mainWidget);
		QLineEdit* const		swapEdit2 = new QLineEdit("", mainWidget);
		
		QLineEdit* const		densityEdit = new QLineEdit("-", mainWidget);
		QLineEdit* const		sparsityEdit = new QLineEdit("-", mainWidget);
		QLineEdit* const		rankEdit = new QLineEdit("Rank 0", mainWidget);
		
		QLineEdit* const		rowOffsetsEdit = new QLineEdit("[ ]", mainWidget);
		QLineEdit* const		columnIndicesEdit = new QLineEdit("[ ]", mainWidget);
		QLineEdit* const		valuesEdit = new QLineEdit("[ ]", mainWidget);
		
		HexRandomGenerator		generator;
		HexSparseMatrix			matrix;
		HexDecomposition		decomp;
		
		inline QString			getDataString(const HexSparseMatrix&) const;
		inline QString			getMatrixString(const HexSparseMatrix&, bool) const;
		inline void			updateEntries(void) const;
	
	private slots:
	
		inline void			decompose(void);
		inline void			downsize(void);
		inline void			insertOne(void);
		inline void			setValue(void);
		inline void			shuffle(void);
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
	
	const auto downsizeButton = new QPushButton("Downsize", QSparseMatrixWindow::mainWidget);
	const auto setValueButton = new QPushButton("Set Value", QSparseMatrixWindow::mainWidget);
	const auto swapRowsButton = new QPushButton("Swap Rows", QSparseMatrixWindow::mainWidget);
	const auto swapColumnsButton = new QPushButton("Swap Columns", QSparseMatrixWindow::mainWidget);
	
	const auto insertButton = new QPushButton("Insert 1", QSparseMatrixWindow::mainWidget);
	const auto shuffleButton = new QPushButton("Shuffle", QSparseMatrixWindow::mainWidget);
	const auto transposeButton = new QPushButton("Transpose", QSparseMatrixWindow::mainWidget);
	const auto decomposeButton = new QPushButton("Decompose", QSparseMatrixWindow::mainWidget);
	
	const auto rowLabel = new QLabel("Row Offsets", QSparseMatrixWindow::mainWidget);
	const auto columnLabel = new QLabel("Column Indices", QSparseMatrixWindow::mainWidget);
	const auto valueLabel = new QLabel("Matrix Values", QSparseMatrixWindow::mainWidget);
	
	const auto layout1 = new QHBoxLayout();
	layout1->addWidget(QSparseMatrixWindow::matrixEdit);
	
	QSparseMatrixWindow::matrixEdit->setFrameStyle(QFrame::NoFrame);
	QSparseMatrixWindow::matrixBox->setLayout(layout1);
	
	const auto layout2 = new QHBoxLayout();
	layout2->addWidget(QSparseMatrixWindow::unitaryEdit);
	
	QSparseMatrixWindow::unitaryEdit->setFrameStyle(QFrame::NoFrame);
	QSparseMatrixWindow::unitaryBox->setLayout(layout2);
	
	const auto layout3 = new QHBoxLayout();
	layout3->addWidget(QSparseMatrixWindow::triangularEdit);
	
	QSparseMatrixWindow::triangularEdit->setFrameStyle(QFrame::NoFrame);
	QSparseMatrixWindow::triangularBox->setLayout(layout3);
	
	const auto layout4 = new QHBoxLayout();
	layout4->addWidget(QSparseMatrixWindow::dataEdit);
	
	QSparseMatrixWindow::dataEdit->setFrameStyle(QFrame::NoFrame);
	QSparseMatrixWindow::dataBox->setLayout(layout4);
	
	const auto grid = new QGridLayout();
	auto rowCount = 0;
	
	grid->addWidget(QSparseMatrixWindow::matrixBox, rowCount, 0, 1, 3);
	grid->addWidget(QSparseMatrixWindow::dataBox, rowCount, 3, 3, 1);
	
	grid->addWidget(QSparseMatrixWindow::unitaryBox, ++rowCount, 0, 1, 3);
	grid->addWidget(QSparseMatrixWindow::triangularBox, ++rowCount, 0, 1, 3);
	
	grid->addWidget(QSparseMatrixWindow::rankEdit, ++rowCount, 0, 1, 1);
	grid->addWidget(QSparseMatrixWindow::sparsityEdit, rowCount, 1, 1, 1);
	grid->addWidget(QSparseMatrixWindow::densityEdit, rowCount, 2, 1, 1);
	grid->addWidget(downsizeButton, rowCount, 3, 1, 1);
	
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
	
	grid->addWidget(transposeButton, ++rowCount, 0, 1, 1);
	grid->addWidget(insertButton, rowCount, 1, 1, 1);
	grid->addWidget(shuffleButton, rowCount, 2, 1, 1);
	grid->addWidget(decomposeButton, rowCount, 3, 1, 1);
	
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
	const auto edits = { QSparseMatrixWindow::rowOffsetsEdit, QSparseMatrixWindow::columnIndicesEdit, QSparseMatrixWindow::valuesEdit, QSparseMatrixWindow::rankEdit, QSparseMatrixWindow::sparsityEdit, QSparseMatrixWindow::densityEdit };
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
	
	QObject::connect(downsizeButton, SIGNAL(clicked(void)), this, SLOT(downsize(void)));
	QObject::connect(setValueButton, SIGNAL(clicked(void)), this, SLOT(setValue(void)));
	QObject::connect(swapRowsButton, SIGNAL(clicked(void)), this, SLOT(swapRows(void)));
	QObject::connect(swapColumnsButton, SIGNAL(clicked(void)), this, SLOT(swapColumns(void)));
	
	QObject::connect(insertButton, SIGNAL(clicked(void)), this, SLOT(insertOne(void)));
	QObject::connect(shuffleButton, SIGNAL(clicked(void)), this, SLOT(shuffle(void)));
	QObject::connect(transposeButton, SIGNAL(clicked(void)), this, SLOT(transpose(void)));
	QObject::connect(decomposeButton, SIGNAL(clicked(void)), this, SLOT(decompose(void)));
	
	const auto s1 = new QShortcut(mainWidget);
	s1->setKeys({ QKeySequence(Qt::Key_Return), QKeySequence(Qt::Key_Enter) });
	
	const auto s2 = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_D), mainWidget);
	const auto s3 = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_I), mainWidget);
	const auto s4 = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_S), mainWidget);
	const auto s5 = new QShortcut(QKeySequence(Qt::CTRL | Qt::Key_T), mainWidget);
	
	QObject::connect(s1, SIGNAL(activated(void)), this, SLOT(setValue(void)));
	QObject::connect(s2, SIGNAL(activated(void)), this, SLOT(decompose(void)));
	QObject::connect(s3, SIGNAL(activated(void)), this, SLOT(insertOne(void)));
	QObject::connect(s4, SIGNAL(activated(void)), this, SLOT(shuffle(void)));
	QObject::connect(s5, SIGNAL(activated(void)), this, SLOT(transpose(void)));
}

QString QSparseMatrixWindow::CellString(bool zero)
{
	const auto theChar = (zero ? '0' : 'x');
	const auto freeSpace = QSparseMatrixWindow::CharactersPerField - 1;
	
	return QString("&nbsp;").repeated(freeSpace/2) + theChar + QString("&nbsp;").repeated(freeSpace - freeSpace/2);
}

void QSparseMatrixWindow::decompose(void)
{
	QSparseMatrixWindow::decomp = QSparseMatrixWindow::matrix.getDecomposition();
	
	const auto unitaryString = QSparseMatrixWindow::getMatrixString(QSparseMatrixWindow::decomp.unitary, false);
	const auto triangularString = QSparseMatrixWindow::getMatrixString(QSparseMatrixWindow::decomp.triangular, true);
	
	QSparseMatrixWindow::unitaryEdit->setHtml(unitaryString);
	QSparseMatrixWindow::triangularEdit->setHtml(triangularString);
	
	const auto rank = decomp.unitary.getHighestColumn() + 1;
	const auto rankString = "Rank " + QString::number(rank);
	
	const auto unitaryDimensionString = QSparseMatrixWindow::decomp.unitary.getDimensionString();
	const auto triangularDimensionString = QSparseMatrixWindow::decomp.triangular.getDimensionString();
	
	const auto unitaryTitle = "Unitary Matrix Q (" + unitaryDimensionString + ')';
	const auto triangularTitle = "Triangular Matrix R (" + triangularDimensionString + ')';
	
	QSparseMatrixWindow::rankEdit->setText(rankString);
	QSparseMatrixWindow::unitaryBox->setTitle(unitaryTitle);
	QSparseMatrixWindow::triangularBox->setTitle(triangularTitle);
}

void QSparseMatrixWindow::downsize(void)
{
	QSparseMatrixWindow::matrix.downsize();
	QSparseMatrixWindow::updateEntries();
}

QString QSparseMatrixWindow::getDataString(const HexSparseMatrix& matrix) const
{
	const auto& rowOffsets = matrix.getRowOffsets();
	const auto& pairs = matrix.getPairs();
	
	auto stopIndex = rowOffsets.cbegin() + 1u;
	auto str = QString();
	auto row = 0;
	
	for (const auto& startIndex : rowOffsets)
	{
		const auto beginning = "a(" + QString::number(row);
		const auto end = pairs.cbegin() + (*stopIndex);
		
		for (auto it = pairs.cbegin() + startIndex; it != end; ++it)
			str += beginning  + ", " + QString::number(it->column) + ") = " + QString::number(it->value) + '\n';
		
		++row;
		++stopIndex;
		
		if (rowOffsets.cend() == stopIndex)
			break;
	}
	
	return str;
}

QString QSparseMatrixWindow::getMatrixString(const HexSparseMatrix& matrix, bool showAllZeroes) const
{
	const auto& rowOffsets = matrix.getRowOffsets();
	const auto& pairs = matrix.getPairs();
	
	const auto numberOfColumns = matrix.getNumberOfColumns();
	const auto highestColumn = matrix.getHighestColumn();
	
	const auto zeroString = QSparseMatrixWindow::CellString(true) + ' ';
	const auto xString = QSparseMatrixWindow::CellString(showAllZeroes) + ' ';
	
	auto stopIndex = rowOffsets.cbegin() + 1u;
	auto str = QString();
	
	for (const auto& startIndex : rowOffsets)
	{
		const auto end = pairs.cbegin() + (*stopIndex);
		auto currentColumn = 0;
		
		for (auto it = pairs.cbegin() + startIndex; it != end; ++it)
		{
			while (currentColumn < it->column)
			{
				str += QSparseMatrixWindow::NumberString(0.) + ' ';
				++currentColumn;
			}
			
			str += QSparseMatrixWindow::NumberString(it->value) + ' ';
			++currentColumn;
		}
		
		while (currentColumn <= highestColumn)
		{
			str += zeroString;
			++currentColumn;
		}
		
		while (currentColumn < numberOfColumns)
		{
			str += xString;
			++currentColumn;
		}
		
		str.chop(1u);
		str += "<br>";
	
		++stopIndex;
		
		if (rowOffsets.cend() == stopIndex)
			break;
	}
	
	return str;
}

void QSparseMatrixWindow::insertOne(void)
{
	const auto success = QSparseMatrixWindow::matrix.insertOne(QSparseMatrixWindow::generator);
	
	if (success)
		QSparseMatrixWindow::updateEntries();
	else
		QMessageBox::critical(this, "Nope", "There is no available field.");
}

QString QSparseMatrixWindow::NumberString(qreal val)
{
	auto str = QString::number(val, 'f', 6);
	const auto pos = str.indexOf('.');
	
	if (pos < 0)
		return "<font color=#FF0000>" + str + "</font>";

	if (pos > QSparseMatrixWindow::CharactersPerField)
		return "<font color=#FF0000>" + str.sliced(0, pos) + "</font>";
	
	while (str.size() > 0 and (str.back() == '.' or str.back() == '0'))
		str.chop(1u);
	
	while (str.size() > QSparseMatrixWindow::CharactersPerField)
		str.chop(1u);
	
	if (str.isEmpty() or str == "-")
	{
		const auto freeSpace = QSparseMatrixWindow::CharactersPerField - 1;
		return QString("&nbsp;").repeated(freeSpace/2) + '0' + QString("&nbsp;").repeated(freeSpace - freeSpace/2);
	}
	
	const auto definitiveSize = str.size();
	const auto freeSpace = QSparseMatrixWindow::CharactersPerField - definitiveSize;
	
	return QString("&nbsp;").repeated(freeSpace/2) + "<font color=#FF0000>" + str + "</font>" + QString("&nbsp;").repeated(freeSpace - freeSpace/2);
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

void QSparseMatrixWindow::shuffle(void)
{
	const auto success = QSparseMatrixWindow::matrix.shuffle(QSparseMatrixWindow::generator);
	
	if (success)
		QSparseMatrixWindow::updateEntries();
	else
		QMessageBox::critical(this, "Nope", "There is nothing to shuffle.");
}

void QSparseMatrixWindow::swapColumns(void)
{
	if (QSparseMatrixWindow::swapEdit1->text().isEmpty() or QSparseMatrixWindow::swapEdit2->text().isEmpty())
		return;
	
	const auto column1 = QSparseMatrixWindow::swapEdit1->text().toInt();
	const auto column2 = QSparseMatrixWindow::swapEdit2->text().toInt();
	
	QSparseMatrixWindow::matrix.swapColumns(column1, column2);
	QSparseMatrixWindow::updateEntries();
}

void QSparseMatrixWindow::swapRows(void)
{
	if (QSparseMatrixWindow::swapEdit1->text().isEmpty() or QSparseMatrixWindow::swapEdit2->text().isEmpty())
		return;
	
	const auto row1 = QSparseMatrixWindow::swapEdit1->text().toInt();
	const auto row2 = QSparseMatrixWindow::swapEdit2->text().toInt();
	
	QSparseMatrixWindow::matrix.swapRows(row1, row2);
	QSparseMatrixWindow::updateEntries();
}

void QSparseMatrixWindow::transpose(void)
{
	QSparseMatrixWindow::matrix.transpose();
	QSparseMatrixWindow::updateEntries();
}

void QSparseMatrixWindow::updateEntries(void) const
{
	const auto dataString = QSparseMatrixWindow::getDataString(QSparseMatrixWindow::matrix);
	const auto matrixString = QSparseMatrixWindow::getMatrixString(QSparseMatrixWindow::matrix, true);
	
	QSparseMatrixWindow::dataEdit->setText(dataString);
	QSparseMatrixWindow::matrixEdit->setHtml(matrixString);
	
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
	const auto densityString = (pairs.empty() ? "-" : "Density " + QString::number(QSparseMatrixWindow::matrix.getDensity()*100., 'f', 2) + '%');
	const auto sparsityString = (pairs.empty() ? "-" : "Sparsity " + QString::number(QSparseMatrixWindow::matrix.getSparsity()*100., 'f', 2) + '%');
	
	QSparseMatrixWindow::rankEdit->setText(rankString);
	QSparseMatrixWindow::densityEdit->setText(densityString);
	QSparseMatrixWindow::sparsityEdit->setText(sparsityString);
	
	if (pairs.empty())
	{
		QSparseMatrixWindow::matrixBox->setTitle("Sparse Matrix A");
		QSparseMatrixWindow::dataBox->setTitle("Non-Zero Values");
	}
	else
	{
		const auto matrixTitle = "Sparse Matrix A (" + QSparseMatrixWindow::matrix.getDimensionString() + ')';
		const auto dataTitle = "Non-Zero Values (" + QString::number(pairs.size()) + ')';
		
		QSparseMatrixWindow::matrixBox->setTitle(matrixTitle);
		QSparseMatrixWindow::dataBox->setTitle(dataTitle);
	}
	
	QSparseMatrixWindow::unitaryBox->setTitle("Unitary Matrix Q");
	QSparseMatrixWindow::triangularBox->setTitle("Triangular Matrix R");
}

#endif
