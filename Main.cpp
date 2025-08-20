// Qt Libraries
#include <QApplication>

// Custom Libraries
#include "QSparseMatrixWindow.hpp"

int main(int argc, char* argv[])
{
	const QString stylesheet = "QGroupBox { border: 1px solid gray; border-radius: 9px; font-size: 12px; font-weight: bold; margin-top: 1.5ex; }"
					"QLineEdit#A   { font-family: monospace; background-color: rgb(204, 204, 255) }"
					"QTextEdit   { font-family: monospace }";
	
	auto app = QApplication(argc, argv);
	auto palette = app.palette();
	
	palette.setColor(QPalette::Window, QColor(255, 255, 255));
	palette.setColor(QPalette::WindowText, QColor(51, 51, 51));
	
	palette.setColor(QPalette::Base, QColor(255, 255, 255));
	palette.setColor(QPalette::AlternateBase, QColor(255, 255, 255));
	
	palette.setColor(QPalette::PlaceholderText, QColor(51, 51, 51));
	palette.setColor(QPalette::Text, QColor(51, 51, 51));
	
	palette.setColor(QPalette::Disabled, QPalette::Button, QColor(170, 170, 170));
	palette.setColor(QPalette::Active, QPalette::Button, QColor(204, 204, 204));
	palette.setColor(QPalette::Inactive, QPalette::Button, QColor(204, 204, 204));
	palette.setColor(QPalette::ButtonText, QColor(51, 51, 51));
	
	palette.setColor(QPalette::Highlight, QColor(142, 45, 197));
	palette.setColor(QPalette::HighlightedText, QColor(255, 255, 255));
	
	app.setStyleSheet(stylesheet);
	app.setPalette(palette);
	
	auto window = QSparseMatrixWindow();
	
	/*
	auto matrix = HexSparseMatrix();
	
	matrix.setValue(1, 1, 5.);
	matrix.setValue(1, 1, 3.);
	matrix.show();
	matrix.setValue(3, 4, -2.);
	matrix.show();
	matrix.setValue(2, 4, -23.);
	matrix.show();
	matrix.setValue(2, 1, 7.);
	matrix.show();
	matrix.setValue(0, 8, 123.);
	matrix.show();
	matrix.setValue(3, 6, -45.32);
	matrix.show();
	matrix.swapColumns(1, 4);
	matrix.show();
	
	matrix.setValue(2, 8, 1.);
	matrix.show();
	matrix.setValue(3, 4, 0.);
	matrix.show();
	matrix.setValue(3, 6, 0.);
	matrix.show();
	matrix.setValue(0, 8, 0.);
	matrix.show();
	matrix.setValue(2, 8, 0.);
	matrix.show();
	matrix.setValue(2, 4, 0.);
	matrix.show();
	matrix.swapRows(1, 0);
	matrix.show();
	matrix.setValue(0, 4, 0.);
	matrix.show();
	*/
	
	window.show();
	return app.exec();
}
