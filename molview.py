from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
from PyQt5.QtGui import QPixmap, QPalette, QColor, QImage, QPainter
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QFileDialog, QHBoxLayout, QLabel, QDialog, QTextEdit, QTabWidget
from PyQt5.QtSvg import QSvgWidget, QSvgRenderer
from PyQt5.QtCore import Qt, QByteArray

class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        self.currentMolecule = Chem.MolFromSmiles("C(CC(=O)NC(CS)C(=O)NCC(=O)O)C(C(=O)O)N")

        # set window title
        self.setWindowTitle("Molecule Viewer")
        self.resize(800, 600)

        # create central widget
        self.central_widget = QWidget(self)

        # create canvas
        self.canvas = QLabel(self.central_widget)
        self.canvas.setAlignment(Qt.AlignCenter)
        self.canvas.setFixedSize(600, 600)

        # set white background
        pal = QPalette()
        pal.setColor(QPalette.Background, QColor(255, 255, 255))
        self.canvas.setAutoFillBackground(True)
        self.canvas.setPalette(pal)

        # create properties tab widget
        self.molecularWeightLabel = QLabel("Average mass:")
        self.monoisotopicMwLabel = QLabel("Monoisotopic mass: ")
        self.properties_tab = QWidget(self.central_widget)
        self.properties_layout = QVBoxLayout(self.properties_tab)
        self.properties_layout.addWidget(self.molecularWeightLabel)
        self.properties_layout.addWidget(self.monoisotopicMwLabel)

        # create tab widget
        self.tab_widget = QTabWidget(self.central_widget)
        self.tab_widget.addTab(self.properties_tab, "Properties")

        # create main layout
        self.main_layout = QHBoxLayout(self.central_widget)
        self.main_layout.addWidget(self.canvas)
        self.main_layout.addWidget(self.tab_widget)

        # set central widget
        self.setCentralWidget(self.central_widget)

        # create menu bar
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("File")
        open_smiles_action = file_menu.addAction("Open SMILES")
        open_smiles_action.triggered.connect(self.open_smiles_dialog)

    def open_smiles_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Open SMILES")

        # create SMILES input box
        smiles_label = QLabel("Enter SMILES code:")
        smiles_textbox = QTextEdit(dialog)
        smiles_textbox.setFixedHeight(60)

        # create submit button
        submit_button = QPushButton("Submit", dialog)
        submit_button.clicked.connect(lambda: self.load_molecule(str(smiles_textbox.toPlainText())))

        # create dialog layout
        dialog_layout = QVBoxLayout(dialog)
        dialog_layout.addWidget(smiles_label)
        dialog_layout.addWidget(smiles_textbox)
        dialog_layout.addWidget(submit_button)

        # show dialog
        dialog.exec_()

    def draw_molecule(self, molecule):
        AllChem.Compute2DCoords(molecule)
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg_data = QByteArray(svg.encode('utf-8'))
        svg_renderer = QSvgRenderer(svg_data)
        img = QImage(self.canvas.size(), QImage.Format_ARGB32)
        img.fill(Qt.white)
        painter = QPainter(img)
        svg_renderer.render(painter)
        painter.end()
        pixmap = QPixmap.fromImage(img)
        self.canvas.setPixmap(pixmap)

    def load_molecule(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        #molecule = Chem.AddHs(molecule)
        self.currentMolecule = molecule
        self.molecularWeightLabel.setText("Average mass: " + str(round(Descriptors.MolWt(molecule), 2)))
        self.monoisotopicMwLabel.setText("Monoisotopic mass: " + str(round(Descriptors.ExactMolWt(molecule), 2)))
        self.draw_molecule(molecule)

if __name__ == "__main__":
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    viewer.load_molecule("C(CC(=O)NC(CS)C(=O)NCC(=O)O)C(C(=O)O)N")
    app.exec_()
