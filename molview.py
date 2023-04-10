from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors
from pyteomics import mass
from PyQt5.QtGui import QPixmap, QPalette, QColor, QImage, QPainter
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QFileDialog, QHBoxLayout, QLabel, QDialog, QTextEdit, QTabWidget
from PyQt5.QtSvg import QSvgWidget, QSvgRenderer
from PyQt5.QtCore import Qt, QByteArray


class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        self.currentMolecule = Chem.MolFromSmiles("CC(C(=O)NC(CCC(=O)NC(CCCC(C(=O)O)N)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)NC(=O)C(C)OC1C(C(OC(C1O)CO)OP(=O)(O)OP(=O)(O)OCC2C(C(C(O2)N3C=CC(=O)NC3=O)O)O)NC(=O)C")

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

        # create properties tab labels for molecular properties
        self.molecularWeightLabel = QLabel("Average mass:")
        self.monoisotopicMwLabel = QLabel("Monoisotopic mass: ")
        self.space1Label = QLabel("")
        self.pve_charge_Label = QLabel("<u>Positive Charge States</u>")
        self.m1_pve_Label = QLabel("[M+H]+: ")
        self.m2_pve_Label = QLabel("[M+2H]<sup>2+</sup>: ")
        self.m3_pve_Label = QLabel("[M+3H]<sup>3+</sup>: ")
        self.space2Label = QLabel("")
        self.nve_charge_Label = QLabel("<u>Negative Charge States</u>")
        self.m1_nve_Label = QLabel("[M-H]-: ")
        self.m2_nve_Label = QLabel("[M-2H]<sup>2-</sup>: ")
        self.m3_nve_Label = QLabel("[M-3H]<sup>3-</sup>: ")

        # set properties of properties tab
        self.properties_tab = QWidget(self.central_widget)
        self.properties_layout = QVBoxLayout(self.properties_tab)
        self.properties_layout.setAlignment(Qt.AlignTop)
        self.properties_layout.setSpacing(0)

        # add labels to layout
        self.properties_layout.addWidget(self.molecularWeightLabel)
        self.properties_layout.addWidget(self.monoisotopicMwLabel)
        self.properties_layout.addWidget(self.space1Label)
        self.properties_layout.addWidget(self.pve_charge_Label)
        self.properties_layout.addWidget(self.m1_pve_Label)
        self.properties_layout.addWidget(self.m2_pve_Label)
        self.properties_layout.addWidget(self.m3_pve_Label)
        self.properties_layout.addWidget(self.space2Label)
        self.properties_layout.addWidget(self.nve_charge_Label)
        self.properties_layout.addWidget(self.m1_nve_Label)
        self.properties_layout.addWidget(self.m2_nve_Label)
        self.properties_layout.addWidget(self.m3_nve_Label)

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

        # create submenu for opening files
        open_menu = QMenu("Open", self)
        file_menu.addMenu(open_menu)

        # create submenu opening options
        open_smiles_action = open_menu.addAction("Paste SMILES")
        open_smiles_action.triggered.connect(self.open_smiles_dialog)
        open_file_action = open_menu.addAction("Open File")
        open_file_action.triggered.connect(self.open_file_dialog)

    def open_smiles_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Open SMILES")

        # create SMILES input box
        smiles_label = QLabel("Enter SMILES code:")
        smiles_textbox = QTextEdit(dialog)
        smiles_textbox.setFixedHeight(60)

        # create submit button
        submit_button = QPushButton("Submit", dialog)
        submit_button.clicked.connect(lambda: (self.load_molecule(Chem.MolFromSmiles(str(smiles_textbox.toPlainText()))), dialog.accept()))

        # create dialog layout
        dialog_layout = QVBoxLayout(dialog)
        dialog_layout.addWidget(smiles_label)
        dialog_layout.addWidget(smiles_textbox)
        dialog_layout.addWidget(submit_button)

        # show dialog
        dialog.exec_()

    def open_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Open Molecule File", "", "Molecule Files (*.mol);;All Files (*)", options=options)
        if file_name:
            with open(file_name, "r") as f:
                molblock = f.read()
            molecule = Chem.MolFromMolBlock(molblock)
            self.load_molecule(molecule)

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

    def load_molecule(self, molecule):
        # clean up and uncharge molecule
        mol = rdMolStandardize.Cleanup(molecule)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)

        self.currentMolecule = mol
        self.molecularWeightLabel.setText("Average mass: " + str(round(Descriptors.MolWt(mol), 2)))
        self.monoisotopicMwLabel.setText("Monoisotopic mass: " + str(round(Descriptors.ExactMolWt(mol), 2)))
        self.m1_pve_Label.setText("[M+H]+: " + str(round(mass.calculate_mass(formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False), ion_type="M", charge = 1), 2)))
        self.m2_pve_Label.setText("[M+2H]<sup>2+</sup>: " + str(round(mass.calculate_mass(formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False), ion_type="M", charge = 2), 2)))
        self.m3_pve_Label.setText("[M+3H]<sup>3+</sup>: " + str(round(mass.calculate_mass(formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False), ion_type="M", charge = 3), 2)))
        self.m1_nve_Label.setText("[M-H]-: " + str(-round(mass.calculate_mass(formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False), ion_type="M", charge = -1), 2)))
        self.m2_nve_Label.setText("[M-2H]<sup>2-</sup>: " + str(-round(mass.calculate_mass(formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False), ion_type="M", charge = -2), 2)))
        self.m3_nve_Label.setText("[M-3H]<sup>3-</sup>: " + str(-round(mass.calculate_mass(formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False), ion_type="M", charge = -3), 2)))
        self.draw_molecule(mol)

if __name__ == "__main__":
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    viewer.load_molecule(viewer.currentMolecule)
    app.exec_()
