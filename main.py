from rdkit import Chem
from rdkit.Chem import rdchem
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

class Molecule(rdchem.Mol):
    def __init__(self, mol):
        rdchem.Mol.__init__(self, mol)
        self.monoisotopicMass = Descriptors.ExactMolWt(mol)
        self.averageMass = Descriptors.MolWt(mol)
        self.formula = Chem.rdMolDescriptors.CalcMolFormula(mol, False)
        self.m1_pve = mass.calculate_mass(formula = self.formula, ion_type="M", charge = 1)
        self.m2_pve = mass.calculate_mass(formula = self.formula, ion_type="M", charge = 2)
        self.m3_pve = mass.calculate_mass(formula = self.formula, ion_type="M", charge = 3)
        self.mNa_pve = mass.calculate_mass(formula = self.formula, charge_carrier="Na+", ion_type="M", charge = 1)
        self.mK_pve = mass.calculate_mass(formula = self.formula, charge_carrier="K+", ion_type="M", charge = 1)
        self.m1_nve = -mass.calculate_mass(formula = self.formula, ion_type="M", charge = -1)
        self.m2_nve = -mass.calculate_mass(formula = self.formula, ion_type="M", charge = -2)
        self.m3_nve = -mass.calculate_mass(formula = self.formula, ion_type="M", charge = -3)


class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        # set window title
        self.setWindowTitle("MolView")
        self.resize(800, 600)

        self.currentMolecule = Molecule(mol=Chem.MolFromSmiles("CCCCCC"))

        # create central widget
        self.central_widget = QWidget(self)

        # create 2D canvas
        self.canvas = QLabel(self.central_widget)
        self.canvas.setAlignment(Qt.AlignTop)
        self.canvas.setFixedSize(600, 600)

        # create layout
        self.display_layout = QVBoxLayout()
        self.display_layout.addWidget(self.canvas)

        # set white background
        pal = QPalette()
        pal.setColor(QPalette.Background, QColor(255, 255, 255))
        self.canvas.setAutoFillBackground(True)
        self.canvas.setPalette(pal)

        # create properties tab labels for molecular properties
        self.molecularWeightLabel = QLabel("<b>Average mass:</b> ")
        self.monoisotopicMwLabel = QLabel("<b>Monoisotopic mass:</b> ")
        self.space1Label = QLabel("")
        self.pve_charge_Label = QLabel("<u><b>Positive Charge States</b></u>")
        self.m1_pve_Label = QLabel("[M+H]<sup>+</sup>: ")
        self.mNa_pve_Label = QLabel("[M+Na]<sup>+<sup/>: ")
        self.mK_pve_label = QLabel("[M+K]<sup>+<sup/>: ")
        self.space2Label = QLabel("")
        self.m2_pve_Label = QLabel("[M+2H]<sup>2+</sup>: ")
        self.m3_pve_Label = QLabel("[M+3H]<sup>3+</sup>: ")
        self.space3Label = QLabel("")
        self.nve_charge_Label = QLabel("<u><b>Negative Charge States</b></u>")
        self.m1_nve_Label = QLabel("[M-H]<sup>-</sup>: ")
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
        self.properties_layout.addWidget(self.mNa_pve_Label)
        self.properties_layout.addWidget(self.mK_pve_label)
        self.properties_layout.addWidget(self.space2Label)
        self.properties_layout.addWidget(self.m2_pve_Label)
        self.properties_layout.addWidget(self.m3_pve_Label)
        self.properties_layout.addWidget(self.space3Label)
        self.properties_layout.addWidget(self.nve_charge_Label)
        self.properties_layout.addWidget(self.m1_nve_Label)
        self.properties_layout.addWidget(self.m2_nve_Label)
        self.properties_layout.addWidget(self.m3_nve_Label)

        # create tab widget
        self.tab_widget = QTabWidget(self.central_widget)
        self.tab_widget.addTab(self.properties_tab, "Properties")

        # create main layout
        self.main_layout = QHBoxLayout(self.central_widget)
        self.main_layout.addLayout(self.display_layout)
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
            newMolecule = Chem.MolFromMolBlock(molblock)
            self.load_molecule(newMolecule)

    def draw_molecule2D(self, molecule):
        AllChem.Compute2DCoords(molecule)
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
        drawer.ClearDrawing()
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

        self.currentMolecule.__init__(mol)
        #self.currentMolecule.setProperties(mol)
        self.molecularWeightLabel.setText("<b>Average mass:</b> " + str(round(self.currentMolecule.averageMass, 4)))
        self.monoisotopicMwLabel.setText("<b>Monoisotopic mass:</b> " + str(round(self.currentMolecule.monoisotopicMass, 4)))
        self.m1_pve_Label.setText("[M+H]<sup>+</sup>: " + str(round(self.currentMolecule.m1_pve, 4)))
        self.m2_pve_Label.setText("[M+2H]<sup>2+</sup>: " + str(round(self.currentMolecule.m2_pve, 4)))
        self.m3_pve_Label.setText("[M+3H]<sup>3+</sup>: " + str(round(self.currentMolecule.m3_pve, 4)))
        self.mNa_pve_Label.setText("[M+Na]<sup>+</sup>: " + str(round(self.currentMolecule.mNa_pve, 4)))
        self.mK_pve_label.setText("[M+K]<sup>+</sup>: " + str(round(self.currentMolecule.mK_pve, 4)))
        self.m1_nve_Label.setText("[M-H]<sup>-</sup>: " + str(round(-1*self.currentMolecule.m1_nve, 4)))
        self.m2_nve_Label.setText("[M-2H]<sup>2-</sup>: " + str(round(-1*self.currentMolecule.m2_nve, 4)))
        self.m3_nve_Label.setText("[M-3H]<sup>3-</sup>: " + str(round(-1*self.currentMolecule.m3_nve, 4)))
        self.draw_molecule2D(self.currentMolecule)

if __name__ == "__main__":
    app = QApplication([])
    viewer = MoleculeViewer()
    viewer.show()
    viewer.load_molecule(Chem.MolFromSmiles("CC(C(=O)NC(CCC(=O)NC(CCCC(C(=O)O)N)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)NC(=O)C(C)OC1C(C(OC(C1O)CO)OP(=O)(O)OP(=O)(O)OCC2C(C(C(O2)N3C=CC(=O)NC3=O)O)O)NC(=O)C"))
    app.exec_()
