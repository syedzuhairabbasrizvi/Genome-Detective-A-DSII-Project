import sys
from PyQt6.QtWidgets import QComboBox, QApplication, QMainWindow, QLabel, QPushButton, QVBoxLayout, QWidget, QMessageBox, QLineEdit
import random
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt
from PyQt6.uic import loadUi
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
import re

class HashMap:
    def __init__(self, size=100):
        self.size = size
        self.map = [None] * self.size

    def _hash(self, key):
        return hash(key) % self.size

    def put(self, key, value):
        index = self._hash(key)
        if self.map[index] is None:
            self.map[index] = [(key, value)]
        else:
            for i, (k, _) in enumerate(self.map[index]):
                if k == key:
                    self.map[index][i] = (key, value)
                    return
            self.map[index].append((key, value))

    def get(self, key):
        index = self._hash(key)
        if self.map[index] is not None:
            for k, v in self.map[index]:
                if k == key:
                    return v
        return []

    def contains_key(self, key):
        index = self._hash(key)
        if self.map[index] is not None:
            for k, _ in self.map[index]:
                if k == key:
                    return True
        return False

    def remove(self, key):
        index = self._hash(key)
        if self.map[index] is not None:
            for i, (k, _) in enumerate(self.map[index]):
                if k == key:
                    del self.map[index][i]
                    return

    def keys(self):
        keys = []
        for bucket in self.map:
            if bucket:
                keys.extend([k for k, _ in bucket])
        return keys

    def values(self):
        values = []
        for bucket in self.map:
            if bucket:
                values.extend([v for _, v in bucket])
        return values

    def items(self):
        items = []
        for bucket in self.map:
            if bucket:
                items.extend(bucket)
        return items

color_map = HashMap()
# Define the codon-to-amino-acid mapping
codon_to_amino_acid_map = HashMap()

codon_to_amino_acid_map.put("TTT", "F")
codon_to_amino_acid_map.put("TTC", "F")
codon_to_amino_acid_map.put("TTA", "L")
codon_to_amino_acid_map.put("TTG", "L")

codon_to_amino_acid_map.put("CTT", "L")
codon_to_amino_acid_map.put("CTC", "L")
codon_to_amino_acid_map.put("CTA", "L")
codon_to_amino_acid_map.put("CTG", "L")

codon_to_amino_acid_map.put("ATT", "I")
codon_to_amino_acid_map.put("ATC", "I")
codon_to_amino_acid_map.put("ATA", "I")
codon_to_amino_acid_map.put("ATG", "M")

codon_to_amino_acid_map.put("GTT", "V")
codon_to_amino_acid_map.put("GTC", "V")
codon_to_amino_acid_map.put("GTA", "V")
codon_to_amino_acid_map.put("GTG", "V")

codon_to_amino_acid_map.put("TCT", "S")
codon_to_amino_acid_map.put("TCC", "S")
codon_to_amino_acid_map.put("TCA", "S")
codon_to_amino_acid_map.put("TCG", "S")

codon_to_amino_acid_map.put("CCT", "P")
codon_to_amino_acid_map.put("CCC", "P")
codon_to_amino_acid_map.put("CCA", "P")
codon_to_amino_acid_map.put("CCG", "P")

codon_to_amino_acid_map.put("ACT", "T")
codon_to_amino_acid_map.put("ACC", "T")
codon_to_amino_acid_map.put("ACA", "T")
codon_to_amino_acid_map.put("ACG", "T")

codon_to_amino_acid_map.put("GCT", "A")
codon_to_amino_acid_map.put("GCC", "A")
codon_to_amino_acid_map.put("GCA", "A")
codon_to_amino_acid_map.put("GCG", "A")

codon_to_amino_acid_map.put("TAT", "Y")
codon_to_amino_acid_map.put("TAC", "Y")
codon_to_amino_acid_map.put("TAA", "*")
codon_to_amino_acid_map.put("TAG", "*")

codon_to_amino_acid_map.put("CAT", "H")
codon_to_amino_acid_map.put("CAC", "H")
codon_to_amino_acid_map.put("CAA", "Q")
codon_to_amino_acid_map.put("CAG", "Q")

codon_to_amino_acid_map.put("AAT", "N")
codon_to_amino_acid_map.put("AAC", "N")
codon_to_amino_acid_map.put("AAA", "K")
codon_to_amino_acid_map.put("AAG", "K")

codon_to_amino_acid_map.put("GAT", "D")
codon_to_amino_acid_map.put("GAC", "D")
codon_to_amino_acid_map.put("GAA", "E")
codon_to_amino_acid_map.put("GAG", "E")

codon_to_amino_acid_map.put("TGT", "C")
codon_to_amino_acid_map.put("TGC", "C")
codon_to_amino_acid_map.put("TGA", "*")
codon_to_amino_acid_map.put("TGG", "W")

codon_to_amino_acid_map.put("CGT", "R")
codon_to_amino_acid_map.put("CGC", "R")
codon_to_amino_acid_map.put("CGA", "R")
codon_to_amino_acid_map.put("CGG", "R")

codon_to_amino_acid_map.put("AGT", "S")
codon_to_amino_acid_map.put("AGC", "S")
codon_to_amino_acid_map.put("AGA", "R")
codon_to_amino_acid_map.put("AGG", "R")

codon_to_amino_acid_map.put("GGT", "G")
codon_to_amino_acid_map.put("GGC", "G")
codon_to_amino_acid_map.put("GGA", "G")
codon_to_amino_acid_map.put("GGG", "G")


class MainPage(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("GenomeDetective")
        self.setGeometry(300, 100, 800, 500)  # Adjust the size of the main page

        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        layout = QVBoxLayout()
        main_widget.setLayout(layout)

        intro_widget = QWidget()
        intro_layout = QVBoxLayout()
        intro_widget.setLayout(intro_layout)
        intro_widget.setStyleSheet("background-color: rgba(255, 255, 255, 0.8); border-radius: 10px; padding: 20px;")

        intro_label = QLabel("<html><b>Welcome to GenomeDetective!</b><br><br>"
                             "GenomeDetective is a powerful software tool designed for gene sequence analysis. <br>"
                             "It provides three main functionalities:<br><br>\n\n"
                             "1. Compare Gene Sequences <br>\n"
                             "2. Disease Detection <br>\n"
                             "3. Sequence Logo <br>\n"
                             "Please select an option below to get started:<br>")
        #intro_label.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Center-align the text
        layout.addWidget(intro_label)
        layout.addSpacing(-50)
        font = QFont()
        font.setPointSize(15)  # Set font size
        font.setFamily("Arial")  # Set font family
        intro_label.setFont(font)

        button_compare = QPushButton("1. Compare Gene Sequences", self)
        button_compare.clicked.connect(self.open_compare)
        layout.addWidget(button_compare)

        button_detection = QPushButton("2. Disease Detection", self)
        button_detection.clicked.connect(self.open_detection)
        layout.addWidget(button_detection)

        button_logo = QPushButton("3. Sequence Logo", self)
        button_logo.clicked.connect(self.open_logo)
        layout.addWidget(button_logo)

        # Apply stylesheet to buttons
        self.setStyleSheet("""
                            QMainWindow {
                background-image: url(finallogo.png);
                background-position: top;
                background-repeat: no-repeat;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: 2px solid #4CAF50;
                border-radius: 10px;
                padding: 10px 20px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)


    def open_compare(self):
        self.compare_window = CompareWindow()
        self.compare_window.show()

    def open_detection(self):
        self.detection_window = DetectionWindow()
        self.detection_window.show()

    def open_logo(self):
        self.logo_window = LogoWindow()
        self.logo_window.show()

# Function to convert nucleotide sequence to amino acid sequence using hash map
def nucleotide_to_amino_acid(sequence, codon_to_amino_acid_map):
    amino_acid_sequence = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = codon_to_amino_acid_map.get(codon)
        if amino_acid != []:
            amino_acid_sequence += amino_acid
        else:
            amino_acid_sequence += "-"  # Handle unknown codons
    return amino_acid_sequence

# Function to calculate amino acid frequency
def calculate_amino_acid_frequency(codon_sequence):
    amino_acid_frequency = {}
    for codon, positions in codon_sequence.items():
        amino_acid_frequency[codon_to_amino_acid_map.get(codon)] = len(positions)
    return amino_acid_frequency

def similarity(seq1, seq2):
    num_codons = len(seq1)
    print(seq1, seq2)
    matches = sum(1 for i in range(num_codons) if seq1[i] == seq2[i])
    print(matches)
    similarity_score = matches / num_codons
    return similarity_score

# def process_sequence(sequence):
#     # Define a string containing all valid nucleotides
#     valid_nucleotides = 'ATCG'
        
#     # Convert the sequence to uppercase and remove spaces
#     sequence = sequence.upper().replace(" ", "")
        
#     # Filter out characters that are not valid nucleotides
#     processed_sequence = ''.join(char for char in sequence if char in valid_nucleotides)
    
#     return processed_sequence

def preprocess_sequence(sequence):
    # Convert to uppercase
    sequence = sequence.upper()
    # Remove spaces and digits
    sequence = re.sub(r'[0-9\s]+', '', sequence)
    return sequence

class CompareWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        loadUi("geneComparator.ui", self)
        self.setWindowTitle("Gene Sequence Comparator")
        self.setGeometry(300, 100, 800, 500) 
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        self.label1 = QLabel("Enter Gene Sequence 1:")
        self.layout.addWidget(self.label1)
        self.input_sequence1 = QLineEdit()
        self.layout.addWidget(self.input_sequence1)

        self.label2 = QLabel("Enter Gene Sequence 2:")
        self.layout.addWidget(self.label2)
        self.input_sequence2 = QLineEdit()
        self.layout.addWidget(self.input_sequence2)

        font = QFont()
        font.setPointSize(15)  # Set font size
        font.setFamily("Sans Serif")  # Set font family
        self.label1.setFont(font)
        self.label2.setFont(font)
        self.input_sequence1.setFont(font)
        self.input_sequence2.setFont(font)

        self.compare_button = QPushButton("Compare Sequences")
        self.compare_button.clicked.connect(self.compare_sequences)
        self.layout.addWidget(self.compare_button)

        # back button
        self.back_button = QPushButton("Back")
        self.back_button.clicked.connect(self.go_back)
        self.layout.addWidget(self.back_button)

        self.setStyleSheet("""
                           QMainWindow {
                background-image: url(DNA2.png);
                background-position: bottom;
                background-repeat: no-repeat;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: 2px solid #4CAF50;
                border-radius: 10px;
                padding: 10px 20px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #45a049;
                border: 2px solid #45a049;
            }
        """)

        self.result_label = QLabel("")
        self.layout.addWidget(self.result_label)

    def go_back(self):
        self.close()

    def compare_sequences(self):
        sequence1 = self.input_sequence1.text()
        sequence2 = self.input_sequence2.text()
        if not sequence1 or not sequence2:
            QMessageBox.warning(self, "Warning", "Please provide an input sequence.")
            return  
        sequence1 = preprocess_sequence(sequence1)
        sequence2 = preprocess_sequence(sequence2)
        if not all(base in ['A', 'T', 'G', 'C'] for base in sequence1) or not all(base in ['A', 'T', 'G', 'C'] for base in sequence2):
            QMessageBox.warning(self, "Warning", "DNA sequence should only contain the bases A, T, G, C!")
            return  
        if len(sequence1) != len(sequence2) or len(sequence1) % 3 != 0:
            QMessageBox.warning(self, "Warning", "Sequences must have the same length and be divisible by 3 (representing codons)")
            return
        # Convert DNA sequences to amino acid sequences using the codon-to-amino-acid map
        amino_acid_sequence1 = nucleotide_to_amino_acid(sequence1, codon_to_amino_acid_map)
        amino_acid_sequence2 = nucleotide_to_amino_acid(sequence2, codon_to_amino_acid_map)

        # Calculate similarity score
        try:
            similarity_score = similarity(amino_acid_sequence1, amino_acid_sequence2)
            self.result_label.setText(f"Similarity Score: {similarity_score}")

        except ValueError as e:
            QMessageBox.warning(self, "Error", str(e))
            self.result_label.clear()

def get_codon_indices(sequence):
    codon_indices_map = HashMap()
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon_indices_map.contains_key(codon):
            codon_indices_map.get(codon).append(i//3)
        else:
            codon_indices_map.put(codon, [i//3])

    return codon_indices_map

def is_color_too_bright(color):
    # Convert color string to RGB tuple
    rgb = mcolors.to_rgba(color)[:3]
    # Calculate luminance using a standard formula (0.2126*R + 0.7152*G + 0.0722*B)
    luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    # Define a threshold for brightness
    brightness_threshold = 0.9
    # Return True if luminance exceeds the threshold, indicating a too bright color
    return luminance > brightness_threshold

class DetectionWindow(QMainWindow):
    def __init__(self):
        super(DetectionWindow, self).__init__()

        self.setWindowTitle("Disease Detection")
        self.setGeometry(300, 100, 800, 500) 

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        layout = QVBoxLayout()
        self.central_widget.setLayout(layout)

        # Title label
        title_label = QLabel("Disease Detection")
        title_font = QFont()
        title_font.setPointSize(20)
        title_font.setBold(True)
        title_label.setFont(title_font)
        layout.addWidget(title_label)

        # DNA sequence input
        sequence_label = QLabel("Enter Gene Sequence:")
        layout.addWidget(sequence_label)
        self.input_sequence = QLineEdit()
        layout.addWidget(self.input_sequence)

        # Disease selection combobox
        disease_label = QLabel("Select Disease:")
        layout.addWidget(disease_label)
        self.comboBox = QComboBox()
        self.populate_combobox()
        layout.addWidget(self.comboBox)

       

        # Detect Disease button
        self.go_button = QPushButton("Detect Disease")
        self.go_button.clicked.connect(self.sequence_checker)
        layout.addWidget(self.go_button)

        # Refresh button
        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.clicked.connect(self.on_refresh_click)
        layout.addWidget(self.refresh_button)

        # back button
        self.back_button = QPushButton("Back")
        self.back_button.clicked.connect(self.go_back)
        layout.addWidget(self.back_button)

         # Result label
        self.result_label = QLabel()
        #self.result_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.result_label.setObjectName("resultLabel")
        layout.addWidget(self.result_label)

        # Apply stylesheet to buttons and labels
        self.setStyleSheet("""
            QLabel {
                font-size: 16px;
                font-weight: bold;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: 2px solid #4CAF50;
                border-radius: 10px;
                padding: 10px 20px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #45a049;
                border: 2px solid #45a049;
            }
            #resultLabel {
                background-color: #FFFF00;
                padding: 5px;
                border-radius: 5px;
            }
        """)

    def go_back(self):
        self.close()
        
    def sequence_checker(self):
        input_sequence = self.input_sequence.text()
        
        if len(input_sequence) < 200:
            QMessageBox.warning(self, "Warning", "Gene sequence is too short!")
            return 
        # Preprocess the input sequence
        processed_sequence = preprocess_sequence(input_sequence)

        if input_sequence == "":
            QMessageBox.warning(self, "Warning", "Please enter a gene sequence!")
            return  
        
        if not all(base in ['A', 'T', 'G', 'C'] for base in processed_sequence) or processed_sequence == "":
            QMessageBox.warning(self, "Warning", "Gene sequence should only contain the bases A, T, G, C!")
            return  
        
        else:
            self.disease_decider()
    
    def populate_combobox(self):
        # Populate the combobox with disease options
        diseases = ["Sickle Cell Anemia", "Huntington's Disease", "Cystic Fibrosis", "Fragile X Syndrome","Duchenne Muscular Dystrophy","Tay-Sachs Disease", "Familial Hypercholesterolemia", "Beta-Thalassemia", "Hemochromatosis", "Alkaptonuria"]
        for disease in diseases:
            self.comboBox.addItem(disease)

    def disease_decider(self):
        # Get selected disease from combobox
        selected_disease = self.comboBox.currentText()

        # Get input sequence from line edit
        input_sequence = self.input_sequence.text()

        # Preprocess the input sequence
        processed_sequence = preprocess_sequence(input_sequence)
        
        if not all(base in ['A', 'T', 'G', 'C'] for base in processed_sequence):
            QMessageBox.warning(self, "Warning", "DNA sequence should only contain the bases A, T, G, C!")
            return  
        # if len(processed_sequence) % 3 != 0:
        #     QMessageBox.warning(self, "Warning", "Sequence must be divisible by 3 (representing codons)")
        #     return
        # Convert the input sequence to amino acid sequence
        codon_indices_map = get_codon_indices(processed_sequence)

        # Call respective functions based on selected disease
        if selected_disease == "Sickle Cell Anemia":
            self.function_for_sickle_cell_anemia(processed_sequence)
        elif selected_disease == "Huntington's Disease":
            self.function_for_huntingtons(codon_indices_map)
        elif selected_disease == "Fragile X Syndrome":
            self.function_for_fragilexsyndrome(codon_indices_map)
        elif selected_disease == "Duchenne Muscular Dystrophy":
            self.function_for_Duchenne_Muscular_Dystrophy(processed_sequence)
        elif selected_disease == "Tay-Sachs Disease":
            self.function_for_Tay_Sachs_Disease(processed_sequence)
        elif selected_disease == "Familial Hypercholesterolemia":
            self.function_for_FH(processed_sequence)
        elif selected_disease == "Beta-Thalassemia":
            self.function_for_Beta_Thalassemia(processed_sequence)
        elif selected_disease == "Cystic Fibrosis":
            self.function_for_cysticfibrosis(codon_indices_map)
        elif selected_disease == "Alkaptonuria":
            self.function_for_alkaptonuria(codon_indices_map)
        elif selected_disease == "Hemochromatosis":
            self.function_for_hemochromatosis(codon_indices_map)

    def function_for_sickle_cell_anemia(self, seq):
        # Check if the 6th codon is 'T'
        if 'T' in seq[5*3:5*3+3]:
            self.result_label.setText("Person has sickle cell anemia")
            return

        self.result_label.setText("Person is healthy")

    def function_for_huntingtons(self, codon_indices_map):
        # Count the number of 'CAG' repeats
        cag_repeats = len(codon_indices_map.get('CAG'))
        if cag_repeats >= 36:  # considering each 'CAG' represents 1 repeat
            self.result_label.setText("Person is at risk of Huntington's disease")
            return

        self.result_label.setText("Person is not at risk of Huntington's disease")

    def function_for_fragilexsyndrome(self, codon_indices_map):
        # Count the number of 'CGG' repeats
        cgg_repeats = len(codon_indices_map.get('CGG'))
        if cgg_repeats >= 200:  # considering each 'CGG' represents 1 repeat
            self.result_label.setText("Person is affected by Fragile X syndrome")
            return

        self.result_label.setText("Person is not affected by Fragile X syndrome")

    def function_for_cysticfibrosis(self, codon_indices_map):
        # Check the codon at position 508
        codon = codon_indices_map.get('TTT')
        codon2 = codon_indices_map.get('TTC')
        codon3 = codon_indices_map.get('TTA')
        # print(codon, codon2, codon3)
        if codon == None and codon2 == None:
            self.result_label.setText("Person does not have cystic fibrosis")
            return
        if 508 in codon or 508 in codon2 or 508 in codon3:
            self.result_label.setText("Person has cystic fibrosis")
            return
        else:
            self.result_label.setText("Person is healthy")
    def function_for_Duchenne_Muscular_Dystrophy(self, seq):
        # Check if the 6th codon is 'T'
        if 'C' == seq[3355]:
            self.result_label.setText("Person has Duchene Muscular Dystrophy")
            return

        self.result_label.setText("Person is healthy")
    def function_for_Tay_Sachs_Disease(self, seq):
        # Check if the 6th codon is 'T'
        if "TATC"== seq[1277:1281]:
            self.result_label.setText("Person has Tay-Sachs Disease")
            return

        self.result_label.setText("Person is healthy")

    def function_for_FH(self, sequence):
        print(sequence)
        if 'A' == sequence[680] or sequence[1476] == 'T': 
            self.result_label.setText("Exonic mutation detected. Person is at risk of Familial Hypercholesterolemia")
        elif sequence[312] == 'A':
            self.result_label.setText("Intronic mutation detected. Person is at risk of Familial Hypercholesterolemia")
        else:
            self.result_label.setText("Person is healthy")

    def function_for_Beta_Thalassemia(self, sequence):
        if 'C' == sequence[5]: 
            self.result_label.setText("Beta+ IVS-I no 6 mutation detected. Person is at risk of Beta-Thalassemia")
        elif sequence[100] == 'T':
            self.result_label.setText("Beta+ 101 C→T mutation detected. Person is at risk of Beta-Thalassemia")
        else:
            self.result_label.setText("Person is healthy")

    def function_for_hemochromatosis(self,codons):
        # C282Y mutation cysteine (C) becomes tyrosine (Y) , G replaced by A
        codon = codons.get('TAT')
        codon2 = codons.get('TAC')
        print(codons.keys)
        # print(codons[281])
        if 281 in codon or 281 in codon2: 
            self.result_label.setText("Person is at risk of Hemochromatosis disease")
        else:
            self.result_label.setText("Person is not at risk of Hemochromatosis disease")
            
    def function_for_alkaptonuria(self,codons):
        codon = codons.get('CGA')
        codon2 = codons.get('AGA')
        codon3 = codons.get('AGG')
        # arginine_codons_with_adenine = ['CGA', 'CGA', 'AGG']
            #the glycine residue at position 537 is replaced with arginine, sub of guanine with cytosine 
        # print(codons[536])
        if 536 in codon or 536 in codon2 or 536 in codon3:
            self.result_label.setText("Person is at risk of Alkaptonuria disease")
        else:
            self.result_label.setText("Person is not at risk of Alkaptonuria disease")



    def on_refresh_click(self):
        # Clear input sequence and result label
        self.input_sequence.clear()
        self.result_label.clear()

# Function to compare the input sequence with disease sequences and find matches
def find_matches(input_sequence, special_sequences):
    matches = []
    for special_sequence in special_sequences:
        index = 0
        while index < len(input_sequence):
            index = input_sequence.find(special_sequence, index)
            if index == -1:
                break
            matches.append((index, index + len(special_sequence), special_sequence))
            index += 1
    return matches

class LogoWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.codon_to_amino_acid_map = codon_to_amino_acid_map
        self.color_map = color_map
        self.setWindowTitle('Sequence Logo')
        self.setGeometry(300, 100, 800, 500) 
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        self.input_label = QLabel('Enter DNA Sequence:')
        self.input_textbox = QLineEdit()
        self.layout.addWidget(self.input_label)
        self.layout.addWidget(self.input_textbox)
        font = QFont()
        font.setPointSize(15)  # Set font size
        font.setFamily("Arial")  # Set font family
        self.input_label.setFont(font)
        self.input_textbox.setFont(font)
        self.frequency_button = QPushButton('Check for Amino Acid Frequency')
        self.frequency_button.clicked.connect(self.check_amino_acid_frequency)
        self.layout.addWidget(self.frequency_button)

        self.disease_button = QPushButton('Check for Important Sequences')
        self.disease_button.clicked.connect(self.check_for_special_sequences)
        self.layout.addWidget(self.disease_button)

        self.back_button = QPushButton("Back")
        self.back_button.clicked.connect(self.go_back)
        self.layout.addWidget(self.back_button)

        self.sequence_logo_label = QLabel()
        self.layout.addWidget(self.sequence_logo_label)

        self.setStyleSheet("""
                           QMainWindow {
                background-image: url(DNAart.png);
                background-position: top;
                background-repeat: no-repeat;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: 2px solid #4CAF50;
                border-radius: 10px;
                padding: 10px 20px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #45a049;
                border: 2px solid #45a049;
            }

        """)
    def go_back(self):
        self.close()

    def check_amino_acid_frequency(self):
        sequence = self.input_textbox.text()
        input_sequence = preprocess_sequence(sequence)
        if not all(base in ['A', 'T', 'G', 'C'] for base in input_sequence):
            QMessageBox.warning(self, "Warning", "DNA sequence should only contain the bases A, T, G, C!")
            return  
        if len(input_sequence) % 3 != 0:
            QMessageBox.warning(self, "Warning", "Sequence must be divisible by 3 (representing codons)")
            return
        codon_indices = get_codon_indices(input_sequence)
        amino_acid_frequency = calculate_amino_acid_frequency(codon_indices)
        self.show_amino_acid_frequency(amino_acid_frequency)

    def show_amino_acid_frequency(self, amino_acid_frequency):
        fig, ax = plt.subplots()
        amino_acids = amino_acid_frequency.keys()
        frequency = amino_acid_frequency.values()
        ax.bar(amino_acids, frequency)
        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Frequency')
        ax.set_title('Amino Acid Frequency')
        plt.show()

    def check_for_special_sequences(self):
        sequence = self.input_textbox.text()
        input_sequence = preprocess_sequence(sequence)
        if not all(base in ['A', 'T', 'G', 'C'] for base in input_sequence):
            QMessageBox.warning(self, "Warning", "DNA sequence should only contain the bases A, T, G, C!")
            return  
        amino_acid_sequence = nucleotide_to_amino_acid(input_sequence, self.codon_to_amino_acid_map)
        # Example disease sequences (amino acid)
        special_sequences_amino_acid = [
            "MLR", "GED", "FPN", "KHV", "YSW",
            "CNG", "RDE", "LIV", "HQT", "SGR",
            "NPG", "TVN", "PAS", "QYC", "DER",
            "WMH", "VLI", "AGP", "CMW", "GAV"
        ]

        # Find matches
        matches = find_matches(amino_acid_sequence, special_sequences_amino_acid)
        if matches:
            self.show_sequence_logo(amino_acid_sequence, matches)
        else:
            self.sequence_logo_label.setText("No matches found.")

    def show_sequence_logo(self, amino_acid_sequence, matches=None):
        fig, ax = plt.subplots()
        positions = range(1, len(amino_acid_sequence) + 1)

        # Plot bars for each position
        for i, aa in enumerate(amino_acid_sequence):
            ax.text(positions[i], 0, aa, ha='center', va='bottom')

        # Highlight matches with randomly chosen colors
        if matches:
            for start, end, disease_sequence in matches:
                color = self.get_disease_color(disease_sequence)
                ax.add_patch(Rectangle((start+0.5, 0), end-start, 1, color=color, alpha=0.7))

            # Create legend for color key
            legend_handles = []
            for disease_sequence, color in self.color_map.items():
                patch = Rectangle((0, 0), 1, 1, color=color, label=disease_sequence)
                legend_handles.append(patch)

            ax.legend(handles=legend_handles, loc='upper left')

        ax.set_xlim(0, len(amino_acid_sequence) + 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('Position')
        ax.set_ylabel('Frequency')
        ax.set_title('Sequence Logo')

        plt.show()

    
    def get_disease_color(self, disease_sequence):
    # Check if color for this disease sequence already exists
        if self.color_map.contains_key(disease_sequence):
            return self.color_map.get(disease_sequence)
        else:
            all_colors = [color for color in mcolors.CSS4_COLORS.values() if not is_color_too_bright(color)]
            # Keep choosing a random color until it's different enough from existing colors
            while True:
                random_color = random.choice(all_colors)
                if random_color not in self.color_map.values():
                    self.color_map.put(disease_sequence, random_color)
                    return random_color

def main():
    app = QApplication(sys.argv)
    main_window = MainPage()
    main_window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
