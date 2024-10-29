This project of mine, "GenomeDetective", concerns bioinformatics, specifically disease prediction.
![image](https://github.com/user-attachments/assets/3c60f076-968c-4f21-a78c-db58a626015b)

Working with Rubab Shah, Yashfeen Zaidi, and Fatima Dossa, I wrote a tool for our CS 201 - Data Structures II project. For data input, we utilized the NCBI database and extracted segments of gene sequences. Think of them as strings thousands of letters long and consisting only of A, C, T, G to represent each nucleotide. The project is written in Python and utilizes the hash map as the primary data structure to represent gene sequences via codon frequency. It is split into three functions: 

-> Compare gene sequences and generate a similarity index based on codon matching
![image](https://github.com/user-attachments/assets/a8e88ba0-3064-41fd-8fc5-c2e139b288bf)
-> Check for mutations leading to diseases; the information was taken from papers in our literature review
![image](https://github.com/user-attachments/assets/6fba9c8d-c444-4f9f-80f9-a5d267f7f8d7)
-> Generate sequence logos via Matplotlib based on amino acid frequency to visualize gene sequences at their most essential, though without order
![image](https://github.com/user-attachments/assets/6fb2bae5-3c4e-445e-859c-31930cbd47fa)
