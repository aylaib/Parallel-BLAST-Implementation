# Parallel BLAST Algorithm Implementation in Python
This project, carried out as part of the "Architecture and Parallel Computing" module, presents the implementation and performance analysis of a parallelized version of the BLAST algorithm (Basic Local Alignment Search Tool).
The main objective is to accelerate the execution of BLAST, an essential algorithm in bioinformatics, by exploiting multi-core processor architectures using Python's `multiprocessing` library.
## 📊 Performance Analysis
A comparative analysis was conducted between the sequential and parallel versions of the algorithm. The results show a **significant speedup of 2.90x** on 8 processes, with an efficiency of 36.22%. The most expensive step, alignment extension, was accelerated by **4.33x**.
### Time and Space Complexity (Sequential Version)
The analysis of the sequential version confirmed a time complexity close to **O(n²)** and a space complexity close to **O(n)**, which is consistent with BLAST theory.
![Complexity Analysis](./results_images/blast_complexity_analysis.png)
*Graphs showing the time complexity (left) and space complexity (right) of the sequential implementation.*
## 🛠️ Parallelization Model
After a profiling analysis (`cProfile`) to identify bottlenecks, a 3-stage pipeline model was designed and implemented:
1.  **Index Creation (25% of cores):** The "subject" sequence is divided into `chunks` and the word index is built in parallel.
2.  **Seed Search (25% of cores):** The "query" sequence is also divided to search for initial matches (seeds) in the index in parallel.
3.  **Alignment Extension (50% of cores):** This is the most computationally expensive step. The found seeds are distributed among multiple processes to extend and score alignments simultaneously.
Communication between processes is managed by a queue system (`multiprocessing.Queue`) to ensure efficient data flow.
## 🚀 Scripts and Usage
The project is structured into several scripts for implementation, analysis and comparison.
- **`src/sequential_blast.py`** : The reference implementation of BLAST, purely sequential.
- **`src/parallel_blast.py`** : The parallel implementation that compares performance with the sequential version.
- **`src/profiling_analysis.py`** : Script used to profile the sequential version and identify bottlenecks.
- **`src/blast-sequentiel-temp.py` & `src/blast-sequentiel-spat.py`** : Scripts dedicated to empirical complexity analysis.
### How to Run
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/YOUR_USERNAME/Parallel-BLAST-Implementation.git
    cd Parallel-BLAST-Implementation
    ```
2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
3.  **Run the main comparison script:**
    ```bash
    python src/parallel_blast.py
    ```
    The script will generate random sequences and display a detailed comparison of execution times between sequential and parallel versions.

    ## 📌 Citation
If you use this project, please cite it as:

Ayoub Laib (2025), *Parallel BLAST Algorithm Implementation in Python*, GitHub repository: https://github.com/aylaib/Parallel-BLAST-Implementation


## 📚 Reference Documents
- **[Complete Project Report](./Rapport_Paralleisation_BLAST.pdf)** : This document contains the detailed presentation of BLAST, complexity analysis, parallelization model, performance results and conclusion.
- **[Project Statement](./Enonce_Projet_ACP.pdf)** : The original specifications.
