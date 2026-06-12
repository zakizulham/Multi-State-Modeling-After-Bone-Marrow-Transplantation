# Multi-State Modeling After Bone Marrow Transplantation: A Continuous-Time Markov Chain Approach

This repository contains the complete, reproducible computational pipeline and LaTeX source documents for our multi-state stochastic analysis on leukemia patient trajectories post-Bone Marrow Transplantation (BMT). This study utilizes a 6-state time-homogeneous Continuous-Time Markov Chain (CTMC) framework to evaluate intermediate transient phases and long-term competing risk absorption allocations.

## 📂 Repository Structure

* `analysis.R` : Core R script containing data restructuring via `msprep`, parametric Maximum Likelihood Estimation (MLE) of the $\mathbf{Q}$-matrix, and matrix exponential transient projections.
* `main.qmd` : Quarto Markdown document rendering the dynamic notebook environment and computational outputs.
* `main.tex` : The final peer-reviewed formatted LaTeX document using the `rho-class` template.
* `rho.cls` : The LaTeX class template file required for document layout compilation.
* `rho.bib` : BibTeX database containing all foundational biostatistical and stochastic references.

## 🛠️ Dependencies and Package Details

Since this project does not employ an isolated `renv` environment, please ensure you have the following packages installed in your global R ecosystem before execution. The pipeline relies on these specific core libraries:

| Package | Purpose in this Project | Key Functions Used |
| :--- | :--- | :--- |
| **mstate** | Multi-state data preparation and long-format conversion | `msprep()`, `transMat()` |
| **survival** | Foundational survival analysis infrastructure | `Surv()` |
| **Matrix** | High-performance matrix operations | `expm()` |
| **tidyverse** | Data manipulation, pipeline piping, and engineering | `dplyr`, `purrr`, `tidyr` |
| **ggplot2** | Continuous state occupancy probability visualizations | `ggplot()`, `geom_line()` |

You can install all required packages simultaneously by running the following command in your R console:
```R
install.packages(c("mstate", "survival", "Matrix", "tidyverse"))
```

## 🚀 How to Reproduce

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/zakizulham/Multi-State-Modeling-After-Bone-Marrow-Transplantation](https://github.com/zakizulham/Multi-State-Modeling-After-Bone-Marrow-Transplantation)
    ```
2.  **Run the Analysis:** Open RStudio, set the repository as your working directory, and execute `analysis.R` or render `main.qmd` to generate the empirical parameters and probability trajectory plots.
3.  **Compile the Paper:** Compile `main.tex` via any standard TeX engine (or Overleaf) with `rho.bib` and `rho.cls` attached to generate the final two-column PDF report.
4. `main.tex` : The final peer-reviewed formatted LaTeX document using the `rho-class` template.
5. `main.pdf` : The compiled, ready-to-read PDF version of the final research paper.
6. `rho.cls`  : The LaTeX class template file required for document layout compilation.

## 📜 Attribution and Citation

If you utilize this analytical framework, mathematical configurations, or code pipelines for your academic research or actuarial applications, please provide credit by citing this repository:

```text
[Muhammad Zaki Zulhamlizar], [Louisa Devina], [Andhika Dzaky Karafathi], & [Muhammad Salman]. (2026). 
Multi-State Modeling After Bone Marrow Transplantation: A Continuous-Time Markov Chain Approach. 
GitHub Repository: [https://github.com/zakizulham/Multi-State-Modeling-After-Bone-Marrow-Transplantation](https://github.com/zakizulham/Multi-State-Modeling-After-Bone-Marrow-Transplantation)
```

## 👥 Contributors

This research project was conducted as part of the foundational Stochastic Models curriculum. We thank all group members for their equal contributions to data validation, mathematical proofs, and computational implementation:

* **[Muhammad Zaki Zulhamlizar]** ([zakizulham205@gmail.com]) - Core Stochastic Modeling & Computation Pipeline.
* **[Muhammad Salman]** - Data Validation & Verification.
* **[Andhika Dzaky Karafathi]** - LaTeX Document Engineering & Formatting.
* **[Louisa Devina]** - Clinical Interpretation & Asset Review.

---
*This project is dedicated to open-science principles and reproducible mathematical modeling.*